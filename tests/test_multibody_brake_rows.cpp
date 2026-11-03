// Brake rows (docs/WHEEL.md, "M3c status"): a dof's friction loss, the reference's brake, as a row of the one soft-constraint
// law and solve, and its exact limit. The reference row against its source's formulas; the reference's own rotor and
// carrier controls (a hinge rotor through reversal, near rest, hold, over capacity, release and capacity changes; a free
// carrier and rotor with an internal drive: equal and opposite moments, momentum, relative-coordinate work and its
// refinement, the brake's moment at the rotor's joint); and the exact brake, which holds zero relative spin below its
// capacity and slips at its capacity above it, where the reference creeps at drive R / B.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include <doctest.h>

#include <sub0math/units.hpp>

#include "sub0physics/multibody.hpp"
#include "test_multibody_models.hpp"

using namespace sub0physics;
using namespace sub0physics::multibody;
namespace f64 = sub0math::f64;

namespace
{

// A hinge rotor of the given axial inertia, no gravity (the reference's fixture: a body on a hinge, diaginertia .4).
Model rotor(double inertia)
{
    using testmodels::diagonalInertia;
    ModelBuilder builder;
    builder.addBody({.name = "rotor",
                     .mass = sub0math::kilograms(1.0),
                     .centreOfMass = {},
                     .inertia = diagonalInertia(inertia, inertia, inertia),
                     .joint = {.type = JointType::Revolute, .axis = {0.0, 0.0, 1.0}}});
    auto model = builder.build();
    REQUIRE(model.has_value());
    return std::move(*model);
}

// A free carrier (inertia 2) and a rotor (0.5) on it, the same axis: the reference's carrier and rotor pair.
Model carrierAndRotor()
{
    using testmodels::diagonalInertia;
    ModelBuilder builder;
    const BodyId carrier = builder.addBody({.name = "carrier",
                                            .mass = sub0math::kilograms(1.0),
                     .centreOfMass = {},
                                            .inertia = diagonalInertia(2.0, 2.0, 2.0),
                                            .joint = {.type = JointType::Revolute, .axis = {0.0, 0.0, 1.0}}});
    builder.addBody({.name = "rotor",
                     .mass = sub0math::kilograms(1.0),
                     .centreOfMass = {},
                     .inertia = diagonalInertia(0.5, 0.5, 0.5),
                     .joint = {.type = JointType::Revolute, .parent = carrier, .axis = {0.0, 0.0, 1.0}}});
    auto model = builder.build();
    REQUIRE(model.has_value());
    return std::move(*model);
}

enum class Brake : std::uint8_t
{
    Reference, ///< The friction-loss row of the reference's law, solved with the mass matrix.
    Exact      ///< exactVelocityRow bounded by the capacity, solved with the step's own matrix.
};

// One brake row per dof, written as a consumer writes it: forward() solves the rows at the start of a step, step()
// completes it with the drive and the rows' generalized force.
struct BrakeStepper
{
    BrakeStepper(const Model& m, Brake brake_)
        : model(m), brake(brake_), n(m.dofCount()), ws(m), solver(m, m.dofCount()), jacobian(n * n, 0.0), kinds(n, SoftRowKind::FrictionLoss),
          reference(n), regularisation(n), force(n), free(n), generalized(n), tau(n), acceleration(n), weight(n)
    {
        Workspace rest(model);
        const State s = State::atRest(model);
        REQUIRE(prepareLinearlyImplicitStep(model, s, rest));
        REQUIRE(factorMassMatrix(model, rest));
        for (std::size_t d = 0; d < n; ++d)
        {
            jacobian[d * n + d] = 1.0;
            const auto w = dofInverseWeight(model, rest, d);
            REQUIRE(w.has_value());
            weight[d] = *w;
        }
    }

    // Solves the rows at the start of a step from `s` under `drive` with capacities `capacity`; returns the rows with a force.
    std::size_t forward(const State& s, std::span<const double> drive, std::span<const double> capacity, double h)
    {
        REQUIRE(prepareLinearlyImplicitStep(model, s, ws));
        for (std::size_t i = 0; i < n; ++i)
        {
            free[i] = ws.elementForce[i] - ws.bias[i] + drive[i];
        }
        std::size_t rows = 0;
        if (brake == Brake::Reference)
        {
            REQUIRE(factorMassMatrix(model, ws));
            REQUIRE(solveMassMatrix(model, ws, free));
            for (std::size_t d = 0; d < n; ++d)
            {
                const SoftConstraintRow row = softConstraintFrictionLossRow({}, s.qd()[d], weight[d], h);
                reference[d] = row.referenceAcceleration;
                regularisation[d] = row.regularisation;
            }
            const auto solved = solver.solve(model, ws, jacobian, reference, regularisation, free, kinds, force, generalized,
                                             {}, DelassusForm::MassSolvePerRow, capacity);
            REQUIRE(solved.has_value());
            rows = solved->activeRows;
            // a = a0 + M^-1 J^T f.
            std::copy(generalized.begin(), generalized.end(), acceleration.begin());
            REQUIRE(solveMassMatrix(model, ws, acceleration));
            for (std::size_t i = 0; i < n; ++i)
            {
                acceleration[i] += free[i];
            }
        }
        else
        {
            REQUIRE(factorPreparedStepMatrix(model, s, f64::Time::fromSi(h), ws));
            REQUIRE(solvePreparedStepMatrix(model, ws, free));
            for (std::size_t d = 0; d < n; ++d)
            {
                const SoftConstraintRow row = exactVelocityRow(s.qd()[d], weight[d], h);
                reference[d] = row.referenceAcceleration;
                regularisation[d] = row.regularisation;
            }
            const auto solved = solver.solveBounded(model, ws, jacobian, reference, regularisation, free, capacity, force,
                                                    generalized);
            REQUIRE(solved.has_value());
            rows = solved->activeRows;
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            tau[i] = drive[i] + generalized[i];
        }
        return rows;
    }

    void step(State& s, std::span<const double> drive, std::span<const double> capacity, double h)
    {
        forward(s, drive, capacity, h);
        REQUIRE(completeLinearlyImplicitStep(model, Forcing{.tau = tau}, s.tick, f64::Time::fromSi(h), s, ws));
    }

    const Model& model;
    Brake brake;
    std::size_t n;
    Workspace ws;
    SoftConstraintSolver solver;
    std::vector<double> jacobian;
    std::vector<SoftRowKind> kinds;
    std::vector<double> reference, regularisation, force, free, generalized, tau, acceleration, weight;
};

} // namespace

TEST_CASE("brake rows: the friction-loss row follows the reference's formulas, and a dof's inverse weight")
{
    // dof_solref (0.02, 1) and dof_solimp (0.9, 0.95, 0.001, 0.5, 2) by default; at zero position d = solimp[0].
    const double h = 1e-4;
    const SoftConstraintRow row = softConstraintFrictionLossRow({}, 2.0, 2.5, h);
    CHECK(row.impedance == doctest::Approx(0.9).epsilon(1e-15));
    CHECK(row.stiffness == 0.0);
    CHECK(row.damping == doctest::Approx(2.0 / (0.95 * 0.02)).epsilon(1e-15));
    CHECK(row.regularisation == doctest::Approx(0.1 * 2.5 / 0.9).epsilon(1e-15));
    CHECK(row.referenceAcceleration == doctest::Approx(-2.0 * 2.0 / (0.95 * 0.02)).epsilon(1e-15));
    // The time constant is raised to 2 h.
    const SoftConstraintRow coarse = softConstraintFrictionLossRow({}, 1.0, 2.5, 0.02);
    CHECK(coarse.damping == doctest::Approx(2.0 / (0.95 * 0.04)).epsilon(1e-15));

    // (M^-1)_ii: the rotor 1 / 0.4; the pair's M = [[2.5, 0.5], [0.5, 0.5]] has inverse [[0.5, -0.5], [-0.5, 2.5]].
    const Model single = rotor(0.4);
    Workspace ws(single);
    REQUIRE(prepareLinearlyImplicitStep(single, State::atRest(single), ws));
    REQUIRE(factorMassMatrix(single, ws));
    CHECK(dofInverseWeight(single, ws, 0).value() == doctest::Approx(2.5).epsilon(1e-14));
    CHECK_FALSE(dofInverseWeight(single, ws, 1).has_value());
    const Model pair = carrierAndRotor();
    Workspace pairWs(pair);
    REQUIRE(prepareLinearlyImplicitStep(pair, State::atRest(pair), pairWs));
    REQUIRE(factorMassMatrix(pair, pairWs));
    CHECK(dofInverseWeight(pair, pairWs, 0).value() == doctest::Approx(0.5).epsilon(1e-14));
    CHECK(dofInverseWeight(pair, pairWs, 1).value() == doctest::Approx(2.5).epsilon(1e-14));
    CHECK_FALSE(dofInverseWeight(pair, ws, 0).has_value());
}

TEST_CASE("brake rows: the reference's rotor controls through reversal, near rest, hold, capacity changes and release")
{
    // BikeLab test/drive-controls.test.mjs, "native rotor brakes ...", on the core: a hinge rotor, inertia 0.4, 1e-4 s.
    const Model model = rotor(0.4);
    BrakeStepper stepper(model, Brake::Reference);
    const double h = 1e-4;
    State s = State::atRest(model);
    const auto observe = [&](double speed, double drive, double capacity) {
        s.qd()[0] = speed;
        const std::array<double, 1> d{drive};
        const std::array<double, 1> c{capacity};
        const std::size_t rows = stepper.forward(s, d, c, h);
        const double torque = stepper.force[0];
        CHECK(stepper.acceleration[0] == doctest::Approx((drive + torque) / 0.4).epsilon(1e-12));
        CHECK(std::abs(torque) <= capacity + 1e-12);
        CHECK(rows == (capacity > 0.0 ? 1U : 0U));
        return torque;
    };
    for (const double sign : {-1.0, 1.0})
    {
        CHECK(std::abs(observe(2.0 * sign, 0.0, 0.8) + 0.8 * sign) <= 1e-8);
        CHECK(observe(1e-8 * sign, 0.0, 0.8) * sign <= 0.0);
        const double hold = observe(0.0, 0.4 * sign, 0.8);
        CHECK(hold * sign < 0.0);
        CHECK(std::abs(hold) <= 0.4);
        CHECK(std::abs(observe(0.0, 2.0 * sign, 0.8) + 0.8 * sign) <= 1e-8);
        CHECK(std::abs(observe(sign, 0.4 * sign, 0.0)) <= 1e-8);
    }
    // Regularised braking approaches rest; released, the rotor reverses under drive.
    observe(1.0, 0.0, 0.8);
    const std::array<double, 1> none{0.0};
    const std::array<double, 1> held{0.8};
    for (int i = 0; i < 6000; ++i)
    {
        stepper.step(s, none, held, h);
    }
    CHECK(std::abs(s.qd()[0]) < 0.01);
    const std::array<double, 1> reverse{-0.4};
    const std::array<double, 1> released{0.0};
    const double before = s.qd()[0];
    for (int i = 0; i < 200; ++i)
    {
        stepper.step(s, reverse, released, h);
    }
    CHECK(s.qd()[0] < before - 0.019);
}

TEST_CASE("brake rows: a free carrier and rotor take equal and opposite drive, with relative-coordinate work and refinement")
{
    // BikeLab test/drive-controls.test.mjs, "free carrier and rotor ...", on the core (implicitfast, here the linearly
    // implicit Euler step: without damping both are the semi-implicit step).
    std::array<double, 2> errors{};
    for (std::size_t pass = 0; pass < 2; ++pass)
    {
        const double h = pass == 0 ? 1e-4 : 5e-5;
        const Model model = carrierAndRotor();
        BrakeStepper stepper(model, Brake::Reference);
        State s = State::atRest(model);
        const std::array<double, 2> drive{0.0, 3.0};
        const std::array<double, 2> free{0.0, 0.0};
        stepper.forward(s, drive, free, h);
        CHECK(stepper.acceleration[0] == doctest::Approx(-1.5).epsilon(1e-12));
        CHECK(stepper.acceleration[1] == doctest::Approx(7.5).epsilon(1e-12));
        // The joints' axial moments: none at the carrier's, the drive at the rotor's.
        Workspace reactions(model);
        REQUIRE(computeReactions(model, s, drive, reactions));
        CHECK(std::abs(reactions.jointReactions[0].moment.z) <= 1e-12);
        CHECK(reactions.jointReactions[1].moment.z == doctest::Approx(3.0).epsilon(1e-12));
        const auto steps = static_cast<int>(std::lround(0.01 / h));
        for (int i = 0; i < steps; ++i)
        {
            stepper.step(s, drive, free, h);
        }
        // Angular momentum about the axis stays zero; the drive's work is on the relative coordinate.
        CHECK(std::abs(2.0 * s.qd()[0] + 0.5 * (s.qd()[0] + s.qd()[1])) <= 1e-12);
        const double work = 3.0 * s.q()[1];
        CHECK(work > 0.0);
        const auto energy = computeEnergy(model, s, reactions);
        REQUIRE(energy.has_value());
        errors[pass] = std::abs(energy->kinetic.si() - work);
        CHECK(std::abs(work - 3.0 * (s.q()[0] + s.q()[1])) > 1e-4);
        // A brake on the rotor's joint: its moment there is the row's generalized force; the carrier's joint carries none.
        const std::array<double, 2> braked{0.0, 4.0};
        stepper.forward(s, free, braked, h);
        REQUIRE(computeReactions(model, s, stepper.tau, reactions));
        CHECK(reactions.jointReactions[1].moment.z == doctest::Approx(stepper.generalized[1]).epsilon(1e-12));
        CHECK(std::abs(stepper.generalized[1]) > 0.0);
        CHECK(std::abs(reactions.jointReactions[0].moment.z) <= 1e-12);
    }
    MESSAGE("work against kinetic energy: " << errors[0] << " at 1e-4 s, " << errors[1] << " at 5e-5 s");
    CHECK(errors[1] < 0.51 * errors[0]);
    CHECK(errors[0] < 2e-5);
}

TEST_CASE("brake rows: the exact brake holds zero spin below capacity and slips at capacity above; the reference creeps")
{
    const Model model = rotor(0.4);
    const double h = 1e-4;
    const std::array<double, 1> drive{0.4};
    const std::array<double, 1> capacity{0.8};
    // From 0.5 rad/s with the drive aiding: sliding at the capacity, (0.4 - 0.8) / 0.4 = -1 rad/s^2, until rest at 0.5 s.
    BrakeStepper exact(model, Brake::Exact);
    State s = State::atRest(model);
    s.qd()[0] = 0.5;
    double worstHold = 0.0;
    double slide = 0.0;
    for (int i = 0; i < 10000; ++i)
    {
        exact.step(s, drive, capacity, h);
        if (i < 4900)
        {
            slide = std::max(slide, std::abs(exact.force[0] + 0.8));
        }
        if (i > 5100)
        {
            worstHold = std::max(worstHold, std::abs(s.qd()[0]));
            CHECK(exact.force[0] == doctest::Approx(-0.4).epsilon(1e-9));
        }
    }
    CHECK(slide == 0.0);
    CHECK(worstHold <= 1e-12);
    // Over capacity it slips at C, either way.
    for (const double sign : {-1.0, 1.0})
    {
        State over = State::atRest(model);
        const std::array<double, 1> strong{2.0 * sign};
        exact.forward(over, strong, capacity, h);
        CHECK(exact.force[0] == -0.8 * sign);
        over.qd()[0] = -2.0 * sign;
        const std::array<double, 1> none{0.0};
        exact.forward(over, none, capacity, h);
        CHECK(exact.force[0] == 0.8 * sign);
        const std::array<double, 1> zero{0.0};
        CHECK(exact.forward(over, strong, zero, h) == 0U);
        CHECK(exact.force[0] == 0.0);
    }
    // The reference's regularised hold creeps at a steady rate: f = -drive inside the box needs B v = drive R.
    BrakeStepper soft(model, Brake::Reference);
    State r = State::atRest(model);
    for (int i = 0; i < 20000; ++i)
    {
        soft.step(r, drive, capacity, h);
    }
    const SoftConstraintRow row = softConstraintFrictionLossRow({}, 0.0, 2.5, h);
    const double creep = 0.4 * row.regularisation / row.damping;
    MESSAGE("reference creep under a held brake: " << r.qd()[0] << " rad/s (drive R / B = " << creep << ")");
    CHECK(r.qd()[0] == doctest::Approx(creep).epsilon(1e-6));
    CHECK(soft.force[0] == doctest::Approx(-0.4).epsilon(1e-9));
}
