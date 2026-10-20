// Constraints (docs/MULTIBODY.md section 6 and the M1c gate): slider-crank closure residuals bounded and not
// drifting, a pendulum against a hard stop that is never overshot, reactions against independent statics, the energy
// ledger with constraint work, and the target-driven levers (achieved against requested, cost against target).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <doctest.h>
#include <ostream> // doctest prints strings (MSVC needs it explicitly)
#include <expected>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

#include <sub0math/units.hpp>

#include "sub0physics/multibody.hpp"
#include "sub0physics/verification/energy_ledger.hpp"
#include "test_multibody_models.hpp"

using namespace sub0physics;
using namespace sub0physics::multibody;
using namespace testmodels;

namespace
{

constexpr double G = 9.81;

sub0math::f64::Time seconds(double s)
{
    return sub0math::f64::Time::fromSi(s);
}

Model sliderCrank(const ConstraintSettings& settings = {}, const SliderCrankSpec& spec = {})
{
    ModelBuilder builder = sliderCrankBuilder(spec);
    builder.setConstraintSettings(settings);
    return *builder.build();
}

// The closure's residual at a state, measured independently of the solver: rod tip minus slider pin.
double closureGap(const Model& model, const State& state, Workspace& ws)
{
    const auto tip = computeSiteKinematics(model, state, *model.findSite("rod tip"), ws).value();
    const auto pin = computeSiteKinematics(model, state, *model.findSite("slider pin"), ws).value();
    const double dx = tip.pose.translation.x.si() - pin.pose.translation.x.si();
    const double dy = tip.pose.translation.y.si() - pin.pose.translation.y.si();
    const double dz = tip.pose.translation.z.si() - pin.pose.translation.z.si();
    return std::max({std::abs(dx), std::abs(dy), std::abs(dz)});
}

// A 1 kg, 1 m link hanging along -y from a revolute joint about z at the origin, with a lower hard stop.
Model pendulumWithStop(double lower, double restitution, const ConstraintSettings& settings = {})
{
    ModelBuilder builder;
    builder.setGravity(gravityAlong(0.0, -G, 0.0));
    const BodyId arm = builder.addBody(planarLink("arm", WorldBody, {}, 1.0, 1.0));
    builder.addConstraint(CoordinateLimitDescription{.name = "stop", .body = arm, .lower = lower, .restitution = restitution});
    builder.setConstraintSettings(settings);
    return *builder.build();
}

double pendulumPotential(double q)
{
    return -G * 0.5 * std::cos(q); // m g y of the centre of mass, m = 1, L = 1
}

} // namespace

TEST_CASE("Slider-crank: closure residual bounded by the requested tolerance and not drifting, every stepper")
{
    const Model model = sliderCrank();
    const SliderCrankSpec spec;
    REQUIRE(model.closureRowCount() == 3U);
    const double requested = model.constraintSettings().closureTolerance;

    const auto run = [&](auto&& step, double h, int steps) {
        Workspace ws(model);
        State s = sliderCrankState(spec, 0.7, 4.0);
        REQUIRE(closureGap(model, s, ws) < 1e-15);
        double worst = 0.0;
        double worstVelocity = 0.0;
        double firstHalf = 0.0;
        double secondHalf = 0.0;
        for (int k = 0; k < steps; ++k)
        {
            REQUIRE(step(s, ws, h));
            const double gap = closureGap(model, s, ws);
            worst = std::max(worst, gap);
            worstVelocity = std::max(worstVelocity, ws.constraintReport.closureVelocityResidual);
            (k < steps / 2 ? firstHalf : secondHalf) = std::max(k < steps / 2 ? firstHalf : secondHalf, gap);
        }
        CHECK(worst <= requested);
        CHECK(worstVelocity < 1e-12);
        CHECK(secondHalf <= std::max(2.0 * firstHalf, requested)); // no drift
        MESSAGE("closure residual: worst " << worst << " m (requested " << requested << "), velocity " << worstVelocity);
        return s;
    };
    SUBCASE("linearly implicit Euler (default)")
    {
        run([&](State& s, Workspace& ws, double h) { return stepLinearlyImplicitEuler(model, {}, s.tick, seconds(h), s, ws); },
            1e-3, 10000);
    }
    SUBCASE("semi-implicit Euler")
    {
        run([&](State& s, Workspace& ws, double h) { return stepSemiImplicitEuler(model, {}, s.tick, seconds(h), s, ws); },
            1e-3, 10000);
    }
    SUBCASE("RK4 reference")
    {
        Workspace rkWs(model);
        Rk4Stepper rk4(model, rkWs);
        run([&](State& s, Workspace& ws, double h) { return rk4.step(model, {}, s.tick, seconds(h), s, ws); }, 1e-3, 3000);
    }
}

TEST_CASE("Slider-crank: without projection the velocity-level closure drifts (the projection is what holds it)")
{
    const Model model = sliderCrank({.maxProjectionIterations = 0});
    Workspace ws(model);
    State s = sliderCrankState({}, 0.7, 4.0);
    for (int k = 0; k < 5000; ++k)
    {
        REQUIRE(stepLinearlyImplicitEuler(model, {}, s.tick, seconds(1e-3), s, ws));
    }
    const double drifted = closureGap(model, s, ws);
    MESSAGE("closure residual after 5 s without projection: " << drifted << " m");
    CHECK(drifted > 1e3 * model.constraintSettings().closureTolerance);
    CHECK(ws.constraintReport.closureResidual == doctest::Approx(drifted).epsilon(1e-6));
}

TEST_CASE("Slider-crank: the energy ledger closes, with the closures' work and the projection's energy reported")
{
    const Model model = sliderCrank();
    const std::array<double, 3> tau{0.05, 0.0, -0.1};
    SUBCASE("RK4: the residual is the projection's (and RK's) own")
    {
        Workspace ws(model);
        Rk4Stepper rk4(model, ws);
        State s = sliderCrankState({}, 0.7, 4.0);
        EnergyLedger ledger(computeEnergy(model, s, ws).value(), sub0math::f64::Energy::fromSi(1.0));
        for (int k = 0; k < 2000; ++k)
        {
            REQUIRE(rk4.step(model, {tau}, s.tick, seconds(1e-3), s, ws));
            const EnergyAccount account = computeEnergy(model, s, ws).value();
            ledger.includeScale(account.kinetic);
            ledger.record(account);
        }
        const double projection = s.projectionEnergy.si();
        MESSAGE("RK4 ledger relative residual " << ledger.relativeResidual() << ", projection energy " << projection
                                                << " J, scale " << ledger.scale().si() << " J");
        CHECK(ledger.relativeResidual() < 1e-6);
        CHECK(std::abs(projection) < 1e-6 * ledger.scale().si());
        // The run's balance residual is accounted for by the projection to within RK's own error.
        const double residual = ledger.residual(computeEnergy(model, s, ws).value()).si();
        CHECK(std::abs(residual - projection) < 1e-6 * ledger.scale().si());
    }
    SUBCASE("implicit Euler: the balance, the closures' work and the projection's energy are all first order")
    {
        // The velocity-level solve uses J at the start of the step and the projection removes what that leaves
        // normal to the manifold at the end, so the closures' impulse work and the projection's energy are O(h)
        // over a run, like the step's own error: halving h halves each channel.
        std::array<std::array<double, 3>, 3> channels{}; // balance residual, closure work, projection energy
        const std::array<double, 3> steps{2e-3, 1e-3, 5e-4};
        for (std::size_t m = 0; m < steps.size(); ++m)
        {
            Workspace ws(model);
            State s = sliderCrankState({}, 0.7, 4.0);
            const EnergyLedger ledger(computeEnergy(model, s, ws).value(), sub0math::f64::Energy::fromSi(1.0));
            const int count = static_cast<int>(2.0 / steps[m]);
            for (int k = 0; k < count; ++k)
            {
                REQUIRE(stepLinearlyImplicitEuler(model, {}, s.tick, seconds(steps[m]), s, ws));
            }
            channels[m] = {ledger.residual(computeEnergy(model, s, ws).value()).si(), s.closureWork.si(),
                           s.projectionEnergy.si()};
            MESSAGE("implicit Euler h " << steps[m] << ": balance residual " << channels[m][0] << " J, closure work "
                                        << channels[m][1] << " J, projection " << channels[m][2] << " J");
        }
        for (std::size_t c = 0; c < 3; ++c)
        {
            CHECK(channels[0][c] / channels[1][c] == doctest::Approx(2.0).epsilon(0.15));
            CHECK(channels[1][c] / channels[2][c] == doctest::Approx(2.0).epsilon(0.15));
        }
    }
}

TEST_CASE("Projection: achieved residual meets the requested target, and cost falls as the target loosens")
{
    std::array<double, 4> targets{1e-5, 1e-8, 1e-11, 1e-14};
    std::array<std::uint32_t, 4> iterations{};
    for (std::size_t t = 0; t < targets.size(); ++t)
    {
        const Model model = sliderCrank({.closureTolerance = targets[t]});
        Workspace ws(model);
        State s = sliderCrankState({}, 0.7, 4.0);
        s.q()[2] += 0.02; // assemble from 2 cm off
        s.q()[1] -= 0.05;
        REQUIRE(projectConstraints(model, s, ws));
        const double achieved = closureGap(model, s, ws);
        iterations[t] = ws.constraintReport.projectionIterations;
        MESSAGE("requested " << targets[t] << ", achieved " << achieved << ", iterations " << iterations[t]);
        CHECK(achieved <= targets[t]);
        CHECK(ws.constraintReport.closureResidual <= targets[t]);
    }
    for (std::size_t t = 1; t < targets.size(); ++t)
    {
        CHECK(iterations[t] >= iterations[t - 1]);
    }
    CHECK(iterations.back() > iterations.front());

    // A target the guard cannot reach is a fault, not a silent miss.
    const Model strict = sliderCrank({.closureTolerance = 1e-14, .maxProjectionIterations = 1});
    Workspace ws(strict);
    State s = sliderCrankState({}, 0.7, 4.0);
    s.q()[2] += 0.02;
    const State before = s;
    const auto refused = projectConstraints(strict, s, ws);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().kind == FaultKind::OutsideDomain);
    CHECK(s == before);
}

TEST_CASE("Hard stop: a pendulum never overshoots its limit, and the impact is located and accounted for")
{
    const double lower = -0.5;
    const double q0 = 1.2;
    const double impactEnergy = pendulumPotential(q0) - pendulumPotential(lower);

    const auto run = [&](double restitution, auto&& step, double h, int steps) {
        const Model model = pendulumWithStop(lower, restitution);
        Workspace ws(model);
        State s = State::fromCoordinates({q0}, {0.0});
        EnergyLedger ledger(computeEnergy(model, s, ws).value(), sub0math::f64::Energy::fromSi(1.0));
        double lowest = q0;
        double impactRate = 0.0;
        double reboundRate = 0.0;
        double dissipatedAtFirst = 0.0;
        bool impacted = false;
        for (int k = 0; k < steps; ++k)
        {
            const double before = s.qd()[0];
            REQUIRE(step(model, s, ws, h));
            lowest = std::min(lowest, s.q()[0]);
            CHECK(s.q()[0] >= lower); // exactly: an accepted state is inside its limits
            if (!impacted && ws.constraintReport.substeps > 1U)
            {
                impacted = true;
                impactRate = before; // within one step's gravity of the rate at contact
                reboundRate = s.qd()[0];
                dissipatedAtFirst = s.dissipated.si();
                CHECK(ws.constraintReport.engagementGap >= 0.0);
                CHECK(ws.constraintReport.engagementGap <= model.constraintSettings().eventTolerance);
                MESSAGE("e " << restitution << ": engagement gap " << ws.constraintReport.engagementGap
                             << " (requested <= " << model.constraintSettings().eventTolerance << "), "
                             << ws.constraintReport.eventIterations << " event trials");
            }
            const EnergyAccount account = computeEnergy(model, s, ws).value();
            ledger.includeScale(account.kinetic);
            ledger.record(account);
        }
        REQUIRE(impacted);
        CHECK(lowest >= lower);
        const double exactImpactRate = -std::sqrt(2.0 * impactEnergy / (1.0 / 3.0));
        MESSAGE("e " << restitution << ": impact rate " << impactRate << " (exact " << exactImpactRate << "), rebound "
                     << reboundRate << ", dissipated at impact " << dissipatedAtFirst << " J (exact "
                     << (1.0 - restitution * restitution) * impactEnergy << "), ledger " << ledger.relativeResidual());
        CHECK(impactRate == doctest::Approx(exactImpactRate).epsilon(0.02));
        CHECK(dissipatedAtFirst == doctest::Approx((1.0 - restitution * restitution) * impactEnergy).epsilon(0.03));
        CHECK(ledger.relativeResidual() < 0.02);
        return reboundRate / impactRate;
    };
    const auto implicitStep = [](const Model& model, State& s, Workspace& ws, double h) {
        return stepLinearlyImplicitEuler(model, {}, s.tick, seconds(h), s, ws);
    };
    const auto semiImplicitStep = [](const Model& model, State& s, Workspace& ws, double h) {
        return stepSemiImplicitEuler(model, {}, s.tick, seconds(h), s, ws);
    };
    SUBCASE("plastic stop (restitution 0), implicit Euler")
    {
        CHECK(std::abs(run(0.0, implicitStep, 1e-3, 3000)) < 0.02);
    }
    SUBCASE("restitution 0.5, implicit Euler")
    {
        CHECK(run(0.5, implicitStep, 1e-3, 3000) == doctest::Approx(-0.5).epsilon(0.03));
    }
    SUBCASE("restitution 0.5, semi-implicit Euler")
    {
        CHECK(run(0.5, semiImplicitStep, 5e-4, 6000) == doctest::Approx(-0.5).epsilon(0.03));
    }
}

TEST_CASE("Event location: achieved engagement gap meets the requested tolerance; cost grows smoothly as it tightens")
{
    std::array<double, 4> targets{1e-3, 1e-6, 1e-9, 1e-12};
    std::array<std::uint32_t, 4> trials{};
    for (std::size_t t = 0; t < targets.size(); ++t)
    {
        const Model model = pendulumWithStop(-0.5, 0.0, {.eventTolerance = targets[t]});
        Workspace ws(model);
        State s = State::fromCoordinates({1.2}, {0.0});
        for (int k = 0; k < 2000 && trials[t] == 0U; ++k)
        {
            REQUIRE(stepLinearlyImplicitEuler(model, {}, s.tick, seconds(1e-3), s, ws));
            if (ws.constraintReport.substeps > 1U)
            {
                trials[t] = ws.constraintReport.eventIterations;
                MESSAGE("requested " << targets[t] << ", achieved gap " << ws.constraintReport.engagementGap
                                     << ", trials " << trials[t]);
                CHECK(ws.constraintReport.engagementGap <= targets[t]);
                CHECK(ws.constraintReport.engagementGap >= 0.0);
            }
        }
        REQUIRE(trials[t] > 0U);
    }
    CHECK(trials.back() >= trials.front());
}

TEST_CASE("Split steps keep integer ticks and determinism")
{
    const Model model = pendulumWithStop(-0.5, 0.3);
    const auto run = [&] {
        Workspace ws(model);
        State s = State::fromCoordinates({1.2}, {0.0});
        for (int k = 0; k < 1500; ++k)
        {
            REQUIRE(stepLinearlyImplicitEuler(model, {}, s.tick, seconds(1e-3), s, ws));
        }
        return s;
    };
    const State a = run();
    const State b = run();
    CHECK(a.tick == 1500U);
    CHECK(a == b);
}

TEST_CASE("Rope (maximum distance limit) is never overshot")
{
    ModelBuilder builder;
    builder.setGravity(gravityAlong(0.0, -G, 0.0));
    const auto point = [](std::string_view name, BodyId parent, sub0math::Vec3d axis) {
        return BodyDescription{.name = name,
                               .mass = sub0math::kilograms(name == "x" ? 1e-3 : 1.0),
                               .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                               .inertia = diagonalInertia(1e-3, 1e-3, 1e-3),
                               .joint = {.type = JointType::Prismatic, .parent = parent, .axis = axis}};
    };
    const BodyId x = builder.addBody(point("x", WorldBody, {1.0, 0.0, 0.0}));
    const BodyId y = builder.addBody(point("y", x, {0.0, 1.0, 0.0}));
    const SiteId anchor = builder.addSite({.name = "anchor", .body = WorldBody, .frame = {}});
    const SiteId bob = builder.addSite({.name = "bob", .body = y, .frame = {}});
    builder.addConstraint(DistanceLimitDescription{.name = "rope", .first = anchor, .second = bob, .maximum = metres(0.5)});
    const Model model = *builder.build();
    Workspace ws(model);
    State s = State::fromCoordinates({0.3, 0.0}, {0.0, 0.0});
    double longest = 0.0;
    int splits = 0;
    for (int k = 0; k < 3000; ++k)
    {
        REQUIRE(stepLinearlyImplicitEuler(model, {}, s.tick, seconds(1e-3), s, ws));
        longest = std::max(longest, std::hypot(s.q()[0], s.q()[1]));
        splits += ws.constraintReport.substeps > 1U ? 1 : 0;
    }
    MESSAGE("longest rope length " << longest << " m (limit 0.5), split steps " << splits);
    CHECK(splits >= 1);
    CHECK(longest <= 0.5 + model.constraintSettings().closureTolerance);
    CHECK(longest > 0.5 - 1e-6);
}

TEST_CASE("Reactions: a loaded slider-crank against its hard stop matches independent statics")
{
    const SliderCrankSpec spec;
    const double theta = 0.6;
    State s = sliderCrankState(spec, theta, 0.0);
    const double x = s.q()[2];
    ModelBuilder builder = sliderCrankBuilder(spec);
    builder.addConstraint(CoordinateLimitDescription{.name = "end stop", .body = BodyId{2}, .upper = x});
    const Model model = *builder.build();
    Workspace ws(model);
    const double torque = 2.0; // on the crank
    const double push = 60.0;  // on the slider, toward the stop
    const std::array<double, 3> tau{torque, 0.0, push};

    // Independent planar statics (free bodies): P is the force of the slider on the rod at the pin.
    const double a = spec.crank;
    const double b = spec.rod;
    const double phi = -std::asin(a * std::sin(theta) / b);
    const double ax = a * std::cos(theta);
    const double ay = a * std::sin(theta);
    const double bx = ax + b * std::cos(phi);
    const double by = ay + b * std::sin(phi);
    const auto cross2 = [](double rx, double ry, double fx, double fy) { return rx * fy - ry * fx; };
    const double w1 = -spec.crankMass * G;
    const double w2 = -spec.rodMass * G;
    const double w3 = -spec.sliderMass * G;
    // Rod about A: cross(B - A, P) + cross(C2 - A, W2) = 0. Crank about O: T + cross(A, P + W2) + cross(C1, W1) = 0.
    const double r11 = -(by - ay);
    const double r12 = bx - ax;
    const double s1 = -cross2(0.5 * (bx - ax), 0.5 * (by - ay), 0.0, w2);
    const double r21 = -ay;
    const double r22 = ax;
    const double s2 = -torque - cross2(ax, ay, 0.0, w2) - cross2(0.5 * ax, 0.5 * ay, 0.0, w1);
    const double det = r11 * r22 - r12 * r21;
    const double px = (s1 * r22 - r12 * s2) / det;
    const double py = (r11 * s2 - s1 * r21) / det;
    const double stopForce = px - push; // on the slider along x; the stop can only push back (negative)
    const double normal = py - w3;
    REQUIRE(stopForce < 0.0);

    REQUIRE(computeReactions(model, s, tau, ws));
    const double scale = std::max({std::abs(px), std::abs(py), push});
    const double tolerance = 1e-9 * scale;
    for (std::size_t i = 0; i < 3; ++i)
    {
        CHECK(std::abs(ws.qdd[i]) < 1e-9); // statically determinate: held still
    }
    const UnilateralLimit& stop = model.limits()[0];
    CHECK(ws.constraintForce[stop.row] == doctest::Approx(-stopForce).epsilon(1e-9));
    // The closure's force on the rod tip is P.
    const ConstraintSiteWrench& onRod = ws.constraintSiteWrenches[0];
    CHECK(std::abs(onRod.force.x - px) < tolerance);
    CHECK(std::abs(onRod.force.y - py) < tolerance);
    // Joint reactions: crank from the world, rod from the crank, slider from the world.
    const JointReaction& crank = ws.jointReactions[0];
    CHECK(std::abs(crank.force.x - (-px)) < tolerance);
    CHECK(std::abs(crank.force.y - (-py - w2 - w1)) < tolerance);
    CHECK(std::abs(crank.moment.z - torque) < tolerance);
    const JointReaction& rod = ws.jointReactions[1];
    CHECK(std::abs(rod.force.x - (-px)) < tolerance);
    CHECK(std::abs(rod.force.y - (-py - w2)) < tolerance);
    CHECK(std::abs(rod.moment.z) < tolerance);
    const JointReaction& slider = ws.jointReactions[2];
    CHECK(std::abs(slider.force.x - (push + stopForce)) < tolerance);
    CHECK(std::abs(slider.force.y - normal) < tolerance);
    CHECK(std::abs(slider.moment.z) < tolerance);
    MESSAGE("statics: P = (" << px << ", " << py << ") N, stop " << stopForce << " N; library closure ("
                             << onRod.force.x << ", " << onRod.force.y << "), stop " << -ws.constraintForce[stop.row]);
}

TEST_CASE("Reactions: each joint's axial component equals its generalized force, in motion too")
{
    const Model model = sliderCrank();
    Workspace ws(model);
    const State s = sliderCrankState({}, 1.1, -3.0);
    const std::array<double, 3> tau{0.7, -0.2, 4.0};
    REQUIRE(computeReactions(model, s, tau, ws));
    for (std::size_t i = 0; i < 3; ++i)
    {
        const JointReaction& r = ws.jointReactions[i];
        // World axis: z for the revolute joints, x for the slider.
        const double axial = model.jointTypes()[i] == JointType::Revolute ? r.moment.z : r.force.x;
        CHECK(axial == doctest::Approx(tau[i]).epsilon(1e-9));
    }
}

TEST_CASE("Frame closure locks a planar three-link chain to the world")
{
    ModelBuilder builder;
    builder.setGravity(gravityAlong(0.0, -G, 0.0));
    BodyId parent = WorldBody;
    for (int k = 0; k < 3; ++k)
    {
        parent = builder.addBody(linkAlongX(k == 0 ? "a" : (k == 1 ? "b" : "c"), parent,
                                            k == 0 ? Pose{} : offset(0.4, 0.0, 0.0), 1.0, 0.4));
    }
    const SiteId tip = builder.addSite({.name = "tip", .body = parent, .frame = offset(0.4, 0.0, 0.0)});
    // Bent, so the six rows have rank three (a straight chain welded at both ends is a singular configuration).
    const std::array<double, 3> bent{0.3, -0.6, 0.3};
    const double reach = 0.4 * (2.0 * std::cos(0.3) + 1.0);
    const SiteId mount = builder.addSite({.name = "mount", .body = WorldBody, .frame = offset(reach, 0.0, 0.0)});
    builder.addConstraint(LoopClosureDescription{.name = "weld", .first = tip, .second = mount, .kind = ClosureKind::Frame});
    const Model model = *builder.build();
    REQUIRE(model.closureRowCount() == 6U);
    Workspace ws(model);
    State s = State::fromCoordinates({bent[0], bent[1], bent[2]}, {0.0, 0.0, 0.0});
    for (int k = 0; k < 1000; ++k)
    {
        REQUIRE(stepLinearlyImplicitEuler(model, {}, s.tick, seconds(1e-3), s, ws));
    }
    for (std::size_t i = 0; i < 3; ++i)
    {
        CHECK(std::abs(s.q()[i] - bent[i]) < 1e-10);
    }
    REQUIRE(computeReactions(model, s, {}, ws));
    // The split of the weight between the pin and the weld is statically indeterminate; their sum is not.
    const double vertical = ws.jointReactions[0].force.y + ws.constraintSiteWrenches[0].force.y;
    CHECK(vertical == doctest::Approx(3.0 * G).epsilon(1e-9));
    CHECK(ws.constraintReport.closureResidual <= model.constraintSettings().closureTolerance);
}

TEST_CASE("Constraint descriptions are validated by the builder")
{
    const auto fails = [](auto&& configure) {
        ModelBuilder builder = sliderCrankBuilder();
        configure(builder);
        const auto built = builder.build();
        return !built.has_value() && built.error().kind == FaultKind::InvalidInput;
    };
    CHECK(fails([](ModelBuilder& b) { b.addConstraint(LoopClosureDescription{.name = "loop", .first = SiteId{0}, .second = SiteId{1}}); }));
    CHECK(fails([](ModelBuilder& b) { b.addConstraint(LoopClosureDescription{.name = "x", .first = SiteId{0}, .second = SiteId{7}}); }));
    CHECK(fails([](ModelBuilder& b) {
        const SiteId other = b.addSite({.name = "other", .body = BodyId{1}, .frame = {}});
        b.addConstraint(LoopClosureDescription{.name = "same body", .first = SiteId{0}, .second = other});
    }));
    CHECK(fails([](ModelBuilder& b) { b.addConstraint(CoordinateLimitDescription{.name = "none", .body = BodyId{0}}); }));
    CHECK(fails([](ModelBuilder& b) {
        b.addConstraint(CoordinateLimitDescription{.name = "crossed", .body = BodyId{0}, .lower = 1.0, .upper = 0.5});
    }));
    CHECK(fails([](ModelBuilder& b) {
        b.addConstraint(CoordinateLimitDescription{.name = "bouncy", .body = BodyId{0}, .lower = 0.0, .restitution = 1.5});
    }));
    CHECK(fails([](ModelBuilder& b) {
        b.addConstraint(DistanceLimitDescription{.name = "rope", .first = SiteId{0}, .second = SiteId{1}, .maximum = metres(-1.0)});
    }));
    CHECK(fails([](ModelBuilder& b) { b.setConstraintSettings({.closureTolerance = 0.0}); }));
    CHECK(fails([](ModelBuilder& b) { b.setConstraintSettings({.solverTolerance = 0.1}); }));
    CHECK_FALSE(fails([](ModelBuilder&) {}));
}

TEST_CASE("Steppers refuse what they cannot honour and stay transactional")
{
    const Model model = pendulumWithStop(-0.5, 0.0);
    Workspace ws(model);
    Rk4Stepper rk4(model, ws);
    State s = State::fromCoordinates({0.2}, {0.0});
    const State before = s;
    const auto refused = rk4.step(model, {}, s.tick, seconds(1e-3), s, ws);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().kind == FaultKind::InvalidInput);
    CHECK(s == before);

    State outside = State::fromCoordinates({-0.6}, {0.0});
    const State outsideBefore = outside;
    const auto violated = stepLinearlyImplicitEuler(model, {}, outside.tick, seconds(1e-3), outside, ws);
    REQUIRE_FALSE(violated);
    CHECK(violated.error().kind == FaultKind::InvalidInput);
    CHECK(outside == outsideBefore);
}
