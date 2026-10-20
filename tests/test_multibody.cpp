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
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <sub0math/linear/quat.hpp>
#include <sub0math/units.hpp>

#include "sub0physics/integrators.hpp"
#include "sub0physics/multibody.hpp"
#include "test_multibody_models.hpp"

using namespace sub0physics;
using namespace sub0physics::multibody;
using namespace testmodels;

namespace
{

// Complete elliptic integral of the first kind K(k) by the arithmetic-geometric mean: K = pi / (2 AGM(1, sqrt(1-k^2))).
// std::comp_ellint_1 is not available in libc++ (macOS), so the test computes it.
double completeEllipticK(double k)
{
    double a = 1.0;
    double b = std::sqrt(1.0 - k * k);
    for (int i = 0; i < 32 && std::abs(a - b) > 1e-16 * a; ++i)
    {
        const double next = 0.5 * (a + b);
        b = std::sqrt(a * b);
        a = next;
    }
    return 3.14159265358979323846 / (2.0 * a);
}


constexpr double G = 9.81;

double element(const std::vector<double>& m, std::size_t n, std::size_t r, std::size_t c)
{
    return m[r * n + c];
}

// A pendulum about z with its centre of mass `length` below the pivot (gravity along -y).
Model pendulum(double mass, double length, double inertiaAboutCentre)
{
    ModelBuilder builder;
    builder.setGravity(gravityAlong(0.0, -G, 0.0));
    builder.addBody({.name = "bob",
                     .mass = sub0math::kilograms(mass),
                     .centreOfMass = {metres(0.0), metres(-length), metres(0.0)},
                     .inertia = diagonalInertia(inertiaAboutCentre, inertiaAboutCentre, inertiaAboutCentre),
                     .joint = {.type = JointType::Revolute, .axis = {0.0, 0.0, 1.0}}});
    return *builder.build();
}

Model doublePendulum()
{
    ModelBuilder builder;
    builder.setGravity(gravityAlong(0.0, -G, 0.0));
    const BodyId upper = builder.addBody(planarLink("upper", WorldBody, {}, 1.0, 1.0));
    builder.addBody(planarLink("lower", upper, offset(0.0, -1.0, 0.0), 0.7, 0.8));
    return *builder.build();
}

// Cart on x carrying two poles about z: the hand-derived reference below.
struct CartPoles
{
    double cart = 2.0;
    std::array<double, 2> mass{0.5, 0.8};
    std::array<double, 2> length{0.6, 0.35}; // pivot to centre of mass
    std::array<double, 2> inertia{0.02, 0.01}; // about the centre of mass, zz
    std::array<double, 2> pivot{-0.2, 0.3};    // along the cart

    [[nodiscard]] Model model() const
    {
        ModelBuilder builder;
        builder.setGravity(gravityAlong(0.0, -G, 0.0));
        const BodyId base = builder.addBody({.name = "cart",
                                             .mass = sub0math::kilograms(cart),
                                             .centreOfMass = {metres(0.1), metres(0.05), metres(0.0)},
                                             .inertia = diagonalInertia(0.3, 0.3, 0.3),
                                             .joint = {.type = JointType::Prismatic, .axis = {1.0, 0.0, 0.0}}});
        for (std::size_t k = 0; k < 2; ++k)
        {
            builder.addBody({.name = k == 0 ? "pole0" : "pole1",
                             .mass = sub0math::kilograms(mass[k]),
                             .centreOfMass = {metres(0.0), metres(-length[k]), metres(0.0)},
                             .inertia = diagonalInertia(inertia[k], inertia[k], inertia[k]),
                             .joint = {.type = JointType::Revolute, .parent = base, .placement = offset(pivot[k], 0.0, 0.0),
                                       .axis = {0.0, 0.0, 1.0}}});
        }
        return *builder.build();
    }

    // Lagrange's equations for q = (x, theta0, theta1).
    [[nodiscard]] std::array<double, 9> massMatrix(const State& s) const
    {
        const double c0 = mass[0] * length[0] * std::cos(s.q()[1]);
        const double c1 = mass[1] * length[1] * std::cos(s.q()[2]);
        return {cart + mass[0] + mass[1], c0, c1,
                c0, inertia[0] + mass[0] * length[0] * length[0], 0.0,
                c1, 0.0, inertia[1] + mass[1] * length[1] * length[1]};
    }

    [[nodiscard]] std::array<double, 3> bias(const State& s) const
    {
        const double s0 = std::sin(s.q()[1]);
        const double s1 = std::sin(s.q()[2]);
        return {-mass[0] * length[0] * s0 * s.qd()[1] * s.qd()[1] - mass[1] * length[1] * s1 * s.qd()[2] * s.qd()[2],
                mass[0] * G * length[0] * s0, mass[1] * G * length[1] * s1};
    }
};

// A general 3D tree: mixed joints, skew axes and placements, full inertia tensors.
Model spatialTree(const std::array<sub0math::Quatd, 6>& frameTurn = {})
{
    using sub0math::Quatd;
    using sub0math::Vec3d;
    const std::array<std::uint32_t, 6> parents{WorldBody.value(), 0, 1, 1, 3, 0};
    const std::array<JointType, 6> types{JointType::Revolute, JointType::Prismatic, JointType::Revolute,
                                         JointType::Revolute, JointType::Prismatic, JointType::Revolute};
    const std::array<Vec3d, 6> axes{Vec3d{0.0, 0.0, 1.0}, Vec3d{1.0, 2.0, -0.5}, Vec3d{0.3, -1.0, 0.2},
                                    Vec3d{1.0, 0.0, 0.0}, Vec3d{0.0, 1.0, 1.0}, Vec3d{-0.4, 0.1, 1.0}};
    ModelBuilder builder;
    builder.setGravity(gravityAlong(0.5, -1.0, -9.6));
    for (std::size_t k = 0; k < 6; ++k)
    {
        const double f = static_cast<double>(k);
        Pose placement{Quatd::fromAxisAngle(sub0math::normalized(Vec3d{1.0, f, 0.5}), sub0math::f64::Angle::fromSi(0.3 * f - 0.4)),
                       {metres(0.1 * f), metres(-0.3), metres(0.05 * f)}};
        sub0math::Vec3d com{0.05 * f, -0.2, 0.1 - 0.03 * f};
        sub0math::Vec3d axis = sub0math::normalized(axes[k]);
        // Principal moments near these diagonal values obey the triangle inequality (a real body).
        RotationalInertia inertia{{sub0math::kilogramMetresSquared(0.05 + 0.004 * f), sub0math::kilogramMetresSquared(0.06),
                                   sub0math::kilogramMetresSquared(0.07 - 0.003 * f)},
                                  {sub0math::kilogramMetresSquared(0.003), sub0math::kilogramMetresSquared(-0.002),
                                   sub0math::kilogramMetresSquared(0.001 * f)}};
        // The same physical body, described in a body frame turned by frameTurn[k] (see the invariance test).
        const Quatd turn = frameTurn[k];
        const Quatd back = sub0math::conjugate(turn);
        const Quatd parentTurn = parents[k] == WorldBody.value() ? Quatd{} : frameTurn[parents[k]];
        placement = Pose{sub0math::conjugate(parentTurn), {}} * placement * Pose{turn, {}};
        axis = rotate(back, axis);
        com = rotate(back, com);
        const sub0math::Mat3d r = sub0math::toMat3(turn);
        sub0math::Mat3d i{};
        i(0, 0) = inertia.moments.x.si();
        i(1, 1) = inertia.moments.y.si();
        i(2, 2) = inertia.moments.z.si();
        i(0, 1) = i(1, 0) = inertia.products.x.si();
        i(0, 2) = i(2, 0) = inertia.products.y.si();
        i(1, 2) = i(2, 1) = inertia.products.z.si();
        const sub0math::Mat3d turned = sub0math::transpose(r) * i * r;
        using sub0math::kilogramMetresSquared;
        inertia = {{kilogramMetresSquared(turned(0, 0)), kilogramMetresSquared(turned(1, 1)), kilogramMetresSquared(turned(2, 2))},
                   {kilogramMetresSquared(turned(0, 1)), kilogramMetresSquared(turned(0, 2)), kilogramMetresSquared(turned(1, 2))}};
        builder.addBody({.name = k == 0 ? "b0" : k == 1 ? "b1" : k == 2 ? "b2" : k == 3 ? "b3" : k == 4 ? "b4" : "b5",
                         .mass = sub0math::kilograms(1.0 + 0.3 * f),
                         .centreOfMass = {metres(com.x), metres(com.y), metres(com.z)},
                         .inertia = inertia,
                         .joint = {.type = types[k], .parent = BodyId{parents[k]}, .placement = placement, .axis = axis}});
    }
    const Pose siteInBody = offset(0.1, -0.2, 0.3);
    builder.addSite({.name = "end", .body = BodyId{4},
                     .frame = Pose{sub0math::conjugate(frameTurn[4]), {}} * siteInBody});
    return builder.build().value();
}

State spatialState(const Model& model)
{
    State s = State::atRest(model);
    for (std::size_t i = 0; i < model.dofCount(); ++i)
    {
        const double f = static_cast<double>(i);
        s.q()[i] = 0.4 * f - 0.7;
        s.qd()[i] = 0.9 - 0.35 * f;
    }
    return s;
}

// Times of upward zero crossings of q[0], linearly interpolated.
std::vector<double> upwardCrossings(const Model& model, State s, double h, double duration)
{
    Workspace ws(model);
    std::vector<double> times;
    const auto steps = static_cast<std::uint64_t>(duration / h);
    for (std::uint64_t k = 0; k < steps; ++k)
    {
        const double before = s.q()[0];
        REQUIRE(stepSemiImplicitEuler(model, {}, s.tick, sub0math::seconds(h), s, ws));
        if (before < 0.0 && s.q()[0] >= 0.0)
        {
            times.push_back((static_cast<double>(k) + before / (before - s.q()[0])) * h);
        }
    }
    return times;
}

} // namespace

TEST_CASE("Single pendulum: exact acceleration, period and bounded energy error")
{
    const double m = 1.3;
    const double l = 0.7;
    const double ic = 0.01;
    const Model model = pendulum(m, l, ic);
    Workspace ws(model);
    State s = State::atRest(model);
    s.q()[0] = 0.3;
    s.qd()[0] = 1.1; // a revolute joint's rate does not change its own acceleration
    REQUIRE(forwardDynamics(model, s, {}, ws));
    const double pivotInertia = ic + m * l * l;
    CHECK(ws.qdd[0] == doctest::Approx(-m * G * l * std::sin(0.3) / pivotInertia).epsilon(1e-13));

    const double h = 1e-4;
    SUBCASE("small amplitude matches 2 pi sqrt(I / m g l)")
    {
        State start = State::atRest(model);
        start.q()[0] = -0.01;
        const auto crossings = upwardCrossings(model, start, h, 6.0);
        REQUIRE(crossings.size() >= 3);
        const double period = (crossings.back() - crossings.front()) / static_cast<double>(crossings.size() - 1);
        const double expected = 2.0 * std::numbers::pi * std::sqrt(pivotInertia / (m * G * l)) * (1.0 + 0.01 * 0.01 / 16.0);
        CHECK(period == doctest::Approx(expected).epsilon(1e-5));
    }
    SUBCASE("large amplitude matches the elliptic-integral period")
    {
        State start = State::atRest(model);
        start.q()[0] = -1.2;
        const auto crossings = upwardCrossings(model, start, h, 8.0);
        REQUIRE(crossings.size() >= 3);
        const double period = (crossings.back() - crossings.front()) / static_cast<double>(crossings.size() - 1);
        const double expected = 4.0 * std::sqrt(pivotInertia / (m * G * l)) * completeEllipticK(std::sin(0.6));
        CHECK(period == doctest::Approx(expected).epsilon(1e-5));
    }
    SUBCASE("energy error stays bounded")
    {
        State state = State::atRest(model);
        state.q()[0] = 1.0;
        const double e0 = computeEnergy(model, state, ws).value().mechanical().si();
        CHECK(e0 == doctest::Approx(-m * G * l * std::cos(1.0)));
        double worst = 0.0;
        for (int k = 0; k < 100000; ++k)
        {
            REQUIRE(stepSemiImplicitEuler(model, {}, state.tick, sub0math::seconds(h), state, ws));
            worst = std::max(worst, std::abs(computeEnergy(model, state, ws).value().mechanical().si() - e0));
        }
        CHECK(state.tick == 100000U);
        CHECK(worst / (m * G * l) < 2e-4);
    }
}

TEST_CASE("Double pendulum conserves energy over long runs")
{
    const Model model = doublePendulum();
    Workspace ws(model);
    State start = State::atRest(model);
    start = State::fromCoordinates({1.2, -0.6}, {0.0, 1.5});
    const double e0 = computeEnergy(model, start, ws).value().mechanical().si();
    const double scale = 1.0 * G * 1.0; // m g l of the upper link

    SUBCASE("symplectic Euler: first-order error that does not grow")
    {
        // Worst energy error over the first and the last 5 s of a 40 s run, at h and h / 2.
        const auto windows = [&](double h) {
            State s = start;
            const auto steps = static_cast<int>(40.0 / h);
            const int window = static_cast<int>(5.0 / h);
            std::array<double, 2> worst{};
            for (int k = 0; k < steps; ++k)
            {
                REQUIRE(stepSemiImplicitEuler(model, {}, s.tick, sub0math::seconds(h), s, ws));
                const double error = std::abs(computeEnergy(model, s, ws).value().mechanical().si() - e0) / scale;
                if (k < window)
                {
                    worst[0] = std::max(worst[0], error);
                }
                else if (k >= steps - window)
                {
                    worst[1] = std::max(worst[1], error);
                }
            }
            return worst;
        };
        const auto coarse = windows(1e-4);
        const auto fine = windows(5e-5);
        MESSAGE("double pendulum symplectic Euler, relative energy error, first / last 5 s of 40 s: h = 0.1 ms "
                << coarse[0] << " / " << coarse[1] << "; h = 0.05 ms " << fine[0] << " / " << fine[1]);
        CHECK(std::max(coarse[0], coarse[1]) < 1e-2);
        // First order in h over the same early trajectory. Later windows differ in h only through chaos.
        // Semi-implicit Euler is not symplectic when M depends on q, so it bounds rather than conserves:
        // conservation to round-off is RK4's job (next subcase).
        CHECK(fine[0] < 0.6 * coarse[0]);
    }
    SUBCASE("RK4 reference stepper, h = 1 ms, 20 s")
    {
        Rk4Stepper rk4(model, ws);
        State s = start;
        for (std::uint64_t tick = 0; tick < 20000U; ++tick)
        {
            REQUIRE(rk4.step(model, {}, tick, sub0math::seconds(1e-3), s, ws));
        }
        CHECK(s.tick == 20000U);
        const double error = std::abs(computeEnergy(model, s, ws).value().mechanical().si() - e0) / scale;
        MESSAGE("double pendulum RK4: relative energy error " << error);
        CHECK(error < 1e-7);
    }
}

TEST_CASE("Slider under gravity follows the exact discrete trajectory")
{
    ModelBuilder builder;
    builder.setGravity(gravityAlong(0.0, -G, 0.0));
    const double angle = std::numbers::pi / 6.0;
    builder.addBody({.name = "slider",
                     .mass = sub0math::kilograms(4.0),
                     .centreOfMass = {metres(0.2), metres(0.1), metres(0.0)},
                     .inertia = diagonalInertia(0.1, 0.2, 0.3),
                     .joint = {.type = JointType::Prismatic,
                               .placement = Pose{sub0math::Quatd::fromAxisAngle({0.0, 0.0, 1.0}, sub0math::f64::Angle::fromSi(angle)), {}},
                               .axis = {1.0, 0.0, 0.0}}});
    const Model model = *builder.build();
    Workspace ws(model);
    State s = State::atRest(model);
    s.q()[0] = 0.5;
    s.qd()[0] = 2.0;
    const double a = -G * std::sin(angle);
    REQUIRE(forwardDynamics(model, s, {}, ws));
    CHECK(ws.qdd[0] == doctest::Approx(a).epsilon(1e-15));

    // Symplectic Euler with constant acceleration: q_N = q0 + N h v0 + a h^2 N (N + 1) / 2.
    const double h = 1e-3;
    const int steps = 1000;
    for (int k = 0; k < steps; ++k)
    {
        REQUIRE(stepSemiImplicitEuler(model, {}, s.tick, sub0math::seconds(h), s, ws));
    }
    const double n = steps;
    CHECK(s.qd()[0] == doctest::Approx(2.0 + a * h * n).epsilon(1e-12));
    CHECK(s.q()[0] == doctest::Approx(0.5 + n * h * 2.0 + a * h * h * n * (n + 1.0) / 2.0).epsilon(1e-12));
}

TEST_CASE("Cart with two poles matches the hand-derived equations of motion")
{
    const CartPoles reference;
    const Model model = reference.model();
    Workspace ws(model);
    State s = State::atRest(model);
    for (const auto& [q, qd] : {std::pair{std::array{0.3, 0.4, -1.1}, std::array{0.7, -2.0, 1.3}},
                                std::pair{std::array{-1.0, 2.5, 0.2}, std::array{-0.3, 0.6, -4.0}}})
    {
        s = State::fromCoordinates(q, qd);
        REQUIRE(computeMassMatrix(model, s, ws));
        REQUIRE(computeBiasForces(model, s, ws));
        const auto m = reference.massMatrix(s);
        const auto c = reference.bias(s);
        for (std::size_t r = 0; r < 3; ++r)
        {
            CHECK(ws.bias[r] == doctest::Approx(c[r]).epsilon(1e-12).scale(1.0));
            for (std::size_t k = 0; k < 3; ++k)
            {
                CHECK(element(ws.massMatrix, 3, r, k) == doctest::Approx(m[r * 3 + k]).epsilon(1e-12).scale(1.0));
            }
        }
        // Forward dynamics solves the reference equations: M qdd = tau - c.
        const std::array<double, 3> tau{5.0, -0.4, 0.25};
        REQUIRE(forwardDynamics(model, s, tau, ws));
        for (std::size_t r = 0; r < 3; ++r)
        {
            double lhs = 0.0;
            for (std::size_t k = 0; k < 3; ++k)
            {
                lhs += m[r * 3 + k] * ws.qdd[k];
            }
            CHECK(lhs == doctest::Approx(tau[r] - c[r]).epsilon(1e-12).scale(1.0));
        }
    }
}

TEST_CASE("Mass matrix: symmetric, positive definite, and equal to RNEA column by column")
{
    for (const Model& model : {spatialTree(), planarTree()})
    {
        const std::size_t n = model.dofCount();
        Workspace ws(model);
        State s = spatialState(model);
        REQUIRE(computeMassMatrix(model, s, ws));
        const std::vector<double> m = ws.massMatrix;
        for (std::size_t r = 0; r < n; ++r)
        {
            for (std::size_t c = 0; c < n; ++c)
            {
                CHECK(element(m, n, r, c) == element(m, n, c, r));
            }
        }
        REQUIRE(factorMassMatrix(model, ws)); // fails unless every pivot is positive
        for (std::size_t r = 0; r < n; ++r)
        {
            CHECK(ws.massFactor[r * n + r] > 0.0);
        }

        // Column j of M is the change in RNEA's torque for a unit acceleration of coordinate j.
        REQUIRE(computeBiasForces(model, s, ws));
        const std::vector<double> bias = ws.bias;
        std::vector<double> qdd(n, 0.0);
        std::vector<double> tau(n, 0.0);
        for (std::size_t j = 0; j < n; ++j)
        {
            qdd.assign(n, 0.0);
            qdd[j] = 1.0;
            REQUIRE(inverseDynamics(model, s, qdd, tau, ws));
            for (std::size_t r = 0; r < n; ++r)
            {
                CHECK(tau[r] - bias[r] == doctest::Approx(element(m, n, r, j)).epsilon(1e-12).scale(1.0));
            }
        }

        // The factor solves M x = b, and forward dynamics inverts inverse dynamics.
        std::vector<double> applied(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            applied[i] = std::sin(static_cast<double>(i) + 0.5);
        }
        REQUIRE(forwardDynamics(model, s, applied, ws));
        const std::vector<double> accelerations = ws.qdd;
        REQUIRE(inverseDynamics(model, s, accelerations, tau, ws));
        for (std::size_t i = 0; i < n; ++i)
        {
            CHECK(tau[i] == doctest::Approx(applied[i]).epsilon(1e-11).scale(1.0));
        }
    }
}

TEST_CASE("Dynamics do not depend on the choice of body frames")
{
    using sub0math::Quatd;
    using sub0math::Vec3d;
    std::array<Quatd, 6> turns{};
    for (std::size_t k = 0; k < turns.size(); ++k)
    {
        turns[k] = Quatd::fromAxisAngle(sub0math::normalized(Vec3d{0.2, 1.0, static_cast<double>(k) - 2.0}),
                                        sub0math::f64::Angle::fromSi(0.5 + 0.7 * static_cast<double>(k)));
    }
    const Model plain = spatialTree();
    const Model turned = spatialTree(turns);
    Workspace a(plain);
    Workspace b(turned);
    const State s = spatialState(plain);
    const std::size_t n = plain.dofCount();

    REQUIRE(forwardDynamics(plain, s, {}, a));
    REQUIRE(forwardDynamics(turned, s, {}, b));
    for (std::size_t i = 0; i < n; ++i)
    {
        CHECK(b.bias[i] == doctest::Approx(a.bias[i]).epsilon(1e-12).scale(1.0));
        CHECK(b.qdd[i] == doctest::Approx(a.qdd[i]).epsilon(1e-11).scale(1.0));
        for (std::size_t j = 0; j < n; ++j)
        {
            CHECK(b.massMatrix[i * n + j] == doctest::Approx(a.massMatrix[i * n + j]).epsilon(1e-12).scale(1.0));
        }
    }
    CHECK(computeEnergy(turned, s, b).value().mechanical().si() == doctest::Approx(computeEnergy(plain, s, a).value().mechanical().si()).epsilon(1e-12));

    const SiteId site = *plain.findSite("end");
    const SiteKinematics ka = computeSiteKinematics(plain, s, site, a).value();
    const SiteKinematics kb = computeSiteKinematics(turned, s, site, b).value();
    CHECK(kb.pose.translation.x.si() == doctest::Approx(ka.pose.translation.x.si()));
    CHECK(kb.pose.translation.y.si() == doctest::Approx(ka.pose.translation.y.si()));
    CHECK(kb.pose.translation.z.si() == doctest::Approx(ka.pose.translation.z.si()));
    CHECK(kb.velocity.z.si() == doctest::Approx(ka.velocity.z.si()));
    CHECK(kb.angularVelocity.x.si() == doctest::Approx(ka.angularVelocity.x.si()));
}

TEST_CASE("Site Jacobian maps rates to the site velocity and matches finite differences")
{
    const Model model = spatialTree();
    const std::size_t n = model.dofCount();
    Workspace ws(model);
    const State s = spatialState(model);
    const SiteId site = *model.findSite("end");
    std::vector<double> j(3 * n);
    REQUIRE(computeSiteJacobian(model, s, site, ws, j));
    CHECK_FALSE(computeSiteJacobian(model, s, site, ws, std::span<double>(j).first(n)));

    const SiteKinematics k = computeSiteKinematics(model, s, site, ws).value();
    const std::array<double, 3> velocity{k.velocity.x.si(), k.velocity.y.si(), k.velocity.z.si()};
    for (std::size_t r = 0; r < 3; ++r)
    {
        double jqd = 0.0;
        for (std::size_t c = 0; c < n; ++c)
        {
            jqd += j[r * n + c] * s.qd()[c];
        }
        CHECK(jqd == doctest::Approx(velocity[r]).epsilon(1e-12).scale(1.0));
    }

    const double step = 1e-6;
    for (std::size_t c = 0; c < n; ++c)
    {
        State plus = s;
        State minus = s;
        plus.q()[c] += step;
        minus.q()[c] -= step;
        const auto p = computeSiteKinematics(model, plus, site, ws).value().pose.translation;
        const auto m = computeSiteKinematics(model, minus, site, ws).value().pose.translation;
        CHECK((p.x.si() - m.x.si()) / (2.0 * step) == doctest::Approx(j[c]).epsilon(1e-7).scale(1.0));
        CHECK((p.y.si() - m.y.si()) / (2.0 * step) == doctest::Approx(j[n + c]).epsilon(1e-7).scale(1.0));
        CHECK((p.z.si() - m.z.si()) / (2.0 * step) == doctest::Approx(j[2 * n + c]).epsilon(1e-7).scale(1.0));
    }
}

TEST_CASE("Model builder validates the description")
{
    const auto body = [](std::string_view name) {
        return BodyDescription{.name = name,
                               .mass = sub0math::kilograms(1.0),
                               .centreOfMass = {metres(0.0), metres(-0.5), metres(0.0)},
                               .inertia = diagonalInertia(0.01, 0.01, 0.01),
                               .joint = {}};
    };
    const auto faultOf = [](const ModelBuilder& builder) {
        const auto built = builder.build();
        REQUIRE_FALSE(built);
        CHECK(built.error().kind == FaultKind::InvalidInput);
        return built.error().detail;
    };

    {
        ModelBuilder ok;
        const BodyId root = ok.addBody(body("root"));
        BodyDescription child = body("child");
        child.joint.parent = root;
        ok.addBody(child);
        ok.addSite({.name = "s", .body = root});
        const auto built = ok.build();
        REQUIRE(built);
        CHECK(built->findBody("child") == BodyId{1});
        CHECK(built->bodyName(BodyId{0}) == "root");
        CHECK(built->parents()[1] == 0U);
        CHECK(built->findSite("s") == SiteId{0});
        CHECK_FALSE(built->findSite("missing"));
    }
    CHECK(faultOf(ModelBuilder{}) == "model has no bodies");
    {
        ModelBuilder b;
        BodyDescription d = body("a");
        d.mass = sub0math::kilograms(0.0);
        b.addBody(d);
        CHECK(faultOf(b) == "body mass must be finite and positive");
    }
    {
        ModelBuilder b;
        BodyDescription d = body("a");
        d.inertia = diagonalInertia(0.01, -0.01, 0.01);
        b.addBody(d);
        CHECK(faultOf(b) == "body inertia must be positive definite");
    }
    {
        // Ixx > Iyy + Izz: positive definite, but no mass distribution has these moments.
        ModelBuilder b;
        BodyDescription d = body("a");
        d.inertia = diagonalInertia(0.05, 0.01, 0.01);
        b.addBody(d);
        CHECK(faultOf(b) == "body inertia violates the triangle inequality (a principal moment exceeds the sum of the other two)");
    }
    {
        // The same violation in a rotated frame: the products hide it from the diagonal.
        ModelBuilder b;
        BodyDescription d = body("a");
        // Principal moments 0.01, 0.01, 0.05: in-plane eigenvalues 0.03 +- 0.02 from Ixy, Izz = 0.01.
        d.inertia = diagonalInertia(0.03, 0.03, 0.01);
        d.inertia.products.x = sub0math::kilogramMetresSquared(0.02);
        b.addBody(d);
        CHECK(faultOf(b) == "body inertia violates the triangle inequality (a principal moment exceeds the sum of the other two)");
    }
    {
        ModelBuilder b;
        b.addBody(body("a"));
        b.addBody(body("a"));
        CHECK(faultOf(b) == "body name is not unique");
    }
    {
        ModelBuilder b;
        BodyDescription d = body("a");
        d.joint.axis = {0.0, 0.0, 2.0};
        b.addBody(d);
        CHECK(faultOf(b) == "joint axis must be a finite unit vector");
    }
    {
        ModelBuilder b;
        BodyDescription d = body("a");
        d.joint.parent = BodyId{0}; // itself
        b.addBody(d);
        CHECK(faultOf(b) == "joint parent must be a body added earlier, or the world");
    }
    {
        ModelBuilder b;
        b.addBody(body("a"));
        b.addSite({.name = "s", .body = BodyId{3}});
        CHECK(faultOf(b) == "site body is not a body of this model, nor the world");
    }
    {
        ModelBuilder b;
        BodyDescription d = body("a");
        d.centreOfMass.x = metres(std::numeric_limits<double>::quiet_NaN());
        b.addBody(d);
        CHECK(faultOf(b) == "body centre of mass and inertia must be finite");
    }
}

TEST_CASE("Stepping is transactional and rejects bad input")
{
    const Model model = pendulum(1.0, 0.5, 0.01);
    Workspace ws(model);
    State s = State::atRest(model);
    s.q()[0] = 0.2;
    const State before = s;
    CHECK_FALSE(stepSemiImplicitEuler(model, {}, s.tick, sub0math::seconds(0.0), s, ws));
    const std::array<double, 2> wrongSize{1.0, 2.0};
    CHECK_FALSE(stepSemiImplicitEuler(model, {wrongSize}, s.tick, sub0math::seconds(1e-3), s, ws));
    s.qd()[0] = std::numeric_limits<double>::infinity();
    const auto r = stepSemiImplicitEuler(model, {}, s.tick, sub0math::seconds(1e-3), s, ws);
    REQUIRE_FALSE(r);
    s.qd()[0] = 0.0;
    CHECK(s == before);
    // The tick must be the state's.
    CHECK_FALSE(stepSemiImplicitEuler(model, {}, s.tick + 1U, sub0math::seconds(1e-3), s, ws));
    CHECK(s == before);
}

TEST_CASE("Every stepper has one shape and is transactional")
{
    const Model model = doublePendulum();
    Workspace ws(model);
    const State start = State::fromCoordinates({0.4, -0.3}, {0.2, 0.1});
    const auto check = [&](auto& stepper) {
        static_assert(Stepper<std::remove_cvref_t<decltype(stepper)>>);
        State s = start;
        REQUIRE(stepper.step(model, {}, 0U, sub0math::seconds(1e-3), s, ws));
        CHECK(s.tick == 1U);
        CHECK(s != start);
        const State before = s;
        const auto rejected = [&](const auto& result) {
            REQUIRE_FALSE(result);
            CHECK(result.error().kind == FaultKind::InvalidInput);
            CHECK(s == before);
        };
        rejected(stepper.step(model, {}, 0U, sub0math::seconds(1e-3), s, ws)); // wrong tick
        rejected(stepper.step(model, {}, 1U, sub0math::seconds(-1e-3), s, ws));
        const std::array<double, 3> wrongSize{};
        rejected(stepper.step(model, {wrongSize}, 1U, sub0math::seconds(1e-3), s, ws));
        s.qd()[1] = std::numeric_limits<double>::infinity();
        const State poisoned = s;
        REQUIRE_FALSE(stepper.step(model, {}, 1U, sub0math::seconds(1e-3), s, ws));
        CHECK(s == poisoned);
    };
    SemiImplicitEuler euler;
    check(euler);
    Rk4Stepper rk4(model, ws);
    check(rk4);
    MidpointStepper midpoint(model, ws);
    check(midpoint);
}

TEST_CASE("RK steppers converge at their order and follow a change of forcing")
{
    const Model model = doublePendulum();
    Workspace ws(model);
    const State start = State::fromCoordinates({0.9, -0.4}, {0.0, 0.8});
    // A torque applied from tick 50 onward; the cached first stage must not survive the change.
    const std::array<double, 2> torque{0.7, -0.2};
    const auto run = [&](auto stepper, double h) {
        State s = start;
        const auto steps = static_cast<std::uint64_t>(std::llround(1.0 / h));
        for (std::uint64_t tick = 0; tick < steps; ++tick)
        {
            const Forcing forcing{tick * 2U < steps ? std::span<const double>{} : std::span<const double>{torque}};
            REQUIRE(stepper.step(model, forcing, tick, sub0math::seconds(h), s, ws));
        }
        return s;
    };
    const State reference = run(Rk4Stepper(model, ws), 1e-4);
    const auto error = [&](const State& s) {
        double sum = 0.0;
        for (std::size_t i = 0; i < s.values.size(); ++i)
        {
            sum += (s.values[i] - reference.values[i]) * (s.values[i] - reference.values[i]);
        }
        return std::sqrt(sum);
    };
    const double rk4Coarse = error(run(Rk4Stepper(model, ws), 0.004));
    const double rk4Fine = error(run(Rk4Stepper(model, ws), 0.002));
    const double midCoarse = error(run(MidpointStepper(model, ws), 0.004));
    const double midFine = error(run(MidpointStepper(model, ws), 0.002));
    MESSAGE("RK4 order " << std::log2(rk4Coarse / rk4Fine) << ", midpoint order " << std::log2(midCoarse / midFine));
    CHECK(std::log2(rk4Coarse / rk4Fine) == doctest::Approx(4.0).epsilon(0.1));
    CHECK(std::log2(midCoarse / midFine) == doctest::Approx(2.0).epsilon(0.1));
}

TEST_CASE("Rk4Stepper matches the fixed-size Rk4 on the same system")
{
    const Model model = doublePendulum();
    Workspace ws(model);
    const State start = State::fromCoordinates({1.2, -0.6}, {0.0, 1.5});
    // The fixed-size integrator over the same adapter, through arrays.
    struct Fixed
    {
        ModelOde ode;
        std::expected<void, Fault> evaluate(double t, const StateVector<4>& x, StateVector<4>& dxdt)
        {
            return ode.evaluate(t, x, dxdt);
        }
    };
    Workspace fixedWs(model);
    Fixed fixed{ModelOde(model, fixedWs)};
    Rk4<4> rk4;
    StateVector<4> x{};
    std::copy(start.values.begin(), start.values.end(), x.begin());
    Rk4Stepper stepper(model, ws);
    State s = start;
    for (std::uint64_t tick = 0; tick < 500U; ++tick)
    {
        REQUIRE(rk4.stepTick(fixed, tick, sub0math::seconds(1e-3), x));
        REQUIRE(stepper.step(model, {}, tick, sub0math::seconds(1e-3), s, ws));
    }
    for (std::size_t i = 0; i < 4; ++i)
    {
        CHECK(s.values[i] == x[i]);
    }
}

TEST_CASE("Entry points reject a workspace built for another model and a mis-sized state")
{
    const Model tree = planarTree();
    // Same number of coordinates, different topology: a size check alone would let this through and
    // leave stale ancestor entries in the mass matrix.
    ModelBuilder chainBuilder;
    chainBuilder.setGravity(gravityAlong(0.0, -G, 0.0));
    BodyId parent = WorldBody;
    for (std::size_t i = 0; i < PlanarTreeDof; ++i)
    {
        parent = chainBuilder.addBody(planarLink("l" + std::to_string(i), parent, offset(0.0, -0.4, 0.0), 1.0, 0.4));
    }
    chainBuilder.addSite({.name = "s", .body = BodyId{0}});
    const Model chain = *chainBuilder.build();
    REQUIRE(chain.dofCount() == tree.dofCount());
    CHECK(chain.id() != tree.id());

    Workspace treeWs(tree);
    Workspace chainWs(chain);
    const State s = State::atRest(tree);
    std::vector<double> tau(tree.dofCount());
    std::vector<double> x(tree.dofCount());
    std::vector<double> jacobian(3 * tree.dofCount());

    const auto rejected = [](const auto& result) {
        REQUIRE_FALSE(result);
        CHECK(result.error().kind == FaultKind::InvalidInput);
    };
    rejected(computeMassMatrix(tree, s, chainWs));
    rejected(computeBiasForces(tree, s, chainWs));
    rejected(forwardDynamics(tree, s, {}, chainWs));
    rejected(inverseDynamics(tree, s, x, tau, chainWs));
    rejected(updateKinematics(tree, s, chainWs));
    rejected(computeSiteKinematics(tree, s, SiteId{0}, chainWs));
    rejected(computeSiteJacobian(tree, s, SiteId{0}, chainWs, jacobian));
    rejected(computeEnergy(tree, s, chainWs));
    rejected(factorMassMatrix(tree, chainWs));
    rejected(solveMassMatrix(tree, chainWs, x));
    State stepped = s;
    rejected(stepSemiImplicitEuler(tree, {}, stepped.tick, sub0math::seconds(1e-3), stepped, chainWs));
    CHECK(stepped == s);

    // A different dof count is caught the same way, and a copy of the model shares its workspace.
    const Model single = pendulum(1.0, 0.5, 0.01);
    Workspace singleWs(single);
    rejected(computeMassMatrix(tree, s, singleWs));
    const Model copy = tree;
    CHECK(computeMassMatrix(copy, s, treeWs));

    // State sizes.
    State shortQ = s;
    shortQ.values.pop_back();
    State shortQd = s;
    shortQd.values.clear();
    for (const State* bad : {&shortQ, &shortQd})
    {
        rejected(computeMassMatrix(tree, *bad, treeWs));
        rejected(computeBiasForces(tree, *bad, treeWs));
        rejected(forwardDynamics(tree, *bad, {}, treeWs));
        rejected(inverseDynamics(tree, *bad, x, tau, treeWs));
        rejected(updateKinematics(tree, *bad, treeWs));
        rejected(computeSiteKinematics(tree, *bad, SiteId{0}, treeWs));
        rejected(computeSiteJacobian(tree, *bad, SiteId{0}, treeWs, jacobian));
        rejected(computeEnergy(tree, *bad, treeWs));
    }
    State badStep = shortQ;
    rejected(stepSemiImplicitEuler(tree, {}, badStep.tick, sub0math::seconds(1e-3), badStep, treeWs));
    CHECK(badStep == shortQ);

    // Sites and span sizes.
    rejected(computeSiteKinematics(tree, s, SiteId{1}, treeWs));
    rejected(computeSiteJacobian(tree, s, SiteId{1}, treeWs, jacobian));
    rejected(solveMassMatrix(tree, treeWs, std::span<double>(x).first(3)));
    CHECK(computeSiteKinematics(tree, s, *tree.findSite("tip"), treeWs));
}

TEST_CASE("Mass matrix entries outside ancestor pairs are zero")
{
    const Model model = planarTree();
    Workspace ws(model);
    for (const double v : ws.massMatrix)
    {
        CHECK(v == 0.0);
    }
    State s = State::atRest(model);
    for (std::size_t i = 0; i < model.dofCount(); ++i)
    {
        s.q()[i] = 0.1 * static_cast<double>(i) - 0.3;
    }
    REQUIRE(computeMassMatrix(model, s, ws));
    const auto parents = model.parents();
    const std::size_t n = model.dofCount();
    const auto isAncestor = [&](std::size_t a, std::size_t b) { // a is an ancestor of b, or equal
        for (std::size_t k = b; k != NoParent; k = parents[k])
        {
            if (k == a) return true;
        }
        return false;
    };
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            if (!isAncestor(i, j) && !isAncestor(j, i))
            {
                CHECK(element(ws.massMatrix, n, i, j) == 0.0);
            }
        }
    }
}

TEST_CASE("Model builder accepts inertias on the edge of the triangle inequality")
{
    const auto build = [](const RotationalInertia& inertia) {
        ModelBuilder b;
        b.addBody({.name = "a",
                   .mass = sub0math::kilograms(1.0),
                   .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                   .inertia = inertia,
                   .joint = {}});
        return b.build();
    };
    // Thin rod along x: Ixx ~ 0, Iyy = Izz. Valid down to a small relative Ixx.
    CHECK(build(diagonalInertia(1e-9, 0.1, 0.1)));
    // Thin disc in the xy plane: Izz = Ixx + Iyy exactly (a planar lamina), in any axis order.
    CHECK(build(diagonalInertia(0.02, 0.02, 0.04)));
    CHECK(build(diagonalInertia(0.04, 0.02, 0.02)));
    CHECK(build(diagonalInertia(0.02, 0.04, 0.02)));
    // Rounding in how a user computed a lamina's moments must not reject it.
    CHECK(build(diagonalInertia(0.02, 0.02, 0.04 * (1.0 + 1e-12))));
    // A product-of-inertia tensor with principal moments 0.01, 0.05, 0.04 (lamina in a turned frame).
    RotationalInertia turned = diagonalInertia(0.03, 0.03, 0.04);
    turned.products.x = sub0math::kilogramMetresSquared(0.02);
    CHECK(build(turned));
    // Just beyond the edge is rejected.
    CHECK_FALSE(build(diagonalInertia(0.02, 0.02, 0.04 * (1.0 + 1e-6))));
    // A sphere and a cube are comfortably inside.
    CHECK(build(diagonalInertia(0.1, 0.1, 0.1)));
}
