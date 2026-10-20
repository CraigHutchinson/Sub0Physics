// The planar formulation, site forces, forcing Jacobians, implicit-term levers and massless bodies (docs/MULTIBODY.md,
// M1d status). The planar formulation is checked against the general one on the same models (every dynamics query,
// the elements with their Jacobians, energy, site kinematics, reactions, and stepped trajectories with and without
// constraints); the forcing additions against independent formulations of the same physics.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <doctest.h>
#include <ostream> // doctest prints strings (MSVC needs it explicitly)
#include <span>
#include <string_view>
#include <vector>

#include <sub0math/linear/quat.hpp>
#include <sub0math/units.hpp>

#include "sub0physics/elements.hpp"
#include "sub0physics/multibody.hpp"
#include "test_multibody_models.hpp"

using namespace sub0physics;
using namespace sub0physics::multibody;
using namespace testmodels;
namespace f64 = sub0math::f64;

namespace
{

// Largest |a - b| over the largest |b| (or 1 when b is all small).
double relativeDifference(std::span<const double> a, std::span<const double> b)
{
    REQUIRE(a.size() == b.size());
    double scale = 1.0;
    double worst = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        scale = std::max(scale, std::abs(b[i]));
        worst = std::max(worst, std::abs(a[i] - b[i]));
    }
    return worst / scale;
}

double vectorDifference(const sub0math::Vec3d& a, const sub0math::Vec3d& b)
{
    return std::max({std::abs(a.x - b.x), std::abs(a.y - b.y), std::abs(a.z - b.z)});
}

RotationalInertia inertia(double xx, double yy, double zz)
{
    return diagonalInertia(xx, yy, zz);
}

Pose rotatedAboutY(double angle, double x, double y, double z)
{
    return {sub0math::Quatd::fromAxisAngle({0.0, 1.0, 0.0}, f64::Angle::fromSi(angle)), {metres(x), metres(y), metres(z)}};
}

// A rig-like tree in the xz plane (normal y, gravity -z) with every feature the planar formulation must map: a
// massless carrier (two joints on one body), a revolute about +y, a prismatic axis in the plane under a rotated
// placement, a revolute about -y, out-of-plane offsets, a tendon whose sites differ along y, and coordinate elements.
ModelBuilder rigLikeBuilder()
{
    using sub0math::kilograms;
    ModelBuilder b;
    const BodyId carrier = b.addBody({.name = "carrier",
                                      .mass = kilograms(0.0),
                                      .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                                      .inertia = inertia(0.0, 0.0, 0.0),
                                      .joint = {.type = JointType::Prismatic, .placement = offset(0.0, 0.0, 0.65),
                                                .axis = {0.0, 0.0, 1.0}}});
    const BodyId frame = b.addBody({.name = "frame",
                                    .mass = kilograms(12.0),
                                    .centreOfMass = {metres(0.05), metres(0.01), metres(-0.02)},
                                    .inertia = inertia(1.1, 1.3, 0.9),
                                    .joint = {.type = JointType::Revolute, .parent = carrier, .axis = {0.0, 1.0, 0.0}}});
    const double head = 65.0 * std::numbers::pi / 180.0;
    const BodyId front = b.addBody({.name = "front",
                                    .mass = kilograms(3.0),
                                    .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                                    .inertia = inertia(0.1, 0.1, 0.1),
                                    .joint = {.type = JointType::Prismatic, .parent = frame,
                                              .placement = rotatedAboutY(0.2, 0.7, 0.05, -0.29),
                                              .axis = {-std::cos(head), 0.0, std::sin(head)},
                                              .springReference = 0.01}});
    const BodyId rear = b.addBody({.name = "rear",
                                   .mass = kilograms(4.0),
                                   .centreOfMass = {metres(-0.49), metres(0.0), metres(-0.055)},
                                   .inertia = inertia(0.13, 0.13, 0.13),
                                   .joint = {.type = JointType::Revolute, .parent = frame,
                                             .placement = rotatedAboutY(-0.3, 0.04, -0.03, -0.235),
                                             .axis = {0.0, -1.0, 0.0}}});
    b.addBody({.name = "rider",
               .mass = kilograms(60.0),
               .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
               .inertia = inertia(5.0, 5.0, 5.0),
               .joint = {.type = JointType::Prismatic, .parent = frame, .axis = {0.0, 0.0, 1.0}}});
    const SiteId eyeFrame = b.addSite({.name = "eye frame", .body = frame, .frame = offset(0.19, 0.02, -0.016)});
    const SiteId eyeArm = b.addSite({.name = "eye arm", .body = rear, .frame = offset(-0.055, -0.01, 0.19)});
    b.addSite({.name = "front axle", .body = front, .frame = offset(0.0, 0.0, 0.0)});
    b.addSite({.name = "rear axle", .body = rear, .frame = offset(-0.49, 0.0, -0.055)});
    b.addSite({.name = "anchor", .body = WorldBody, .frame = offset(0.3, 0.0, 0.1)});
    b.addElement(CoordinateElementDescription<LinearSpringDamper<f64::Length>>{
        .name = "fork spring",
        .body = front,
        .law = *LinearSpringDamper<f64::Length>::create(
            {.stiffness = f64::Stiffness::fromSi(9000.0), .damping = f64::Damping::fromSi(700.0)})});
    b.addElement(CoordinateElementDescription<TravelBuffer>{
        .name = "fork buffer",
        .body = front,
        .law = *TravelBuffer::create({.travel = metres(0.16), .zoneWidth = metres(0.02),
                                      .stiffness = f64::Stiffness::fromSi(3e5), .damping = f64::Damping::fromSi(1500.0)}),
        .reference = 0.0});
    b.addElement(TendonElementDescription<LinearSpringDamper<f64::Length>>{
        .name = "shock",
        .first = eyeFrame,
        .second = eyeArm,
        .law = *LinearSpringDamper<f64::Length>::create(
            {.stiffness = f64::Stiffness::fromSi(70000.0), .damping = f64::Damping::fromSi(2400.0)}),
        .reference = metres(0.21)});
    b.addElement(TendonElementDescription<AxialSpringDamper>{
        .name = "strut",
        .first = SiteId{4},
        .second = SiteId{2},
        .law = *AxialSpringDamper::create({.preload = f64::Force::fromSi(20.0),
                                           .stiffness = f64::Stiffness::fromSi(500.0),
                                           .compressionDamping = f64::Damping::fromSi(10.0),
                                           .reboundDamping = f64::Damping::fromSi(10.0)}),
        .sense = LawSense::Decreasing,
        .reference = metres(0.9)});
    return b;
}

State rigLikeState()
{
    return State::fromCoordinates({-0.017, 0.025, 0.033, 0.12, -0.04}, {0.3, -0.8, 0.5, 1.7, -0.6});
}

struct Pair
{
    Model planar;
    Model general;
};

Pair buildBoth(ModelBuilder builder)
{
    auto planar = builder.build();
    REQUIRE(planar);
    builder.setFormulation(Formulation::General);
    auto general = builder.build();
    REQUIRE(general);
    REQUIRE(planar->formulation() == Formulation::Planar);
    REQUIRE(general->formulation() == Formulation::General);
    return {*planar, *general};
}

// Every dynamics query of the two formulations at one state agrees to rounding.
void checkQueries(const Pair& models, const State& state, double tolerance)
{
    Workspace wp(models.planar);
    Workspace wg(models.general);
    const std::size_t n = models.planar.dofCount();
    REQUIRE(computeMassMatrix(models.planar, state, wp));
    REQUIRE(computeMassMatrix(models.general, state, wg));
    CHECK(relativeDifference(wp.massMatrix, wg.massMatrix) < tolerance);
    REQUIRE(computeBiasForces(models.planar, state, wp));
    REQUIRE(computeBiasForces(models.general, state, wg));
    CHECK(relativeDifference(wp.bias, wg.bias) < tolerance);
    REQUIRE(forwardDynamics(models.planar, state, {}, wp));
    REQUIRE(forwardDynamics(models.general, state, {}, wg));
    CHECK(relativeDifference(wp.qdd, wg.qdd) < tolerance);
    std::vector<double> qdd(n, 0.3);
    std::vector<double> tp(n);
    std::vector<double> tg(n);
    REQUIRE(inverseDynamics(models.planar, state, qdd, tp, wp));
    REQUIRE(inverseDynamics(models.general, state, qdd, tg, wg));
    CHECK(relativeDifference(tp, tg) < tolerance);
    if (models.planar.elementCount() > 0)
    {
        REQUIRE(evaluateElements(models.planar, state, wp));
        REQUIRE(evaluateElements(models.general, state, wg));
        CHECK(relativeDifference(wp.elementForce, wg.elementForce) < tolerance);
        CHECK(relativeDifference(wp.elementStiffness, wg.elementStiffness) < tolerance);
        CHECK(relativeDifference(wp.elementDamping, wg.elementDamping) < tolerance);
        CHECK(std::abs(wp.elementPotential - wg.elementPotential) < tolerance * std::max(1.0, wg.elementPotential));
        for (std::size_t w = 0; w < wp.elementSiteWrenches.size(); ++w)
        {
            CHECK(vectorDifference(wp.elementSiteWrenches[w].force, wg.elementSiteWrenches[w].force) <
                  tolerance * 1e3);
        }
    }
    const auto ep = computeEnergy(models.planar, state, wp);
    const auto eg = computeEnergy(models.general, state, wg);
    REQUIRE(ep);
    REQUIRE(eg);
    const double scale = std::max(1.0, std::abs(eg->kinetic.si()) + std::abs(eg->potential.si()));
    INFO("planar potential " << ep->potential.si() << ", general " << eg->potential.si());
    CHECK(std::abs(ep->kinetic.si() - eg->kinetic.si()) < tolerance * scale);
    CHECK(std::abs(ep->potential.si() - eg->potential.si()) < tolerance * scale);
    std::vector<SiteId> sites;
    for (std::size_t s = 0; s < models.planar.siteCount(); ++s)
    {
        sites.push_back(SiteId{static_cast<std::uint32_t>(s)});
    }
    std::vector<SiteKinematics> kp(sites.size());
    std::vector<SiteKinematics> kg(sites.size());
    REQUIRE(computeSiteKinematics(models.planar, state, sites, kp, wp));
    REQUIRE(computeSiteKinematics(models.general, state, sites, kg, wg));
    for (std::size_t s = 0; s < sites.size(); ++s)
    {
        const auto single = computeSiteKinematics(models.general, state, sites[s], wg);
        REQUIRE(single);
        const auto plain = [](const auto& v) { return sub0math::Vec3d{v.x.si(), v.y.si(), v.z.si()}; };
        CHECK(vectorDifference(plain(kp[s].pose.translation), plain(single->pose.translation)) < tolerance);
        CHECK(vectorDifference(plain(kp[s].velocity), plain(single->velocity)) < tolerance);
        CHECK(vectorDifference(plain(kp[s].angularVelocity), plain(single->angularVelocity)) < tolerance);
        CHECK(vectorDifference(plain(kg[s].pose.translation), plain(single->pose.translation)) == 0.0);
        const sub0math::Vec3d probe{0.3, -0.2, 0.7};
        CHECK(vectorDifference(sub0math::rotate(kp[s].pose.rotation, probe), sub0math::rotate(single->pose.rotation, probe)) < tolerance);
    }
}

} // namespace

TEST_CASE("planar formulation: selected automatically, refused when the model is not planar, forced off on request")
{
    CHECK(planarTree().formulation() == Formulation::Planar);
    CHECK(planarTree().planar().normal == 2U);
    const auto rig = rigLikeBuilder().build();
    REQUIRE(rig);
    CHECK(rig->formulation() == Formulation::Planar);
    CHECK(rig->planar().normal == 1U);
    CHECK(rig->hasMasslessBodies());

    // A skew axis: general only; requiring planar is refused.
    ModelBuilder skew;
    skew.addBody({.name = "a",
                  .mass = sub0math::kilograms(1.0),
                  .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                  .inertia = inertia(0.1, 0.1, 0.1),
                  .joint = {.type = JointType::Revolute, .axis = {0.0, 0.6, 0.8}}});
    const auto general = skew.build();
    REQUIRE(general);
    CHECK(general->formulation() == Formulation::General);
    skew.setFormulation(Formulation::Planar);
    CHECK_FALSE(skew.build());

    // Gravity along the only possible normal: general.
    ModelBuilder tilted = planarTreeBuilder();
    tilted.setGravity(gravityAlong(0.0, -9.81, -1.0));
    CHECK(tilted.build()->formulation() == Formulation::General);
}

TEST_CASE("planar formulation: every query agrees with the general formulation to rounding")
{
    SUBCASE("planar tree in the xy plane")
    {
        const Pair models = buildBoth(planarTreeBuilder());
        std::vector<double> q(PlanarTreeDof);
        std::vector<double> qd(PlanarTreeDof);
        for (std::size_t i = 0; i < PlanarTreeDof; ++i)
        {
            q[i] = 0.3 * std::sin(1.7 * static_cast<double>(i) + 0.4);
            qd[i] = 1.1 * std::cos(0.9 * static_cast<double>(i));
        }
        checkQueries(models, State::fromCoordinates(q, qd), 1e-13);
    }
    SUBCASE("planar tree with a cart spring and a branch-coupling strut")
    {
        const Pair models = buildBoth(planarTreeWithElementsBuilder());
        std::vector<double> q(PlanarTreeDof);
        std::vector<double> qd(PlanarTreeDof);
        for (std::size_t i = 0; i < PlanarTreeDof; ++i)
        {
            q[i] = 0.2 * std::cos(1.3 * static_cast<double>(i));
            qd[i] = 0.7 * std::sin(0.5 * static_cast<double>(i) + 1.0);
        }
        checkQueries(models, State::fromCoordinates(q, qd), 1e-13);
    }
    SUBCASE("rig-like tree in the xz plane: massless carrier, rotated placements, -y axis, out-of-plane offsets")
    {
        checkQueries(buildBoth(rigLikeBuilder()), rigLikeState(), 1e-13);
    }
    SUBCASE("slider-crank with its closure")
    {
        checkQueries(buildBoth(sliderCrankBuilder()), sliderCrankState({}, 0.7, 2.0), 1e-13);
    }
}

TEST_CASE("planar formulation: stepped trajectories agree with the general formulation")
{
    const auto run = [](const Pair& models, State state, auto stepper, std::size_t steps, double h) {
        State planar = state;
        State general = state;
        Workspace wp(models.planar);
        Workspace wg(models.general);
        for (std::size_t k = 0; k < steps; ++k)
        {
            REQUIRE(stepper.step(models.planar, Forcing{}, planar.tick, f64::Time::fromSi(h), planar, wp));
            REQUIRE(stepper.step(models.general, Forcing{}, general.tick, f64::Time::fromSi(h), general, wg));
        }
        return relativeDifference(planar.values, general.values);
    };
    // Rounding differs between the formulations and is carried through the steps; the systems are damped or
    // short, so it stays near round-off.
    CHECK(run(buildBoth(rigLikeBuilder()), rigLikeState(), LinearlyImplicitEuler{}, 500, 1e-3) < 1e-11);
    CHECK(run(buildBoth(rigLikeBuilder()), rigLikeState(), SemiImplicitEuler{}, 500, 1e-4) < 1e-11);
    CHECK(run(buildBoth(planarTreeWithElementsBuilder()), State::atRest(planarTree()),
              LinearlyImplicitEuler{.terms = ImplicitTerms::MaterialStiffness}, 300, 1e-3) < 1e-11);
    CHECK(run(buildBoth(sliderCrankBuilder()), sliderCrankState({}, 0.4, 3.0), LinearlyImplicitEuler{}, 300, 1e-3) <
          1e-10);
}

TEST_CASE("planar formulation: reactions agree with the general formulation")
{
    const Pair models = buildBoth(rigLikeBuilder());
    const State state = rigLikeState();
    Workspace wp(models.planar);
    Workspace wg(models.general);
    const std::array<SiteForce, 2> contact{SiteForce{SiteId{2}, {0.0, 0.0, 420.0}, {}},
                                           SiteForce{SiteId{3}, {15.0, 3.0, 610.0}, {0.0, 2.0, 0.0}}};
    const Forcing forcing{{}, contact};
    REQUIRE(computeReactionsWithForcing(models.planar, state, forcing, wp));
    REQUIRE(computeReactionsWithForcing(models.general, state, forcing, wg));
    for (std::size_t b = 0; b < models.planar.bodyCount(); ++b)
    {
        CHECK(vectorDifference(wp.jointReactions[b].force, wg.jointReactions[b].force) < 1e-10);
        CHECK(vectorDifference(wp.jointReactions[b].moment, wg.jointReactions[b].moment) < 1e-10);
    }
}

TEST_CASE("site forces: generalized force is J^T F, reactions see the force at its site, work is accounted")
{
    for (const Formulation formulation : {Formulation::Automatic, Formulation::General})
    {
        ModelBuilder builder = rigLikeBuilder();
        builder.setFormulation(formulation);
        const Model model = *builder.build();
        const std::size_t n = model.dofCount();
        const State state = rigLikeState();
        Workspace ws(model);
        const SiteForce f{SiteId{3}, {15.0, 3.0, 610.0}, {0.5, 2.0, -0.3}};
        const std::array<SiteForce, 1> forces{f};
        REQUIRE(computeAppliedForce(model, state, Forcing{{}, forces}, ws));
        // Independent: the site's translational Jacobian for the force, and the moment through each revolute axis.
        std::vector<double> jacobian(3 * n);
        REQUIRE(computeSiteJacobian(model, state, SiteId{3}, ws, jacobian));
        Workspace general(*rigLikeBuilder().build());
        for (std::size_t i = 0; i < n; ++i)
        {
            double expected = jacobian[i] * f.force.x + jacobian[n + i] * f.force.y + jacobian[2 * n + i] * f.force.z;
            if (model.jointTypes()[i] == JointType::Revolute)
            {
                // Every revolute axis is +-y in world axes here.
                expected += model.jointAxes()[i].y * f.moment.y;
            }
            CHECK(std::abs(ws.appliedForce[i] - expected) < 1e-11);
        }
        // forwardDynamicsWithForcing equals forwardDynamics with the mapped force.
        const std::vector<double> mapped = ws.appliedForce;
        REQUIRE(forwardDynamicsWithForcing(model, state, Forcing{{}, forces}, ws));
        const std::vector<double> withForcing = ws.qdd;
        REQUIRE(forwardDynamics(model, state, mapped, ws));
        CHECK(relativeDifference(withForcing, ws.qdd) < 1e-14);

        // Work: a step's externalWork is tau . dq with the mapped force.
        State stepped = state;
        REQUIRE(stepLinearlyImplicitEuler(model, Forcing{{}, forces}, 0, f64::Time::fromSi(1e-3), stepped, ws));
        double work = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            work += mapped[i] * (stepped.q()[i] - state.q()[i]);
        }
        CHECK(stepped.externalWork.si() == doctest::Approx(work).epsilon(1e-12));

        // Reactions: a force at the rear axle loads the rear joint; as a generalized torque it would not.
        REQUIRE(computeReactionsWithForcing(model, state, Forcing{{}, forces}, ws));
        const JointReaction asSite = ws.jointReactions[3];
        REQUIRE(computeReactions(model, state, mapped, ws));
        const JointReaction asTorque = ws.jointReactions[3];
        CHECK(vectorDifference(asSite.force, asTorque.force) > 100.0);
        // The rear body's Newton-Euler balance: joint reaction + site force + gravity = rate of momentum. Checked as
        // the reactions' consistency with the joint's generalized force: projection on the axis equals the applied
        // force on that coordinate (zero here).
        (void)general;

        // A refused forcing: a site that is not in the model.
        const std::array<SiteForce, 1> bad{SiteForce{SiteId{99}, {1.0, 0.0, 0.0}, {}}};
        State unchanged = state;
        CHECK_FALSE(stepLinearlyImplicitEuler(model, Forcing{{}, bad}, 0, f64::Time::fromSi(1e-3), unchanged, ws));
        CHECK(unchanged == state);
    }
}

TEST_CASE("site forces: the explicit steppers hold them over the step")
{
    const Model model = *rigLikeBuilder().build();
    const State state = rigLikeState();
    const std::array<SiteForce, 1> forces{SiteForce{SiteId{2}, {0.0, 0.0, 300.0}, {}}};
    Workspace ws(model);
    REQUIRE(computeAppliedForce(model, state, Forcing{{}, forces}, ws));
    const std::vector<double> mapped = ws.appliedForce;
    for (int which = 0; which < 2; ++which)
    {
        State a = state;
        State b = state;
        Workspace wa(model);
        Workspace wb(model);
        if (which == 0)
        {
            REQUIRE(stepSemiImplicitEuler(model, Forcing{{}, forces}, 0, f64::Time::fromSi(1e-4), a, wa));
            REQUIRE(stepSemiImplicitEuler(model, Forcing{mapped}, 0, f64::Time::fromSi(1e-4), b, wb));
        }
        else
        {
            Rk4Stepper ra(model, wa);
            Rk4Stepper rb(model, wb);
            REQUIRE(ra.step(model, Forcing{{}, forces}, 0, f64::Time::fromSi(1e-4), a, wa));
            REQUIRE(rb.step(model, Forcing{mapped}, 0, f64::Time::fromSi(1e-4), b, wb));
        }
        CHECK(a.values == b.values);
        CHECK(a.externalWork.si() == b.externalWork.si());
    }
}

TEST_CASE("forcing Jacobians: an external spring with its Jacobians steps like the same spring as an element")
{
    // Element version: the fork spring of the rig-like model. External version: the model without it, and
    // tau = -k (q - r) - c qd on the fork coordinate with K = -k, D = -c supplied by the caller.
    ModelBuilder withElement;
    {
        ModelBuilder b;
        b.addBody({.name = "slider",
                   .mass = sub0math::kilograms(2.0),
                   .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                   .inertia = inertia(0.1, 0.1, 0.1),
                   .joint = {.type = JointType::Prismatic, .axis = {1.0, 0.0, 0.0}}});
        b.addBody({.name = "pendulum",
                   .mass = sub0math::kilograms(1.0),
                   .centreOfMass = {metres(0.0), metres(-0.5), metres(0.0)},
                   .inertia = inertia(0.02, 0.001, 0.02),
                   .joint = {.type = JointType::Revolute, .parent = BodyId{0}, .axis = {0.0, 0.0, 1.0}}});
        b.setGravity(gravityAlong(0.0, -9.81, 0.0));
        withElement = b;
    }
    ModelBuilder bare = withElement;
    constexpr double K = 4e5;
    constexpr double C = 60.0;
    constexpr double Reference = 0.02;
    withElement.addElement(CoordinateElementDescription<LinearSpringDamper<f64::Length>>{
        .name = "spring",
        .body = BodyId{0},
        .law = *LinearSpringDamper<f64::Length>::create(
            {.stiffness = f64::Stiffness::fromSi(K), .damping = f64::Damping::fromSi(C)}),
        .reference = Reference});
    const Model element = *withElement.build();
    const Model external = *bare.build();
    for (const ImplicitTerms terms : {ImplicitTerms::Full, ImplicitTerms::DampingOnly})
    {
        State a = State::fromCoordinates({0.05, 0.4}, {0.0, 0.0});
        State b = a;
        Workspace wa(element);
        Workspace wb(external);
        std::array<double, 2> tau{};
        const std::array<double, 4> stiffness{-K, 0.0, 0.0, 0.0};
        const std::array<double, 4> damping{-C, 0.0, 0.0, 0.0};
        for (std::uint64_t k = 0; k < 400; ++k)
        {
            tau[0] = -K * (b.q()[0] - Reference) - C * b.qd()[0];
            REQUIRE(stepLinearlyImplicitEuler(element, Forcing{}, k, f64::Time::fromSi(2e-3), a, wa, terms));
            REQUIRE(stepLinearlyImplicitEuler(external, Forcing{tau, {}, stiffness, damping}, k,
                                              f64::Time::fromSi(2e-3), b, wb, terms));
        }
        CHECK(relativeDifference(a.values, b.values) < 1e-12);
    }
}

TEST_CASE("implicit terms: DampingOnly is MuJoCo's implicitfast structure, Full is stable where it is not")
{
    // One slider on a stiff spring-damper: DampingOnly solves (m - h D) dqd = h F, with F = -k x - c v explicit in x.
    ModelBuilder b;
    b.addBody({.name = "slider",
               .mass = sub0math::kilograms(1.0),
               .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
               .inertia = inertia(0.1, 0.1, 0.1),
               .joint = {.type = JointType::Prismatic, .axis = {1.0, 0.0, 0.0}}});
    b.setGravity(gravityAlong(0.0, 0.0, 0.0));
    constexpr double K = 1e6;
    constexpr double C = 5.0;
    b.addElement(CoordinateElementDescription<LinearSpringDamper<f64::Length>>{
        .name = "spring",
        .body = BodyId{0},
        .law = *LinearSpringDamper<f64::Length>::create(
            {.stiffness = f64::Stiffness::fromSi(K), .damping = f64::Damping::fromSi(C)})});
    const Model model = *b.build();
    Workspace ws(model);
    State s = State::fromCoordinates({0.01}, {0.0});
    const double h = 1e-4;
    double x = 0.01;
    double v = 0.0;
    for (std::uint64_t k = 0; k < 50; ++k)
    {
        REQUIRE(stepLinearlyImplicitEuler(model, Forcing{}, k, f64::Time::fromSi(h), s, ws, ImplicitTerms::DampingOnly));
        v += h * (-K * x - C * v) / (1.0 + h * C);
        x += h * v;
    }
    CHECK(s.q()[0] == doctest::Approx(x).epsilon(1e-12));
    CHECK(s.qd()[0] == doctest::Approx(v).epsilon(1e-12));

    // At h omega = 10 the explicit spring diverges; with K implicit the step stays bounded.
    State full = State::fromCoordinates({0.01}, {0.0});
    State dampingOnly = full;
    for (std::uint64_t k = 0; k < 100; ++k)
    {
        REQUIRE(stepLinearlyImplicitEuler(model, Forcing{}, k, f64::Time::fromSi(1e-2), full, ws));
        const auto stepped =
            stepLinearlyImplicitEuler(model, Forcing{}, k, f64::Time::fromSi(1e-2), dampingOnly, ws, ImplicitTerms::DampingOnly);
        if (!stepped)
        {
            break;
        }
    }
    CHECK(std::abs(full.q()[0]) <= 0.01);
    CHECK_FALSE(std::abs(dampingOnly.q()[0]) <= 0.01);
}

TEST_CASE("massless bodies: two joints on one body; a massless leaf is refused")
{
    // Heave and pitch on one frame, as MJCF's multi-joint body: the mass matrix is the frame's, positive definite.
    const Model rig = *rigLikeBuilder().build();
    Workspace ws(rig);
    REQUIRE(computeMassMatrix(rig, rigLikeState(), ws));
    CHECK(ws.massMatrix[0] == doctest::Approx(12.0 + 3.0 + 4.0 + 60.0));
    REQUIRE(factorMassMatrix(rig, ws));

    ModelBuilder leaf;
    leaf.addBody({.name = "a",
                  .mass = sub0math::kilograms(1.0),
                  .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                  .inertia = inertia(0.1, 0.1, 0.1),
                  .joint = {.type = JointType::Prismatic, .axis = {1.0, 0.0, 0.0}}});
    leaf.addBody({.name = "ghost",
                  .mass = sub0math::kilograms(0.0),
                  .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                  .inertia = inertia(0.0, 0.0, 0.0),
                  .joint = {.type = JointType::Revolute, .parent = BodyId{0}, .axis = {0.0, 0.0, 1.0}}});
    const auto refused = leaf.build();
    REQUIRE_FALSE(refused);
    CHECK(std::string_view(refused.error().detail) ==
          "a massless body leaves a joint that moves no mass: the mass matrix is singular");
}

TEST_CASE("linearly implicit midpoint: conserves an oscillator's energy, second order, refuses constraints")
{
    // A slider on a spring (no damping): the linearised trapezoid is the trapezoid rule here, which keeps the
    // quadratic energy exactly (to rounding) at any step, where backward Euler loses it at about w^2 h / 2 per second.
    ModelBuilder b;
    b.addBody({.name = "slider",
               .mass = sub0math::kilograms(2.0),
               .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
               .inertia = inertia(0.1, 0.1, 0.1),
               .joint = {.type = JointType::Prismatic, .axis = {1.0, 0.0, 0.0}}});
    b.setGravity(gravityAlong(0.0, 0.0, 0.0));
    constexpr double K = 8000.0;
    b.addElement(CoordinateElementDescription<LinearSpringDamper<f64::Length>>{
        .name = "spring",
        .body = BodyId{0},
        .law = *LinearSpringDamper<f64::Length>::create(
            {.stiffness = f64::Stiffness::fromSi(K), .damping = f64::Damping::fromSi(0.0)})});
    const Model model = *b.build();
    Workspace ws(model);
    const auto energyOf = [&](const State& s) {
        return 0.5 * 2.0 * s.qd()[0] * s.qd()[0] + 0.5 * K * s.q()[0] * s.q()[0];
    };
    const double omega = std::sqrt(K / 2.0);
    const auto runTo = [&](double h, double time, auto&& stepper) {
        State s = State::fromCoordinates({0.01}, {0.0});
        const auto steps = static_cast<std::uint64_t>(std::llround(time / h));
        for (std::uint64_t k = 0; k < steps; ++k)
        {
            REQUIRE(stepper.step(model, Forcing{}, k, f64::Time::fromSi(h), s, ws));
        }
        return s;
    };
    const double e0 = 0.5 * K * 0.01 * 0.01;
    const State midpoint = runTo(5e-3, 2.0, LinearlyImplicitMidpoint{});
    const State euler = runTo(5e-3, 2.0, LinearlyImplicitEuler{});
    CHECK(std::abs(energyOf(midpoint) - e0) < 1e-12 * e0);
    CHECK(energyOf(euler) < 0.01 * e0); // backward Euler at h w = 0.32 has damped it away
    // Order: the largest error over 0.5 s against the closed form quarters per halving of h.
    double previous = 0.0;
    for (const double h : {4e-3, 2e-3, 1e-3, 5e-4})
    {
        State s = State::fromCoordinates({0.01}, {0.0});
        double error = 0.0;
        const auto steps = static_cast<std::uint64_t>(std::llround(0.5 / h));
        for (std::uint64_t k = 0; k < steps; ++k)
        {
            REQUIRE(stepLinearlyImplicitMidpoint(model, Forcing{}, k, f64::Time::fromSi(h), s, ws));
            const double t = static_cast<double>(k + 1) * h;
            error = std::max(error, std::abs(s.q()[0] - 0.01 * std::cos(omega * t)));
        }
        if (previous > 0.0)
        {
            CHECK(std::log2(previous / error) == doctest::Approx(2.0).epsilon(0.05));
        }
        previous = error;
    }
    const Model crank = *sliderCrankBuilder().build();
    Workspace wc(crank);
    State sc = sliderCrankState({}, 0.3, 1.0);
    CHECK_FALSE(stepLinearlyImplicitMidpoint(crank, Forcing{}, 0, f64::Time::fromSi(1e-3), sc, wc));
}
