// Force elements and the linearly implicit Euler step (docs/MULTIBODY.md sections 4 and 5, the M1b gate):
// element Jacobians against finite differences, the tendon geometry against a closed form, stiff stability
// against explicit RK, the energy ledger, and the measured convergence order.

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
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <sub0math/linear/quat.hpp>
#include <sub0math/units.hpp>

#include "sub0physics/elements.hpp"
#include "sub0physics/multibody.hpp"
#include "sub0physics/verification/convergence.hpp"
#include "sub0physics/verification/energy_ledger.hpp"
#include "test_multibody_models.hpp"

using namespace sub0physics;
using namespace sub0physics::multibody;
using namespace testmodels;
namespace f64 = sub0math::f64;

namespace
{

using LinearLaw = LinearSpringDamper<f64::Length>;
using TorsionLaw = LinearSpringDamper<f64::Angle>;

LinearLaw linearLaw(double k, double c)
{
    return *LinearLaw::create({.stiffness = f64::Stiffness::fromSi(k), .damping = f64::Damping::fromSi(c)});
}

TorsionLaw torsionLaw(double k, double c)
{
    return *TorsionLaw::create({.stiffness = TorsionLaw::Spec::StiffnessType::fromSi(k),
                                .damping = TorsionLaw::Spec::DampingType::fromSi(c)});
}

AxialSpringDamper axialLaw(double preload, double k, double compression, double rebound)
{
    return *AxialSpringDamper::create({.preload = f64::Force::fromSi(preload),
                                       .stiffness = f64::Stiffness::fromSi(k),
                                       .compressionDamping = f64::Damping::fromSi(compression),
                                       .reboundDamping = f64::Damping::fromSi(rebound)});
}

BodyDescription spatialBody(std::string_view name, BodyId parent, const Pose& placement, JointType type,
                            sub0math::Vec3d axis, double mass)
{
    return {.name = name,
            .mass = f64::Mass::fromSi(mass),
            .centreOfMass = {metres(0.05), metres(-0.1), metres(0.02)},
            .inertia = diagonalInertia(0.02 * mass, 0.03 * mass, 0.04 * mass),
            .joint = {.type = type, .parent = parent, .placement = placement, .axis = sub0math::normalized(axis)}};
}

Pose turned(double x, double y, double z, sub0math::Vec3d axis, double angle)
{
    return {sub0math::Quatd::fromAxisAngle(sub0math::normalized(axis), f64::Angle::fromSi(angle)),
            {metres(x), metres(y), metres(z)}};
}

// A spatial tree with mixed joints, rotated placements and two branches, so the tendon geometry is
// exercised in 3D: a slider, an arm on it (two revolute joints with skew axes and a prismatic tip), and a
// branch of two revolute joints.
ModelBuilder spatialChainBuilder()
{
    using sub0math::Vec3d;
    ModelBuilder b;
    b.setGravity(gravityAlong(0.3, -9.81, 0.5));
    const BodyId slide = b.addBody(spatialBody("slide", WorldBody, {}, JointType::Prismatic, Vec3d{0.6, 0.8, 0.0}, 2.0));
    const BodyId arm = b.addBody(spatialBody("arm", slide, turned(0.1, 0.2, 0.0, Vec3d{1.0, 0.2, 0.3}, 0.4),
                                             JointType::Revolute, Vec3d{0.0, 0.0, 1.0}, 1.2));
    const BodyId forearm = b.addBody(spatialBody("forearm", arm, turned(0.4, 0.0, 0.1, Vec3d{0.0, 1.0, 0.0}, -0.3),
                                                 JointType::Revolute, Vec3d{0.0, 0.6, 0.8}, 0.9));
    const BodyId tip = b.addBody(spatialBody("tip", forearm, offset(0.3, 0.0, 0.0), JointType::Prismatic,
                                             Vec3d{0.0, 0.0, 1.0}, 0.5));
    const BodyId branch = b.addBody(spatialBody("branch", slide, turned(-0.2, 0.0, 0.0, Vec3d{0.0, 0.0, 1.0}, 0.7),
                                                JointType::Revolute, Vec3d{1.0, 0.0, 0.0}, 0.8));
    const BodyId twig = b.addBody(spatialBody("twig", branch, offset(0.0, -0.3, 0.0), JointType::Revolute,
                                              Vec3d{0.0, 1.0, 0.0}, 0.6));
    b.addSite({.name = "a", .body = tip, .frame = offset(0.1, 0.05, 0.0)});
    b.addSite({.name = "b", .body = twig, .frame = offset(0.0, -0.2, 0.1)});
    b.addSite({.name = "anchor", .body = WorldBody, .frame = offset(0.5, 1.0, -0.3)});
    b.addSite({.name = "c", .body = forearm, .frame = offset(0.2, 0.0, 0.05)});
    return b;
}

State spatialChainState(double scale = 1.0)
{
    return State::fromCoordinates({0.05, 0.4, -0.6, 0.03, 0.3, -0.5},
                                  {0.3 * scale, -0.7 * scale, 0.9 * scale, 0.2 * scale, -0.4 * scale, 0.6 * scale});
}

double siteDistance(const Model& model, const State& s, std::string_view first, std::string_view second)
{
    Workspace ws(model);
    const auto p = computeSiteKinematics(model, s, *model.findSite(first), ws).value().pose.translation;
    const auto q = computeSiteKinematics(model, s, *model.findSite(second), ws).value().pose.translation;
    const double dx = p.x.si() - q.x.si();
    const double dy = p.y.si() - q.y.si();
    const double dz = p.z.si() - q.z.si();
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

struct SpatialOptions
{
    bool coordinateElements = true;
    bool damping = true;
};

// The spatial chain with every law: four tendons (one to a world anchor) and five coordinate elements.
// References are set from the test state so that every law sits inside a branch (not on a kink).
Model spatialChainWithElements(SpatialOptions options = {})
{
    ModelBuilder b = spatialChainBuilder();
    const Model bare = *spatialChainBuilder().build();
    const State s = spatialChainState();
    const double ab = siteDistance(bare, s, "a", "b");
    const double cb = siteDistance(bare, s, "c", "b");
    const double anchorA = siteDistance(bare, s, "anchor", "a");
    const double c = options.damping ? 1.0 : 0.0;
    b.addElement(TendonElementDescription<AxialSpringDamper>{.name = "strut",
                                                             .first = *bare.findSite("a"),
                                                             .second = *bare.findSite("b"),
                                                             .law = axialLaw(60.0, 2500.0, 40.0 * c, 70.0 * c),
                                                             .sense = LawSense::Decreasing,
                                                             .reference = metres(ab + 0.03)});
    b.addElement(TendonElementDescription<LinearLaw>{.name = "cable",
                                                     .first = *bare.findSite("anchor"),
                                                     .second = *bare.findSite("c"),
                                                     .law = linearLaw(800.0, 15.0 * c),
                                                     .reference = metres(0.5)});
    b.addElement(TendonElementDescription<UnilateralSpringDamper>{
        .name = "bump",
        .first = *bare.findSite("c"),
        .second = *bare.findSite("b"),
        .law = *UnilateralSpringDamper::create({.stiffness = f64::Stiffness::fromSi(1.0e4), .damping = f64::Damping::fromSi(50.0 * c)}),
        .sense = LawSense::Decreasing,
        .reference = metres(cb + 0.02)});
    b.addElement(TendonElementDescription<TravelBuffer>{
        .name = "buffer",
        .first = *bare.findSite("anchor"),
        .second = *bare.findSite("a"),
        .law = *TravelBuffer::create({.travel = metres(0.2), .zoneWidth = metres(0.05),
                                      .stiffness = f64::Stiffness::fromSi(2.0e4), .damping = f64::Damping::fromSi(300.0 * c)}),
        .sense = LawSense::Decreasing,
        .reference = metres(anchorA + 0.17)});
    if (options.coordinateElements)
    {
        b.addElement(CoordinateElementDescription<TravelBuffer>{
            .name = "slide buffer",
            .body = *bare.findBody("slide"),
            .law = *TravelBuffer::create({.travel = metres(0.3), .zoneWidth = metres(0.08),
                                          .stiffness = f64::Stiffness::fromSi(5.0e4), .damping = f64::Damping::fromSi(400.0 * c)}),
            .reference = 0.01});
        b.addElement(CoordinateElementDescription<TorsionLaw>{.name = "arm torsion",
                                                              .body = *bare.findBody("arm"),
                                                              .law = torsionLaw(30.0, 0.8 * c),
                                                              .sense = LawSense::Decreasing,
                                                              .reference = 0.3});
        b.addElement(CoordinateElementDescription<AxialSpringDamper>{.name = "tip spring",
                                                                     .body = *bare.findBody("tip"),
                                                                     .law = axialLaw(5.0, 900.0, 12.0 * c, 20.0 * c)});
        b.addElement(CoordinateElementDescription<LinearLaw>{.name = "tip damper",
                                                             .body = *bare.findBody("tip"),
                                                             .law = linearLaw(0.0, 7.0 * c)});
        b.addElement(CoordinateElementDescription<UnilateralSpringDamper>{
            .name = "slide stop",
            .body = *bare.findBody("slide"),
            .law = *UnilateralSpringDamper::create({.stiffness = f64::Stiffness::fromSi(3.0e4), .damping = f64::Damping::fromSi(90.0 * c)}),
            .sense = LawSense::Decreasing,
            .reference = 0.08});
    }
    return *b.build();
}

std::vector<double> elementForceAt(const Model& model, const State& s, Workspace& ws)
{
    REQUIRE(evaluateElements(model, s, ws, ElementJacobians::Skip));
    return ws.elementForce;
}

struct FiniteDifferenceErrors
{
    double stiffness = 0.0; ///< Largest |K - K_fd| / (max |K| + 1).
    double damping = 0.0;   ///< Largest |D - D_fd| / (max |D| + 1).
};

// Compares elementStiffness and elementDamping with central differences of elementForce.
FiniteDifferenceErrors jacobianErrors(const Model& model, const State& s)
{
    const std::size_t n = model.dofCount();
    Workspace ws(model);
    REQUIRE(evaluateElements(model, s, ws, ElementJacobians::Compute));
    const std::vector<double> k = ws.elementStiffness;
    const std::vector<double> d = ws.elementDamping;
    const auto scale = [](const std::vector<double>& m) {
        double largest = 0.0;
        for (const double v : m)
        {
            largest = std::max(largest, std::abs(v));
        }
        return largest + 1.0;
    };
    FiniteDifferenceErrors errors;
    const double step = 1e-6;
    for (std::size_t j = 0; j < n; ++j)
    {
        for (const bool rate : {false, true})
        {
            State plus = s;
            State minus = s;
            (rate ? plus.qd() : plus.q())[j] += step;
            (rate ? minus.qd() : minus.q())[j] -= step;
            const std::vector<double> up = elementForceAt(model, plus, ws);
            const std::vector<double> down = elementForceAt(model, minus, ws);
            for (std::size_t i = 0; i < n; ++i)
            {
                const double fd = (up[i] - down[i]) / (2.0 * step);
                if (rate)
                {
                    errors.damping = std::max(errors.damping, std::abs(d[i * n + j] - fd) / scale(d));
                }
                else
                {
                    errors.stiffness = std::max(errors.stiffness, std::abs(k[i * n + j] - fd) / scale(k));
                }
            }
        }
    }
    return errors;
}

double balanceOf(const Model& model, const State& s, Workspace& ws)
{
    return computeEnergy(model, s, ws).value().balance().si();
}

// A slider of mass m along x with a coordinate spring-damper (k, c), no gravity along the slide.
Model stiffSlider(double mass, double k, double c)
{
    ModelBuilder b;
    b.setGravity(gravityAlong(0.0, -9.81, 0.0));
    const BodyId slider = b.addBody({.name = "slider",
                                     .mass = f64::Mass::fromSi(mass),
                                     .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                                     .inertia = diagonalInertia(0.01, 0.01, 0.01),
                                     .joint = {.type = JointType::Prismatic, .axis = {1.0, 0.0, 0.0}}});
    b.addElement(CoordinateElementDescription<LinearLaw>{.name = "spring", .body = slider, .law = linearLaw(k, c)});
    return *b.build();
}

// A pendulum about z (bob at radius r below the pivot) with a tendon from the bob to a world anchor at
// (0, -a, 0): l(theta) = sqrt(a^2 + r^2 - 2 a r cos theta).
template<typename L>
Model anchoredPendulum(double r, double a, const L& law, LawSense sense, double reference)
{
    ModelBuilder b;
    b.setGravity(gravityAlong(0.0, -9.81, 0.0));
    const BodyId bob = b.addBody({.name = "bob",
                                  .mass = f64::Mass::fromSi(1.0),
                                  .centreOfMass = {metres(0.0), metres(-r), metres(0.0)},
                                  .inertia = diagonalInertia(0.01, 0.01, 0.01),
                                  .joint = {.type = JointType::Revolute, .axis = {0.0, 0.0, 1.0}}});
    const SiteId end = b.addSite({.name = "bob", .body = bob, .frame = offset(0.0, -r, 0.0)});
    const SiteId anchor = b.addSite({.name = "anchor", .body = WorldBody, .frame = offset(0.0, -a, 0.0)});
    b.addElement(TendonElementDescription<L>{.name = "tendon", .first = end, .second = anchor, .law = law,
                                             .sense = sense, .reference = metres(reference)});
    return *b.build();
}

} // namespace

// ── Description and builder ───────────────────────────────────────────────────

TEST_CASE("Model builder compiles and validates force elements")
{
    const Model model = spatialChainWithElements();
    CHECK(model.elementCount() == 9U);
    CHECK(model.tendonCount() == 4U);
    CHECK_FALSE(model.elementsFollowTree()); // the strut joins two branches
    CHECK(model.findElement("strut") == ElementId{0});
    CHECK(model.findElement("slide stop") == ElementId{8});
    CHECK_FALSE(model.findElement("nothing"));
    CHECK(model.siteBody(*model.findSite("anchor")) == WorldBody);
    CHECK(planarTreeWithElements().elementCount() == 2U);
    CHECK(anchoredPendulum(1.0, 2.0, linearLaw(1.0, 0.0), LawSense::Increasing, 1.0).elementsFollowTree());

    const auto faultOf = [](ModelBuilder& b) { return std::string_view(b.build().error().detail); };
    {
        ModelBuilder b = spatialChainBuilder(); // a law on an angle on a prismatic joint
        b.addElement(CoordinateElementDescription<TorsionLaw>{.name = "e", .body = BodyId{0}, .law = torsionLaw(1.0, 0.0)});
        CHECK(faultOf(b).starts_with("coordinate element law does not suit the joint"));
    }
    {
        ModelBuilder b = spatialChainBuilder(); // a law on a length on a revolute joint
        b.addElement(CoordinateElementDescription<LinearLaw>{.name = "e", .body = BodyId{1}, .law = linearLaw(1.0, 0.0)});
        CHECK(faultOf(b).starts_with("coordinate element law does not suit the joint"));
    }
    {
        ModelBuilder b = spatialChainBuilder();
        b.addElement(CoordinateElementDescription<LinearLaw>{.name = "e", .body = BodyId{17}, .law = linearLaw(1.0, 0.0)});
        CHECK(faultOf(b) == "coordinate element body is not a body of this model");
    }
    {
        ModelBuilder b = spatialChainBuilder();
        b.addElement(CoordinateElementDescription<LinearLaw>{
            .name = "e", .body = BodyId{0}, .law = linearLaw(1.0, 0.0), .reference = std::numeric_limits<double>::infinity()});
        CHECK(faultOf(b) == "coordinate element reference must be finite");
    }
    {
        ModelBuilder b = spatialChainBuilder();
        b.addElement(CoordinateElementDescription<LinearLaw>{.name = "e", .body = BodyId{0}, .law = linearLaw(1.0, 0.0)});
        b.addElement(CoordinateElementDescription<LinearLaw>{.name = "e", .body = BodyId{3}, .law = linearLaw(1.0, 0.0)});
        CHECK(faultOf(b) == "element name is not unique");
    }
    {
        ModelBuilder b = spatialChainBuilder();
        b.addElement(CoordinateElementDescription<LinearLaw>{.name = "", .body = BodyId{0}, .law = linearLaw(1.0, 0.0)});
        CHECK(faultOf(b) == "element name is empty");
    }
    {
        ModelBuilder b = spatialChainBuilder();
        b.addElement(TendonElementDescription<LinearLaw>{.name = "t", .first = SiteId{0}, .second = SiteId{9}, .law = linearLaw(1.0, 0.0)});
        CHECK(faultOf(b) == "tendon site is not a site of this model");
    }
    {
        ModelBuilder b = spatialChainBuilder();
        b.addSite({.name = "a2", .body = BodyId{3}});
        b.addElement(TendonElementDescription<LinearLaw>{.name = "t", .first = SiteId{0}, .second = SiteId{4}, .law = linearLaw(1.0, 0.0)});
        CHECK(faultOf(b) == "tendon sites must be on different bodies (or one on the world)");
    }
    {
        ModelBuilder b = spatialChainBuilder();
        b.addElement(TendonElementDescription<LinearLaw>{
            .name = "t", .first = SiteId{0}, .second = SiteId{1}, .law = linearLaw(1.0, 0.0), .reference = metres(-0.1)});
        CHECK(faultOf(b) == "tendon reference length must be finite and not negative");
    }
}

TEST_CASE("Linear spring-damper law on an angle: exact response and power balance")
{
    const TorsionLaw law = torsionLaw(40.0, 2.5);
    const LawResponse<f64::Angle> r = law.respond(f64::Angle::fromSi(0.2), f64::AngularVelocity::fromSi(-1.5));
    CHECK(r.force.si() == doctest::Approx(-40.0 * 0.2 + 2.5 * 1.5));
    CHECK(r.dForceDPosition.si() == -40.0);
    CHECK(r.dForceDVelocity.si() == -2.5);
    CHECK(r.potential.si() == doctest::Approx(0.5 * 40.0 * 0.04));
    CHECK(r.dissipationPower.si() == doctest::Approx(2.5 * 2.25));
    // force * qd = -dV/dt - P, with dV/dt = k q qd.
    CHECK(r.force.si() * -1.5 == doctest::Approx(-(40.0 * 0.2 * -1.5) - r.dissipationPower.si()));
    CHECK_FALSE(TorsionLaw::create({.stiffness = TorsionLaw::Spec::StiffnessType::fromSi(-1.0),
                                    .damping = TorsionLaw::Spec::DampingType::fromSi(0.0)}));
}

// ── Element evaluation ────────────────────────────────────────────────────────

TEST_CASE("Every element's Jacobians match finite differences")
{
    const Model model = spatialChainWithElements();
    double worstStiffness = 0.0;
    double worstDamping = 0.0;
    for (const double scale : {1.0, -0.6, 2.0})
    {
        const FiniteDifferenceErrors errors = jacobianErrors(model, spatialChainState(scale));
        worstStiffness = std::max(worstStiffness, errors.stiffness);
        worstDamping = std::max(worstDamping, errors.damping);
    }
    const FiniteDifferenceErrors tree = jacobianErrors(planarTreeWithElements(), [] {
        State s = State::atRest(planarTree());
        for (std::size_t i = 0; i < PlanarTreeDof; ++i)
        {
            s.q()[i] = 0.2 * std::sin(1.3 * static_cast<double>(i) + 0.4);
            s.qd()[i] = 0.7 * std::cos(0.9 * static_cast<double>(i));
        }
        return s;
    }());
    MESSAGE("element Jacobians vs central differences (h = 1e-6), relative to the largest entry: spatial chain dtau/dq "
            << worstStiffness << ", dtau/dqd " << worstDamping << "; planar tree " << tree.stiffness << ", " << tree.damping);
    CHECK(worstStiffness < 1e-7);
    CHECK(worstDamping < 1e-7);
    CHECK(tree.stiffness < 1e-7);
    CHECK(tree.damping < 1e-7);
}

TEST_CASE("Element forces are minus the gradient of their potential at rest rates")
{
    // With zero rates every damper is idle, so tau_e = -dV/dq: checks J_l and the signs independently of K.
    const Model model = spatialChainWithElements();
    State s = spatialChainState(0.0);
    Workspace ws(model);
    const std::vector<double> tau = elementForceAt(model, s, ws);
    const double step = 1e-6;
    for (std::size_t j = 0; j < model.dofCount(); ++j)
    {
        State plus = s;
        State minus = s;
        plus.q()[j] += step;
        minus.q()[j] -= step;
        REQUIRE(evaluateElements(model, plus, ws, ElementJacobians::Skip));
        const double up = ws.elementPotential;
        REQUIRE(evaluateElements(model, minus, ws, ElementJacobians::Skip));
        const double down = ws.elementPotential;
        CHECK(tau[j] == doctest::Approx(-(up - down) / (2.0 * step)).epsilon(1e-6).scale(1.0));
    }
}

TEST_CASE("Tendon geometric stiffness matches the closed form")
{
    // A pendulum (r = 0.8) with a linear tendon to a world anchor at a = 1.5 below the pivot.
    const double r = 0.8;
    const double a = 1.5;
    const double k = 700.0;
    const double reference = 0.9;
    const Model model = anchoredPendulum(r, a, linearLaw(k, 0.0), LawSense::Increasing, reference);
    Workspace ws(model);
    for (const double theta : {-2.5, -0.7, 0.3, 1.1, 2.9})
    {
        const State s = State::fromCoordinates({theta}, {0.4});
        REQUIRE(evaluateElements(model, s, ws));
        const double l = std::sqrt(a * a + r * r - 2.0 * a * r * std::cos(theta));
        const double dl = a * r * std::sin(theta) / l;
        const double d2l = a * r * std::cos(theta) / l - dl * dl / l;
        const double force = -k * (l - reference); // conjugate to l
        CHECK(ws.elementForce[0] == doctest::Approx(force * dl).epsilon(1e-12));
        CHECK(ws.elementStiffness[0] == doctest::Approx(-k * dl * dl + force * d2l).epsilon(1e-12));
        CHECK(ws.elementPotential == doctest::Approx(0.5 * k * (l - reference) * (l - reference)).epsilon(1e-12));
    }

    // A pure preload (no rate, no damping): the whole stiffness is the geometric term F d2l/dq2, which an
    // implementation without it would report as zero.
    const Model preloaded = spatialChainWithElements({.coordinateElements = false, .damping = false});
    ModelBuilder strutOnly = spatialChainBuilder();
    strutOnly.addElement(TendonElementDescription<AxialSpringDamper>{.name = "preload",
                                                                     .first = SiteId{0},
                                                                     .second = SiteId{1},
                                                                     .law = axialLaw(500.0, 0.0, 0.0, 0.0),
                                                                     .sense = LawSense::Decreasing});
    const Model pure = *strutOnly.build();
    const State s = spatialChainState();
    Workspace pw(pure);
    REQUIRE(evaluateElements(pure, s, pw));
    double largest = 0.0;
    for (const double v : pw.elementStiffness)
    {
        largest = std::max(largest, std::abs(v));
    }
    MESSAGE("pure-preload tendon: largest geometric stiffness entry " << largest << " N m/rad");
    CHECK(largest > 10.0);
    CHECK(jacobianErrors(pure, s).stiffness < 1e-8);
    CHECK(jacobianErrors(preloaded, s).stiffness < 1e-8);
}

TEST_CASE("Element site wrenches reproduce the generalized force")
{
    // tau_e = sum over tendon ends of J_site^T f: the loads reported at the sites are the forces that moved
    // the bodies.
    const Model model = spatialChainWithElements({.coordinateElements = false});
    const std::size_t n = model.dofCount();
    const State s = spatialChainState();
    Workspace ws(model);
    REQUIRE(evaluateElements(model, s, ws));
    const std::vector<double> tau = ws.elementForce;
    const std::vector<ElementSiteWrench> wrenches = ws.elementSiteWrenches;
    REQUIRE(wrenches.size() == 8U);
    std::vector<double> sum(n, 0.0);
    std::vector<double> jacobian(3 * n);
    for (const ElementSiteWrench& w : wrenches)
    {
        CHECK(w.element.valid());
        CHECK(w.moment.x == 0.0);
        REQUIRE(computeSiteJacobian(model, s, w.site, ws, jacobian));
        for (std::size_t c = 0; c < n; ++c)
        {
            sum[c] += jacobian[c] * w.force.x + jacobian[n + c] * w.force.y + jacobian[2 * n + c] * w.force.z;
        }
    }
    for (std::size_t c = 0; c < n; ++c)
    {
        CHECK(sum[c] == doctest::Approx(tau[c]).epsilon(1e-12).scale(1.0));
    }
    // The two ends of a tendon carry equal and opposite forces.
    CHECK(wrenches[0].force.x == doctest::Approx(-wrenches[1].force.x));
}

TEST_CASE("Forward and inverse dynamics and energy include the elements")
{
    const Model model = spatialChainWithElements();
    const Model bare = *spatialChainBuilder().build();
    const State s = spatialChainState();
    Workspace ws(model);
    Workspace bw(bare);
    const std::vector<double> applied{0.3, -0.2, 0.1, 0.4, 0.0, -0.1};
    REQUIRE(forwardDynamics(model, s, applied, ws));
    std::vector<double> total = applied;
    for (std::size_t i = 0; i < total.size(); ++i)
    {
        total[i] += ws.elementForce[i];
    }
    REQUIRE(forwardDynamics(bare, s, total, bw));
    std::vector<double> back(applied.size());
    REQUIRE(inverseDynamics(model, s, ws.qdd, back, ws));
    for (std::size_t i = 0; i < applied.size(); ++i)
    {
        CHECK(ws.qdd[i] == doctest::Approx(bw.qdd[i]).epsilon(1e-12).scale(1.0));
        CHECK(back[i] == doctest::Approx(applied[i]).epsilon(1e-9).scale(1.0));
    }
    const EnergyAccount withElements = computeEnergy(model, s, ws).value();
    const EnergyAccount without = computeEnergy(bare, s, bw).value();
    REQUIRE(evaluateElements(model, s, ws));
    CHECK(withElements.potential.si() == doctest::Approx(without.potential.si() + ws.elementPotential).epsilon(1e-14));
    CHECK(withElements.kinetic == without.kinetic);
}

TEST_CASE("A tendon whose sites coincide is outside its domain, and steps refuse it")
{
    ModelBuilder b;
    b.setGravity(gravityAlong(0.0, 0.0, 0.0));
    const BodyId slider = b.addBody({.name = "slider",
                                     .mass = f64::Mass::fromSi(1.0),
                                     .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                                     .inertia = diagonalInertia(0.01, 0.01, 0.01),
                                     .joint = {.type = JointType::Prismatic, .axis = {1.0, 0.0, 0.0}}});
    const SiteId end = b.addSite({.name = "end", .body = slider});
    const SiteId anchor = b.addSite({.name = "anchor", .body = WorldBody, .frame = offset(0.5, 0.0, 0.0)});
    b.addElement(TendonElementDescription<LinearLaw>{.name = "t", .first = end, .second = anchor, .law = linearLaw(10.0, 0.0)});
    const Model model = *b.build();
    Workspace ws(model);
    State s = State::fromCoordinates({0.5}, {0.0});
    const State before = s;
    const auto refused = stepLinearlyImplicitEuler(model, {}, 0U, sub0math::seconds(1e-3), s, ws);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().kind == FaultKind::OutsideDomain);
    CHECK(s == before);
    CHECK_FALSE(stepSemiImplicitEuler(model, {}, 0U, sub0math::seconds(1e-3), s, ws));
    CHECK(s == before);
}

// ── The linearly implicit Euler step ──────────────────────────────────────────

TEST_CASE("Linearly implicit Euler is a transactional Stepper")
{
    static_assert(Stepper<LinearlyImplicitEuler>);
    const Model model = planarTreeWithElements();
    Workspace ws(model);
    LinearlyImplicitEuler stepper;
    State s = State::atRest(model);
    s.q()[2] = 0.3;
    REQUIRE(stepper.step(model, {}, 0U, sub0math::seconds(1e-3), s, ws));
    CHECK(s.tick == 1U);
    const State before = s;
    const auto rejected = [&](const std::expected<void, Fault>& result, FaultKind kind) {
        REQUIRE_FALSE(result);
        CHECK(result.error().kind == kind);
        CHECK(s == before);
    };
    rejected(stepper.step(model, {}, 0U, sub0math::seconds(1e-3), s, ws), FaultKind::InvalidInput); // wrong tick
    rejected(stepper.step(model, {}, 1U, sub0math::seconds(0.0), s, ws), FaultKind::InvalidInput);
    const std::array<double, 3> wrongSize{};
    rejected(stepper.step(model, {wrongSize}, 1U, sub0math::seconds(1e-3), s, ws), FaultKind::InvalidInput);
    const Model other = planarTreeWithElements();
    rejected(stepper.step(other, {}, 1U, sub0math::seconds(1e-3), s, ws), FaultKind::InvalidInput);
    s.qd()[4] = std::numeric_limits<double>::infinity();
    const State poisoned = s;
    REQUIRE_FALSE(stepper.step(model, {}, 1U, sub0math::seconds(1e-3), s, ws));
    CHECK(s == poisoned);
}

TEST_CASE("Linearly implicit Euler follows the exact discrete map of a linear oscillator")
{
    // Slider m with spring k and damper c, and a constant applied force f: the step is linear, so it is
    // exactly backward Euler, (m + h c + h^2 k) dv = h (f - k x - c v - h k v).
    const double m = 2.0;
    const double k = 5.0e4;
    const double c = 30.0;
    const double f = 3.0;
    const double h = 2e-3;
    const Model model = stiffSlider(m, k, c);
    Workspace ws(model);
    State s = State::fromCoordinates({0.02}, {-0.5});
    double x = 0.02;
    double v = -0.5;
    const std::array<double, 1> tau{f};
    for (std::uint64_t tick = 0; tick < 200U; ++tick)
    {
        REQUIRE(stepLinearlyImplicitEuler(model, {tau}, tick, sub0math::seconds(h), s, ws));
        const double dv = h * (f - k * x - c * v - h * k * v) / (m + h * c + h * h * k);
        v += dv;
        x += h * v;
        CHECK(s.q()[0] == doctest::Approx(x).epsilon(1e-12).scale(1.0));
        CHECK(s.qd()[0] == doctest::Approx(v).epsilon(1e-12).scale(1.0));
    }
    CHECK_FALSE(ws.geometricStiffnessDropped);
}

TEST_CASE("Linearly implicit Euler solves the assembled system on both factorisation paths")
{
    // Independent check of the dense LDL^T (branches coupled) and the tree LTDL (elements on one chain):
    // assemble M - h D - h^2 K~ and the right-hand side from the workspace and solve by Gaussian elimination.
    const auto check = [](const Model& model, const State& start, double h) {
        const std::size_t n = model.dofCount();
        Workspace ws(model);
        State s = start;
        std::vector<double> tau(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            tau[i] = 0.05 * static_cast<double>(i + 1);
        }
        REQUIRE(stepLinearlyImplicitEuler(model, {tau}, 0U, sub0math::seconds(h), s, ws));
        REQUIRE_FALSE(ws.geometricStiffnessDropped);
        Workspace fresh(model);
        REQUIRE(computeMassMatrix(model, start, fresh));
        REQUIRE(computeBiasForces(model, start, fresh));
        REQUIRE(evaluateElements(model, start, fresh));
        std::vector<double> a(n * (n + 1));
        for (std::size_t i = 0; i < n; ++i)
        {
            double stiffnessRate = 0.0;
            for (std::size_t j = 0; j < n; ++j)
            {
                const double k = 0.5 * (fresh.elementStiffness[i * n + j] + fresh.elementStiffness[j * n + i]);
                a[i * (n + 1) + j] = fresh.massMatrix[i * n + j] - h * fresh.elementDamping[i * n + j] - h * h * k;
                stiffnessRate += k * start.qd()[j];
            }
            a[i * (n + 1) + n] = h * (tau[i] + fresh.elementForce[i] - fresh.bias[i] + h * stiffnessRate);
        }
        for (std::size_t col = 0; col < n; ++col)
        {
            std::size_t pivot = col;
            for (std::size_t r = col + 1; r < n; ++r)
            {
                if (std::abs(a[r * (n + 1) + col]) > std::abs(a[pivot * (n + 1) + col]))
                {
                    pivot = r;
                }
            }
            for (std::size_t c = 0; c <= n; ++c)
            {
                std::swap(a[col * (n + 1) + c], a[pivot * (n + 1) + c]);
            }
            for (std::size_t r = col + 1; r < n; ++r)
            {
                const double factor = a[r * (n + 1) + col] / a[col * (n + 1) + col];
                for (std::size_t c = col; c <= n; ++c)
                {
                    a[r * (n + 1) + c] -= factor * a[col * (n + 1) + c];
                }
            }
        }
        std::vector<double> dqd(n);
        for (std::size_t r = n; r-- > 0;)
        {
            double sum = a[r * (n + 1) + n];
            for (std::size_t c = r + 1; c < n; ++c)
            {
                sum -= a[r * (n + 1) + c] * dqd[c];
            }
            dqd[r] = sum / a[r * (n + 1) + r];
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            CHECK(s.qd()[i] == doctest::Approx(start.qd()[i] + dqd[i]).epsilon(1e-11).scale(1.0));
            CHECK(s.q()[i] == doctest::Approx(start.q()[i] + h * (start.qd()[i] + dqd[i])).epsilon(1e-11).scale(1.0));
        }
    };
    const Model coupled = planarTreeWithElements({.strutStiffness = 2.0e5});
    REQUIRE_FALSE(coupled.elementsFollowTree());
    State start = State::atRest(coupled);
    for (std::size_t i = 0; i < PlanarTreeDof; ++i)
    {
        start.q()[i] = 0.1 * std::sin(static_cast<double>(i));
        start.qd()[i] = 0.3 * std::cos(static_cast<double>(i));
    }
    check(coupled, start, 4e-3);
    const Model chain = spatialChainWithElements(); // dense too; skew axes
    check(chain, spatialChainState(), 2e-3);
    const Model onePath = anchoredPendulum(0.8, 1.5, linearLaw(2.0e4, 40.0), LawSense::Increasing, 0.9);
    REQUIRE(onePath.elementsFollowTree());
    check(onePath, State::fromCoordinates({0.4}, {1.0}), 5e-3);
}

TEST_CASE("Stiff spring-damper: implicit Euler is stable at a step where explicit RK4 blows up")
{
    // m = 1 kg on k = 1e6 N/m (omega = 1000 rad/s, period 6.3 ms), c = 10 N s/m: RK4's stability limit is
    // h omega ~ 2.8, so h = 10 ms (h omega = 10) is far outside it.
    const Model model = stiffSlider(1.0, 1.0e6, 10.0);
    const double h = 1e-2;
    const State start = State::fromCoordinates({0.01}, {0.0});
    Workspace ws(model);
    const double initial = balanceOf(model, start, ws);

    Rk4Stepper rk4(model, ws);
    State explicitState = start;
    bool blewUp = false;
    double explicitPeak = 0.0;
    for (std::uint64_t tick = 0; tick < 100U && !blewUp; ++tick)
    {
        if (!rk4.step(model, {}, tick, sub0math::seconds(h), explicitState, ws))
        {
            blewUp = true;
            break;
        }
        explicitPeak = std::max(explicitPeak, computeEnergy(model, explicitState, ws).value().mechanical().si());
        blewUp = !(explicitPeak < 1e6 * initial);
    }

    State s = start;
    double peakBalance = initial;
    for (std::uint64_t tick = 0; tick < 100U; ++tick)
    {
        REQUIRE(stepLinearlyImplicitEuler(model, {}, tick, sub0math::seconds(h), s, ws));
        peakBalance = std::max(peakBalance, balanceOf(model, s, ws));
    }
    const EnergyAccount end = computeEnergy(model, s, ws).value();
    MESSAGE("stiff slider, h = 10 ms: RK4 " << std::string(blewUp ? "blew up" : "stayed bounded") << " (peak energy "
            << explicitPeak << " J from " << initial << " J); implicit Euler: energy " << end.mechanical().si()
            << " J, dissipated " << end.dissipated.si() << " J, peak balance " << peakBalance << " J");
    CHECK(blewUp);
    // Starting at rest, the left-rectangle dissipation telescopes to at most the true discrete dissipation,
    // so the balance never rises above its start (backward Euler's own damping only lowers it).
    CHECK(peakBalance <= initial * (1.0 + 1e-12));
    CHECK(end.mechanical().si() < 1e-3 * initial);
    CHECK(std::abs(s.q()[0]) < 1e-4);
}

TEST_CASE("Stiff coupled tree: implicit Euler is stable where RK4 blows up")
{
    // The planar tree with a 1e6 N/m strut between branches and a 1e6 N/m cart spring, stepped at 5 ms.
    const Model model = planarTreeWithElements({.cartStiffness = 1.0e6, .cartDamping = 50.0, .strutStiffness = 1.0e6,
                                                .strutDamping = 200.0});
    State start = State::atRest(model);
    start.q()[0] = 0.02; // the cart spring loaded (reference 0.05 m)
    start.q()[4] = 0.3;  // the middle branch, which the strut does not touch, swung
    start.qd()[1] = 0.5;
    const double h = 5e-3;
    Workspace ws(model);
    const double initial = computeEnergy(model, start, ws).value().mechanical().si();

    Rk4Stepper rk4(model, ws);
    State explicitState = start;
    bool blewUp = false;
    for (std::uint64_t tick = 0; tick < 200U && !blewUp; ++tick)
    {
        blewUp = !rk4.step(model, {}, tick, sub0math::seconds(h), explicitState, ws) ||
                 !(std::abs(computeEnergy(model, explicitState, ws).value().mechanical().si()) < 1e6 * std::abs(initial) + 1e3);
    }

    State s = start;
    const double startBalance = balanceOf(model, s, ws);
    double drift = 0.0;
    int fallbacks = 0;
    for (std::uint64_t tick = 0; tick < 400U; ++tick)
    {
        const auto stepped = stepLinearlyImplicitEuler(model, {}, tick, sub0math::seconds(h), s, ws);
        if (!stepped)
        {
            MESSAGE("tick " << tick << ": " << std::string(stepped.error().detail) << " dropped " << ws.geometricStiffnessDropped);
        }
        REQUIRE(stepped);
        fallbacks += ws.geometricStiffnessDropped ? 1 : 0;
        drift = std::max(drift, balanceOf(model, s, ws) - startBalance);
    }
    MESSAGE("stiff tree, h = 5 ms: RK4 " << std::string(blewUp ? "blew up" : "stayed bounded") << "; implicit Euler balance rose at most "
            << drift << " J from " << startBalance << " J, dissipated " << s.dissipated.si() << " J; geometric-stiffness fallbacks "
            << fallbacks << " of 400 steps");
    CHECK(blewUp);
    CHECK(drift < 0.05 * std::abs(startBalance) + 1e-6);
}

TEST_CASE("An indefinite stiffness falls back to the damping-only matrix")
{
    // A strut pushing a pendulum's bob away from an anchor beyond it: at theta = 0 the length is minimal, so a
    // pushing force is destabilising (positive dtau/dq = F a r / l). With F = 1e5 N and h = 10 ms,
    // M - h^2 K is indefinite.
    const Model model = anchoredPendulum(1.0, 3.0, axialLaw(1.0e5, 0.0, 0.0, 0.0), LawSense::Decreasing, 2.0);
    Workspace ws(model);
    State s = State::fromCoordinates({0.01}, {0.0});
    REQUIRE(stepLinearlyImplicitEuler(model, {}, 0U, sub0math::seconds(1e-2), s, ws));
    CHECK(ws.geometricStiffnessDropped);
    CHECK(s.qd()[0] > 0.0); // pushed further from the unstable equilibrium
    State gentle = State::fromCoordinates({0.01}, {0.0});
    REQUIRE(stepLinearlyImplicitEuler(model, {}, 0U, sub0math::seconds(1e-4), gentle, ws));
    CHECK_FALSE(ws.geometricStiffnessDropped);
}

TEST_CASE("Energy ledger closes over a run with elements and applied forces")
{
    // The planar tree with both elements, an applied force on the cart and torques on two links, for 1 s.
    const Model model = planarTreeWithElements();
    State start = State::atRest(model);
    for (std::size_t i = 0; i < PlanarTreeDof; ++i)
    {
        start.q()[i] = 0.25 * std::sin(1.7 * static_cast<double>(i) + 0.3);
    }
    std::vector<double> tau(PlanarTreeDof, 0.0);
    tau[0] = 4.0;
    tau[2] = 0.3;
    tau[8] = -0.2;
    const auto run = [&](auto stepper, double h) {
        Workspace ws(model);
        State s = start;
        const EnergyAccount initial = computeEnergy(model, s, ws).value();
        EnergyLedger ledger(initial, initial.mechanical());
        const auto steps = static_cast<std::uint64_t>(std::llround(1.0 / h));
        for (std::uint64_t tick = 0; tick < steps; ++tick)
        {
            REQUIRE(stepper.step(model, {tau}, tick, sub0math::seconds(h), s, ws));
            const EnergyAccount account = computeEnergy(model, s, ws).value();
            ledger.includeScale(account.kinetic);
            ledger.includeScale(account.potential);
            ledger.includeScale(account.dissipated);
            ledger.includeScale(account.externalWork);
            ledger.record(account);
        }
        CHECK(s.dissipated.si() > 0.0);
        CHECK(s.externalWork.si() != 0.0);
        return ledger.relativeResidual();
    };
    Workspace rkWorkspace(model);
    const double rk4 = run(Rk4Stepper(model, rkWorkspace), 1e-3);
    std::array<double, 3> implicit{};
    std::array<double, 3> semi{};
    for (std::size_t k = 0; k < 3; ++k)
    {
        const double h = 2e-3 / static_cast<double>(1U << k);
        implicit[k] = run(LinearlyImplicitEuler{}, h);
        semi[k] = run(SemiImplicitEuler{}, h);
    }
    MESSAGE("energy ledger relative residual over 1 s: RK4 (h = 1 ms) " << rk4 << "; implicit Euler (h = 2, 1, 0.5 ms) "
            << implicit[0] << ", " << implicit[1] << ", " << implicit[2] << " (order "
            << observedOrder(implicit[1], implicit[2], 2.0) << "); semi-implicit Euler " << semi[0] << ", " << semi[1]
            << ", " << semi[2]);
    CHECK(rk4 < 1e-6);
    CHECK(implicit[2] < 0.05);
    CHECK(observedOrder(implicit[1], implicit[2], 2.0) == doctest::Approx(1.0).epsilon(0.2));
    CHECK(semi[2] < 0.05);
}

TEST_CASE("Linearly implicit Euler converges at first order")
{
    const Model model = planarTreeWithElements();
    State start = State::atRest(model);
    for (std::size_t i = 0; i < PlanarTreeDof; ++i)
    {
        start.q()[i] = 0.3 * std::sin(1.1 * static_cast<double>(i) + 0.5);
        start.qd()[i] = 0.4 * std::cos(0.8 * static_cast<double>(i));
    }
    const double horizon = 0.5;
    const auto run = [&](auto stepper, double h) {
        Workspace ws(model);
        State s = start;
        const auto steps = static_cast<std::uint64_t>(std::llround(horizon / h));
        for (std::uint64_t tick = 0; tick < steps; ++tick)
        {
            REQUIRE(stepper.step(model, {}, tick, sub0math::seconds(h), s, ws));
        }
        return s;
    };
    Workspace refWs(model);
    const State reference = run(Rk4Stepper(model, refWs), 2e-5);
    const auto error = [&](const State& s) {
        double sum = 0.0;
        for (std::size_t i = 0; i < s.values.size(); ++i)
        {
            sum += (s.values[i] - reference.values[i]) * (s.values[i] - reference.values[i]);
        }
        return std::sqrt(sum);
    };
    std::array<double, 4> errors{};
    for (std::size_t k = 0; k < errors.size(); ++k)
    {
        errors[k] = error(run(LinearlyImplicitEuler{}, 2e-3 / static_cast<double>(1U << k)));
    }
    MESSAGE("implicit Euler error at t = 0.5 s for h = 2, 1, 0.5, 0.25 ms: " << errors[0] << ", " << errors[1] << ", "
            << errors[2] << ", " << errors[3] << "; observed orders " << observedOrder(errors[0], errors[1], 2.0) << ", "
            << observedOrder(errors[1], errors[2], 2.0) << ", " << observedOrder(errors[2], errors[3], 2.0));
    CHECK(observedOrder(errors[1], errors[2], 2.0) == doctest::Approx(1.0).epsilon(0.1));
    CHECK(observedOrder(errors[2], errors[3], 2.0) == doctest::Approx(1.0).epsilon(0.1));
}
