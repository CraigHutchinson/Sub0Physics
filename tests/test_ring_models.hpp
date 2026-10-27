#pragma once

// Models with deformable rings (docs/WHEEL.md section 2.3), shared by the ring tests, the allocation gate and the
// benchmark: a small planar tree in the xz plane (a frame on two slides and a pitch hinge, carrying two hubs), each hub
// with a ring, built either with ring substructures (condensed) or with the same description expanded into ordinary
// bodies, coordinate elements and tendons (the oracle).

#include <array>
#include <cstddef>
#include <random>
#include <string>
#include <vector>

#include <sub0math/units.hpp>

#include "sub0physics/multibody.hpp"
#include "test_multibody_models.hpp"

namespace ringmodels
{

using namespace sub0physics::multibody;
using testmodels::diagonalInertia;
using testmodels::gravityAlong;
using testmodels::metres;
using testmodels::offset;

/** The tree's coordinates: travel, heave, pitch, rear spin, front spin. */
inline constexpr std::size_t TreeDof = 5;

/** Ring parameters (the reference's defaults, docs/WHEEL.md section 1). */
struct RingParameters
{
    std::size_t nodes = 16;
    double radius = 0.352;      ///< Centreline radius [m].
    double skin = 0.008;        ///< [m].
    double mass = 1.0;          ///< Ring mass [kg].
    double foundation = 3e5;    ///< [N/m^2].
    double extension = 500.0;   ///< EA [N].
    double viscosity = 0.5;     ///< [N s].
    double gravity = 9.80665;   ///< Along -z [m/s^2].
};

inline UniformRingSpec ringSpec(std::string_view name, BodyId host, const RingParameters& p)
{
    using namespace sub0math;
    return {.name = name,
            .host = host,
            .nodes = p.nodes,
            .centrelineRadius = metres(p.radius),
            .skinRadius = metres(p.skin),
            .ringMass = kilograms(p.mass),
            .foundationStiffness = f64::Pressure::fromSi(p.foundation),
            .extensionStiffness = f64::Force::fromSi(p.extension),
            .edgeViscosity = f64::Impulse::fromSi(p.viscosity)};
}

/** Adds the tree; returns the two hubs (rear, front). */
inline std::array<BodyId, 2> addTree(ModelBuilder& builder, double gravity)
{
    using sub0math::kilograms;
    builder.setGravity(gravityAlong(0.0, 0.0, -gravity));
    const BodyDescription massless{};
    BodyDescription travel = massless;
    travel.name = "travel";
    travel.joint = {.type = JointType::Prismatic, .axis = {1.0, 0.0, 0.0}};
    const BodyId t = builder.addBody(travel);
    BodyDescription heave = massless;
    heave.name = "heave";
    heave.joint = {.type = JointType::Prismatic, .parent = t, .axis = {0.0, 0.0, 1.0}};
    const BodyId h = builder.addBody(heave);
    const BodyId frame = builder.addBody({.name = "frame",
                                          .mass = kilograms(20.0),
                                          .centreOfMass = {metres(0.1), metres(0.0), metres(0.2)},
                                          .inertia = diagonalInertia(1.0, 1.5, 1.0),
                                          .joint = {.type = JointType::Revolute, .parent = h, .axis = {0.0, 1.0, 0.0}}});
    const auto hub = [&](std::string_view name, double x) {
        return builder.addBody({.name = name,
                                .mass = kilograms(1.5),
                                .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                                .inertia = diagonalInertia(0.04, 0.04, 0.04),
                                .joint = {.type = JointType::Revolute,
                                          .parent = frame,
                                          .placement = offset(x, 0.0, -0.3),
                                          .axis = {0.0, 1.0, 0.0}}});
    };
    return {hub("rearHub", -0.55), hub("frontHub", 0.6)};
}

/** Expands a ring description into bodies (a massless carrier on a_0, the node on a_1), foundation coordinate
 * elements and edge tendons: the oracle. Coordinates come out in the ring's order. */
inline void addExpandedRing(ModelBuilder& builder, const DeformableRingDescription& d, const std::string& prefix)
{
    using sub0math::f64::Mass;
    const std::size_t n = d.restPositions.size();
    std::vector<SiteId> sites;
    // Names must outlive the builder's copy: the builder copies them, so locals suffice.
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::string carrierName = prefix + "_c" + std::to_string(i);
        const std::string nodeName = prefix + "_n" + std::to_string(i);
        const auto& r = d.restPositions[i];
        BodyDescription carrier{};
        carrier.name = carrierName;
        carrier.joint = {.type = JointType::Prismatic,
                         .parent = d.host,
                         .placement = {{}, {r.x, r.y, r.z}},
                         .axis = d.displacementAxes[0]};
        const BodyId c = builder.addBody(carrier);
        const double j = d.nodeInertia[i].si();
        const BodyId node = builder.addBody({.name = nodeName,
                                             .mass = d.nodeMass[i],
                                             .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                                             .inertia = diagonalInertia(j, j, j),
                                             .joint = {.type = JointType::Prismatic,
                                                       .parent = c,
                                                       .axis = d.displacementAxes[1]}});
        builder.addElement(CoordinateElementDescription<sub0physics::LinearSpringDamper<sub0math::f64::Length>>{
            .name = carrierName + "_k", .body = c, .law = d.foundation[i], .reference = 0.0});
        builder.addElement(CoordinateElementDescription<sub0physics::LinearSpringDamper<sub0math::f64::Length>>{
            .name = nodeName + "_k", .body = node, .law = d.foundation[i], .reference = 0.0});
        sites.push_back(builder.addSite({.name = nodeName + "_s", .body = node}));
    }
    for (std::size_t e = 0; e < n; ++e)
    {
        const std::string name = prefix + "_e" + std::to_string(e);
        builder.addElement(TendonElementDescription<sub0physics::LinearSpringDamper<sub0math::f64::Length>>{
            .name = name,
            .first = sites[e],
            .second = sites[(e + 1) % n],
            .law = d.edgeLaws[e],
            .sense = LawSense::Increasing,
            .reference = d.edgeRestLength[e]});
    }
}

/** The tree with two rings: condensed substructures, or the expanded oracle. */
inline Model ringRig(const RingParameters& p, bool expanded, bool secondRing = true)
{
    ModelBuilder builder;
    const auto hubs = addTree(builder, p.gravity);
    const auto rear = uniformRing(ringSpec("rear", hubs[0], p));
    const auto front = uniformRing(ringSpec("front", hubs[1], p));
    if (expanded)
    {
        addExpandedRing(builder, *rear, "rear");
        if (secondRing)
        {
            addExpandedRing(builder, *front, "front");
        }
    }
    else
    {
        builder.addRing(*rear);
        if (secondRing)
        {
            builder.addRing(*front);
        }
    }
    return *builder.build();
}

/** A deterministic state: the tree moving, the ring nodes displaced by a few millimetres and moving. */
inline State randomState(const Model& model, unsigned seed, double displacement = 0.004)
{
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    State s = State::atRest(model);
    const std::size_t n = model.dofCount();
    auto q = s.q();
    auto qd = s.qd();
    q[0] = 0.2 * u(rng);
    q[1] = 0.4 + 0.05 * u(rng);
    q[2] = 0.3 * u(rng);
    q[3] = 3.0 * u(rng);
    q[4] = 3.0 * u(rng);
    qd[0] = 2.0 * u(rng);
    qd[1] = 0.5 * u(rng);
    qd[2] = 0.5 * u(rng);
    qd[3] = 5.0 * u(rng);
    qd[4] = 5.0 * u(rng);
    for (std::size_t i = TreeDof; i < n; ++i)
    {
        q[i] = displacement * u(rng);
        qd[i] = 0.2 * u(rng);
    }
    return s;
}

} // namespace ringmodels
