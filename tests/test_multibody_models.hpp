#pragma once

// Small multibody models shared by the multibody tests, the allocation gate and the benchmark.

#include <array>
#include <cmath>
#include <cstddef>
#include <string>

#include <sub0math/linear/quat.hpp>
#include <sub0math/units.hpp>

#include "sub0physics/multibody.hpp"

namespace testmodels
{

using namespace sub0physics::multibody;
using Pose = sub0math::RigidTransform<double, sub0math::f64::Length>;

inline sub0math::f64::Length metres(double v)
{
    return sub0math::f64::Length::fromSi(v);
}

inline sub0math::Vec3<sub0math::f64::Acceleration> gravityAlong(double x, double y, double z)
{
    using A = sub0math::f64::Acceleration;
    return {A::fromSi(x), A::fromSi(y), A::fromSi(z)};
}

inline RotationalInertia diagonalInertia(double xx, double yy, double zz)
{
    using sub0math::kilogramMetresSquared;
    return {{kilogramMetresSquared(xx), kilogramMetresSquared(yy), kilogramMetresSquared(zz)}};
}

inline Pose offset(double x, double y, double z)
{
    return {{}, {metres(x), metres(y), metres(z)}};
}

/** Planar link hanging along -y from its joint: centre of mass at (0, -length / 2, 0). */
inline BodyDescription planarLink(std::string_view name, BodyId parent, const Pose& placement, double mass,
                                  double length)
{
    const double izz = mass * length * length / 12.0;
    return {.name = name,
            .mass = sub0math::kilograms(mass),
            .centreOfMass = {metres(0.0), metres(-0.5 * length), metres(0.0)},
            .inertia = diagonalInertia(izz, 0.1 * izz + 1e-6, izz),
            .joint = {.type = JointType::Revolute, .parent = parent, .placement = placement, .axis = {0.0, 0.0, 1.0}}};
}

/** Number of coordinates of planarTree(). */
inline constexpr std::size_t PlanarTreeDof = 10;

/** Parent of each body of planarTree() (NoParent for the cart). */
inline constexpr std::array<std::size_t, PlanarTreeDof> PlanarTreeParents{NoParent, 0, 1, 2, 0, 4, 5, 0, 7, 8};

/** Link length of the revolute bodies of planarTree(). */
inline constexpr double PlanarTreeLink = 0.4;

/**
 * A fixed 10-DOF planar tree in the xy plane, gravity along -y: a cart sliding along x carrying three
 * hanging branches of three links each, every link a revolute joint about z.
 */
inline ModelBuilder planarTreeBuilder()
{
    ModelBuilder builder;
    builder.setGravity(gravityAlong(0.0, -9.81, 0.0));
    const BodyId cart = builder.addBody({.name = "cart",
                                         .mass = sub0math::kilograms(3.0),
                                         .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                                         .inertia = diagonalInertia(0.05, 0.05, 0.05),
                                         .joint = {.type = JointType::Prismatic, .axis = {1.0, 0.0, 0.0}}});
    const std::array<double, 3> roots{-0.3, 0.0, 0.35};
    for (std::size_t b = 0; b < 3; ++b)
    {
        BodyId parent = cart;
        for (std::size_t k = 0; k < 3; ++k)
        {
            const std::string name = "link" + std::to_string(b) + std::to_string(k);
            const Pose placement = k == 0 ? offset(roots[b], 0.0, 0.0) : offset(0.0, -PlanarTreeLink, 0.0);
            parent = builder.addBody(planarLink(name, parent, placement, 0.8 - 0.1 * static_cast<double>(k), PlanarTreeLink));
        }
    }
    builder.addSite({.name = "tip", .body = BodyId{9}, .frame = offset(0.0, -PlanarTreeLink, 0.0)});
    return builder;
}

/** The fixed 10-DOF planar tree; see planarTreeBuilder. */
inline Model planarTree()
{
    return *planarTreeBuilder().build();
}

/** Parameters of the two elements of planarTreeWithElements. */
struct PlanarTreeElements
{
    double cartStiffness = 2000.0;   ///< Coordinate spring on the cart [N/m].
    double cartDamping = 20.0;       ///< Its damping [N s/m].
    double cartReference = 0.05;     ///< Its relaxed cart position [m].
    double strutPreload = 40.0;      ///< Tendon preload, spreading [N].
    double strutStiffness = 3000.0;  ///< Tendon rate [N/m].
    double strutDamping = 30.0;      ///< Tendon damping, both ways [N s/m].
    double strutReference = 0.65;    ///< Tendon length at zero compression [m].
};

/**
 * The planar tree with two force elements, as a builder: a linear spring-damper on the cart's coordinate, and a preloaded
 * strut (an axial spring-damper on compression) between the tips of the first and the last branch, which
 * couples the branches.
 */
inline ModelBuilder planarTreeWithElementsBuilder(const PlanarTreeElements& p = {})
{
    using namespace sub0math;
    ModelBuilder builder = planarTreeBuilder();
    const SiteId first = builder.addSite({.name = "tip0", .body = BodyId{3}, .frame = offset(0.0, -PlanarTreeLink, 0.0)});
    const SiteId last{0}; // "tip", on the last branch
    builder.addElement(CoordinateElementDescription<sub0physics::LinearSpringDamper<f64::Length>>{
        .name = "cart spring",
        .body = BodyId{0},
        .law = *sub0physics::LinearSpringDamper<f64::Length>::create(
            {.stiffness = f64::Stiffness::fromSi(p.cartStiffness), .damping = f64::Damping::fromSi(p.cartDamping)}),
        .reference = p.cartReference});
    builder.addElement(TendonElementDescription<sub0physics::AxialSpringDamper>{
        .name = "strut",
        .first = first,
        .second = last,
        .law = *sub0physics::AxialSpringDamper::create({.preload = f64::Force::fromSi(p.strutPreload),
                                                        .stiffness = f64::Stiffness::fromSi(p.strutStiffness),
                                                        .compressionDamping = f64::Damping::fromSi(p.strutDamping),
                                                        .reboundDamping = f64::Damping::fromSi(p.strutDamping)}),
        .sense = LawSense::Decreasing,
        .reference = metres(p.strutReference)});
    return builder;
}

/** The planar tree with its two elements; see planarTreeWithElementsBuilder. */
inline Model planarTreeWithElements(const PlanarTreeElements& p = {})
{
    return *planarTreeWithElementsBuilder(p).build();
}

/** Geometry and masses of sliderCrank(). */
struct SliderCrankSpec
{
    double crank = 0.1;       ///< Crank length a [m].
    double rod = 0.3;         ///< Rod length b [m].
    double crankMass = 1.0;   ///< [kg]
    double rodMass = 0.5;     ///< [kg]
    double sliderMass = 2.0;  ///< [kg]
};

/** Planar link along +x from its joint: centre of mass at (length / 2, 0, 0). */
inline BodyDescription linkAlongX(std::string_view name, BodyId parent, const Pose& placement, double mass, double length)
{
    const double izz = mass * length * length / 12.0;
    return {.name = name,
            .mass = sub0math::kilograms(mass),
            .centreOfMass = {metres(0.5 * length), metres(0.0), metres(0.0)},
            .inertia = diagonalInertia(0.1 * izz + 1e-6, izz, izz),
            .joint = {.type = JointType::Revolute, .parent = parent, .placement = placement, .axis = {0.0, 0.0, 1.0}}};
}

/**
 * A planar slider-crank in the xy plane, gravity along -y: a crank (body 0) turning about z at the origin, a rod
 * (body 1) pinned to the crank's tip, and a slider (body 2) moving along the world x axis. The loop is closed by a
 * point closure between the rod's tip (site "rod tip") and the slider's origin (site "slider pin"); its z row is
 * dependent. Coordinates: crank angle, rod angle relative to the crank, slider position.
 */
inline ModelBuilder sliderCrankBuilder(const SliderCrankSpec& spec = {})
{
    ModelBuilder builder;
    builder.setGravity(gravityAlong(0.0, -9.81, 0.0));
    const BodyId crank = builder.addBody(linkAlongX("crank", WorldBody, {}, spec.crankMass, spec.crank));
    const BodyId rod = builder.addBody(linkAlongX("rod", crank, offset(spec.crank, 0.0, 0.0), spec.rodMass, spec.rod));
    const double slider = 0.01 * spec.sliderMass;
    const BodyId block = builder.addBody({.name = "slider",
                                          .mass = sub0math::kilograms(spec.sliderMass),
                                          .centreOfMass = {metres(0.0), metres(0.0), metres(0.0)},
                                          .inertia = diagonalInertia(slider, slider, slider),
                                          .joint = {.type = JointType::Prismatic, .axis = {1.0, 0.0, 0.0}}});
    const SiteId tip = builder.addSite({.name = "rod tip", .body = rod, .frame = offset(spec.rod, 0.0, 0.0)});
    const SiteId pin = builder.addSite({.name = "slider pin", .body = block, .frame = {}});
    builder.addConstraint(LoopClosureDescription{.name = "loop", .first = tip, .second = pin});
    return builder;
}

/**
 * A state of sliderCrank() on the closure manifold: crank angle and rate, the rod and slider from the closed form.
 *
 * @param theta Crank angle [rad].
 * @param rate Crank rate [rad/s].
 */
inline State sliderCrankState(const SliderCrankSpec& spec, double theta, double rate)
{
    const double a = spec.crank;
    const double b = spec.rod;
    const double phi = -std::asin(a * std::sin(theta) / b); // rod's world angle
    const double phiRate = -a * std::cos(theta) * rate / (b * std::cos(phi));
    const double x = a * std::cos(theta) + b * std::cos(phi);
    const double xRate = -a * std::sin(theta) * rate - b * std::sin(phi) * phiRate;
    return State::fromCoordinates({theta, phi - theta, x}, {rate, phiRate - rate, xRate});
}

} // namespace testmodels
