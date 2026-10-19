#pragma once

// Small multibody models shared by the multibody tests, the allocation gate and the benchmark.

#include <array>
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
inline Model planarTree()
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
    return *builder.build();
}

} // namespace testmodels
