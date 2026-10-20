#pragma once

#include <cmath>
#include <cstddef>
#include <expected>

#include <sub0math/config.hpp>
#include <sub0math/linear/mat.hpp>
#include <sub0math/linear/plucker_transform.hpp>
#include <sub0math/linear/spatial_force.hpp>
#include <sub0math/linear/spatial_inertia.hpp>
#include <sub0math/linear/spatial_motion.hpp>
#include <sub0math/linear/vec.hpp>

#include "sub0physics/fault.hpp"
#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/state.hpp"
#include "sub0physics/multibody/workspace.hpp"

// Joint transforms shared by the kinematics and dynamics translation units.

namespace sub0physics::multibody::detail
{

// The one precondition check every public entry point makes: the workspace belongs to the model and
// the state has one entry per coordinate. A few comparisons; no allocation.
[[nodiscard]] SUB0MATH_FORCE_INLINE std::expected<void, Fault> checkInputs(const Model& model, const State& state,
                                                                          const Workspace& ws) noexcept
{
    if (!ws.isBuiltFor(model))
    {
        return invalidInput("workspace was not built for this model");
    }
    if (!state.isSizedFor(model))
    {
        return invalidInput("state is not sized for the model");
    }
    return {};
}

// Transform from the parent's frame to the body's frame at coordinate q: X_J(q) X_T.
// Revolute: X_J = rot(R(a, q)^T), so X_up = (E_J E_T, r_T). Prismatic: X_J = xlt(a q), so
// X_up = (E_T, r_T + E_T^T a q).
SUB0MATH_FORCE_INLINE sub0math::PluckerTransformd parentToBody(JointType type, const sub0math::Vec3d& a,
                                                               const sub0math::PluckerTransformd& tree, bool unrotated,
                                                               double q) noexcept
{
    if (type == JointType::Revolute)
    {
        // E_J = R^T = c 1 - s [a]x + (1 - c) a a^T (Rodrigues, transposed).
        const double c = std::cos(q);
        const double s = std::sin(q);
        const double t = 1.0 - c;
        sub0math::Mat3d e{};
        e(0, 0) = c + t * a.x * a.x;
        e(1, 1) = c + t * a.y * a.y;
        e(2, 2) = c + t * a.z * a.z;
        e(0, 1) = t * a.x * a.y + s * a.z;
        e(1, 0) = t * a.x * a.y - s * a.z;
        e(0, 2) = t * a.x * a.z - s * a.y;
        e(2, 0) = t * a.x * a.z + s * a.y;
        e(1, 2) = t * a.y * a.z + s * a.x;
        e(2, 1) = t * a.y * a.z - s * a.x;
        // PERF: skipping the identity placement product saves 27 multiplies per joint.
        return {unrotated ? e : e * tree.rotation, tree.translation};
    }
    const sub0math::Vec3d slide = a * q;
    if (unrotated)
    {
        return {tree.rotation, tree.translation + slide};
    }
    return {tree.rotation, tree.translation + sub0math::Vec3d{dot(tree.rotation.columns[0], slide),
                                                              dot(tree.rotation.columns[1], slide),
                                                              dot(tree.rotation.columns[2], slide)}};
}

// S . f, using that S has only an angular (revolute) or only a linear (prismatic) part.
SUB0MATH_FORCE_INLINE double subspaceDot(JointType type, const sub0math::Vec3d& a, const sub0math::SpatialForced& f) noexcept
{
    return type == JointType::Revolute ? dot(a, f.angular) : dot(a, f.linear);
}

// I S, the same way: {I_o a, -h x a} for a revolute joint, {h x a, m a} for a prismatic one.
SUB0MATH_FORCE_INLINE sub0math::SpatialForced inertiaTimesSubspace(JointType type, const sub0math::Vec3d& a,
                                                                   const sub0math::SpatialInertiad& i) noexcept
{
    if (type == JointType::Revolute)
    {
        return {i.rotational * a, cross(a, i.firstMoment)};
    }
    return {cross(i.firstMoment, a), a * i.mass};
}

// Fills parentToBody and, when requested, bodyVelocity for every body.
template<bool WithVelocity>
SUB0MATH_FORCE_INLINE void updateTransforms(const Model& model, const State& state, Workspace& ws) noexcept
{
    const auto parents = model.parents();
    const auto types = model.jointTypes();
    const auto axes = model.jointAxes();
    const auto tree = model.treeTransforms();
    const auto subspace = model.motionSubspaces();
    const auto unrotated = model.placementUnrotated();
    for (std::size_t i = 0; i < parents.size(); ++i)
    {
        const double q = state.q[i];
        ws.parentToBody[i] = parentToBody(types[i], axes[i], tree[i], unrotated[i] != 0U, q);
        if constexpr (WithVelocity)
        {
            const sub0math::SpatialMotiond joint = subspace[i] * state.qd[i];
            ws.bodyVelocity[i] =
                parents[i] == NoParent ? joint : ws.parentToBody[i].apply(ws.bodyVelocity[parents[i]]) + joint;
        }
    }
}

} // namespace sub0physics::multibody::detail
