#pragma once

#include <cmath>
#include <cstddef>
#include <expected>
#include <span>

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
#include "sub0physics/multibody/stepper.hpp"
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
        sub0math::Mat3d e{};
        if (a.x == 0.0 && a.y == 0.0)
        {
            // PERF: z-axis branch beats general Rodrigues (about 4%).
            const double sz = s * a.z;
            e(0, 0) = c;
            e(1, 1) = c;
            e(2, 2) = 1.0;
            e(0, 1) = sz;
            e(1, 0) = -sz;
            return {unrotated ? e : e * tree.rotation, tree.translation};
        }
        const double t = 1.0 - c;
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

// E^T v: a body-axes vector in world axes, for a world-to-body rotation E, without forming the transpose.
SUB0MATH_FORCE_INLINE sub0math::Vec3d bodyToWorld(const sub0math::Mat3d& e, const sub0math::Vec3d& v) noexcept
{
    return {dot(e.columns[0], v), dot(e.columns[1], v), dot(e.columns[2], v)};
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
    const std::size_t n = parents.size();
    const double* q = state.values.data();
    const double* qd = q + n;
    for (std::size_t i = 0; i < n; ++i)
    {
        ws.parentToBody[i] = parentToBody(types[i], axes[i], tree[i], unrotated[i] != 0U, q[i]);
        if constexpr (WithVelocity)
        {
            const sub0math::SpatialMotiond joint = subspace[i] * qd[i];
            ws.bodyVelocity[i] =
                parents[i] == NoParent ? joint : ws.parentToBody[i].apply(ws.bodyVelocity[parents[i]]) + joint;
        }
    }
}

// Composes each body's world transform from the parent-to-body transforms already in the workspace.
SUB0MATH_FORCE_INLINE void composeWorldTransforms(const Model& model, Workspace& ws) noexcept
{
    const auto parents = model.parents();
    for (std::size_t i = 0; i < parents.size(); ++i)
    {
        const std::size_t p = parents[i];
        ws.worldToBody[i] = p == NoParent ? ws.parentToBody[i] : ws.parentToBody[i] * ws.worldToBody[p];
    }
}

// Sums the model's force elements at a state into the workspace (elementForce, the energy terms, the site
// wrenches and, when asked, elementStiffness and elementDamping, with or without the tendons' geometric
// terms). Needs worldToBody at the state when the model has tendons. OutsideDomain when a tendon's sites
// coincide.
[[nodiscard]] std::expected<void, Fault> evaluateElementsUnchecked(const Model& model, const State& state,
                                                                   Workspace& ws, bool jacobians,
                                                                   bool geometricStiffness = true) noexcept;

// Kinetic and potential energy (gravity and elements) [J].
struct EnergyTerms
{
    double kinetic;
    double potential;
};

// The energy of a state whose transforms, velocities and world transforms are in the workspace (the formula of
// computeEnergy, which calls it). Evaluates the elements' energy terms.
[[nodiscard]] std::expected<EnergyTerms, Fault> energyTerms(const Model& model, const State& state,
                                                            Workspace& workspace) noexcept;

// The pieces of a step that the constraint solver (multibody_constraints.cpp) shares with the steppers
// (multibody_dynamics.cpp), so the formulas exist once.

// Transforms, velocities, world transforms, bias forces (RNEA at zero qdd, so bodyAcceleration holds the bias and
// gravity), the mass matrix and the elements, with Jacobians when asked.
[[nodiscard]] std::expected<void, Fault> prepareStepDynamics(const Model& model, const State& state, Workspace& ws,
                                                             bool elementJacobians) noexcept;

// Forms and factorises the step matrix (M - h D - h^2 K~ when implicit, else M) and its right-hand side
// h (tau + tau_e - c [+ h K~ qd]) in implicitRhs. Returns whether the factor is dense.
[[nodiscard]] std::expected<bool, Fault> formStepMatrix(const Model& model, const State& state,
                                                        std::span<const double> tau, const double* elementForce,
                                                        Workspace& ws, double dt, bool implicit,
                                                        bool allowFallback) noexcept;

// x := A^-1 x with the factor from formStepMatrix.
void solveStepMatrix(const Model& model, const Workspace& ws, bool dense, double* x) noexcept;

// x := M^-1 x with massFactor.
void solveMassFactor(const Model& model, const Workspace& ws, double* x) noexcept;

// RNEA with accelerations and external body wrenches (body frame): bodyForce ends as each joint's transmitted wrench,
// tau as its generalized force. Updates the transforms, velocities and world transforms.
void transmittedForces(const Model& model, const State& state, const double* qdd, double* tau,
                       const sub0math::SpatialForced* external, Workspace& ws) noexcept;

// The constrained step (multibody_constraints.cpp): semi-implicit (implicit = false) or linearly implicit Euler with
// the model's closures and limits. The caller has checked the inputs.
[[nodiscard]] std::expected<void, Fault> stepConstrained(const Model& model, const Forcing& forcing, double dt,
                                                         State& state, Workspace& ws, bool implicit) noexcept;

// Adds the constraint forces to the free accelerations in ws.qdd (acceleration level, closures and engaged limits),
// after forwardDynamics' free solve: massFactor, bias accelerations and kinematics at `state` in the workspace.
[[nodiscard]] std::expected<void, Fault> constrainAccelerations(const Model& model, const State& state,
                                                                Workspace& ws) noexcept;

} // namespace sub0physics::multibody::detail
