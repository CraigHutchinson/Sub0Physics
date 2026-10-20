#pragma once

#include <expected>
#include <span>

#include "sub0physics/fault.hpp"
#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/state.hpp"
#include "sub0physics/multibody/workspace.hpp"

// Rigid-body dynamics of a tree, M(q) qdd + c(q, qd) = tau (docs/MULTIBODY.md section 3), in
// Featherstone's formulation: RNEA for c and inverse dynamics, CRBA for M, and a factorisation of M
// that follows the tree so it creates no fill-in. None of these allocate. Every public entry point checks
// that the workspace was built for the model and that the state is sized for it, and returns InvalidInput
// if not.

namespace sub0physics::multibody
{

/**
 * Computes the joint-space mass matrix with the composite-rigid-body algorithm.
 *
 * Writes `massMatrix` (symmetric, positive definite for a valid model) and the kinematics.
 *
 * @param model The model.
 * @param state Coordinates (rates are not used).
 * @param workspace Workspace built for the model.
 * @return Nothing, or InvalidInput if the workspace was built for another model or the state is not
 *         sized for this one.
 */
[[nodiscard]] std::expected<void, Fault> computeMassMatrix(const Model& model, const State& state,
                                                           Workspace& workspace) noexcept;

/**
 * Computes the bias forces c(q, qd): Coriolis, centrifugal and gravity, with the recursive Newton-Euler
 * algorithm at zero acceleration.
 *
 * Writes `bias` and the kinematics.
 *
 * @param model The model.
 * @param state Coordinates and rates.
 * @param workspace Workspace built for the model.
 * @return Nothing, or InvalidInput if the workspace was built for another model or the state is not
 *         sized for this one.
 */
[[nodiscard]] std::expected<void, Fault> computeBiasForces(const Model& model, const State& state,
                                                           Workspace& workspace) noexcept;

/**
 * Computes the generalized forces that produce given accelerations: tau = M qdd + c.
 *
 * @param model The model.
 * @param state Coordinates and rates.
 * @param qdd Generalized accelerations, dofCount entries.
 * @param tau Output, dofCount entries.
 * @param workspace Workspace built for the model.
 * @return Nothing, or InvalidInput if a span has the wrong size, the workspace was built for another
 *         model or the state is not sized for this one.
 */
[[nodiscard]] std::expected<void, Fault> inverseDynamics(const Model& model, const State& state,
                                                         std::span<const double> qdd, std::span<double> tau,
                                                         Workspace& workspace) noexcept;

/**
 * Factorises `massMatrix` as M = L^T D L, with L unit lower triangular, eliminating from the leaves to
 * the root (Featherstone's LTDL). Along a tree this ordering creates no fill-in: L is non-zero only
 * where a coordinate's body is an ancestor of another's, so cost follows the tree's depth, not n^3.
 *
 * Writes `massFactor`.
 *
 * @param model The model whose `massMatrix` is in the workspace.
 * @param workspace Workspace holding M from computeMassMatrix.
 * @return Nothing, or InvalidInput if the workspace was built for another model or M is not positive
 *         definite.
 */
[[nodiscard]] std::expected<void, Fault> factorMassMatrix(const Model& model, Workspace& workspace) noexcept;

/**
 * Solves M x = b in place with the factor from factorMassMatrix.
 *
 * @param model The model.
 * @param workspace Workspace holding the factor.
 * @param x On entry b, on exit M^-1 b; dofCount entries.
 * @return Nothing, or InvalidInput if `x` has the wrong size or the workspace was built for another model.
 */
[[nodiscard]] std::expected<void, Fault> solveMassMatrix(const Model& model, const Workspace& workspace,
                                                         std::span<double> x) noexcept;

/**
 * Computes the generalized accelerations qdd = M^-1 (tau - c).
 *
 * Writes `qdd`, and on the way `massMatrix`, `massFactor`, `bias` and the kinematics.
 *
 * @param model The model.
 * @param state Coordinates and rates.
 * @param tau Applied generalized forces, dofCount entries; empty for none.
 * @param workspace Workspace built for the model.
 * @return Nothing, or a fault: InvalidInput for a wrong-sized span or state, a workspace built for
 *         another model, or a mass matrix that is not positive definite.
 */
[[nodiscard]] std::expected<void, Fault> forwardDynamics(const Model& model, const State& state,
                                                         std::span<const double> tau, Workspace& workspace) noexcept;

} // namespace sub0physics::multibody
