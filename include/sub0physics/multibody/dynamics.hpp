#pragma once

#include <cstdint>
#include <expected>
#include <span>

#include "sub0physics/fault.hpp"
#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/state.hpp"
#include "sub0physics/multibody/workspace.hpp"

// Rigid-body dynamics of a tree, M(q) qdd + c(q, qd) = tau + tau_e(q, qd) (docs/MULTIBODY.md section 3), where
// tau is applied and tau_e comes from the model's force elements (section 4), in
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

/** Whether evaluateElements also computes the elements' Jacobians. */
enum class ElementJacobians : std::uint8_t
{
    Skip,   ///< Force, energy, power and site wrenches only.
    Compute ///< Also `elementStiffness` (dtau/dq) and `elementDamping` (dtau/dqd).
};

/**
 * Evaluates the model's force elements at a state (docs/MULTIBODY.md section 4).
 *
 * Writes `elementForce`, `elementPotential`, `elementDissipationPower`, `elementSiteWrenches` and, with
 * ElementJacobians::Compute, `elementStiffness` and `elementDamping` as dense `dofCount x dofCount` sums
 * (each element adds only on the coordinates it touches). Also writes the kinematics when the model has
 * tendons. Never allocates.
 *
 * @param model The model.
 * @param state Coordinates and rates.
 * @param workspace Workspace built for the model.
 * @param jacobians Whether to compute the Jacobians.
 * @return Nothing, or a fault: InvalidInput for a workspace built for another model or a mis-sized state,
 *         OutsideDomain when a tendon's two sites coincide.
 */
[[nodiscard]] std::expected<void, Fault> evaluateElements(const Model& model, const State& state, Workspace& workspace,
                                                          ElementJacobians jacobians = ElementJacobians::Compute) noexcept;

/**
 * Computes the applied generalized forces that produce given accelerations: tau = M qdd + c - tau_e, so
 * that forwardDynamics with this tau returns qdd. tau_e is the model's force elements (zero without any).
 *
 * @param model The model.
 * @param state Coordinates and rates.
 * @param qdd Generalized accelerations, dofCount entries.
 * @param tau Output, dofCount entries.
 * @param workspace Workspace built for the model.
 * @return Nothing, or InvalidInput if a span has the wrong size, the workspace was built for another
 *         model or the state is not sized for this one; OutsideDomain from an element.
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
 * Computes the generalized accelerations qdd = M^-1 (tau + tau_e - c + J^T lambda), where tau_e is the model's
 * force elements and lambda the forces of its constraints (docs/MULTIBODY.md section 6).
 *
 * Constraints act at acceleration level here: every loop closure (`J qdd + Jdot qd = 0`) and every limit that is
 * engaged (gap within ConstraintSettings::eventTolerance) and not separating, as a small LCP. An impact needs an
 * impulse, which only a step gives: the steppers solve constraints at velocity level instead. The forces are
 * written to `constraintForce` in the model's constraint layout.
 *
 * Writes `qdd`, and on the way `massMatrix`, `massFactor`, `bias`, the kinematics and, when the model has
 * elements, their force, energy terms and site wrenches (not their Jacobians).
 *
 * @param model The model.
 * @param state Coordinates and rates.
 * @param tau Applied generalized forces, dofCount entries; empty for none.
 * @param workspace Workspace built for the model.
 * @return Nothing, or a fault: InvalidInput for a wrong-sized span or state, a workspace built for
 *         another model, or a mass matrix that is not positive definite; OutsideDomain from an element.
 */
[[nodiscard]] std::expected<void, Fault> forwardDynamics(const Model& model, const State& state,
                                                         std::span<const double> tau, Workspace& workspace) noexcept;

/**
 * Projects a state onto its model's constraints (docs/MULTIBODY.md section 6): positions onto the loop closures and
 * the violated limits by Gauss-Newton in the mass metric, `dq = -M^-1 J^T (J M^-1 J^T)^-1 phi`, until the residual
 * meets ConstraintSettings::closureTolerance; coordinate limits still violated are then set to their bound; and
 * velocities onto the closures' tangent space. Use it to assemble a mechanism from an approximate configuration;
 * the steppers project every accepted step themselves.
 *
 * The achieved residuals and iterations are written to `constraintReport`; the energy the projection changed is
 * added to `state.projectionEnergy` (when ConstraintSettings::measureProjectionEnergy).
 *
 * @param model The model.
 * @param state The state; projected in place only on success.
 * @param workspace Workspace built for the model.
 * @return Nothing, or a fault: InvalidInput for a workspace built for another model or a mis-sized state;
 *         OutsideDomain when the tolerance is not met within ConstraintSettings::maxProjectionIterations, or a
 *         distance limit's sites coincide.
 */
[[nodiscard]] std::expected<void, Fault> projectConstraints(const Model& model, State& state,
                                                            Workspace& workspace) noexcept;

/**
 * Computes the joint reactions and constraint site loads at a state (docs/MULTIBODY.md section 7): the constrained
 * accelerations from forwardDynamics, then a recursive Newton-Euler pass that includes every force: gravity, the
 * applied forces, the elements (coordinate elements in their joints, tendons as wrenches at their sites) and the
 * constraint forces (closures and distance limits as wrenches at their sites, coordinate limits in their joints).
 * No force category is left out.
 *
 * Writes `jointReactions` (the wrench each joint transmits from the parent to its body, world axes, about the joint
 * origin), `constraintSiteWrenches`, `constraintForce`, `qdd` and, through forwardDynamics, the element site
 * wrenches. A joint reaction's projection on the joint's motion axis equals the applied, coordinate-element and
 * coordinate-limit force on that coordinate.
 *
 * @param model The model.
 * @param state Coordinates and rates.
 * @param tau Applied generalized forces, dofCount entries; empty for none.
 * @param workspace Workspace built for the model.
 * @return Nothing, or the forwardDynamics fault.
 */
[[nodiscard]] std::expected<void, Fault> computeReactions(const Model& model, const State& state,
                                                          std::span<const double> tau, Workspace& workspace) noexcept;

} // namespace sub0physics::multibody
