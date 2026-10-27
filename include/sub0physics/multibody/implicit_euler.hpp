#pragma once

#include <cstdint>
#include <expected>
#include <span>

#include <sub0math/units/quantities.hpp>

#include "sub0physics/fault.hpp"
#include "sub0physics/multibody/kinematics.hpp"
#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/state.hpp"
#include "sub0physics/multibody/stepper.hpp"
#include "sub0physics/multibody/workspace.hpp"

namespace sub0physics::multibody
{

/**
 * Advances a state by one linearly implicit Euler step on velocities (docs/MULTIBODY.md section 5): the
 * default stepper, stable for stiff force elements because they enter through their own Jacobians.
 *
 * With `K = dtau_e/dq` and `D = dtau_e/dqd` summed over the model's elements, it solves
 *
 *     (M - h D - h^2 K~) dqd = h (tau + tau_e - c + h K~ qd)
 *
 * once (no iteration), then sets `qd += dqd` and `q += h qd` (both joint types have one coordinate, so the
 * manifold update is that sum). This is backward Euler linearized at the start of the step. `K~` is the
 * symmetric part of `K`: `K` is exact but not symmetric for a damped tendon, and the step is a W-method, which
 * keeps first order with any Jacobian, so only stability depends on the choice. `M` is evaluated at the start
 * of the step and not differentiated, like MuJoCo's `implicitfast`.
 *
 * **Factorisation.** The matrix is symmetric. While every tendon's ends lie on one chain
 * (Model::elementsFollowTree) it keeps the mass matrix's tree sparsity and Featherstone's LTDL factorises it
 * without fill-in; a tendon between branches couples them, and a dense `LDL^T` (n^3/6, n <= 40) is used.
 * **Indefinite matrix.** Every law here is passive (`dForceDPosition <= 0`, `dForceDVelocity <= 0`), so
 * `M - h D - h^2 K` is positive definite with the material terms alone. A tendon's geometric term
 * `F d2l/dq2` (and its rate transport) is indefinite in general: a compressed strut softens transverse
 * motion. When it outgrows `M / h^2` (an instability faster than the step can follow) the factorisation meets
 * a pivot that is not positive; the step then re-evaluates the elements with the material stiffness only
 * (still first order, and stiff elements stay implicit), sets `Workspace::geometricStiffnessDropped`, and
 * faults only if that fails too.
 *
 * **Energy.** The elements' dissipated power at the start of the step, times h, is added to
 * `state.dissipated` (first order, like the step), and the applied forces' work `tau . (q1 - q0)`, exact
 * for forces constant over the step, to `state.externalWork`. Implicit Euler also damps numerically: the
 * energy balance falls by O(h) per unit time, which the ledger shows as a negative residual.
 *
 * Transactional: on a fault the state is unchanged. Never allocates: every matrix is in the workspace.
 *
 * @param model The model.
 * @param forcing Applied generalized forces, constant over the step.
 * @param tick The tick the state is at; must equal `state.tick`.
 * @param h Step size; finite and positive.
 * @param state State at tick; replaced by the state at tick + 1 only on success.
 * @param workspace Workspace built for the model; on success `qdd` holds `dqd / h` and the element and
 *        implicit-matrix fields hold the start of the step.
 * @param terms Which Jacobians are implicit; the default, Full, is described above. The forcing's Jacobians follow
 *        the same choice (its stiffness is dropped with DampingOnly).
 * @return Nothing, or the fault: InvalidInput (step size, tick, sizes, a workspace built for another model,
 *         a matrix that is not positive definite), OutsideDomain (an element) or NonFiniteState.
 */
[[nodiscard]] std::expected<void, Fault> stepLinearlyImplicitEuler(const Model& model, const Forcing& forcing,
                                                                   std::uint64_t tick, sub0math::f64::Time h,
                                                                   State& state, Workspace& workspace,
                                                                   ImplicitTerms terms = ImplicitTerms::Full) noexcept;

/** The linearly implicit Euler step as a Stepper; its one lever is which Jacobians are implicit. */
struct LinearlyImplicitEuler
{
    ImplicitTerms terms = ImplicitTerms::Full; ///< Which Jacobians are implicit.

    /** Integrator order. */
    static constexpr int Order = 1;

    /**
     * Advances one tick; see stepLinearlyImplicitEuler.
     *
     * @param model The model.
     * @param forcing Applied generalized forces.
     * @param tick The tick the state is at; must equal `state.tick`.
     * @param h Step size; finite and positive.
     * @param state State at tick; replaced by the state at tick + 1 only on success.
     * @param workspace Workspace built for the model.
     * @return Nothing, or the fault.
     */
    [[nodiscard]] std::expected<void, Fault> step(const Model& model, const Forcing& forcing, std::uint64_t tick,
                                                  sub0math::f64::Time h, State& state,
                                                  Workspace& workspace) const noexcept
    {
        return stepLinearlyImplicitEuler(model, forcing, tick, h, state, workspace, terms);
    }
};

static_assert(Stepper<LinearlyImplicitEuler>);

/**
 * Advances a state by one linearly implicit midpoint step: the energy-better model beside linearly implicit Euler
 * (docs/MULTIBODY.md, M1d status; docs/LOD.md, reference and improved models).
 *
 * The same linearisation with weight 1/2, the trapezoid rule on velocities:
 *
 *     (M - h/2 D - h^2/4 K~) dqd = h (tau + tau_e - c + h/2 K~ qd),   q += h (qd + dqd / 2),   qd += dqd
 *
 * It is A-stable like the Euler form but adds no numerical damping to a linear oscillator (backward Euler damps a
 * mode of frequency w at about w^2 h / 2 per second), and it is second order for constant M. `M` and `c` are taken
 * at the start of the step, so a configuration-dependent mass matrix keeps first-order terms of size h times the
 * rate of change of M. Loop closures and limits are not supported yet (InvalidInput): its constrained form is a
 * RATTLE-type scheme, still open. Energy bookkeeping, transactional behaviour and allocation are as
 * stepLinearlyImplicitEuler.
 *
 * @param model The model; without constraints.
 * @param forcing Applied generalized forces and site forces, constant over the step.
 * @param tick The tick the state is at; must equal `state.tick`.
 * @param h Step size; finite and positive.
 * @param state State at tick; replaced by the state at tick + 1 only on success.
 * @param workspace Workspace built for the model.
 * @param terms Which Jacobians are implicit.
 * @return Nothing, or the fault, as stepLinearlyImplicitEuler; InvalidInput for a model with constraints.
 */
[[nodiscard]] std::expected<void, Fault> stepLinearlyImplicitMidpoint(const Model& model, const Forcing& forcing,
                                                                      std::uint64_t tick, sub0math::f64::Time h,
                                                                      State& state, Workspace& workspace,
                                                                      ImplicitTerms terms = ImplicitTerms::Full) noexcept;

/** The linearly implicit midpoint step as a Stepper. */
struct LinearlyImplicitMidpoint
{
    /** Integrator order (for constant M). */
    static constexpr int Order = 2;

    ImplicitTerms terms = ImplicitTerms::Full; ///< Which Jacobians are implicit.

    /**
     * Advances one tick; see stepLinearlyImplicitMidpoint.
     *
     * @param model The model.
     * @param forcing Applied forces.
     * @param tick The tick the state is at; must equal `state.tick`.
     * @param h Step size; finite and positive.
     * @param state State at tick; replaced by the state at tick + 1 only on success.
     * @param workspace Workspace built for the model.
     * @return Nothing, or the fault.
     */
    [[nodiscard]] std::expected<void, Fault> step(const Model& model, const Forcing& forcing, std::uint64_t tick,
                                                  sub0math::f64::Time h, State& state,
                                                  Workspace& workspace) const noexcept
    {
        return stepLinearlyImplicitMidpoint(model, forcing, tick, h, state, workspace, terms);
    }
};

static_assert(Stepper<LinearlyImplicitMidpoint>);

/** The weight of a linearly implicit step: Euler (1) or the midpoint rule (1/2). */
enum class ImplicitScheme : std::uint8_t
{
    Euler,   ///< stepLinearlyImplicitEuler.
    Midpoint ///< stepLinearlyImplicitMidpoint.
};

/**
 * Prepares a linearly implicit step at a state, the first half of a two-phase step: transforms, velocities, bias forces,
 * the mass matrix and the elements with their Jacobians, as the one-call steppers compute them at the start of the
 * step. Between this and completeLinearlyImplicitStep the caller may read the prepared state with
 * computePreparedSiteKinematics (sites and energy, at no further kinematics cost) to form its forcing, as a contact
 * model does; nothing else may use the workspace, and any other entry point that writes it withdraws the preparation.
 *
 * prepare then complete gives exactly the one-call step (bit for bit).
 *
 * @param model The model; without constraints (the constrained steps need their own preparation).
 * @param state The state the step will start from.
 * @param workspace Workspace built for the model.
 * @param terms Which Jacobians will be implicit.
 * @return Nothing, or InvalidInput (sizes, another model's workspace, a model with constraints) or an element fault.
 */
[[nodiscard]] std::expected<void, Fault> prepareLinearlyImplicitStep(const Model& model, const State& state,
                                                                     Workspace& workspace,
                                                                     ImplicitTerms terms = ImplicitTerms::Full) noexcept;

/**
 * Reads site kinematics, and optionally the energy account, from a prepared step (prepareLinearlyImplicitStep), with
 * no further kinematics pass. Leaves the preparation in place.
 *
 * @param model The model.
 * @param state The state the step was prepared at.
 * @param sites Sites of the model.
 * @param out One result per site.
 * @param workspace The prepared workspace.
 * @param energy When given, receives computeEnergy's account of the state.
 * @return Nothing, or InvalidInput when the workspace holds no preparation at this state, or for a site or output size.
 */
[[nodiscard]] std::expected<void, Fault> computePreparedSiteKinematics(const Model& model, const State& state,
                                                                       std::span<const SiteId> sites,
                                                                       std::span<SiteKinematics> out,
                                                                       const Workspace& workspace,
                                                                       EnergyAccount* energy = nullptr) noexcept;

/**
 * Reads site positions and velocities (no orientations, the cheapest site query) from a prepared step, and optionally
 * the energy account; see computePreparedSiteKinematics.
 *
 * @param model The model.
 * @param state The state the step was prepared at.
 * @param sites Sites of the model.
 * @param out One result per site.
 * @param workspace The prepared workspace.
 * @param energy When given, receives computeEnergy's account of the state.
 * @return Nothing, or InvalidInput when the workspace holds no preparation at this state, or for a site or output size.
 */
[[nodiscard]] std::expected<void, Fault> computePreparedSiteMotion(const Model& model, const State& state,
                                                                   std::span<const SiteId> sites,
                                                                   std::span<SiteMotion> out, const Workspace& workspace,
                                                                   EnergyAccount* energy = nullptr) noexcept;

/** What completing a prepared step does to the preparation. */
enum class PreparedStepUse : std::uint8_t
{
    Consume, ///< The preparation is used up (the state it was made at is replaced).
    Keep     ///< The preparation stays valid for the starting state, so the step can be completed again from a copy of
             ///< it, with another step size or scheme (step doubling, an embedded pair): it costs one solve, not a
             ///< second preparation. If a completion fell back to the material stiffness (geometric term dropped,
             ///< Workspace::geometricStiffnessDropped), later completions from the same preparation use it too.
};

/**
 * Completes a prepared step (prepareLinearlyImplicitStep): forms and solves the step with the forcing and commits it,
 * as stepLinearlyImplicitEuler or stepLinearlyImplicitMidpoint would from the same state. By default the preparation
 * is used up; with PreparedStepUse::Keep it stays for further completions from the same starting values (each one bit
 * for bit what a fresh preparation would give, except after a geometric-stiffness fallback, see PreparedStepUse).
 *
 * @param model The model.
 * @param forcing Applied forces and site forces, constant over the step.
 * @param tick The tick the state is at; must equal `state.tick`.
 * @param h Step size; finite and positive.
 * @param state The state the step was prepared at; replaced by the state at tick + 1 only on success.
 * @param workspace The prepared workspace.
 * @param scheme Euler or midpoint.
 * @param use Whether the preparation is used up (the default) or kept for another completion.
 * @return Nothing, or InvalidInput when the workspace holds no preparation at this state, or a stepper fault.
 */
[[nodiscard]] std::expected<void, Fault>
completeLinearlyImplicitStep(const Model& model, const Forcing& forcing, std::uint64_t tick, sub0math::f64::Time h,
                             State& state, Workspace& workspace, ImplicitScheme scheme = ImplicitScheme::Euler,
                             PreparedStepUse use = PreparedStepUse::Consume) noexcept;

} // namespace sub0physics::multibody
