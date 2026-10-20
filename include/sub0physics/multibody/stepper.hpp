#pragma once

#include <concepts>
#include <cstdint>
#include <expected>
#include <span>

#include <sub0math/units/quantities.hpp>

#include "sub0physics/fault.hpp"
#include "sub0physics/multibody/model.hpp"
#include "sub0physics/multibody/state.hpp"
#include "sub0physics/multibody/workspace.hpp"

namespace sub0physics::multibody
{

/**
 * What drives a model over one step besides its own dynamics.
 *
 * For now, applied generalized forces held constant over the step. Force elements join it in M1b
 * (docs/MULTIBODY.md section 4) without changing the stepper signature.
 */
struct Forcing
{
    std::span<const double> tau; ///< Applied generalized forces, dofCount entries (N or N m); empty for none.
};

/**
 * A multibody stepper (docs/MULTIBODY.md, "M1b contracts"): it advances a state over one fixed tick.
 *
 * `step(model, forcing, tick, h, state, workspace)` advances `state` from time `tick * h` to
 * `(tick + 1) * h` and sets `state.tick = tick + 1`. Time is derived from the tick inside the stepper,
 * never passed separately. It is **transactional**: on a fault the state is unchanged. A stepper never
 * allocates while stepping; any scratch is sized when it is built. `FidelitySelector` can switch between
 * steppers because they share this shape.
 *
 * @tparam S The stepper type.
 */
template<typename S>
concept Stepper = requires(S& stepper, const Model& model, const Forcing& forcing, std::uint64_t tick,
                           sub0math::f64::Time h, State& state, Workspace& workspace) {
    { stepper.step(model, forcing, tick, h, state, workspace) } -> std::same_as<std::expected<void, Fault>>;
};

/**
 * Advances a state by one explicit semi-implicit Euler step: `qd += h qdd(q, qd)`, then `q += h qd`
 * with the new rate, and `tick += 1`.
 *
 * This is semi-implicit (symplectic-style) Euler, but it is **not symplectic** for a configuration-
 * dependent mass matrix `M(q)`: the scheme is symplectic only for separable Hamiltonians, and a
 * multibody's kinetic energy depends on `q`. Energy is therefore bounded to first order in `h` (it
 * drifts with `h`, it is not conserved to round-off). Where long-horizon energy conservation matters use
 * RungeKuttaStepper, or the linearly implicit integrator (M1b) for stiff elements.
 *
 * First order, explicit and cheap: one forward-dynamics evaluation per step. The step is transactional:
 * on a fault the state is unchanged.
 *
 * @param model The model.
 * @param forcing Applied generalized forces.
 * @param tick The tick the state is at; must equal `state.tick`.
 * @param h Step size; finite and positive.
 * @param state State at tick; replaced by the state at tick + 1 only on success.
 * @param workspace Workspace built for the model; holds the dynamics at the start of the step.
 * @return Nothing, or the fault: InvalidInput (step size, tick, sizes, a workspace built for another
 *         model, mass matrix) or NonFiniteState.
 */
[[nodiscard]] std::expected<void, Fault> stepSemiImplicitEuler(const Model& model, const Forcing& forcing,
                                                               std::uint64_t tick, sub0math::f64::Time h,
                                                               State& state, Workspace& workspace) noexcept;

/** The semi-implicit Euler step as a Stepper; stateless. */
struct SemiImplicitEuler
{
    /** Integrator order. */
    static constexpr int Order = 1;

    /**
     * Advances one tick; see stepSemiImplicitEuler.
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
        return stepSemiImplicitEuler(model, forcing, tick, h, state, workspace);
    }
};

static_assert(Stepper<SemiImplicitEuler>);

} // namespace sub0physics::multibody
