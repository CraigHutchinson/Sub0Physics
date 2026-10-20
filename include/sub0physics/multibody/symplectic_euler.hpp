#pragma once

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
 * Advances a state by one explicit semi-implicit Euler step: `qd += h qdd(q, qd)`, then `q += h qd`
 * with the new rate, and `tick += 1`.
 *
 * This is semi-implicit (symplectic-style) Euler, but it is **not symplectic** for a configuration-
 * dependent mass matrix `M(q)`: the scheme is symplectic only for separable Hamiltonians, and a
 * multibody's kinetic energy depends on `q`. Energy is therefore bounded to first order in `h` (it
 * drifts with `h`, it is not conserved to round-off). Do not use it where long-horizon energy
 * conservation matters; use `Rk4` through an `OdeSystem` adapter, or the linearly implicit integrator
 * (M1b) for stiff elements.
 *
 * First order, explicit and cheap: one forward-dynamics evaluation per step. It is the reference
 * stepping for the M1a analytic gates. The step is transactional: on a fault the state is unchanged.
 *
 * @param model The model.
 * @param h Step size; finite and positive.
 * @param tau Applied generalized forces, dofCount entries; empty for none.
 * @param state State at tick; replaced by the state at tick + 1 only on success.
 * @param workspace Workspace built for the model.
 * @return Nothing, or the fault: InvalidInput (step size, sizes, a workspace built for another model, mass
 *         matrix) or NonFiniteState.
 */
[[nodiscard]] std::expected<void, Fault> stepSymplecticEuler(const Model& model, sub0math::f64::Time h,
                                                             std::span<const double> tau, State& state,
                                                             Workspace& workspace) noexcept;

} // namespace sub0physics::multibody
