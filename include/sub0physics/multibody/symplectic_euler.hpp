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
 * Advances a state by one explicit semi-implicit (symplectic) Euler step: `qd += h qdd(q, qd)`, then
 * `q += h qd` with the new rate, and `tick += 1`.
 *
 * First order, explicit and cheap: one forward-dynamics evaluation per step. It is the reference
 * stepping for the M1a analytic gates; stiff force elements need the linearly implicit integrator
 * (M1b). The step is transactional: on a fault the state is unchanged.
 *
 * @param model The model.
 * @param h Step size; finite and positive.
 * @param tau Applied generalized forces, dofCount entries; empty for none.
 * @param state State at tick; replaced by the state at tick + 1 only on success.
 * @param workspace Workspace sized for the model.
 * @return Nothing, or the fault: InvalidInput (step size, sizes, mass matrix) or NonFiniteState.
 */
[[nodiscard]] std::expected<void, Fault> stepSymplecticEuler(const Model& model, sub0math::f64::Time h,
                                                             std::span<const double> tau, State& state,
                                                             Workspace& workspace) noexcept;

} // namespace sub0physics::multibody
