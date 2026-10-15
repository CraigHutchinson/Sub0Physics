#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <expected>

#include "sub0physics/fault.hpp"

namespace sub0physics
{

/**
 * Fixed-size state of an ordinary differential equation.
 *
 * @tparam N Number of state variables.
 */
template<std::size_t N>
using StateVector = std::array<double, N>;

/**
 * A first-order system dx/dt = f(t, x) that an integrator can advance.
 *
 * `evaluate` writes f(t, x) into its output and may refuse a state outside the model's domain. A
 * system may cache observations of the last evaluated state; integrators finish every successful
 * step by evaluating the committed endpoint, so those observations then describe the committed state.
 *
 * @tparam S The system type.
 * @tparam N Number of state variables.
 */
template<typename S, std::size_t N>
concept OdeSystem = requires(S& system, double t, const StateVector<N>& x, StateVector<N>& dxdt) {
    { system.evaluate(t, x, dxdt) } -> std::same_as<std::expected<void, Fault>>;
};

} // namespace sub0physics
