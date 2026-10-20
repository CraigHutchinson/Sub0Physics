#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <expected>
#include <span>
#include <type_traits>
#include <vector>

#include "sub0physics/fault.hpp"

namespace sub0physics
{

/**
 * State of an ordinary differential equation.
 *
 * A fixed size (`std::array`) for small systems whose size is known at compile time, which the
 * integrators unroll; `std::dynamic_extent` selects a `std::vector`, sized once by its owner, for
 * systems sized at run time (a multibody model).
 *
 * @tparam N Number of state variables, or `std::dynamic_extent`.
 */
template<std::size_t N>
using StateVector = std::conditional_t<N == std::dynamic_extent, std::vector<double>, std::array<double, N>>;

/**
 * A first-order system dx/dt = f(t, x) that an integrator can advance.
 *
 * `evaluate` writes f(t, x) into its output and may refuse a state outside the model's domain. A
 * system may cache observations of the last evaluated state; integrators finish every successful
 * step by evaluating the committed endpoint, so those observations then describe the committed state.
 *
 * A fixed-size system takes `const StateVector<N>&` and `StateVector<N>&`. A dynamic-size system
 * (`N == std::dynamic_extent`) takes `std::span<const double>` and `std::span<double>` of the size
 * its integrator was built with. Time is plain seconds: systems work in SI doubles internally.
 *
 * @tparam S The system type.
 * @tparam N Number of state variables, or `std::dynamic_extent`.
 */
template<typename S, std::size_t N>
concept OdeSystem =
    (N != std::dynamic_extent &&
     requires(S& system, double t, const StateVector<N>& x, StateVector<N>& dxdt) {
         { system.evaluate(t, x, dxdt) } -> std::same_as<std::expected<void, Fault>>;
     }) ||
    (N == std::dynamic_extent && requires(S& system, double t, std::span<const double> x, std::span<double> dxdt) {
        { system.evaluate(t, x, dxdt) } -> std::same_as<std::expected<void, Fault>>;
    });

/**
 * Jacobian of a fixed-size system, row-major: `jacobian[i][j]` is the partial derivative of `dx_i/dt` with
 * respect to `x_j`.
 *
 * @tparam N Number of state variables.
 */
template<std::size_t N>
using JacobianMatrix = std::array<std::array<double, N>, N>;

/**
 * A fixed-size system that also provides its Jacobian `d f / d x` at `(t, x)`, for linearly implicit methods.
 *
 * The Jacobian may be approximate (any matrix keeps a W-method at its order); exact stiff terms are what make
 * large steps stable. Like `evaluate`, `jacobian` may refuse a state outside the model's domain.
 *
 * @tparam S The system type.
 * @tparam N Number of state variables.
 */
template<typename S, std::size_t N>
concept OdeSystemWithJacobian =
    OdeSystem<S, N> && N != std::dynamic_extent &&
    requires(S& system, double t, const StateVector<N>& x, JacobianMatrix<N>& jacobian) {
        { system.jacobian(t, x, jacobian) } -> std::same_as<std::expected<void, Fault>>;
    };

} // namespace sub0physics
