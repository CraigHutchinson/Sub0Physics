#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>

#include "sub0physics/fault.hpp"
#include "sub0physics/integrator/ode_system.hpp"

namespace sub0physics
{

/** The explicit Runge-Kutta tableaux this library provides. */
enum class RungeKuttaMethod : std::uint8_t
{
    Midpoint,  ///< Explicit midpoint: second order, two evaluations per step.
    ClassicRk4 ///< Classic fourth-order method: four evaluations per step.
};

/**
 * Transactional fixed-step explicit Runge-Kutta integrator that never allocates.
 *
 * A step either commits a finite endpoint that the system accepted, or leaves the caller's state
 * untouched and the system observing that unchanged state. The derivative at a committed endpoint is
 * reused as the first stage of the next step when the caller passes the same time and state, so the
 * steady-state cost is one evaluation per stage.
 *
 * @tparam N Number of state variables.
 * @tparam Method The tableau.
 * @note Not thread-safe; use one instance per integrated system.
 */
template<std::size_t N, RungeKuttaMethod Method>
class ExplicitRungeKutta
{
public:
    /** Global order of accuracy of the method. */
    static constexpr int Order = Method == RungeKuttaMethod::ClassicRk4 ? 4 : 2;

    /**
     * Advances `state` from `t` to `t + h`.
     *
     * @tparam S The system type.
     * @param system The system to integrate; it may refuse any trial or endpoint state.
     * @param t Time at the start of the step.
     * @param h Step size; finite and positive.
     * @param state State at `t`; replaced by the state at `t + h` only on success.
     * @return Nothing on success; otherwise the first fault, with `state` unchanged.
     */
    template<OdeSystem<N> S>
    [[nodiscard]] std::expected<void, Fault> step(S& system, double t, double h, StateVector<N>& state)
    {
        if (!std::isfinite(t) || !std::isfinite(h) || !(h > 0.0))
        {
            return std::unexpected(Fault{FaultKind::InvalidInput, "step needs a finite time and positive step"});
        }
        auto result = advance(system, t, h, state);
        if (!result)
        {
            // Leave the system observing the last valid state. The trial's fault is the one to report.
            cacheValid_ = system.evaluate(t, state, k1_).has_value();
            cachedTime_ = t;
            cachedState_ = state;
            return result;
        }
        state = next_;
        cachedTime_ = t + h;
        cachedState_ = state;
        k1_ = endDerivative_;
        cacheValid_ = true;
        return {};
    }

    /** Forgets the cached first stage; call after changing the system's parameters in place. */
    void invalidate() noexcept
    {
        cacheValid_ = false;
    }

private:
    template<typename S>
    std::expected<void, Fault> advance(S& system, double t, double h, const StateVector<N>& x)
    {
        if (!cacheValid_ || cachedTime_ != t || cachedState_ != x)
        {
            if (auto r = system.evaluate(t, x, k1_); !r)
            {
                return r;
            }
        }
        const double half = 0.5 * h;
        stage(x, k1_, half);
        if (auto r = system.evaluate(t + half, trial_, k2_); !r)
        {
            return r;
        }
        if constexpr (Method == RungeKuttaMethod::Midpoint)
        {
            for (std::size_t i = 0; i < N; ++i)
            {
                next_[i] = x[i] + h * k2_[i];
            }
        }
        else
        {
            stage(x, k2_, half);
            if (auto r = system.evaluate(t + half, trial_, k3_); !r)
            {
                return r;
            }
            stage(x, k3_, h);
            if (auto r = system.evaluate(t + h, trial_, k4_); !r)
            {
                return r;
            }
            for (std::size_t i = 0; i < N; ++i)
            {
                next_[i] = x[i] + h * (k1_[i] + 2.0 * k2_[i] + 2.0 * k3_[i] + k4_[i]) / 6.0;
            }
        }
        for (const double value : next_)
        {
            if (!std::isfinite(value))
            {
                return std::unexpected(Fault{FaultKind::NonFiniteState, "nonfinite physical state"});
            }
        }
        return system.evaluate(t + h, next_, endDerivative_);
    }

    void stage(const StateVector<N>& x, const StateVector<N>& k, double scale) noexcept
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            trial_[i] = x[i] + scale * k[i];
        }
    }

    StateVector<N> k1_{};
    StateVector<N> k2_{};
    StateVector<N> k3_{};
    StateVector<N> k4_{};
    StateVector<N> trial_{};
    StateVector<N> next_{};
    StateVector<N> endDerivative_{};
    StateVector<N> cachedState_{};
    double cachedTime_ = 0.0;
    bool cacheValid_ = false;
};

/**
 * Classic fourth-order Runge-Kutta; the reference tier.
 *
 * @tparam N Number of state variables.
 */
template<std::size_t N>
using Rk4 = ExplicitRungeKutta<N, RungeKuttaMethod::ClassicRk4>;

/**
 * Explicit midpoint; half the cost of Rk4, and exact while acceleration is constant.
 *
 * @tparam N Number of state variables.
 */
template<std::size_t N>
using Midpoint = ExplicitRungeKutta<N, RungeKuttaMethod::Midpoint>;

} // namespace sub0physics
