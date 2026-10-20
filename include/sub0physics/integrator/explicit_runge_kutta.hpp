#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <tuple>
#include <type_traits>

#include <sub0math/units/quantities.hpp>

#include "sub0physics/fault.hpp"
#include "sub0physics/integrator/ode_system.hpp"

// Inlines a stepper's whole call tree (the system's evaluate included) into one function.
#if defined(__GNUC__) || defined(__clang__)
#define SUB0PHYSICS_FLATTEN __attribute__((flatten))
#else
#define SUB0PHYSICS_FLATTEN
#endif

namespace sub0physics
{

/** The explicit Runge-Kutta tableaux this library provides. */
enum class RungeKuttaMethod : std::uint8_t
{
    Midpoint,  ///< Explicit midpoint: second order, two evaluations per step.
    ClassicRk4 ///< Classic fourth-order method: four evaluations per step.
};

/**
 * Transactional fixed-step explicit Runge-Kutta integrator that never allocates while stepping.
 *
 * A step either commits a finite endpoint that the system accepted, or leaves the caller's state
 * untouched and the system observing that unchanged state. The derivative at a committed endpoint is
 * reused as the first stage of the next step when the caller passes the same time and state, so the
 * steady-state cost is one evaluation per stage.
 *
 * Two forms share this one algorithm:
 *  - **fixed size** (`N` a number): the state is a `StateVector<N>` (`std::array`); small station
 *    systems use it, and the compiler unrolls its loops;
 *  - **dynamic size** (`N == std::dynamic_extent`): the state is a `std::span<double>` of the size
 *    given at construction, where the scratch is allocated once; it drives systems sized at run time,
 *    such as a multibody model (multibody::RungeKuttaStepper).
 *
 * @tparam N Number of state variables, or `std::dynamic_extent`.
 * @tparam Method The tableau.
 * @note Not thread-safe; use one instance per integrated system.
 */
template<std::size_t N, RungeKuttaMethod Method>
class ExplicitRungeKutta
{
    static constexpr bool Dynamic = N == std::dynamic_extent;

public:
    /** Global order of accuracy of the method. */
    static constexpr int Order = Method == RungeKuttaMethod::ClassicRk4 ? 4 : 2;

    /** The caller's state: the array itself for the fixed size, a span over it for the dynamic size. */
    using StateRef = std::conditional_t<Dynamic, std::span<double>, StateVector<N>&>;

    /** Builds a fixed-size integrator. */
    ExplicitRungeKutta() noexcept
        requires(!Dynamic)
    = default;

    /**
     * Builds a dynamic-size integrator; allocates its scratch once.
     *
     * @param size Number of state variables of every state it will step.
     */
    explicit ExplicitRungeKutta(std::size_t size)
        requires Dynamic
        : k1_(size),
          scratch_{StateVector<N>(size), StateVector<N>(size), StateVector<N>(size), StateVector<N>(size),
                   StateVector<N>(size), StateVector<N>(size)},
          cachedState_(size)
    {
    }

    /** @return Number of state variables. */
    [[nodiscard]] std::size_t size() const noexcept
    {
        if constexpr (Dynamic)
        {
            return k1_.size();
        }
        else
        {
            return N;
        }
    }

    /**
     * Advances `state` from `t` to `t + h`.
     *
     * @tparam S The system type.
     * @param system The system to integrate; it may refuse any trial or endpoint state.
     * @param t Time at the start of the step.
     * @param h Step size; finite and positive.
     * @param state State at `t`; replaced by the state at `t + h` only on success.
     * @return Nothing on success; otherwise the first fault, with `state` unchanged. InvalidInput for a
     *         bad time or step, or a dynamic state of the wrong size.
     */
    template<OdeSystem<N> S>
    [[nodiscard]] std::expected<void, Fault> step(S& system, sub0math::f64::Time t, sub0math::f64::Time h,
                                                  StateRef state)
    {
        return stepBetween(system, t.si(), h.si(), t.si() + h.si(), state);
    }

    /**
     * Advances `state` over fixed-step tick `tick`, from `tick * h` to `(tick + 1) * h`.
     *
     * Times come from the integer tick, never from accumulating `h`, so they do not drift and the
     * endpoint of one tick is exactly the start of the next. That keeps runs reproducible and lets
     * every step reuse the previous endpoint derivative.
     *
     * @tparam S The system type.
     * @param system The system to integrate; it may refuse any trial or endpoint state.
     * @param tick Index of the step; the state is at `tick * h`.
     * @param h Step size; finite and positive.
     * @param state State at `tick * h`; replaced by the state at `(tick + 1) * h` only on success.
     * @return Nothing on success; otherwise the first fault, with `state` unchanged. InvalidInput for a
     *         bad step, or a dynamic state of the wrong size.
     */
    // PERF: tick times keep the endpoint-derivative cache hit (31% missed).
    template<OdeSystem<N> S>
    [[nodiscard]] std::expected<void, Fault> stepTick(S& system, std::uint64_t tick, sub0math::f64::Time h,
                                                      StateRef state)
    {
        const double dt = h.si();
        return stepBetween(system, static_cast<double>(tick) * dt, dt, static_cast<double>(tick + 1U) * dt, state);
    }

    /** Forgets the cached first stage; call after changing the system's parameters in place. */
    void invalidate() noexcept
    {
        cacheValid_ = false;
    }

private:
    using ConstStateRef = std::conditional_t<Dynamic, std::span<const double>, const StateVector<N>&>;

    // The stages of one step. The fixed size keeps them local to the step, so the compiler can hold them in
    // registers; the dynamic size keeps one set, allocated at construction.
    struct Scratch
    {
        StateVector<N> k2;
        StateVector<N> k3;
        StateVector<N> k4;
        StateVector<N> trial;
        StateVector<N> next;
        StateVector<N> end;
    };

    // PERF: flattened: inlined evaluate 43 vs 52 ns per RK4 step.
    template<typename S>
    SUB0PHYSICS_FLATTEN std::expected<void, Fault> stepBetween(S& system, double t, double h, double end, StateRef state)
    {
        if (!std::isfinite(t) || !std::isfinite(h) || !(h > 0.0) || !std::isfinite(end))
        {
            return std::unexpected(Fault{FaultKind::InvalidInput, "step needs a finite time and positive step"});
        }
        if constexpr (Dynamic)
        {
            if (state.size() != size())
            {
                return std::unexpected(Fault{FaultKind::InvalidInput, "state size differs from the integrator's"});
            }
            return stepWith(system, t, h, end, state, scratch_);
        }
        else
        {
            // PERF: local stages avoid store-forwarding stalls (Midpoint 1.3x to 1.0x).
            Scratch local;
            return stepWith(system, t, h, end, state, local);
        }
    }

    template<typename S>
    std::expected<void, Fault> stepWith(S& system, double t, double h, double end, StateRef state, Scratch& s)
    {
        auto result = advance(system, t, h, end, state, s);
        if (!result)
        {
            // Leave the system observing the last valid state. The trial's fault is the one to report.
            cacheValid_ = system.evaluate(t, state, k1_).has_value();
            cachedTime_ = t;
            assign(cachedState_, state);
            return result;
        }
        assign(state, s.next);
        cachedTime_ = end;
        assign(cachedState_, s.next);
        if constexpr (Dynamic)
        {
            k1_.swap(s.end); // the old first stage is scratch now
        }
        else
        {
            k1_ = s.end;
        }
        cacheValid_ = true;
        return {};
    }

    template<typename S>
    std::expected<void, Fault> advance(S& system, double t, double h, double end, ConstStateRef x, Scratch& s)
    {
        const std::size_t n = size();
        if (!cacheValid_ || cachedTime_ != t || !sameState(x))
        {
            if (auto r = system.evaluate(t, x, k1_); !r)
            {
                return r;
            }
        }
        const double half = 0.5 * h;
        stage(x, k1_, half, s.trial);
        if (auto r = system.evaluate(t + half, s.trial, s.k2); !r)
        {
            return r;
        }
        if constexpr (Method == RungeKuttaMethod::Midpoint)
        {
            for (std::size_t i = 0; i < n; ++i)
            {
                s.next[i] = x[i] + h * s.k2[i];
            }
        }
        else
        {
            stage(x, s.k2, half, s.trial);
            if (auto r = system.evaluate(t + half, s.trial, s.k3); !r)
            {
                return r;
            }
            stage(x, s.k3, h, s.trial);
            if (auto r = system.evaluate(end, s.trial, s.k4); !r)
            {
                return r;
            }
            for (std::size_t i = 0; i < n; ++i)
            {
                s.next[i] = x[i] + h * (k1_[i] + 2.0 * s.k2[i] + 2.0 * s.k3[i] + s.k4[i]) / 6.0;
            }
        }
        for (const double value : s.next)
        {
            if (!std::isfinite(value))
            {
                return std::unexpected(Fault{FaultKind::NonFiniteState, "nonfinite physical state"});
            }
        }
        return system.evaluate(end, s.next, s.end);
    }

    void stage(ConstStateRef x, const StateVector<N>& k, double scale, StateVector<N>& trial) const noexcept
    {
        for (std::size_t i = 0; i < size(); ++i)
        {
            trial[i] = x[i] + scale * k[i];
        }
    }

    [[nodiscard]] bool sameState(ConstStateRef x) const noexcept
    {
        if constexpr (Dynamic)
        {
            return std::equal(x.begin(), x.end(), cachedState_.begin());
        }
        else
        {
            return cachedState_ == x;
        }
    }

    template<typename To, typename From>
    static void assign(To& to, const From& from) noexcept
    {
        if constexpr (Dynamic)
        {
            std::copy(from.begin(), from.end(), to.begin());
        }
        else
        {
            to = from;
        }
    }

    StateVector<N> k1_{};
    std::conditional_t<Dynamic, Scratch, std::tuple<>> scratch_{}; // the dynamic size's stages
    StateVector<N> cachedState_{};
    double cachedTime_ = 0.0;
    bool cacheValid_ = false;
};

/**
 * Classic fourth-order Runge-Kutta; the reference tier.
 *
 * @tparam N Number of state variables, or `std::dynamic_extent`.
 */
template<std::size_t N>
using Rk4 = ExplicitRungeKutta<N, RungeKuttaMethod::ClassicRk4>;

/**
 * Explicit midpoint; half the cost of Rk4, and exact while acceleration is constant.
 *
 * @tparam N Number of state variables, or `std::dynamic_extent`.
 */
template<std::size_t N>
using Midpoint = ExplicitRungeKutta<N, RungeKuttaMethod::Midpoint>;

} // namespace sub0physics
