#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <type_traits>
#include <utility>

#include "sub0physics/fault.hpp"
#include "sub0physics/integrator/ode_system.hpp"

namespace sub0physics
{

/**
 * Per-component error tolerance of an adaptive step: component i may err by
 * `absolute[i] + relative * max(|x0[i]|, |x1[i]|)` per step. An infinite `absolute[i]` leaves the component out of
 * the control (a quadrature that nothing depends on, such as accumulated work).
 *
 * @tparam N Number of state variables.
 */
template<std::size_t N>
struct ErrorTolerance
{
    StateVector<N> absolute{}; ///< Absolute tolerance per component, in that component's unit.
    double relative = 0.0;     ///< Relative tolerance, applied to every controlled component.
};

/**
 * Largest ratio of a step's error estimate to its tolerance (the max norm): at most 1 means the step meets it.
 *
 * @tparam N Number of state variables.
 * @param error The stepper's local error estimate.
 * @param x0 State at the start of the step.
 * @param x1 Proposed state at its end.
 * @param tolerance The tolerance.
 * @return The ratio; infinity when an estimate is not finite.
 */
template<std::size_t N>
[[nodiscard]] double errorRatio(const StateVector<N>& error, const StateVector<N>& x0, const StateVector<N>& x1,
                                const ErrorTolerance<N>& tolerance) noexcept
{
    double ratio = 0.0;
    for (std::size_t i = 0; i < N; ++i)
    {
        const double absolute = tolerance.absolute[i];
        if (absolute == std::numeric_limits<double>::infinity())
        {
            continue;
        }
        const double allowed = absolute + tolerance.relative * std::max(std::abs(x0[i]), std::abs(x1[i]));
        const double r = std::abs(error[i]) / allowed;
        if (!(r <= std::numeric_limits<double>::max()))
        {
            return std::numeric_limits<double>::infinity();
        }
        ratio = std::max(ratio, r);
    }
    return ratio;
}

/**
 * Deterministic step-size controller (elementary, with safety factor and limits): the next step is
 * `h * clamp(safety * ratio^(-1/(q+1)), minFactor, maxFactor)`, then clamped to `[minStep, maxStep]`, with `q` the
 * error estimate's order. It depends only on the ratio and the step, never on wall-clock time, so a replay
 * reproduces every step.
 */
struct StepSizeController
{
    double minStep = 1e-7;   ///< Smallest step, s; a step this small is accepted whatever its error.
    double maxStep = 1e-2;   ///< Largest step, s.
    double safety = 0.9;     ///< Factor below the predicted optimal step.
    double minFactor = 0.2;  ///< Largest shrink per step.
    double maxFactor = 5.0;  ///< Largest growth per step.

    /**
     * @param h The step just attempted, s.
     * @param ratio Its error ratio (errorRatio).
     * @param errorOrder The error estimate's order q.
     * @return The next step to try, s.
     */
    [[nodiscard]] double next(double h, double ratio, int errorOrder) const noexcept
    {
        double factor = maxFactor;
        if (ratio > 0.0)
        {
            factor = ratio == std::numeric_limits<double>::infinity()
                         ? minFactor
                         : std::clamp(safety * std::pow(ratio, -1.0 / static_cast<double>(errorOrder + 1)), minFactor,
                                      maxFactor);
        }
        return std::clamp(h * factor, minStep, maxStep);
    }
};

/** What an adaptive advance spent. */
struct AdaptiveCounters
{
    std::uint64_t accepted = 0U;      ///< Accepted steps.
    std::uint64_t rejected = 0U;      ///< Attempts rejected for error or refused by the system, then retried smaller.
    std::uint64_t eventRetries = 0U;  ///< Attempts shortened to end on an event (a sign change of the event function).
};

/** The event function of an advance without events. */
struct NoEvent
{
};

namespace detail
{

// What an event function returns: one value, or an array of them (several switches).
template<typename Event, std::size_t N>
struct EventValueOf
{
    using type = std::remove_cvref_t<decltype(std::declval<Event&>()(0.0, std::declval<const StateVector<N>&>()))>;
};

template<std::size_t N>
struct EventValueOf<NoEvent, N>
{
    using type = double;
};

// Fraction of the step at which the value changes sign (secant), or -1 without a change.
[[nodiscard]] inline double crossingFraction(double start, double end) noexcept
{
    if ((start > 0.0) == (end > 0.0) || start == end)
    {
        return -1.0;
    }
    return std::clamp(start / (start - end), 0.0, 1.0);
}

// The earliest sign change over several event functions.
template<std::size_t M>
[[nodiscard]] double crossingFraction(const std::array<double, M>& start, const std::array<double, M>& end) noexcept
{
    double earliest = -1.0;
    for (std::size_t i = 0; i < M; ++i)
    {
        const double f = crossingFraction(start[i], end[i]);
        if (f >= 0.0 && (earliest < 0.0 || f < earliest))
        {
            earliest = f;
        }
    }
    return earliest;
}

} // namespace detail

/**
 * A stepper that proposes a step with an error estimate and commits it on request: EmbeddedRungeKutta and
 * LinearlyImplicitEuler.
 */
template<typename T, typename S, std::size_t N>
concept ErrorEstimatingStepper = requires(T& stepper, S& system, double t, double h, StateVector<N>& x) {
    { stepper.attempt(system, t, h, std::as_const(x)) } -> std::same_as<std::expected<void, Fault>>;
    { stepper.proposal() } -> std::same_as<const StateVector<N>&>;
    { stepper.error() } -> std::same_as<const StateVector<N>&>;
    stepper.accept(x);
    { T::ErrorOrder } -> std::convertible_to<int>;
};

/**
 * Advances `x` from `t0` to exactly `t1` with error-controlled sub-steps: the continuous step-size lever driven to
 * a tolerance target (docs/LOD.md, "Levers").
 *
 * The outer interval is a fixed tick of the caller (a recording interval or a frame), so the caller's grid and
 * replay are unchanged; inside it, steps shrink where the estimate is large and grow where it is small. The last
 * sub-step is shortened to land on `t1`; `h` carries the controller's step to the next tick. A system refusal (a
 * fault) on an attempt is treated as an unbounded error and retried smaller, until `minStep`, where the fault is
 * returned and the system observes `(t, x)` again.
 *
 * **Events.** A local error estimate cannot see a switch in the model (a contact closing or releasing, a damper
 * changing branch): the derivative has a kink there, and a step across it loses the method's order however small
 * its estimate. With an event function `g(t, x)` whose sign marks the branch, a step whose end has the other sign
 * is retried to end just past the zero (secant on the step length), until the zero lies within `eventTolerance` of
 * the step's end, so every step integrates one smooth branch and the switch is located to that time. An event
 * function may return a `std::array<double, M>` for several switches; the earliest crossing in a step decides.
 *
 * Decisions depend only on the state, the configuration and the times, so the result is deterministic.
 *
 * @tparam N Number of state variables.
 * @param stepper The stepper.
 * @param system The system.
 * @param t0 Start time, s.
 * @param t1 End time, s; greater than t0.
 * @param x State at `t0`; the state at `t1` on success, the last accepted state on failure.
 * @param h In: the step to try first; out: the controller's step for the next call. Positive.
 * @param tolerance Error tolerance per step.
 * @param controller Step-size rules.
 * @param counters Accumulates accepted, rejected and event-shortened attempts.
 * @param onAccepted Called as `onAccepted(t, x)` after each accepted sub-step (`std::expected<void, Fault>`); a
 *        fault stops the advance and is returned.
 * @param stoppedAt Receives the time of the last accepted state.
 * @param event `double event(double t, const StateVector<N>& x)` (or returning `std::array<double, M>`), or NoEvent.
 * @param eventTolerance Time within which an event is located, s.
 * @return Nothing, or the first fault.
 */
template<std::size_t N, typename Stepper, typename S, typename OnAccepted, typename Event = NoEvent>
    requires ErrorEstimatingStepper<Stepper, S, N>
[[nodiscard]] std::expected<void, Fault> advanceAdaptively(Stepper& stepper, S& system, double t0, double t1,
                                                           StateVector<N>& x, double& h,
                                                           const ErrorTolerance<N>& tolerance,
                                                           const StepSizeController& controller,
                                                           AdaptiveCounters& counters, OnAccepted&& onAccepted,
                                                           double& stoppedAt, Event&& event = Event{},
                                                           double eventTolerance = 1e-6)
{
    constexpr bool HasEvent = !std::is_same_v<std::remove_cvref_t<Event>, NoEvent>;
    double t = t0;
    stoppedAt = t0;
    using EventValue = typename detail::EventValueOf<std::remove_cvref_t<Event>, N>::type;
    EventValue eventAtStart{};
    if constexpr (HasEvent)
    {
        eventAtStart = event(t, std::as_const(x));
    }
    double eventStep = 0.0; // when positive, the next attempt's length, chosen to end just past an event
    while (t < t1)
    {
        const double wanted = eventStep > 0.0 ? eventStep : h;
        const bool shortenedForEvent = eventStep > 0.0;
        eventStep = 0.0;
        // A remainder shorter than a millionth of the step joins it, so no sliver step follows.
        const bool reachesEnd = t + wanted >= t1 || t1 - (t + wanted) <= 1e-6 * wanted;
        const double end = reachesEnd ? t1 : t + wanted;
        const double step = end - t;
        const auto attempted = stepper.attempt(system, t, step, std::as_const(x));
        double ratio = std::numeric_limits<double>::infinity();
        if (attempted)
        {
            ratio = errorRatio(stepper.error(), x, stepper.proposal(), tolerance);
        }
        // The requested step, since end - t can round just above the floor.
        const bool atMinimum = wanted <= controller.minStep || step <= controller.minStep;
        if (!attempted && atMinimum)
        {
            StateVector<N> scratch{};
            (void)system.evaluate(t, x, scratch); // observe the last accepted state again
            return std::unexpected(attempted.error());
        }
        if (attempted && (ratio <= 1.0 || atMinimum))
        {
            EventValue eventAtEnd{};
            if constexpr (HasEvent)
            {
                eventAtEnd = event(end, stepper.proposal());
                const double fraction = detail::crossingFraction(eventAtStart, eventAtEnd);
                if (fraction >= 0.0 && step > eventTolerance)
                {
                    const double crossing = step * fraction;
                    const double retry = crossing + 0.5 * eventTolerance;
                    if (step - crossing > eventTolerance && retry < step)
                    {
                        ++counters.eventRetries;
                        eventStep = std::max(retry, controller.minStep);
                        continue; // the step-size controller's h is untouched
                    }
                }
            }
            stepper.accept(x);
            t = end;
            stoppedAt = t;
            eventAtStart = eventAtEnd;
            ++counters.accepted;
            const double proposed = controller.next(step, ratio, Stepper::ErrorOrder);
            // A step shortened to land on t1 or on an event says nothing about longer ones: keep h unless it must shrink.
            h = (reachesEnd || shortenedForEvent) && proposed >= step ? h : proposed;
            if (auto r = onAccepted(t, std::as_const(x)); !r)
            {
                return r;
            }
            continue;
        }
        ++counters.rejected;
        h = controller.next(step, ratio, Stepper::ErrorOrder); // below step: the ratio exceeds 1
    }
    return {};
}

} // namespace sub0physics
