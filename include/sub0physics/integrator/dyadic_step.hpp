#pragma once

#include <cstdint>
#include <expected>
#include <type_traits>
#include <utility>

#include "sub0physics/fault.hpp"

// Error-controlled steps on a dyadic lattice of a base tick (docs/LOD.md, "Levers"; the consumer's time-base decision
// D4: an error-controlled step is `base * 2^k`, never crosses a command tick, and lands exactly on scheduled ticks).
//
// Time is an integer count of fine units, `base / 2^subdivision` each, so every step boundary is exact and a replay
// does not depend on floating-point accumulation of step sizes. A step of level k is `2^(k + subdivision)` fine units
// and starts only at a multiple of its own length (aligned): a step never crosses a multiple of a longer allowed step,
// so every tick of the coarsest lattice (and any boundary the caller names) is landed on.

namespace sub0physics
{

/** The lattice: a base tick split into `2^subdivision` fine units. */
struct DyadicClock
{
    double base = 1.0;    ///< The base tick [s] (a scenario's command tick); positive.
    int subdivision = 0;  ///< S >= 0: the finest step is `base / 2^S`. At most 30.

    /** @return The fine units in one base tick, `2^S`. */
    [[nodiscard]] constexpr std::uint64_t unitsPerTick() const noexcept
    {
        return std::uint64_t{1} << static_cast<unsigned>(subdivision);
    }

    /**
     * @param level k, with `-subdivision <= k` and `k + subdivision <= 62`.
     * @return The step's length in fine units, `2^(k + S)`.
     */
    [[nodiscard]] constexpr std::uint64_t units(int level) const noexcept
    {
        return std::uint64_t{1} << static_cast<unsigned>(level + subdivision);
    }

    /**
     * @param fine A position in fine units.
     * @return Its time [s]: `fine * (base / 2^S)`, one multiplication by an exact power-of-two scaling of the base, so
     *         it is the same number on every platform.
     */
    [[nodiscard]] constexpr double time(std::uint64_t fine) const noexcept
    {
        return static_cast<double>(fine) * (base / static_cast<double>(unitsPerTick()));
    }

    /**
     * @param level k.
     * @return The step's length [s], `base * 2^k` (exact scaling).
     */
    [[nodiscard]] constexpr double step(int level) const noexcept
    {
        return static_cast<double>(units(level)) * (base / static_cast<double>(unitsPerTick()));
    }
};

/**
 * The step-level controller: the dyadic counterpart of StepSizeController. Levels change by whole factors of two, so
 * the decision compares the error ratio with fixed thresholds and evaluates no power function: deterministic on every
 * platform.
 *
 * With an error estimate of order q (a local error that scales as h^(q+1)), the predicted ratio after changing the step
 * by 2^d is `ratio * 2^(d (q + 1))`; a step is kept below `safety^(q + 1)` of the tolerance.
 */
struct DyadicStepController
{
    int minLevel = -4;    ///< Smallest step, `base * 2^minLevel`: an attempt this small is accepted whatever its error or event (the event-location floor). At least `-subdivision`.
    int maxLevel = 5;     ///< Largest step, `base * 2^maxLevel`.
    double safety = 0.8;  ///< Margin below the tolerance the next step aims for (in (0, 1]).
    int maxShrink = 3;    ///< Most levels dropped after one rejection.

    /**
     * The level to try after an attempt at `level` with error ratio `ratio` (error over tolerance).
     *
     * @param level The level attempted.
     * @param ratio Its error ratio; infinity for a failed attempt.
     * @param order q, the error estimate's order.
     * @return `level + 1` when the doubled step is predicted to meet the target, `level` when the step itself does,
     *         else the smallest drop (at most maxShrink) predicted to meet it; clamped to [minLevel, maxLevel].
     */
    [[nodiscard]] int next(int level, double ratio, int order) const noexcept
    {
        // safety^(q+1) and 2^(q+1), by multiplication.
        double target = 1.0;
        double perLevel = 1.0;
        for (int i = 0; i <= order; ++i)
        {
            target *= safety;
            perLevel *= 2.0;
        }
        int k = level;
        if (ratio * perLevel <= target)
        {
            k = level + 1;
        }
        else if (!(ratio <= 1.0))
        {
            double predicted = ratio;
            int drop = 0;
            while (drop < maxShrink && !(predicted <= target))
            {
                predicted /= perLevel;
                ++drop;
            }
            k = level - (drop == 0 ? 1 : drop);
        }
        return k < minLevel ? minLevel : (k > maxLevel ? maxLevel : k);
    }

    /**
     * The largest level not above `wanted` whose step starts aligned at `position` and ends at or before `boundary`.
     *
     * @param clock The lattice.
     * @param position The start, in fine units.
     * @param boundary The next tick that must not be crossed (a command tick, the end), in fine units; after position.
     * @param wanted The controller's level.
     * @return The level; minLevel when nothing larger fits (the caller's boundaries are on the minLevel lattice).
     */
    [[nodiscard]] int admissible(const DyadicClock& clock, std::uint64_t position, std::uint64_t boundary,
                                 int wanted) const noexcept
    {
        int k = wanted > maxLevel ? maxLevel : wanted;
        while (k > minLevel)
        {
            const std::uint64_t length = clock.units(k);
            if (position % length == 0 && boundary - position >= length)
            {
                break;
            }
            --k;
        }
        return k;
    }
};

/** What one attempt of a dyadic system reports. */
struct DyadicAttempt
{
    /** How the attempt ended. */
    enum class Outcome : std::uint8_t
    {
        Proposed, ///< A step was formed; `errorRatio` and `event` describe it, and `accept` commits it.
        Failed,   ///< The step faulted (`fault`); the system's state is unchanged. Retried smaller, above the floor.
        Stop      ///< The system refuses to go on (a guard over the proposed state): the advance stops, nothing is committed.
    };

    Outcome outcome = Outcome::Proposed;
    double errorRatio = 0.0; ///< The step's error estimate over its tolerance (max norm): at most 1 meets it.
    bool event = false;      ///< The step crosses a switch of the model (a contact onset or release) or could not be certified free of one: bisected down to the floor.
    Fault fault{};           ///< When Failed.
};

/** What a dyadic advance spent. */
struct DyadicCounters
{
    std::uint64_t accepted = 0U;      ///< Accepted steps.
    std::uint64_t rejected = 0U;      ///< Attempts rejected for their error, retried smaller.
    std::uint64_t eventRetries = 0U;  ///< Attempts bisected because they held an event.
    std::uint64_t failures = 0U;      ///< Attempts that faulted and were retried smaller.
    std::uint64_t floorAccepts = 0U;  ///< Attempts accepted at the floor although their error exceeded the target or they held an event (an event located to the floor step).
};

/** How a dyadic advance ended. */
enum class DyadicStop : std::uint8_t
{
    Reached,  ///< Landed on the boundary.
    Refused,  ///< The system answered Stop; `position` is the last accepted state.
    Callback  ///< `onAccepted` asked to stop after an accepted step.
};

/**
 * A system advanced by advanceDyadically. `attempt(h, level)` proposes a step of h seconds from its committed state
 * (with its error estimate and whether it holds an event); `accept()` commits the last proposal.
 */
template<typename S>
concept DyadicSystem = requires(S& system, double h, int level) {
    { system.attempt(h, level) } -> std::same_as<DyadicAttempt>;
    system.accept();
};

/**
 * Advances a system from `position` to exactly `boundary` (fine units) with error-controlled steps on the dyadic
 * lattice: every attempt is `base * 2^k`, aligned, never crosses `boundary`, and the last lands on it.
 *
 * - **Error control.** A proposal whose error ratio exceeds 1 is retried at the level the controller predicts, until
 *   the floor (minLevel), where it is accepted (counted as a floor accept).
 * - **Events.** A proposal holding an event (a switch of the force law inside the step, or a contact it cannot certify
 *   absent) is retried one level smaller: bisection on the lattice. A step without the event is accepted; the one
 *   holding it shrinks until the floor, so the switch is located to the floor step and every other step integrates
 *   one smooth branch.
 * - **Faults** are retried one level smaller; at the floor the fault is returned.
 *
 * Decisions depend only on the attempts' reports, the configuration and the position, so a replay reproduces every
 * step. Nothing allocates.
 *
 * @param system The system.
 * @param clock The lattice.
 * @param controller The level controller.
 * @param order The error estimate's order q (2 for step doubling of a second-order step).
 * @param position In: the start (fine units); out: the last accepted state's position.
 * @param boundary The end (fine units), greater than position, on the minLevel lattice.
 * @param level In: the level to try first; out: the controller's level for the next call.
 * @param counters Accumulates what the advance spent.
 * @param onAccepted Called as `onAccepted(position)` after each accepted step; returns true to go on.
 * @return How it stopped, or the fault of a failed attempt at the floor.
 */
template<typename System, typename OnAccepted>
    requires DyadicSystem<System>
[[nodiscard]] std::expected<DyadicStop, Fault> advanceDyadically(System& system, const DyadicClock& clock,
                                                                 const DyadicStepController& controller, int order,
                                                                 std::uint64_t& position, std::uint64_t boundary,
                                                                 int& level, DyadicCounters& counters,
                                                                 OnAccepted&& onAccepted)
{
    while (position < boundary)
    {
        const int k = controller.admissible(clock, position, boundary, level);
        const DyadicAttempt a = system.attempt(clock.step(k), k);
        const bool floor = k <= controller.minLevel;
        if (a.outcome == DyadicAttempt::Outcome::Stop)
        {
            return DyadicStop::Refused;
        }
        if (a.outcome == DyadicAttempt::Outcome::Failed)
        {
            if (floor)
            {
                return std::unexpected(a.fault);
            }
            ++counters.failures;
            level = k - 1;
            continue;
        }
        if (a.event && !floor)
        {
            ++counters.eventRetries;
            level = k - 1;
            continue;
        }
        if (!(a.errorRatio <= 1.0) && !floor)
        {
            ++counters.rejected;
            level = controller.next(k, a.errorRatio, order);
            continue;
        }
        if (a.event || !(a.errorRatio <= 1.0))
        {
            ++counters.floorAccepts;
        }
        system.accept();
        position += clock.units(k);
        ++counters.accepted;
        // After a located switch the step regrows from the floor, one level per aligned step.
        level = a.event ? controller.minLevel : controller.next(k, a.errorRatio, order);
        if (!onAccepted(position))
        {
            return DyadicStop::Callback;
        }
    }
    return DyadicStop::Reached;
}

} // namespace sub0physics
