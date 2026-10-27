#include <array>
#include <cmath>
#include <cstdint>
#include <doctest.h>
#include <vector>

#include "sub0physics/integrator/dyadic_step.hpp"

// The dyadic lattice and its advance (docs/LOD.md, "Levers": the step-size lever on a lattice of the base tick).

using namespace sub0physics;

namespace
{

// A portable generator (splitmix64): the same numbers on every platform, unlike the standard distributions.
std::uint64_t splitmix(std::uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30U)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27U)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31U);
}

double unit(std::uint64_t x)
{
    return static_cast<double>(splitmix(x) >> 11U) * (1.0 / 9007199254740992.0);
}

// A system that reports a pseudo-random error ratio per (position, level) and an event inside chosen fine intervals.
struct Synthetic
{
    DyadicClock clock;
    std::uint64_t position = 0;
    std::uint64_t proposed = 0;
    std::vector<std::uint64_t> events; // a switch at the start of each of these fine units
    std::vector<std::array<std::uint64_t, 2>> steps;

    DyadicAttempt attempt(double h, int level)
    {
        const std::uint64_t length = clock.units(level);
        CHECK(clock.step(level) == h);
        CHECK(position % length == 0);
        proposed = position + length;
        DyadicAttempt a;
        a.errorRatio = 1.1 * unit(position * 131U + static_cast<std::uint64_t>(level + 64));
        for (const std::uint64_t e : events)
        {
            a.event = a.event || (e > position && e <= proposed - 1U);
        }
        return a;
    }

    void accept()
    {
        steps.push_back({position, proposed});
        position = proposed;
    }
};

static_assert(DyadicSystem<Synthetic>);

// x'' = -w^2 x (w = 2) by the trapezoid rule (exact linear solve), with step doubling as the error estimate.
struct Oscillator
{
    double omega = 2.0;
    std::array<double, 2> x{1.0, 0.0};
    std::array<double, 2> next{};
    double tolerance = 1e-6;
    std::uint64_t solves = 0;

    [[nodiscard]] std::array<double, 2> trapezoid(const std::array<double, 2>& s, double h)
    {
        ++solves;
        const double a = 0.25 * omega * omega * h * h;
        // [1, -h/2; w^2 h/2, 1] x1 = [1, h/2; -w^2 h/2, 1] x0
        const double r0 = s[0] + 0.5 * h * s[1];
        const double r1 = s[1] - 0.5 * h * omega * omega * s[0];
        const double det = 1.0 + a;
        return {(r0 + 0.5 * h * r1) / det, (r1 - 0.5 * h * omega * omega * r0) / det};
    }

    DyadicAttempt attempt(double h, int /*level*/)
    {
        const auto coarse = trapezoid(x, h);
        next = trapezoid(trapezoid(x, 0.5 * h), 0.5 * h);
        DyadicAttempt a;
        a.errorRatio = std::max(std::abs(next[0] - coarse[0]), std::abs(next[1] - coarse[1]) / omega) / tolerance;
        return a;
    }

    void accept()
    {
        x = next;
    }
};

} // namespace

TEST_CASE("dyadic clock: exact times and lengths on the lattice")
{
    const DyadicClock clock{.base = 3.125e-5, .subdivision = 6};
    CHECK(clock.unitsPerTick() == 64U);
    CHECK(clock.units(0) == 64U);
    CHECK(clock.units(-6) == 1U);
    CHECK(clock.units(5) == 2048U);
    CHECK(clock.step(0) == 3.125e-5);
    CHECK(clock.step(5) == 1e-3);
    CHECK(clock.time(64U * 32000U) == 1.0); // 32000 base ticks are one second, exactly
}

TEST_CASE("dyadic controller: thresholds instead of powers, clamped")
{
    const DyadicStepController c{.minLevel = -4, .maxLevel = 5, .safety = 0.8, .maxShrink = 3};
    // Order 2: grow when ratio * 8 <= 0.512, keep when ratio <= 1, else drop until ratio / 8^d <= 0.512.
    CHECK(c.next(0, 0.06, 2) == 1);
    CHECK(c.next(0, 0.07, 2) == 0);
    CHECK(c.next(0, 1.0, 2) == 0);
    CHECK(c.next(0, 1.1, 2) == -1);
    CHECK(c.next(0, 5.0, 2) == -2);
    CHECK(c.next(0, 1e9, 2) == -3);
    CHECK(c.next(-4, 1e9, 2) == -4);
    CHECK(c.next(5, 0.0, 2) == 5);
    const DyadicClock clock{.base = 1.0, .subdivision = 4};
    // Alignment and boundaries: from 3 units, only the 1-unit level (-4) is aligned.
    CHECK(c.admissible(clock, 3, 100, 5) == -4);
    CHECK(c.admissible(clock, 4, 100, 5) == -2);
    CHECK(c.admissible(clock, 32, 100, 5) == 1);
    CHECK(c.admissible(clock, 32, 40, 5) == -1); // 8 units fit before the boundary
}

TEST_CASE("dyadic advance: aligned steps, boundaries landed on and never crossed, events bisected to the floor")
{
    Synthetic system;
    system.clock = {.base = 1.0, .subdivision = 6};
    // Switches inside fine units (not on the coarse lattice).
    system.events = {1000U, 1001U, 5000U, 12345U};
    const DyadicStepController controller{.minLevel = -6, .maxLevel = 4};
    // Command ticks: every 3 base ticks, so 4-tick steps never fit across them.
    std::uint64_t position = 0;
    int level = 0;
    DyadicCounters counters;
    for (std::uint64_t boundary = 3U * 64U; boundary <= 300U * 64U; boundary += 3U * 64U)
    {
        const auto stop = advanceDyadically(system, system.clock, controller, 2, position, boundary, level, counters,
                                            [](std::uint64_t) { return true; });
        REQUIRE(stop.has_value());
        CHECK(*stop == DyadicStop::Reached);
        CHECK(position == boundary);
    }
    std::uint64_t at = 0;
    for (const auto& s : system.steps)
    {
        CHECK(s[0] == at);
        at = s[1];
        const std::uint64_t length = s[1] - s[0];
        CHECK((length & (length - 1U)) == 0U); // a power of two
        CHECK(s[0] % length == 0U);            // aligned
        CHECK(s[0] / (3U * 64U) == (s[1] - 1U) / (3U * 64U)); // inside one command interval
        for (const std::uint64_t e : system.events)
        {
            if (e > s[0] && e <= s[1] - 1U)
            {
                CHECK(length == 1U); // an event is only ever inside a floor step
            }
        }
    }
    MESSAGE("accepted " << counters.accepted << ", rejected " << counters.rejected << ", event retries "
                        << counters.eventRetries << ", floor accepts " << counters.floorAccepts);
    CHECK(counters.eventRetries > 0U);
    CHECK(counters.floorAccepts >= system.events.size());
    // A replay takes exactly the same steps.
    Synthetic again;
    again.clock = system.clock;
    again.events = system.events;
    position = 0;
    level = 0;
    DyadicCounters replay;
    for (std::uint64_t boundary = 3U * 64U; boundary <= 300U * 64U; boundary += 3U * 64U)
    {
        REQUIRE(advanceDyadically(again, again.clock, controller, 2, position, boundary, level, replay,
                                  [](std::uint64_t) { return true; })
                    .has_value());
    }
    CHECK(again.steps == system.steps);
}

TEST_CASE("dyadic advance with step doubling meets its tolerance, and cost falls as it loosens")
{
    const DyadicClock clock{.base = 0.01, .subdivision = 8};
    const DyadicStepController controller{.minLevel = -8, .maxLevel = 4};
    double previousError = 0.0;
    std::uint64_t previousSolves = ~std::uint64_t{0};
    for (const double tolerance : {1e-10, 1e-8, 1e-6, 1e-4})
    {
        Oscillator system;
        system.tolerance = tolerance;
        std::uint64_t position = 0;
        int level = 0;
        DyadicCounters counters;
        const std::uint64_t end = 1000U * clock.unitsPerTick(); // 10 s
        REQUIRE(advanceDyadically(system, clock, controller, 2, position, end, level, counters,
                                  [](std::uint64_t) { return true; })
                    .has_value());
        const double t = clock.time(position);
        const double error = std::abs(system.x[0] - std::cos(2.0 * t));
        MESSAGE("tolerance " << tolerance << ": error " << error << ", solves " << system.solves << ", accepted "
                             << counters.accepted << ", rejected " << counters.rejected);
        // Error per step is controlled: the global error stays below the accepted steps times the target.
        CHECK(error < static_cast<double>(counters.accepted) * tolerance);
        CHECK(error >= previousError);
        CHECK(system.solves < previousSolves);
        previousError = error;
        previousSolves = system.solves;
    }
}
