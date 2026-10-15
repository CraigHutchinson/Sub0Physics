#include <doctest.h>
#include <limits>
#include <vector>

#include "sub0physics/lod.hpp"

using namespace sub0physics;

namespace
{

constexpr double Infinity = std::numeric_limits<double>::infinity();

FidelitySelector threeTier(std::uint64_t dwell)
{
    auto selector = FidelitySelector::create(
        {
            {{"coarse", 1U, 1.0}, 0.2, 0.3},
            {{"medium", 1U, 2.0}, 0.5, 0.7},
            {{"detailed", 1U, 8.0}, Infinity, Infinity},
        },
        dwell, 0U);
    REQUIRE(selector);
    return *std::move(selector);
}

} // namespace

TEST_CASE("Selector validates its bands")
{
    CHECK_FALSE(FidelitySelector::create({}, 0U, 0U));
    CHECK_FALSE(FidelitySelector::create({{{"a", 1U, 1.0}, 0.5, 0.4}}, 0U, 0U));
    CHECK_FALSE(FidelitySelector::create({{{"a", 1U, 2.0}, 0.1, 0.2}, {{"b", 1U, 1.0}, 0.3, 0.4}}, 0U, 0U));
    CHECK_FALSE(FidelitySelector::create({{{"a", 1U, 1.0}, 0.1, 0.2}}, 0U, 1U));
}

TEST_CASE("Refinement is immediate and picks the cheapest valid model")
{
    FidelitySelector selector = threeTier(100U);
    const FidelityDecision stay = selector.decide(0.25, 1U);
    CHECK(stay.reason == SwitchReason::Stay);

    const FidelityDecision refine = selector.decide(0.6, 2U);
    CHECK(refine.reason == SwitchReason::Refine);
    CHECK(refine.to == 1U);
    REQUIRE(selector.commit(refine, 2U));
    CHECK(selector.active() == 1U);

    const FidelityDecision jump = selector.decide(5.0, 3U);
    CHECK(jump.to == 2U);
}

TEST_CASE("Coarsening waits for the dwell and respects hysteresis")
{
    FidelitySelector selector = threeTier(10U);
    REQUIRE(selector.commit(selector.decide(1.0, 0U), 0U));
    REQUIRE(selector.active() == 2U);

    CHECK(selector.decide(0.1, 5U).reason == SwitchReason::Stay);
    const FidelityDecision coarsen = selector.decide(0.1, 10U);
    CHECK(coarsen.reason == SwitchReason::Coarsen);
    CHECK(coarsen.to == 0U);

    // 0.6 is inside medium's envelope but above its entry threshold, so it cannot be entered.
    CHECK(selector.decide(0.6, 20U).reason == SwitchReason::Stay);
    CHECK(selector.decide(0.45, 20U).to == 1U);
}

TEST_CASE("No valid model keeps the active one and says so")
{
    auto selector = FidelitySelector::create({{{"only", 1U, 1.0}, 0.5, 1.0}}, 0U, 0U);
    REQUIRE(selector);
    const FidelityDecision d = selector->decide(2.0, 1U);
    CHECK(d.reason == SwitchReason::NoValidModel);
    CHECK_FALSE(d.isSwitch());
    CHECK(selector->decide(std::numeric_limits<double>::quiet_NaN(), 1U).reason == SwitchReason::NoValidModel);
}

TEST_CASE("A stale decision cannot be committed")
{
    FidelitySelector selector = threeTier(0U);
    const FidelityDecision first = selector.decide(0.6, 1U);
    const FidelityDecision second = selector.decide(5.0, 1U);
    REQUIRE(selector.commit(first, 1U));
    CHECK_FALSE(selector.commit(second, 1U));
    CHECK(selector.active() == 1U);
    CHECK(selector.switchCount() == 1U);
}
