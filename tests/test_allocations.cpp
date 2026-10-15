// Counts every heap allocation in this executable, to pin the "no allocation while stepping" contract.
#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <doctest.h>
#include <expected>
#include <new>
#include <vector>

#include "sub0physics/sub0physics.hpp"
#include "test_profiles.hpp"

namespace
{

std::atomic<std::size_t> allocations{0};

struct AllocationCount
{
    std::size_t start = allocations.load();

    [[nodiscard]] std::size_t since() const noexcept
    {
        return allocations.load() - start;
    }
};

} // namespace

void* operator new(std::size_t size)
{
    ++allocations;
    if (void* p = std::malloc(size == 0U ? 1U : size))
    {
        return p;
    }
    throw std::bad_alloc();
}

void operator delete(void* p) noexcept
{
    std::free(p);
}

void operator delete(void* p, std::size_t) noexcept
{
    std::free(p);
}

using namespace sub0physics;

namespace
{

struct Contact
{
    UnilateralSpringDamper tyre{sub0math::newtonsPerMetre(90000.0), sub0math::newtonSecondsPerMetre(400.0)};

    std::expected<void, Fault> evaluate(double, const StateVector<2>& x, StateVector<2>& dxdt) const
    {
        dxdt[0] = x[1];
        dxdt[1] = tyre.respond(sub0math::metres(0.36 - x[0]), sub0math::metresPerSecond(-x[1])).force.si() / 45.0 - StandardGravity;
        return {};
    }
};

} // namespace

TEST_CASE("Integrator steps never allocate")
{
    Contact system;
    Rk4<2> rk4;
    Midpoint<2> midpoint;
    StateVector<2> x{0.4, 0.0};
    const AllocationCount count;
    for (std::uint64_t tick = 0; tick < 10000U; ++tick)
    {
        REQUIRE(((tick % 2U == 0U) ? rk4.stepTick(system, tick, 0.0005, x) : midpoint.stepTick(system, tick, 0.0005, x)));
    }
    CHECK(count.since() == 0U);
}

TEST_CASE("Recording, ledgers, statistics and selection never allocate after construction")
{
    ChannelTable<4> table(1000);
    EnergyLedger ledger(1.0, 1.0);
    ErrorStatistics stats(1.0);
    auto selector = FidelitySelector::create({{{"cheap", 1U, 1.0}, 0.5, 1.0}, {{"full", 1U, 2.0}, 10.0, 10.0}}, 3U, 0U);
    REQUIRE(selector);
    const AllocationCount count;
    for (std::size_t i = 0; i < 1000U; ++i)
    {
        const double v = static_cast<double>(i);
        REQUIRE(table.append({v, v, v, v}));
        ledger.record(v, 0.0, 0.0);
        ledger.includeScale(v);
        stats.add(v * 1e-6);
        const auto d = selector->decide(static_cast<double>(i % 3U), i);
        REQUIRE(selector->commit(d, i));
        (void)evaluateDampedOscillator({sub0math::kilograms(1.0), sub0math::newtonsPerMetre(100.0), sub0math::newtonSecondsPerMetre(1.0), sub0math::metres(0.1), {}}, sub0math::seconds(v * 1e-3));
    }
    CHECK(count.since() == 0U);
}

TEST_CASE("A clearance query allocates only its witness list")
{
    auto profile = PolylineProfile::create(test_profiles::metresPoints({{-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.1}, {0.2, 0.1}, {0.2, 0.0}, {1.0, 0.0}}));
    REQUIRE(profile);
    const AllocationCount count;
    const CircleClearance c = measureCircleClearance(*profile, {{test_profiles::at(-0.1, 0.25), sub0math::metres(0.1)}, sub0math::metres(1e-10)});
    CHECK(c.witnesses.size() == 1U);
    CHECK(count.since() <= 1U);
    const AllocationCount refused;
    CHECK(measureCircleClearance(*profile, {{test_profiles::at(0.9, 1.0), sub0math::metres(0.2)}, sub0math::metres(1e-10)}).status == ClearanceStatus::OutsideDomain);
    CHECK(refused.since() == 0U);
}
