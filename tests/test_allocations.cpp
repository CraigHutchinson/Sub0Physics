// Counts every heap allocation in this executable, to pin the "no allocation while stepping" contract.
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <doctest.h>
#include <expected>
#include <new>
#include <span>
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
    UnilateralSpringDamper tyre =
        *UnilateralSpringDamper::create({.stiffness = sub0math::newtonsPerMetre(90000.0), .damping = sub0math::newtonSecondsPerMetre(400.0)});

    std::expected<void, Fault> evaluate(double, const StateVector<2>& x, StateVector<2>& dxdt) const
    {
        dxdt[0] = x[1];
        // The generalized force on the penetration 0.36 - x[0] acts on x[0] with the opposite sign.
        dxdt[1] = -tyre.respond(sub0math::metres(0.36 - x[0]), sub0math::metresPerSecond(-x[1])).force.si() / 45.0 - StandardGravity;
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
        REQUIRE(((tick % 2U == 0U) ? rk4.stepTick(system, tick, sub0math::seconds(0.0005), x) : midpoint.stepTick(system, tick, sub0math::seconds(0.0005), x)));
    }
    CHECK(count.since() == 0U);
}

TEST_CASE("Recording, ledgers, statistics and selection never allocate after construction")
{
    ChannelTable<4> table(1000);
    EnergyLedger ledger(EnergyAccount{.kinetic = sub0math::joules(1.0)}, sub0math::joules(1.0));
    ErrorStatistics stats(1.0);
    auto selector = FidelitySelector::create({{{"cheap", 1U, 1.0}, 0.5, 1.0}, {{"full", 1U, 2.0}, 10.0, 10.0}}, 3U, 0U);
    REQUIRE(selector);
    const AllocationCount count;
    for (std::size_t i = 0; i < 1000U; ++i)
    {
        const double v = static_cast<double>(i);
        REQUIRE(table.append({v, v, v, v}));
        ledger.record(EnergyAccount{.kinetic = sub0math::joules(v)});
        ledger.includeScale(sub0math::joules(v));
        stats.add(v * 1e-6);
        const auto d = selector->decide(static_cast<double>(i % 3U), i);
        REQUIRE(selector->commit(d, i));
        (void)evaluateDampedOscillator({sub0math::kilograms(1.0), sub0math::newtonsPerMetre(100.0), sub0math::newtonSecondsPerMetre(1.0), sub0math::metres(0.1), {}}, sub0math::seconds(v * 1e-3));
    }
    CHECK(count.since() == 0U);
}

TEST_CASE("A clearance query never allocates, whatever the witness buffer holds")
{
    auto profile = PolylineProfile::create(test_profiles::metresPoints({{-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.1}, {0.2, 0.1}, {0.2, 0.0}, {1.0, 0.0}}));
    REQUIRE(profile);
    const auto valley = PolylineProfile::create(test_profiles::metresPoints({{-1.0, 1.0}, {0.0, 0.0}, {1.0, 1.0}}));
    REQUIRE(valley);
    std::array<Witness2, DefaultWitnessCapacity> buffer{};
    const AllocationCount count;
    const CircleClearance c = measureCircleClearance(*profile, {{test_profiles::at(-0.1, 0.25), sub0math::metres(0.1)}, sub0math::metres(1e-10)}, buffer);
    CHECK(c.witnessCount == 1U);
    // Two ties into a buffer of one: overflow is a flag, not a reallocation.
    const CircleClearance tied = measureCircleClearance(*valley, {{test_profiles::at(0.0, 0.5), sub0math::metres(0.1)}, sub0math::metres(1e-10)},
                                                       std::span<Witness2>(buffer).first(1));
    CHECK(tied.overflowed);
    CHECK(measureCircleClearance(*profile, {{test_profiles::at(0.9, 1.0), sub0math::metres(0.2)}, sub0math::metres(1e-10)}, buffer).status == ClearanceStatus::OutsideDomain);
    CHECK(count.since() == 0U);
}

TEST_CASE("Prepared profile queries never allocate once the ray origin has been sized")
{
    using sub0math::metres;
    const auto step = PolylineProfile::create(test_profiles::metresPoints({{-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.1}, {0.2, 0.1}, {0.2, 0.0}, {1.0, 0.0}}));
    REQUIRE(step);
    const PreparedProfile prepared(*step);
    RayOrigin origin;
    prepared.bindRayOrigin(test_profiles::at(0.1, 0.4), origin); // sizes the arrays
    std::array<Witness2, DefaultWitnessCapacity> buffer{};
    std::array<sub0math::Point2<sub0math::f64::Length>, 4> contour{{test_profiles::at(-0.3, 0.05), test_profiles::at(0.3, 0.05),
                                                                    test_profiles::at(0.3, -0.05), test_profiles::at(-0.3, -0.05)}};
    const AllocationCount count;
    prepared.bindRayOrigin(test_profiles::at(0.1, 0.4), origin);
    CHECK(prepared.castRay(origin, sub0math::Vec2d{0.0, -1.0}, 1e-10).status == RayStatus::Hit);
    CHECK(prepared.castRay(test_profiles::at(0.1, 0.4), sub0math::Vec2d{0.3, -0.9}, 1e-10).status == RayStatus::Hit);
    CHECK(prepared.sweepCircle(test_profiles::at(-0.5, 0.15), {metres(1.0), metres(0.0)}, metres(0.1)).status == CircleSweepStatus::Contact);
    CHECK(prepared.onSurface(test_profiles::at(0.1, 0.1), 1e-10));
    CHECK(prepared.contourIntrusion(contour).status == ContourStatus::Intrudes);
    CHECK(measureCircleClearance(prepared, {{test_profiles::at(-0.1, 0.25), metres(0.1)}, metres(1e-10)}, buffer).witnessCount == 1U);
    CHECK(count.since() == 0U);
}

TEST_CASE("Convex terrain contact and clearance queries never allocate")
{
    using sub0math::metres;
    const auto step = PolylineProfile::create(test_profiles::metresPoints({{-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.1}, {0.2, 0.1}, {0.2, 0.0}, {1.0, 0.0}}));
    REQUIRE(step);
    const ConvexTerrain terrain(*step);
    std::array<sub0math::Point2<sub0math::f64::Length>, 4> ring{{test_profiles::at(-0.05, 0.105), test_profiles::at(0.05, 0.105),
                                                                 test_profiles::at(0.05, 0.3), test_profiles::at(-0.05, 0.3)}};
    std::array<TerrainContact, 8> buffer{};
    const StepMotionBound motion{test_profiles::at(0.0, 0.2), metres(0.001), 0.01, metres(0.0005)};
    const AllocationCount count;
    CHECK(terrain.ringContacts(ring, metres(0.008), metres(0.001), buffer).count > 0U);
    CHECK(terrain.contacts({{ring[0], ring[1]}, metres(0.008)}, metres(0.0), std::span<TerrainContact>{}).overflowed);
    CHECK(terrain.certifyRingClearance(ring, metres(0.008), motion, metres(0.1)).status == TerrainQueryStatus::Available);
    CHECK(count.since() == 0U);
}

TEST_CASE("The bounded periodic solver allocates nothing once the solution has capacity")
{
    std::vector<double> lower(512, 0.0);
    for (std::size_t i = 200; i < 260; ++i)
    {
        lower[i] = 0.015;
    }
    BoundedPeriodicSolution solution;
    solution.reserve(512);
    const AllocationCount count;
    CHECK(solveBoundedPeriodicTridiagonal(lower, 300000.0 * 0.0044, 30.0 / 0.0044, solution));
    CHECK(solveBoundedPeriodicTridiagonal(lower, 300000.0 * 0.0044, 200.0 / 0.0044, solution));
    CHECK(solution.converged);
    CHECK(count.since() == 0U);
}

TEST_CASE("Warm-started and active-set chain solves allocate nothing once the solution has capacity")
{
    std::vector<double> lower(512, 0.0);
    for (std::size_t i = 200; i < 260; ++i)
    {
        lower[i] = 0.015;
    }
    BoundedPeriodicSolution solution;
    solution.reserveActiveSet(512);
    const AllocationCount count;
    CHECK(solveBoundedPeriodicTridiagonalActiveSet(lower, 300000.0 * 0.0044, 30.0 / 0.0044, solution));
    CHECK(solveBoundedPeriodicTridiagonalActiveSet(lower, 300000.0 * 0.0044, 30.0 / 0.0044, solution, {.warmStart = true}));
    CHECK(solveBoundedPeriodicTridiagonal(lower, 300000.0 * 0.0044, 30.0 / 0.0044, solution, {.warmStart = true}));
    CHECK(solution.converged);
    CHECK(count.since() == 0U);
}

TEST_CASE("Adaptive and linearly implicit steps never allocate")
{
    Contact system;
    EmbeddedRungeKutta<2, DormandPrince54> stepper;
    StateVector<2> x{0.4, 0.0};
    const ErrorTolerance<2> tolerance{{1e-6, 1e-4}, 0.0};
    AdaptiveCounters counters;
    double h = 1e-3;
    double stopped = 0.0;
    const AllocationCount count;
    for (int tick = 0; tick < 100; ++tick)
    {
        REQUIRE(advanceAdaptively<2>(stepper, system, tick * 0.01, (tick + 1) * 0.01, x, h, tolerance, StepSizeController{},
                                     counters, [](double, const StateVector<2>&) -> std::expected<void, Fault> { return {}; },
                                     stopped));
    }
    CHECK(count.since() == 0U);
}

TEST_CASE("Corridor rasterising and grid lookups never allocate")
{
    using namespace sub0physics::terrain;
    const auto grid = GridSpec::create({.finestSpacing = sub0math::f64::Length::fromSi(0.02),
                                        .cellsPerSide = 40,
                                        .registration = Registration::CellCentre,
                                        .anchor = {sub0math::f64::Length{}, sub0math::f64::Length::fromSi(-1.6)},
                                        .pages = {10, 4}});
    REQUIRE(grid);
    std::array<std::uint64_t, 512> words{};
    std::size_t cells = 0;
    const AllocationCount count;
    for (int k = 0; k < 100; ++k)
    {
        const double x = 1.0 + 0.01 * k;
        const auto corridor = CorridorBuilder::create({.from = {sub0math::f64::Length::fromSi(x), sub0math::f64::Length::fromSi(0.1), {}},
                                                       .to = {sub0math::f64::Length::fromSi(x + 0.2), sub0math::f64::Length::fromSi(0.3),
                                                              sub0math::f64::Angle::fromSi(0.2)},
                                                       .footprint = {sub0math::f64::Length::fromSi(0.38), sub0math::f64::Length::fromSi(0.05),
                                                                     sub0math::f64::Length::fromSi(0.1)}});
        REQUIRE(corridor);
        const Raster r = corridor->rasterise(grid->cellGrid(0), words);
        REQUIRE(r.status == RasterStatus::Complete);
        cells += r.cells;
        const TerrainPoint p{{sub0math::f64::Length::fromSi(x), sub0math::f64::Length::fromSi(0.1), sub0math::f64::Length{}}};
        cells += grid->cellOf(p, 0).known ? 1U : 0U;
    }
    CHECK(cells > 0U);
    CHECK(count.since() == 0U);
}

TEST_CASE("Terrain page lookups never allocate")
{
    using namespace sub0physics::terrain;
    constexpr std::size_t count = static_cast<std::size_t>(HeightPage::Nodes) * HeightPage::Nodes;
    std::vector<sub0math::f64::Length> heights(count);
    for (std::size_t n = 0; n < count; ++n) heights[n] = sub0math::metres(1.0 + 0.001 * static_cast<double>(n % HeightPage::Nodes));
    const std::vector<std::uint8_t> known(count, 1);
    const auto page = buildHeightPage(*GridSpec::create(sub0math::metres(0.01), 2), TileKey{0, 0, 0}, SourceRevision{1}, heights, known);
    REQUIRE(page);
    const AllocationCount lookups;
    float sum = 0.0F;
    for (int k = 0; k < 1000; ++k)
    {
        if (const auto z = page->heightAt(sub0math::Length::fromSi(0.002F * static_cast<float>(k)), sub0math::Length::fromSi(1.0F))) sum += z.height.si();
        sum += static_cast<float>(page->blockRange(k % 6, 0, 0).hi);
    }
    CHECK(sum > 0.0F);
    CHECK(lookups.since() == 0U);
}

TEST_CASE("Tile cache, extraction and ray casts never allocate once the buffers exist")
{
    using namespace sub0physics::terrain;
    struct Tile
    {
        int value = 0;
    };
    struct Source final : TileSource<Tile>
    {
        LoadResult load(const TileKey& key, Tile& slot) noexcept override
        {
            slot.value = key.i;
            return {LoadState::Ready, {}};
        }
    } source;
    auto cache = SourceTileCache<Tile>::create({.residentTiles = 4, .pendingLoads = 4, .extraSlots = 2}, &source);
    REQUIRE(cache);

    const auto grid = GridSpec::create({.finestSpacing = sub0math::f64::Length::fromSi(0.01), .pages = {1, 1}});
    REQUIRE(grid);
    constexpr std::size_t nodes = static_cast<std::size_t>(HeightPage::Nodes) * HeightPage::Nodes;
    std::vector<sub0math::f64::Length> heights(nodes);
    for (std::size_t n = 0; n < nodes; ++n) heights[n] = sub0math::metres(0.2 + 0.0004 * static_cast<double>(n % HeightPage::Nodes));
    const std::vector<std::uint8_t> known(nodes, 1);
    const auto page = buildHeightPage(*grid, TileKey{0, 0, 0}, SourceRevision{1}, heights, known);
    REQUIRE(page);
    const std::array<const HeightPage*, 1> pointers{&*page};
    const PageWindow window{{{0, 0}, {0, 0}}, 0, pointers};
    HeightExtract out;
    std::vector<sub0math::Length> xs(256, sub0math::Length::fromSi(0.5F)), ys(256, sub0math::Length::fromSi(0.7F)), hs(256);
    std::vector<std::uint8_t> flags(256);
    const AllocationCount count;
    std::size_t sum = 0;
    for (int k = 0; k < 50; ++k)
    {
        const TileKey key{k % 7, 0, 0};
        if (cache->access(key).state == TileAccess::Started)
        {
            cache->pump([](const LoadEvent&) {});
        }
        HeightRequest request;
        request.cells = {{10, 10}, {90 + k, 60}};
        request.stride = (k % 2 == 0) ? 1U : 4U;
        sum += extractHeights(*grid, window, request, out).deliveredCells;
        const auto hit = castRay(out, {{sub0math::Length::fromSi(0.3F), sub0math::Length::fromSi(0.3F), sub0math::Length::fromSi(1.0F)}, {0.0F, 0.0F, -1.0F}},
                                 sub0math::Length::fromSi(2.0F));
        sum += hit.hit ? 1U : 0U;
        sum += page->heightsAt(xs, ys, hs, flags);
    }
    CHECK(sum > 0U);
    CHECK(count.since() == 0U);
}

TEST_CASE("Travel buffer creation, refusal and evaluation never allocate")
{
    const AllocationCount count;
    double sum = 0.0;
    const auto buffer = TravelBuffer::create({.travel = sub0math::metres(0.065), .zoneWidth = sub0math::metres(0.008), .stiffness = sub0math::newtonsPerMetre(1.0e6), .damping = sub0math::newtonSecondsPerMetre(3000.0)});
    REQUIRE(buffer);
    CHECK_FALSE(TravelBuffer::create({.travel = sub0math::metres(0.01), .zoneWidth = sub0math::metres(0.008), .stiffness = sub0math::newtonsPerMetre(1.0e6), .damping = sub0math::newtonSecondsPerMetre(0.0)}));
    for (int k = 0; k < 1000; ++k)
    {
        sum += buffer->respond(sub0math::metres(0.0001 * static_cast<double>(k % 700)), sub0math::metresPerSecond(0.01 * static_cast<double>(k % 200 - 100))).force.si();
    }
    CHECK(sum != 0.0);
    CHECK(count.since() == 0U);
}
