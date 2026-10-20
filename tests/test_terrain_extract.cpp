// Source tile cache, column extraction and height extraction, checked against the tier-C analytic surfaces and
// against brute-force or hand-written controls.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <doctest.h>
#include <map>
#include <optional>
#include <random>
#include <vector>

#include "sub0physics/terrain.hpp"
#include "sub0physics/terrain/testing.hpp"

using namespace sub0physics;
using namespace sub0physics::terrain;
namespace f64 = sub0math::f64;
using f64::Length;
using LengthF = sub0math::Length;

namespace
{

// ---- Cache ---------------------------------------------------------------------------------------------------

struct IntTile
{
    int value = -1;
    std::int32_t i = -1;
};

class ScriptedSource final : public TileSource<IntTile>
{
public:
    LoadResult load(const TileKey& key, IntTile& slot) noexcept override
    {
        ++asked;
        if (failing.contains(key.i))
        {
            return {LoadState::Failed, "boom"};
        }
        if (!open)
        {
            return {LoadState::Waiting, {}};
        }
        slot.value = key.i * 10 + key.j;
        slot.i = key.i;
        return {LoadState::Ready, {}};
    }

    bool open = false;
    int asked = 0;
    std::map<std::int32_t, bool> failing;
};

std::vector<TileKey> residentKeys(const SourceTileCache<IntTile>& cache)
{
    std::array<TileKey, 64> buffer{};
    const auto n = cache.resident(buffer);
    return {buffer.begin(), buffer.begin() + static_cast<std::ptrdiff_t>(n)};
}

// ---- Surfaces ------------------------------------------------------------------------------------------------

const GridSpec& pageGrid()
{
    static const GridSpec grid = *GridSpec::create({.finestSpacing = Length::fromSi(0.01), .pages = {2, 2}});
    return grid;
}

template<typename S>
HeightPage bakePage(const S& surface, TileKey key, const std::vector<std::pair<int, int>>& unknown = {})
{
    const auto n = static_cast<std::size_t>(HeightPage::Nodes);
    std::vector<Length> heights(n * n);
    std::vector<std::uint8_t> known(n * n);
    testing::sampleSurface(pageGrid(), key, surface, heights, known);
    for (const auto& [i, j] : unknown)
    {
        known[static_cast<std::size_t>(j) * n + static_cast<std::size_t>(i)] = 0U;
    }
    auto page = buildHeightPage(pageGrid(), key, SourceRevision{1}, heights, known);
    REQUIRE(page);
    return std::move(*page);
}

// A 2 x 2 window of pages of one surface.
struct Quad
{
    std::vector<HeightPage> pages;
    std::vector<const HeightPage*> pointers;

    template<typename S>
    explicit Quad(const S& surface, const std::vector<std::pair<int, int>>& unknownInFirst = {})
    {
        pages.reserve(4);
        for (std::int32_t j = 0; j < 2; ++j)
        {
            for (std::int32_t i = 0; i < 2; ++i)
            {
                pages.push_back(bakePage(surface, {i, j, 0}, i == 0 && j == 0 ? unknownInFirst : std::vector<std::pair<int, int>>{}));
            }
        }
        for (const HeightPage& p : pages)
        {
            pointers.push_back(&p);
        }
    }

    [[nodiscard]] PageWindow window() const
    {
        return {{{0, 0}, {1, 1}}, 0, pointers};
    }
};

double terrainZ(const HeightExtract& e, std::size_t node)
{
    return e.origin().z.si() + static_cast<double>(e.z()[node].si());
}

} // namespace

TEST_CASE("Cache: de-duplicated loads, hits and misses, key-ordered residency")
{
    ScriptedSource source;
    auto cache = SourceTileCache<IntTile>::create({.residentTiles = 3, .pendingLoads = 4, .extraSlots = 2, .tileBytes = 100}, &source);
    REQUIRE(cache);
    CHECK(cache->access({1, 0, 0}).state == TileAccess::Started);
    CHECK(cache->access({1, 0, 0}).state == TileAccess::Loading);
    CHECK(cache->access({0, 2, 0}).state == TileAccess::Started);
    CHECK(cache->stats() == CacheStats{1, 2, 0, 0});
    CHECK(source.asked == 0); // nothing is asked outside pump
    CHECK_FALSE(cache->pump([](const LoadEvent&) {})); // not available yet
    CHECK(source.asked == 2);
    source.open = true;
    std::vector<TileKey> loaded;
    CHECK(cache->pump([&](const LoadEvent& e) { loaded.push_back(e.key); }));
    CHECK(loaded == std::vector<TileKey>{{1, 0, 0}, {0, 2, 0}}); // oldest first
    CHECK(cache->bytes() == 200U);
    CHECK(residentKeys(*cache) == std::vector<TileKey>{{0, 2, 0}, {1, 0, 0}});
    const auto hit = cache->access({1, 0, 0});
    REQUIRE(hit.state == TileAccess::Resident);
    CHECK(cache->tile(hit.slot).value == 10);
    CHECK(cache->stats() == CacheStats{2, 2, 2, 0});
}

TEST_CASE("Cache: least recently accessed tile is evicted, pins keep storage valid, budget is visible")
{
    ScriptedSource source;
    source.open = true;
    auto cache = SourceTileCache<IntTile>::create({.residentTiles = 2, .pendingLoads = 2, .extraSlots = 1}, &source);
    REQUIRE(cache);
    const auto ignore = [](const LoadEvent&) {};
    for (std::int32_t i = 0; i < 2; ++i)
    {
        (void)cache->access({i, 0, 0});
    }
    cache->pump(ignore);
    // Touch tile 0 so tile 1 is the oldest.
    const auto zero = cache->access({0, 0, 0});
    REQUIRE(zero.state == TileAccess::Resident);
    auto one = cache->access({1, 0, 0});
    REQUIRE(one.state == TileAccess::Resident);
    auto lease = SourceTileCache<IntTile>::Lease(*cache, one.slot);
    (void)cache->access({0, 0, 0}); // 1 is now least recently used
    CHECK(cache->access({2, 0, 0}).state == TileAccess::Started);
    CHECK(cache->access({3, 0, 0}).state == TileAccess::Started);
    CHECK(cache->access({4, 0, 0}).state == TileAccess::Refused); // two loads in flight
    cache->pump(ignore);
    CHECK(residentKeys(*cache) == std::vector<TileKey>{{2, 0, 0}, {3, 0, 0}}); // 0 and 1 evicted, 0 least recently used
    CHECK(cache->stats().evictions == 2U);
    CHECK(lease.tile().value == 10); // evicted but pinned: storage intact
    CHECK(cache->residentSlot({1, 0, 0}) < 0);
    lease.release();
}

TEST_CASE("Cache: failed loads are reported once and leave no trace; events arrive before the next load")
{
    ScriptedSource source;
    source.open = true;
    source.failing[7] = true;
    auto cache = SourceTileCache<IntTile>::create({.residentTiles = 1, .pendingLoads = 4, .extraSlots = 1}, &source);
    REQUIRE(cache);
    (void)cache->access({7, 0, 0});
    (void)cache->access({1, 0, 0});
    (void)cache->access({2, 0, 0});
    std::vector<std::pair<std::int32_t, bool>> events;
    std::int32_t pinned = -1;
    cache->pump([&](const LoadEvent& e) {
        events.emplace_back(e.key.i, e.failed);
        if (!e.failed && e.key.i == 1)
        {
            pinned = e.slot;
            cache->pin(e.slot); // pinned before the next load can evict it
        }
    });
    CHECK(events == std::vector<std::pair<std::int32_t, bool>>{{7, true}, {1, false}, {2, false}});
    CHECK(cache->pendingCount() == 0U);
    CHECK(cache->tile(pinned).i == 1);
    CHECK(residentKeys(*cache) == std::vector<TileKey>{{2, 0, 0}});
    cache->unpin(pinned);
    CHECK(cache->stats().loads == 2U);
    CHECK_FALSE(SourceTileCache<IntTile>::create({.residentTiles = 0}, &source));
    CHECK_FALSE(SourceTileCache<IntTile>::create({}, nullptr));
}

TEST_CASE("Cache: a random access sequence matches a straightforward LRU model")
{
    ScriptedSource source;
    source.open = true;
    constexpr std::uint32_t Limit = 5;
    auto cache = SourceTileCache<IntTile>::create({.residentTiles = Limit, .pendingLoads = 4}, &source);
    REQUIRE(cache);
    std::vector<std::int32_t> model; // oldest first
    std::mt19937 rng(7);
    std::uint64_t evictions = 0;
    for (int step = 0; step < 2000; ++step)
    {
        const auto i = static_cast<std::int32_t>(rng() % 9U);
        const auto r = cache->access({i, 0, 0});
        const auto it = std::ranges::find(model, i);
        if (it != model.end())
        {
            REQUIRE(r.state == TileAccess::Resident);
            model.erase(it);
            model.push_back(i);
        }
        else
        {
            REQUIRE(r.state == TileAccess::Started);
            cache->pump([](const LoadEvent&) {});
            if (model.size() == Limit)
            {
                model.erase(model.begin());
                ++evictions;
            }
            model.push_back(i);
        }
    }
    std::vector<std::int32_t> actual;
    for (const TileKey& k : residentKeys(*cache))
    {
        actual.push_back(k.i);
    }
    std::ranges::sort(model);
    CHECK(actual == model);
    CHECK(cache->stats().evictions == evictions);
}

// ---- Column extraction ---------------------------------------------------------------------------------------

namespace
{

// A cell-centre tile of a catalogue surface, with a second span under the bridge deck.
struct SurfaceColumns
{
    std::int32_t side = 0;
    std::vector<double> heights, lower2;
    std::vector<sub0math::Vec3d> normals;
    std::vector<std::uint8_t> valid;
    std::vector<std::uint8_t> twoSpans;

    [[nodiscard]] bool known(std::size_t c) const { return valid[c] != 0U; }
    [[nodiscard]] Length top(std::size_t c) const { return Length::fromSi(heights[c]); }
    [[nodiscard]] sub0math::Vec3d normal(std::size_t c) const { return normals[c]; }
    [[nodiscard]] std::size_t spanCount(std::size_t c) const { return twoSpans[c] != 0U ? 2U : 1U; }
    [[nodiscard]] ColumnSpan span(std::size_t c, std::size_t k) const
    {
        if (twoSpans[c] == 0U)
        {
            return {Length::fromSi(heights[c]), Length::fromSi(-1.0), 1};
        }
        return k == 0 ? ColumnSpan{Length::fromSi(heights[c]), Length::fromSi(lower2[c]), 2}
                      : ColumnSpan{Length::fromSi(lower2[c]), Length::fromSi(-1.0), 1};
    }
};

using TopAdapter = SurfaceColumns;
static_assert(ColumnTile<SurfaceColumns>);

TopAdapter bakeColumns(const GridSpec& grid, const TileKey& key, const testing::Bridge& bridge,
                       std::optional<std::size_t> unknownCell = std::nullopt)
{
    TopAdapter t;
    const std::int32_t n = grid.cellsPerSide();
    t.side = n;
    for (std::int32_t j = 0; j < n; ++j)
    {
        for (std::int32_t i = 0; i < n; ++i)
        {
            const auto p = grid.samplePosition(key, i, j);
            const double x = p.x.si(), y = p.y.si();
            t.heights.push_back(bridge.height(x, y));
            t.normals.push_back(bridge.normal(x, y));
            const auto spans = bridge.spansAt(x, y);
            t.twoSpans.push_back(spans.count == 2U ? 1U : 0U);
            t.lower2.push_back(spans.spans[1].lower);
            t.valid.push_back(1U);
        }
    }
    if (unknownCell)
    {
        t.valid[*unknownCell] = 0U;
    }
    return t;
}

} // namespace

TEST_CASE("Column extract: a cell-centre block across page seams reproduces the analytic columns")
{
    const auto grid = *GridSpec::create({.finestSpacing = Length::fromSi(0.02),
                                         .cellsPerSide = 40,
                                         .registration = Registration::CellCentre,
                                         .anchor = {Length::fromSi(0.0), Length::fromSi(-1.6)},
                                         .pages = {10, 4}});
    const testing::Bridge bridge{.x0 = 1.1, .x1 = 1.5, .deckTop = 0.6, .thickness = 0.1};
    std::map<std::pair<int, int>, TopAdapter> tiles;
    for (int j = 0; j < 4; ++j)
    {
        for (int i = 0; i < 10; ++i)
        {
            tiles.emplace(std::pair{i, j}, bakeColumns(grid, {i, j, 0}, bridge));
        }
    }
    const auto lookup = [&](const TileKey& k) -> const TopAdapter* { return &tiles.at({k.i, k.j}); };
    ColumnExtract out;
    const sub0math::Aabb<Length, 2> box{{Length::fromSi(0.9), Length::fromSi(-0.3)}, {Length::fromSi(1.7), Length::fromSi(0.3)}};
    const auto cells = grid.cellsTouching(box, 0);
    const auto r = extractColumns<TopAdapter>(grid, 0, cells, lookup, out);
    REQUIRE(r.status == ColumnStatus::Complete);
    REQUIRE(out.ids().size() == static_cast<std::size_t>(cells.count()));
    REQUIRE(out.offset().size() == out.ids().size() + 1U);
    const std::uint64_t across = 400;
    std::size_t spans = 0;
    for (std::size_t c = 0; c < out.ids().size(); ++c)
    {
        const auto i = static_cast<double>(out.ids()[c] % across), j = static_cast<double>(out.ids()[c] / across);
        CHECK(out.x()[c].si() == (i + 0.5) * 0.02);
        CHECK(out.y()[c].si() == -1.6 + (j + 0.5) * 0.02);
        CHECK(out.top()[c].si() == bridge.height(out.x()[c].si(), out.y()[c].si()));
        const auto expect = bridge.spansAt(out.x()[c].si(), out.y()[c].si());
        const std::size_t count = out.offset()[c + 1] - out.offset()[c];
        REQUIRE(count == expect.count);
        spans += count;
    }
    CHECK(out.upper().size() == spans);
    CHECK(out.material().size() == spans);
    // Budget: a block over the column capacity is refused, not truncated.
    ColumnExtract small(ColumnLimits{.columns = 10, .spans = 20});
    CHECK(extractColumns<TopAdapter>(grid, 0, cells, lookup, small).status == ColumnStatus::ColumnBudget);
    CHECK(small.ids().empty());
    ColumnExtract fewSpans(ColumnLimits{.columns = 4096, .spans = 5});
    CHECK(extractColumns<TopAdapter>(grid, 0, cells, lookup, fewSpans).status == ColumnStatus::SpanBudget);
}

TEST_CASE("Column extract: unauthored cells make the result Partial and are never delivered as flat")
{
    const auto grid = *GridSpec::create({.finestSpacing = Length::fromSi(0.02),
                                         .cellsPerSide = 40,
                                         .registration = Registration::CellCentre,
                                         .pages = {2, 2}});
    const testing::Bridge bridge{};
    TopAdapter holey = bakeColumns(grid, {0, 0, 0}, bridge, 5U * 40U + 7U);
    TopAdapter full = bakeColumns(grid, {1, 0, 0}, bridge);
    const auto lookup = [&](const TileKey& k) -> const TopAdapter* { return k.i == 0 ? &holey : (k.j == 0 ? &full : nullptr); };
    ColumnExtract out;
    const sub0math::GridRange<2> block{{2, 3}, {20, 9}};
    const auto r = extractColumns<TopAdapter>(grid, 0, block, lookup, out);
    REQUIRE(r.status == ColumnStatus::Partial);
    REQUIRE(out.missing().size() == 1U);
    CHECK(out.missing()[0] == 5U * 80U + 7U);
    CHECK(out.ids().empty());
    const sub0math::GridRange<2> needsAbsentPage{{70, 70}, {75, 75}};
    CHECK(extractColumns<TopAdapter>(grid, 0, needsAbsentPage, lookup, out).status == ColumnStatus::MissingTile);
}

// ---- Height extraction ---------------------------------------------------------------------------------------

TEST_CASE("Height extract: a plane across four pages is exact and the seams agree")
{
    const testing::Plane plane{0.30, 0.10, -0.05};
    Quad quad(plane);
    HeightExtract out;
    HeightRequest request;
    request.cells = {{230, 240}, {290, 300}}; // straddles both seams
    const auto r = extractHeights(pageGrid(), quad.window(), request, out);
    REQUIRE(r.status == HeightStatus::Complete);
    CHECK(r.deliveredCells == 61U * 61U);
    CHECK(out.triangles().size() == 2U * r.deliveredCells);
    const auto [cx, cy] = out.cells();
    REQUIRE(cx == 61);
    const double h = 0.01;
    for (std::int32_t b = 0; b <= cy; ++b)
    {
        for (std::int32_t a = 0; a <= cx; ++a)
        {
            const auto node = static_cast<std::size_t>(b * (cx + 1) + a);
            const double x = pageGrid().anchor().x.si() + h * (230 + a), y = pageGrid().anchor().y.si() + h * (240 + b);
            CHECK(out.nodeKnown()[node] == 1U);
            CHECK(std::abs(terrainZ(out, node) - plane.height(x, y)) < 2e-4); // 0.1 mm steps
        }
    }
    // Triangles: every vertex height is the node height, and the normal is the plane's.
    const auto n = plane.normal(0, 0);
    for (std::size_t t = 0; t < out.triangles().size(); ++t)
    {
        CHECK(out.triangleNormal()[t].x == doctest::Approx(static_cast<float>(n.x)).epsilon(5e-3));
        CHECK(out.triangleNormal()[t].z == doctest::Approx(static_cast<float>(n.z)).epsilon(5e-3));
    }
    for (std::size_t c = 0; c < out.cellState().size(); ++c)
    {
        REQUIRE(out.cellState()[c] == CellState::Delivered);
        const auto t0 = out.triangles()[out.cellTriangle()[c]];
        const auto a = static_cast<std::int32_t>(c) % cx, b = static_cast<std::int32_t>(c) / cx;
        const auto node = [&](std::int32_t u, std::int32_t v) { return static_cast<std::size_t>((b + v) * (cx + 1) + a + u); };
        CHECK(t0.a.z.si() == out.z()[node(0, 0)].si());
        CHECK(t0.b.z.si() == out.z()[node(1, 0)].si());
        CHECK(t0.c.z.si() == out.z()[node(1, 1)].si());
    }
}

TEST_CASE("Height extract: unknown nodes are never interpolated across, and a mask limits what is asked")
{
    Quad quad(testing::Plane{0.2, 0.0, 0.0}, {{10, 10}});
    HeightExtract out;
    HeightRequest request;
    request.cells = {{5, 5}, {20, 20}};
    const auto r = extractHeights(pageGrid(), quad.window(), request, out);
    REQUIRE(r.status == HeightStatus::Partial);
    CHECK(r.unknownCells == 4U); // the four cells around node (10, 10)
    CHECK(r.deliveredCells == 16U * 16U - 4U);
    CHECK(out.cellState()[(10 - 5) * 16 + (10 - 5)] == CellState::Unknown);
    CHECK(out.cellState()[(9 - 5) * 16 + (9 - 5)] == CellState::Unknown);
    CHECK(out.cellState()[(8 - 5) * 16 + (8 - 5)] == CellState::Delivered);
    CHECK(out.cellTriangle()[(10 - 5) * 16 + (10 - 5)] == HeightExtract::NoTriangle);
    CHECK(out.nodeKnown()[(10 - 5) * 17 + (10 - 5)] == 0U);

    // A mask that asks for a cell range away from the hole.
    std::array<std::uint64_t, 8> words{};
    CellMask mask{{{12, 12}, {14, 13}}, words};
    words[0] = 0b111111U;
    request.mask = &mask;
    const auto m = extractHeights(pageGrid(), quad.window(), request, out);
    CHECK(m.status == HeightStatus::Complete);
    CHECK(m.deliveredCells == 6U);
    CHECK(out.cellState()[0] == CellState::NotRequested);
}

TEST_CASE("Height extract: a missing page and the grid edge are Unknown; budgets and bad requests are refused")
{
    Quad quad(testing::Plane{});
    quad.pointers[1] = nullptr; // page (1, 0) absent
    HeightExtract out;
    HeightRequest request;
    request.cells = {{250, 10}, {260, 20}};
    const auto r = extractHeights(pageGrid(), quad.window(), request, out);
    CHECK(r.status == HeightStatus::Partial);
    CHECK(r.unknownCells == 11U * 5U); // cells 256..260 are on the absent page
    request.cells = {{0, 0}, {511, 511}};
    HeightExtract tiny(HeightExtractLimits{.nodes = 100, .triangles = 100});
    CHECK(extractHeights(pageGrid(), quad.window(), request, tiny).status == HeightStatus::Budget);
    request.stride = 3;
    CHECK(extractHeights(pageGrid(), quad.window(), request, out).status == HeightStatus::InvalidRequest);
    request.stride = 1;
    request.cells = {{0, 0}, {512, 3}};
    CHECK(extractHeights(pageGrid(), quad.window(), request, out).status == HeightStatus::InvalidRequest);
}

TEST_CASE("Height extract: a step is a ramp across its cell, never a riser; the achieved spacing is reported")
{
    const testing::Step step{0.6, 0.1};
    Quad quad(step);
    HeightExtract out;
    HeightRequest request;
    request.cells = {{55, 20}, {66, 22}};
    REQUIRE(extractHeights(pageGrid(), quad.window(), request, out).status == HeightStatus::Complete);
    // The cell holding the riser interpolates between the two node heights: a ramp, which is why the riser is a
    // feature (docs/TERRAIN.md), and the extract's normal there is tilted, not horizontal.
    const std::size_t riserCell = 4U; // node 59 (x = 0.59) is 0, node 60 (x = 0.60) is the rise: cell 59 -> a = 4
    const auto t = out.cellTriangle()[riserCell];
    CHECK(std::abs(out.triangleNormal()[t].x) > 0.9F);
    CHECK(out.spacing().si() == doctest::Approx(0.01F));
    request.stride = 4;
    request.cells = {{56, 20}, {67, 23}};
    HeightExtract coarse;
    REQUIRE(extractHeights(pageGrid(), quad.window(), request, coarse).status == HeightStatus::Complete);
    CHECK(coarse.spacing().si() == doctest::Approx(0.04F));
    CHECK(coarse.stride() == 4U);
    CHECK(coarse.firstCell() == std::array<std::int32_t, 2>{56, 20});
}

TEST_CASE("Height extract: coarse cells need every finer node, and a coarse surface of a plane is the plane")
{
    const testing::Plane plane{0.1, 0.2, 0.1};
    Quad quad(plane, {{3, 3}});
    HeightExtract out;
    HeightRequest request;
    request.cells = {{0, 0}, {15, 15}};
    request.stride = 4;
    const auto r = extractHeights(pageGrid(), quad.window(), request, out);
    CHECK(r.status == HeightStatus::Partial);
    CHECK(r.unknownCells == 1U); // only the coarse cell that covers node (3, 3)
    CHECK(out.cellState()[0] == CellState::Unknown);
    for (std::size_t t = 2; t < out.triangles().size(); ++t)
    {
        const auto& tri = out.triangles()[t];
        const double x = pageGrid().anchor().x.si() + static_cast<double>(tri.a.x.si());
        const double y = pageGrid().anchor().y.si() + static_cast<double>(tri.a.y.si());
        CHECK(std::abs(static_cast<double>(tri.a.z.si()) + out.origin().z.si() - plane.height(x, y)) < 2e-4);
    }
}

TEST_CASE("Height extract: corridor coverage delivers exactly the touched cells")
{
    const testing::Ramp ramp{.x0 = 1.0, .length = 1.0, .rise = 0.3};
    Quad quad(ramp);
    const auto corridor = CorridorBuilder::create(
        {{Length::fromSi(0.8), Length::fromSi(0.8), f64::Angle::fromSi(0.3)},
         {Length::fromSi(1.4), Length::fromSi(1.0), f64::Angle::fromSi(0.5)},
         {Length::fromSi(0.1), Length::fromSi(0.04), Length::fromSi(0.05)},
         8});
    REQUIRE(corridor);
    HeightExtract out;
    const auto pages = pageGrid().pagesTouching(corridor->bounds(), 0);
    REQUIRE_FALSE(pages.empty());
    const auto r = extractCorridor(pageGrid(), *corridor, quad.window(), 0, 1, out);
    REQUIRE(r.status == HeightStatus::Complete);
    const auto raster = corridor->rasterise(pageGrid().cellGrid(0), out.maskWords());
    std::size_t delivered = 0;
    const auto [cx, cy] = out.cells();
    for (std::int32_t b = 0; b < cy; ++b)
    {
        for (std::int32_t a = 0; a < cx; ++a)
        {
            const bool touched = raster.mask.contains(out.firstCell()[0] + a, out.firstCell()[1] + b);
            const auto s = out.cellState()[static_cast<std::size_t>(b * cx + a)];
            CHECK((s == CellState::Delivered) == touched);
            delivered += touched ? 1U : 0U;
        }
    }
    CHECK(delivered == r.deliveredCells);
    CHECK(r.deliveredCells < static_cast<std::size_t>(cx * cy)); // a diagonal corridor is not its bounding box
}

TEST_CASE("Height extract: ray casts agree with a brute-force test of every triangle")
{
    const testing::Berm berm{.y0 = 0.5, .radius = 1.5, .maxAngle = 0.5};
    Quad quad(berm);
    HeightExtract out;
    HeightRequest request;
    request.cells = {{40, 40}, {140, 120}};
    REQUIRE(extractHeights(pageGrid(), quad.window(), request, out).status == HeightStatus::Complete);
    std::mt19937 rng(11);
    std::uniform_real_distribution<float> px(-0.2F, 1.2F), py(-0.2F, 1.0F), az(0.2F, 0.9F);
    std::uniform_real_distribution<float> angle(0.0F, 6.2831853F), pitch(-1.3F, 0.4F);
    int hits = 0;
    for (int k = 0; k < 600; ++k)
    {
        const float pitchAngle = pitch(rng);
        const sub0math::Ray<LengthF, 3> ray{
            {LengthF::fromSi(px(rng)), LengthF::fromSi(py(rng)), LengthF::fromSi(az(rng))},
            {std::cos(pitchAngle) * std::cos(angle(rng)), std::cos(pitchAngle) * std::sin(angle(rng)), std::sin(pitchAngle)}};
        auto unit = ray;
        const float len = std::sqrt(unit.direction.x * unit.direction.x + unit.direction.y * unit.direction.y + unit.direction.z * unit.direction.z);
        unit.direction = {unit.direction.x / len, unit.direction.y / len, unit.direction.z / len};
        const auto fast = castRay(out, unit, LengthF::fromSi(3.0F));
        const auto slow = castRayBruteForce(out, unit, LengthF::fromSi(3.0F));
        REQUIRE(fast.hit == slow.hit);
        if (fast.hit)
        {
            ++hits;
            CHECK(fast.distance.si() == doctest::Approx(slow.distance.si()).epsilon(1e-4));
            CHECK(fast.witness.kind == WitnessKind::Face);
            const auto& n = fast.witness.normal;
            CHECK(n[0] * unit.direction.x + n[1] * unit.direction.y + n[2] * unit.direction.z <= 0.0); // faces the origin
        }
    }
    CHECK(hits > 100);
}

TEST_CASE("Batch height lookup gives exactly what heightAt gives")
{
    const testing::Berm berm{};
    const auto page = bakePage(berm, {0, 0, 0}, {{20, 20}});
    std::mt19937 rng(3);
    std::uniform_real_distribution<float> u(-0.1F, 2.7F);
    std::vector<LengthF> xs, ys, hs(1024);
    std::vector<std::uint8_t> known(1024);
    for (int k = 0; k < 1024; ++k)
    {
        xs.push_back(LengthF::fromSi(u(rng)));
        ys.push_back(LengthF::fromSi(u(rng)));
    }
    std::size_t expected = 0;
    const auto count = page.heightsAt(xs, ys, hs, known);
    for (std::size_t k = 0; k < xs.size(); ++k)
    {
        const auto s = page.heightAt(xs[k], ys[k]);
        REQUIRE((known[k] == 1U) == s.known);
        if (s.known)
        {
            CHECK(hs[k].si() == s.height.si());
            ++expected;
        }
    }
    CHECK(count == expected);
    std::vector<LengthF> shortOut(10);
    CHECK(page.heightsAt(xs, ys, shortOut, known) == 0U);
}

TEST_CASE("GridSpec domain and covers describe the authored rectangle")
{
    const auto grid = *GridSpec::create({.finestSpacing = Length::fromSi(0.02),
                                         .cellsPerSide = 40,
                                         .registration = Registration::CellCentre,
                                         .anchor = {Length::fromSi(0.0), Length::fromSi(-1.6)},
                                         .pages = {10, 4}});
    CHECK(grid.domain().hi.x.si() == 8.0);
    CHECK(grid.domain().hi.y.si() == 1.6);
    using Box = sub0math::Aabb<Length, 2>;
    CHECK(grid.covers(Box{{Length::fromSi(0.0), Length::fromSi(-1.6)}, {Length::fromSi(8.0), Length::fromSi(1.6)}}));
    CHECK_FALSE(grid.covers(Box{{Length::fromSi(-0.001), Length::fromSi(0.0)}, {Length::fromSi(1.0), Length::fromSi(1.0)}}));
    CHECK_FALSE(grid.covers(Box{{Length::fromSi(std::nan("")), Length::fromSi(0.0)}, {Length::fromSi(1.0), Length::fromSi(1.0)}}));
}

TEST_CASE("ColumnBuffer refuses to grow past its capacity and never truncates silently")
{
    ColumnBuffer<int> buffer;
    CHECK_FALSE(buffer.pushBack(1)); // nothing reserved
    buffer.reserve(3);
    CHECK(buffer.pushBack(1));
    CHECK(buffer.resize(3));
    buffer[2] = 9;
    CHECK_FALSE(buffer.pushBack(4));
    CHECK_FALSE(buffer.resize(4));
    CHECK(buffer.size() == 3U);
    const ColumnBuffer<int> copy = buffer;
    CHECK(copy.size() == 3U);
    CHECK(copy[2] == 9);
    CHECK(copy.capacity() == 3U);
}
