#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <doctest.h>
#include <numbers>
#include <random>
#include <vector>

#include "sub0physics/terrain.hpp"

using namespace sub0physics;
using namespace sub0physics::terrain;
namespace f64 = sub0math::f64;
using f64::Angle;
using f64::Length;

namespace
{

// The R0a source shape: 2 cm cells over x in [0, 8], y in [-1.6, 1.6].
const sub0math::RegularGrid<Length, 2> Cells{{Length::fromSi(0.0), Length::fromSi(-1.6)}, Length::fromSi(0.02), {401, 161}};

// R0a's fixed envelope (revision "fixed-1"): radius 0.36, width 0.06, margin 0.02, lookbehind 0.10.
constexpr double Radius = 0.36, Width = 0.06, Margin = 0.02, Behind = 0.10;
constexpr double HalfLength = Radius + Margin, HalfWidth = Width / 2.0 + Margin;

struct Pose2
{
    double x, y, yaw;
};

PlanPose plan(const Pose2& p)
{
    return {Length::fromSi(p.x), Length::fromSi(p.y), Angle::fromSi(p.yaw)};
}

CorridorSpec specOf(const Pose2& a, const Pose2& b, std::uint32_t k)
{
    return {plan(a), plan(b), {Length::fromSi(HalfLength), Length::fromSi(HalfWidth), Length::fromSi(Behind)}, k};
}

// ---- The oracles, written from the definitions in BikeLab's terrain-request (R0a) ----------------------------

struct Box
{
    double xMin, xMax, yMin, yMax;
};

// R0a's footprint corners in its order (rear right, front right, front left, rear left).
std::array<std::array<double, 2>, 4> corners(const Pose2& p)
{
    const double c = std::cos(p.yaw), s = std::sin(p.yaw);
    const std::array<std::array<double, 2>, 4> local{{{-HalfLength, -HalfWidth}, {HalfLength, -HalfWidth}, {HalfLength, HalfWidth}, {-HalfLength, HalfWidth}}};
    std::array<std::array<double, 2>, 4> out{};
    for (std::size_t i = 0; i < 4; ++i)
    {
        out[i] = {p.x + c * local[i][0] - s * local[i][1], p.y + s * local[i][0] + c * local[i][1]};
    }
    return out;
}

// R0a's sweptBox: the box of both poses' corners and their lookbehind copies, padded by R |delta yaw| / 2 where
// R is the corner radius.
Box sweptBox(const Pose2& a, const Pose2& b)
{
    const double pad = std::hypot(HalfLength, HalfWidth) * std::abs(b.yaw - a.yaw) / 2.0;
    std::vector<std::array<double, 2>> points;
    for (const auto& c : corners(a)) points.push_back(c);
    for (const auto& c : corners(b)) points.push_back(c);
    const double dx = b.x - a.x, dy = b.y - a.y, d = std::hypot(dx, dy);
    const std::size_t n = points.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        points.push_back(d != 0.0 ? std::array<double, 2>{points[i][0] - Behind * dx / d, points[i][1] - Behind * dy / d}
                                  : std::array<double, 2>{points[i][0] - Behind * std::cos(a.yaw), points[i][1] - Behind * std::sin(a.yaw)});
    }
    Box box{points[0][0], points[0][0], points[0][1], points[0][1]};
    for (const auto& [x, y] : points)
    {
        box = {std::min(box.xMin, x), std::max(box.xMax, x), std::min(box.yMin, y), std::max(box.yMax, y)};
    }
    return {box.xMin - pad, box.xMax + pad, box.yMin - pad, box.yMax + pad};
}

// Cells of a closed range on one axis, as R0a counts them.
struct CellRange
{
    int lo, hi;
};

CellRange axisCells(double lo, double hi, double origin, int cells)
{
    return {static_cast<int>(std::max(0.0, std::ceil((lo - origin) / 0.02) - 1.0)),
            static_cast<int>(std::min(static_cast<double>(cells - 1), std::floor((hi - origin) / 0.02)))};
}

// Every point of the footprint at pose p (corners, edge midpoints and interior points), and the lookbehind copies.
void footprintPoints(const Pose2& p, double shiftX, double shiftY, std::vector<std::array<double, 2>>& out)
{
    const double c = std::cos(p.yaw), s = std::sin(p.yaw);
    for (int iu = 0; iu <= 6; ++iu)
    {
        for (int iv = 0; iv <= 2; ++iv)
        {
            const double u = HalfLength * (-1.0 + iu / 3.0), v = HalfWidth * (-1.0 + iv);
            const double x = p.x + c * u - s * v, y = p.y + s * u + c * v;
            out.push_back({x, y});
            out.push_back({x + shiftX, y + shiftY});
        }
    }
}

Pose2 lerp(const Pose2& a, const Pose2& b, double t)
{
    return {a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, a.yaw + (b.yaw - a.yaw) * t};
}

struct Trial
{
    Pose2 a, b;
    std::uint32_t k;
};

Trial randomTrial(std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> px(1.0, 7.0), py(-1.0, 1.0), heading(-std::numbers::pi, std::numbers::pi);
    std::uniform_real_distribution<double> step(0.0, 0.6), turn(-1.5, 1.5);
    std::uniform_int_distribution<int> choose(0, 4);
    const Pose2 a{px(rng), py(rng), heading(rng)};
    const double dir = heading(rng);
    const double d = choose(rng) == 0 ? 0.0 : step(rng); // sometimes stationary
    const Pose2 b{a.x + d * std::cos(dir), a.y + d * std::sin(dir), choose(rng) == 0 ? a.yaw : a.yaw + turn(rng)};
    const std::array<std::uint32_t, 5> ks{1, 2, 4, 8, 16};
    return {a, b, ks[static_cast<std::size_t>(choose(rng))]};
}

} // namespace

TEST_CASE("Corridor specs are validated")
{
    const Pose2 a{1.0, 0.0, 0.0}, b{1.5, 0.0, 0.1};
    CHECK(CorridorBuilder::create(specOf(a, b, 8)).has_value());
    CHECK_FALSE(CorridorBuilder::create(specOf(a, b, 0)).has_value());
    CHECK_FALSE(CorridorBuilder::create(specOf(a, b, CorridorBuilder::MaxSubIntervals + 1)).has_value());
    CHECK(CorridorBuilder::create(specOf(a, b, CorridorBuilder::MaxSubIntervals)).has_value());
    CHECK_FALSE(CorridorBuilder::create(specOf({std::nan(""), 0.0, 0.0}, b, 8)).has_value());
    CHECK_FALSE(CorridorBuilder::create(specOf(a, {HUGE_VAL, 0.0, 0.0}, 8)).has_value());
    CHECK_FALSE(CorridorBuilder::create(specOf(a, {0.0, 0.0, 1e9}, 8)).has_value());
    auto spec = specOf(a, b, 8);
    spec.footprint.halfLength = Length{};
    CHECK_FALSE(CorridorBuilder::create(spec).has_value());
    spec = specOf(a, b, 8);
    spec.footprint.halfWidth = Length::fromSi(-0.1);
    CHECK_FALSE(CorridorBuilder::create(spec).has_value());
    spec = specOf(a, b, 8);
    spec.footprint.behind = Length::fromSi(-0.01);
    CHECK_FALSE(CorridorBuilder::create(spec).has_value());
    spec.footprint.behind = Length::fromSi(std::nan(""));
    const auto refused = CorridorBuilder::create(spec);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().kind == FaultKind::InvalidInput);
}

TEST_CASE("A stationary wheel gets the footprint and its lookbehind, along the start heading")
{
    const auto corridor = CorridorBuilder::create(specOf({4.0, 0.0, 0.0}, {4.0, 0.0, 0.0}, 1));
    REQUIRE(corridor);
    REQUIRE(corridor->boxes().size() == 1U);
    const CorridorBox& box = corridor->boxes()[0];
    CHECK(box.centre.x.si() == doctest::Approx(4.0 - Behind / 2.0));
    CHECK(box.centre.y.si() == doctest::Approx(0.0));
    CHECK(box.halfExtents.x.si() == doctest::Approx(HalfLength + Behind / 2.0).epsilon(1e-6));
    CHECK(box.halfExtents.y.si() == doctest::Approx(HalfWidth).epsilon(1e-6));
}

TEST_CASE("Sub-intervals tighten a turning corridor")
{
    const Pose2 a{3.0, 0.0, 0.0}, b{3.3, 0.2, 1.2};
    const auto coarse = CorridorBuilder::create(specOf(a, b, 1));
    const auto fine = CorridorBuilder::create(specOf(a, b, 16));
    REQUIRE(coarse);
    REQUIRE(fine);
    std::array<std::uint64_t, 4096> words{};
    const Raster rc = coarse->rasterise(Cells, words);
    REQUIRE(rc.status == RasterStatus::Complete);
    const std::size_t coarseCells = rc.cells;
    const Raster rf = fine->rasterise(Cells, words);
    REQUIRE(rf.status == RasterStatus::Complete);
    CHECK(rf.cells < coarseCells);
}

TEST_CASE("The corridor covers every cell the footprint touches along the motion (R0a's 1001-instant check)")
{
    std::mt19937_64 rng(20261020U);
    std::vector<std::uint64_t> words(4096);
    std::vector<std::array<double, 2>> points;
    std::size_t checkedCells = 0;
    for (int trial = 0; trial < 600; ++trial)
    {
        const Trial t = randomTrial(rng);
        const auto corridor = CorridorBuilder::create(specOf(t.a, t.b, t.k));
        REQUIRE(corridor);
        const Raster raster = corridor->rasterise(Cells, words);
        REQUIRE(raster.status == RasterStatus::Complete);

        // The lookbehind: opposite the displacement, along the start heading when stationary.
        const double dx = t.b.x - t.a.x, dy = t.b.y - t.a.y, d = std::hypot(dx, dy);
        const double shiftX = d != 0.0 ? -Behind * dx / d : -Behind * std::cos(t.a.yaw);
        const double shiftY = d != 0.0 ? -Behind * dy / d : -Behind * std::sin(t.a.yaw);
        points.clear();
        for (int n = 0; n <= 1000; ++n)
        {
            footprintPoints(lerp(t.a, t.b, n / 1000.0), shiftX, shiftY, points);
        }
        for (const auto& [x, y] : points)
        {
            const CellRange cx = axisCells(x, x, 0.0, 400), cy = axisCells(y, y, -1.6, 160);
            for (int j = cy.lo; j <= cy.hi; ++j)
            {
                for (int i = cx.lo; i <= cx.hi; ++i)
                {
                    INFO("trial " << trial << " k " << t.k << " cell " << i << "," << j);
                    REQUIRE(raster.mask.contains(i, j));
                    ++checkedCells;
                }
            }
        }
    }
    CHECK(checkedCells > 1000000U);
}

TEST_CASE("R0a's swept box contains the footprint; the corridor is the tighter cover")
{
    std::mt19937_64 rng(31U);
    std::vector<std::uint64_t> words(4096);
    std::vector<std::array<double, 2>> points;
    double corridorCells = 0.0, boxCells = 0.0;
    for (int trial = 0; trial < 600; ++trial)
    {
        const Trial t = randomTrial(rng);
        const Box swept = sweptBox(t.a, t.b);
        const auto corridor = CorridorBuilder::create(specOf(t.a, t.b, 8));
        REQUIRE(corridor);
        const Raster raster = corridor->rasterise(Cells, words);
        REQUIRE(raster.status == RasterStatus::Complete);

        // The oracle itself: every footprint point of the true motion lies in the swept box.
        const double dx = t.b.x - t.a.x, dy = t.b.y - t.a.y, d = std::hypot(dx, dy);
        points.clear();
        for (int n = 0; n <= 200; ++n)
        {
            footprintPoints(lerp(t.a, t.b, n / 200.0), d != 0.0 ? -Behind * dx / d : -Behind * std::cos(t.a.yaw),
                            d != 0.0 ? -Behind * dy / d : -Behind * std::sin(t.a.yaw), points);
        }
        for (const auto& [x, y] : points)
        {
            REQUIRE(x >= swept.xMin - 1e-12);
            REQUIRE(x <= swept.xMax + 1e-12);
            REQUIRE(y >= swept.yMin - 1e-12);
            REQUIRE(y <= swept.yMax + 1e-12);
        }

        // Both are covers of the same motion; the oriented boxes need fewer cells than the axis-aligned one.
        const CellRange sx = axisCells(swept.xMin, swept.xMax, 0.0, 400), sy = axisCells(swept.yMin, swept.yMax, -1.6, 160);
        boxCells += static_cast<double>((sx.hi - sx.lo + 1) * (sy.hi - sy.lo + 1));
        corridorCells += static_cast<double>(raster.cells);
    }
    INFO("corridor " << corridorCells << " cells, swept box " << boxCells << " cells");
    CHECK(corridorCells < boxCells);
}

TEST_CASE("Rasterising agrees with a brute-force overlap of every cell against every box")
{
    std::mt19937_64 rng(5U);
    std::vector<std::uint64_t> words(4096);
    for (int trial = 0; trial < 50; ++trial)
    {
        const Trial t = randomTrial(rng);
        const auto corridor = CorridorBuilder::create(specOf(t.a, t.b, t.k));
        REQUIRE(corridor);
        const Raster raster = corridor->rasterise(Cells, words);
        REQUIRE(raster.status == RasterStatus::Complete);
        std::size_t expected = 0;
        for (std::int32_t j = 0; j < 160; ++j)
        {
            for (std::int32_t i = 0; i < 400; ++i)
            {
                const auto lo = Cells.nodeAt({i, j});
                const sub0math::Aabb<Length, 2> cell{lo, Cells.nodeAt({i + 1, j + 1})};
                bool any = false;
                for (const CorridorBox& box : corridor->boxes())
                {
                    any = any || sub0math::overlaps(box, cell);
                }
                expected += any ? 1U : 0U;
                CHECK(raster.mask.contains(i, j) == any);
            }
        }
        CHECK(raster.cells == expected);
    }
}

TEST_CASE("Rasterising reports a missed grid or a small buffer, and writes nothing then")
{
    const auto corridor = CorridorBuilder::create(specOf({4.0, 0.0, 0.0}, {4.4, 0.1, 0.2}, 4));
    REQUIRE(corridor);
    std::array<std::uint64_t, 1> tiny{0xDEADBEEFU};
    const Raster small = corridor->rasterise(Cells, tiny);
    CHECK(small.status == RasterStatus::BufferTooSmall);
    CHECK(tiny[0] == 0xDEADBEEFU);
    CHECK(small.cells == 0U);
    CHECK(CorridorBuilder::wordsFor(small.mask.range) > 1U);

    const auto away = CorridorBuilder::create(specOf({40.0, 0.0, 0.0}, {40.4, 0.0, 0.0}, 2));
    REQUIRE(away);
    std::array<std::uint64_t, 4> words{};
    const Raster missed = away->rasterise(Cells, words);
    CHECK(missed.status == RasterStatus::OutsideGrid);
    CHECK(missed.mask.range.empty());
    CHECK_FALSE(missed.mask.contains(0, 0));

    // Partly outside: clipped to the grid, still complete for the part inside.
    const auto edge = CorridorBuilder::create(specOf({7.9, 0.0, 0.0}, {8.4, 0.0, 0.0}, 2));
    REQUIRE(edge);
    std::vector<std::uint64_t> big(4096);
    const Raster clipped = edge->rasterise(Cells, big);
    CHECK(clipped.status == RasterStatus::Complete);
    CHECK(clipped.mask.range.hi[0] == 399);
}

TEST_CASE("Rasterising works on the cell grid of a GridSpec at any level")
{
    const auto grid = *GridSpec::create({.finestSpacing = Length::fromSi(0.02),
                                         .levels = 3,
                                         .cellsPerSide = 40,
                                         .registration = Registration::CellCentre,
                                         .anchor = {Length::fromSi(0.0), Length::fromSi(-1.6)},
                                         .pages = {12, 4}});
    const auto corridor = CorridorBuilder::create(specOf({4.0, 0.0, 0.3}, {4.5, 0.2, 0.6}, 8));
    REQUIRE(corridor);
    std::vector<std::uint64_t> words(4096);
    std::size_t previous = 0;
    for (std::uint8_t level = 0; level < 3; ++level)
    {
        const Raster r = corridor->rasterise(grid.cellGrid(level), words);
        REQUIRE(r.status == RasterStatus::Complete);
        if (level > 0)
        {
            CHECK(r.cells < previous); // coarser cells, fewer of them
        }
        previous = r.cells;
        // Every page that owns a cell of the mask is among the pages the corridor's bounds touch.
        const auto pages = grid.pagesTouching(corridor->bounds(), level);
        REQUIRE_FALSE(pages.empty());
        for (std::int32_t j = r.mask.range.lo[1]; j <= r.mask.range.hi[1]; ++j)
        {
            for (std::int32_t i = r.mask.range.lo[0]; i <= r.mask.range.hi[0]; ++i)
            {
                if (!r.mask.contains(i, j)) continue;
                CHECK(i / 40 >= pages.lo[0]);
                CHECK(i / 40 <= pages.hi[0]);
                CHECK(j / 40 >= pages.lo[1]);
                CHECK(j / 40 <= pages.hi[1]);
            }
        }
    }
}
