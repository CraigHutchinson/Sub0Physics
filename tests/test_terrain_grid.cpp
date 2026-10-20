#include <algorithm>
#include <cmath>
#include <cstdint>
#include <doctest.h>
#include <random>
#include <vector>

#include "sub0physics/terrain.hpp"

using namespace sub0physics;
using namespace sub0physics::terrain;
namespace f64 = sub0math::f64;
using f64::Length;

namespace
{

TerrainPoint at(double x, double y)
{
    return TerrainPoint{{Length::fromSi(x), Length::fromSi(y), Length{}}};
}

sub0math::Aabb<Length, 2> box(double x0, double y0, double x1, double y1)
{
    return {{Length::fromSi(x0), Length::fromSi(y0)}, {Length::fromSi(x1), Length::fromSi(y1)}};
}

// A source shaped like the R0a fixture: one level, 40-cell pages of 2 cm cells, sampled at cell centres, over
// x in [0, 8] and y in [-1.6, 1.6].
GridSpec r0aShaped()
{
    return *GridSpec::create({.finestSpacing = Length::fromSi(0.02),
                              .levels = 1,
                              .cellsPerSide = 40,
                              .registration = Registration::CellCentre,
                              .anchor = {Length::fromSi(0.0), Length::fromSi(-1.6)},
                              .pages = {10, 4}});
}

// The R0a oracle, from its definition: cell i of an axis is [origin + i h, origin + (i + 1) h); a closed box
// takes ceil(lo) - 1 to floor(hi) in grid units, limited to the cells.
struct AxisOracle
{
    double origin;
    double spacing;
    int cells;

    [[nodiscard]] std::pair<int, int> closed(double lo, double hi) const
    {
        const double a = std::max(0.0, std::ceil((lo - origin) / spacing) - 1.0);
        const double b = std::min(static_cast<double>(cells - 1), std::floor((hi - origin) / spacing));
        return {static_cast<int>(a), static_cast<int>(b)};
    }

    [[nodiscard]] int owner(double p) const
    {
        return static_cast<int>(std::min(static_cast<double>(cells - 1), std::floor((p - origin) / spacing)));
    }
};

} // namespace

TEST_CASE("An R0a-shaped source is expressible exactly")
{
    const GridSpec grid = r0aShaped();
    CHECK(grid.cellsPerSide() == 40);
    CHECK(grid.samplesPerSide() == 40);
    CHECK(grid.registration() == Registration::CellCentre);
    CHECK(grid.pageCount(0, 0) == 10);
    CHECK(grid.pageCount(0, 1) == 4);
    CHECK(grid.pageSide(0).si() == doctest::Approx(0.8));

    // Samples sit where the fixture evaluates its surface: (ix + 0.5) * 0.02 and -1.6 + (iy + 0.5) * 0.02.
    for (int ty = 0; ty < 4; ++ty)
    {
        for (int tx = 0; tx < 10; ++tx)
        {
            for (const int i : {0, 1, 17, 39})
            {
                for (const int j : {0, 5, 38, 39})
                {
                    const auto p = grid.samplePosition({tx, ty, 0}, i, j);
                    CHECK(p.x.si() == (static_cast<double>(tx * 40 + i) + 0.5) * 0.02);
                    CHECK(p.y.si() == -1.6 + (static_cast<double>(ty * 40 + j) + 0.5) * 0.02);
                }
            }
        }
    }

    // Ownership agrees with the fixture's cell for random points, for points on the lattice, and on the far edge.
    const AxisOracle ox{0.0, 0.02, 400}, oy{-1.6, 0.02, 160};
    std::mt19937_64 rng(1020U);
    std::uniform_real_distribution<double> ux(0.0, 8.0), uy(-1.6, 1.6);
    std::vector<std::pair<double, double>> points;
    for (int k = 0; k < 20000; ++k) points.emplace_back(ux(rng), uy(rng));
    for (int ix = 0; ix <= 400; ix += 3)
        for (int iy = 0; iy <= 160; iy += 7) points.emplace_back(0.02 * ix, -1.6 + 0.02 * iy);
    points.emplace_back(8.0, 1.6);
    for (const auto& [x, y] : points)
    {
        const CellLookup c = grid.cellOf(at(x, y), 0);
        REQUIRE(c.known);
        CHECK(c.page.i * 40 + c.i == ox.owner(x));
        CHECK(c.page.j * 40 + c.j == oy.owner(y));
        CHECK(c.page.level == 0U);
    }
    CHECK(grid.cellOf(at(8.0, 1.6), 0).page == TileKey{9, 3, 0}); // the far edge belongs to the last page
    CHECK_FALSE(grid.cellOf(at(8.001, 0.0), 0));
    CHECK_FALSE(grid.cellOf(at(-0.001, 0.0), 0));
    CHECK_FALSE(grid.cellOf(at(1.0, -1.601), 0));
}

TEST_CASE("Closed coverage equals the fixture's cell rule, and pages follow from cells")
{
    const GridSpec grid = r0aShaped();
    const AxisOracle ox{0.0, 0.02, 400}, oy{-1.6, 0.02, 160};
    std::mt19937_64 rng(77U);
    std::uniform_real_distribution<double> ux(0.0, 8.0), uy(-1.6, 1.6), size(0.0, 0.9);
    std::uniform_int_distribution<int> lattice(0, 3);
    for (int k = 0; k < 20000; ++k)
    {
        double x0 = ux(rng), y0 = uy(rng);
        double x1 = std::min(8.0, x0 + size(rng)), y1 = std::min(1.6, y0 + size(rng));
        if (lattice(rng) == 0) // edges exactly on cell boundaries
        {
            x0 = 0.02 * std::floor(x0 / 0.02);
            x1 = 0.02 * std::ceil(x1 / 0.02);
        }
        if (lattice(rng) == 0)
        {
            y0 = -1.6 + 0.02 * std::floor((y0 + 1.6) / 0.02);
            y1 = std::min(1.6, -1.6 + 0.02 * std::ceil((y1 + 1.6) / 0.02));
        }
        const auto [lx, hx] = ox.closed(x0, x1);
        const auto [ly, hy] = oy.closed(y0, y1);
        const auto cells = grid.cellsTouching(box(x0, y0, x1, y1), 0);
        CHECK(cells.lo[0] == lx);
        CHECK(cells.hi[0] == hx);
        CHECK(cells.lo[1] == ly);
        CHECK(cells.hi[1] == hy);
        const auto pages = grid.pagesTouching(box(x0, y0, x1, y1), 0);
        CHECK(pages.lo[0] == lx / 40);
        CHECK(pages.hi[0] == hx / 40);
        CHECK(pages.lo[1] == ly / 40);
        CHECK(pages.hi[1] == hy / 40);
    }
    // A box that touches a page edge exactly takes the pages on both sides.
    const auto both = grid.pagesTouching(box(0.7, 0.0, 0.8, 0.1), 0);
    CHECK(both.lo[0] == 0);
    CHECK(both.hi[0] == 1);
    // Misses and non-finite corners give empty ranges, never a converted integer.
    CHECK(grid.cellsTouching(box(9.0, 0.0, 10.0, 1.0), 0).empty());
    CHECK(grid.cellsTouching(box(std::nan(""), 0.0, 1.0, 1.0), 0).empty());
    CHECK(grid.pagesTouching(box(0.0, 0.0, HUGE_VAL, 1.0), 0).empty());
    CHECK(grid.pagesTouching(box(0.0, 0.0, 1.0, 1.0), 1).empty()); // level outside the grid
}

TEST_CASE("Registration changes where samples are, not which cell owns a point")
{
    const auto make = [](Registration r) {
        return *GridSpec::create({.finestSpacing = Length::fromSi(0.1), .cellsPerSide = 16, .registration = r, .pages = {4, 4}});
    };
    const GridSpec nodes = make(Registration::NodeIsPoint);
    const GridSpec centres = make(Registration::CellCentre);
    CHECK(nodes.samplesPerSide() == 17);
    CHECK(centres.samplesPerSide() == 16);
    const TileKey key{2, 1, 0};
    const auto corner = nodes.pageCorner(key);
    CHECK(nodes.samplePosition(key, 0, 0).x == corner.x);
    CHECK(nodes.samplePosition(key, 16, 16).x.si() == doctest::Approx(corner.x.si() + 1.6));
    CHECK(centres.samplePosition(key, 0, 0).x.si() == doctest::Approx(corner.x.si() + 0.05));
    CHECK(centres.samplePosition(key, 15, 15).y.si() == doctest::Approx(corner.y.si() + 1.55));
    // A cell-centre sample is inside the cell that owns it.
    for (int i = 0; i < 16; i += 5)
    {
        const auto p = centres.samplePosition(key, i, 3);
        const CellLookup c = centres.cellOf(TerrainPoint{p}, 0);
        REQUIRE(c.known);
        CHECK(c.page == key);
        CHECK(c.i == i);
        CHECK(c.j == 3);
    }
    for (const double v : {0.0, 0.049, 0.05, 0.1, 0.55, 1.599})
    {
        const CellLookup a = nodes.cellOf(at(v, v), 0), b = centres.cellOf(at(v, v), 0);
        CHECK(a.page == b.page);
        CHECK(a.i == b.i);
        CHECK(a.j == b.j);
    }
}

TEST_CASE("Levels double the spacing, and a coarse page is the parent of its fine pages")
{
    const auto grid = *GridSpec::create({.finestSpacing = Length::fromSi(0.5),
                                         .levels = 4,
                                         .cellsPerSide = 8,
                                         .anchor = {Length::fromSi(-3.0), Length::fromSi(7.0)},
                                         .pages = {13, 6}});
    CHECK(grid.pageCount(0, 0) == 13);
    CHECK(grid.pageCount(1, 0) == 7); // rounds up: the coarse page may extend past the authored rectangle
    CHECK(grid.pageCount(2, 0) == 4);
    CHECK(grid.pageCount(3, 0) == 2);
    CHECK(grid.pageCount(3, 1) == 1);
    CHECK(grid.spacing(2).si() == 2.0);
    CHECK(grid.pageSide(2).si() == 16.0);
    std::mt19937_64 rng(3U);
    std::uniform_real_distribution<double> ux(-3.0, -3.0 + 13 * 4.0 - 1e-9), uy(7.0, 7.0 + 6 * 4.0 - 1e-9);
    for (int k = 0; k < 5000; ++k)
    {
        const auto p = at(ux(rng), uy(rng));
        const PageLookup fine = grid.pageOf(p, 0);
        REQUIRE(fine);
        for (std::uint8_t level = 1; level < 4; ++level)
        {
            const PageLookup coarse = grid.pageOf(p, level);
            REQUIRE(coarse);
            CHECK(coarse.key.i == fine.key.i >> level);
            CHECK(coarse.key.j == fine.key.j >> level);
        }
    }
    // A point beyond the last fine page but inside the rounded-up coarse page has a coarse page only.
    const auto edge = at(-3.0 + 13 * 4.0 + 1.0, 8.0);
    CHECK_FALSE(grid.pageOf(edge, 0));
    CHECK(grid.pageOf(edge, 1).key.i == 6);
}

TEST_CASE("A page corner is owned by its own page, and a point inside it by that page")
{
    for (const GridSpec& grid : {r0aShaped(), *GridSpec::create(Length::fromSi(0.01), 3),
                                 *GridSpec::create({.finestSpacing = Length::fromSi(0.37),
                                                    .levels = 3,
                                                    .anchor = {Length::fromSi(-123.456), Length::fromSi(98.7654)},
                                                    .pages = {300, 200}})})
    {
        std::mt19937_64 rng(11U);
        for (std::uint8_t level = 0; level < grid.levels(); ++level)
        {
            std::uniform_int_distribution<int> ki(0, grid.pageCount(level, 0) - 1), kj(0, grid.pageCount(level, 1) - 1);
            for (int k = 0; k < 500; ++k)
            {
                const TileKey key{ki(rng), kj(rng), level};
                const auto corner = grid.pageCorner(key);
                const double side = grid.pageSide(level).si();
                const PageLookup inside = grid.pageOf(at(corner.x.si() + 0.37 * side, corner.y.si() + 0.62 * side), level);
                REQUIRE(inside);
                CHECK(inside.key == key);
                // A corner is a shared edge: it belongs to the page on the positive side, up to the rounding of the
                // corner's own position (here: within one page of the key, never beyond).
                const PageLookup exact = grid.pageOf(TerrainPoint{corner}, level);
                REQUIRE(exact);
                CHECK(std::abs(exact.key.i - key.i) <= 1);
                CHECK(std::abs(exact.key.j - key.j) <= 1);
            }
        }
    }
}

TEST_CASE("Grid layouts are validated")
{
    const auto make = [](GridLayout layout) { return GridSpec::create(layout); };
    const GridLayout base{.finestSpacing = Length::fromSi(0.02), .cellsPerSide = 40, .pages = {10, 4}};
    CHECK(make(base).has_value());
    auto bad = base;
    bad.cellsPerSide = 0;
    CHECK_FALSE(make(bad).has_value());
    bad = base;
    bad.cellsPerSide = GridSpec::MaxCellsPerSide + 1;
    CHECK_FALSE(make(bad).has_value());
    bad = base;
    bad.pages = {0, 4};
    CHECK_FALSE(make(bad).has_value());
    bad = base;
    bad.pages = {1 << 25, 1 << 25}; // 2^25 * 40 cells exceeds 2^30
    CHECK_FALSE(make(bad).has_value());
    bad = base;
    bad.anchor = {Length::fromSi(std::nan("")), Length{}};
    CHECK_FALSE(make(bad).has_value());
    bad = base;
    bad.anchor = {Length{}, Length::fromSi(HUGE_VAL)};
    CHECK_FALSE(make(bad).has_value());
    bad = base;
    bad.registration = static_cast<Registration>(7);
    CHECK_FALSE(make(bad).has_value());
    bad = base;
    bad.levels = 17;
    CHECK_FALSE(make(bad).has_value());
    bad = base;
    bad.finestSpacing = Length::fromSi(std::nan(""));
    CHECK_FALSE(make(bad).has_value());
    CHECK(make(bad).error().kind == FaultKind::InvalidInput);

    // Defaults: 256 cells per page, node registration; the two-argument form is the same grid, centred.
    const GridLayout defaults{.finestSpacing = Length::fromSi(0.01)};
    const auto d = make(defaults);
    REQUIRE(d);
    CHECK(d->cellsPerSide() == 256);
    CHECK(d->registration() == Registration::NodeIsPoint);
    CHECK(d->samplesPerSide() == 257);
    CHECK(GridSpec::create(Length::fromSi(0.01), 1)->cellsPerSide() == GridSpec::DefaultCellsPerSide);
}

TEST_CASE("Height pages need the 256-cell node-registered layout")
{
    const std::size_t count = static_cast<std::size_t>(HeightPage::Nodes) * HeightPage::Nodes;
    const std::vector<Length> heights(count, Length::fromSi(1.0));
    const std::vector<std::uint8_t> known(count, 1);
    const auto good = *GridSpec::create(Length::fromSi(0.01), 2);
    CHECK(buildHeightPage(good, {32768, 32768, 0}, SourceRevision{1}, heights, known));
    CHECK_FALSE(buildHeightPage(good, {-1, 0, 0}, SourceRevision{1}, heights, known));    // before the extent
    CHECK_FALSE(buildHeightPage(good, {65536, 0, 0}, SourceRevision{1}, heights, known)); // after it
    CHECK(buildHeightPage(good, {32767, 0, 1}, SourceRevision{1}, heights, known));       // a coarse page in range
    CHECK_FALSE(buildHeightPage(good, {32768, 0, 1}, SourceRevision{1}, heights, known)); // coarse pages are half as many
    CHECK_FALSE(buildHeightPage(r0aShaped(), {0, 0, 0}, SourceRevision{1}, heights, known));
    const auto centred = *GridSpec::create({.finestSpacing = Length::fromSi(0.01), .registration = Registration::CellCentre});
    const auto refused = buildHeightPage(centred, {1, 1, 0}, SourceRevision{1}, heights, known);
    REQUIRE_FALSE(refused);
    CHECK(refused.error().kind == FaultKind::InvalidInput);
}
