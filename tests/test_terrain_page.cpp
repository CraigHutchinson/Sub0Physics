#include <cmath>
#include <concepts>
#include <cstdint>
#include <doctest.h>
#include <type_traits>
#include <vector>

#include "sub0physics/terrain.hpp"

using namespace sub0physics;
using namespace sub0physics::terrain;
using sub0math::Length;
namespace f64 = sub0math::f64;

namespace
{

template<typename From, typename To>
concept ConvertsTo = std::convertible_to<From, To>;

constexpr std::size_t NodeCount = static_cast<std::size_t>(HeightPage::Nodes) * HeightPage::Nodes;

struct Samples
{
    std::vector<f64::Length> heights = std::vector<f64::Length>(NodeCount);
    std::vector<std::uint8_t> known = std::vector<std::uint8_t>(NodeCount, 1);
};

template<typename F>
Samples sample(const GridSpec& grid, const TileKey& key, F&& height)
{
    Samples s;
    const auto corner = grid.pageCorner(key);
    const double h = grid.spacing(key.level).si();
    for (std::int32_t j = 0; j < HeightPage::Nodes; ++j)
        for (std::int32_t i = 0; i < HeightPage::Nodes; ++i)
            s.heights[static_cast<std::size_t>(j) * HeightPage::Nodes + static_cast<std::size_t>(i)] =
                f64::Length::fromSi(height(corner.x.si() + i * h, corner.y.si() + j * h));
    return s;
}

const GridSpec Grid = *GridSpec::create(f64::Length::fromSi(0.01), 4);

} // namespace

// Frames cannot be mixed: a terrain point is not a local point, and height queries take local coordinates.
static_assert(!ConvertsTo<TerrainPoint, LocalPoint>);
static_assert(!ConvertsTo<LocalPoint, TerrainPoint>);
static_assert(!ConvertsTo<sub0math::Point3<f64::Length>, TerrainPoint>);

TEST_CASE("Grid spec: spacing per level, page lookup and corners")
{
    CHECK(Grid.spacing(0).si() == doctest::Approx(0.01));
    CHECK(Grid.spacing(3).si() == doctest::Approx(0.08));
    CHECK(Grid.pageSide(0).si() == doctest::Approx(2.56));
    const TerrainPoint p{{f64::Length::fromSi(-0.1), f64::Length::fromSi(5.2), f64::Length{}}};
    REQUIRE(Grid.pageOf(p, 0));
    CHECK(Grid.pageOf(p, 0).key == TileKey{-1, 2, 0});
    REQUIRE(Grid.pageOf(p, 2));
    CHECK(Grid.pageOf(p, 2).key == TileKey{-1, 0, 2});
    CHECK(Grid.pageCorner({-1, 2, 0}).x.si() == doctest::Approx(-2.56));
}

TEST_CASE("Grid spec: create validates levels and spacing")
{
    const auto make = [](double spacing, int levels) {
        return GridSpec::create(f64::Length::fromSi(spacing), static_cast<std::uint8_t>(levels));
    };
    CHECK(make(0.01, 1).has_value());
    CHECK(make(0.01, 16).has_value());
    CHECK(make(GridSpec::MinSpacing, 4).has_value());
    CHECK(make(GridSpec::MaxSpacing, 16).has_value());
    CHECK_FALSE(make(0.01, 0).has_value());
    CHECK_FALSE(make(0.01, 17).has_value());
    CHECK_FALSE(make(0.01, 255).has_value());
    CHECK_FALSE(make(0.0, 4).has_value());
    CHECK_FALSE(make(-0.01, 4).has_value());
    CHECK_FALSE(make(std::nan(""), 4).has_value());
    CHECK_FALSE(make(HUGE_VAL, 4).has_value());
    CHECK_FALSE(make(GridSpec::MinSpacing / 2, 4).has_value());
    CHECK_FALSE(make(GridSpec::MaxSpacing * 2, 4).has_value());
    const auto bad = make(0.01, 0);
    CHECK(bad.error().kind == FaultKind::InvalidInput);
    // The only way to a grid is create: no aggregate or default construction.
    static_assert(!std::is_default_constructible_v<GridSpec>);
    static_assert(!std::is_constructible_v<GridSpec, f64::Length, std::uint8_t>);
    CHECK(make(0.01, 16)->levels() == 16);
    CHECK(Grid.finestSpacing().si() == 0.01);
}

TEST_CASE("Grid spec: page lookup refuses what it cannot index")
{
    const auto at = [](double x, double y) {
        return TerrainPoint{{f64::Length::fromSi(x), f64::Length::fromSi(y), f64::Length{}}};
    };
    const double nan = std::nan("");
    CHECK(Grid.pageOf(at(1.0, 1.0), 0));
    CHECK(Grid.pageOf(at(1.0, 1.0), 3));
    CHECK_FALSE(Grid.pageOf(at(1.0, 1.0), 4));   // level outside the grid
    CHECK_FALSE(Grid.pageOf(at(1.0, 1.0), 255));
    CHECK_FALSE(Grid.pageOf(at(nan, 1.0), 0));
    CHECK_FALSE(Grid.pageOf(at(1.0, nan), 0));
    CHECK_FALSE(Grid.pageOf(at(HUGE_VAL, 1.0), 0));
    CHECK_FALSE(Grid.pageOf(at(1.0, -HUGE_VAL), 0));
    CHECK_FALSE(Grid.pageOf(at(1e300, 1.0), 0));  // index beyond 32 bits
    CHECK_FALSE(Grid.pageOf(at(-1e300, 1.0), 0));
    const double side = Grid.pageSide(0).si();
    CHECK_FALSE(Grid.pageOf(at(side * 2147483648.5, 0.0), 0)); // one past INT32_MAX
    CHECK(Grid.pageOf(at(side * 2147483647.5, 0.0), 0).key.i == 2147483647);
    CHECK(Grid.pageOf(at(side * -2147483647.5, 0.0), 0).key.i == -2147483647 - 1);
    CHECK_FALSE(Grid.pageOf(at(side * -2147483648.5, 0.0), 0));
}

TEST_CASE("Page frames keep positions exact far from the terrain origin")
{
    // 1000 km out: float alone would round to 6 cm; the frame keeps micrometres.
    const TileKey key{390625, -390625, 0};
    const auto corner = Grid.pageCorner(key);
    const PageFrame frame{key, {corner.x, corner.y, f64::Length::fromSi(812.0)}};
    const TerrainPoint p{{corner.x + f64::Length::fromSi(1.234567), corner.y + f64::Length::fromSi(0.000321),
                          f64::Length::fromSi(812.0345)}};
    const LocalPoint local = toLocal(frame, p);
    CHECK(local.x.si() == doctest::Approx(1.234567).epsilon(1e-7));
    CHECK(local.z.si() == doctest::Approx(0.0345).epsilon(1e-6));
    const TerrainPoint back = toTerrain(frame, local);
    CHECK(std::abs((back.value().x - p.value().x).si()) < 1e-6);
    CHECK(std::abs((back.value().y - p.value().y).si()) < 1e-6);
    CHECK(std::abs((back.value().z - p.value().z).si()) < 1e-6);
}

TEST_CASE("Height pages reproduce a plane exactly and quantise to 0.1 mm")
{
    const TileKey key{3, -2, 0};
    const auto plane = [](double x, double y) { return 100.0 + 0.2 * x - 0.1 * y; };
    const Samples s = sample(Grid, key, plane);
    const auto page = buildHeightPage(Grid, key, SourceRevision{7}, s.heights, s.known);
    REQUIRE(page);
    CHECK(page->revision() == SourceRevision{7});
    const auto& frame = page->frame();
    double worst = 0.0;
    for (int k = 0; k < 200; ++k)
    {
        const float x = 0.0123F * static_cast<float>(k);
        const float y = 2.56F - 0.0117F * static_cast<float>(k);
        const auto z = page->heightAt(Length::fromSi(x), Length::fromSi(y));
        REQUIRE(z);
        const double terrainZ = frame.origin.z.si() + static_cast<double>(z.height.si());
        const double expected = plane(frame.origin.x.si() + x, frame.origin.y.si() + y);
        worst = std::max(worst, std::abs(terrainZ - expected));
    }
    CHECK(worst <= 0.5e-4 + 1e-6); // half a quantum, plus float rounding of a ~0.5 m local height
}

TEST_CASE("Unknown nodes are never interpolated across")
{
    const TileKey key{0, 0, 0};
    Samples s = sample(Grid, key, [](double, double) { return 1.0; });
    s.known[static_cast<std::size_t>(10) * HeightPage::Nodes + 10] = 0; // node (10, 10)
    const auto page = buildHeightPage(Grid, key, SourceRevision{1}, s.heights, s.known);
    REQUIRE(page);
    CHECK_FALSE(page->known(10, 10));
    // Both triangles touching node (10,10) from cell (9,9) and cell (10,10) are refused.
    CHECK_FALSE(page->heightAt(Length::fromSi(0.0995F), Length::fromSi(0.0991F)));
    CHECK_FALSE(page->heightAt(Length::fromSi(0.1005F), Length::fromSi(0.1001F)));
    CHECK(page->heightAt(Length::fromSi(0.2F), Length::fromSi(0.2F)));
    CHECK_FALSE(page->heightAt(Length::fromSi(-0.001F), Length::fromSi(0.5F))); // outside the page
}

TEST_CASE("The min/max pyramid bounds every node, and coarse blocks bound their children")
{
    const TileKey key{0, 0, 1};
    Samples s = sample(Grid, key, [](double x, double y) { return 0.3 * std::sin(3.0 * x) * std::cos(2.0 * y); });
    s.known[0] = 0;
    const auto page = buildHeightPage(Grid, key, SourceRevision{1}, s.heights, s.known);
    REQUIRE(page);
    std::size_t violations = 0;
    for (std::int32_t j = 0; j < HeightPage::Nodes; ++j)
        for (std::int32_t i = 0; i < HeightPage::Nodes; ++i)
        {
            if (!page->known(i, j)) continue;
            const float z = page->nodeHeight(i, j).si();
            for (std::int32_t level = 0; level < HeightPage::PyramidLevels; ++level)
            {
                const std::int32_t span = HeightPage::BlockCells << level;
                // A node on a block edge belongs to both neighbouring blocks; check the one below it.
                const std::int32_t bi = std::min(i / span, (32 >> level) - 1);
                const std::int32_t bj = std::min(j / span, (32 >> level) - 1);
                const HeightRange r = page->blockRange(level, bi, bj);
                if (!r.anyKnown || HeightPage::toHeight(r.lo).si() > z || HeightPage::toHeight(r.hi).si() < z) ++violations;
            }
        }
    CHECK(violations == 0);
    CHECK(page->blockRange(0, 0, 0).anyUnknown);
    CHECK(page->blockRange(HeightPage::PyramidLevels - 1, 0, 0).anyUnknown);
    CHECK_FALSE(page->blockRange(0, 5, 5).anyUnknown);
}

TEST_CASE("Invalid samples are refused")
{
    const TileKey key{0, 0, 0};
    Samples s = sample(Grid, key, [](double, double) { return 0.0; });
    CHECK(!buildHeightPage(Grid, TileKey{0, 0, 9}, SourceRevision{1}, s.heights, s.known));        // level
    CHECK(!buildHeightPage(Grid, key, SourceRevision{1}, std::span(s.heights).first(10), s.known)); // count
    s.heights[5] = f64::Length::fromSi(std::nan(""));
    CHECK(!buildHeightPage(Grid, key, SourceRevision{1}, s.heights, s.known));                      // non-finite
    s.known[5] = 0;
    CHECK(buildHeightPage(Grid, key, SourceRevision{1}, s.heights, s.known));                       // ignored when unknown
    s.heights[6] = f64::Length::fromSi(7.0);
    const auto tall = buildHeightPage(Grid, key, SourceRevision{1}, s.heights, s.known);            // 7 m range
    REQUIRE_FALSE(tall);
    CHECK(tall.error().kind == FaultKind::InvalidInput);
}
