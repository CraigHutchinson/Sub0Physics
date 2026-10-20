#include <algorithm>
#include <cmath>
#include <cstdint>
#include <doctest.h>
#include <ostream>
#include <random>
#include <string_view>
#include <vector>

#include "sub0physics/terrain.hpp"
#include "sub0physics/terrain/testing.hpp"

using namespace sub0physics::terrain;
using namespace sub0physics::terrain::testing;
namespace f64 = sub0math::f64;

namespace
{

double length(const sub0math::Vec3d& v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

} // namespace

TEST_CASE("The catalogue has one name per surface")
{
    const auto surfaces = defaultSurfaces();
    std::vector<std::string_view> seen;
    for (std::size_t i = 0; i < surfaces.size(); ++i)
    {
        CHECK(nameOf(surfaces[i]) == SurfaceNames[i]);
        CHECK(std::find(seen.begin(), seen.end(), SurfaceNames[i]) == seen.end());
        seen.push_back(SurfaceNames[i]);
    }
    CHECK(seen.size() == 7U);
    CHECK(std::find(seen.begin(), seen.end(), "Square") == seen.end()); // Step, not Square
    CHECK(std::find(seen.begin(), seen.end(), "Step") != seen.end());
}

TEST_CASE("Catalogue normals are unit length and agree with the height field's gradient")
{
    std::mt19937_64 rng(20261020U);
    std::uniform_real_distribution<double> xs(0.0, 3.5), ys(-0.5, 0.8);
    constexpr double Eps = 1e-6;
    for (const AnySurface& surface : defaultSurfaces())
    {
        INFO(nameOf(surface));
        std::size_t compared = 0;
        for (int trial = 0; trial < 4000; ++trial)
        {
            const double x = xs(rng), y = ys(rng);
            const auto n = normalOf(surface, x, y);
            REQUIRE(length(n) == doctest::Approx(1.0).epsilon(1e-12));
            CHECK(n.z > 0.0); // points out of the solid
            // One-sided differences on both sides; only smooth points have a gradient to compare with.
            const double h = heightOf(surface, x, y);
            const double fx = (heightOf(surface, x + Eps, y) - h) / Eps, bx = (h - heightOf(surface, x - Eps, y)) / Eps;
            const double fy = (heightOf(surface, x, y + Eps) - h) / Eps, by = (h - heightOf(surface, x, y - Eps)) / Eps;
            if (std::abs(fx - bx) > 1e-4 || std::abs(fy - by) > 1e-4 || std::abs(fx) > 10.0 || std::abs(fy) > 10.0)
            {
                continue;
            }
            const double gx = 0.5 * (fx + bx), gy = 0.5 * (fy + by);
            const double g = std::sqrt(1.0 + gx * gx + gy * gy);
            CHECK(n.x == doctest::Approx(-gx / g).epsilon(1e-4).scale(1.0));
            CHECK(n.y == doctest::Approx(-gy / g).epsilon(1e-4).scale(1.0));
            CHECK(n.z == doctest::Approx(1.0 / g).epsilon(1e-4));
            ++compared;
        }
        CHECK(compared > 1000U);
    }
}

TEST_CASE("Catalogue surfaces have the documented shapes")
{
    const Ramp ramp{1.0, 2.0, 0.5};
    CHECK(ramp.height(0.0, 0.0) == 0.0);
    CHECK(ramp.height(2.0, 0.0) == doctest::Approx(0.25));
    CHECK(ramp.height(9.0, 0.0) == 0.5);
    CHECK(ramp.normal(2.0, 0.0).x == doctest::Approx(-0.25 / std::sqrt(1.0 + 0.0625)));

    const Step step{1.0, 0.1};
    CHECK(step.height(0.999, 0.0) == 0.0);
    CHECK(step.height(1.0, 0.0) == 0.1); // the point belongs to the piece that starts there
    CHECK(Step::riserNormal().x == -1.0);

    const RoundedBlock block{1.0, 1.4, -0.2, 0.2, 0.15, 0.03};
    CHECK(block.height(1.2, 0.0) == 0.15);
    CHECK(block.height(0.5, 0.0) == 0.0);
    CHECK(block.height(1.0 + 1e-6, 0.0) == doctest::Approx(0.15 - 0.03).epsilon(1e-2)); // the wall's top: the rounding ends vertical
    CHECK(block.height(0.9999, 0.0) == 0.0);
    const auto edge = block.normal(1.0 + 1e-6, 0.0);
    CHECK(edge.x == doctest::Approx(-1.0).epsilon(1e-3));
    CHECK(edge.z < 0.01);
    CHECK(block.normal(1.03, 0.0).z == doctest::Approx(1.0));

    const Bridge bridge;
    CHECK(bridge.height(1.5, 0.0) == bridge.deckTop);
    CHECK(bridge.spansAt(1.5, 0.0).count == 2U);
    CHECK(bridge.spansAt(1.5, 0.0).spans[1].lower == doctest::Approx(0.5));
    CHECK(bridge.spansAt(0.5, 0.0).count == 1U);
    CHECK(Bridge::undersideNormal().z == -1.0);

    const Berm berm{0.0, 2.0, 0.6};
    CHECK(berm.height(0.0, -1.0) == 0.0);
    CHECK(berm.height(0.0, 0.5) == doctest::Approx(2.0 - std::sqrt(4.0 - 0.25)));
    const double arcEnd = 2.0 * std::sin(0.6);
    CHECK(berm.height(0.0, arcEnd + 0.1) == doctest::Approx(2.0 * (1.0 - std::cos(0.6)) + 0.1 * std::tan(0.6)));
    CHECK(berm.normal(0.0, arcEnd + 0.1).y == doctest::Approx(-std::sin(0.6)));

    const Rut rut{0.0, 0.1, 0.05};
    CHECK(rut.height(0.0, 0.0) == doctest::Approx(-0.05));
    CHECK(rut.height(0.0, 0.1) == 0.0);
    CHECK(rut.normal(0.0, 0.0).z == doctest::Approx(1.0));
}

TEST_CASE("A height page reproduces the catalogue to within its interpolation and quantisation bounds")
{
    const auto grid = *GridSpec::create(f64::Length::fromSi(0.01), 1);
    const TileKey key{32768, 32768, 0}; // corner at the terrain origin: x and y from 0 to 2.56 m
    const auto build = [&](const auto& surface) {
        std::vector<f64::Length> heights(static_cast<std::size_t>(HeightPage::Nodes) * HeightPage::Nodes);
        std::vector<std::uint8_t> known(heights.size());
        sampleSurface(grid, key, surface, heights, known);
        return buildHeightPage(grid, key, SourceRevision{1}, heights, known);
    };
    std::mt19937_64 rng(5U);
    std::uniform_real_distribution<double> u(0.0, 2.56);
    const auto worst = [&](const auto& surface, const HeightPage& page) {
        double w = 0.0;
        for (int k = 0; k < 4000; ++k)
        {
            const double x = u(rng), y = u(rng);
            const auto z = page.heightAt(sub0math::Length::fromSi(static_cast<float>(x)), sub0math::Length::fromSi(static_cast<float>(y)));
            REQUIRE(z);
            w = std::max(w, std::abs(page.frame().origin.z.si() + static_cast<double>(z.height.si()) - surface.height(x, y)));
        }
        return w;
    };
    constexpr double Quantum = 0.5e-4 + 1e-5; // half a step, plus float rounding of the local height
    {
        const Plane plane{0.3, 0.1, -0.05};
        const auto page = build(plane);
        REQUIRE(page);
        CHECK(worst(plane, *page) <= Quantum); // linear: exact up to quantisation
    }
    {
        const Ramp ramp{1.0, 1.0, 0.4}; // both kinks on nodes: piecewise linear, exact up to quantisation
        const auto page = build(ramp);
        REQUIRE(page);
        CHECK(worst(ramp, *page) <= Quantum);
    }
    {
        const Rut rut{1.28, 0.2, 0.05}; // parabola: error at most h^2/8 times the curvature 2 d / w^2
        const auto page = build(rut);
        REQUIRE(page);
        CHECK(worst(rut, *page) <= 1e-4 / 8.0 * (2.0 * 0.05 / 0.04) + Quantum);
    }
    {
        const Step step{1.0, 0.1};
        const auto page = build(step);
        REQUIRE(page);
        // Away from the riser's cell the page is exact; inside it the fixed triangulation interpolates.
        for (const double x : {0.2, 0.9, 1.01, 1.7})
        {
            const auto z = page->heightAt(sub0math::Length::fromSi(static_cast<float>(x)), sub0math::Length::fromSi(1.0F));
            REQUIRE(z);
            CHECK(std::abs(page->frame().origin.z.si() + static_cast<double>(z.height.si()) - step.height(x, 1.0)) <= Quantum);
        }
        const auto between = page->heightAt(sub0math::Length::fromSi(0.995F), sub0math::Length::fromSi(1.0F));
        REQUIRE(between);
        CHECK(between.height.si() > 0.0F);
        CHECK(between.height.si() < 0.1F);
    }
}
