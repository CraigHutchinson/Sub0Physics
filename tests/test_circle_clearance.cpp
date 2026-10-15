#include <algorithm>
#include <cmath>
#include <doctest.h>
#include <limits>
#include <numbers>
#include <random>
#include <vector>

#include "sub0physics/geometry/circle_clearance.hpp"

using namespace sub0physics;

namespace
{

constexpr double Tie = 1e-10;

PolylineProfile squareStep()
{
    return *PolylineProfile::create({{-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.1}, {0.2, 0.1}, {0.2, 0.0}, {1.0, 0.0}});
}

} // namespace

TEST_CASE("Clearance above flat ground is the height minus the radius")
{
    const auto flat = *PolylineProfile::create({{-1.0, 0.0}, {1.0, 0.0}});
    const CircleClearance c = measureCircleClearance(flat, {{0.3, 0.5}, 0.2, Tie});
    REQUIRE(c.status == ClearanceStatus::Available);
    CHECK(c.distance == doctest::Approx(0.5));
    CHECK(c.gap == doctest::Approx(0.3));
    REQUIRE(c.witnesses.size() == 1U);
    CHECK(c.witnesses[0].kind == FeatureKind::Segment);
    CHECK(c.witnesses[0].normal.y == doctest::Approx(1.0));
}

TEST_CASE("A square edge is approached through its corner")
{
    const CircleClearance c = measureCircleClearance(squareStep(), {{-0.1, 0.25}, 0.1, Tie});
    REQUIRE(c.status == ClearanceStatus::Available);
    REQUIRE(c.witnesses.size() == 1U);
    CHECK(c.witnesses[0].kind == FeatureKind::Vertex);
    CHECK(c.witnesses[0].index == 2U);
    CHECK(c.distance == doctest::Approx(std::hypot(0.1, 0.15)));
    CHECK(c.witnesses[0].normal.x == doctest::Approx(-0.1 / std::hypot(0.1, 0.15)));
}

TEST_CASE("Equidistant features are all reported once")
{
    const auto valley = *PolylineProfile::create({{-1.0, 1.0}, {0.0, 0.0}, {1.0, 1.0}});
    const CircleClearance c = measureCircleClearance(valley, {{0.0, 0.5}, 0.1, Tie});
    REQUIRE(c.status == ClearanceStatus::Available);
    CHECK(c.distance == doctest::Approx(0.5 / std::numbers::sqrt2));
    REQUIRE(c.witnesses.size() == 2U);
    CHECK(c.witnesses[0].kind == FeatureKind::Segment);
    CHECK(c.witnesses[0].index == 0U);
    CHECK(c.witnesses[1].index == 1U);
    CHECK(c.witnesses[0].normal.x == doctest::Approx(-c.witnesses[1].normal.x));
}

TEST_CASE("Clearance refuses undefined queries")
{
    const PolylineProfile profile = squareStep();
    CHECK(measureCircleClearance(profile, {{0.9, 1.0}, 0.2, Tie}).status == ClearanceStatus::OutsideDomain);
    CHECK(measureCircleClearance(profile, {{0.1, 0.1}, 0.05, Tie}).status == ClearanceStatus::OnBoundary);
    CHECK(measureCircleClearance(profile, {{0.1, 0.05}, 0.01, Tie}).status == ClearanceStatus::CentreInside);
}

TEST_CASE("Pruned search agrees with brute force on random terrain")
{
    std::mt19937_64 rng(20261015U);
    std::uniform_real_distribution<double> step(0.0, 0.05);
    std::uniform_real_distribution<double> rise(-0.04, 0.04);
    std::bernoulli_distribution vertical(0.1);
    std::vector<Point2> points{{-5.0, 0.0}};
    while (points.back().x < 5.0)
    {
        const Point2 last = points.back();
        const bool face = vertical(rng) && (points.size() < 2U || points[points.size() - 2U].x != last.x);
        points.push_back(face ? Point2{last.x, last.y + 0.05} : Point2{last.x + 0.001 + step(rng), last.y + rise(rng)});
    }
    const auto profile = *PolylineProfile::create(points);

    std::uniform_real_distribution<double> x(-4.0, 4.0);
    std::uniform_real_distribution<double> lift(0.0, 1.0);
    for (int trial = 0; trial < 2000; ++trial)
    {
        const double cx = x(rng);
        const Point2 centre{cx, profile.surfaceHeight(cx) + 1e-3 + lift(rng)};
        double brute = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i + 1U < points.size(); ++i)
        {
            const Point2 a = points[i];
            const Point2 b = points[i + 1U];
            const double dx = b.x - a.x;
            const double dy = b.y - a.y;
            const double t = std::clamp(((centre.x - a.x) * dx + (centre.y - a.y) * dy) / (dx * dx + dy * dy), 0.0, 1.0);
            brute = std::min(brute, std::hypot(centre.x - (a.x + t * dx), centre.y - (a.y + t * dy)));
        }
        const CircleClearance c = measureCircleClearance(profile, {centre, 0.5, Tie});
        REQUIRE(c.status == ClearanceStatus::Available);
        CHECK(c.distance == brute);
        CHECK_FALSE(c.witnesses.empty());
    }
}
