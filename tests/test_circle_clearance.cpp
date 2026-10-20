#include <algorithm>
#include <array>
#include <cmath>
#include <doctest.h>
#include <limits>
#include <numbers>
#include <random>
#include <span>
#include <sub0math/units.hpp>
#include <utility>
#include <vector>

#include "sub0physics/geometry/circle_clearance.hpp"
#include "test_profiles.hpp"

using namespace sub0physics;
using sub0math::metres;
using test_profiles::at;
using test_profiles::metresPoints;

namespace
{

const sub0math::f64::Length Tie = metres(1e-10);

PolylineProfile squareStep()
{
    return *PolylineProfile::create(metresPoints({{-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.1}, {0.2, 0.1}, {0.2, 0.0}, {1.0, 0.0}}));
}

// A clearance with the buffer its witnesses were written to.
struct Measured : CircleClearance
{
    std::array<Witness2, DefaultWitnessCapacity> buffer{};

    [[nodiscard]] std::span<const Witness2> witnesses() const noexcept
    {
        return std::span<const Witness2>(buffer).first(witnessCount);
    }
};

Measured query(const PolylineProfile& profile, double cx, double cy, double radius, std::size_t capacity = DefaultWitnessCapacity)
{
    Measured m;
    static_cast<CircleClearance&>(m) = measureCircleClearance(profile, {{at(cx, cy), metres(radius)}, Tie},
                                                              std::span<Witness2>(m.buffer).first(capacity));
    return m;
}

} // namespace

TEST_CASE("Clearance above flat ground is the height minus the radius")
{
    const auto flat = *PolylineProfile::create(metresPoints({{-1.0, 0.0}, {1.0, 0.0}}));
    const Measured c = query(flat, 0.3, 0.5, 0.2);
    REQUIRE(c.status == ClearanceStatus::Available);
    CHECK(c.distance.si() == doctest::Approx(0.5));
    CHECK(c.gap.si() == doctest::Approx(0.3));
    REQUIRE(c.witnesses().size() == 1U);
    CHECK(c.witnesses()[0].kind == WitnessKind::Edge);
    CHECK(c.witnesses()[0].normal.y == doctest::Approx(1.0));
}

TEST_CASE("A square edge is approached through its corner")
{
    const Measured c = query(squareStep(), -0.1, 0.25, 0.1);
    REQUIRE(c.status == ClearanceStatus::Available);
    REQUIRE(c.witnesses().size() == 1U);
    CHECK(c.witnesses()[0].kind == WitnessKind::Vertex);
    CHECK(c.witnesses()[0].index == 2U);
    CHECK(c.distance.si() == doctest::Approx(std::hypot(0.1, 0.15)));
    CHECK(c.witnesses()[0].normal.x == doctest::Approx(-0.1 / std::hypot(0.1, 0.15)));
}

TEST_CASE("Equidistant features are all reported once")
{
    const auto valley = *PolylineProfile::create(metresPoints({{-1.0, 1.0}, {0.0, 0.0}, {1.0, 1.0}}));
    const Measured c = query(valley, 0.0, 0.5, 0.1);
    REQUIRE(c.status == ClearanceStatus::Available);
    CHECK(c.distance.si() == doctest::Approx(0.5 / std::numbers::sqrt2));
    REQUIRE(c.witnesses().size() == 2U);
    CHECK(c.witnesses()[0].kind == WitnessKind::Edge);
    CHECK(c.witnesses()[0].index == 0U);
    CHECK(c.witnesses()[1].index == 1U);
    CHECK(c.witnesses()[0].normal.x == doctest::Approx(-c.witnesses()[1].normal.x));
}

TEST_CASE("Clearance refuses undefined queries")
{
    const PolylineProfile profile = squareStep();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK(query(profile, 0.9, 1.0, 0.2).status == ClearanceStatus::OutsideDomain);
    CHECK(query(profile, 0.1, 1e200, 0.05).status == ClearanceStatus::OutsideDomain);
    CHECK(query(profile, 0.1, nan, 0.05).status == ClearanceStatus::OutsideDomain);
    CHECK(query(profile, 0.1, 1.0, nan).status == ClearanceStatus::OutsideDomain);
    CHECK(query(profile, 0.1, 0.1, 0.05).status == ClearanceStatus::OnBoundary);
    CHECK(query(profile, 0.1, 0.05, 0.01).status == ClearanceStatus::CentreInside);
}

TEST_CASE("Pruned search agrees with brute force on random terrain")
{
    std::mt19937_64 rng(20261015U);
    std::uniform_real_distribution<double> step(0.0, 0.05);
    std::uniform_real_distribution<double> rise(-0.04, 0.04);
    std::bernoulli_distribution vertical(0.1);
    std::vector<std::pair<double, double>> raw{{-5.0, 0.0}};
    while (raw.back().first < 5.0)
    {
        const auto [lx, ly] = raw.back();
        const bool face = vertical(rng) && (raw.size() < 2U || raw[raw.size() - 2U].first != lx);
        raw.push_back(face ? std::pair{lx, ly + 0.05} : std::pair{lx + 0.001 + step(rng), ly + rise(rng)});
    }
    std::vector<ProfilePoint> points;
    for (const auto& [px, py] : raw)
    {
        points.push_back(at(px, py));
    }
    const auto profile = *PolylineProfile::create(points);

    std::uniform_real_distribution<double> x(-4.0, 4.0);
    std::uniform_real_distribution<double> lift(0.0, 1.0);
    for (int trial = 0; trial < 2000; ++trial)
    {
        const double cx = x(rng);
        const double cy = profile.surfaceHeight(metres(cx)).si() + 1e-3 + lift(rng);
        double brute = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i + 1U < raw.size(); ++i)
        {
            const auto [ax, ay] = raw[i];
            const auto [bx, by] = raw[i + 1U];
            const double dx = bx - ax;
            const double dy = by - ay;
            const double t = std::clamp(((cx - ax) * dx + (cy - ay) * dy) / (dx * dx + dy * dy), 0.0, 1.0);
            brute = std::min(brute, std::hypot(cx - (ax + t * dx), cy - (ay + t * dy)));
        }
        const Measured c = query(profile, cx, cy, 0.5);
        REQUIRE(c.status == ClearanceStatus::Available);
        CHECK(c.distance.si() == doctest::Approx(brute).epsilon(1e-14));
        CHECK_FALSE(c.witnesses().empty());
    }
}

TEST_CASE("Non-finite or out-of-range query values are outside the domain, never available")
{
    const auto flat = *PolylineProfile::create(metresPoints({{-1.0, 0.0}, {1.0, 0.0}}));
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const auto status = [&](double cx, double cy, double radius, double tie) {
        return measureCircleClearance(flat, {{at(cx, cy), metres(radius)}, metres(tie)}, {}).status;
    };
    REQUIRE(status(0.0, 0.5, 0.1, 1e-10) == ClearanceStatus::Available);
    CHECK(status(nan, 0.5, 0.1, 1e-10) == ClearanceStatus::OutsideDomain);
    CHECK(status(inf, 0.5, 0.1, 1e-10) == ClearanceStatus::OutsideDomain);
    CHECK(status(-inf, 0.5, 0.1, 1e-10) == ClearanceStatus::OutsideDomain);
    CHECK(status(1e101, 0.5, 0.1, 1e-10) == ClearanceStatus::OutsideDomain);
    CHECK(status(0.0, nan, 0.1, 1e-10) == ClearanceStatus::OutsideDomain);
    CHECK(status(0.0, inf, 0.1, 1e-10) == ClearanceStatus::OutsideDomain);
    CHECK(status(0.0, 0.5, nan, 1e-10) == ClearanceStatus::OutsideDomain);
    CHECK(status(0.0, 0.5, inf, 1e-10) == ClearanceStatus::OutsideDomain);
    CHECK(status(0.0, 0.5, 0.1, nan) == ClearanceStatus::OutsideDomain);
    CHECK(status(0.0, 0.5, 0.1, inf) == ClearanceStatus::OutsideDomain);
    CHECK(status(0.0, 0.5, 0.1, -1e-10) == ClearanceStatus::OutsideDomain);
    CHECK(status(0.0, 0.5, 0.1, 0.0) == ClearanceStatus::Available); // zero tolerance is allowed
}

TEST_CASE("A full witness buffer keeps the first ties and says it overflowed")
{
    const auto valley = *PolylineProfile::create(metresPoints({{-1.0, 1.0}, {0.0, 0.0}, {1.0, 1.0}}));
    const Measured all = query(valley, 0.0, 0.5, 0.1);
    REQUIRE(all.status == ClearanceStatus::Available);
    CHECK_FALSE(all.overflowed);
    REQUIRE(all.witnessCount == 2U);

    const Measured one = query(valley, 0.0, 0.5, 0.1, 1U);
    CHECK(one.overflowed);
    REQUIRE(one.witnessCount == 1U);
    CHECK(one.witnesses()[0].index == all.witnesses()[0].index); // the first tie, in profile order
    CHECK(one.distance.si() == all.distance.si());                // the measurement itself is unaffected
    CHECK(one.buffer[1].index == 0U);                              // nothing written past the capacity

    const Measured none = query(valley, 0.0, 0.5, 0.1, 0U);       // distance only
    CHECK(none.witnessCount == 0U);
    CHECK(none.overflowed);
    CHECK(none.gap.si() == all.gap.si());
}

TEST_CASE("Witness kinds: a corner is a vertex, a face interior is an edge, never a blend")
{
    const Measured corner = query(squareStep(), -0.1, 0.25, 0.1);
    REQUIRE(corner.witnessCount == 1U);
    CHECK(corner.witnesses()[0].kind == WitnessKind::Vertex);
    const Measured face = query(squareStep(), 0.1, 0.5, 0.1);
    REQUIRE(face.witnessCount == 1U);
    CHECK(face.witnesses()[0].kind == WitnessKind::Edge);
    CHECK(face.witnesses()[0].normal.y == doctest::Approx(1.0));
    // A witness normal is a unit vector pointing from the surface point to the centre.
    for (const Witness2& w : {corner.witnesses()[0], face.witnesses()[0]})
    {
        CHECK(std::hypot(w.normal.x, w.normal.y) == doctest::Approx(1.0));
    }
}
