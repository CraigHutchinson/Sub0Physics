#include <algorithm>
#include <array>
#include <cmath>
#include <doctest.h>
#include <limits>
#include <numbers>
#include <random>
#include <span>
#include <sub0math/units.hpp>
#include <vector>

#include "sub0physics/geometry.hpp"
#include "test_profiles.hpp"

using namespace sub0physics;
using sub0math::metres;
using sub0math::Vec2d;
using sub0math::f64::Length;
using test_profiles::at;
using test_profiles::metresPoints;

namespace
{

constexpr double Infinity = std::numeric_limits<double>::infinity();

// ---- Straightforward references: the formulations PreparedProfile replaces, written plainly over the points ----

// Scaled, compensated norm (V8's Math.hypot for two arguments): the arithmetic the prepared normals use.
double scaledNorm(double a, double b)
{
    const double largest = std::max(std::abs(a), std::abs(b));
    if (largest == 0.0)
    {
        return 0.0;
    }
    double sum = 0.0;
    double compensation = 0.0;
    for (const double v : {a, b})
    {
        const double n = std::abs(v) / largest;
        const double summand = n * n - compensation;
        const double preliminary = sum + summand;
        compensation = (preliminary - sum) - summand;
        sum = preliminary;
    }
    return std::sqrt(sum) * largest;
}

struct Raw
{
    std::vector<double> x, y;

    explicit Raw(const PolylineProfile& profile)
    {
        for (const auto& p : profile.points())
        {
            x.push_back(p.x.si());
            y.push_back(p.y.si());
        }
    }

    [[nodiscard]] std::size_t segments() const
    {
        return x.size() - 1U;
    }
};

struct RefRay
{
    bool hit = false;
    bool ambiguous = false;
    double distance = Infinity;
    std::size_t feature = 0;
};

RefRay referenceRay(const Raw& line, double cx, double cz, double ox, double oz, double tie)
{
    RefRay out;
    double hitNx = 0.0;
    double hitNz = 0.0;
    for (std::size_t j = 0; j < line.segments(); ++j)
    {
        const double dx = line.x[j + 1U] - line.x[j];
        const double dz = line.y[j + 1U] - line.y[j];
        const double length = scaledNorm(dx, dz);
        const double nx = -dz / length;
        const double nz = dx / length;
        const double ax = line.x[j] - cx;
        const double az = line.y[j] - cz;
        if (nx * ox + nz * oz >= -1e-12)
        {
            continue;
        }
        const double den = ox * dz - oz * dx;
        const double r = (ax * dz - az * dx) / den;
        const double t = (ax * oz - az * ox) / den;
        if (r <= 0.0 || t < -tie || t > 1.0 + tie || r > out.distance + tie)
        {
            continue;
        }
        if (std::abs(r - out.distance) < tie)
        {
            if (scaledNorm(hitNx - nx, hitNz - nz) > 1e-6)
            {
                out.ambiguous = true;
            }
            continue;
        }
        out.distance = r;
        out.feature = j;
        out.hit = true;
        out.ambiguous = false;
        hitNx = nx;
        hitNz = nz;
    }
    return out;
}

// A different formulation of the same first hit, for independence: intersect the ray with each segment by solving
// ray(s) = a + u (b - a) for (s, u) with Cramer's rule, and keep the nearest hit on a segment that faces the origin.
double independentRayDistance(const Raw& line, double cx, double cz, double ox, double oz)
{
    double best = Infinity;
    for (std::size_t j = 0; j < line.segments(); ++j)
    {
        const double ex = line.x[j + 1U] - line.x[j];
        const double ez = line.y[j + 1U] - line.y[j];
        // Solid is below the segment direction (left normal points up for forward x): front-facing means the origin
        // is on the left of the segment and the ray heads towards its line.
        const double side = ex * (cz - line.y[j]) - ez * (cx - line.x[j]);
        const double det = ox * (-ez) - oz * (-ex);
        if (side <= 0.0 || std::abs(det) < 1e-14)
        {
            continue;
        }
        const double qx = line.x[j] - cx;
        const double qz = line.y[j] - cz;
        const double s = (qx * (-ez) - qz * (-ex)) / det;
        const double u = (ox * qz - oz * qx) / det;
        if (s > 0.0 && u >= -1e-9 && u <= 1.0 + 1e-9)
        {
            best = std::min(best, s);
        }
    }
    return best;
}

bool referenceOnSurface(const Raw& line, double px, double pz, double tol)
{
    for (std::size_t j = 0; j < line.segments(); ++j)
    {
        const double dx = line.x[j + 1U] - line.x[j];
        const double dz = line.y[j + 1U] - line.y[j];
        const double t = ((px - line.x[j]) * dx + (pz - line.y[j]) * dz) / (dx * dx + dz * dz);
        if (t >= -tol && t <= 1.0 + tol && std::abs((px - line.x[j]) * dz - (pz - line.y[j]) * dx) / scaledNorm(dx, dz) < tol)
        {
            return true;
        }
    }
    return false;
}

// The original contour check: every chord scanned against every vertex, sorted, every sample against every segment.
double referenceIntrusion(const PolylineProfile& profile, const Raw& line, std::span<const sub0math::Point2<Length>> outline,
                          double tol)
{
    double maximum = 0.0;
    for (std::size_t i = 0; i < outline.size(); ++i)
    {
        const auto a = outline[i];
        const auto b = outline[(i + 1U) % outline.size()];
        const double ax = a.x.si();
        const double az = a.y.si();
        const double dx = b.x.si() - ax;
        const double dz = b.y.si() - az;
        std::vector<double> cuts{0.0, 1.0};
        if (std::abs(dx) > 1e-14)
        {
            for (const double vertexX : line.x)
            {
                const double t = (vertexX - ax) / dx;
                if (t > 0.0 && t < 1.0)
                {
                    cuts.push_back(t);
                }
            }
        }
        std::sort(cuts.begin(), cuts.end());
        for (std::size_t j = 1; j < cuts.size(); ++j)
        {
            for (const double f : {1e-8, 0.5, 1.0 - 1e-8})
            {
                const double t = cuts[j - 1U] + (cuts[j] - cuts[j - 1U]) * f;
                const double x = ax + t * dx;
                const double z = az + t * dz;
                if (!referenceOnSurface(line, x, z, tol))
                {
                    maximum = std::max(maximum, profile.surfaceHeight(Length::fromSi(x)).si() - z);
                }
            }
        }
    }
    return std::max(0.0, maximum);
}

struct RefSweep
{
    bool contact = false;
    double fraction = 1.0;
};

// The swept-circle roots over every vertex and segment, with no window.
RefSweep referenceSweep(const Raw& line, double cx, double cz, double dx, double dz, double r, double tol)
{
    double hit = Infinity;
    const auto candidate = [&](double t) {
        if (t >= -tol && t <= 1.0 + tol)
        {
            hit = std::min(hit, std::max(0.0, std::min(1.0, t)));
        }
    };
    for (std::size_t v = 0; v < line.x.size(); ++v)
    {
        const double x = cx - line.x[v];
        const double z = cz - line.y[v];
        const double a = dx * dx + dz * dz;
        if (a == 0.0)
        {
            break;
        }
        const double b = 2.0 * (x * dx + z * dz);
        const double c = x * x + z * z - r * r;
        const double discriminant = b * b - 4.0 * a * c;
        const double scale = b * b + std::abs(4.0 * a * c);
        if (discriminant < -32.0 * std::numeric_limits<double>::epsilon() * scale)
        {
            continue;
        }
        const double root = std::sqrt(std::max(0.0, discriminant));
        const double q = -0.5 * (b + (b >= 0.0 ? root : -root));
        if (q == 0.0)
        {
            candidate(-b / (2.0 * a));
        }
        else
        {
            candidate(q / a);
            candidate(c / q);
        }
    }
    for (std::size_t i = 0; i < line.segments(); ++i)
    {
        const double sx = line.x[i + 1U] - line.x[i];
        const double sz = line.y[i + 1U] - line.y[i];
        const double length = scaledNorm(sx, sz);
        const double nx = -sz / length;
        const double nz = sx / length;
        const double d0 = (cx - line.x[i]) * nx + (cz - line.y[i]) * nz;
        const double dv = dx * nx + dz * nz;
        if (dv == 0.0)
        {
            continue;
        }
        for (const double sign : {-1.0, 1.0})
        {
            const double t = (sign * r - d0) / dv;
            if (t < -tol || t > 1.0 + tol)
            {
                continue;
            }
            const double projection = ((cx + t * dx - line.x[i]) * sx + (cz + t * dz - line.y[i]) * sz) / (length * length);
            if (projection >= 0.0 && projection <= 1.0)
            {
                candidate(t);
            }
        }
    }
    return {hit < Infinity, hit < Infinity ? hit : 1.0};
}

// ---- Random profiles ----

// A height field with random spacing, random heights and occasional vertical faces (never reversing).
PolylineProfile randomProfile(std::mt19937_64& rng, std::size_t points, bool allowVertical)
{
    std::uniform_real_distribution<double> step(0.002, 0.15);
    std::uniform_real_distribution<double> height(-0.15, 0.25);
    std::bernoulli_distribution vertical(0.2);
    std::vector<ProfilePoint> out;
    double x = -1.0;
    double y = height(rng);
    out.push_back(at(x, y));
    bool lastVertical = false;
    while (out.size() < points)
    {
        if (allowVertical && !lastVertical && vertical(rng))
        {
            double ny = height(rng);
            while (std::abs(ny - y) < 0.01)
            {
                ny = height(rng);
            }
            y = ny;
            lastVertical = true;
        }
        else
        {
            x += step(rng);
            y = height(rng);
            lastVertical = false;
        }
        out.push_back(at(x, y));
    }
    auto profile = PolylineProfile::create(std::move(out));
    REQUIRE(profile.has_value());
    return *profile;
}

Length m(double v)
{
    return metres(v);
}

} // namespace

TEST_CASE("PreparedProfile keeps the geometry of the profile it was built from")
{
    const auto profile = *PolylineProfile::create(metresPoints({{-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.1}, {0.2, 0.1}, {0.2, 0.0}, {1.0, 0.0}}));
    const PreparedProfile prepared(profile);
    REQUIRE(prepared.vertexCount() == 6U);
    REQUIRE(prepared.segmentCount() == 5U);
    CHECK(prepared.minX().si() == -1.0);
    CHECK(prepared.maxX().si() == 1.0);
    // Segment 1 is the riser: it faces back (-x) with the solid below and behind it.
    CHECK(prepared.segmentNormalX()[1] == -1.0);
    CHECK(prepared.segmentNormalY()[1] == 0.0);
    // Segment 2 is the top, facing up.
    CHECK(prepared.segmentNormalX()[2] == 0.0);
    CHECK(prepared.segmentNormalY()[2] == 1.0);
    CHECK(prepared.segmentSumSquares()[2] == 0.2 * 0.2);
}

TEST_CASE("Prepared surface height is bit-identical to the profile's")
{
    std::mt19937_64 rng(11);
    for (int round = 0; round < 20; ++round)
    {
        const auto profile = randomProfile(rng, 5U + static_cast<std::size_t>(round) * 3U, true);
        const PreparedProfile prepared(profile);
        std::uniform_real_distribution<double> xs(profile.minX().si() - 0.1, profile.maxX().si() + 0.1);
        for (int k = 0; k < 300; ++k)
        {
            const double x = xs(rng);
            CHECK(prepared.surfaceHeight(m(x)).si() == profile.surfaceHeight(m(x)).si());
        }
        for (const auto& p : profile.points())
        {
            CHECK(prepared.surfaceHeight(p.x).si() == profile.surfaceHeight(p.x).si());
        }
    }
    const auto flat = *PolylineProfile::create(metresPoints({{-1.0, 0.0}, {1.0, 0.0}}));
    CHECK(PreparedProfile(flat).surfaceHeight(m(std::nan(""))).si() == -Infinity);
}

TEST_CASE("A ray from above flat ground hits it at the height")
{
    const auto flat = *PolylineProfile::create(metresPoints({{-1.0, 0.0}, {1.0, 0.0}}));
    const PreparedProfile prepared(flat);
    const RayHit down = prepared.castRay(at(0.2, 0.5), Vec2d{0.0, -1.0}, 1e-10);
    REQUIRE(down.status == RayStatus::Hit);
    CHECK(down.distance == doctest::Approx(0.5));
    CHECK_FALSE(down.ambiguous);
    CHECK(down.witness.kind == WitnessKind::Edge);
    CHECK(down.witness.index == 0U);
    CHECK(down.witness.normal.y == 1.0);
    CHECK(down.witness.point.y.si() == doctest::Approx(0.0));
    CHECK(down.witness.point.x.si() == doctest::Approx(0.2));

    const RayHit up = prepared.castRay(at(0.2, 0.5), Vec2d{0.0, 1.0}, 1e-10);
    CHECK(up.status == RayStatus::Miss);
    CHECK(up.distance == Infinity);
    const RayHit sideways = prepared.castRay(at(0.2, 0.5), Vec2d{1.0, 0.0}, 1e-10);
    CHECK(sideways.status == RayStatus::Miss);
    // Past the end of the profile: the parameter is outside the segment.
    CHECK(prepared.castRay(at(2.0, 0.5), Vec2d{0.0, -1.0}, 1e-10).status == RayStatus::Miss);
}

TEST_CASE("A ray into a valley corner is ambiguous; one into a face is not")
{
    // A V: two faces meet at (0, 0) with normals leaning apart.
    const auto valley = *PolylineProfile::create(metresPoints({{-1.0, 1.0}, {0.0, 0.0}, {1.0, 1.0}}));
    const PreparedProfile prepared(valley);
    const RayHit corner = prepared.castRay(at(0.0, 0.5), Vec2d{0.0, -1.0}, 1e-10);
    REQUIRE(corner.status == RayStatus::Hit);
    CHECK(corner.distance == doctest::Approx(0.5));
    CHECK(corner.ambiguous);
    const RayHit face = prepared.castRay(at(0.3, 0.5), Vec2d{0.0, -1.0}, 1e-10);
    REQUIRE(face.status == RayStatus::Hit);
    CHECK_FALSE(face.ambiguous);
    CHECK(face.witness.index == 1U);
}

TEST_CASE("A ray cast from a bound origin equals the one-off cast and the reference")
{
    std::mt19937_64 rng(5);
    RayOrigin bound;
    for (int round = 0; round < 30; ++round)
    {
        const auto profile = randomProfile(rng, 4U + static_cast<std::size_t>(round) * 2U, true);
        const Raw raw(profile);
        const PreparedProfile prepared(profile);
        std::uniform_real_distribution<double> ox(-0.5, 0.5);
        std::uniform_real_distribution<double> oz(0.35, 0.6);
        const double cx = ox(rng);
        const double cz = oz(rng);
        prepared.bindRayOrigin(at(cx, cz), bound);
        for (int k = 0; k < 200; ++k)
        {
            const double angle = -std::numbers::pi + (static_cast<double>(k) + 0.5) * 2.0 * std::numbers::pi / 200.0;
            const Vec2d direction{std::sin(angle), -std::cos(angle)};
            const RayHit viaBound = prepared.castRay(bound, direction, 1e-10);
            const RayHit oneOff = prepared.castRay(at(cx, cz), direction, 1e-10);
            const RefRay ref = referenceRay(raw, cx, cz, direction.x, direction.y, 1e-10);
            REQUIRE((viaBound.status == RayStatus::Hit) == ref.hit);
            CHECK(viaBound.status == oneOff.status);
            CHECK(viaBound.distance == ref.distance);
            CHECK(oneOff.distance == ref.distance);
            CHECK(viaBound.ambiguous == ref.ambiguous);
            if (ref.hit)
            {
                CHECK(viaBound.witness.index == ref.feature);
                CHECK(oneOff.witness.index == ref.feature);
                // Independent formulation of the same first hit.
                CHECK(independentRayDistance(raw, cx, cz, direction.x, direction.y) == doctest::Approx(ref.distance).epsilon(1e-9));
            }
            else
            {
                CHECK(independentRayDistance(raw, cx, cz, direction.x, direction.y) == Infinity);
            }
        }
    }
}

TEST_CASE("Ray casts reject invalid queries")
{
    const auto flat = *PolylineProfile::create(metresPoints({{-1.0, 0.0}, {1.0, 0.0}}));
    const PreparedProfile prepared(flat);
    CHECK(prepared.castRay(at(0.0, 0.5), Vec2d{0.0, -1.0}, -1.0).status == RayStatus::InvalidQuery);
    CHECK(prepared.castRay(at(0.0, 0.5), Vec2d{0.0, -1.0}, std::nan("")).status == RayStatus::InvalidQuery);
    CHECK(prepared.castRay(at(std::nan(""), 0.5), Vec2d{0.0, -1.0}, 1e-10).status == RayStatus::InvalidQuery);
    CHECK(prepared.castRay(at(0.0, 0.5), Vec2d{Infinity, -1.0}, 1e-10).status == RayStatus::InvalidQuery);
    RayOrigin unbound;
    CHECK(prepared.castRay(unbound, Vec2d{0.0, -1.0}, 1e-10).status == RayStatus::InvalidQuery);
}

TEST_CASE("Sweeping a circle finds the first touch of a flat ground")
{
    const auto flat = *PolylineProfile::create(metresPoints({{-1.0, 0.0}, {1.0, 0.0}}));
    const PreparedProfile prepared(flat);
    const CircleSweep down = prepared.sweepCircle(at(0.0, 0.5), {m(0.0), m(-0.4)}, m(0.2));
    REQUIRE(down.status == CircleSweepStatus::Contact);
    CHECK(down.fraction == doctest::Approx(0.75));
    CHECK(down.witness.kind == WitnessKind::Edge);
    CHECK(down.witness.normal.y == 1.0);
    CHECK(down.witness.point.y.si() == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(prepared.sweepCircle(at(0.0, 0.5), {m(0.5), m(0.0)}, m(0.2)).status == CircleSweepStatus::Clear);
    CHECK(prepared.sweepCircle(at(0.0, 0.5), {m(0.0), m(0.0)}, m(0.2)).status == CircleSweepStatus::Clear);
    CHECK(prepared.sweepCircle(at(0.0, 0.5), {m(0.0), m(-0.4)}, m(-0.2)).status == CircleSweepStatus::InvalidQuery);
    CHECK(prepared.sweepCircle(at(0.0, 0.5), {m(0.0), m(std::nan(""))}, m(0.2)).status == CircleSweepStatus::InvalidQuery);
}

TEST_CASE("A circle swept past a square corner touches the vertex")
{
    const auto step = *PolylineProfile::create(metresPoints({{-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.1}, {0.2, 0.1}, {0.2, 0.0}, {1.0, 0.0}}));
    const PreparedProfile prepared(step);
    // Moving right with the centre 0.05 above the block's top (radius 0.1): the circle meets the top-left corner
    // (0, 0.1) before the riser face, when the centre is sqrt(0.1^2 - 0.05^2) to its left.
    const CircleSweep sweep = prepared.sweepCircle(at(-0.5, 0.15), {m(1.0), m(0.0)}, m(0.1));
    REQUIRE(sweep.status == CircleSweepStatus::Contact);
    const double offset = std::sqrt(0.1 * 0.1 - 0.05 * 0.05);
    CHECK(sweep.fraction == doctest::Approx(0.5 - offset));
    CHECK(sweep.witness.kind == WitnessKind::Vertex);
    CHECK(sweep.witness.index == 2U);
    CHECK(sweep.witness.normal.x == doctest::Approx(-offset / 0.1));
    CHECK(sweep.witness.normal.y == doctest::Approx(0.5));
    // At the riser's height the face itself is met first.
    const CircleSweep face = prepared.sweepCircle(at(-0.5, 0.05), {m(1.0), m(0.0)}, m(0.1));
    REQUIRE(face.status == CircleSweepStatus::Contact);
    CHECK(face.fraction == doctest::Approx(0.4));
    CHECK(face.witness.kind == WitnessKind::Edge);
    CHECK(face.witness.index == 1U);
    CHECK(face.witness.normal.x == -1.0);
}

TEST_CASE("The windowed sweep equals the full scan and a sampled search")
{
    std::mt19937_64 rng(23);
    int contacts = 0;
    int clear = 0;
    for (int round = 0; round < 40; ++round)
    {
        const auto profile = randomProfile(rng, 6U + static_cast<std::size_t>(round) * 4U, round % 2 == 0);
        const Raw raw(profile);
        const PreparedProfile prepared(profile);
        const double radius = 0.04 + 0.002 * round;
        std::uniform_real_distribution<double> sx(profile.minX().si() + radius + 0.01, profile.maxX().si() - radius - 0.01);
        std::uniform_real_distribution<double> dxs(-0.8, 0.8);
        std::uniform_real_distribution<double> dzs(-0.5, 0.2);
        for (int k = 0; k < 40; ++k)
        {
            const double cx = sx(rng);
            const double surface = profile.surfaceHeight(m(cx)).si();
            const double cz = surface + radius + 0.02 + 0.2 * std::generate_canonical<double, 53>(rng);
            const double dx = dxs(rng);
            const double dz = dzs(rng);
            std::array<Witness2, DefaultWitnessCapacity> buffer;
            const auto startClearance = measureCircleClearance(profile, {{at(cx, cz), m(radius)}, m(1e-10)}, buffer);
            if (startClearance.status != ClearanceStatus::Available || startClearance.gap.si() <= 1e-6)
            {
                continue;
            }
            const double endX = cx + dx;
            if (endX - radius < profile.minX().si() || endX + radius > profile.maxX().si())
            {
                continue;
            }
            const CircleSweep sweep = prepared.sweepCircle(at(cx, cz), {m(dx), m(dz)}, m(radius));
            const RefSweep ref = referenceSweep(raw, cx, cz, dx, dz, radius, 1e-10);
            REQUIRE((sweep.status == CircleSweepStatus::Contact) == ref.contact);
            CHECK(sweep.fraction == ref.fraction);
            // Sampled search: the clearance stays positive before the contact and is zero at it.
            constexpr int Samples = 2000;
            double firstTouch = 1.0;
            bool touched = false;
            for (int s = 0; s <= Samples; ++s)
            {
                const double f = static_cast<double>(s) / Samples;
                const auto c = measureCircleClearance(profile, {{at(cx + f * dx, cz + f * dz), m(radius)}, m(1e-10)}, buffer);
                if (c.status != ClearanceStatus::Available || c.gap.si() <= 0.0)
                {
                    firstTouch = f;
                    touched = true;
                    break;
                }
            }
            if (sweep.status == CircleSweepStatus::Contact)
            {
                ++contacts;
                const auto at0 = measureCircleClearance(profile, {{at(cx + sweep.fraction * dx, cz + sweep.fraction * dz), m(radius)}, m(1e-10)}, buffer);
                if (at0.status == ClearanceStatus::Available)
                {
                    CHECK(std::abs(at0.gap.si()) < 1e-9);
                }
                // The sampled search can only be late, by less than a sample, unless it tunnelled through a thin feature.
                if (touched)
                {
                    CHECK(firstTouch + 1e-12 >= sweep.fraction - 1e-9);
                    CHECK(firstTouch - sweep.fraction < 1.5 / Samples + 1e-9);
                }
            }
            else
            {
                ++clear;
                CHECK_FALSE(touched);
            }
        }
    }
    CHECK(contacts > 100);
    CHECK(clear > 50);
}

TEST_CASE("onSurface agrees with the full scan, on and off the profile")
{
    std::mt19937_64 rng(31);
    int on = 0;
    int off = 0;
    for (int round = 0; round < 20; ++round)
    {
        const auto profile = randomProfile(rng, 8U + static_cast<std::size_t>(round) * 4U, true);
        const Raw raw(profile);
        const PreparedProfile prepared(profile);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        std::uniform_real_distribution<double> offset(-3e-10, 3e-10);
        for (int k = 0; k < 400; ++k)
        {
            const std::size_t seg = static_cast<std::size_t>(unit(rng) * static_cast<double>(raw.segments()));
            const double t = k % 7 == 0 ? (k % 2 == 0 ? -1e-11 : 1.0 + 1e-11) : unit(rng);
            const double px = raw.x[seg] + t * (raw.x[seg + 1U] - raw.x[seg]) + (k % 3 == 0 ? offset(rng) : 0.0);
            const double pz = raw.y[seg] + t * (raw.y[seg + 1U] - raw.y[seg]) + (k % 3 == 1 ? offset(rng) : 0.0);
            const bool expected = referenceOnSurface(raw, px, pz, 1e-10);
            CHECK(prepared.onSurface(at(px, pz), 1e-10) == expected);
            (expected ? on : off)++;
        }
        // Far from the profile.
        CHECK_FALSE(prepared.onSurface(at(0.0, 5.0), 1e-10));
    }
    CHECK(on > 500);
    CHECK(off > 500);
    const auto flat = *PolylineProfile::create(metresPoints({{-1.0, 0.0}, {1.0, 0.0}}));
    const PreparedProfile prepared(flat);
    CHECK(prepared.onSurface(at(0.3, 0.0), 1e-10));
    CHECK_FALSE(prepared.onSurface(at(0.3, 1e-9), 1e-10));
    CHECK_FALSE(prepared.onSurface(at(std::nan(""), 0.0), 1e-10));
    CHECK_FALSE(prepared.onSurface(at(0.3, 0.0), -1.0));
}

TEST_CASE("onSurface finds a point on a profile with one very long segment")
{
    // The overshoot allowed beyond a segment's end is a fraction of its length, so a long segment reaches far in x.
    const auto longProfile = *PolylineProfile::create(metresPoints({{0.0, 0.0}, {1.0, 0.0}, {1.0001, 0.0001}, {1000.0, 0.0001}}));
    const Raw raw(longProfile);
    const PreparedProfile prepared(longProfile);
    for (const double x : {-5e-8, -1e-8, 1.0 + 5e-11, 1000.0 + 5e-8, 1000.0 + 5e-6})
    {
        CHECK(prepared.onSurface(at(x, 0.0001), 1e-10) == referenceOnSurface(raw, x, 0.0001, 1e-10));
        CHECK(prepared.onSurface(at(x, 0.0), 1e-10) == referenceOnSurface(raw, x, 0.0, 1e-10));
    }
}

TEST_CASE("Contour intrusion on known shapes")
{
    const auto flat = *PolylineProfile::create(metresPoints({{-1.0, 0.0}, {1.0, 0.0}}));
    const PreparedProfile prepared(flat);
    const std::array<sub0math::Point2<Length>, 3> triangle{{at(-0.1, 0.02), at(0.1, 0.02), at(0.0, -0.03)}};
    const ContourIntrusion cut = prepared.contourIntrusion(triangle);
    REQUIRE(cut.status == ContourStatus::Intrudes);
    CHECK(cut.depth.si() == doctest::Approx(0.03).epsilon(1e-6));
    CHECK(cut.witness.kind == WitnessKind::Edge);
    CHECK(cut.witness.index == 0U);
    CHECK(cut.witness.normal.y == 1.0);
    CHECK(cut.witness.point.y.si() == 0.0);

    const std::array<sub0math::Point2<Length>, 3> above{{at(-0.1, 0.02), at(0.1, 0.02), at(0.0, 0.1)}};
    const ContourIntrusion clear = prepared.contourIntrusion(above);
    CHECK(clear.status == ContourStatus::Clear);
    CHECK(clear.depth.si() == 0.0);

    // A contour lying along the surface is not an intrusion.
    const std::array<sub0math::Point2<Length>, 3> along{{at(-0.1, 0.0), at(0.1, 0.0), at(0.0, 0.05)}};
    CHECK(prepared.contourIntrusion(along).status == ContourStatus::Clear);
    CHECK(prepared.contourIntrusion({}).status == ContourStatus::Clear);

    const std::array<sub0math::Point2<Length>, 3> bad{{at(-0.1, 0.02), at(std::nan(""), 0.02), at(0.0, 0.1)}};
    CHECK(prepared.contourIntrusion(bad).status == ContourStatus::InvalidQuery);
    CHECK(prepared.contourIntrusion(above, -1.0).status == ContourStatus::InvalidQuery);
}

TEST_CASE("A chord cutting a square corner is caught by the one-sided samples")
{
    const auto step = *PolylineProfile::create(metresPoints({{-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.1}, {0.2, 0.1}, {0.2, 0.0}, {1.0, 0.0}}));
    const PreparedProfile prepared(step);
    // The chord from (-0.05, 0.12) to (0.05, 0.02) passes through the block's top-left corner region.
    const std::array<sub0math::Point2<Length>, 3> contour{{at(-0.05, 0.12), at(0.05, 0.02), at(-0.05, 0.02)}};
    const ContourIntrusion hit = prepared.contourIntrusion(contour);
    CHECK(hit.status == ContourStatus::Intrudes);
    const Raw raw(step);
    CHECK(hit.depth.si() == referenceIntrusion(step, raw, contour, 1e-10));
}

TEST_CASE("Contour intrusion equals the original full-scan algorithm bit for bit")
{
    std::mt19937_64 rng(47);
    int intruding = 0;
    int clear = 0;
    for (int round = 0; round < 40; ++round)
    {
        const auto profile = randomProfile(rng, 5U + static_cast<std::size_t>(round) * 3U, round % 3 != 0);
        const Raw raw(profile);
        const PreparedProfile prepared(profile);
        std::uniform_real_distribution<double> cx(-0.5, 0.5);
        std::uniform_real_distribution<double> cz(0.0, 0.3);
        std::uniform_real_distribution<double> radius(0.02, 0.4);
        for (int k = 0; k < 15; ++k)
        {
            const std::size_t n = 8U + static_cast<std::size_t>(k) * 17U;
            const double centreX = cx(rng);
            const double centreZ = cz(rng);
            const double rad = radius(rng);
            std::vector<sub0math::Point2<Length>> outline;
            for (std::size_t i = 0; i < n; ++i)
            {
                const double angle = static_cast<double>(i) * 2.0 * std::numbers::pi / static_cast<double>(n);
                outline.push_back(at(centreX + rad * std::cos(angle), centreZ + rad * std::sin(angle)));
            }
            const double expected = referenceIntrusion(profile, raw, outline, 1e-10);
            const ContourIntrusion result = prepared.contourIntrusion(outline);
            CHECK(result.depth.si() == expected);
            CHECK((result.status == ContourStatus::Intrudes) == (expected > 0.0));
            (expected > 0.0 ? intruding : clear)++;
            if (result.status == ContourStatus::Intrudes)
            {
                // The witness is the surface straight above the deepest sample, on the segment that set the height.
                CHECK(result.witness.kind == WitnessKind::Edge);
                CHECK(result.witness.point.y.si() == profile.surfaceHeight(result.witness.point.x).si());
            }
        }
    }
    CHECK(intruding > 100);
    CHECK(clear > 100);
}

TEST_CASE("The prepared clearance measures the same circle as the unprepared one")
{
    std::mt19937_64 rng(59);
    int available = 0;
    for (int round = 0; round < 30; ++round)
    {
        const auto profile = randomProfile(rng, 4U + static_cast<std::size_t>(round) * 5U, round % 2 == 0);
        const PreparedProfile prepared(profile);
        std::uniform_real_distribution<double> xs(profile.minX().si() - 0.05, profile.maxX().si() + 0.05);
        std::uniform_real_distribution<double> ys(-0.4, 0.6);
        std::uniform_real_distribution<double> radius(0.01, 0.5);
        for (int k = 0; k < 150; ++k)
        {
            const CircleQuery query{{at(xs(rng), ys(rng)), m(radius(rng))}, m(1e-10)};
            std::array<Witness2, DefaultWitnessCapacity> a;
            std::array<Witness2, DefaultWitnessCapacity> b;
            const CircleClearance plain = measureCircleClearance(profile, query, a);
            const CircleClearance fast = measureCircleClearance(prepared, query, b);
            REQUIRE(fast.status == plain.status);
            if (plain.status != ClearanceStatus::Available)
            {
                continue;
            }
            ++available;
            CHECK(fast.distance.si() == plain.distance.si());
            CHECK(fast.gap.si() == plain.gap.si());
            REQUIRE(fast.witnessCount == plain.witnessCount);
            CHECK(fast.overflowed == plain.overflowed);
            for (std::size_t i = 0; i < plain.witnessCount; ++i)
            {
                CHECK(b[i].kind == a[i].kind);
                CHECK(b[i].index == a[i].index);
                CHECK(b[i].point.x.si() == a[i].point.x.si());
                CHECK(b[i].point.y.si() == a[i].point.y.si());
                CHECK(b[i].normal.x == a[i].normal.x);
                CHECK(b[i].normal.y == a[i].normal.y);
            }
        }
    }
    CHECK(available > 1000);
}

TEST_CASE("The prepared clearance keeps the unprepared statuses, ties and buffer rules")
{
    const auto valley = *PolylineProfile::create(metresPoints({{-1.0, 1.0}, {0.0, 0.0}, {1.0, 1.0}}));
    const PreparedProfile prepared(valley);
    std::array<Witness2, DefaultWitnessCapacity> buffer;
    const CircleClearance tie = measureCircleClearance(prepared, {{at(0.0, 0.5), m(0.1)}, m(1e-10)}, buffer);
    REQUIRE(tie.status == ClearanceStatus::Available);
    CHECK(tie.witnessCount == 2U);
    CHECK(tie.distance.si() == doctest::Approx(0.5 / std::numbers::sqrt2));
    const CircleClearance small = measureCircleClearance(prepared, {{at(0.0, 0.5), m(0.1)}, m(1e-10)}, std::span<Witness2>(buffer).first(1));
    CHECK(small.witnessCount == 1U);
    CHECK(small.overflowed);
    const CircleClearance none = measureCircleClearance(prepared, {{at(0.0, 0.5), m(0.1)}, m(1e-10)}, std::span<Witness2>{});
    CHECK(none.status == ClearanceStatus::Available);
    CHECK(none.witnessCount == 0U);
    CHECK(measureCircleClearance(prepared, {{at(0.0, -0.5), m(0.1)}, m(1e-10)}, buffer).status == ClearanceStatus::CentreInside);
    CHECK(measureCircleClearance(prepared, {{at(0.0, 0.0), m(0.1)}, m(1e-10)}, buffer).status == ClearanceStatus::OnBoundary);
    CHECK(measureCircleClearance(prepared, {{at(0.95, 0.5), m(0.1)}, m(1e-10)}, buffer).status == ClearanceStatus::OutsideDomain);
    CHECK(measureCircleClearance(prepared, {{at(std::nan(""), 0.5), m(0.1)}, m(1e-10)}, buffer).status == ClearanceStatus::OutsideDomain);
    CHECK(measureCircleClearance(prepared, {{at(0.0, 0.5), m(0.1)}, m(-1.0)}, buffer).status == ClearanceStatus::OutsideDomain);
    CHECK(measureCircleClearance(prepared, {{at(0.0, 0.5), m(0.0)}, m(1e-10)}, buffer).status == ClearanceStatus::OutsideDomain);
}
