#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <doctest.h>
#include <limits>
#include <numbers>
#include <random>
#include <span>
#include <sub0math/units.hpp>
#include <utility>
#include <vector>

#include "sub0physics/geometry.hpp"
#include "test_profiles.hpp"

using namespace sub0physics;
using sub0math::metres;
using sub0math::Point2d;
using sub0math::Vec2d;
using sub0math::f64::Length;
using test_profiles::at;

namespace
{

using LD = long double;
using Node = sub0math::Point2<Length>;

// ---- The reference vehicle experiments' synthetic profiles (trail fixtures) ----

std::vector<std::pair<double, double>> crown(int segments = 32)
{
    // A 15 mm polynomial crown over 0.8 m: 16 h t^2 (1 - t)^2.
    std::vector<std::pair<double, double>> p{{-3.0, 0.0}};
    for (int i = 0; i <= segments; ++i)
    {
        const double t = static_cast<double>(i) / segments;
        p.emplace_back(0.6 + 0.8 * t, 0.015 * 16 * t * t * (1 - t) * (1 - t));
    }
    p.emplace_back(5.0, 0.0);
    return p;
}

std::vector<std::pair<double, double>> squareRock()
{
    return {{-3.0, 0.0}, {0.6, 0.0}, {0.6, 0.04}, {0.9, 0.04}, {0.9, 0.0}, {5.0, 0.0}};
}

std::vector<std::pair<double, double>> facetedRock()
{
    return {{-3.0, 0.0}, {0.6, 0.0}, {0.68, 0.045}, {0.85, 0.06}, {1.0, 0.015}, {1.08, 0.0}, {5.0, 0.0}};
}

std::vector<std::pair<double, double>> squareObstacle()
{
    return {{-3.0, 0.0}, {0.6, 0.0}, {0.6, 0.2}, {0.9, 0.2}, {0.9, 0.0}, {5.0, 0.0}};
}

PolylineProfile profileOf(const std::vector<std::pair<double, double>>& points)
{
    std::vector<ProfilePoint> v;
    for (const auto& [x, y] : points)
    {
        v.push_back(at(x, y));
    }
    auto p = PolylineProfile::create(std::move(v));
    REQUIRE(p);
    return *p;
}

// Splits every segment into `parts` pieces: the same physical terrain, authored more finely.
std::vector<std::pair<double, double>> subdivide(const std::vector<std::pair<double, double>>& points, int parts)
{
    std::vector<std::pair<double, double>> out{points.front()};
    for (std::size_t i = 1; i < points.size(); ++i)
    {
        const auto [ax, ay] = points[i - 1];
        const auto [bx, by] = points[i];
        for (int k = 1; k < parts; ++k)
        {
            // Exactly representable fractions keep the inserted points on the segment's line to rounding.
            const double t = static_cast<double>(k) / parts;
            out.emplace_back(ax + (bx - ax) * t, ay + (by - ay) * t);
        }
        out.push_back(points[i]);
    }
    return out;
}

// A casing ring of `n` nodes on centreline radius `r` about `hub`, with radial dents from a seed.
std::vector<Node> ring(Point2d hub, double r, std::size_t n, double dent, std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<Node> nodes(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        const double angle = 2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(n);
        const double radius = r + dent * u(rng);
        nodes[i] = {metres(hub.x + radius * std::cos(angle)), metres(hub.y + radius * std::sin(angle))};
    }
    return nodes;
}

// ---- Independent references in long double ----

// Signed distance of a segment from a convex polygon: the largest separation over every direction (exact for
// convex bodies, in contact or not), by a dense angle scan and golden refinement.
LD signedDistanceReference(Point2d a, Point2d b, const TerrainPiece& p)
{
    const auto separation = [&](LD angle) {
        const LD nx = std::cos(angle);
        const LD ny = std::sin(angle);
        LD hp = -std::numeric_limits<LD>::infinity();
        for (const auto& v : p.vertices)
        {
            hp = std::max(hp, v.x * nx + v.y * ny);
        }
        return std::min(a.x * nx + a.y * ny, b.x * nx + b.y * ny) - hp;
    };
    constexpr int Scan = 4000;
    const LD twoPi = 2 * std::numbers::pi_v<LD>;
    int best = 0;
    LD bestValue = separation(0);
    for (int i = 1; i < Scan; ++i)
    {
        const LD v = separation(twoPi * i / Scan);
        if (v > bestValue)
        {
            bestValue = v;
            best = i;
        }
    }
    LD lo = twoPi * (best - 1) / Scan;
    LD hi = twoPi * (best + 1) / Scan;
    const LD g = (std::sqrt(LD{5}) - 1) / 2;
    for (int it = 0; it < 120; ++it)
    {
        const LD m1 = hi - g * (hi - lo);
        const LD m2 = lo + g * (hi - lo);
        if (separation(m1) > separation(m2))
        {
            hi = m2;
        }
        else
        {
            lo = m1;
        }
    }
    return std::max(bestValue, separation((lo + hi) / 2));
}

// Distance from a segment to an authored (uncanonicalised) polyline: per authored segment, the point-segment distance
// along the axis is convex, so a golden search per pair is exact.
LD polylineDistanceReference(Point2d a, Point2d b, const std::vector<std::pair<double, double>>& points)
{
    LD best = std::numeric_limits<LD>::infinity();
    for (std::size_t s = 1; s < points.size(); ++s)
    {
        const LD px = points[s - 1].first;
        const LD py = points[s - 1].second;
        const LD dx = points[s].first - px;
        const LD dy = points[s].second - py;
        const auto f = [&](LD t) {
            const LD x = a.x + t * (static_cast<LD>(b.x) - a.x);
            const LD y = a.y + t * (static_cast<LD>(b.y) - a.y);
            const LD u = std::clamp<LD>(((x - px) * dx + (y - py) * dy) / (dx * dx + dy * dy), 0, 1);
            const LD ex = px + u * dx - x;
            const LD ey = py + u * dy - y;
            return ex * ex + ey * ey;
        };
        LD lo = 0;
        LD hi = 1;
        const LD g = (std::sqrt(LD{5}) - 1) / 2;
        for (int it = 0; it < 120; ++it)
        {
            const LD m1 = hi - g * (hi - lo);
            const LD m2 = lo + g * (hi - lo);
            if (f(m1) < f(m2))
            {
                hi = m2;
            }
            else
            {
                lo = m1;
            }
        }
        best = std::min({best, f(0), f(1), f((lo + hi) / 2)});
    }
    return std::sqrt(best);
}

Point2d si(const Node& p)
{
    return {p.x.si(), p.y.si()};
}

std::uint64_t mix(std::uint64_t h, double v)
{
    h ^= std::bit_cast<std::uint64_t>(v);
    h *= 0x100000001b3ULL;
    return h ^ (h >> 29);
}

} // namespace

TEST_CASE("Convex terrain: canonical pieces follow the prism construction")
{
    const ConvexTerrain square(profileOf(squareRock()));
    // Vertical faces have no piece: three pieces, floor one metre under the lowest vertex.
    REQUIRE(square.pieces().size() == 3U);
    CHECK(square.floor() == -1.0);
    const TerrainPiece& top = square.pieces()[1];
    CHECK(top.vertices[3] == Point2d{0.6, 0.04});
    CHECK(top.vertices[2] == Point2d{0.9, 0.04});
    CHECK(top.normals[2] == Vec2d{0.0, 1.0});
    CHECK(top.normals[3] == Vec2d{-1.0, 0.0}); // the rock's left face is the piece's side
    CHECK(square.pieceSegments()[1] == 2U);

    // Collinear subdivision is merged away; a reversal or a turn is not.
    const ConvexTerrain fine(profileOf(subdivide(squareRock(), 4)));
    CHECK(std::ranges::equal(fine.vertexX(), square.vertexX()));
    CHECK(std::ranges::equal(fine.vertexY(), square.vertexY()));
    REQUIRE(fine.pieces().size() == square.pieces().size());
    for (std::size_t i = 0; i < fine.pieces().size(); ++i)
    {
        CHECK(fine.pieces()[i] == square.pieces()[i]);
    }

    // The crown keeps its facets (distinct slopes); the flat run before it merges with the crown's first point.
    const ConvexTerrain bump(profileOf(crown()));
    CHECK(bump.vertexX().size() == crown().size()); // every facet turns, so nothing merges
    CHECK(bump.pieces().size() == crown().size() - 1U);
}

TEST_CASE("Convex terrain: contacts match the direction search per piece, and the polyline distance overall")
{
    const std::array fixtures{crown(), squareRock(), facetedRock(), squareObstacle()};
    std::mt19937_64 rng(2026);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::array<TerrainContact, 64> buffer{};
    std::size_t compared = 0;
    std::size_t vertexWitnesses = 0;
    for (const auto& authored : fixtures)
    {
        const ConvexTerrain terrain(profileOf(authored));
        for (int k = 0; k < 40; ++k)
        {
            // A 64-node ring, 8 mm skin, hub over the feature, pressed a few millimetres into the surface.
            const double hx = 0.4 + 0.7 * u(rng);
            double surface = -1.0;
            for (std::size_t s = 1; s < authored.size(); ++s)
            {
                if (authored[s - 1].first <= hx && hx <= authored[s].first)
                {
                    surface = std::max({surface, authored[s - 1].second, authored[s].second});
                }
            }
            const Point2d hub{hx, std::max(surface, 0.0) + 0.338 - 0.006 * u(rng) + (k % 4 == 0 ? 0.03 : 0.0)};
            const std::vector<Node> nodes = ring(hub, 0.33, 64, 0.002, rng);
            const TerrainContactSet set = terrain.ringContacts(nodes, metres(0.008), metres(0.002), buffer);
            REQUIRE(set.status == TerrainQueryStatus::Available);
            REQUIRE_FALSE(set.overflowed);
            for (std::size_t c = 0; c < set.count; ++c)
            {
                const TerrainContact& tc = buffer[c];
                const Point2d a = si(nodes[tc.edge]);
                const Point2d b = si(nodes[(tc.edge + 1) % nodes.size()]);
                const auto reference = static_cast<double>(signedDistanceReference(a, b, terrain.pieces()[tc.piece])) - 0.008;
                CHECK(std::abs(tc.gap.si() - reference) < 1e-9);
                CHECK(tc.gap.si() <= 0.002);
                ++compared;
                // Multiplicity: at most one contact per (edge, piece).
                if (c > 0)
                {
                    CHECK((buffer[c - 1].edge != tc.edge || buffer[c - 1].piece != tc.piece));
                }
                if (tc.witness.kind == WitnessKind::Vertex && !tc.axisCrosses)
                {
                    ++vertexWitnesses;
                    const auto cone = sub0math::vertexNormalCone(terrain.pieces()[tc.piece], tc.witness.index);
                    CHECK(cone.contains(tc.witness.normal, 1e-12));
                }
            }
            // The smallest gap is the distance to the authored polyline (the union of the pieces) less the skin.
            for (std::size_t e = 0; e < nodes.size(); ++e)
            {
                const TerrainContactSet one = terrain.contacts({{nodes[e], nodes[(e + 1) % nodes.size()]}, metres(0.008)},
                                                               metres(0.05), std::span<TerrainContact>{});
                // Exact when within the margin (every nearer piece is in the window); beyond it, both exceed the margin.
                const auto polyline = polylineDistanceReference(si(nodes[e]), si(nodes[(e + 1) % nodes.size()]), authored);
                if (one.minimumGap.si() <= -0.008)
                {
                    continue; // The axis is inside the solid: the polyline distance is unsigned there.
                }
                if (static_cast<double>(polyline) - 0.008 > 0.05)
                {
                    CHECK(one.minimumGap.si() > 0.05);
                }
                else
                {
                    CHECK(std::abs(one.minimumGap.si() - (static_cast<double>(polyline) - 0.008)) < 1e-9);
                }
            }
        }
    }
    CHECK(compared > 300U);
    CHECK(vertexWitnesses > 20U);
}

TEST_CASE("Convex terrain: a square edge gives the corner vertex, inside the terrain's corner cone, at the corner")
{
    const ConvexTerrain terrain(profileOf(squareRock()));
    bool convex = false;
    const auto cone = terrain.cornerCone(2, convex); // the rock's top left corner (0.6, 0.04)
    CHECK(convex);
    CHECK(cone.from == Vec2d{0.0, 1.0});
    CHECK(cone.to == Vec2d{-1.0, 0.0});
    bool concave = true;
    (void)terrain.cornerCone(1, concave); // the foot of the face
    CHECK_FALSE(concave);

    // Edges approaching the corner from every direction inside its cone: the witness is the corner, never a face.
    std::array<TerrainContact, 8> buffer{};
    for (int k = 1; k < 30; ++k)
    {
        const double angle = std::numbers::pi / 2.0 + (std::numbers::pi / 2.0) * k / 30.0;
        const Vec2d d{std::cos(angle), std::sin(angle)};
        const Point2d centre{0.6 + 0.009 * d.x, 0.04 + 0.009 * d.y};
        const Vec2d along{-d.y, d.x};
        const Node a{metres(centre.x - 0.01 * along.x), metres(centre.y - 0.01 * along.y)};
        const Node b{metres(centre.x + 0.01 * along.x), metres(centre.y + 0.01 * along.y)};
        const TerrainContactSet set = terrain.contacts({{a, b}, metres(0.008)}, metres(0.002), buffer);
        REQUIRE(set.count >= 1U);
        for (std::size_t c = 0; c < set.count; ++c)
        {
            const TerrainContact& tc = buffer[c];
            CHECK(tc.witness.kind == WitnessKind::Vertex);
            CHECK(tc.witness.point.x.si() == 0.6);
            CHECK(tc.witness.point.y.si() == 0.04);
            CHECK(cone.contains(tc.witness.normal, 1e-12));
            // The true approach direction, well inside the 0.005 rad budget of the reference.
            CHECK(std::abs(tc.witness.normal.x * d.y - tc.witness.normal.y * d.x) < 1e-12);
            CHECK(tc.gap.si() == doctest::Approx(0.001).epsilon(1e-9));
        }
    }
}

TEST_CASE("Convex terrain: collinear subdivision gives identical contacts, bit for bit")
{
    const std::array fixtures{squareRock(), facetedRock(), squareObstacle()};
    std::mt19937_64 rng(5);
    std::array<TerrainContact, 64> coarseOut{};
    std::array<TerrainContact, 64> fineOut{};
    for (const auto& authored : fixtures)
    {
        const ConvexTerrain coarse(profileOf(authored));
        const ConvexTerrain fine(profileOf(subdivide(authored, 8)));
        for (int k = 0; k < 20; ++k)
        {
            const std::vector<Node> nodes = ring({0.5 + 0.03 * k, 0.335 + (k % 3 == 0 ? 0.04 : 0.0)}, 0.33, 64, 0.001, rng);
            const auto a = coarse.ringContacts(nodes, metres(0.008), metres(0.002), coarseOut);
            const auto b = fine.ringContacts(nodes, metres(0.008), metres(0.002), fineOut);
            REQUIRE(a.count == b.count);
            for (std::size_t c = 0; c < a.count; ++c)
            {
                CHECK(coarseOut[c].gap == fineOut[c].gap);
                CHECK(coarseOut[c].witness.point == fineOut[c].witness.point);
                CHECK(coarseOut[c].witness.normal == fineOut[c].witness.normal);
                CHECK(coarseOut[c].witness.kind == fineOut[c].witness.kind);
            }
        }
    }
}

TEST_CASE("Convex terrain: buffers, overflow and invalid queries")
{
    const ConvexTerrain terrain(profileOf(squareRock()));
    const CasingEdge flat{{at(0.0, 0.005), at(1.5, 0.005)}, metres(0.008)}; // touches every piece
    std::array<TerrainContact, 2> two{};
    const auto set = terrain.contacts(flat, metres(0.0), two, 7);
    CHECK(set.found == 3U);
    CHECK(set.count == 2U);
    CHECK(set.overflowed);
    CHECK(two[0].edge == 7U);
    CHECK(two[0].piece == 0U);
    CHECK(two[1].piece == 1U);
    CHECK(set.minimumGap.si() < 0.0);
    CHECK(terrain.contacts({{at(0.0, std::nan("")), at(1.0, 0.0)}, metres(0.008)}, metres(0.0), two).status ==
          TerrainQueryStatus::InvalidQuery);
    CHECK(terrain.contacts(flat, metres(-1.0), two).status == TerrainQueryStatus::InvalidQuery);
    // Far above: nothing is measured at all.
    const auto clear = terrain.contacts({{at(0.0, 2.0), at(0.1, 2.0)}, metres(0.008)}, metres(0.01), two);
    CHECK(clear.pairsTested == 0U);
    CHECK(clear.count == 0U);
}

TEST_CASE("Clearance bound: conservative under translation, rotation and deformation; the rotation term is needed")
{
    // Property test: random steps of a ring over every fixture. The true gap at 1001 instants never falls below the
    // certified bound. A bound without the rotation term is violated by some of the same steps.
    const std::array fixtures{crown(), squareRock(), facetedRock(), squareObstacle()};
    std::mt19937_64 rng(77);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::size_t certified = 0;
    std::size_t violatedWithoutRotation = 0;
    for (const auto& authored : fixtures)
    {
        const ConvexTerrain terrain(profileOf(authored));
        for (int k = 0; k < 30; ++k)
        {
            constexpr std::size_t N = 32;
            const double r = 0.33;
            const Point2d c0{0.3 + 0.8 * u(rng), 0.36 + 0.03 * u(rng)};
            // A curved hub path (quadratic), a spin of up to 0.3 rad and node displacements in the hub frame.
            const Vec2d v{0.03 * (u(rng) - 0.5), -0.03 * u(rng)};
            const Vec2d acc{0.01 * (u(rng) - 0.5), 0.01 * (u(rng) - 0.5)};
            const double spin = 0.3 * (u(rng) - 0.5);
            const double theta0 = u(rng);
            std::array<Vec2d, N> q0{};
            std::array<Vec2d, N> dq{};
            for (std::size_t i = 0; i < N; ++i)
            {
                const double angle = 2.0 * std::numbers::pi * static_cast<double>(i) / N;
                q0[i] = {r * std::cos(angle), r * std::sin(angle)};
                dq[i] = {0.002 * (u(rng) - 0.5), 0.002 * (u(rng) - 0.5)};
            }
            const auto nodesAt = [&](double t) {
                std::vector<Node> nodes(N);
                const Point2d c{c0.x + v.x * t + acc.x * t * t, c0.y + v.y * t + acc.y * t * t};
                const double th = theta0 + spin * t;
                for (std::size_t i = 0; i < N; ++i)
                {
                    const Vec2d q = q0[i] + dq[i] * t;
                    nodes[i] = {metres(c.x + std::cos(th) * q.x - std::sin(th) * q.y),
                                metres(c.y + std::sin(th) * q.x + std::cos(th) * q.y)};
                }
                return nodes;
            };
            double travel = 0.0;
            double nodeTravel = 0.0;
            for (int s = 0; s <= 1000; ++s)
            {
                const double t = s / 1000.0;
                travel = std::max(travel, std::hypot(v.x * t + acc.x * t * t, v.y * t + acc.y * t * t));
            }
            for (std::size_t i = 0; i < N; ++i)
            {
                nodeTravel = std::max(nodeTravel, std::sqrt(lengthSquared(dq[i])));
            }
            const std::vector<Node> start = nodesAt(0.0);
            // The extremal hub excursion is sampled; pad by a micrometre for the samples between.
            const StepMotionBound motion{{metres(c0.x), metres(c0.y)}, metres(travel + 1e-6), std::abs(spin), metres(nodeTravel)};
            const StepMotionBound noRotation{motion.hubCentre, motion.hubTravel, 0.0, motion.nodeTravel};
            const ClearanceBound bound = terrain.certifyRingClearance(start, metres(0.008), motion, metres(0.5));
            const ClearanceBound naive = terrain.certifyRingClearance(start, metres(0.008), noRotation, metres(0.5));
            REQUIRE(bound.status == TerrainQueryStatus::Available);
            double trueMinimum = std::numeric_limits<double>::infinity();
            for (int s = 0; s <= 1000; ++s)
            {
                const auto set = terrain.ringContacts(nodesAt(s / 1000.0), metres(0.008), metres(0.0), std::span<TerrainContact>{});
                trueMinimum = std::min(trueMinimum, set.minimumGap.si());
            }
            // The bound is capped at the horizon; below it, it must hold.
            CHECK(bound.clearance.si() <= trueMinimum + 1e-12);
            if (bound.certified)
            {
                ++certified;
                CHECK(trueMinimum > 0.0);
            }
            if (naive.clearance.si() > trueMinimum + 1e-12)
            {
                ++violatedWithoutRotation;
            }
        }
    }
    CHECK(certified > 5U);
    CHECK(violatedWithoutRotation > 0U);
}

TEST_CASE("Convex terrain: results are deterministic across compilers (digest of a fixed batch)")
{
    const ConvexTerrain terrain(profileOf(facetedRock()));
    std::mt19937_64 rng(99);
    std::array<TerrainContact, 64> buffer{};
    std::uint64_t h = 0xcbf29ce484222325ULL;
    for (int k = 0; k < 50; ++k)
    {
        // The generator's integers are portable; its real distributions are not, so build doubles from bits.
        const auto unit = [&] { return static_cast<double>(rng() >> 11) * 0x1.0p-53; };
        std::vector<Node> nodes(64);
        for (std::size_t i = 0; i < nodes.size(); ++i)
        {
            nodes[i] = {metres(0.4 + 0.8 * unit()), metres(0.07 * unit() - 0.005)};
        }
        const auto set = terrain.ringContacts(nodes, metres(0.008), metres(0.002), buffer);
        for (std::size_t c = 0; c < set.count; ++c)
        {
            h = mix(h, buffer[c].gap.si());
            h = mix(h, buffer[c].witness.point.x.si());
            h = mix(h, buffer[c].witness.point.y.si());
            h = mix(h, buffer[c].witness.normal.x);
            h = mix(h, buffer[c].witness.normal.y);
            h = mix(h, static_cast<double>(buffer[c].witness.index));
        }
    }
    MESSAGE("digest " << h);
    CHECK(h == 8356751711351349667ULL); // GCC 16 and Clang 22, -march=x86-64-v3 (FMA available, contraction off)
}
