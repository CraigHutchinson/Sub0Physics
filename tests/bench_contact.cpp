// Casing contact geometry (M2b): the library's ring query paired with a hand-written twin that does the same job on
// raw arrays: the canonical vertices, a binary-searched x window, the height cull and the closed-form capsule-slab
// distance with its witness, the slab corners built in registers from the vertex arrays (no stored polygons). The
// twin is the reference; an API gap to it is recorded in docs/PERFORMANCE.md. Outputs are checked equal first.
#include <algorithm>
#include <array>
#include <bench_runner.hpp>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <numbers>
#include <sub0math/units.hpp>
#include <vector>

#include "sub0physics/geometry.hpp"

using namespace sub0physics;

namespace
{

struct TwinContact
{
    std::size_t edge, piece;
    double gap, px, py, nx, ny;
    std::size_t index;
    bool vertex;
};

struct Twin
{
    std::vector<double> x, y, right, top, normalX, normalY, left;
    std::vector<std::size_t> segment;
    double floor = 0.0;
    double highest = -1e300;

    explicit Twin(const ConvexTerrain& terrain)
    {
        x.assign(terrain.vertexX().begin(), terrain.vertexX().end());
        y.assign(terrain.vertexY().begin(), terrain.vertexY().end());
        floor = terrain.floor();
        for (std::size_t s = 0; s + 1 < x.size(); ++s)
        {
            if (x[s + 1] == x[s])
            {
                continue;
            }
            const sub0math::Vec2d n = terrain.segmentNormal(s);
            left.push_back(x[s]);
            right.push_back(x[s + 1]);
            top.push_back(std::max(y[s], y[s + 1]));
            highest = std::max(highest, top.back());
            normalX.push_back(n.x);
            normalY.push_back(n.y);
            segment.push_back(s);
        }
    }

    // One capsule against one slab, the closed form of `sub0math::contact` written out on scalars.
    TwinContact slab(std::size_t p, double ax, double ay, double bx, double by, double radius) const noexcept
    {
        const std::size_t s = segment[p];
        const double vx[4] = {x[s], x[s + 1], x[s + 1], x[s]};
        const double vy[4] = {floor, floor, y[s + 1], y[s]};
        const double nx[4] = {0.0, 1.0, normalX[p], -1.0};
        const double ny[4] = {-1.0, 0.0, normalY[p], 0.0};
        const double abx = bx - ax;
        const double aby = by - ay;
        const double ab2 = abx * abx + aby * aby;
        double oa[4], ob[4];
        double edgeGap = -1e300;
        std::size_t edgeAxis = 0;
        for (std::size_t i = 0; i < 4; ++i)
        {
            oa[i] = nx[i] * (ax - vx[i]) + ny[i] * (ay - vy[i]);
            ob[i] = nx[i] * (bx - vx[i]) + ny[i] * (by - vy[i]);
            const double g = std::min(oa[i], ob[i]);
            if (g > edgeGap)
            {
                edgeGap = g;
                edgeAxis = i;
            }
        }
        double lo = 1e300, hi = -1e300;
        std::size_t loV = 0, hiV = 0;
        for (std::size_t k = 0; k < 4; ++k)
        {
            const double o = -aby * (vx[k] - ax) + abx * (vy[k] - ay);
            if (o < lo)
            {
                lo = o;
                loV = k;
            }
            if (o > hi)
            {
                hi = o;
                hiV = k;
            }
        }
        const bool misses = ab2 > 0.0 && (lo > 0.0 || hi < 0.0);
        if (edgeGap > 0.0 || misses)
        {
            double best = 1e300, sx = 0.0, sy = 0.0, wx = 0.0, wy = 0.0;
            std::size_t index = 0;
            bool edge = false;
            const auto end = [&](double ex, double ey, const double* o) {
                for (std::size_t i = 0; i < 4; ++i)
                {
                    if (!(o[i] > 0.0))
                    {
                        continue;
                    }
                    const std::size_t j = (i + 1) & 3U;
                    const double dx = vx[j] - vx[i];
                    const double dy = vy[j] - vy[i];
                    const double t = std::clamp(((ex - vx[i]) * dx + (ey - vy[i]) * dy) / (dx * dx + dy * dy), 0.0, 1.0);
                    const double qx = t == 0.0 ? vx[i] : (t == 1.0 ? vx[j] : vx[i] + dx * t);
                    const double qy = t == 0.0 ? vy[i] : (t == 1.0 ? vy[j] : vy[i] + dy * t);
                    const double d2 = (ex - qx) * (ex - qx) + (ey - qy) * (ey - qy);
                    if (d2 < best)
                    {
                        best = d2;
                        index = t == 1.0 ? j : i;
                        edge = t != 0.0 && t != 1.0;
                        sx = ex;
                        sy = ey;
                        wx = qx;
                        wy = qy;
                    }
                }
            };
            end(ax, ay, oa);
            end(bx, by, ob);
            if (misses)
            {
                const std::size_t v = lo > 0.0 ? loV : hiV;
                const double t = ((vx[v] - ax) * abx + (vy[v] - ay) * aby) / ab2;
                if (t > 0.0 && t < 1.0)
                {
                    const double ox = ax + abx * t;
                    const double oy = ay + aby * t;
                    const double d2 = (ox - vx[v]) * (ox - vx[v]) + (oy - vy[v]) * (oy - vy[v]);
                    if (d2 < best)
                    {
                        best = d2;
                        index = v;
                        edge = false;
                        sx = ox;
                        sy = oy;
                        wx = vx[v];
                        wy = vy[v];
                    }
                }
            }
            const double d = std::sqrt(best);
            return {0, p, d - radius, wx, wy, edge ? nx[index] : (sx - wx) / d, edge ? ny[index] : (sy - wy) / d, index, !edge};
        }
        const double inv = 1.0 / std::sqrt(ab2);
        const double seg = std::max(-hi * inv, lo * inv);
        if (seg > edgeGap)
        {
            const bool plus = -hi * inv >= lo * inv;
            const std::size_t v = plus ? hiV : loV;
            const double sign = plus ? 1.0 : -1.0;
            return {0, p, seg - radius, vx[v], vy[v], -aby * inv * sign, abx * inv * sign, v, true};
        }
        const bool endA = oa[edgeAxis] <= ob[edgeAxis];
        return {0, p, edgeGap - radius, (endA ? ax : bx) - nx[edgeAxis] * edgeGap, (endA ? ay : by) - ny[edgeAxis] * edgeGap,
                nx[edgeAxis], ny[edgeAxis], edgeAxis, false};
    }

    std::size_t ring(const std::vector<double>& nx, const std::vector<double>& ny, double skin, double margin,
                     std::array<TwinContact, 64>& out) const noexcept
    {
        std::size_t count = 0;
        const std::size_t n = nx.size();
        const double pad = skin + margin;
        for (std::size_t e = 0; e < n; ++e)
        {
            const std::size_t f = e + 1 == n ? 0 : e + 1;
            const double ax = nx[e], ay = ny[e], bx = nx[f], by = ny[f];
            const double low = std::min(ay, by) - pad;
            if (low > highest)
            {
                continue; // above every slab
            }
            const double xlo = std::min(ax, bx) - pad;
            const double xhi = std::max(ax, bx) + pad;
            std::size_t p = static_cast<std::size_t>(std::lower_bound(right.begin(), right.end(), xlo) - right.begin());
            for (; p < left.size() && left[p] <= xhi; ++p)
            {
                if (low > top[p])
                {
                    continue;
                }
                TwinContact c = slab(p, ax, ay, bx, by, skin);
                if (c.gap <= margin && count < out.size())
                {
                    c.edge = e;
                    out[count++] = c;
                }
            }
        }
        return count;
    }
};

} // namespace

void benchContact(sub0bench::BenchRunner& runner)
{
    // The 40 mm square rock of the reference trail fixtures, a 64-node ring (0.33 m centreline, 8 mm skin) pressed
    // about 3 mm onto the rock's left corner and top.
    std::vector<ProfilePoint> points;
    for (const auto& [px, py] : std::array<std::array<double, 2>, 6>{{{-3.0, 0.0}, {0.6, 0.0}, {0.6, 0.04}, {0.9, 0.04}, {0.9, 0.0}, {5.0, 0.0}}})
    {
        points.push_back({sub0math::metres(px), sub0math::metres(py)});
    }
    const auto profile = PolylineProfile::create(std::move(points));
    const ConvexTerrain terrain(*profile);
    const Twin twin(terrain);
    constexpr std::size_t N = 64;
    std::vector<sub0math::Point2<sub0math::f64::Length>> nodes(N);
    std::vector<double> nx(N), ny(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        const double angle = 2.0 * std::numbers::pi * static_cast<double>(i) / N;
        nx[i] = 0.62 + 0.33 * std::cos(angle);
        ny[i] = 0.375 + 0.33 * std::sin(angle);
        nodes[i] = {sub0math::metres(nx[i]), sub0math::metres(ny[i])};
    }
    std::array<TerrainContact, 64> libraryOut{};
    std::array<TwinContact, 64> twinOut{};
    const auto set = terrain.ringContacts(nodes, sub0math::metres(0.008), sub0math::metres(0.001), libraryOut);
    const std::size_t count = twin.ring(nx, ny, 0.008, 0.001, twinOut);
    if (set.count != count || count == 0)
    {
        std::fprintf(stderr, "contact bench: twin and library disagree on the count (%zu vs %zu)\n", count, set.count);
        std::abort();
    }
    for (std::size_t c = 0; c < count; ++c)
    {
        const TerrainContact& l = libraryOut[c];
        const TwinContact& t = twinOut[c];
        if (l.edge != t.edge || l.piece != t.piece || l.gap.si() != t.gap || l.witness.point.x.si() != t.px ||
            l.witness.point.y.si() != t.py || l.witness.normal.x != t.nx || l.witness.normal.y != t.ny ||
            l.witness.index != t.index || (l.witness.kind == WitnessKind::Vertex) != t.vertex)
        {
            std::fprintf(stderr, "contact bench: twin and library disagree at contact %zu\n", c);
            std::abort();
        }
    }
    runner.run("ring contacts 64 edges, square rock, hand-written twin", [&] {
        ankerl::nanobench::doNotOptimizeAway(twin.ring(nx, ny, 0.008, 0.001, twinOut));
    });
    runner.run("ring contacts 64 edges, square rock, ConvexTerrain::ringContacts", [&] {
        ankerl::nanobench::doNotOptimizeAway(
            terrain.ringContacts(nodes, sub0math::metres(0.008), sub0math::metres(0.001), libraryOut));
    });
    const StepMotionBound motion{{sub0math::metres(0.62), sub0math::metres(0.375)}, sub0math::metres(0.005), 0.02,
                                 sub0math::metres(0.001)};
    runner.run("ring clearance bound 64 edges, ConvexTerrain::certifyRingClearance", [&] {
        ankerl::nanobench::doNotOptimizeAway(
            terrain.certifyRingClearance(nodes, sub0math::metres(0.008), motion, sub0math::metres(0.05)));
    });
}
