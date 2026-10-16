// Terrain representation decision experiment D0: see experiments/terrain_d0/README.md.

#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "terrain_d0.hpp"

using namespace d0;

namespace
{

volatile float sink = 0.0f;

template<typename F>
double bestSeconds(int repeats, F&& f)
{
    double best = std::numeric_limits<double>::infinity();
    for (int r = 0; r < repeats; ++r)
    {
        const auto t0 = std::chrono::steady_clock::now();
        f();
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double>(t1 - t0).count());
    }
    return best;
}

// ── Scenes ──────────────────────────────────────────────────────────────────────────────────────

Scene makeBerm()
{
    const double wallY = 0.3, height = 0.6, tan70 = std::tan(70.0 * Pi / 180.0), lipDepth = 0.03, lipThick = 0.02;
    const auto wallAt = [=](double z) { return wallY + z / tan70; };
    const double lipEdge = wallAt(height) - lipDepth;
    const std::vector<V2d> poly = {{-2, -1}, {3, -1}, {3, height}, {lipEdge, height}, {lipEdge, height - lipThick},
                                   {wallAt(height - lipThick), height - lipThick}, {wallY, 0}, {-2, 0}};
    Scene s;
    s.name = "Berm 70° with 3 cm lip";
    s.sdf = [poly](V3d p) { return sdPolygon(poly, {p.y, p.z}); };
    s.ground = [=](double, double y) { return y < wallY ? 0.0 : std::min(height, (y - wallY) * tan70); };
    s.meshes.push_back(makePrismX({{lipEdge, height - lipThick}, {wallAt(height - lipThick), height - lipThick},
                                   {wallAt(height), height}, {lipEdge, height}},
                                  -1.0, 6.0));
    s.domain = {{0, -1, -0.3}, {5, 1, 0.8}};
    s.queryRegion = {{1, 0.1, -0.03}, {4, 0.6, 0.65}};
    return s;
}

Scene makeRuts()
{
    const double width = 0.08, depth = 0.06;
    const std::vector<double> centres = {0.12, -0.12}; // descending, for polygon order
    const auto top = [=](double y) {
        for (double c : centres)
            if (std::abs(y - c) < width / 2) return -depth * 0.5 * (1.0 + std::cos(2.0 * Pi * (y - c) / width));
        return 0.0;
    };
    std::vector<V2d> poly = {{-2, -1}, {2, -1}, {2, 0}};
    for (double c : centres)
        for (int k = 0; k <= 320; ++k)
        {
            const double y = c + width / 2 - k * (width / 320.0);
            poly.push_back({y, top(y)});
        }
    poly.push_back({-2, 0});
    Scene s;
    s.name = "Two ruts 8 cm wide, 6 cm deep";
    s.sdf = [poly](V3d p) { return sdPolygon(poly, {p.y, p.z}); };
    s.ground = [=](double, double y) { return top(y); };
    s.domain = {{0, -1, -0.3}, {5, 1, 0.3}};
    s.queryRegion = {{1, -0.3, -0.07}, {4, 0.3, 0.04}};
    return s;
}

Scene makeRocks()
{
    Scene s;
    s.name = "Rock garden: 2 boxes + mesh rock, half buried";
    s.boxes.push_back(makeObb({2.0, -0.1, 0.04}, {0.15, 0.10, 0.12}, 0.35, 0.0));
    s.boxes.push_back(makeObb({2.8, 0.15, 0.0}, {0.12, 0.12, 0.10}, -0.6, 0.3));
    s.meshes.push_back(makeEllipsoidRock({3.4, -0.05, 0.0}, {0.18, 0.14, 0.11}));
    auto boxes = s.boxes;
    auto meshes = s.meshes;
    s.sdf = [boxes, meshes](V3d p) {
        double d = p.z;
        for (auto& b : boxes) d = std::min(d, sdObb(b, p).d);
        for (auto& m : meshes) d = std::min(d, sdConvexMesh(m, p).d);
        return d;
    };
    s.ground = [](double, double) { return 0.0; };
    s.domain = {{0, -1, -0.3}, {5, 1, 0.4}};
    s.queryRegion = {{1.6, -0.45, -0.02}, {3.7, 0.45, 0.25}};
    return s;
}

// ── Sampling and statistics ─────────────────────────────────────────────────────────────────────

struct Sample
{
    V3d p;
    double ref;
    V3d n;
    bool steep;
};

std::vector<Sample> samplePoints(const std::function<double(V3d)>& sdf, const Box3& r, std::size_t count,
                                 std::uint32_t seed, const std::function<bool(V3d)>& filter = {})
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> ux(r.lo.x, r.hi.x), uy(r.lo.y, r.hi.y), uz(r.lo.z, r.hi.z);
    std::vector<Sample> out;
    for (std::size_t tries = 0; out.size() < count && tries < count * 2000; ++tries)
    {
        const V3d p{ux(rng), uy(rng), uz(rng)};
        if (filter && !filter(p)) continue;
        const double d = sdf(p);
        if (d < 0.002 || d > 0.03) continue;
        const V3d n = referenceNormal(sdf, p);
        out.push_back({p, d, n, n.z < std::cos(20.0 * Pi / 180.0)});
    }
    return out;
}

struct Stats
{
    double rms = 0, p99 = 0, max = 0;
};

Stats stats(std::vector<double> v)
{
    Stats s;
    if (v.empty()) return s;
    double sum = 0;
    for (double x : v) sum += x * x;
    std::sort(v.begin(), v.end());
    s.rms = std::sqrt(sum / static_cast<double>(v.size()));
    s.p99 = v[std::min(v.size() - 1, static_cast<std::size_t>(0.99 * static_cast<double>(v.size())))];
    s.max = v.back();
    return s;
}

struct QHit
{
    float d;
    V3f n;
    bool known;
};

struct Accuracy
{
    Stats gapMm, normalDeg, steepGapMm;
    std::size_t unknown = 0, samples = 0, steep = 0;
    double queryNs = 0;
};

template<typename Q>
Accuracy measure(const std::vector<Sample>& samples, Q&& q)
{
    Accuracy a;
    std::vector<double> g, n, sg;
    for (const auto& s : samples)
    {
        const QHit h = q(cast<float>(s.p));
        if (!h.known)
        {
            ++a.unknown;
            continue;
        }
        const double e = std::abs(static_cast<double>(h.d) - s.ref) * 1000.0;
        g.push_back(e);
        if (s.steep) sg.push_back(e);
        if (s.ref <= 0.005) // near the surface; farther out the medial axis makes the nearest normal ambiguous
            n.push_back(std::acos(std::clamp(dot(cast<double>(h.n), s.n), -1.0, 1.0)) * 180.0 / Pi);
    }
    a.gapMm = stats(g);
    a.normalDeg = stats(n);
    a.steepGapMm = stats(sg);
    a.samples = samples.size();
    a.steep = sg.size();
    std::vector<V3f> pts;
    for (const auto& s : samples) pts.push_back(cast<float>(s.p));
    const double sec = bestSeconds(5, [&] {
        float acc = 0;
        for (const auto& p : pts) acc += q(p).d;
        sink = acc;
    });
    a.queryNs = sec / static_cast<double>(pts.size()) * 1e9;
    return a;
}

std::string fmt(double v, int prec = 2)
{
    std::ostringstream o;
    o.setf(std::ios::fixed);
    o.precision(prec);
    o << v;
    return o.str();
}

std::string bytesText(double b)
{
    if (b >= 1024.0 * 1024.0 * 1024.0) return fmt(b / (1024.0 * 1024.0 * 1024.0), 1) + " GiB";
    if (b >= 1024.0 * 1024.0) return fmt(b / (1024.0 * 1024.0), 1) + " MiB";
    if (b >= 1024.0) return fmt(b / 1024.0, 1) + " KiB";
    return fmt(b, 0) + " B";
}

std::size_t featureBytes(const Scene& s)
{
    std::size_t b = s.boxes.size() * sizeof(ObbT<float>);
    for (auto& m : s.meshes) b += m.v.size() * 12 + m.tri.size() * (12 + 16);
    return b;
}

double area(const Box3& b) { return (b.hi.x - b.lo.x) * (b.hi.y - b.lo.y); }

std::unique_ptr<HeightField> bakeField(const Box3& domain, float spacing, const std::function<double(double, double)>& top)
{
    auto hf = std::make_unique<HeightField>(domain.lo.x, domain.lo.y, domain.hi.x, domain.hi.y, spacing);
    hf->bake(top);
    return hf;
}

float bandFor(float voxel) { return std::max(4.0f * voxel, 0.05f); }

} // namespace

int main(int argc, char** argv)
{
    const std::string outPath = argc > 1 ? argv[1] : "terrain_d0_results.md";
    std::ostringstream md;
    md << "# Terrain D0 results\n\n";
    md << "Generated by `experiments/terrain_d0`. Errors are against exact analytic references, at points 2–30 mm "
          "outside the solid. Steep = reference normal more than 20° from vertical. Normal errors use only points 2–5 mm from "
          "the surface, away from medial axes where the nearest surface switches. Times are best of 5 on this "
          "host.\n\n";

    const std::vector<Scene> scenes = {makeBerm(), makeRuts(), makeRocks()};
    const std::vector<float> spacings = {0.01f, 0.02f};

    // ── 1. Accuracy, memory and query cost ──────────────────────────────────────────────────────
    md << "## 1. Accuracy, memory and query cost\n\n";
    md << "| Fixture | Representation | Spacing | Bytes/m² | Gap RMS mm | Gap p99 mm | Gap max mm | Steep gap p99 mm "
          "| Normal p99 ° (2–5 mm) | Normal max ° (2–5 mm) | Query ns | Unknown |\n";
    md << "|---|---|---|---|---|---|---|---|---|---|---|---|\n";
    for (const auto& scene : scenes)
    {
        const auto samples = samplePoints(scene.sdf, scene.queryRegion, 20000, 1234);
        std::cerr << scene.name << ": " << samples.size() << " samples\n";
        for (float h : spacings)
        {
            const double a = area(scene.domain);
            const bool hasFeatures = !scene.boxes.empty() || !scene.meshes.empty();
            const auto row = [&](const std::string& rep, double bytes, const Accuracy& acc) {
                md << "| " << scene.name << " | " << rep << " | " << fmt(h * 100.0, 0) << " cm | " << bytesText(bytes / a)
                   << " | " << fmt(acc.gapMm.rms) << " | " << fmt(acc.gapMm.p99) << " | " << fmt(acc.gapMm.max) << " | "
                   << fmt(acc.steepGapMm.p99) << " | " << fmt(acc.normalDeg.p99, 1) << " | " << fmt(acc.normalDeg.max, 1)
                   << " | " << fmt(acc.queryNs, 0) << " | " << acc.unknown << " |\n";
            };

            if (hasFeatures)
            {
                // B only: everything, features included, baked into heights along vertical rays.
                auto all = bakeField(scene.domain, h, [&](double x, double y) {
                    return topByRay(scene.sdf, x, y, scene.domain.hi.z, scene.domain.lo.z);
                });
                row("B (all baked into heights)", static_cast<double>(all->baseBytes()),
                    measure(samples, [&](V3f p) { auto r = all->query(p); return QHit{r.d, r.n, true}; }));
            }
            auto ground = bakeField(scene.domain, h, scene.ground);
            HeightFieldWithFeatures be{ground.get(), {}, {}};
            for (auto& b : scene.boxes) be.boxes.push_back(castObb<float>(b));
            for (auto& m : scene.meshes) be.meshes.push_back(castMesh<float>(m));
            row(hasFeatures ? "B+E (heights + exact features)" : "B (heights)",
                static_cast<double>(ground->baseBytes() + featureBytes(scene)),
                measure(samples, [&](V3f p) { auto r = be.query(p); return QHit{r.d, r.n, true}; }));

            SdfBricks bricks(h, bandFor(h));
            bricks.bake(scene.sdf, scene.domain);
            row("H (SDF bricks, band " + fmt(bricks.band() * 100.0, 0) + " cm)", static_cast<double>(bricks.bytes()),
                measure(samples, [&](V3f p) { auto r = bricks.query(p); return QHit{r.d, r.n, r.known}; }));
        }
    }

    // ── 2. Memory per km of trail ───────────────────────────────────────────────────────────────
    md << "\n## 2. Memory per km of trail (4 m wide, ruts fixture as typical ground)\n\n";
    md << "| Spacing | B heights | H SDF bricks |\n|---|---|---|\n";
    for (float h : spacings)
    {
        const Scene& ruts = scenes[1];
        auto hf = bakeField(ruts.domain, h, ruts.ground);
        SdfBricks bricks(h, bandFor(h));
        bricks.bake(ruts.sdf, ruts.domain);
        const double a = area(ruts.domain);
        md << "| " << fmt(h * 100.0, 0) << " cm | " << bytesText(static_cast<double>(hf->baseBytes()) / a * 4000.0) << " | "
           << bytesText(static_cast<double>(bricks.bytes()) / a * 4000.0) << " |\n";
    }

    // ── 3. R0 corridor extraction ───────────────────────────────────────────────────────────────
    md << "\n## 3. R0 corridor extraction (ruts fixture; x 1.52–3.63 m)\n\n";
    md << "| Spacing | Corridor width | Columns | B: heights + triangle normals µs | H: surface by sphere tracing + normals µs |\n";
    md << "|---|---|---|---|---|\n";
    for (float h : spacings)
    {
        const Scene& ruts = scenes[1];
        auto hf = bakeField(ruts.domain, h, ruts.ground);
        SdfBricks bricks(h, bandFor(h));
        bricks.bake(ruts.sdf, ruts.domain);
        for (double width : {0.1, 0.3})
        {
            const double x0 = 1.52, x1 = 3.63, y0 = -width / 2, y1 = width / 2;
            const int i0 = static_cast<int>(std::lround(x0 / h)), i1 = static_cast<int>(std::lround(x1 / h));
            const int j0 = static_cast<int>(std::lround((y0 + 1.0) / h)), j1 = static_cast<int>(std::lround((y1 + 1.0) / h));
            const std::size_t cols = static_cast<std::size_t>((i1 - i0 + 1) * (j1 - j0 + 1));
            std::vector<float> topZ(cols), nx(cols * 2), ny(cols * 2), nz(cols * 2);
            const double tb = bestSeconds(5, [&] {
                std::size_t k = 0;
                for (int j = j0; j <= j1; ++j)
                    for (int i = i0; i <= i1; ++i, ++k)
                    {
                        const float za = hf->z(i, j), zb = hf->z(i + 1, j), zc = hf->z(i, j + 1), zd = hf->z(i + 1, j + 1);
                        topZ[k] = za;
                        const V3f n0 = normalized(V3f{-(zb - za) / h, -(zd - zb) / h, 1.0f});
                        const V3f n1 = normalized(V3f{-(zd - zc) / h, -(zc - za) / h, 1.0f});
                        nx[2 * k] = n0.x; ny[2 * k] = n0.y; nz[2 * k] = n0.z;
                        nx[2 * k + 1] = n1.x; ny[2 * k + 1] = n1.y; nz[2 * k + 1] = n1.z;
                    }
                sink = topZ[cols / 2] + nz[cols];
            });
            const double th = bestSeconds(5, [&] {
                std::size_t k = 0;
                for (int j = j0; j <= j1; ++j)
                {
                    float zStart = static_cast<float>(ruts.domain.hi.z);
                    for (int i = i0; i <= i1; ++i, ++k)
                    {
                        const float x = static_cast<float>(i) * h, y = static_cast<float>(j) * h - 1.0f;
                        float z = zStart;
                        SdfBricks::Hit q{};
                        for (int it = 0; it < 128; ++it)
                        {
                            q = bricks.query({x, y, z});
                            const float d = q.known ? q.d : bricks.band();
                            if (d < 1e-4f) break;
                            z -= std::max(d, 1e-4f);
                        }
                        topZ[k] = z;
                        nx[k] = q.n.x; ny[k] = q.n.y; nz[k] = q.n.z;
                        zStart = z + bricks.band(); // coherent start for the next column
                    }
                }
                sink = topZ[cols / 2] + nz[cols / 2];
            });
            md << "| " << fmt(h * 100.0, 0) << " cm | " << fmt(width, 1) << " m | " << cols << " | " << fmt(tb * 1e6, 1)
               << " | " << fmt(th * 1e6, 1) << " |\n";
        }
    }

    // ── 4. Applying a rut: K (delta on heights) vs H (CSG on bricks) ────────────────────────────
    md << "\n## 4. Applying a rut (3 m long, 6 cm wide, 3 cm deep) on flat ground\n\n";
    md << "| Spacing | Representation | Edit µs | Bytes added | Gap p99 mm after | Gap max mm after | Normal p99 ° after (2–5 mm) | Query ns after |\n";
    md << "|---|---|---|---|---|---|---|---|\n";
    {
        const Scene& ruts = scenes[1];
        const V3d s0{1.0, 0.35, 0.0}, s1{4.0, 0.35, 0.0};
        const double r = 0.03;
        const auto segDist2D = [&](double x, double y) {
            const double t = std::clamp((x - s0.x) / (s1.x - s0.x), 0.0, 1.0);
            const double px = s0.x + (s1.x - s0.x) * t;
            return std::hypot(x - px, y - s0.y);
        };
        const auto capsule = [&](V3d p) {
            const double t = std::clamp((p.x - s0.x) / (s1.x - s0.x), 0.0, 1.0);
            return length(p - V3d{s0.x + (s1.x - s0.x) * t, s0.y, s0.z}) - r;
        };
        // Exact reference away from the end caps: flat ground with a semicircular groove, as a fine polygon.
        std::vector<V2d> groove = {{-2, -1}, {2, -1}, {2, 0}};
        for (int k = 0; k <= 720; ++k)
        {
            const double a = Pi * k / 720.0; // from the +y rim, under the axis, to the -y rim
            groove.push_back({s0.y + r * std::cos(a), -r * std::sin(a)});
        }
        groove.push_back({-2, 0});
        const auto after = [&](V3d p) { return sdPolygon(groove, {p.y, p.z}); };
        const auto samples = samplePoints(after, {{1.5, 0.28, -0.04}, {3.5, 0.42, 0.02}}, 10000, 99);
        for (float h : spacings)
        {
            // K
            double editK = std::numeric_limits<double>::infinity();
            std::unique_ptr<HeightField> hf;
            for (int rep = 0; rep < 5; ++rep)
            {
                hf = bakeField(ruts.domain, h, ruts.ground);
                const auto t0 = std::chrono::steady_clock::now();
                const int i0 = static_cast<int>(std::floor((s0.x - r) / h)), i1 = static_cast<int>(std::ceil((s1.x + r) / h));
                const int j0 = static_cast<int>(std::floor((s0.y - r + 1.0) / h)), j1 = static_cast<int>(std::ceil((s0.y + r + 1.0) / h));
                for (int j = j0; j <= j1; ++j)
                    for (int i = i0; i <= i1; ++i)
                    {
                        const double rho = segDist2D(hf->nodeX(i), hf->nodeY(j));
                        if (rho < r) hf->lowerTo(i, j, static_cast<float>(-std::sqrt(r * r - rho * rho)));
                    }
                editK = std::min(editK, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
            }
            const auto accK = measure(samples, [&](V3f p) { auto q = hf->query(p); return QHit{q.d, q.n, true}; });
            md << "| " << fmt(h * 100.0, 0) << " cm | B + K delta | " << fmt(editK * 1e6, 1) << " | "
               << bytesText(static_cast<double>(hf->deltaBytes())) << " | " << fmt(accK.gapMm.p99) << " | "
               << fmt(accK.gapMm.max) << " | " << fmt(accK.normalDeg.p99, 1) << " | " << fmt(accK.queryNs, 0) << " |\n";

            // H
            SdfBricks base(h, bandFor(h));
            base.bake(ruts.sdf, ruts.domain);
            double editH = std::numeric_limits<double>::infinity();
            SdfBricks edited = base;
            std::size_t allocated = 0;
            const Box3 region{{s0.x - r - 0.1, s0.y - r - 0.1, -r - 0.1}, {s1.x + r + 0.1, s0.y + r + 0.1, 0.1}};
            for (int rep = 0; rep < 5; ++rep)
            {
                edited = base;
                const auto t0 = std::chrono::steady_clock::now();
                allocated = edited.subtract(capsule, region, [&](V3d c) { return ruts.sdf(c) < 0.0; });
                editH = std::min(editH, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
            }
            const auto accH = measure(samples, [&](V3f p) { auto q = edited.query(p); return QHit{q.d, q.n, q.known}; });
            md << "| " << fmt(h * 100.0, 0) << " cm | H CSG subtract | " << fmt(editH * 1e6, 1) << " | "
               << bytesText(static_cast<double>(edited.bytes() - base.bytes())) << " (" << allocated << " new bricks) | "
               << fmt(accH.gapMm.p99) << " | " << fmt(accH.gapMm.max) << " | " << fmt(accH.normalDeg.p99, 1) << " | "
               << fmt(accH.queryNs, 0) << " |\n";
        }
    }

    // ── 5. Seams between levels ─────────────────────────────────────────────────────────────────
    md << "\n## 5. Seam between a 1 cm page and a 2 cm page\n\n";
    {
        const auto f = [](double x, double y) {
            return 0.05 * std::sin(2 * Pi * x / 0.37) * std::cos(2 * Pi * y / 0.53) + 0.02 * std::sin(2 * Pi * (x + y) / 0.11);
        };
        double crack = 0, coarseInterior = 0;
        for (int k = 1; k < 100; k += 2)
        {
            const double y = k * 0.01;
            crack = std::max(crack, std::abs(f(1.0, y) - 0.5 * (f(1.0, y - 0.01) + f(1.0, y + 0.01))));
        }
        for (int k = 0; k < 50; ++k)
            for (int m = 0; m < 50; ++m)
            {
                const double x = 1.01 + 0.02 * k, y = 0.01 + 0.02 * m; // cell centres of the coarse page
                const double interp = 0.25 * (f(x - 0.01, y - 0.01) + f(x + 0.01, y - 0.01) + f(x - 0.01, y + 0.01) +
                                              f(x + 0.01, y + 0.01));
                coarseInterior = std::max(coarseInterior, std::abs(interp - f(x, y)));
            }
        md << "| Measure | mm |\n|---|---|\n";
        md << "| Crack height along the seam, unstitched (max) | " << fmt(crack * 1000.0) << " |\n";
        md << "| Crack height, stitched (fine edge nodes on the coarse edge) | 0.00 (by construction) |\n";
        md << "| Largest stitch adjustment of a fine edge node | " << fmt(crack * 1000.0) << " |\n";
        md << "| Coarse page's own interpolation error (max), for comparison | " << fmt(coarseInterior * 1000.0) << " |\n";
    }

    // ── 6. Carving a trench through the half-buried rock garden ─────────────────────────────────
    md << "\n## 6. Carving an 8 cm trench through the half-buried rock garden\n\n";
    md << "Errors at points within 30 mm of a rock after the carve: the newly exposed faces.\n\n";
    md << "| Spacing | Representation | Edit ms | Gap p99 mm | Gap max mm | Normal p99 ° (2–5 mm) | Unknown |\n|---|---|---|---|---|---|---|\n";
    {
        const Scene& rocks = scenes[2];
        const Box3 trench{{1.5, -0.5, -0.08}, {4.0, 0.5, 1.0}};
        const auto sdTrench = [&](V3d p) {
            const V3d c = (trench.lo + trench.hi) * 0.5, e = (trench.hi - trench.lo) * 0.5;
            const V3d q{std::abs(p.x - c.x) - e.x, std::abs(p.y - c.y) - e.y, std::abs(p.z - c.z) - e.z};
            const V3d o{std::max(q.x, 0.0), std::max(q.y, 0.0), std::max(q.z, 0.0)};
            return length(o) + std::min(std::max(q.x, std::max(q.y, q.z)), 0.0);
        };
        // Exact soil after the carve: {z < -0.08} united with {z < 0 outside the trench footprint}.
        const auto soil = [&](V3d p) {
            const double dx = std::max(trench.lo.x - p.x, p.x - trench.hi.x);
            const double dy = std::max(trench.lo.y - p.y, p.y - trench.hi.y);
            const double toLeaveFootprint = (dx > 0 || dy > 0) ? 0.0 : -std::max(dx, dy);
            const bool inside = p.z < -0.08 || (p.z < 0 && toLeaveFootprint == 0.0);
            const double d = std::min(std::hypot(std::max(p.z, 0.0), toLeaveFootprint), p.z + 0.08);
            return inside ? -1e-3 : d; // accuracy samples are all outside
        };
        const auto after = [&](V3d p) { return std::min(soil(p), featureSdf(rocks, p)); };
        const auto samples = samplePoints(after, {{1.6, -0.45, -0.08}, {3.7, 0.45, 0.06}}, 10000, 77,
                                          [&](V3d p) { return featureSdf(rocks, p) < 0.03; });
        const auto inTrench = [&](double x, double y) {
            return x >= trench.lo.x && x <= trench.hi.x && y >= trench.lo.y && y <= trench.hi.y;
        };
        for (float h : spacings)
        {
            const auto rowOut = [&](const std::string& rep, double ms, const Accuracy& a) {
                md << "| " << fmt(h * 100.0, 0) << " cm | " << rep << " | " << fmt(ms, 2) << " | " << fmt(a.gapMm.p99) << " | "
                   << fmt(a.gapMm.max) << " | " << fmt(a.normalDeg.p99, 1) << " | " << a.unknown << " |\n";
            };
            // B+E(+K): lower the ground; features untouched and therefore exact.
            std::unique_ptr<HeightField> hf;
            double tBE = std::numeric_limits<double>::infinity();
            for (int rep = 0; rep < 5; ++rep)
            {
                hf = bakeField(rocks.domain, h, rocks.ground);
                const auto t0 = std::chrono::steady_clock::now();
                for (int j = 0; j < hf->ny(); ++j)
                    for (int i = 0; i < hf->nx(); ++i)
                        if (inTrench(hf->nodeX(i), hf->nodeY(j))) hf->lowerTo(i, j, -0.08f);
                tBE = std::min(tBE, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
            }
            HeightFieldWithFeatures be{hf.get(), {}, {}};
            for (auto& b : rocks.boxes) be.boxes.push_back(castObb<float>(b));
            for (auto& m : rocks.meshes) be.meshes.push_back(castMesh<float>(m));
            rowOut("B+E+K: lower heights, features untouched", tBE * 1e3,
                   measure(samples, [&](V3f p) { auto r = be.query(p); return QHit{r.d, r.n, true}; }));

            // B only: rocks baked into heights; rebake the carved region.
            const auto tB0 = std::chrono::steady_clock::now();
            auto all = bakeField(rocks.domain, h, [&](double x, double y) {
                return inTrench(x, y) ? topByRay(after, x, y, rocks.domain.hi.z, rocks.domain.lo.z)
                                      : topByRay(rocks.sdf, x, y, rocks.domain.hi.z, rocks.domain.lo.z);
            });
            const double tB = std::chrono::duration<double>(std::chrono::steady_clock::now() - tB0).count();
            rowOut("B only: rebake heights (rocks were baked in)", tB * 1e3,
                   measure(samples, [&](V3f p) { auto r = all->query(p); return QHit{r.d, r.n, true}; }));

            // H naive: one field for soil and rock; subtracting the trench removes rock too.
            SdfBricks base(h, bandFor(h));
            base.bake(rocks.sdf, rocks.domain);
            SdfBricks naive = base;
            const Box3 region{{trench.lo.x - 0.1, trench.lo.y - 0.1, -0.3}, {trench.hi.x + 0.1, trench.hi.y + 0.1, 0.4}};
            const auto tN0 = std::chrono::steady_clock::now();
            naive.subtract(sdTrench, region, [&](V3d c) { return rocks.sdf(c) < 0.0; });
            const double tN = std::chrono::duration<double>(std::chrono::steady_clock::now() - tN0).count();
            rowOut("H: CSG subtract on a single field", tN * 1e3,
                   measure(samples, [&](V3f p) { auto r = naive.query(p); return QHit{r.d, r.n, r.known}; }));

            // H rebake: resample affected bricks from separated sources (soil carved, rock kept).
            SdfBricks rebaked = base;
            const auto tR0 = std::chrono::steady_clock::now();
            {
                const double bs = SdfBricks::Cells * static_cast<double>(h);
                for (int k = static_cast<int>(std::floor(region.lo.z / bs)); k <= static_cast<int>(std::floor(region.hi.z / bs)); ++k)
                    for (int j = static_cast<int>(std::floor(region.lo.y / bs)); j <= static_cast<int>(std::floor(region.hi.y / bs)); ++j)
                        for (int i = static_cast<int>(std::floor(region.lo.x / bs)); i <= static_cast<int>(std::floor(region.hi.x / bs)); ++i)
                            rebaked.bakeBrick(after, i, j, k, false);
            }
            const double tR = std::chrono::duration<double>(std::chrono::steady_clock::now() - tR0).count();
            rowOut("H: rebake bricks from separate soil and rock sources", tR * 1e3,
                   measure(samples, [&](V3f p) { auto r = rebaked.query(p); return QHit{r.d, r.n, r.known}; }));
        }
    }

    // ── 7. Editing: carve a 10 m curved trail, paint a 2 m stroke ───────────────────────────────
    md << "\n## 7. Editing at 1 cm: carve a 10 m curved trail (1.4 m wide, 10 cm deep), paint a 2 m stroke\n\n";
    {
        const Box3 domain{{0, -1.5, -0.3}, {12, 1.5, 0.3}};
        const float h = 0.01f;
        const auto profile = [](double s) {
            const double a = std::abs(s);
            return a <= 0.6 ? -0.1 : (a < 0.7 ? -0.1 + (a - 0.6) : 0.0);
        };
        const std::vector<V2d> air = {{-0.7, 0.5}, {-0.7, 0.0}, {-0.6, -0.1}, {0.6, -0.1}, {0.7, 0.0}, {0.7, 0.5}};
        const auto radial = [](double x, double y) { return std::hypot(x - 6.0, y + 7.0) - 8.0; };

        double tCarveB = std::numeric_limits<double>::infinity();
        std::unique_ptr<HeightField> hf;
        for (int rep = 0; rep < 3; ++rep)
        {
            hf = bakeField(domain, h, [](double, double) { return 0.0; });
            const auto t0 = std::chrono::steady_clock::now();
            for (int j = 0; j < hf->ny(); ++j)
                for (int i = 0; i < hf->nx(); ++i)
                {
                    const double s = radial(hf->nodeX(i), hf->nodeY(j));
                    if (std::abs(s) < 0.7) hf->lowerTo(i, j, static_cast<float>(profile(s)));
                }
            tCarveB = std::min(tCarveB, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        }

        SdfBricks base(h, bandFor(h));
        base.bake([](V3d p) { return p.z; }, domain);
        double tCarveH = std::numeric_limits<double>::infinity();
        std::size_t newBricks = 0, baseBricks = base.brickCount();
        for (int rep = 0; rep < 3; ++rep)
        {
            SdfBricks e = base;
            const auto t0 = std::chrono::steady_clock::now();
            newBricks = e.subtract([&](V3d p) { return sdPolygon(air, {radial(p.x, p.y), p.z}); },
                                   {{0, -1.5, -0.3}, {12, 1.5, 0.1}}, [](V3d c) { return c.z < 0.0; });
            tCarveH = std::min(tCarveH, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        }

        std::vector<std::uint8_t> channel(static_cast<std::size_t>(hf->nx()) * static_cast<std::size_t>(hf->ny()), 0);
        const double tPaint = bestSeconds(5, [&] {
            const double r = 0.25;
            const int i0 = static_cast<int>((5.0 - r) / h), i1 = static_cast<int>((7.0 + r) / h);
            const int j0 = static_cast<int>((-r + 1.5) / h), j1 = static_cast<int>((r + 1.5) / h);
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i)
                {
                    const double x = hf->nodeX(i), y = hf->nodeY(j);
                    const double px = std::clamp(x, 5.0, 7.0);
                    const double rho = std::hypot(x - px, y);
                    if (rho >= r) continue;
                    const double t = rho / r, w = 1.0 - t * t * (3.0 - 2.0 * t);
                    auto& c = channel[static_cast<std::size_t>(j) * static_cast<std::size_t>(hf->nx()) + static_cast<std::size_t>(i)];
                    c = static_cast<std::uint8_t>(std::lround(c + (200.0 - c) * w));
                }
            sink = channel[channel.size() / 2];
        });

        md << "| Operation | Representation | Time | Notes |\n|---|---|---|---|\n";
        md << "| Carve trail | B + K delta | " << fmt(tCarveB * 1e3, 1) << " ms | " << hf->nx() * hf->ny()
           << " nodes visited; delta " << bytesText(static_cast<double>(hf->deltaBytes())) << " |\n";
        md << "| Carve trail | H CSG subtract | " << fmt(tCarveH * 1e3, 1) << " ms | " << baseBricks << " bricks before, "
           << newBricks << " allocated below the old band |\n";
        md << "| Paint 2 m stroke, 0.25 m radius | B grid u8 channel | " << fmt(tPaint * 1e6, 1)
           << " µs | H has no surface grid to paint; channels would live on B pages anyway |\n";
    }

    std::ofstream(outPath) << md.str();
    std::cout << md.str();
    return 0;
}
