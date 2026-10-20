// Terrain extraction (G3): each library call paired with hand-written code for the same job, on raw arrays with
// the page shape, tile layout and output arrays known at compile time. The twin is the reference; an API gap to it
// is recorded in docs/PERFORMANCE.md. Every pair checks that both sides produce the same output before timing.
#include <algorithm>
#include <array>
#include <bench_runner.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "sub0physics/terrain.hpp"

using namespace sub0physics::terrain;

namespace
{

[[noreturn]] void mismatch(const char* what)
{
    std::fprintf(stderr, "terrain bench: twin and library disagree on %s\n", what);
    std::abort();
}

// ---- Batch height lookup -------------------------------------------------------------------------------------

constexpr int Cells = 256;
constexpr int Nodes = Cells + 1;

struct RawPage
{
    std::vector<std::uint16_t> heights = std::vector<std::uint16_t>(static_cast<std::size_t>(Nodes) * Nodes);
    std::vector<std::uint64_t> known = std::vector<std::uint64_t>((static_cast<std::size_t>(Nodes) * Nodes + 63) / 64, ~std::uint64_t{0});
    std::vector<sub0math::f64::Length> authored = std::vector<sub0math::f64::Length>(heights.size());
    std::vector<std::uint8_t> flags = std::vector<std::uint8_t>(heights.size(), 1U);
    double originZ = 0.0;

    RawPage()
    {
        for (int j = 0; j < Nodes; ++j)
        {
            for (int i = 0; i < Nodes; ++i)
            {
                const double z = 0.2 + 0.05 * std::sin(0.11 * i) * std::cos(0.07 * j);
                authored[static_cast<std::size_t>(j) * Nodes + static_cast<std::size_t>(i)] = sub0math::f64::Length::fromSi(z);
            }
        }
    }
};

const GridSpec& grid256()
{
    static const GridSpec grid = *GridSpec::create({.finestSpacing = sub0math::f64::Length::fromSi(0.01), .pages = {1, 1}});
    return grid;
}

void benchBatchHeight(sub0bench::BenchRunner& runner)
{
    RawPage raw;
    const auto page = *buildHeightPage(grid256(), {0, 0, 0}, SourceRevision{1}, raw.authored, raw.flags);
    // The twin reads the page's own quantised heights, so both sides see identical data.
    for (int j = 0; j < Nodes; ++j)
    {
        const auto row = page.nodeRow(j);
        std::copy(row.begin(), row.end(), raw.heights.begin() + static_cast<std::ptrdiff_t>(j) * Nodes);
    }
    std::vector<float> xs(1024), ys(1024), twinOut(1024);
    std::vector<std::uint8_t> twinKnown(1024);
    std::vector<sub0math::Length> lx(1024), ly(1024), libOut(1024);
    std::vector<std::uint8_t> libKnown(1024);
    for (std::size_t k = 0; k < xs.size(); ++k)
    {
        xs[k] = 0.0025F * static_cast<float>(k % 997);
        ys[k] = 2.55F - 0.0023F * static_cast<float>((k * 7) % 1013);
        lx[k] = sub0math::Length::fromSi(xs[k]);
        ly[k] = sub0math::Length::fromSi(ys[k]);
    }
    const auto twin = [&] {
        const std::uint16_t* h = raw.heights.data();
        const std::uint64_t* mask = raw.known.data();
        std::size_t count = 0;
        for (std::size_t k = 0; k < xs.size(); ++k)
        {
            const float px = xs[k] * 100.0F, py = ys[k] * 100.0F;
            twinKnown[k] = 0U;
            if (!(px >= 0.0F && py >= 0.0F && px <= 256.0F && py <= 256.0F)) continue;
            const int i = std::min(static_cast<int>(px), Cells - 1), j = std::min(static_cast<int>(py), Cells - 1);
            const float fx = px - static_cast<float>(i), fy = py - static_cast<float>(j);
            const std::size_t a = static_cast<std::size_t>(j) * Nodes + static_cast<std::size_t>(i);
            const std::size_t d = a + Nodes + 1U;
            const auto bit = [mask](std::size_t n) { return ((mask[n >> 6] >> (n & 63)) & 1U) != 0U; };
            float z;
            if (fx >= fy)
            {
                const std::size_t b = a + 1U;
                if (!bit(a) || !bit(b) || !bit(d)) continue;
                const float za = h[a], zb = h[b], zd = h[d];
                z = (za + (zb - za) * fx + (zd - zb) * fy) * 1e-4F;
            }
            else
            {
                const std::size_t c = a + Nodes;
                if (!bit(a) || !bit(d) || !bit(c)) continue;
                const float za = h[a], zd = h[d], zc = h[c];
                z = (za + (zd - zc) * fx + (zc - za) * fy) * 1e-4F;
            }
            twinOut[k] = z;
            twinKnown[k] = 1U;
            ++count;
        }
        return count;
    };
    if (twin() != page.heightsAt(lx, ly, libOut, libKnown))
    {
        mismatch("batch height count");
    }
    for (std::size_t k = 0; k < xs.size(); ++k)
    {
        if (twinKnown[k] != libKnown[k] || (twinKnown[k] != 0U && twinOut[k] != libOut[k].si())) mismatch("batch height");
    }
    runner.run("page heights batch x1024, hand-written raw arrays", [&] { ankerl::nanobench::doNotOptimizeAway(twin()); });
    runner.run("page heights batch x1024, HeightPage::heightsAt",
               [&] { ankerl::nanobench::doNotOptimizeAway(page.heightsAt(lx, ly, libOut, libKnown)); });
}

// ---- Height extraction ---------------------------------------------------------------------------------------

// The job: a 128 x 64-cell block of one page into node heights, cell states and triangles with normals.
void benchHeightExtract(sub0bench::BenchRunner& runner)
{
    RawPage raw;
    const auto page = *buildHeightPage(grid256(), {0, 0, 0}, SourceRevision{1}, raw.authored, raw.flags);
    for (int j = 0; j < Nodes; ++j)
    {
        const auto row = page.nodeRow(j);
        std::copy(row.begin(), row.end(), raw.heights.begin() + static_cast<std::ptrdiff_t>(j) * Nodes);
    }
    const std::array<const HeightPage*, 1> pointers{&page};
    const PageWindow window{{{0, 0}, {0, 0}}, 0, pointers};
    HeightRequest request;
    request.cells = {{20, 30}, {147, 93}};
    HeightExtract out;
    const auto r = extractHeights(grid256(), window, request, out);
    if (r.status != HeightStatus::Complete) mismatch("extract status");

    constexpr int cx = 128, cy = 64, lo0 = 20, lo1 = 30;
    std::vector<float> nodeZ((cx + 1) * (cy + 1));
    std::vector<std::uint8_t> nodeKnown(nodeZ.size()), cellState(cx * cy);
    std::vector<std::uint32_t> cellTriangle(cx * cy), triangleCell(2 * cx * cy);
    std::vector<float> triangles(2 * cx * cy * 9), normals(2 * cx * cy * 3);
    const float step = 1e-4F, hh = static_cast<float>(grid256().spacing(0).si());
    const auto twin = [&] {
        const std::uint16_t* h = raw.heights.data();
        const std::uint64_t* mask = raw.known.data();
        const auto bit = [mask](std::size_t n) { return ((mask[n >> 6] >> (n & 63)) & 1U) != 0U; };
        for (int b = 0; b <= cy; ++b)
        {
            const std::size_t row = static_cast<std::size_t>(lo1 + b) * Nodes + lo0;
            for (int a = 0; a <= cx; ++a)
            {
                const std::size_t at = static_cast<std::size_t>(b) * (cx + 1) + static_cast<std::size_t>(a);
                const bool k = bit(row + static_cast<std::size_t>(a));
                nodeKnown[at] = k ? 1U : 0U;
                nodeZ[at] = k ? static_cast<float>(h[row + static_cast<std::size_t>(a)]) * step : 0.0F;
            }
        }
        std::size_t tri = 0;
        for (int b = 0; b < cy; ++b)
        {
            for (int a = 0; a < cx; ++a)
            {
                const std::size_t cell = static_cast<std::size_t>(b) * cx + static_cast<std::size_t>(a);
                const std::size_t n0 = static_cast<std::size_t>(b) * (cx + 1) + static_cast<std::size_t>(a), n1 = n0 + (cx + 1);
                if (!nodeKnown[n0] || !nodeKnown[n0 + 1] || !nodeKnown[n1] || !nodeKnown[n1 + 1])
                {
                    cellState[cell] = 2U;
                    cellTriangle[cell] = 0xFFFFFFFFU;
                    continue;
                }
                cellState[cell] = 1U;
                cellTriangle[cell] = static_cast<std::uint32_t>(tri);
                const float za = nodeZ[n0], zb = nodeZ[n0 + 1], zc = nodeZ[n1], zd = nodeZ[n1 + 1];
                const float x0 = static_cast<float>(a) * hh, x1 = static_cast<float>(a + 1) * hh;
                const float y0 = static_cast<float>(b) * hh, y1 = static_cast<float>(b + 1) * hh;
                float* t = &triangles[tri * 9];
                t[0] = x0; t[1] = y0; t[2] = za; t[3] = x1; t[4] = y0; t[5] = zb; t[6] = x1; t[7] = y1; t[8] = zd;
                t[9] = x0; t[10] = y0; t[11] = za; t[12] = x1; t[13] = y1; t[14] = zd; t[15] = x0; t[16] = y1; t[17] = zc;
                float* nrm = &normals[tri * 3];
                const float n0x = -(zb - za), n0y = zb - zd, i0 = 1.0F / std::sqrt(n0x * n0x + n0y * n0y + hh * hh);
                nrm[0] = n0x * i0; nrm[1] = n0y * i0; nrm[2] = hh * i0;
                const float n1x = zc - zd, n1y = za - zc, i1 = 1.0F / std::sqrt(n1x * n1x + n1y * n1y + hh * hh);
                nrm[3] = n1x * i1; nrm[4] = n1y * i1; nrm[5] = hh * i1;
                triangleCell[tri] = triangleCell[tri + 1] = static_cast<std::uint32_t>(cell);
                tri += 2;
            }
        }
        return tri;
    };
    if (twin() != out.triangles().size()) mismatch("extract triangle count");
    for (std::size_t n = 0; n < nodeZ.size(); ++n)
    {
        if (nodeZ[n] != out.z()[n].si() || nodeKnown[n] != out.nodeKnown()[n]) mismatch("extract nodes");
    }
    for (std::size_t t = 0; t < out.triangles().size(); ++t)
    {
        const auto& lt = out.triangles()[t];
        const float* tt = &triangles[t * 9];
        if (tt[0] != lt.a.x.si() || tt[2] != lt.a.z.si() || tt[5] != lt.b.z.si() || tt[8] != lt.c.z.si() || tt[7] != lt.c.y.si()) mismatch("extract triangles");
        if (normals[t * 3] != out.triangleNormal()[t].x || normals[t * 3 + 2] != out.triangleNormal()[t].z) mismatch("extract normals");
    }
    runner.bench().minEpochIterations(200);
    runner.run("height extract 128x64 cells, hand-written raw arrays", [&] { ankerl::nanobench::doNotOptimizeAway(twin()); });
    runner.run("height extract 128x64 cells, extractHeights", [&] {
        ankerl::nanobench::doNotOptimizeAway(extractHeights(grid256(), window, request, out).deliveredCells);
    });
    // The same block at stride 4: the resolution lever, recorded (no twin).
    HeightRequest coarse = request;
    coarse.stride = 4;
    runner.run("height extract 128x64 cells at stride 4, extractHeights", [&] {
        ankerl::nanobench::doNotOptimizeAway(extractHeights(grid256(), window, coarse, out).deliveredCells);
    });
    // Rays against the 8192-cell patch: the cell walk against testing every triangle (no twin; the pair is information).
    (void)extractHeights(grid256(), window, request, out);
    const sub0math::Ray<sub0math::Length, 3> ray{{sub0math::Length::fromSi(0.05F), sub0math::Length::fromSi(0.1F), sub0math::Length::fromSi(1.0F)},
                                                  {0.9F, 0.2F, -0.3835F}};
    runner.bench().minEpochIterations(2000);
    runner.run("ray cast against 8192 cells, castRay", [&] {
        ankerl::nanobench::doNotOptimizeAway(castRay(out, ray, sub0math::Length::fromSi(3.0F)).hit);
    });
    runner.bench().minEpochIterations(20);
    runner.run("ray cast against 8192 cells, every triangle", [&] {
        ankerl::nanobench::doNotOptimizeAway(castRayBruteForce(out, ray, sub0math::Length::fromSi(3.0F)).hit);
    });
}

// ---- Column extraction (R0a shape) ---------------------------------------------------------------------------

struct BenchTile
{
    static constexpr std::size_t Cells = 1600;
    std::array<double, Cells * 7> values{};
    std::array<std::uint8_t, Cells * 3> materials{};
    std::array<std::uint8_t, Cells> count{};
    std::array<std::uint8_t, Cells> valid{};

    [[nodiscard]] bool known(std::size_t c) const { return valid[c] != 0U; }
    [[nodiscard]] sub0math::f64::Length top(std::size_t c) const { return sub0math::f64::Length::fromSi(values[7 * c]); }
    [[nodiscard]] sub0math::Vec3d normal(std::size_t c) const { return {values[7 * c + 1], values[7 * c + 2], values[7 * c + 3]}; }
    [[nodiscard]] std::size_t spanCount(std::size_t c) const { return count[c]; }
    [[nodiscard]] ColumnSpan span(std::size_t c, std::size_t k) const
    {
        const double* v = &values[7 * c];
        return {sub0math::f64::Length::fromSi(k != 0U ? v[3 + k] : v[0]), sub0math::f64::Length::fromSi(v[4 + k]), materials[3 * c + k]};
    }
};

void benchColumns(sub0bench::BenchRunner& runner)
{
    const auto grid = *GridSpec::create({.finestSpacing = sub0math::f64::Length::fromSi(0.02),
                                         .cellsPerSide = 40,
                                         .registration = Registration::CellCentre,
                                         .anchor = {sub0math::f64::Length{}, sub0math::f64::Length::fromSi(-1.6)},
                                         .pages = {10, 4}});
    std::vector<BenchTile> tiles(40);
    for (std::size_t t = 0; t < tiles.size(); ++t)
    {
        for (std::size_t c = 0; c < BenchTile::Cells; ++c)
        {
            double* v = &tiles[t].values[7 * c];
            v[0] = 0.02 + 0.001 * static_cast<double>((c + t) % 7);
            v[1] = 0.0; v[2] = 0.0; v[3] = 1.0;
            v[4] = -0.01; v[5] = -0.2; v[6] = 0.0;
            tiles[t].count[c] = (c % 5 == 0) ? 2U : 1U;
            tiles[t].materials[3 * c] = 1; tiles[t].materials[3 * c + 1] = 2; tiles[t].materials[3 * c + 2] = 0;
            tiles[t].valid[c] = 1U;
        }
    }
    const sub0math::GridRange<2> block{{67, 58}, {90, 71}}; // 24 x 14 cells over four tiles
    const auto lookup = [&](const TileKey& k) -> const BenchTile* { return &tiles[static_cast<std::size_t>(k.j * 10 + k.i)]; };
    ColumnExtract out;
    if (extractColumns<BenchTile>(grid, 0, block, lookup, out).status != ColumnStatus::Complete) mismatch("column status");

    const std::size_t n = 24 * 14;
    std::vector<std::uint32_t> ids(n), offsets(n + 1);
    std::vector<double> xs(n), ys(n), tops(n), normals(n * 3), uppers(n * 3), lowers(n * 3);
    std::vector<std::uint8_t> mats(n * 3);
    const auto twin = [&] {
        constexpr int T = 40;
        std::uint32_t missing = 0, spans = 0;
        for (int y = 58; y <= 71; ++y)
        {
            for (int x = 67; x <= 90;)
            {
                const int end = std::min(90, (x / T + 1) * T - 1);
                const BenchTile& t = tiles[static_cast<std::size_t>((y / T) * 10 + x / T)];
                std::size_t index = static_cast<std::size_t>((y % T) * T + x % T);
                for (int i = x; i <= end; ++i, ++index)
                {
                    missing += t.valid[index] == 0U ? 1U : 0U;
                    spans += t.count[index];
                }
                x = end + 1;
            }
        }
        if (missing != 0U || spans > 12288U) return std::size_t{0};
        std::size_t j = 0;
        std::uint32_t s = 0;
        for (int y = 58; y <= 71; ++y)
        {
            const double cy = -1.6 + 0.02 * (y + 0.5);
            for (int x = 67; x <= 90;)
            {
                const int end = std::min(90, (x / T + 1) * T - 1);
                const BenchTile& t = tiles[static_cast<std::size_t>((y / T) * 10 + x / T)];
                std::size_t index = static_cast<std::size_t>((y % T) * T + x % T);
                for (int i = x; i <= end; ++i, ++index, ++j)
                {
                    const double* v = &t.values[7 * index];
                    ids[j] = static_cast<std::uint32_t>(y * 400 + i);
                    xs[j] = 0.0 + 0.02 * (i + 0.5);
                    ys[j] = cy;
                    tops[j] = v[0];
                    normals[3 * j] = v[1]; normals[3 * j + 1] = v[2]; normals[3 * j + 2] = v[3];
                    offsets[j] = s;
                    for (std::uint8_t k = 0; k < t.count[index]; ++k, ++s)
                    {
                        uppers[s] = k != 0U ? v[3 + k] : v[0];
                        lowers[s] = v[4 + k];
                        mats[s] = t.materials[3 * index + k];
                    }
                }
                x = end + 1;
            }
        }
        offsets[n] = s;
        return j;
    };
    if (twin() != n) mismatch("column count");
    for (std::size_t c = 0; c < n; ++c)
    {
        if (ids[c] != out.ids()[c] || xs[c] != out.x()[c].si() || ys[c] != out.y()[c].si() || tops[c] != out.top()[c].si() ||
            offsets[c] != out.offset()[c])
            mismatch("columns");
    }
    for (std::size_t s = 0; s < out.upper().size(); ++s)
    {
        if (uppers[s] != out.upper()[s].si() || lowers[s] != out.lower()[s].si() || mats[s] != out.material()[s]) mismatch("spans");
    }
    runner.bench().minEpochIterations(2000);
    runner.run("column extract 336 columns, 4 tiles, hand-written", [&] { ankerl::nanobench::doNotOptimizeAway(twin()); });
    runner.run("column extract 336 columns, 4 tiles, extractColumns", [&] {
        ankerl::nanobench::doNotOptimizeAway(extractColumns<BenchTile>(grid, 0, block, lookup, out).columns);
    });
}

// ---- Tile cache ----------------------------------------------------------------------------------------------

struct EmptyTile
{
    int value = 0;
};

class InstantSource final : public TileSource<EmptyTile>
{
public:
    LoadResult load(const TileKey& key, EmptyTile& slot) noexcept override
    {
        slot.value = key.i;
        return {LoadState::Ready, {}};
    }
};

// The twin keeps resident keys in an array, oldest first, as R0a's cache does: a linear find and a rotate.
struct ArrayLru
{
    static constexpr std::size_t Limit = 32;
    std::array<std::int32_t, Limit> keys{};
    std::size_t count = 0;
    std::uint64_t hits = 0, misses = 0, evictions = 0;

    void access(std::int32_t key) noexcept
    {
        const auto end = keys.begin() + static_cast<std::ptrdiff_t>(count);
        const auto it = std::find(keys.begin(), end, key);
        if (it != end)
        {
            ++hits;
            std::rotate(it, it + 1, end);
            return;
        }
        ++misses;
        if (count == Limit)
        {
            std::move(keys.begin() + 1, keys.end(), keys.begin());
            --count;
            ++evictions;
        }
        keys[count++] = key;
    }
};

void benchCache(sub0bench::BenchRunner& runner)
{
    // 1024 accesses over 48 tiles with a moving working set, so there are hits, misses and evictions.
    std::vector<std::int32_t> sequence(1024);
    std::uint32_t state = 12345;
    for (std::size_t k = 0; k < sequence.size(); ++k)
    {
        state = state * 1664525U + 1013904223U;
        sequence[k] = static_cast<std::int32_t>((k / 16 + ((state >> 16) % 24)) % 48);
    }
    InstantSource source;
    auto cache = *SourceTileCache<EmptyTile>::create({.residentTiles = 32, .pendingLoads = 8, .extraSlots = 4}, &source);
    ArrayLru twin;
    for (int round = 0; round < 3; ++round)
    {
        for (const std::int32_t key : sequence)
        {
            twin.access(key);
            if (cache.access({key, 0, 0}).state == TileAccess::Started)
            {
                cache.pump([](const LoadEvent&) {});
            }
        }
    }
    if (twin.hits != cache.stats().hits || twin.misses != cache.stats().misses || twin.evictions != cache.stats().evictions) mismatch("cache counters");
    runner.bench().minEpochIterations(2000);
    runner.run("tile cache x1024 accesses, hand-written array LRU", [&] {
        for (const std::int32_t key : sequence) twin.access(key);
        ankerl::nanobench::doNotOptimizeAway(twin.hits);
    });
    runner.run("tile cache x1024 accesses, SourceTileCache", [&] {
        for (const std::int32_t key : sequence)
        {
            if (cache.access({key, 0, 0}).state == TileAccess::Started)
            {
                cache.pump([](const LoadEvent&) {});
            }
        }
        ankerl::nanobench::doNotOptimizeAway(cache.stats().hits);
    });
}

// ---- Corridor extraction -------------------------------------------------------------------------------------

void benchCorridor(sub0bench::BenchRunner& runner)
{
    RawPage raw;
    const auto page = *buildHeightPage(grid256(), {0, 0, 0}, SourceRevision{1}, raw.authored, raw.flags);
    const std::array<const HeightPage*, 1> pointers{&page};
    const PageWindow window{{{0, 0}, {0, 0}}, 0, pointers};
    const auto corridor = *CorridorBuilder::create({{sub0math::f64::Length::fromSi(0.3), sub0math::f64::Length::fromSi(0.5), {}},
                                                    {sub0math::f64::Length::fromSi(0.5), sub0math::f64::Length::fromSi(0.6), sub0math::f64::Angle::fromSi(0.3)},
                                                    {sub0math::f64::Length::fromSi(0.38), sub0math::f64::Length::fromSi(0.08), sub0math::f64::Length::fromSi(0.1)},
                                                    8});
    HeightExtract out;
    if (extractCorridor(grid256(), corridor, window, 0, 1, out).status != HeightStatus::Complete) mismatch("corridor status");
    runner.bench().minEpochIterations(200);
    runner.run("corridor extract, 8 boxes, 1 cm, extractCorridor", [&] {
        ankerl::nanobench::doNotOptimizeAway(extractCorridor(grid256(), corridor, window, 0, 1, out).deliveredCells);
    });
}

} // namespace

void benchTerrain(sub0bench::BenchRunner& runner)
{
    benchBatchHeight(runner);
    benchHeightExtract(runner);
    benchColumns(runner);
    benchCache(runner);
    benchCorridor(runner);
}
