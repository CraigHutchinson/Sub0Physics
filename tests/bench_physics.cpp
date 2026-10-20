#define ANKERL_NANOBENCH_IMPLEMENT
#include <algorithm>
#include <array>
#include <bench_runner.hpp>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <expected>
#include <limits>
#include <span>
#include <sub0math/units.hpp>
#include <vector>

#include "sub0physics/sub0physics.hpp"

using namespace sub0physics;

void benchMultibody(sub0bench::BenchRunner& runner); // bench_multibody.cpp

namespace
{

// A loaded wheel on a clipped contact: the shape of BikeLab's first station.
struct ContactOscillator
{
    UnilateralSpringDamper contact =
        *UnilateralSpringDamper::create({.stiffness = sub0math::newtonsPerMetre(90000.0), .damping = sub0math::newtonSecondsPerMetre(400.0)});
    double mass = 45.0;
    double radius = 0.36;

    std::expected<void, Fault> evaluate(double, const StateVector<2>& x, StateVector<2>& dxdt) const
    {
        const ContactResponse r = contact.respond(sub0math::metres(radius - x[0]), sub0math::metresPerSecond(-x[1]));
        dxdt[0] = x[1];
        dxdt[1] = -r.force.si() / mass - StandardGravity; // force on the penetration radius - x[0]
        return {};
    }
};

template<typename Integrator>
void benchStep(sub0bench::BenchRunner& runner, const char* name)
{
    ContactOscillator system;
    Integrator integrator;
    StateVector<2> x{0.36, 0.0};
    std::uint64_t tick = 0;
    runner.run(name, [&] {
        if (!integrator.stepTick(system, tick++, sub0math::seconds(0.0005), x))
        {
            x = {0.36, 0.0};
        }
        ankerl::nanobench::doNotOptimizeAway(x);
    });
}

// Height lookup on a terrain page, against the fastest hand-written form of the same job: 16-bit heights,
// a known-node bit mask and the fixed-diagonal triangulation, read from raw arrays.
void benchPageHeight(sub0bench::BenchRunner& runner)
{
    using namespace sub0physics::terrain;
    constexpr int nodes = HeightPage::Nodes;
    constexpr std::size_t count = static_cast<std::size_t>(nodes) * nodes;
    const GridSpec grid = *GridSpec::create(sub0math::f64::Length::fromSi(0.01), 4);
    std::vector<sub0math::f64::Length> heights(count);
    std::vector<std::uint8_t> known(count, 1);
    std::vector<std::uint16_t> raw(count);
    std::vector<std::uint64_t> mask((count + 63) / 64, ~std::uint64_t{0});
    for (int j = 0; j < nodes; ++j)
        for (int i = 0; i < nodes; ++i)
        {
            const double z = 0.2 + 0.05 * std::sin(0.11 * i) * std::cos(0.07 * j);
            const auto n = static_cast<std::size_t>(j) * nodes + static_cast<std::size_t>(i);
            heights[n] = sub0math::f64::Length::fromSi(z);
            raw[n] = static_cast<std::uint16_t>(std::lround((z - 0.15) / 1e-4));
        }
    const auto page = *buildHeightPage(grid, {0, 0, 0}, SourceRevision{1}, heights, known);
    std::vector<float> xs(1024), ys(1024);
    for (std::size_t k = 0; k < xs.size(); ++k)
    {
        xs[k] = 0.0025F * static_cast<float>(k % 997);
        ys[k] = 2.55F - 0.0023F * static_cast<float>((k * 7) % 1013);
    }

    runner.run("page height x1024, hand-written raw arrays", [&] {
        float sum = 0.0F;
        const auto bit = [&](int i, int j) {
            const auto n = static_cast<std::size_t>(j) * nodes + static_cast<std::size_t>(i);
            return ((mask[n >> 6] >> (n & 63)) & 1U) != 0U;
        };
        const auto at = [&](int i, int j) { return static_cast<float>(raw[static_cast<std::size_t>(j) * nodes + static_cast<std::size_t>(i)]); };
        for (std::size_t k = 0; k < xs.size(); ++k)
        {
            const float lx = xs[k] * 100.0F, ly = ys[k] * 100.0F;
            if (!(lx >= 0.0F && ly >= 0.0F && lx <= 256.0F && ly <= 256.0F)) continue;
            const int i = std::min(static_cast<int>(lx), 255), j = std::min(static_cast<int>(ly), 255);
            const float fx = lx - static_cast<float>(i), fy = ly - static_cast<float>(j);
            if (fx >= fy)
            {
                if (!bit(i, j) || !bit(i + 1, j) || !bit(i + 1, j + 1)) continue;
                sum += (at(i, j) + (at(i + 1, j) - at(i, j)) * fx + (at(i + 1, j + 1) - at(i + 1, j)) * fy) * 1e-4F;
            }
            else
            {
                if (!bit(i, j) || !bit(i + 1, j + 1) || !bit(i, j + 1)) continue;
                sum += (at(i, j) + (at(i + 1, j + 1) - at(i, j + 1)) * fx + (at(i, j + 1) - at(i, j)) * fy) * 1e-4F;
            }
        }
        ankerl::nanobench::doNotOptimizeAway(sum);
    });
    runner.run("page height x1024, HeightPage::heightAt", [&] {
        float sum = 0.0F;
        for (std::size_t k = 0; k < xs.size(); ++k)
        {
            if (const auto z = page.heightAt(sub0math::Length::fromSi(xs[k]), sub0math::Length::fromSi(ys[k]))) sum += z.height.si();
        }
        ankerl::nanobench::doNotOptimizeAway(sum);
    });
}

// Page height with the cell count read at run time, against the same hand-written kernel with 256 as a
// compile-time constant: what a runtime-sized HeightPage would cost on its lookup (the decision in docs/TERRAIN.md).
void benchPageSize(sub0bench::BenchRunner& runner)
{
    constexpr int fixedCells = 256;
    std::vector<std::uint16_t> raw(static_cast<std::size_t>(fixedCells + 1) * (fixedCells + 1));
    std::vector<std::uint64_t> mask((raw.size() + 63) / 64, ~std::uint64_t{0});
    for (std::size_t n = 0; n < raw.size(); ++n) raw[n] = static_cast<std::uint16_t>(2000 + (n * 37) % 500);
    std::vector<float> xs(1024), ys(1024);
    for (std::size_t k = 0; k < xs.size(); ++k)
    {
        xs[k] = 0.0025F * static_cast<float>(k % 997);
        ys[k] = 2.55F - 0.0023F * static_cast<float>((k * 7) % 1013);
    }
    static volatile int runtimeCells = fixedCells;
    const auto kernel = [&](auto cellsOf) {
        const int cells = cellsOf();
        const int nodes = cells + 1;
        const float inverse = 100.0F;
        float sum = 0.0F;
        const auto bit = [&](int i, int j) {
            const auto n = static_cast<std::size_t>(j) * static_cast<std::size_t>(nodes) + static_cast<std::size_t>(i);
            return ((mask[n >> 6] >> (n & 63)) & 1U) != 0U;
        };
        const auto at = [&](int i, int j) {
            return static_cast<float>(raw[static_cast<std::size_t>(j) * static_cast<std::size_t>(nodes) + static_cast<std::size_t>(i)]);
        };
        for (std::size_t k = 0; k < xs.size(); ++k)
        {
            const float lx = xs[k] * inverse, ly = ys[k] * inverse;
            if (!(lx >= 0.0F && ly >= 0.0F && lx <= static_cast<float>(cells) && ly <= static_cast<float>(cells))) continue;
            const int i = std::min(static_cast<int>(lx), cells - 1), j = std::min(static_cast<int>(ly), cells - 1);
            const float fx = lx - static_cast<float>(i), fy = ly - static_cast<float>(j);
            if (fx >= fy)
            {
                if (!bit(i, j) || !bit(i + 1, j) || !bit(i + 1, j + 1)) continue;
                sum += (at(i, j) + (at(i + 1, j) - at(i, j)) * fx + (at(i + 1, j + 1) - at(i + 1, j)) * fy) * 1e-4F;
            }
            else
            {
                if (!bit(i, j) || !bit(i + 1, j + 1) || !bit(i, j + 1)) continue;
                sum += (at(i, j) + (at(i + 1, j + 1) - at(i, j + 1)) * fx + (at(i, j + 1) - at(i, j)) * fy) * 1e-4F;
            }
        }
        return sum;
    };
    runner.run("page size x1024, 256 cells compile-time", [&] {
        ankerl::nanobench::doNotOptimizeAway(kernel([] { return fixedCells; }));
    });
    runner.run("page size x1024, 256 cells at run time", [&] {
        ankerl::nanobench::doNotOptimizeAway(kernel([] { return static_cast<int>(runtimeCells); }));
    });
}

// Circle clearance on a 2048-point profile over 256 queries, against a hand-written twin that does the same job
// on prepared structure-of-arrays data (per-segment extents and reciprocal squared lengths): a window from the
// segment under the centre, four independent minima, ties into a fixed array. The twin is the reference; the
// API's gap to it is recorded in docs/PERFORMANCE.md.
struct ClearanceTwin
{
    std::vector<double> x, y, dx, dy, inverseLength2;
    struct Out
    {
        int status = 0; // 1 = available
        double distance = 0.0, gap = 0.0;
        std::size_t witnesses = 0;
        std::array<double, 8> px{}, py{}, nx{}, ny{};
        std::array<std::size_t, 8> index{};
        std::array<bool, 8> vertex{};
    };

    explicit ClearanceTwin(std::span<const ProfilePoint> points)
    {
        for (const ProfilePoint& p : points)
        {
            x.push_back(p.x.si());
            y.push_back(p.y.si());
        }
        for (std::size_t i = 0; i + 1U < x.size(); ++i)
        {
            dx.push_back(x[i + 1U] - x[i]);
            dy.push_back(y[i + 1U] - y[i]);
            inverseLength2.push_back(1.0 / (dx[i] * dx[i] + dy[i] * dy[i]));
        }
    }

    Out measure(double cx, double cy, double radius) const noexcept
    {
        constexpr double Tie = 1e-10;
        Out out;
        if (cx - radius < x.front() || cx + radius > x.back())
        {
            return out;
        }
        const std::size_t segments = x.size() - 1U;
        const double* const xs = x.data();
        const auto it = std::lower_bound(xs + 1, xs + x.size(), cx, [](double px, double v) { return px < v; });
        const std::size_t start = std::min(static_cast<std::size_t>(it - xs), x.size() - 1U) - 1U;
        const auto squared = [&](std::size_t i) {
            const double ax = cx - x[i];
            const double ay = cy - y[i];
            const double t = std::clamp((ax * dx[i] + ay * dy[i]) * inverseLength2[i], 0.0, 1.0);
            const double ex = ax - dx[i] * t;
            const double ey = ay - dy[i] * t;
            return ex * ex + ey * ey;
        };
        // PERF: window from centre segment, 4 minima; faster than pruned walk.
        const double bound = std::sqrt(squared(start));
        const std::size_t lo = static_cast<std::size_t>(std::lower_bound(xs, xs + x.size(), cx - bound) - xs);
        const std::size_t hi = static_cast<std::size_t>(std::upper_bound(xs, xs + x.size(), cx + bound) - xs);
        std::size_t i = lo == 0U ? 0U : lo - 1U;
        const std::size_t end = std::min(hi, segments);
        double m0 = std::numeric_limits<double>::infinity(), m1 = m0, m2 = m0, m3 = m0;
        for (; i + 4U <= end; i += 4U)
        {
            m0 = std::min(m0, squared(i));
            m1 = std::min(m1, squared(i + 1U));
            m2 = std::min(m2, squared(i + 2U));
            m3 = std::min(m3, squared(i + 3U));
        }
        for (; i < end; ++i)
        {
            m0 = std::min(m0, squared(i));
        }
        const double distance = std::sqrt(std::min(std::min(m0, m1), std::min(m2, m3)));
        if (distance <= Tie)
        {
            return out;
        }
        double surface = -std::numeric_limits<double>::infinity();
        {
            const auto first = std::lower_bound(xs, xs + x.size(), cx, [](double px, double v) { return px < v; });
            std::size_t k = first == xs ? 0U : static_cast<std::size_t>(first - xs) - 1U;
            for (; k + 1U < x.size() && xs[k] <= cx; ++k)
            {
                if (cx > xs[k + 1U]) continue;
                const double h = dx[k] == 0.0 ? std::max(y[k], y[k + 1U]) : y[k] + (cx - xs[k]) * dy[k] / dx[k];
                surface = std::max(surface, h);
            }
        }
        if (cy < surface)
        {
            return out;
        }
        const double reach = distance + Tie;
        const double reach2 = reach * reach;
        std::size_t first = start;
        while (first > 0U && cx - xs[first] <= reach) --first;
        for (std::size_t k = first; k < segments && xs[k] - cx <= reach; ++k)
        {
            const double ax = cx - x[k], ay = cy - y[k];
            const double t = std::clamp((ax * dx[k] + ay * dy[k]) * inverseLength2[k], 0.0, 1.0);
            const double px = x[k] + dx[k] * t, py = y[k] + dy[k] * t;
            const double ex = cx - px, ey = cy - py;
            const double d2 = ex * ex + ey * ey;
            if (d2 > reach2) continue;
            const bool vertex = t == 0.0 || t == 1.0;
            const std::size_t index = t == 1.0 ? k + 1U : k;
            if (out.witnesses != 0U && out.witnesses <= out.index.size() && out.vertex[out.witnesses - 1U] == vertex &&
                out.index[out.witnesses - 1U] == index)
            {
                continue;
            }
            if (out.witnesses < out.index.size())
            {
                const double d = std::sqrt(d2);
                out.px[out.witnesses] = px;
                out.py[out.witnesses] = py;
                out.nx[out.witnesses] = ex / d;
                out.ny[out.witnesses] = ey / d;
                out.vertex[out.witnesses] = vertex;
                out.index[out.witnesses] = index;
            }
            ++out.witnesses;
        }
        out.status = 1;
        out.distance = distance;
        out.gap = distance - radius;
        return out;
    }
};

void benchClearance(sub0bench::BenchRunner& runner, const PolylineProfile& profile)
{
    const ClearanceTwin twin(profile.points());
    std::vector<double> cx(256), cy(256);
    for (std::size_t k = 0; k < cx.size(); ++k)
    {
        cx[k] = -9.0 + 18.0 * static_cast<double>((k * 97) % 256) / 256.0;
        cy[k] = profile.surfaceHeight(sub0math::metres(cx[k])).si() + 0.05 + 0.6 * static_cast<double>((k * 61) % 256) / 256.0;
    }
    const sub0math::f64::Length tie = sub0math::metres(1e-10), radius = sub0math::metres(0.36);
    // The twin must agree with the API before its time means anything.
    std::array<Witness2, DefaultWitnessCapacity> witnesses;
    for (std::size_t k = 0; k < cx.size(); ++k)
    {
        const CircleClearance api = measureCircleClearance(profile, {{{sub0math::metres(cx[k]), sub0math::metres(cy[k])}, radius}, tie}, witnesses);
        const ClearanceTwin::Out ref = twin.measure(cx[k], cy[k], radius.si());
        const bool same = api.status == ClearanceStatus::Available && ref.status == 1 &&
                          std::abs(api.distance.si() - ref.distance) <= 1e-12 && api.witnessCount == std::min(ref.witnesses, ref.index.size()) &&
                          witnesses[0].index == ref.index[0];
        if (!same)
        {
            std::fprintf(stderr, "clearance twin disagrees with the API at query %zu\n", k);
            std::abort();
        }
    }

    runner.run("clearance x256, hand-written twin", [&] {
        double sum = 0.0;
        for (std::size_t k = 0; k < cx.size(); ++k)
        {
            const ClearanceTwin::Out o = twin.measure(cx[k], cy[k], radius.si());
            sum += o.gap + o.nx[0] + static_cast<double>(o.witnesses);
        }
        ankerl::nanobench::doNotOptimizeAway(sum);
    });
    runner.run("clearance x256, measureCircleClearance", [&] {
        double sum = 0.0;
        std::array<Witness2, DefaultWitnessCapacity> w;
        for (std::size_t k = 0; k < cx.size(); ++k)
        {
            const CircleClearance c = measureCircleClearance(profile, {{{sub0math::metres(cx[k]), sub0math::metres(cy[k])}, radius}, tie}, w);
            sum += c.gap.si() + w[0].normal.x + static_cast<double>(c.witnessCount);
        }
        ankerl::nanobench::doNotOptimizeAway(sum);
    });
}

// Page and cell lookup through GridSpec (RegularGrid ownership), against hand-written floor and range tests on the
// same layout, and the corridor rasteriser for a typical wheel step (no twin: it is new; the number is recorded).
void benchGrid(sub0bench::BenchRunner& runner)
{
    using namespace sub0physics::terrain;
    const auto grid = *GridSpec::create({.finestSpacing = sub0math::f64::Length::fromSi(0.02),
                                         .cellsPerSide = 40,
                                         .registration = Registration::CellCentre,
                                         .anchor = {sub0math::f64::Length{}, sub0math::f64::Length::fromSi(-1.6)},
                                         .pages = {10, 4}});
    std::vector<double> xs(1024), ys(1024);
    for (std::size_t k = 0; k < xs.size(); ++k)
    {
        xs[k] = 8.0 * static_cast<double>((k * 389) % 1024) / 1024.0;
        ys[k] = -1.6 + 3.2 * static_cast<double>((k * 613) % 1024) / 1024.0;
    }
    runner.run("grid cellOf x1024, hand-written floor and range", [&] {
        long sum = 0;
        for (std::size_t k = 0; k < xs.size(); ++k)
        {
            const double ux = (xs[k] - 0.0) / 0.02, uy = (ys[k] + 1.6) / 0.02;
            if (!(ux >= 0.0 && uy >= 0.0 && ux <= 400.0 && uy <= 160.0)) continue;
            const int ix = std::min(static_cast<int>(std::floor(ux)), 399), iy = std::min(static_cast<int>(std::floor(uy)), 159);
            sum += (ix / 40) * 7 + (iy / 40) * 3 + ix % 40 + iy % 40;
        }
        ankerl::nanobench::doNotOptimizeAway(sum);
    });
    runner.run("grid cellOf x1024, GridSpec::cellOf", [&] {
        long sum = 0;
        for (std::size_t k = 0; k < xs.size(); ++k)
        {
            const TerrainPoint p{{sub0math::f64::Length::fromSi(xs[k]), sub0math::f64::Length::fromSi(ys[k]), sub0math::f64::Length{}}};
            const CellLookup c = grid.cellOf(p, 0);
            if (c.known) sum += c.page.i * 7 + c.page.j * 3 + c.i + c.j;
        }
        ankerl::nanobench::doNotOptimizeAway(sum);
    });
    const auto corridor = *CorridorBuilder::create(
        {.from = {sub0math::f64::Length::fromSi(4.0), sub0math::f64::Length::fromSi(0.0), sub0math::f64::Angle::fromSi(0.1)},
         .to = {sub0math::f64::Length::fromSi(4.1), sub0math::f64::Length::fromSi(0.02), sub0math::f64::Angle::fromSi(0.12)},
         .footprint = {sub0math::f64::Length::fromSi(0.38), sub0math::f64::Length::fromSi(0.05), sub0math::f64::Length::fromSi(0.1)}});
    std::array<std::uint64_t, 256> words;
    runner.run("corridor rasterise, 8 boxes, 2 cm cells", [&] {
        ankerl::nanobench::doNotOptimizeAway(corridor.rasterise(grid.cellGrid(0), words));
    });
    runner.run("corridor create, 8 boxes", [&] {
        ankerl::nanobench::doNotOptimizeAway(CorridorBuilder::create(
            {.from = {sub0math::f64::Length::fromSi(4.0), sub0math::f64::Length::fromSi(0.0), sub0math::f64::Angle::fromSi(0.1)},
             .to = {sub0math::f64::Length::fromSi(4.1), sub0math::f64::Length::fromSi(0.02), sub0math::f64::Angle::fromSi(0.12)},
             .footprint = {sub0math::f64::Length::fromSi(0.38), sub0math::f64::Length::fromSi(0.05), sub0math::f64::Length::fromSi(0.1)}}));
    });
}

// Travel buffer over 1024 (q, v) samples spanning both zones and the free zone, against plain doubles.
void benchTravelBuffer(sub0bench::BenchRunner& runner)
{
    constexpr double travel = 0.065, width = 0.008, stiffness = 1.0e6, damping = 3000.0;
    const auto buffer = *TravelBuffer::create({.travel = sub0math::metres(travel),
                                              .zoneWidth = sub0math::metres(width),
                                              .stiffness = sub0math::newtonsPerMetre(stiffness),
                                              .damping = sub0math::newtonSecondsPerMetre(damping)});
    std::vector<double> qs(1024), vs(1024);
    for (std::size_t k = 0; k < qs.size(); ++k)
    {
        qs[k] = -0.004 + 0.073 * static_cast<double>((k * 389) % 1024) / 1024.0;
        vs[k] = -2.0 + 4.0 * static_cast<double>((k * 613) % 1024) / 1024.0;
    }

    runner.run("travel buffer x1024, hand-written scalars", [&] {
        double sum = 0.0;
        const double invWidth = 1.0 / width;
        for (std::size_t k = 0; k < qs.size(); ++k)
        {
            const double q = qs[k], v = vs[k];
            const double top = std::max(0.0, width - q), bottom = std::max(0.0, q - (travel - width));
            const bool closingTop = v < 0.0;
            const double speed = closingTop ? -v : v;
            const double topForce = (stiffness * top + (closingTop ? damping * speed : 0.0)) * top * invWidth;
            const double bottomForce = (stiffness * bottom + (closingTop ? 0.0 : damping * speed)) * bottom * invWidth;
            const double depth = closingTop ? top : bottom;
            const double slope = (top > 0.0 ? 2.0 * stiffness * top + (closingTop ? damping * speed : 0.0) : 0.0) +
                                 (bottom > 0.0 ? 2.0 * stiffness * bottom + (closingTop ? 0.0 : damping * speed) : 0.0);
            sum += (topForce - bottomForce) + stiffness * (top * top * top + bottom * bottom * bottom) * invWidth / 3.0 +
                   damping * speed * speed * depth * invWidth + topForce + bottomForce - slope * invWidth - damping * depth * invWidth;
        }
        ankerl::nanobench::doNotOptimizeAway(sum);
    });
    runner.run("travel buffer x1024, TravelBuffer::respond", [&] {
        double sum = 0.0;
        for (std::size_t k = 0; k < qs.size(); ++k)
        {
            const TravelBufferResponse r = buffer.respond(sub0math::f64::Length::fromSi(qs[k]), sub0math::f64::Velocity::fromSi(vs[k]));
            sum += r.force.si() + r.potential.si() + r.dissipationPower.si() + r.topForce.si() + r.bottomForce.si() +
                   r.dForceDPosition.si() + r.dForceDVelocity.si();
        }
        ankerl::nanobench::doNotOptimizeAway(sum);
    });
}

} // namespace

int main(int argc, char** argv)
{
    sub0bench::BenchRunner runner(argc, argv, "Sub0Physics");
    runner.bench().minEpochIterations(20000);
    benchStep<Rk4<2>>(runner, "Rk4 step, contact oscillator");
    benchStep<Midpoint<2>>(runner, "Midpoint step, contact oscillator");

    std::vector<ProfilePoint> points;
    for (int i = 0; i <= 2047; ++i)
    {
        const double x = -10.0 + 20.0 * i / 2047.0;
        points.push_back({sub0math::metres(x), sub0math::metres(0.05 * std::sin(3.0 * x))});
    }
    const auto profile = *PolylineProfile::create(points);
    runner.bench().minEpochIterations(5000);
    runner.run("Circle clearance, 2048-point profile", [&] {
        std::array<Witness2, DefaultWitnessCapacity> witnesses;
        auto c = measureCircleClearance(profile, {{{sub0math::metres(0.3), sub0math::metres(0.5)}, sub0math::metres(0.36)}, sub0math::metres(1e-10)}, witnesses);
        ankerl::nanobench::doNotOptimizeAway(c);
        ankerl::nanobench::doNotOptimizeAway(witnesses);
    });
    benchClearance(runner, profile);
    runner.bench().minEpochIterations(20000);
    runner.run("Surface height, 2048-point profile", [&] {
        ankerl::nanobench::doNotOptimizeAway(profile.surfaceHeight(sub0math::metres(0.3)));
    });
    runner.run("Damped oscillator closed form", [&] {
        ankerl::nanobench::doNotOptimizeAway(evaluateDampedOscillator({sub0math::kilograms(45.0), sub0math::newtonsPerMetre(9000.0), sub0math::newtonSecondsPerMetre(500.0), sub0math::metres(0.03), sub0math::metresPerSecond(0.4)}, sub0math::seconds(0.731)));
    });
    benchPageHeight(runner);
    benchPageSize(runner);
    benchGrid(runner);
    benchTravelBuffer(runner);
    benchMultibody(runner);
    return runner.finish();
}
