#define ANKERL_NANOBENCH_IMPLEMENT
#include <algorithm>
#include <bench_runner.hpp>
#include <cmath>
#include <cstdint>
#include <expected>
#include <sub0math/units.hpp>
#include <vector>

#include "sub0physics/sub0physics.hpp"

using namespace sub0physics;

namespace
{

// A loaded wheel on a clipped contact: the shape of BikeLab's first station.
struct ContactOscillator
{
    UnilateralSpringDamper contact{sub0math::newtonsPerMetre(90000.0), sub0math::newtonSecondsPerMetre(400.0)};
    double mass = 45.0;
    double radius = 0.36;

    std::expected<void, Fault> evaluate(double, const StateVector<2>& x, StateVector<2>& dxdt) const
    {
        const ContactResponse r = contact.respond(sub0math::metres(radius - x[0]), sub0math::metresPerSecond(-x[1]));
        dxdt[0] = x[1];
        dxdt[1] = r.force.si() / mass - StandardGravity;
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
        if (!integrator.stepTick(system, tick++, 0.0005, x))
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
    const GridSpec grid{sub0math::f64::Length::fromSi(0.01), 4};
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

// Travel buffer over 1024 (q, v) samples spanning both zones and the free zone, against plain doubles.
void benchTravelBuffer(sub0bench::BenchRunner& runner)
{
    constexpr double travel = 0.065, width = 0.008, stiffness = 1.0e6, damping = 3000.0;
    const auto buffer = *TravelBuffer::create(sub0math::metres(travel), sub0math::metres(width),
                                              sub0math::newtonsPerMetre(stiffness), sub0math::newtonSecondsPerMetre(damping));
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
            sum += r.force.si() + r.potential.si() + r.lossPower.si() + r.topForce.si() + r.bottomForce.si() +
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
        auto c = measureCircleClearance(profile, {{{sub0math::metres(0.3), sub0math::metres(0.5)}, sub0math::metres(0.36)}, sub0math::metres(1e-10)});
        ankerl::nanobench::doNotOptimizeAway(c);
    });
    runner.bench().minEpochIterations(20000);
    runner.run("Surface height, 2048-point profile", [&] {
        ankerl::nanobench::doNotOptimizeAway(profile.surfaceHeight(sub0math::metres(0.3)));
    });
    runner.run("Damped oscillator closed form", [&] {
        ankerl::nanobench::doNotOptimizeAway(evaluateDampedOscillator({sub0math::kilograms(45.0), sub0math::newtonsPerMetre(9000.0), sub0math::newtonSecondsPerMetre(500.0), sub0math::metres(0.03), sub0math::metresPerSecond(0.4)}, sub0math::seconds(0.731)));
    });
    benchPageHeight(runner);
    benchTravelBuffer(runner);
    return runner.finish();
}
