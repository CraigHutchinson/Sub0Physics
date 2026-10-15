#define ANKERL_NANOBENCH_IMPLEMENT
#include <bench_runner.hpp>
#include <cmath>
#include <cstdint>
#include <expected>
#include <vector>

#include "sub0physics/sub0physics.hpp"

using namespace sub0physics;

namespace
{

// A loaded wheel on a clipped contact: the shape of BikeLab's first station.
struct ContactOscillator
{
    UnilateralSpringDamper contact{90000.0, 400.0};
    double mass = 45.0;
    double radius = 0.36;

    std::expected<void, Fault> evaluate(double, const StateVector<2>& x, StateVector<2>& dxdt) const
    {
        const ContactResponse r = contact.respond(radius - x[0], -x[1]);
        dxdt[0] = x[1];
        dxdt[1] = r.force / mass - StandardGravity;
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

} // namespace

int main(int argc, char** argv)
{
    sub0bench::BenchRunner runner(argc, argv, "Sub0Physics");
    runner.bench().minEpochIterations(20000);
    benchStep<Rk4<2>>(runner, "Rk4 step, contact oscillator");
    benchStep<Midpoint<2>>(runner, "Midpoint step, contact oscillator");

    std::vector<Point2> points;
    for (int i = 0; i <= 2047; ++i)
    {
        const double x = -10.0 + 20.0 * i / 2047.0;
        points.push_back({x, 0.05 * std::sin(3.0 * x)});
    }
    const auto profile = *PolylineProfile::create(points);
    runner.bench().minEpochIterations(200);
    runner.run("Circle clearance, 2048-point profile", [&] {
        auto c = measureCircleClearance(profile, {{0.3, 0.5}, 0.36, 1e-10});
        ankerl::nanobench::doNotOptimizeAway(c);
    });
    runner.bench().minEpochIterations(20000);
    runner.run("Surface height, 2048-point profile", [&] {
        ankerl::nanobench::doNotOptimizeAway(profile.surfaceHeight(0.3));
    });
    runner.run("Damped oscillator closed form", [&] {
        ankerl::nanobench::doNotOptimizeAway(evaluateDampedOscillator({45.0, 9000.0, 500.0, 0.03, 0.4}, 0.731));
    });
    return runner.finish();
}
