#define ANKERL_NANOBENCH_IMPLEMENT
#include <cmath>
#include <expected>
#include <nanobench.h>
#include <numbers>
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
void benchStep(ankerl::nanobench::Bench& bench, const char* name)
{
    ContactOscillator system;
    Integrator integrator;
    StateVector<2> x{0.36, 0.0};
    double t = 0.0;
    const double h = 0.0005;
    bench.run(name, [&] {
        if (!integrator.step(system, t, h, x))
        {
            x = {0.36, 0.0};
        }
        t += h;
        ankerl::nanobench::doNotOptimizeAway(x);
    });
}

} // namespace

int main()
{
    ankerl::nanobench::Bench bench;
    bench.title("Sub0Physics").relative(true).warmup(1000).minEpochIterations(20000);
    benchStep<Rk4<2>>(bench, "Rk4 step, contact oscillator");
    benchStep<Midpoint<2>>(bench, "Midpoint step, contact oscillator");

    std::vector<Point2> points;
    for (int i = 0; i <= 2047; ++i)
    {
        const double x = -10.0 + 20.0 * i / 2047.0;
        points.push_back({x, 0.05 * std::sin(3.0 * x)});
    }
    const auto profile = *PolylineProfile::create(points);
    bench.relative(false).minEpochIterations(200);
    bench.run("Circle clearance, 2048-point profile", [&] {
        auto c = measureCircleClearance(profile, {{0.3, 0.5}, 0.36, 1e-10});
        ankerl::nanobench::doNotOptimizeAway(c);
    });
    bench.minEpochIterations(20000);
    bench.run("Surface height, 2048-point profile", [&] {
        ankerl::nanobench::doNotOptimizeAway(profile.surfaceHeight(0.3));
    });
}
