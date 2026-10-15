#include <cmath>
#include <doctest.h>

#include "sub0physics/integrators.hpp"
#include "sub0physics/references.hpp"

using namespace sub0physics;

namespace
{

struct LinearOscillator
{
    DampedOscillator p;

    std::expected<void, Fault> evaluate(double, const StateVector<2>& x, StateVector<2>& dxdt) const
    {
        dxdt[0] = x[1];
        dxdt[1] = -(p.stiffness * x[0] + p.damping * x[1]) / p.mass;
        return {};
    }
};

// Integrates with a fine RK4 step: an independent check of the closed form.
OscillatorState integrate(const DampedOscillator& p, double t)
{
    LinearOscillator system{p};
    Rk4<2> rk4;
    StateVector<2> x{p.initialDisplacement, p.initialVelocity};
    const int steps = 20000;
    const double h = t / steps;
    for (int i = 0; i < steps; ++i)
    {
        REQUIRE(rk4.step(system, i * h, h, x));
    }
    return {x[0], x[1]};
}

void checkAgainstIntegration(double damping)
{
    const DampedOscillator p{45.0, 9000.0, damping, 0.02, -0.3};
    for (const double t : {0.0, 1e-6, 0.05, 0.4, 1.5})
    {
        const OscillatorState exact = evaluateDampedOscillator(p, t);
        const OscillatorState numeric = t == 0.0 ? OscillatorState{p.initialDisplacement, p.initialVelocity} : integrate(p, t);
        CHECK(exact.displacement == doctest::Approx(numeric.displacement).epsilon(1e-7).scale(0.02));
        CHECK(exact.velocity == doctest::Approx(numeric.velocity).epsilon(1e-7).scale(0.3));
    }
}

} // namespace

TEST_CASE("Closed form matches integration in every damping regime")
{
    const double critical = 2.0 * std::sqrt(9000.0 * 45.0);
    SUBCASE("undamped") { checkAgainstIntegration(0.0); }
    SUBCASE("under-damped") { checkAgainstIntegration(350.0); }
    SUBCASE("critical") { checkAgainstIntegration(critical); }
    SUBCASE("just over-damped") { checkAgainstIntegration(critical * (1.0 + 1e-9)); }
    SUBCASE("over-damped") { checkAgainstIntegration(critical * 1.5); }
    SUBCASE("strongly over-damped") { checkAgainstIntegration(critical * 20.0); }
}

TEST_CASE("Closed form starts at the initial conditions")
{
    const OscillatorState s = evaluateDampedOscillator({30.0, 22000.0, 500.0, -0.01, 0.7}, 0.0);
    CHECK(s.displacement == doctest::Approx(-0.01));
    CHECK(s.velocity == doctest::Approx(0.7));
}

TEST_CASE("Strongly over-damped decay stays finite for long times")
{
    const OscillatorState s = evaluateDampedOscillator({1.0, 1000.0, 8000.0, 0.1, 0.0}, 1000.0);
    CHECK(std::isfinite(s.displacement));
    CHECK(std::isfinite(s.velocity));
}
