#include <cmath>
#include <doctest.h>
#include <sub0math/units.hpp>

#include "sub0physics/integrators.hpp"
#include "sub0physics/references.hpp"

using namespace sub0physics;
using namespace sub0math;

namespace
{

// Integrator state is plain SI; the model converts its typed parameters once, at its edge.
struct LinearOscillator
{
    double k;
    double c;
    double m;

    std::expected<void, Fault> evaluate(double, const StateVector<2>& x, StateVector<2>& dxdt) const
    {
        dxdt[0] = x[1];
        dxdt[1] = -(k * x[0] + c * x[1]) / m;
        return {};
    }
};

// Integrates with a fine RK4 step: an independent check of the closed form.
OscillatorState integrate(const DampedOscillator& p, double t)
{
    LinearOscillator system{p.stiffness.si(), p.damping.si(), p.mass.si()};
    Rk4<2> rk4;
    StateVector<2> x{p.initialDisplacement.si(), p.initialVelocity.si()};
    const int steps = 20000;
    const double h = t / steps;
    for (int i = 0; i < steps; ++i)
    {
        REQUIRE(rk4.step(system, sub0math::seconds(i * h), sub0math::seconds(h), x));
    }
    return {f64::Length::fromSi(x[0]), f64::Velocity::fromSi(x[1])};
}

DampedOscillator oscillator(double mass, double stiffness, double damping, double e0, double u0)
{
    return {kilograms(mass), newtonsPerMetre(stiffness), newtonSecondsPerMetre(damping), metres(e0),
            metresPerSecond(u0)};
}

void checkAgainstIntegration(double damping)
{
    const DampedOscillator p = oscillator(45.0, 9000.0, damping, 0.02, -0.3);
    for (const double t : {0.0, 1e-6, 0.05, 0.4, 1.5})
    {
        const OscillatorState exact = evaluateDampedOscillator(p, seconds(t));
        const OscillatorState numeric = t == 0.0 ? OscillatorState{p.initialDisplacement, p.initialVelocity} : integrate(p, t);
        CHECK(exact.displacement.si() == doctest::Approx(numeric.displacement.si()).epsilon(1e-7).scale(0.02));
        CHECK(exact.velocity.si() == doctest::Approx(numeric.velocity.si()).epsilon(1e-7).scale(0.3));
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
    const OscillatorState s = evaluateDampedOscillator(oscillator(30.0, 22000.0, 500.0, -0.01, 0.7), seconds(0.0));
    CHECK(s.displacement.si() == doctest::Approx(-0.01));
    CHECK(s.velocity.si() == doctest::Approx(0.7));
}

TEST_CASE("Strongly over-damped decay stays finite for long times")
{
    const OscillatorState s = evaluateDampedOscillator(oscillator(1.0, 1000.0, 8000.0, 0.1, 0.0), seconds(1000.0));
    CHECK(isFinite(s.displacement));
    CHECK(isFinite(s.velocity));
}
