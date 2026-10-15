#include <cmath>

#include "sub0physics/reference/damped_oscillator.hpp"

namespace sub0physics
{

namespace
{

OscillatorState state(double displacement, double velocity) noexcept
{
    return {sub0math::f64::Length::fromSi(displacement), sub0math::f64::Velocity::fromSi(velocity)};
}

} // namespace

OscillatorState evaluateDampedOscillator(const DampedOscillator& oscillator, sub0math::f64::Time time) noexcept
{
    // Closed form in SI numbers; the interface carries the units.
    const double t = time.si();
    const double e0 = oscillator.initialDisplacement.si();
    const double u0 = oscillator.initialVelocity.si();
    const double alpha = oscillator.damping.si() / (2.0 * oscillator.mass.si());
    const double w2 = oscillator.stiffness.si() / oscillator.mass.si();
    const double d = alpha * alpha - w2;
    const double b = u0 + alpha * e0;

    if (d <= 0.0)
    {
        const double w = std::sqrt(-d);
        const double x = w * t;
        const double decay = std::exp(-alpha * t);
        // sin(wt)/w -> t as w -> 0 (critical damping).
        const double s = std::abs(x) < 1e-4 ? t * (1.0 - x * x / 6.0 + x * x * x * x / 120.0) : std::sin(x) / w;
        const double c = std::cos(x);
        return state(decay * (e0 * c + b * s), decay * (u0 * c - (alpha * u0 + w2 * e0) * s));
    }

    const double beta = std::sqrt(d);
    const double x = beta * t;
    if (x < 1.0)
    {
        const double decay = std::exp(-alpha * t);
        const double c = std::cosh(x);
        const double s = x < 1e-4 ? t * (1.0 + x * x / 6.0 + x * x * x * x / 120.0) : std::sinh(x) / beta;
        return state(decay * (e0 * c + b * s), decay * (u0 * c - (alpha * u0 + w2 * e0) * s));
    }

    // -alpha + beta cancels badly and cosh(beta t) overflows; use the two decay rates directly.
    const double slow = -w2 / (alpha + beta);
    const double fast = -alpha - beta;
    const double a = (u0 - fast * e0) / (2.0 * beta);
    const double bFast = e0 - a;
    const double es = std::exp(slow * t);
    const double ef = std::exp(fast * t);
    return state(a * es + bFast * ef, a * slow * es + bFast * fast * ef);
}

} // namespace sub0physics
