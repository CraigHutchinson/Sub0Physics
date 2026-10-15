#pragma once

namespace sub0physics
{

/** Parameters of m e'' + c e' + k e = 0 with its initial conditions. */
struct DampedOscillator
{
    double mass;                ///< m, kg; positive.
    double stiffness;           ///< k, N/m; positive.
    double damping;             ///< c, N s/m; not negative.
    double initialDisplacement; ///< e(0), m.
    double initialVelocity;     ///< e'(0), m/s.
};

/** Displacement and velocity of an oscillator at one time. */
struct OscillatorState
{
    double displacement; ///< e(t), m.
    double velocity;     ///< e'(t), m/s.
};

/**
 * Evaluates the exact solution of a damped linear oscillator.
 *
 * Covers under-, critically and over-damped cases. Small-argument series avoid the 0/0 forms near
 * critical damping, and the strongly over-damped branch is written in its two decay rates to avoid
 * cancellation and overflow. Independent of any integrator, so it can judge one.
 *
 * @param oscillator Parameters; the caller guarantees positive mass and stiffness, non-negative
 *        damping and finite values.
 * @param t Time since the initial state, s; not negative.
 * @return The state at `t`.
 */
[[nodiscard]] OscillatorState evaluateDampedOscillator(const DampedOscillator& oscillator, double t) noexcept;

} // namespace sub0physics
