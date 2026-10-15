#pragma once

#include <sub0math/units/quantities.hpp>

namespace sub0physics
{

/** Parameters of m e'' + c e' + k e = 0 for a translational displacement e, with initial conditions. */
struct DampedOscillator
{
    sub0math::f64::Mass mass;                ///< m; positive.
    sub0math::f64::Stiffness stiffness;      ///< k; positive.
    sub0math::f64::Damping damping;          ///< c; not negative.
    sub0math::f64::Length initialDisplacement; ///< e(0).
    sub0math::f64::Velocity initialVelocity;   ///< e'(0).
};

/** Displacement and velocity of an oscillator at one time. */
struct OscillatorState
{
    sub0math::f64::Length displacement; ///< e(t).
    sub0math::f64::Velocity velocity;   ///< e'(t).
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
 * @param t Time since the initial state; not negative.
 * @return The state at `t`.
 */
[[nodiscard]] OscillatorState evaluateDampedOscillator(const DampedOscillator& oscillator, sub0math::f64::Time t) noexcept;

} // namespace sub0physics
