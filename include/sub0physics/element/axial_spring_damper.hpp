#pragma once

#include <sub0math/units/quantities.hpp>

namespace sub0physics
{

/** Forces and power of an axial element at one instant. Signs: positive force spreads the element. */
struct AxialResponse
{
    sub0math::f64::Force springForce;     ///< Preload plus linear spring force.
    sub0math::f64::Force damperForce;     ///< Viscous force; same sign as the compression rate.
    sub0math::f64::Force force;           ///< Total spreading force.
    sub0math::f64::Energy storedEnergy;   ///< Elastic energy relative to zero compression.
    sub0math::f64::Power dissipationRate; ///< Power removed by the damper; never negative.
};

/**
 * Preloaded linear spring in parallel with a direction-dependent viscous damper.
 *
 * Compression `q` grows as the element shortens; `u = dq/dt`. Compression damping applies while
 * `u >= 0` and rebound damping while `u < 0`. The element is defined for any `q`; travel limits belong
 * to the caller, which knows its stroke. Initialise with designated initialisers so the two damping
 * coefficients cannot be swapped: `{.preload = ..., .stiffness = ..., .compressionDamping = ...,
 * .reboundDamping = ...}`.
 */
struct AxialSpringDamper
{
    sub0math::f64::Force preload;              ///< Spring force at zero compression.
    sub0math::f64::Stiffness stiffness;        ///< Spring rate.
    sub0math::f64::Damping compressionDamping; ///< Damping coefficient while compressing.
    sub0math::f64::Damping reboundDamping;     ///< Damping coefficient while extending.

    /**
     * Evaluates the element.
     *
     * @param compression Compression `q`.
     * @param compressionRate Compression rate `u`.
     * @return Forces, stored energy and dissipated power at this state.
     */
    [[nodiscard]] constexpr AxialResponse respond(sub0math::f64::Length compression,
                                                  sub0math::f64::Velocity compressionRate) const noexcept
    {
        const sub0math::f64::Damping damping =
            compressionRate >= sub0math::f64::Velocity{} ? compressionDamping : reboundDamping;
        const sub0math::f64::Force spring = preload + stiffness * compression;
        const sub0math::f64::Force damper = damping * compressionRate;
        return AxialResponse{
            spring,
            damper,
            spring + damper,
            preload * compression + 0.5 * stiffness * compression * compression,
            damping * compressionRate * compressionRate,
        };
    }
};

} // namespace sub0physics
