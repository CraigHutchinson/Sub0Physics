#pragma once

namespace sub0physics
{

/** Forces and power of an axial element at one instant. Signs: positive force spreads the element. */
struct AxialResponse
{
    double springForce;     ///< Preload plus linear spring force, N.
    double damperForce;     ///< Viscous force, N; same sign as the compression rate.
    double force;           ///< Total spreading force, N.
    double storedEnergy;    ///< Elastic energy relative to zero compression, J.
    double dissipationRate; ///< Power removed by the damper, W; never negative.
};

/**
 * Preloaded linear spring in parallel with a direction-dependent viscous damper.
 *
 * Compression `q` grows as the element shortens; `u = dq/dt`. Compression damping applies while
 * `u >= 0` and rebound damping while `u < 0`. The element is defined for any `q`; travel limits belong
 * to the caller, which knows its stroke.
 */
struct AxialSpringDamper
{
    double preload;            ///< Spring force at zero compression, N.
    double stiffness;          ///< Spring rate, N/m.
    double compressionDamping; ///< Damping coefficient while compressing, N s/m.
    double reboundDamping;     ///< Damping coefficient while extending, N s/m.

    /**
     * Evaluates the element.
     *
     * @param compression Compression `q`, m.
     * @param compressionRate Compression rate `u`, m/s.
     * @return Forces, stored energy and dissipated power at this state.
     */
    [[nodiscard]] constexpr AxialResponse respond(double compression, double compressionRate) const noexcept
    {
        const double damping = compressionRate >= 0.0 ? compressionDamping : reboundDamping;
        const double spring = preload + stiffness * compression;
        const double damper = damping * compressionRate;
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
