#pragma once

namespace sub0physics
{

/** Response of a compression-only contact at one instant. */
struct ContactResponse
{
    double compression;     ///< Penetration clipped at zero, m.
    double force;           ///< Normal force, N; never negative (no adhesion).
    double storedEnergy;    ///< Elastic energy, J.
    double dissipationRate; ///< Power not returned by the spring, W.
};

/**
 * Spring and damper that can push but never pull: a penalty contact surrogate.
 *
 * When the damper would make the total negative the force is clipped to zero; the power the spring
 * then gives up is reported as dissipation, so that stored energy, dissipation and external work
 * still balance.
 */
struct UnilateralSpringDamper
{
    double stiffness; ///< Contact stiffness, N/m.
    double damping;   ///< Contact damping, N s/m.

    /**
     * Evaluates the contact.
     *
     * @param penetration Overlap, m; zero or negative means separated.
     * @param penetrationRate Rate of change of overlap, m/s.
     * @return The clipped response; all zero when separated.
     */
    [[nodiscard]] constexpr ContactResponse respond(double penetration, double penetrationRate) const noexcept
    {
        if (!(penetration > 0.0))
        {
            return ContactResponse{0.0, 0.0, 0.0, 0.0};
        }
        const double elastic = stiffness * penetration;
        const double total = elastic + damping * penetrationRate;
        const double force = total > 0.0 ? total : 0.0;
        return ContactResponse{
            penetration,
            force,
            0.5 * stiffness * penetration * penetration,
            (force - elastic) * penetrationRate,
        };
    }
};

} // namespace sub0physics
