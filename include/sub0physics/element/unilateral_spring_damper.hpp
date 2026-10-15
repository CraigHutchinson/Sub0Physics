#pragma once

#include <sub0math/units/quantities.hpp>

namespace sub0physics
{

/** Response of a compression-only contact at one instant. */
struct ContactResponse
{
    sub0math::f64::Length compression;    ///< Penetration clipped at zero.
    sub0math::f64::Force force;           ///< Normal force; never negative (no adhesion).
    sub0math::f64::Energy storedEnergy;   ///< Elastic energy.
    sub0math::f64::Power dissipationRate; ///< Power not returned by the spring.
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
    sub0math::f64::Stiffness stiffness; ///< Contact stiffness.
    sub0math::f64::Damping damping;     ///< Contact damping.

    /**
     * Evaluates the contact.
     *
     * @param penetration Overlap; zero or negative means separated.
     * @param penetrationRate Rate of change of overlap.
     * @return The clipped response; all zero when separated.
     */
    [[nodiscard]] constexpr ContactResponse respond(sub0math::f64::Length penetration,
                                                    sub0math::f64::Velocity penetrationRate) const noexcept
    {
        if (!(penetration > sub0math::f64::Length{}))
        {
            return ContactResponse{};
        }
        const sub0math::f64::Force elastic = stiffness * penetration;
        const sub0math::f64::Force total = elastic + damping * penetrationRate;
        const sub0math::f64::Force force = total > sub0math::f64::Force{} ? total : sub0math::f64::Force{};
        return ContactResponse{
            penetration,
            force,
            0.5 * stiffness * penetration * penetration,
            (force - elastic) * penetrationRate,
        };
    }
};

} // namespace sub0physics
