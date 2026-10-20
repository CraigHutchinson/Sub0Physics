#pragma once

#include <cmath>
#include <expected>

#include <sub0math/config.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/element/law_response.hpp"
#include "sub0physics/fault.hpp"

namespace sub0physics
{

/**
 * Response of a UnilateralSpringDamper: the common law response plus the contact's own view.
 *
 * `force` is the generalized force on the penetration coordinate, so it is never positive:
 * `force = -normalForce`.
 */
struct ContactResponse : LawResponse<sub0math::f64::Length>
{
    sub0math::f64::Length compression; ///< Penetration clipped at zero.
    sub0math::f64::Force normalForce;  ///< Contact force pushing the bodies apart; never negative (no adhesion).
};

/** The parameters of a UnilateralSpringDamper: `{.stiffness = ..., .damping = ...}`. */
struct UnilateralSpringDamperSpec
{
    sub0math::f64::Stiffness stiffness; ///< Contact stiffness [N/m]; finite and non-negative.
    sub0math::f64::Damping damping;     ///< Contact damping [N s/m]; finite and non-negative.
};

/**
 * Spring and damper that can push but never pull: a penalty contact surrogate, as a law on the
 * penetration coordinate.
 *
 * Penetration `q` is the overlap (zero or negative means separated) and `u = dq/dt`. While `q > 0` the
 * normal force is `max(0, k q + c u)`. When the damper would make it negative the force is clipped to
 * zero; the power the spring then gives up is reported as dissipation, so that
 * `force u = -dV/dt - dissipationPower` still holds with `V = k q^2 / 2`.
 *
 * The Jacobians are exact piecewise: `-k` and `-c` while the contact pushes, zero when separated
 * (`q <= 0`) or clipped (`k q + c u <= 0`).
 */
class UnilateralSpringDamper
{
public:
    /**
     * Validates and builds the law.
     *
     * @param spec The parameters; see UnilateralSpringDamperSpec.
     * @return The law, or an `InvalidInput` fault for a non-finite or negative value.
     */
    [[nodiscard]] static constexpr std::expected<UnilateralSpringDamper, Fault>
    create(const UnilateralSpringDamperSpec& spec) noexcept
    {
        const double k = spec.stiffness.si();
        const double c = spec.damping.si();
        if (!std::isfinite(k) || !(k >= 0.0) || !std::isfinite(c) || !(c >= 0.0))
        {
            return invalidInput("unilateral spring-damper needs finite non-negative stiffness and damping");
        }
        return UnilateralSpringDamper{spec};
    }

    /** @return The parameters it was built with. */
    [[nodiscard]] constexpr const UnilateralSpringDamperSpec& spec() const noexcept
    {
        return spec_;
    }

    /**
     * Evaluates the contact.
     *
     * @param penetration Overlap `q` [m]; zero or negative means separated.
     * @param penetrationRate Rate of change of overlap `u` [m/s].
     * @return The clipped response; all zero when separated.
     */
    [[nodiscard]] SUB0MATH_FORCE_INLINE constexpr ContactResponse respond(sub0math::f64::Length penetration,
                                                    sub0math::f64::Velocity penetrationRate) const noexcept
    {
        if (!(penetration > sub0math::f64::Length{}))
        {
            return ContactResponse{};
        }
        const sub0math::f64::Force elastic = spec_.stiffness * penetration;
        const sub0math::f64::Force total = elastic + spec_.damping * penetrationRate;
        const bool pushing = total > sub0math::f64::Force{};
        const sub0math::f64::Force normal = pushing ? total : sub0math::f64::Force{};
        return ContactResponse{
            {
                -normal,
                pushing ? -spec_.stiffness : sub0math::f64::Stiffness{},
                pushing ? -spec_.damping : sub0math::f64::Damping{},
                0.5 * spec_.stiffness * penetration * penetration,
                (normal - elastic) * penetrationRate,
            },
            penetration,
            normal,
        };
    }

private:
    constexpr explicit UnilateralSpringDamper(const UnilateralSpringDamperSpec& spec) noexcept : spec_(spec)
    {
    }

    UnilateralSpringDamperSpec spec_;
};

static_assert(Law<UnilateralSpringDamper, sub0math::f64::Length>);

} // namespace sub0physics
