#pragma once

#include <cmath>
#include <expected>

#include <sub0math/units/quantities.hpp>

#include "sub0physics/element/law_response.hpp"
#include "sub0physics/fault.hpp"

namespace sub0physics
{

/**
 * Response of an AxialSpringDamper: the common law response plus the spring and damper parts.
 *
 * `force` is the generalized force on the compression coordinate, so it is negative while the element
 * pushes its ends apart: `force = -(springForce + damperForce)`.
 */
struct AxialResponse : LawResponse<sub0math::f64::Length>
{
    sub0math::f64::Force springForce; ///< Spreading force of the spring: preload plus linear spring force.
    sub0math::f64::Force damperForce; ///< Spreading force of the damper; same sign as the compression rate.
};

/**
 * The parameters of an AxialSpringDamper. Initialise with designated initialisers so the two damping
 * coefficients cannot be swapped: `{.preload = ..., .stiffness = ..., .compressionDamping = ...,
 * .reboundDamping = ...}`.
 */
struct AxialSpringDamperSpec
{
    sub0math::f64::Force preload;              ///< Spring force at zero compression [N]; finite.
    sub0math::f64::Stiffness stiffness;        ///< Spring rate [N/m]; finite and non-negative.
    sub0math::f64::Damping compressionDamping; ///< Damping while compressing [N s/m]; finite and non-negative.
    sub0math::f64::Damping reboundDamping;     ///< Damping while extending [N s/m]; finite and non-negative.
};

/**
 * Preloaded linear spring in parallel with a direction-dependent viscous damper: a law on a compression
 * coordinate.
 *
 * Compression `q` grows as the element shortens; `u = dq/dt`. Compression damping applies while
 * `u >= 0` and rebound damping while `u < 0`. The element is defined for any `q`; travel limits belong
 * to the caller, which knows its stroke.
 *
 * `force` is the generalized force on `q` (negative while the element spreads), `potential` is
 * `preload q + k q^2 / 2` and `dissipationPower` the damper power, so
 * `force u = -dV/dt - dissipationPower`. The Jacobians are exact: `dForceDPosition = -k`, and
 * `dForceDVelocity` is minus the damping of the branch the rate selects (compression at `u = 0`).
 */
class AxialSpringDamper
{
public:
    /**
     * Validates and builds the law.
     *
     * @param spec The parameters; see AxialSpringDamperSpec.
     * @return The law, or an `InvalidInput` fault for a non-finite value or a negative rate or damping.
     */
    [[nodiscard]] static constexpr std::expected<AxialSpringDamper, Fault>
    create(const AxialSpringDamperSpec& spec) noexcept
    {
        const double p = spec.preload.si();
        const double k = spec.stiffness.si();
        const double cc = spec.compressionDamping.si();
        const double cr = spec.reboundDamping.si();
        if (!std::isfinite(p) || !std::isfinite(k) || !(k >= 0.0))
        {
            return invalidInput("axial spring-damper needs a finite preload and a finite non-negative stiffness");
        }
        if (!std::isfinite(cc) || !(cc >= 0.0) || !std::isfinite(cr) || !(cr >= 0.0))
        {
            return invalidInput("axial spring-damper needs finite non-negative damping coefficients");
        }
        return AxialSpringDamper{spec};
    }

    /** @return The parameters it was built with. */
    [[nodiscard]] constexpr const AxialSpringDamperSpec& spec() const noexcept
    {
        return spec_;
    }

    /**
     * Evaluates the law.
     *
     * @param compression Compression `q` [m].
     * @param compressionRate Compression rate `u` [m/s].
     * @return Generalized force, Jacobians, stored energy, dissipated power and the two spreading parts.
     */
    [[nodiscard]] constexpr AxialResponse respond(sub0math::f64::Length compression,
                                                  sub0math::f64::Velocity compressionRate) const noexcept
    {
        const sub0math::f64::Damping damping =
            compressionRate >= sub0math::f64::Velocity{} ? spec_.compressionDamping : spec_.reboundDamping;
        const sub0math::f64::Force spring = spec_.preload + spec_.stiffness * compression;
        const sub0math::f64::Force damper = damping * compressionRate;
        AxialResponse r{};
        r.force = -(spring + damper);
        r.dForceDPosition = -spec_.stiffness;
        r.dForceDVelocity = -damping;
        r.potential = spec_.preload * compression + 0.5 * spec_.stiffness * compression * compression;
        r.dissipationPower = damping * compressionRate * compressionRate;
        r.springForce = spring;
        r.damperForce = damper;
        return r;
    }

private:
    constexpr explicit AxialSpringDamper(const AxialSpringDamperSpec& spec) noexcept : spec_(spec)
    {
    }

    AxialSpringDamperSpec spec_;
};

static_assert(Law<AxialSpringDamper, sub0math::f64::Length>);

} // namespace sub0physics
