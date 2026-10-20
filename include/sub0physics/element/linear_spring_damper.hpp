#pragma once

#include <cmath>
#include <expected>

#include <sub0math/units/quantities.hpp>

#include "sub0physics/element/law_response.hpp"
#include "sub0physics/fault.hpp"

namespace sub0physics
{

/**
 * The parameters of a LinearSpringDamper on a coordinate: `{.stiffness = ..., .damping = ...}`.
 *
 * @tparam Coordinate `f64::Length` (stiffness in N/m, damping in N s/m) or `f64::Angle` (N m/rad, N m s/rad).
 */
template<typename Coordinate>
struct LinearSpringDamperSpec
{
    using StiffnessType = decltype(GeneralizedForce<Coordinate>{} / Coordinate{});             ///< Force per unit coordinate.
    using DampingType = decltype(GeneralizedForce<Coordinate>{} / CoordinateRate<Coordinate>{}); ///< Force per unit rate.

    StiffnessType stiffness; ///< Spring rate; finite and non-negative.
    DampingType damping;     ///< Viscous damping; finite and non-negative.
};

/**
 * A linear spring in parallel with a linear viscous damper, relaxed at coordinate zero: the law of a
 * coordinate spring-damper (docs/MULTIBODY.md section 4), on a linear or a rotational coordinate.
 *
 * `force = -k q - c qd`, `potential = k q^2 / 2` and `dissipationPower = c qd^2`, so
 * `force qd = -dV/dt - dissipationPower`. The Jacobians are exact: `-k` and `-c`. A reference other than
 * zero belongs to the element that binds the law (multibody::CoordinateElementDescription).
 *
 * @tparam Coordinate `f64::Length` or `f64::Angle`.
 */
template<typename Coordinate>
class LinearSpringDamper
{
public:
    /** The parameters. */
    using Spec = LinearSpringDamperSpec<Coordinate>;

    /**
     * Validates and builds the law.
     *
     * @param spec The parameters.
     * @return The law, or an `InvalidInput` fault for a non-finite or negative value.
     */
    [[nodiscard]] static constexpr std::expected<LinearSpringDamper, Fault> create(const Spec& spec) noexcept
    {
        const double k = spec.stiffness.si();
        const double c = spec.damping.si();
        if (!std::isfinite(k) || !(k >= 0.0) || !std::isfinite(c) || !(c >= 0.0))
        {
            return invalidInput("linear spring-damper needs finite non-negative stiffness and damping");
        }
        return LinearSpringDamper{spec};
    }

    /** @return The parameters it was built with. */
    [[nodiscard]] constexpr const Spec& spec() const noexcept
    {
        return spec_;
    }

    /**
     * Evaluates the law.
     *
     * @param q Coordinate, relative to where the spring is relaxed (m or rad).
     * @param qd Coordinate rate (m/s or rad/s).
     * @return Generalized force on `q` (N or N m), the exact Jacobians, stored energy and dissipated power.
     */
    [[nodiscard]] constexpr LawResponse<Coordinate> respond(Coordinate q, CoordinateRate<Coordinate> qd) const noexcept
    {
        const double k = spec_.stiffness.si();
        const double c = spec_.damping.si();
        const double x = q.si();
        const double v = qd.si();
        LawResponse<Coordinate> r{};
        r.force = GeneralizedForce<Coordinate>::fromSi(-k * x - c * v);
        r.dForceDPosition = decltype(r.dForceDPosition)::fromSi(-k);
        r.dForceDVelocity = decltype(r.dForceDVelocity)::fromSi(-c);
        r.potential = sub0math::f64::Energy::fromSi(0.5 * k * x * x);
        r.dissipationPower = sub0math::f64::Power::fromSi(c * v * v);
        return r;
    }

private:
    constexpr explicit LinearSpringDamper(const Spec& spec) noexcept : spec_(spec)
    {
    }

    Spec spec_;
};

static_assert(Law<LinearSpringDamper<sub0math::f64::Length>, sub0math::f64::Length>);
static_assert(Law<LinearSpringDamper<sub0math::f64::Angle>, sub0math::f64::Angle>);

} // namespace sub0physics
