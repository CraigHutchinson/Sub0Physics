#pragma once

#include <concepts>
#include <type_traits>

#include <sub0math/units/quantities.hpp>

namespace sub0physics
{

/**
 * The quantities that belong to a law's coordinate.
 *
 * Specialised for `f64::Length` (a linear coordinate: rate `Velocity`, generalized force `Force`) and
 * `f64::Angle` (a rotational coordinate: rate `AngularVelocity`, generalized force `Torque`).
 *
 * @tparam Coordinate The coordinate quantity.
 */
template<typename Coordinate>
struct CoordinateTraits;

/** A linear coordinate, in metres. */
template<>
struct CoordinateTraits<sub0math::f64::Length>
{
    using Rate = sub0math::f64::Velocity; ///< Coordinate rate.
    using Force = sub0math::f64::Force;   ///< Generalized force: power per unit rate.
};

/** A rotational coordinate, in radians. */
template<>
struct CoordinateTraits<sub0math::f64::Angle>
{
    using Rate = sub0math::f64::AngularVelocity; ///< Coordinate rate.
    using Force = sub0math::f64::Torque;          ///< Generalized force: power per unit rate.
};

/**
 * The rate of a coordinate: `Velocity` for `Length`, `AngularVelocity` for `Angle`.
 *
 * @tparam Coordinate The coordinate quantity.
 */
template<typename Coordinate>
using CoordinateRate = typename CoordinateTraits<Coordinate>::Rate;

/**
 * The generalized force on a coordinate: `Force` for `Length`, `Torque` for `Angle`.
 *
 * @tparam Coordinate The coordinate quantity.
 */
template<typename Coordinate>
using GeneralizedForce = typename CoordinateTraits<Coordinate>::Force;

/**
 * What a law answers for one coordinate state: its generalized force, the exact partial derivatives of
 * that force, its stored energy and its dissipated power (docs/MULTIBODY.md, "M1b contracts").
 *
 * **Sign.** `force` is the generalized force on the law's own coordinate `q`: positive pushes toward
 * larger `q`. For every law, `force * qd = -dV/dt - dissipationPower`, where `V` is `potential`. A spring
 * on a compression coordinate therefore reports a negative force when compressed.
 *
 * A law with extra diagnostic outputs derives its response from this one, so the common fields keep
 * their names and an element can read any law's response through this type.
 *
 * @tparam Coordinate `f64::Length` (linear) or `f64::Angle` (rotational).
 */
template<typename Coordinate>
struct LawResponse
{
    using ForceType = GeneralizedForce<Coordinate>; ///< `Force` or `Torque`.
    using RateType = CoordinateRate<Coordinate>;    ///< `Velocity` or `AngularVelocity`.

    ForceType force;                                         ///< Generalized force on `q`; positive toward larger `q`.
    decltype(ForceType{} / Coordinate{}) dForceDPosition;    ///< Exact partial derivative of `force` with respect to `q`.
    decltype(ForceType{} / RateType{}) dForceDVelocity;      ///< Exact partial derivative of `force` with respect to `qd`.
    sub0math::f64::Energy potential;                         ///< Stored energy relative to the law's reference.
    sub0math::f64::Power dissipationPower;                   ///< Power removed; never negative.
};

static_assert(std::same_as<decltype(LawResponse<sub0math::f64::Length>{}.dForceDPosition), sub0math::f64::Stiffness>);
static_assert(std::same_as<decltype(LawResponse<sub0math::f64::Length>{}.dForceDVelocity), sub0math::f64::Damping>);

/**
 * A scalar constitutive law of one coordinate and its rate: a value type answering `respond(q, qd)` with
 * a LawResponse (or a type derived from it that adds diagnostics). Laws are built by a validating
 * `static create(const Spec&) -> std::expected<Law, Fault>`.
 *
 * @tparam L The law type.
 * @tparam Coordinate The coordinate quantity.
 */
template<typename L, typename Coordinate>
concept Law = requires(const L law, Coordinate q, CoordinateRate<Coordinate> qd) {
    { law.respond(q, qd) } noexcept;
    requires std::derived_from<std::remove_cvref_t<decltype(law.respond(q, qd))>, LawResponse<Coordinate>>;
};

} // namespace sub0physics
