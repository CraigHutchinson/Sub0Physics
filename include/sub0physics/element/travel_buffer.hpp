#pragma once

#include <cmath>
#include <expected>
#include <sub0math/config.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/fault.hpp"

namespace sub0physics
{

/**
 * Forces, energy, power and Jacobians of a travel buffer at one instant.
 *
 * Sign convention: `force` is the generalized force on the compression coordinate `q`. The top buffer
 * pushes toward larger `q` (positive), the bottom buffer toward smaller `q` (negative).
 */
struct TravelBufferResponse
{
    sub0math::f64::Force force;               ///< Generalized force on `q`: `topForce - bottomForce`.
    sub0math::f64::Energy potential;          ///< Elastic energy of the two buffers; never negative.
    sub0math::f64::Power lossPower;           ///< Power removed by the damping; never negative.
    sub0math::f64::Force topForce;            ///< Top buffer (near `q = 0`) force magnitude; never negative.
    sub0math::f64::Force bottomForce;         ///< Bottom buffer (near `q = travel`) force magnitude; never negative.
    sub0math::f64::Stiffness dForceDPosition; ///< Exact partial derivative of `force` with respect to `q`.
    sub0math::f64::Damping dForceDVelocity;   ///< Exact partial derivative of `force` with respect to `u`.
};

/**
 * Progressive compliant end-stops inside a finite stroke.
 *
 * Compression `q` runs from 0 (fully extended) to `travel`. A top buffer of depth `width` engages for
 * `q < width` and a bottom buffer of the same depth for `q > travel - width`; the zones never
 * overlap and the free zone between them exerts nothing. With penetration depth `d` into a buffer the
 * spring force is `k d^2 / width` (progressive, not constant tangent stiffness) and the damper adds
 * `c s d / width`, where `s` is the speed of entering the buffer. Damping acts only while entering;
 * leaving a buffer is undamped.
 *
 * `potential` is the exact integral of the spring part, so `-d(potential)/dq` equals the spring part
 * of `force`. `lossPower` is the damper power, so `force * u = -d(potential)/dt - lossPower`.
 *
 * Jacobians are exact piecewise partial derivatives of `force`. At the kinks the derivative of the
 * branch the input selects is used, consistently:
 *  - the top buffer is active for `q < width` and the bottom buffer for `q > travel - width`
 *    (strict), so a buffer contributes nothing to either Jacobian entry at its edge;
 *  - `u < 0` enters the top buffer and `u >= 0` the bottom one, so at `u = 0` the bottom branch is
 *    the active one.
 * Only the damping term has a one-sided slope at a `q` kink; the spring term is continuous in slope.
 *
 * Non-finite inputs give non-finite outputs; the caller owns that check, as for the other elements.
 */
class TravelBuffer
{
public:
    /**
     * Validates and builds a buffer.
     *
     * @param travel Total stroke [m]; finite.
     * @param width Depth of each buffer zone [m]; finite, positive, and `2 * width <= travel`.
     * @param stiffness Spring scale `k` [N/m]; finite and non-negative.
     * @param damping Damping scale `c` [N s/m]; finite and non-negative.
     * @return The buffer, or an `InvalidInput` fault.
     */
    [[nodiscard]] static constexpr std::expected<TravelBuffer, Fault> create(sub0math::f64::Length travel,
                                                                             sub0math::f64::Length width,
                                                                             sub0math::f64::Stiffness stiffness,
                                                                             sub0math::f64::Damping damping) noexcept
    {
        const double t = travel.si();
        const double w = width.si();
        const double k = stiffness.si();
        const double c = damping.si();
        if (!std::isfinite(t) || !std::isfinite(w) || !(w > 0.0))
        {
            return std::unexpected(
                Fault{FaultKind::InvalidInput, "travel buffer needs a finite travel and a positive width"});
        }
        if (!(2.0 * w <= t))
        {
            return std::unexpected(Fault{FaultKind::InvalidInput, "travel buffer zones overlap: need 2 * width <= travel"});
        }
        if (!std::isfinite(k) || !(k >= 0.0) || !std::isfinite(c) || !(c >= 0.0))
        {
            return std::unexpected(
                Fault{FaultKind::InvalidInput, "travel buffer needs finite non-negative stiffness and damping"});
        }
        return TravelBuffer{t, w, k, c};
    }

    /** @return The stroke. */
    [[nodiscard]] constexpr sub0math::f64::Length travel() const noexcept
    {
        return sub0math::f64::Length::fromSi(travel_);
    }

    /** @return The depth of each buffer zone. */
    [[nodiscard]] constexpr sub0math::f64::Length width() const noexcept
    {
        return sub0math::f64::Length::fromSi(width_);
    }

    /** @return The spring scale `k`. */
    [[nodiscard]] constexpr sub0math::f64::Stiffness stiffness() const noexcept
    {
        return sub0math::f64::Stiffness::fromSi(stiffness_);
    }

    /** @return The damping scale `c`. */
    [[nodiscard]] constexpr sub0math::f64::Damping damping() const noexcept
    {
        return sub0math::f64::Damping::fromSi(damping_);
    }

    /**
     * Evaluates the buffer.
     *
     * @param position Compression `q` [m]; 0 is fully extended, `travel` fully compressed.
     * @param velocity Compression rate `u = dq/dt` [m/s].
     * @return Force, energy, loss power, per-buffer forces and the two Jacobian entries.
     */
    [[nodiscard]] SUB0MATH_FORCE_INLINE TravelBufferResponse respond(sub0math::f64::Length position,
                                                                     sub0math::f64::Velocity velocity) const noexcept
    {
        const double q = position.si();
        const double v = velocity.si();
        const double k = stiffness_;
        const double c = damping_;
        const double invWidth = 1.0 / width_;

        // PERF: fmax measured equal to compare-select; kept.
        const double top = std::fmax(0.0, width_ - q);
        const double bottom = std::fmax(0.0, q - (travel_ - width_));
        // A buffer is damped only while its own side closes: u < 0 closes the top, u >= 0 the bottom.
        const bool closingTop = v < 0.0;
        const double speed = closingTop ? -v : v;
        const double topDamping = closingTop ? c * speed : 0.0;
        const double bottomDamping = closingTop ? 0.0 : c * speed;
        const double closingDepth = closingTop ? top : bottom;

        const double topForce = (k * top + topDamping) * top * invWidth;
        const double bottomForce = (k * bottom + bottomDamping) * bottom * invWidth;

        // d(topForce)/dq = -(2 k top + topDamping) / width; d(bottomForce)/dq = +(2 k bottom + bottomDamping) / width,
        // each only while its buffer is engaged (depth > 0).
        const double slopeTop = top > 0.0 ? 2.0 * k * top + topDamping : 0.0;
        const double slopeBottom = bottom > 0.0 ? 2.0 * k * bottom + bottomDamping : 0.0;

        return TravelBufferResponse{
            sub0math::f64::Force::fromSi(topForce - bottomForce),
            sub0math::f64::Energy::fromSi(k * (top * top * top + bottom * bottom * bottom) * invWidth / 3.0),
            sub0math::f64::Power::fromSi(c * speed * speed * closingDepth * invWidth),
            sub0math::f64::Force::fromSi(topForce),
            sub0math::f64::Force::fromSi(bottomForce),
            sub0math::f64::Stiffness::fromSi(-(slopeTop + slopeBottom) * invWidth),
            sub0math::f64::Damping::fromSi(-c * closingDepth * invWidth),
        };
    }

private:
    constexpr TravelBuffer(double travel, double width, double stiffness, double damping) noexcept
        : travel_(travel), width_(width), stiffness_(stiffness), damping_(damping)
    {
    }

    double travel_;    ///< Stroke [m].
    double width_;     ///< Zone depth [m].
    double stiffness_; ///< Spring scale [N/m].
    double damping_;   ///< Damping scale [N s/m].
};

} // namespace sub0physics
