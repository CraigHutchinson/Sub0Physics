#pragma once

#include <algorithm>
#include <cmath>
#include <expected>
#include <sub0math/config.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/element/law_response.hpp"
#include "sub0physics/fault.hpp"

namespace sub0physics
{

/**
 * Response of a TravelBuffer: the common law response plus the force of each buffer.
 *
 * Sign convention: `force` is the generalized force on the compression coordinate `q`. The top buffer
 * pushes toward larger `q` (positive), the bottom buffer toward smaller `q` (negative), so
 * `force = topForce - bottomForce`. `potential` is the elastic energy of the two buffers and
 * `dissipationPower` the damper power, both never negative.
 */
struct TravelBufferResponse : LawResponse<sub0math::f64::Length>
{
    sub0math::f64::Force topForce;    ///< Top buffer (near `q = 0`) force magnitude; never negative.
    sub0math::f64::Force bottomForce; ///< Bottom buffer (near `q = travel`) force magnitude; never negative.
};

/**
 * The parameters of a TravelBuffer. Initialise with designated initialisers, so the two lengths (stroke
 * and zone depth) and the two scales cannot be swapped:
 * `{.travel = ..., .zoneWidth = ..., .stiffness = ..., .damping = ...}`.
 */
struct TravelBufferSpec
{
    sub0math::f64::Length travel;      ///< Total stroke [m]; finite.
    sub0math::f64::Length zoneWidth;   ///< Depth of each buffer zone [m]; finite, positive, and `2 * zoneWidth < travel`.
    sub0math::f64::Stiffness stiffness; ///< Spring scale `k` [N/m]; finite and non-negative.
    sub0math::f64::Damping damping;    ///< Damping scale `c` [N s/m]; finite and non-negative.
};

/**
 * Progressive compliant end-stops inside a finite stroke.
 *
 * Compression `q` runs from 0 (fully extended) to `travel`. A top buffer of depth `width` engages for
 * `q < width` and a bottom buffer of the same depth for `q > travel - width`; the zones never
 * overlap or touch (`2 * width < travel`, the BikeLab rule) and the free zone between them exerts nothing. With penetration depth `d` into a buffer the
 * spring force is `k d^2 / width` (progressive, not constant tangent stiffness) and the damper adds
 * `c s d / width`, where `s` is the speed of entering the buffer. Damping acts only while entering;
 * leaving a buffer is undamped.
 *
 * `potential` is the exact integral of the spring part, so `-d(potential)/dq` equals the spring part
 * of `force`. `dissipationPower` is the damper power, so `force * u = -d(potential)/dt - dissipationPower`.
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
     * @param spec The parameters; see TravelBufferSpec. Zones that touch (`2 * zoneWidth == travel`) or
     *        overlap are refused, as is any non-finite value.
     * @return The buffer, or an `InvalidInput` fault.
     */
    [[nodiscard]] static constexpr std::expected<TravelBuffer, Fault> create(const TravelBufferSpec& spec) noexcept
    {
        const double t = spec.travel.si();
        const double w = spec.zoneWidth.si();
        const double k = spec.stiffness.si();
        const double c = spec.damping.si();
        if (!std::isfinite(t) || !std::isfinite(w) || !(w > 0.0))
        {
            return invalidInput("travel buffer needs a finite travel and a positive zone width");
        }
        if (!(2.0 * w < t))
        {
            return invalidInput("travel buffer zones overlap or touch: need 2 * zoneWidth < travel");
        }
        if (!std::isfinite(k) || !(k >= 0.0) || !std::isfinite(c) || !(c >= 0.0))
        {
            return invalidInput("travel buffer needs finite non-negative stiffness and damping");
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
     * @return Force, the two Jacobian entries, energy, dissipated power and the per-buffer forces.
     */
    [[nodiscard]] SUB0MATH_FORCE_INLINE TravelBufferResponse respond(sub0math::f64::Length position,
                                                                     sub0math::f64::Velocity velocity) const noexcept
    {
        const double q = position.si();
        const double v = velocity.si();
        const double k = stiffness_;
        const double c = damping_;
        const double invWidth = 1.0 / width_;

        // PERF: std::max (maxsd) beats std::fmax (a libcall on MinGW GCC).
        const double top = std::max(0.0, width_ - q);
        const double bottom = std::max(0.0, q - (travel_ - width_));
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
            {
                sub0math::f64::Force::fromSi(topForce - bottomForce),
                sub0math::f64::Stiffness::fromSi(-(slopeTop + slopeBottom) * invWidth),
                sub0math::f64::Damping::fromSi(-c * closingDepth * invWidth),
                sub0math::f64::Energy::fromSi(k * (top * top * top + bottom * bottom * bottom) * invWidth / 3.0),
                sub0math::f64::Power::fromSi(c * speed * speed * closingDepth * invWidth),
            },
            sub0math::f64::Force::fromSi(topForce),
            sub0math::f64::Force::fromSi(bottomForce),
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

static_assert(Law<TravelBuffer, sub0math::f64::Length>);

} // namespace sub0physics
