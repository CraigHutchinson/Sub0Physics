#pragma once

#include <expected>
#include <span>
#include <sub0math/geometry/point.hpp>
#include <sub0math/units/quantities.hpp>
#include <sub0math/units/quantity.hpp>
#include <vector>

#include "sub0physics/fault.hpp"

namespace sub0physics
{

/** A profile vertex: x along the ground, y up. */
using ProfilePoint = sub0math::Point2<sub0math::f64::Length>;

/**
 * A rigid ground profile: an ordered polyline with solid below it.
 *
 * x never decreases, so the profile is a height field that also allows vertical faces (square edges)
 * but no overhangs. Validity is established once by `create`; every query then trusts it.
 */
class PolylineProfile
{
public:
    /** Shortest segment accepted; shorter ones make normals meaningless. */
    static constexpr sub0math::f64::Length MinimumSegmentLength = sub0math::f64::Length::fromSi(1e-6);
    /** Largest coordinate magnitude accepted; keeps squared distances finite. */
    static constexpr sub0math::f64::Length CoordinateLimit = sub0math::f64::Length::fromSi(1e100);

    /**
     * Validates points and builds a profile.
     *
     * Rejects fewer than two points, non-finite coordinates or ones beyond CoordinateLimit, decreasing x, segments shorter than
     * MinimumSegmentLength, consecutive vertical segments that reverse direction, and a profile with
     * no horizontal extent.
     *
     * @param points Vertices in order along x.
     * @return The profile, or an InvalidInput fault naming the first violated rule.
     */
    [[nodiscard]] static std::expected<PolylineProfile, Fault> create(std::vector<ProfilePoint> points);

    /** @return The vertices in order. */
    [[nodiscard]] std::span<const ProfilePoint> points() const noexcept
    {
        return points_;
    }

    /** @return Smallest x covered. */
    [[nodiscard]] sub0math::f64::Length minX() const noexcept
    {
        return points_.front().x;
    }

    /** @return Largest x covered. */
    [[nodiscard]] sub0math::f64::Length maxX() const noexcept
    {
        return points_.back().x;
    }

    /**
     * Finds the top of the solid at one horizontal position.
     *
     * Where a vertical face passes through `x`, the top of the face is returned. Runs in
     * O(log n + k) for k segments touching `x`.
     *
     * @param x Horizontal position; inside [minX(), maxX()].
     * @return Surface height; negative infinity outside the profile.
     */
    [[nodiscard]] sub0math::f64::Length surfaceHeight(sub0math::f64::Length x) const noexcept;

private:
    explicit PolylineProfile(std::vector<ProfilePoint> points) noexcept;

    std::vector<ProfilePoint> points_;
};

} // namespace sub0physics
