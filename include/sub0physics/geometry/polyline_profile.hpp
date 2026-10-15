#pragma once

#include <expected>
#include <span>
#include <vector>

#include "sub0physics/fault.hpp"
#include "sub0physics/geometry/point2.hpp"

namespace sub0physics
{

/**
 * A rigid ground profile: an ordered polyline with solid below it.
 *
 * x never decreases, so the profile is a height field that also allows vertical faces (square edges)
 * but no overhangs. Validity is established once by `create`; every query then trusts it.
 */
class PolylineProfile
{
public:
    /** Shortest segment accepted, m; shorter ones make normals meaningless. */
    static constexpr double MinimumSegmentLength = 1e-6;

    /**
     * Validates points and builds a profile.
     *
     * Rejects fewer than two points, non-finite coordinates, decreasing x, segments shorter than
     * MinimumSegmentLength, consecutive vertical segments that reverse direction, and a profile with
     * no horizontal extent.
     *
     * @param points Vertices in order along x.
     * @return The profile, or an InvalidInput fault naming the first violated rule.
     */
    [[nodiscard]] static std::expected<PolylineProfile, Fault> create(std::vector<Point2> points);

    /** @return The vertices in order. */
    [[nodiscard]] std::span<const Point2> points() const noexcept
    {
        return points_;
    }

    /** @return Smallest x covered, m. */
    [[nodiscard]] double minX() const noexcept
    {
        return points_.front().x;
    }

    /** @return Largest x covered, m. */
    [[nodiscard]] double maxX() const noexcept
    {
        return points_.back().x;
    }

    /**
     * Finds the top of the solid at one horizontal position.
     *
     * Where a vertical face passes through `x`, the top of the face is returned. Runs in
     * O(log n + k) for k segments touching `x`.
     *
     * @param x Horizontal position, m; inside [minX(), maxX()].
     * @return Surface height, m; negative infinity outside the profile.
     */
    [[nodiscard]] double surfaceHeight(double x) const noexcept;

private:
    explicit PolylineProfile(std::vector<Point2> points) noexcept;

    std::vector<Point2> points_;
};

} // namespace sub0physics
