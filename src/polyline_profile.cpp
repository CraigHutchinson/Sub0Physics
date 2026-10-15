#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include "sub0physics/geometry/polyline_profile.hpp"

namespace sub0physics
{

namespace
{

std::unexpected<Fault> invalid(std::string_view detail) noexcept
{
    return std::unexpected(Fault{FaultKind::InvalidInput, detail});
}

} // namespace

std::expected<PolylineProfile, Fault> PolylineProfile::create(std::vector<Point2> points)
{
    if (points.size() < 2U)
    {
        return invalid("profile needs at least two points");
    }
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        const Point2 p = points[i];
        if (!(std::abs(p.x) <= CoordinateLimit) || !(std::abs(p.y) <= CoordinateLimit))
        {
            return invalid("profile coordinates must be finite and within the coordinate limit");
        }
        if (i == 0U)
        {
            continue;
        }
        const double dx = p.x - points[i - 1U].x;
        const double dy = p.y - points[i - 1U].y;
        if (dx < 0.0)
        {
            return invalid("profile x must not decrease; overhangs are unsupported");
        }
        if (std::hypot(dx, dy) < MinimumSegmentLength)
        {
            return invalid("profile segments must be at least 1 micrometre long");
        }
        if (i > 1U && dx == 0.0 && points[i - 1U].x == points[i - 2U].x
            && dy * (points[i - 1U].y - points[i - 2U].y) < 0.0)
        {
            return invalid("vertical profile segments must not reverse or overlap");
        }
    }
    if (points.back().x - points.front().x < MinimumSegmentLength)
    {
        return invalid("profile must span a horizontal domain");
    }
    return PolylineProfile(std::move(points));
}

PolylineProfile::PolylineProfile(std::vector<Point2> points) noexcept
    : points_(std::move(points))
{
}

double PolylineProfile::surfaceHeight(double x) const noexcept
{
    double surface = -std::numeric_limits<double>::infinity();
    // Segments touching x are contiguous because x never decreases: start at the last vertex left of x.
    const auto first = std::lower_bound(points_.begin(), points_.end(), x,
                                        [](const Point2& p, double value) { return p.x < value; });
    std::size_t i = first == points_.begin() ? 0U : static_cast<std::size_t>(first - points_.begin()) - 1U;
    for (; i + 1U < points_.size() && points_[i].x <= x; ++i)
    {
        const Point2 a = points_[i];
        const Point2 b = points_[i + 1U];
        if (x > b.x)
        {
            continue;
        }
        const double dx = b.x - a.x;
        const double height = dx == 0.0 ? std::max(a.y, b.y) : a.y + (x - a.x) * (b.y - a.y) / dx;
        surface = std::max(surface, height);
    }
    return surface;
}

} // namespace sub0physics
