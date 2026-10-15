#include <algorithm>
#include <cmath>
#include <limits>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/functions.hpp>
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

using sub0math::f64::Length;

std::expected<PolylineProfile, Fault> PolylineProfile::create(std::vector<ProfilePoint> points)
{
    if (points.size() < 2U)
    {
        return invalid("profile needs at least two points");
    }
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        const ProfilePoint p = points[i];
        if (!(sub0math::abs(p.x) <= CoordinateLimit) || !(sub0math::abs(p.y) <= CoordinateLimit))
        {
            return invalid("profile coordinates must be finite and within the coordinate limit");
        }
        if (i == 0U)
        {
            continue;
        }
        const Length dx = p.x - points[i - 1U].x;
        const Length dy = p.y - points[i - 1U].y;
        if (dx < Length{})
        {
            return invalid("profile x must not decrease; overhangs are unsupported");
        }
        if (sub0math::length(p - points[i - 1U]) < MinimumSegmentLength)
        {
            return invalid("profile segments must be at least 1 micrometre long");
        }
        if (i > 1U && dx == Length{} && points[i - 1U].x == points[i - 2U].x
            && dy * (points[i - 1U].y - points[i - 2U].y) < sub0math::f64::Area{})
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

PolylineProfile::PolylineProfile(std::vector<ProfilePoint> points) noexcept
    : points_(std::move(points))
{
}

Length PolylineProfile::surfaceHeight(Length x) const noexcept
{
    Length surface = Length::fromSi(-std::numeric_limits<double>::infinity());
    // Segments touching x are contiguous because x never decreases: start at the last vertex left of x.
    const auto first = std::lower_bound(points_.begin(), points_.end(), x,
                                        [](const ProfilePoint& p, Length value) { return p.x < value; });
    std::size_t i = first == points_.begin() ? 0U : static_cast<std::size_t>(first - points_.begin()) - 1U;
    for (; i + 1U < points_.size() && points_[i].x <= x; ++i)
    {
        const ProfilePoint a = points_[i];
        const ProfilePoint b = points_[i + 1U];
        if (x > b.x)
        {
            continue;
        }
        const Length dx = b.x - a.x;
        const Length height = dx == Length{} ? std::max(a.y, b.y) : a.y + (x - a.x) * (b.y - a.y) / dx;
        surface = std::max(surface, height);
    }
    return surface;
}

} // namespace sub0physics
