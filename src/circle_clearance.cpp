#include <algorithm>
#include <limits>
#include <span>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/functions.hpp>

#include "sub0physics/geometry/circle_clearance.hpp"

namespace sub0physics
{

namespace
{

using sub0math::f64::Area;
using sub0math::f64::Length;

struct Nearest
{
    ProfilePoint point;
    Area distanceSquared;
    double t;
};

Nearest nearestOnSegment(const ProfilePoint& a, const ProfilePoint& b, const ProfilePoint& c) noexcept
{
    const sub0math::Vec2<Length> ab = b - a;
    const sub0math::Vec2<Length> ac = c - a;
    const double t = std::clamp(static_cast<double>(dot(ac, ab) / dot(ab, ab)), 0.0, 1.0);
    const ProfilePoint p{a.x + ab.x * t, a.y + ab.y * t};
    // Squared distance: profiles and queries are bounded by PolylineProfile::CoordinateLimit, so no
    // overflow guard is needed, and square roots are taken only when the nearest improves.
    return Nearest{p, lengthSquared(c - p), t};
}

// Index of the first segment whose right end is at or beyond x.
std::size_t segmentAt(std::span<const ProfilePoint> points, Length x) noexcept
{
    const auto it = std::lower_bound(points.begin() + 1, points.end(), x,
                                     [](const ProfilePoint& p, Length value) { return p.x < value; });
    return std::min(static_cast<std::size_t>(it - points.begin()), points.size() - 1U) - 1U;
}

} // namespace

CircleClearance measureCircleClearance(const PolylineProfile& profile, const CircleQuery& query)
{
    const ProfilePoint c = query.centre;
    const Length r = query.radius;
    CircleClearance result{ClearanceStatus::Available, Length{}, Length{}, {}};
    if (!(sub0math::abs(c.y) <= PolylineProfile::CoordinateLimit) || !(r > Length{}) || c.x - r < profile.minX()
        || c.x + r > profile.maxX())
    {
        result.status = ClearanceStatus::OutsideDomain;
        return result;
    }

    // x never decreases, so segments farther away horizontally than the best distance so far cannot
    // be nearer. Search outward from the segment under the centre and stop at that bound.
    const auto points = profile.points();
    const std::size_t segments = points.size() - 1U;
    const std::size_t start = segmentAt(points, c.x);
    Area best = Area::fromSi(std::numeric_limits<double>::infinity());
    Length distance = Length::fromSi(std::numeric_limits<double>::infinity());
    for (std::size_t i = start; i < segments && points[i].x - c.x <= distance; ++i)
    {
        if (const Area d2 = nearestOnSegment(points[i], points[i + 1U], c).distanceSquared; d2 < best)
        {
            best = d2;
            distance = sqrt(best);
        }
    }
    for (std::size_t i = start; i-- > 0U && c.x - points[i + 1U].x <= distance;)
    {
        if (const Area d2 = nearestOnSegment(points[i], points[i + 1U], c).distanceSquared; d2 < best)
        {
            best = d2;
            distance = sqrt(best);
        }
    }
    if (distance <= query.tieTolerance)
    {
        result.status = ClearanceStatus::OnBoundary;
        return result;
    }
    if (c.y < profile.surfaceHeight(c.x))
    {
        result.status = ClearanceStatus::CentreInside;
        return result;
    }

    // Second pass collects ties. A vertex is reached from its two adjacent segments in turn, so
    // comparing with the previous witness removes every duplicate.
    const Length reach = distance + query.tieTolerance;
    const Area reachSquared = reach * reach;
    std::size_t first = start;
    while (first > 0U && c.x - points[first].x <= reach)
    {
        --first;
    }
    for (std::size_t i = first; i < segments && points[i].x - c.x <= reach; ++i)
    {
        const Nearest n = nearestOnSegment(points[i], points[i + 1U], c);
        if (n.distanceSquared > reachSquared)
        {
            continue;
        }
        const bool atStart = n.t == 0.0;
        const bool atEnd = n.t == 1.0;
        const FeatureKind kind = atStart || atEnd ? FeatureKind::Vertex : FeatureKind::Segment;
        const std::size_t index = atEnd ? i + 1U : i;
        if (!result.witnesses.empty() && result.witnesses.back().kind == kind && result.witnesses.back().index == index)
        {
            continue;
        }
        const Length d = sqrt(n.distanceSquared);
        result.witnesses.push_back(ClearanceWitness{
            n.point,
            sub0math::Vec2d{(c.x - n.point.x) / d, (c.y - n.point.y) / d},
            kind,
            index,
        });
    }
    result.distance = distance;
    result.gap = distance - r;
    return result;
}

} // namespace sub0physics
