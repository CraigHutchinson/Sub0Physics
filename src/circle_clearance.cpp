#include <algorithm>
#include <cmath>
#include <limits>
#include <span>

#include "sub0physics/geometry/circle_clearance.hpp"

namespace sub0physics
{

namespace
{

struct Nearest
{
    Point2 point;
    double distanceSquared;
    double t;
};

Nearest nearestOnSegment(Point2 a, Point2 b, Point2 c) noexcept
{
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double t = std::clamp(((c.x - a.x) * dx + (c.y - a.y) * dy) / (dx * dx + dy * dy), 0.0, 1.0);
    const Point2 p{a.x + t * dx, a.y + t * dy};
    const double dxc = c.x - p.x;
    const double dyc = c.y - p.y;
    // Squared distance: profiles and queries are bounded by PolylineProfile::CoordinateLimit, so no
    // std::hypot overflow guard is needed, and square roots are taken only when the nearest improves.
    return Nearest{p, dxc * dxc + dyc * dyc, t};
}

// Index of the first segment whose right end is at or beyond x.
std::size_t segmentAt(std::span<const Point2> points, double x) noexcept
{
    const auto it = std::lower_bound(points.begin() + 1, points.end(), x,
                                     [](const Point2& p, double value) { return p.x < value; });
    return std::min(static_cast<std::size_t>(it - points.begin()), points.size() - 1U) - 1U;
}

} // namespace

CircleClearance measureCircleClearance(const PolylineProfile& profile, const CircleQuery& query)
{
    const Point2 c = query.centre;
    const double r = query.radius;
    CircleClearance result{ClearanceStatus::Available, 0.0, 0.0, {}};
    if (!(std::abs(c.y) <= PolylineProfile::CoordinateLimit) || !(r > 0.0) || c.x - r < profile.minX()
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
    double best = std::numeric_limits<double>::infinity();
    double distance = best;
    for (std::size_t i = start; i < segments && points[i].x - c.x <= distance; ++i)
    {
        if (const double d2 = nearestOnSegment(points[i], points[i + 1U], c).distanceSquared; d2 < best)
        {
            best = d2;
            distance = std::sqrt(best);
        }
    }
    for (std::size_t i = start; i-- > 0U && c.x - points[i + 1U].x <= distance;)
    {
        if (const double d2 = nearestOnSegment(points[i], points[i + 1U], c).distanceSquared; d2 < best)
        {
            best = d2;
            distance = std::sqrt(best);
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
    const double reach = distance + query.tieTolerance;
    const double reachSquared = reach * reach;
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
        const double d = std::sqrt(n.distanceSquared);
        result.witnesses.push_back(ClearanceWitness{
            n.point,
            Point2{(c.x - n.point.x) / d, (c.y - n.point.y) / d},
            kind,
            index,
        });
    }
    result.distance = distance;
    result.gap = distance - r;
    return result;
}

} // namespace sub0physics
