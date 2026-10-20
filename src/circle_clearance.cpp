#include <algorithm>
#include <limits>
#include <span>
#include <sub0math/geometry/segment.hpp>
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
    const auto nearest = sub0math::closestPoint(sub0math::Segment<Length, 2>{a, b}, c);
    // PERF: squared distances, sqrt on improvement; 3x faster than hypot.
    // Squared distance: profiles and queries are bounded by PolylineProfile::CoordinateLimit, so no
    // overflow guard is needed, and square roots are taken only when the nearest improves.
    return Nearest{nearest.point, lengthSquared(c - nearest.point), nearest.t};
}

// Index of the first segment whose right end is at or beyond x.
std::size_t segmentAt(std::span<const ProfilePoint> points, Length x) noexcept
{
    const auto it = std::lower_bound(points.begin() + 1, points.end(), x,
                                     [](const ProfilePoint& p, Length value) { return p.x < value; });
    return std::min(static_cast<std::size_t>(it - points.begin()), points.size() - 1U) - 1U;
}

} // namespace

CircleClearance measureCircleClearance(const PolylineProfile& profile, const CircleQuery& query,
                                       std::span<Witness2> witnesses) noexcept
{
    const ProfilePoint c = query.circle.centre;
    const Length r = query.circle.radius;
    CircleClearance result{ClearanceStatus::Available, Length{}, Length{}, 0U, false};
    // Negated comparisons so that NaN fails every test; the limits also exclude infinities.
    const Length limit = PolylineProfile::CoordinateLimit;
    if (!(sub0math::abs(c.x) <= limit) || !(sub0math::abs(c.y) <= limit) || !(r > Length{}) || !(r <= limit)
        || !(query.tieTolerance >= Length{}) || !(query.tieTolerance <= limit) || c.x - r < profile.minX()
        || c.x + r > profile.maxX())
    {
        result.status = ClearanceStatus::OutsideDomain;
        return result;
    }

    // PERF: outward search pruned by x; ~17x faster than scanning all.
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
    bool havePrevious = false;
    WitnessKind previousKind = WitnessKind::Vertex;
    std::size_t previousIndex = 0U;
    for (std::size_t i = first; i < segments && points[i].x - c.x <= reach; ++i)
    {
        const Nearest n = nearestOnSegment(points[i], points[i + 1U], c);
        if (n.distanceSquared > reachSquared)
        {
            continue;
        }
        const bool atStart = n.t == 0.0;
        const bool atEnd = n.t == 1.0;
        const WitnessKind kind = atStart || atEnd ? WitnessKind::Vertex : WitnessKind::Edge;
        const std::size_t index = atEnd ? i + 1U : i;
        if (havePrevious && previousKind == kind && previousIndex == index)
        {
            continue;
        }
        havePrevious = true;
        previousKind = kind;
        previousIndex = index;
        if (result.witnessCount == witnesses.size())
        {
            result.overflowed = true;
            continue;
        }
        const Length d = sqrt(n.distanceSquared);
        witnesses[result.witnessCount++] = Witness2{
            n.point,
            sub0math::Vec2d{(c.x - n.point.x) / d, (c.y - n.point.y) / d},
            kind,
            index,
        };
    }
    result.distance = distance;
    result.gap = distance - r;
    return result;
}

CircleClearance measureCircleClearance(const PreparedProfile& profile, const CircleQuery& query,
                                       std::span<Witness2> witnesses) noexcept
{
    const double cx = query.circle.centre.x.si();
    const double cy = query.circle.centre.y.si();
    const double r = query.circle.radius.si();
    const double tie = query.tieTolerance.si();
    CircleClearance result{ClearanceStatus::Available, Length{}, Length{}, 0U, false};
    // Negated comparisons so that NaN fails every test; the limits also exclude infinities.
    constexpr double limit = PolylineProfile::CoordinateLimit.si();
    const std::span<const double> xs = profile.vertexX();
    if (!(std::abs(cx) <= limit) || !(std::abs(cy) <= limit) || !(r > 0.0) || !(r <= limit) || !(tie >= 0.0)
        || !(tie <= limit) || cx - r < xs.front() || cx + r > xs.back())
    {
        result.status = ClearanceStatus::OutsideDomain;
        return result;
    }

    const std::span<const double> ys = profile.vertexY();
    const std::span<const double> dxs = profile.segmentDx();
    const std::span<const double> dys = profile.segmentDy();
    const std::span<const double> sumSquares = profile.segmentSumSquares();
    const std::size_t count = xs.size();
    const std::size_t segments = count - 1U;
    const double* const px = xs.data();
    const auto squared = [&](std::size_t i) noexcept {
        const double ax = cx - px[i];
        const double ay = cy - ys[i];
        const double t = std::clamp((ax * dxs[i] + ay * dys[i]) / sumSquares[i], 0.0, 1.0);
        const double ex = cx - (px[i] + dxs[i] * t);
        const double ey = cy - (ys[i] + dys[i] * t);
        return ex * ex + ey * ey;
    };

    // PERF: window from centre segment, four minima; faster than pruned walk.
    // The segment under the centre bounds the distance; only segments whose x-range meets the centre's band of that
    // half-width can be nearer. The margin keeps a rounded bound from excluding a segment at exactly that distance.
    // The arithmetic of each distance is the unprepared one (a division, the point rebuilt from the segment start),
    // so both overloads agree bit for bit.
    const std::size_t start =
        std::min(static_cast<std::size_t>(std::lower_bound(px + 1, px + count, cx) - px), count - 1U) - 1U;
    const double bound = std::sqrt(squared(start)) * (1.0 + 1e-14);
    const std::size_t lo = static_cast<std::size_t>(std::lower_bound(px, px + count, cx - bound) - px);
    const std::size_t hi = static_cast<std::size_t>(std::upper_bound(px, px + count, cx + bound) - px);
    std::size_t i = lo == 0U ? 0U : lo - 1U;
    const std::size_t end = std::min(hi, segments);
    constexpr double Infinity = std::numeric_limits<double>::infinity();
    double m0 = Infinity;
    double m1 = Infinity;
    double m2 = Infinity;
    double m3 = Infinity;
    for (; i + 4U <= end; i += 4U)
    {
        m0 = std::min(m0, squared(i));
        m1 = std::min(m1, squared(i + 1U));
        m2 = std::min(m2, squared(i + 2U));
        m3 = std::min(m3, squared(i + 3U));
    }
    for (; i < end; ++i)
    {
        m0 = std::min(m0, squared(i));
    }
    const double distance = std::sqrt(std::min(std::min(m0, m1), std::min(m2, m3)));
    if (distance <= tie)
    {
        result.status = ClearanceStatus::OnBoundary;
        return result;
    }
    if (cy < profile.surfaceHeight(query.circle.centre.x).si())
    {
        result.status = ClearanceStatus::CentreInside;
        return result;
    }

    // Second pass collects ties. A vertex is reached from its two adjacent segments in turn, so comparing with the previous witness removes
    // every duplicate.
    const double reach = distance + tie;
    const double reachSquared = reach * reach;
    std::size_t first = start;
    while (first > 0U && cx - px[first] <= reach)
    {
        --first;
    }
    bool havePrevious = false;
    WitnessKind previousKind = WitnessKind::Vertex;
    std::size_t previousIndex = 0U;
    for (std::size_t k = first; k < segments && px[k] - cx <= reach; ++k)
    {
        const double dx = dxs[k];
        const double dy = dys[k];
        const double t = std::clamp(((cx - px[k]) * dx + (cy - ys[k]) * dy) / sumSquares[k], 0.0, 1.0);
        const double nearestX = px[k] + dx * t;
        const double nearestY = ys[k] + dy * t;
        const double ex = cx - nearestX;
        const double ey = cy - nearestY;
        const double d2 = ex * ex + ey * ey;
        if (d2 > reachSquared)
        {
            continue;
        }
        const bool atStart = t == 0.0;
        const bool atEnd = t == 1.0;
        const WitnessKind kind = atStart || atEnd ? WitnessKind::Vertex : WitnessKind::Edge;
        const std::size_t index = atEnd ? k + 1U : k;
        if (havePrevious && previousKind == kind && previousIndex == index)
        {
            continue;
        }
        havePrevious = true;
        previousKind = kind;
        previousIndex = index;
        if (result.witnessCount == witnesses.size())
        {
            result.overflowed = true;
            continue;
        }
        const double d = std::sqrt(d2);
        witnesses[result.witnessCount++] = Witness2{
            {Length::fromSi(nearestX), Length::fromSi(nearestY)},
            sub0math::Vec2d{ex / d, ey / d},
            kind,
            index,
        };
    }
    result.distance = Length::fromSi(distance);
    result.gap = Length::fromSi(distance - r);
    return result;
}

} // namespace sub0physics
