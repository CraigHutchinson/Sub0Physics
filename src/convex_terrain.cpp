#include "sub0physics/geometry/convex_terrain.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <cmath>
#include <limits>
#include <sub0math/config.hpp>
#include <sub0math/geometry/capsule_contact.hpp>

#include "sub0physics/geometry/norm.hpp"

namespace sub0physics
{

namespace
{

using sub0math::f64::Length;
using Capsule2d = sub0math::Capsule<double, 2>;

constexpr double Infinity = std::numeric_limits<double>::infinity();
constexpr std::uint64_t ExponentBits = 0x7ff0000000000000ULL; // All set: an infinity or a NaN.

Capsule2d toSi(const CasingEdge& edge) noexcept
{
    return {{{edge.axis.a.x.si(), edge.axis.a.y.si()}, {edge.axis.b.x.si(), edge.axis.b.y.si()}}, edge.radius.si()};
}

bool finite(const Capsule2d& edge) noexcept
{
    return std::isfinite(edge.axis.a.x) && std::isfinite(edge.axis.a.y) && std::isfinite(edge.axis.b.x) &&
           std::isfinite(edge.axis.b.y) && std::isfinite(edge.radius) && edge.radius >= 0.0;
}

bool nonNegative(double v) noexcept
{
    return std::isfinite(v) && v >= 0.0;
}

bool validMotion(const StepMotionBound& motion) noexcept
{
    return std::isfinite(motion.hubCentre.x.si()) && std::isfinite(motion.hubCentre.y.si()) &&
           nonNegative(motion.hubTravel.si()) && nonNegative(motion.hubRotation) && nonNegative(motion.nodeTravel.si());
}

} // namespace

ConvexTerrain::ConvexTerrain(const PolylineProfile& profile)
{
    // Collinear merging, as the reference: drop the middle vertex while the turn is negligible and not a reversal.
    for (const ProfilePoint& p : profile.points())
    {
        const double cx = p.x.si();
        const double cy = p.y.si();
        while (x_.size() > 1)
        {
            const std::size_t n = x_.size();
            const double ux = x_[n - 1] - x_[n - 2];
            const double uy = y_[n - 1] - y_[n - 2];
            const double vx = cx - x_[n - 1];
            const double vy = cy - y_[n - 1];
            if (std::abs(ux * vy - uy * vx) > CollinearSine * compensatedNorm(ux, uy) * compensatedNorm(vx, vy) ||
                ux * vx + uy * vy <= 0.0)
            {
                break;
            }
            x_.pop_back();
            y_.pop_back();
        }
        x_.push_back(cx);
        y_.push_back(cy);
    }

    floor_ = *std::min_element(y_.begin(), y_.end()) - FloorDepth;
    for (std::size_t s = 0; s + 1 < x_.size(); ++s)
    {
        const double ax = x_[s];
        const double ay = y_[s];
        const double bx = x_[s + 1];
        const double by = y_[s + 1];
        if (bx == ax)
        {
            continue; // A vertical face: the side of the neighbouring piece that rises above it.
        }
        const sub0math::Vec2d top = segmentNormal(s);
        // Counter-clockwise from the bottom left; edge normals: floor, right side, top (reversed segment), left side.
        pieces_.push_back(TerrainPiece{{sub0math::Point2d{ax, floor_}, sub0math::Point2d{bx, floor_},
                                        sub0math::Point2d{bx, by}, sub0math::Point2d{ax, ay}},
                                       {sub0math::Vec2d{0.0, -1.0}, sub0math::Vec2d{1.0, 0.0}, top,
                                        sub0math::Vec2d{-1.0, 0.0}}});
        pieceRight_.push_back(bx);
        pieceTop_.push_back(std::max(ay, by));
        highest_ = std::max(highest_, std::max(ay, by));
        segment_.push_back(s);
    }
}

sub0math::Vec2d ConvexTerrain::segmentNormal(std::size_t segment) const noexcept
{
    const double dx = x_[segment + 1] - x_[segment];
    const double dy = y_[segment + 1] - y_[segment];
    const double norm = compensatedNorm(dx, dy);
    return {-dy / norm, dx / norm};
}

sub0math::NormalCone<double> ConvexTerrain::cornerCone(std::size_t vertex, bool& convex) const noexcept
{
    const sub0math::Vec2d before = vertex == 0 ? sub0math::Vec2d{-1.0, 0.0} : segmentNormal(vertex - 1);
    const sub0math::Vec2d after = vertex + 1 >= x_.size() ? sub0math::Vec2d{1.0, 0.0} : segmentNormal(vertex);
    convex = after.x * before.y - after.y * before.x > 0.0;
    return {after, before};
}

void ConvexTerrain::window(double lo, double hi, std::size_t& first, std::size_t& last) const noexcept
{
    first = static_cast<std::size_t>(std::lower_bound(pieceRight_.begin(), pieceRight_.end(), lo) - pieceRight_.begin());
    last = first;
    while (last < pieces_.size() && pieces_[last].vertices[0].x <= hi)
    {
        ++last;
    }
}

// PERF: scalar arguments; a Capsule by reference stalled store forwarding.
SUB0MATH_FORCE_INLINE void ConvexTerrain::edgeContacts(double ax, double ay, double bx, double by, double radius, double margin,
                                 std::span<TerrainContact> out, std::size_t edgeIndex, TerrainContactSet& set,
                                 double& minimum) const noexcept
{
    const double pad = radius + margin;
    const double low = std::min(ay, by) - pad;
    if (low > highest_)
    {
        return; // Above every piece by more than the pad.
    }
    const double high = std::max(ay, by) + pad;
    std::size_t first = 0;
    std::size_t last = 0;
    window(std::min(ax, bx) - pad, std::max(ax, bx) + pad, first, last);
    const Capsule2d e{{{ax, ay}, {bx, by}}, radius};
    for (std::size_t i = first; i < last; ++i)
    {
        if (low > pieceTop_[i] || high < floor_)
        {
            continue; // Every point of the piece is further than the pad, vertically.
        }
        ++set.pairsTested;
        const sub0math::CapsuleContact<double, 2> c = sub0math::contact(e, pieces_[i]);
        minimum = std::min(minimum, c.gap);
        if (!(c.gap <= margin))
        {
            continue;
        }
        ++set.found;
        if (set.count == out.size())
        {
            set.overflowed = true;
            continue;
        }
        out[set.count++] = {edgeIndex,
                            i,
                            segment_[i],
                            Length::fromSi(c.gap),
                            c.axisT,
                            c.axisCrosses,
                            {Length::fromSi(c.capsulePoint.x), Length::fromSi(c.capsulePoint.y)},
                            {{Length::fromSi(c.witness.point.x), Length::fromSi(c.witness.point.y)},
                             c.witness.normal,
                             c.witness.kind,
                             c.witness.index}};
    }
}

TerrainContactSet ConvexTerrain::contacts(const CasingEdge& edge, Length margin, std::span<TerrainContact> out,
                                          std::size_t edgeIndex) const noexcept
{
    TerrainContactSet set{TerrainQueryStatus::Available, 0, 0, false, 0, Length::fromSi(Infinity)};
    const Capsule2d e = toSi(edge);
    if (!finite(e) || !nonNegative(margin.si()))
    {
        set.status = TerrainQueryStatus::InvalidQuery;
        return set;
    }
    double minimum = Infinity;
    edgeContacts(e.axis.a.x, e.axis.a.y, e.axis.b.x, e.axis.b.y, e.radius, margin.si(), out, edgeIndex, set, minimum);
    set.minimumGap = Length::fromSi(minimum);
    return set;
}

TerrainContactSet ConvexTerrain::ringContacts(std::span<const sub0math::Point2<Length>> nodes, Length skin,
                                              Length margin, std::span<TerrainContact> out) const noexcept
{
    TerrainContactSet set{TerrainQueryStatus::Available, 0, 0, false, 0, Length::fromSi(Infinity)};
    // One pass over the ring: exponent bits all set is an infinity or a NaN.
    std::uint64_t nonFinite = 0;
    for (const auto& node : nodes)
    {
        nonFinite |= (~std::bit_cast<std::uint64_t>(node.x.si()) & ExponentBits) == 0 ? 1U : 0U;
        nonFinite |= (~std::bit_cast<std::uint64_t>(node.y.si()) & ExponentBits) == 0 ? 1U : 0U;
    }
    if (nodes.size() < 2 || nonFinite != 0 || !nonNegative(skin.si()) || !nonNegative(margin.si()))
    {
        set.status = TerrainQueryStatus::InvalidQuery;
        return set;
    }
    double minimum = Infinity;
    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        const std::size_t next = i + 1 == nodes.size() ? 0 : i + 1;
        edgeContacts(nodes[i].x.si(), nodes[i].y.si(), nodes[next].x.si(), nodes[next].y.si(), skin.si(), margin.si(), out,
                     i, set, minimum);
    }
    set.minimumGap = Length::fromSi(minimum);
    return set;
}

SUB0MATH_FORCE_INLINE double ConvexTerrain::edgeClearance(double ax, double ay, double bx, double by, double radius,
                                                         const StepMotionBound& motion, double horizon,
                                                         double& reach) const noexcept
{
    // Largest distance of any point of the edge from the hub, then the movement bound of the class comment.
    const double hx = motion.hubCentre.x.si();
    const double hy = motion.hubCentre.y.si();
    const double lever = std::sqrt(std::max((ax - hx) * (ax - hx) + (ay - hy) * (ay - hy),
                                            (bx - hx) * (bx - hx) + (by - hy) * (by - hy)));
    reach = motion.hubTravel.si() + motion.nodeTravel.si() + lever * std::min(motion.hubRotation, 2.0);
    // Pieces beyond the pad have a gap above reach + horizon, so they cannot lower the capped bound.
    return edgeGap(ax, ay, bx, by, radius, reach + horizon) - reach;
}

SUB0MATH_FORCE_INLINE double ConvexTerrain::edgeGap(double ax, double ay, double bx, double by, double radius,
                                                   double cap) const noexcept
{
    const double pad = radius + cap;
    const double low = std::min(ay, by) - pad;
    double minimum = cap;
    if (low <= highest_)
    {
        std::size_t first = 0;
        std::size_t last = 0;
        window(std::min(ax, bx) - pad, std::max(ax, bx) + pad, first, last);
        const Capsule2d e{{{ax, ay}, {bx, by}}, radius};
        for (std::size_t i = first; i < last; ++i)
        {
            if (low > pieceTop_[i])
            {
                continue;
            }
            minimum = std::min(minimum, sub0math::contact(e, pieces_[i]).gap);
        }
    }
    return minimum;
}

ClearanceBound ConvexTerrain::certifyClearance(const CasingEdge& edge, const StepMotionBound& motion,
                                               Length horizon) const noexcept
{
    const Capsule2d e = toSi(edge);
    if (!finite(e) || !validMotion(motion) || !nonNegative(horizon.si()))
    {
        return {TerrainQueryStatus::InvalidQuery, Length::fromSi(0.0), Length::fromSi(0.0), false};
    }
    double reach = 0.0;
    const double clearance =
        edgeClearance(e.axis.a.x, e.axis.a.y, e.axis.b.x, e.axis.b.y, e.radius, motion, horizon.si(), reach);
    return {TerrainQueryStatus::Available, Length::fromSi(clearance), Length::fromSi(reach), clearance > 0.0};
}

ClearanceBound ConvexTerrain::certifyRingClearance(std::span<const sub0math::Point2<Length>> nodes, Length skin,
                                                   const StepMotionBound& motion, Length horizon) const noexcept
{
    std::uint64_t nonFinite = 0;
    for (const auto& node : nodes)
    {
        nonFinite |= (~std::bit_cast<std::uint64_t>(node.x.si()) & ExponentBits) == 0 ? 1U : 0U;
        nonFinite |= (~std::bit_cast<std::uint64_t>(node.y.si()) & ExponentBits) == 0 ? 1U : 0U;
    }
    if (nodes.size() < 2 || nonFinite != 0 || !nonNegative(skin.si()) || !validMotion(motion) ||
        !nonNegative(horizon.si()))
    {
        return {TerrainQueryStatus::InvalidQuery, Length::fromSi(0.0), Length::fromSi(0.0), false};
    }
    double clearance = horizon.si();
    double largestReach = 0.0;
    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        const std::size_t next = i + 1 == nodes.size() ? 0 : i + 1;
        double reach = 0.0;
        clearance = std::min(clearance, edgeClearance(nodes[i].x.si(), nodes[i].y.si(), nodes[next].x.si(),
                                                      nodes[next].y.si(), skin.si(), motion, horizon.si(), reach));
        largestReach = std::max(largestReach, reach);
    }
    return {TerrainQueryStatus::Available, Length::fromSi(clearance), Length::fromSi(largestReach), clearance > 0.0};
}

ClearanceBound ConvexTerrain::certifySweptRingClearance(std::span<const sub0math::Point2<Length>> start,
                                                        std::span<const sub0math::Point2<Length>> end, Length skin,
                                                        const SweptStepMotion& motion, Length horizon,
                                                        std::span<double> edgeBounds) const noexcept
{
    std::uint64_t nonFinite = 0;
    for (const auto& node : start)
    {
        nonFinite |= (~std::bit_cast<std::uint64_t>(node.x.si()) & ExponentBits) == 0 ? 1U : 0U;
        nonFinite |= (~std::bit_cast<std::uint64_t>(node.y.si()) & ExponentBits) == 0 ? 1U : 0U;
    }
    for (const auto& node : end)
    {
        nonFinite |= (~std::bit_cast<std::uint64_t>(node.x.si()) & ExponentBits) == 0 ? 1U : 0U;
        nonFinite |= (~std::bit_cast<std::uint64_t>(node.y.si()) & ExponentBits) == 0 ? 1U : 0U;
    }
    const double cx0 = motion.hubStart.x.si();
    const double cy0 = motion.hubStart.y.si();
    const double cx1 = motion.hubEnd.x.si();
    const double cy1 = motion.hubEnd.y.si();
    const double travel = std::sqrt((cx1 - cx0) * (cx1 - cx0) + (cy1 - cy0) * (cy1 - cy0));
    const double rotation = motion.hubRotation;
    const double nodeTravel = motion.nodeTravel.si();
    if (start.size() < 2 || end.size() != start.size() || nonFinite != 0 || !nonNegative(skin.si()) ||
        !std::isfinite(travel) || !nonNegative(rotation) || !nonNegative(nodeTravel) || !nonNegative(horizon.si()) ||
        (!edgeBounds.empty() && edgeBounds.size() != start.size()))
    {
        return {TerrainQueryStatus::InvalidQuery, Length::fromSi(0.0), Length::fromSi(0.0), false};
    }
    const double h = horizon.si();
    double clearance = h;
    double largestReach = 0.0;
    const std::size_t n = start.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        const std::size_t next = i + 1 == n ? 0 : i + 1;
        const double ax0 = start[i].x.si();
        const double ay0 = start[i].y.si();
        const double bx0 = start[next].x.si();
        const double by0 = start[next].y.si();
        const double ax1 = end[i].x.si();
        const double ay1 = end[i].y.si();
        const double bx1 = end[next].x.si();
        const double by1 = end[next].y.si();
        // The lever from each end's own hub centre (a point's hub-frame position at that end).
        const double start2 = std::max((ax0 - cx0) * (ax0 - cx0) + (ay0 - cy0) * (ay0 - cy0),
                                       (bx0 - cx0) * (bx0 - cx0) + (by0 - cy0) * (by0 - cy0));
        const double end2 = std::max((ax1 - cx1) * (ax1 - cx1) + (ay1 - cy1) * (ay1 - cy1),
                                     (bx1 - cx1) * (bx1 - cx1) + (by1 - cy1) * (by1 - cy1));
        const double lever = std::sqrt(std::max(start2, end2));
        const double reach = travel + nodeTravel + lever * rotation;
        const double cap = reach + h;
        const double g0 = edgeGap(ax0, ay0, bx0, by0, skin.si(), cap);
        const double g1 = edgeGap(ax1, ay1, bx1, by1, skin.si(), cap);
        const double bound = std::min(h, 0.5 * (g0 + g1 - reach));
        if (!edgeBounds.empty())
        {
            edgeBounds[i] = bound;
        }
        clearance = std::min(clearance, bound);
        largestReach = std::max(largestReach, reach);
    }
    return {TerrainQueryStatus::Available, Length::fromSi(clearance), Length::fromSi(largestReach), clearance > 0.0};
}

} // namespace sub0physics
