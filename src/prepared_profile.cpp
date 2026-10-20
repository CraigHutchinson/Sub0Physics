#include <algorithm>
#include <cmath>
#include <limits>

#include "sub0physics/geometry/prepared_profile.hpp"

namespace sub0physics
{

using sub0math::Point2;
using sub0math::Vec2;
using sub0math::Vec2d;
using sub0math::f64::Length;

namespace
{

constexpr double Infinity = std::numeric_limits<double>::infinity();

// Length of (a, b) by scaling with the larger magnitude and summing the squares with Kahan compensation: the
// arithmetic of V8's Math.hypot. It does not depend on the C library, and ported goldens rely on its rounding.
double scaledNorm(double a, double b) noexcept
{
    const double largest = std::max(std::abs(a), std::abs(b));
    if (largest == 0.0)
    {
        return 0.0;
    }
    double sum = 0.0;
    double compensation = 0.0;
    for (const double v : {a, b})
    {
        const double scaled = std::abs(v) / largest;
        const double summand = scaled * scaled - compensation;
        const double preliminary = sum + summand;
        compensation = (preliminary - sum) - summand;
        sum = preliminary;
    }
    return std::sqrt(sum) * largest;
}

bool finite(double a, double b) noexcept
{
    return std::isfinite(a) && std::isfinite(b);
}

std::size_t lowerIndex(const std::vector<double>& xs, double value) noexcept
{
    return static_cast<std::size_t>(std::lower_bound(xs.begin(), xs.end(), value) - xs.begin());
}

std::size_t upperIndex(const std::vector<double>& xs, double value) noexcept
{
    return static_cast<std::size_t>(std::upper_bound(xs.begin(), xs.end(), value) - xs.begin());
}

} // namespace

PreparedProfile::PreparedProfile(const PolylineProfile& profile)
{
    const auto points = profile.points();
    const std::size_t segments = points.size() - 1U;
    x_.reserve(points.size());
    y_.reserve(points.size());
    for (const ProfilePoint& p : points)
    {
        x_.push_back(p.x.si());
        y_.push_back(p.y.si());
    }
    for (auto* array : {&dx_, &dy_, &sumSquares_, &norm_, &normSquared_, &nx_, &ny_})
    {
        array->reserve(segments);
    }
    for (std::size_t i = 0; i < segments; ++i)
    {
        const double dx = x_[i + 1U] - x_[i];
        const double dy = y_[i + 1U] - y_[i];
        const double sumSquares = dx * dx + dy * dy;
        const double norm = scaledNorm(dx, dy);
        dx_.push_back(dx);
        dy_.push_back(dy);
        sumSquares_.push_back(sumSquares);
        norm_.push_back(norm);
        normSquared_.push_back(norm * norm);
        nx_.push_back(-dy / norm);
        ny_.push_back(dx / norm);
        longestSegment_ = std::max(longestSegment_, norm);
    }
}

double PreparedProfile::surfaceAt(double x, std::size_t* feature) const noexcept
{
    // Same arithmetic as PolylineProfile::surfaceHeight; segments touching x are contiguous because x never decreases.
    double surface = -Infinity;
    std::size_t i = lowerIndex(x_, x);
    i = i == 0U ? 0U : i - 1U;
    for (; i + 1U < x_.size() && x_[i] <= x; ++i)
    {
        if (x > x_[i + 1U])
        {
            continue;
        }
        const double run = dx_[i];
        const double height = run == 0.0 ? std::max(y_[i], y_[i + 1U]) : y_[i] + (x - x_[i]) * dy_[i] / run;
        if (height > surface && feature != nullptr)
        {
            *feature = i;
        }
        surface = std::max(surface, height);
    }
    return surface;
}

Length PreparedProfile::surfaceHeight(Length x) const noexcept
{
    return Length::fromSi(surfaceAt(x.si(), nullptr));
}

void PreparedProfile::bindRayOrigin(Point2<Length> origin, RayOrigin& bound) const
{
    const std::size_t segments = segmentCount();
    const double ox = origin.x.si();
    const double oz = origin.y.si();
    bound.origin = origin;
    bound.ax.resize(segments);
    bound.az.resize(segments);
    bound.cross.resize(segments);
    for (std::size_t j = 0; j < segments; ++j)
    {
        const double ax = x_[j] - ox;
        const double az = y_[j] - oz;
        bound.ax[j] = ax;
        bound.az[j] = az;
        bound.cross[j] = ax * dy_[j] - az * dx_[j];
    }
}

template<typename RelativeAt>
RayHit PreparedProfile::castFrom(Point2<Length> origin, Vec2d direction, double tieTolerance,
                                 const RelativeAt& relativeAt) const noexcept
{
    const double ox = direction.x;
    const double oz = direction.y;
    RayHit hit{RayStatus::Miss, false, Infinity, Witness2{}};
    if (!finite(ox, oz) || !finite(origin.x.si(), origin.y.si()) || !(tieTolerance >= 0.0) || !std::isfinite(tieTolerance))
    {
        hit.status = RayStatus::InvalidQuery;
        return hit;
    }
    const std::size_t segments = segmentCount();
    const double* const sx = dx_.data();
    const double* const sz = dy_.data();
    const double* const snx = nx_.data();
    const double* const snz = ny_.data();
    const double upper = 1.0 + tieTolerance;
    double distance = Infinity;
    std::size_t feature = 0;
    bool ambiguous = false;
    bool found = false;
    for (std::size_t j = 0; j < segments; ++j)
    {
        const double incidence = snx[j] * ox + snz[j] * oz;
        if (incidence >= -GrazingIncidence)
        {
            continue;
        }
        const double den = ox * sz[j] - oz * sx[j];
        const Relative rel = relativeAt(j);
        const double r = rel.cross / den;
        const double t = (rel.ax * oz - rel.az * ox) / den;
        if (r <= 0.0 || t < -tieTolerance || t > upper || r > distance + tieTolerance)
        {
            continue;
        }
        if (std::abs(r - distance) < tieTolerance)
        {
            if (scaledNorm(snx[feature] - snx[j], snz[feature] - snz[j]) > DistinctNormal)
            {
                ambiguous = true;
            }
            continue;
        }
        distance = r;
        feature = j;
        ambiguous = false;
        found = true;
    }
    if (found)
    {
        hit.status = RayStatus::Hit;
        hit.ambiguous = ambiguous;
        hit.distance = distance;
        hit.witness = Witness2{{Length::fromSi(origin.x.si() + distance * ox), Length::fromSi(origin.y.si() + distance * oz)},
                               Vec2d{snx[feature], snz[feature]},
                               WitnessKind::Edge,
                               feature};
    }
    return hit;
}

RayHit PreparedProfile::castRay(const RayOrigin& bound, Vec2d direction, double tieTolerance) const noexcept
{
    if (bound.ax.size() != segmentCount())
    {
        return RayHit{RayStatus::InvalidQuery, false, Infinity, Witness2{}};
    }
    const double* const ax = bound.ax.data();
    const double* const az = bound.az.data();
    const double* const cross = bound.cross.data();
    return castFrom(bound.origin, direction, tieTolerance,
                    [=](std::size_t j) noexcept { return Relative{ax[j], az[j], cross[j]}; });
}

RayHit PreparedProfile::castRay(Point2<Length> origin, Vec2d direction, double tieTolerance) const noexcept
{
    const double ox = origin.x.si();
    const double oz = origin.y.si();
    return castFrom(origin, direction, tieTolerance, [&](std::size_t j) noexcept {
        const double ax = x_[j] - ox;
        const double az = y_[j] - oz;
        return Relative{ax, az, ax * dy_[j] - az * dx_[j]};
    });
}

CircleSweep PreparedProfile::sweepCircle(Point2<Length> start, Vec2<Length> delta, Length radiusLength,
                                         double fractionTolerance) const noexcept
{
    const double cx = start.x.si();
    const double cz = start.y.si();
    const double dx = delta.x.si();
    const double dz = delta.y.si();
    const double radius = radiusLength.si();
    CircleSweep sweep{CircleSweepStatus::InvalidQuery, 0.0, Witness2{}};
    if (!finite(cx, cz) || !finite(dx, dz) || !(radius > 0.0) || !std::isfinite(radius) || !(fractionTolerance >= 0.0)
        || !std::isfinite(fractionTolerance))
    {
        return sweep;
    }

    // Features beyond this x-window are farther than the radius from every point of the path (with room for the
    // root tolerance and for the fraction slack, which lets the centre overshoot the path ends).
    const double reach = radius * (1.0 + 1e-9) + 1e-9 * (1.0 + std::abs(dx) + std::abs(dz)) + fractionTolerance * std::abs(dx);
    const std::size_t first = lowerIndex(x_, std::min(cx, cx + dx) - reach);
    const std::size_t past = upperIndex(x_, std::max(cx, cx + dx) + reach);
    const double upper = 1.0 + fractionTolerance;

    double hit = Infinity;
    WitnessKind hitKind = WitnessKind::Vertex;
    std::size_t hitIndex = 0;
    double hitSign = 0.0;
    const auto candidate = [&](double t, WitnessKind kind, std::size_t index, double sign) noexcept {
        if (t >= -fractionTolerance && t <= upper)
        {
            const double clamped = std::max(0.0, std::min(1.0, t));
            if (clamped < hit)
            {
                hit = clamped;
                hitKind = kind;
                hitIndex = index;
                hitSign = sign;
            }
        }
    };

    // Vertex caps: roots t of |(x, z) + t (dx, dz)|^2 = r^2. The q form avoids cancellation.
    const double a = dx * dx + dz * dz;
    if (a != 0.0)
    {
        const double fourA = 4.0 * a;
        const double radiusSquared = radius * radius;
        constexpr double RootSlack = -32.0 * std::numeric_limits<double>::epsilon();
        for (std::size_t v = first; v < past; ++v)
        {
            const double x = cx - x_[v];
            const double z = cz - y_[v];
            const double b = 2.0 * (x * dx + z * dz);
            const double c = x * x + z * z - radiusSquared;
            const double discriminant = b * b - fourA * c;
            const double scale = b * b + std::abs(fourA * c);
            if (discriminant < RootSlack * scale)
            {
                continue;
            }
            const double root = std::sqrt(std::max(0.0, discriminant));
            const double q = -0.5 * (b + (b >= 0.0 ? root : -root));
            if (q == 0.0)
            {
                candidate(-b / (2.0 * a), WitnessKind::Vertex, v, 0.0);
            }
            else
            {
                candidate(q / a, WitnessKind::Vertex, v, 0.0);
                candidate(c / q, WitnessKind::Vertex, v, 0.0);
            }
        }
    }

    // Segment faces: the centre is `radius` from the segment's line on either side, with its foot on the segment.
    const std::size_t firstSegment = first == 0U ? 0U : first - 1U;
    const std::size_t endSegment = std::min(past, segmentCount());
    for (std::size_t k = firstSegment; k < endSegment; ++k)
    {
        const double nx = nx_[k];
        const double nz = ny_[k];
        const double dv = dx * nx + dz * nz;
        if (dv == 0.0)
        {
            continue;
        }
        const double ax = x_[k];
        const double az = y_[k];
        const double d0 = (cx - ax) * nx + (cz - az) * nz;
        for (const double sign : {-1.0, 1.0})
        {
            const double t = (sign * radius - d0) / dv;
            if (t < -fractionTolerance || t > upper)
            {
                continue;
            }
            const double projection = ((cx + t * dx - ax) * dx_[k] + (cz + t * dz - az) * dy_[k]) / normSquared_[k];
            if (projection >= 0.0 && projection <= 1.0)
            {
                candidate(t, WitnessKind::Edge, k, sign);
            }
        }
    }

    if (!(hit < Infinity))
    {
        sweep.status = CircleSweepStatus::Clear;
        sweep.fraction = 1.0;
        return sweep;
    }
    sweep.status = CircleSweepStatus::Contact;
    sweep.fraction = hit;
    const double centreX = cx + hit * dx;
    const double centreZ = cz + hit * dz;
    if (hitKind == WitnessKind::Vertex)
    {
        const double ex = centreX - x_[hitIndex];
        const double ez = centreZ - y_[hitIndex];
        const double length = std::sqrt(ex * ex + ez * ez);
        sweep.witness = Witness2{{Length::fromSi(x_[hitIndex]), Length::fromSi(y_[hitIndex])},
                                 length > 0.0 ? Vec2d{ex / length, ez / length} : Vec2d{0.0, 1.0},
                                 WitnessKind::Vertex,
                                 hitIndex};
    }
    else
    {
        const double nx = hitSign * nx_[hitIndex];
        const double nz = hitSign * ny_[hitIndex];
        sweep.witness = Witness2{{Length::fromSi(centreX - radius * nx), Length::fromSi(centreZ - radius * nz)},
                                 Vec2d{nx, nz},
                                 WitnessKind::Edge,
                                 hitIndex};
    }
    return sweep;
}

bool PreparedProfile::onSurfaceWithin(double px, double pz, double tolerance, double slack) const noexcept
{
    // Segments whose x-range, widened by `slack`, reaches px: those from the vertex before the first at or beyond
    // px - slack up to the last vertex at or before px + slack.
    const std::size_t lo = lowerIndex(x_, px - slack);
    const std::size_t hi = upperIndex(x_, px + slack);
    const std::size_t end = std::min(hi, segmentCount());
    for (std::size_t k = lo == 0U ? 0U : lo - 1U; k < end; ++k)
    {
        const double dx = dx_[k];
        const double dz = dy_[k];
        const double t = ((px - x_[k]) * dx + (pz - y_[k]) * dz) / sumSquares_[k];
        if (t >= -tolerance && t <= 1.0 + tolerance
            && std::abs((px - x_[k]) * dz - (pz - y_[k]) * dx) / norm_[k] < tolerance)
        {
            return true;
        }
    }
    return false;
}

// x reach of a point within `tolerance` of a segment, measured from the segment's x-range: the overshoot along the
// segment is a fraction of its length, so the longest segment bounds it. Padded so rounding cannot shrink it.
double PreparedProfile::reachFor(double tolerance) const noexcept
{
    return tolerance * (longestSegment_ * (1.0 + 1e-6) + 1.0) + 1e-12;
}

bool PreparedProfile::onSurface(Point2<Length> point, double tolerance) const noexcept
{
    const double px = point.x.si();
    const double pz = point.y.si();
    if (!finite(px, pz) || !(tolerance >= 0.0) || !std::isfinite(tolerance))
    {
        return false;
    }
    return onSurfaceWithin(px, pz, tolerance, reachFor(tolerance));
}

ContourIntrusion PreparedProfile::contourIntrusion(std::span<const Point2<Length>> outline, double tolerance) const noexcept
{
    ContourIntrusion result{ContourStatus::InvalidQuery, Length{}, Witness2{}};
    if (!(tolerance >= 0.0) || !std::isfinite(tolerance))
    {
        return result;
    }
    for (const Point2<Length>& p : outline)
    {
        if (!finite(p.x.si(), p.y.si()))
        {
            return result;
        }
    }
    const double slack = reachFor(tolerance);
    double maximum = 0.0;
    double deepestX = 0.0;
    for (std::size_t i = 0; i < outline.size(); ++i)
    {
        const Point2<Length> a = outline[i];
        const Point2<Length> b = outline[(i + 1U) % outline.size()];
        const double ax = a.x.si();
        const double az = a.y.si();
        const double bx = b.x.si();
        const double dx = bx - ax;
        const double dz = b.y.si() - az;
        const auto samplePart = [&](double from, double to) noexcept {
            for (const double f : {SampleEdge, 0.5, 1.0 - SampleEdge})
            {
                const double t = from + (to - from) * f;
                const double x = ax + t * dx;
                const double z = az + t * dz;
                if (!onSurfaceWithin(x, z, tolerance, slack))
                {
                    const double depth = surfaceAt(x, nullptr) - z;
                    if (depth > maximum)
                    {
                        maximum = depth;
                        deepestX = x;
                    }
                }
            }
        };
        // The parts of the chord between the terrain vertices' x, in order along the chord. The vertex x never
        // decreases, so the cuts come out sorted without a scratch array: ascending for a chord running
        // forward, descending for one running back.
        double previous = 0.0;
        if (std::abs(dx) > VerticalChord)
        {
            if (dx > 0.0)
            {
                for (std::size_t v = upperIndex(x_, ax), end = lowerIndex(x_, bx); v < end; ++v)
                {
                    const double t = (x_[v] - ax) / dx;
                    if (t > 0.0 && t < 1.0)
                    {
                        samplePart(previous, t);
                        previous = t;
                    }
                }
            }
            else
            {
                for (std::size_t v = lowerIndex(x_, ax), end = upperIndex(x_, bx); v-- > end;)
                {
                    const double t = (x_[v] - ax) / dx;
                    if (t > 0.0 && t < 1.0)
                    {
                        samplePart(previous, t);
                        previous = t;
                    }
                }
            }
        }
        samplePart(previous, 1.0);
    }
    if (!(maximum > 0.0))
    {
        result.status = ContourStatus::Clear;
        return result;
    }
    std::size_t feature = 0;
    const double surface = surfaceAt(deepestX, &feature);
    result.status = ContourStatus::Intrudes;
    result.depth = Length::fromSi(maximum);
    result.witness = Witness2{{Length::fromSi(deepestX), Length::fromSi(surface)},
                              Vec2d{nx_[feature], ny_[feature]},
                              WitnessKind::Edge,
                              feature};
    return result;
}

} // namespace sub0physics
