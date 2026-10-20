#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <sub0math/geometry/point.hpp>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>
#include <vector>

#include "sub0physics/geometry/polyline_profile.hpp"
#include "sub0physics/geometry/witness.hpp"

namespace sub0physics
{

/** Outcome of a ray cast. */
enum class RayStatus : std::uint8_t
{
    Hit,         ///< The ray meets a front-facing segment; the other fields are valid.
    Miss,        ///< No front-facing segment is met; `distance` is infinite.
    InvalidQuery ///< A value is non-finite, the tolerance is negative, or the bound origin was built for a profile of another size.
};

/**
 * First hit of a ray on a profile. Fields other than `status` are meaningful only when Hit, except
 * `distance`, which is infinite for a Miss.
 */
struct RayHit
{
    RayStatus status;   ///< Whether the remaining fields are valid.
    bool ambiguous;     ///< Another segment with a distinctly different normal ties for first hit (a corner or a tangent ray).
    double distance;    ///< Along the ray, in multiples of the direction's length (metres for a unit direction).
    Witness2 witness;   ///< The segment hit: `Edge`, `index` its first vertex, `normal` its outward unit normal (solid below).
};

/**
 * The profile's segments expressed relative to one ray origin, so many rays from the same origin share the
 * subtraction and the cross product that do not depend on the direction. Reuse one object across solves: it
 * allocates only when the profile has more segments than it has held before.
 */
struct RayOrigin
{
    sub0math::Point2<sub0math::f64::Length> origin{}; ///< The origin the arrays were built for.
    std::vector<double> ax;                           ///< Segment start x minus origin x, m.
    std::vector<double> az;                           ///< Segment start y minus origin y, m.
    std::vector<double> cross;                        ///< `ax * dy - az * dx` per segment, m^2.
};

/** Outcome of a circle sweep. */
enum class CircleSweepStatus : std::uint8_t
{
    Contact,     ///< The circle touches the profile at `fraction`.
    Clear,       ///< No contact along the whole path; `fraction` is 1.
    InvalidQuery ///< A value is non-finite, the radius is not positive or the tolerance is negative.
};

/** First touch of a circle swept along a straight path. */
struct CircleSweep
{
    CircleSweepStatus status; ///< Whether the remaining fields are valid.
    double fraction;          ///< Path fraction of the first touch, 0 to 1; not time. For Clear: 1.
    Witness2 witness;         ///< Touching feature (Vertex, or Edge for a segment face); `normal` points to the circle centre.
};

/** Outcome of a contour check. */
enum class ContourStatus : std::uint8_t
{
    Clear,       ///< The contour does not enter the solid; `depth` is zero.
    Intrudes,    ///< The contour cuts into the solid; `depth` and `witness` are valid.
    InvalidQuery ///< A coordinate is non-finite or the tolerance is negative.
};

/** Deepest intrusion of a closed contour into the solid below a profile. */
struct ContourIntrusion
{
    ContourStatus status;     ///< Whether the remaining fields are valid.
    sub0math::f64::Length depth; ///< Largest depth below the surface at a sampled contour point, m; zero when Clear.
    Witness2 witness;         ///< Surface feature straight above the deepest sample (Edge, or Vertex at a vertex x).
};

/**
 * A profile prepared for repeated geometric queries.
 *
 * Built once from a validated `PolylineProfile`, it keeps structure-of-arrays vertices and, per segment,
 * the extents, squared lengths and the outward unit normal (solid below). Every query is
 * `noexcept` and allocates nothing; queries that report several features write them to caller-owned storage.
 * The vertex x never decreases, so queries narrow their work to the horizontal window that can matter by
 * binary search instead of scanning every segment.
 *
 * Exactness. Each query uses the same operand order as the straightforward formulation it replaces
 * (origin subtracted first, `x * dy - z * dx` cross products, divisions where the formula divides), so results
 * are reproducible bit for bit across compilers when floating-point contraction is off, which this library
 * enforces. Unit normals use a scaled, compensated norm (the arithmetic of V8's `Math.hypot`, which the
 * ported goldens rely on), not `std::hypot`, so they do not depend on the C library. Projection parameters are
 * divided by the squared length, not multiplied by its reciprocal: the reciprocal measured about 10% faster in
 * the clearance query but changes the last bit of distances and parameters, and with them witness classification
 * at a vertex, so the prepared and unprepared clearance queries would no longer agree bit for bit.
 *
 * The object copies the profile, so it does not refer to it afterwards.
 */
class PreparedProfile
{
public:
    /** Ray and contour rule: a segment is front-facing when `normal . direction` is below minus this. */
    static constexpr double GrazingIncidence = 1e-12;
    /** Ray rule: tied hits are ambiguous when their normals differ by more than this (unit vectors). */
    static constexpr double DistinctNormal = 1e-6;
    /** Contour rule: chords with a shorter x extent are not split at terrain vertices, m. */
    static constexpr double VerticalChord = 1e-14;
    /** Contour rule: each part of a chord is sampled at these fractions, which catch vertical faces from either side. */
    static constexpr double SampleEdge = 1e-8;

    /**
     * Prepares a profile. Allocates the arrays.
     *
     * @param profile A profile validated by `PolylineProfile::create`.
     */
    explicit PreparedProfile(const PolylineProfile& profile);

    /** @return Vertices; segments number one fewer. */
    [[nodiscard]] std::size_t vertexCount() const noexcept
    {
        return x_.size();
    }

    /** @return Number of segments. */
    [[nodiscard]] std::size_t segmentCount() const noexcept
    {
        return dx_.size();
    }

    /** @return Smallest x covered. */
    [[nodiscard]] sub0math::f64::Length minX() const noexcept
    {
        return sub0math::f64::Length::fromSi(x_.front());
    }

    /** @return Largest x covered. */
    [[nodiscard]] sub0math::f64::Length maxX() const noexcept
    {
        return sub0math::f64::Length::fromSi(x_.back());
    }

    /** @return Vertex x in order, m (never decreasing). */
    [[nodiscard]] std::span<const double> vertexX() const noexcept
    {
        return x_;
    }

    /** @return Vertex height (profile y) in order, m. */
    [[nodiscard]] std::span<const double> vertexY() const noexcept
    {
        return y_;
    }

    /**
     * The top of the solid at one horizontal position; the same arithmetic as `PolylineProfile::surfaceHeight`.
     *
     * @param x Horizontal position.
     * @return Surface height; negative infinity outside the profile or for a NaN x.
     */
    [[nodiscard]] sub0math::f64::Length surfaceHeight(sub0math::f64::Length x) const noexcept;

    /**
     * Expresses the segments relative to an origin for repeated `castRay` calls.
     *
     * @param origin Where the rays start.
     * @param bound Receives the arrays; its previous contents are replaced and its capacity is reused.
     */
    void bindRayOrigin(sub0math::Point2<sub0math::f64::Length> origin, RayOrigin& bound) const;

    /**
     * Casts a ray from a bound origin and reports the nearest front-facing segment.
     *
     * Segments facing away (`normal . direction >= -GrazingIncidence`) are skipped, the hit parameter along the
     * segment may overshoot its ends by `tieTolerance`, and a later hit within `tieTolerance` of the best replaces
     * nothing but marks the result ambiguous when its normal differs by more than `DistinctNormal`. The scan is
     * over every segment, in profile order, so ties resolve to the first.
     *
     * @param bound The origin, from `bindRayOrigin` on this profile (a different segment count is an InvalidQuery;
     *        binding to one profile and casting on another of equal size is the caller's error).
     * @param direction Ray direction; unit length for the incidence and ambiguity tolerances to mean what they say.
     * @param tieTolerance Distance ties and segment-parameter overshoot (a fraction of the segment); not negative.
     * @return The hit, a Miss, or InvalidQuery.
     */
    [[nodiscard]] RayHit castRay(const RayOrigin& bound, sub0math::Vec2d direction, double tieTolerance) const noexcept;

    /**
     * Casts one ray without a bound origin; the same result as binding the origin and casting, for a single ray.
     *
     * @param origin Where the ray starts.
     * @param direction Ray direction, as above.
     * @param tieTolerance As above.
     * @return The hit, a Miss, or InvalidQuery.
     */
    [[nodiscard]] RayHit castRay(sub0math::Point2<sub0math::f64::Length> origin, sub0math::Vec2d direction,
                                 double tieTolerance) const noexcept;

    /**
     * Finds the first touch of a circle moved along a straight path, independently of any sampling.
     *
     * Solves a quadratic for each vertex cap and a linear equation for each segment face, over the vertices and
     * segments within the horizontal window the circle can reach; no other feature can touch. The earliest root
     * in `[-fractionTolerance, 1 + fractionTolerance]` wins, clamped to `[0, 1]`; the first found wins a tie.
     * A circle already overlapping the profile at the start is reported by whichever roots exist; check the
     * starting clearance first.
     *
     * @param start Circle centre at the start.
     * @param delta Centre displacement over the path (not a velocity).
     * @param radius Circle radius; positive.
     * @param fractionTolerance Path-fraction slack on the root range; not negative.
     * @return Contact with its fraction and witness, Clear, or InvalidQuery.
     */
    [[nodiscard]] CircleSweep sweepCircle(sub0math::Point2<sub0math::f64::Length> start,
                                          sub0math::Vec2<sub0math::f64::Length> delta,
                                          sub0math::f64::Length radius,
                                          double fractionTolerance = 1e-10) const noexcept;

    /**
     * Tests whether a point lies on the profile.
     *
     * True when some segment has the point's projection within `tolerance` (a fraction of the segment) of its
     * ends and the point within `tolerance` metres of its line.
     *
     * @param point The point.
     * @param tolerance Both slacks, as above; not negative.
     * @return True on the surface; false for a non-finite point or tolerance.
     */
    [[nodiscard]] bool onSurface(sub0math::Point2<sub0math::f64::Length> point, double tolerance) const noexcept;

    /**
     * Finds how deeply a closed contour cuts into the solid below the profile.
     *
     * Each chord (the last point joins the first) is split at the x of the vertices it spans, and both ends and
     * the middle of every part are tested (as `SampleEdge` fractions, one-sided so vertical faces are caught).
     * A sample on the surface (`onSurface` with `tolerance`) is not an intrusion; the others contribute their
     * depth below `surfaceHeight`. No full-profile scan is made per sample, and nothing is allocated.
     *
     * @param outline Contour points, in order; any number, including none.
     * @param tolerance On-surface tolerance, as `onSurface`.
     * @return Clear, Intrudes with the deepest sample and the surface feature above it, or InvalidQuery.
     */
    [[nodiscard]] ContourIntrusion contourIntrusion(std::span<const sub0math::Point2<sub0math::f64::Length>> outline,
                                                    double tolerance = 1e-10) const noexcept;

    /** @return Segment extents `vertex[i + 1] - vertex[i]` in x, m. */
    [[nodiscard]] std::span<const double> segmentDx() const noexcept
    {
        return dx_;
    }

    /** @return Segment extents in y, m. */
    [[nodiscard]] std::span<const double> segmentDy() const noexcept
    {
        return dy_;
    }

    /** @return `dx * dx + dy * dy` per segment, m^2. */
    [[nodiscard]] std::span<const double> segmentSumSquares() const noexcept
    {
        return sumSquares_;
    }

    /** @return Outward unit normal x per segment (solid below), computed with the scaled norm. */
    [[nodiscard]] std::span<const double> segmentNormalX() const noexcept
    {
        return nx_;
    }

    /** @return Outward unit normal y per segment. */
    [[nodiscard]] std::span<const double> segmentNormalY() const noexcept
    {
        return ny_;
    }

private:
    /** The three origin-relative values of one segment. */
    struct Relative
    {
        double ax, az, cross;
    };

    template<typename RelativeAt>
    [[nodiscard]] RayHit castFrom(sub0math::Point2<sub0math::f64::Length> origin, sub0math::Vec2d direction,
                                  double tieTolerance, const RelativeAt& relativeAt) const noexcept;
    [[nodiscard]] double reachFor(double tolerance) const noexcept;
    [[nodiscard]] bool onSurfaceWithin(double px, double pz, double tolerance, double slack) const noexcept;
    [[nodiscard]] double surfaceAt(double x, std::size_t* feature) const noexcept;

    std::vector<double> x_, y_;              ///< Vertices, m.
    std::vector<double> dx_, dy_;            ///< Segment extents `vertex[i + 1] - vertex[i]`, m.
    std::vector<double> sumSquares_;         ///< `dx * dx + dy * dy`, m^2.
    std::vector<double> norm_;               ///< Segment length by the scaled, compensated norm, m.
    std::vector<double> normSquared_;        ///< `norm * norm`, m^2.
    std::vector<double> nx_, ny_;            ///< Outward unit normal `(-dy, dx) / norm` (solid below).
    double longestSegment_ = 0.0;            ///< Largest `norm`, m: bounds the x reach of a parameter overshoot.
};

} // namespace sub0physics
