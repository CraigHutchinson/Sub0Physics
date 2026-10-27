#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <sub0math/geometry/capsule.hpp>
#include <sub0math/geometry/convex_polygon.hpp>
#include <sub0math/geometry/point.hpp>
#include <sub0math/units/quantities.hpp>
#include <vector>

#include "sub0physics/geometry/polyline_profile.hpp"
#include "sub0physics/geometry/witness.hpp"

namespace sub0physics
{

/** A casing edge: a segment between two nodes, thickened by the skin radius. */
using CasingEdge = sub0math::Capsule<sub0math::f64::Length, 2>;

/** One convex piece of a terrain: four corners, counter-clockwise, with outward edge normals. */
using TerrainPiece = sub0math::ConvexPolygon<double, 4>;

/** Whether a terrain query has an answer. */
enum class TerrainQueryStatus : std::uint8_t
{
    Available,   ///< The other fields are valid.
    InvalidQuery ///< A coordinate is non-finite, or a radius, margin, bound or horizon is negative or non-finite.
};

/**
 * One contact between a casing edge and one convex piece of the terrain.
 *
 * The witness is on the piece. Its `index` is the piece-local feature: vertices 0 (bottom left), 1 (bottom right),
 * 2 (top right: canonical vertex `segment + 1`), 3 (top left: canonical vertex `segment`); edges 0 (bottom), 1 (right
 * side), 2 (top: canonical segment `segment`), 3 (left side). Side faces are real faces: they are the terrain's
 * vertical faces where a neighbour is lower.
 */
struct TerrainContact
{
    std::size_t edge;                                ///< The caller's index of the casing edge.
    std::size_t piece;                               ///< Index of the convex piece.
    std::size_t segment;                             ///< Canonical segment (first vertex) on top of the piece.
    sub0math::f64::Length gap;                       ///< Skin surface to piece: negative when overlapping, m.
    double axisT;                                    ///< Position along the edge of the nearest (or deepest) axis point: 0 at `a`, 1 at `b`.
    bool axisCrosses;                                ///< The edge's axis itself touches or enters the piece (deep penetration).
    sub0math::Point2<sub0math::f64::Length> capsulePoint; ///< The skin point facing the piece.
    Witness2 witness;                                ///< The piece's feature; normal from the piece towards the edge.
};

/**
 * Outcome of a contact query. The contacts themselves go to the caller's buffer, in edge order then piece order.
 */
struct TerrainContactSet
{
    TerrainQueryStatus status;           ///< Whether the remaining fields are valid.
    std::size_t count;                   ///< Contacts written to the buffer.
    std::size_t found;                   ///< Pairs within the margin; more than `count` when the buffer overflowed.
    bool overflowed;                     ///< More pairs than the buffer holds: the first ones were kept.
    std::size_t pairsTested;             ///< (Edge, piece) pairs measured exactly after the window and height cull.
    sub0math::f64::Length minimumGap;    ///< Smallest gap over the pairs measured: exact when within the margin, else only known to exceed it; infinite when none was measured.
};

/**
 * A bound on how far any point of a casing can move during one step: the hub carrying the ring translates and
 * rotates, and the nodes move relative to the hub. Each value bounds the whole step, not only its end.
 */
struct StepMotionBound
{
    sub0math::Point2<sub0math::f64::Length> hubCentre; ///< Hub centre at the start of the step (the rotation's pivot).
    sub0math::f64::Length hubTravel;                   ///< Bound on `|c(t) - c(0)|` over the step, m.
    double hubRotation;                                ///< Bound on `|theta(t) - theta(0)|` over the step, rad.
    sub0math::f64::Length nodeTravel;                  ///< Bound on a node's movement in the hub frame, `|q(t) - q(0)|`, m.
};

/**
 * The motion of a casing over one step between two known configurations, along the step's straight interpolant in
 * the generalized coordinates (the hub's centre and angle and the nodes' hub-frame displacements each change linearly):
 * every point of the casing then moves at most `tau R` from where it starts and `(1 - tau) R` from where it ends, with
 * `R = |hubEnd - hubStart| + nodeTravel + lever |hubRotation|` (the rotation term uncapped, so it stays proportional).
 */
struct SweptStepMotion
{
    sub0math::Point2<sub0math::f64::Length> hubStart; ///< Hub centre at the start, m.
    sub0math::Point2<sub0math::f64::Length> hubEnd;   ///< Hub centre at the end, m.
    double hubRotation;                               ///< |theta(1) - theta(0)|, rad.
    sub0math::f64::Length nodeTravel;                 ///< Largest |w_i(1) - w_i(0)| of a node in the hub frame, m.
};

/**
 * A certified lower bound on a casing's gap to the terrain over a whole step.
 *
 * `clearance` is at most `horizon` (pieces further away are not measured). When positive, no point of the moving casing
 * can touch the terrain during the step: no contact can be missed or tunnelled through.
 */
struct ClearanceBound
{
    TerrainQueryStatus status;        ///< Whether the remaining fields are valid.
    sub0math::f64::Length clearance;  ///< Lower bound on the gap at every instant of the step, capped at `horizon`, m.
    sub0math::f64::Length reach;      ///< The largest movement of any casing point over the step, m.
    bool certified;                   ///< `clearance > 0`: the step cannot bring the casing into contact.
};

/**
 * A terrain profile split into canonical convex pieces for exact contact with casing edges.
 *
 * Built once from a validated `PolylineProfile` (the constructor allocates; nothing else does). The construction is
 * the convex prism decomposition of the reference vehicle experiments:
 *
 * - **Collinear runs are merged first.** A vertex is dropped while its two neighbouring segments point the same way
 *   (`|u x v| <= 1e-12 |u| |v|` and `u . v > 0`, the norms in the compensated arithmetic of `compensatedNorm`), so
 *   subdividing an authored line adds no support and changes no contact.
 * - **One piece per non-vertical canonical segment:** the quadrilateral from the segment down to a common floor one
 *   metre below the lowest vertex, with vertical sides. A vertical segment has no piece of its own; it is the side
 *   face of the neighbour that rises above it, so vertical faces are real faces with normal +-x.
 *
 * Contact multiplicity follows the reference: one contact per (casing edge, piece) pair whose gap is within the
 * margin, independent of the authored subdivision. Witnesses are exact features (`sub0math::contact`), never averaged
 * normals; at a convex terrain corner both pieces that meet there may report the vertex, each with the true
 * approach direction.
 */
class ConvexTerrain
{
public:
    /** Collinearity test for merging: the sine of the turn, relative to the segment lengths. */
    static constexpr double CollinearSine = 1e-12;
    /** The floor of every piece lies this far below the lowest canonical vertex, m. */
    static constexpr double FloorDepth = 1.0;

    /**
     * Canonicalises a profile and builds its pieces. Allocates.
     *
     * @param profile A profile validated by `PolylineProfile::create`.
     */
    explicit ConvexTerrain(const PolylineProfile& profile);

    /** @return The canonical vertices' x, m (never decreasing). */
    [[nodiscard]] std::span<const double> vertexX() const noexcept
    {
        return x_;
    }

    /** @return The canonical vertices' heights, m. */
    [[nodiscard]] std::span<const double> vertexY() const noexcept
    {
        return y_;
    }

    /** @return The convex pieces, in x order. */
    [[nodiscard]] std::span<const TerrainPiece> pieces() const noexcept
    {
        return pieces_;
    }

    /** @return For each piece, the canonical segment on its top (index of its first vertex). */
    [[nodiscard]] std::span<const std::size_t> pieceSegments() const noexcept
    {
        return segment_;
    }

    /** @return The common floor of the pieces, m. */
    [[nodiscard]] double floor() const noexcept
    {
        return floor_;
    }

    /**
     * The outward normal of a canonical segment (solid below): `(-dy, dx) / |d|`, the norm compensated.
     *
     * @param segment Segment index, less than `vertexX().size() - 1`.
     * @return The unit normal; +-x for a vertical segment.
     */
    [[nodiscard]] sub0math::Vec2d segmentNormal(std::size_t segment) const noexcept;

    /**
     * The admissible normals at a canonical vertex of the terrain: the cone between its two segments' normals (at the
     * ends of the profile, the end piece's vertical side). Valid only where the terrain turns down (a convex corner);
     * at a straight or concave vertex there is no exposed corner and `convex` is false.
     *
     * @param vertex Canonical vertex index.
     * @param convex Receives whether the terrain is convex at the vertex.
     * @return The cone, counter-clockwise from the following segment's normal to the preceding one's.
     */
    [[nodiscard]] sub0math::NormalCone<double> cornerCone(std::size_t vertex, bool& convex) const noexcept;

    /**
     * Finds every piece within `margin` of one casing edge and writes one contact per piece.
     *
     * Pieces are culled by their x extent (binary search on the sorted pieces) and their top height, then measured
     * exactly. Nothing is allocated.
     *
     * @param edge The edge; its radius is the skin, not negative.
     * @param margin Activation distance: pairs whose gap is at most this are reported; not negative.
     * @param out The caller's buffer.
     * @param edgeIndex Copied into each contact's `edge`.
     * @return Counts and the smallest gap measured, or InvalidQuery.
     */
    [[nodiscard]] TerrainContactSet contacts(const CasingEdge& edge, sub0math::f64::Length margin,
                                             std::span<TerrainContact> out, std::size_t edgeIndex = 0) const noexcept;

    /**
     * Finds the contacts of a closed ring of casing edges: edge `i` joins node `i` to node `i + 1`, and the last joins
     * the last node to the first.
     *
     * @param nodes Node positions, in ring order; at least two.
     * @param skin Edge radius, not negative.
     * @param margin Activation distance, not negative.
     * @param out The caller's buffer, filled in edge order then piece order.
     * @return Counts over all edges and the smallest gap measured, or InvalidQuery.
     */
    [[nodiscard]] TerrainContactSet ringContacts(std::span<const sub0math::Point2<sub0math::f64::Length>> nodes,
                                                 sub0math::f64::Length skin, sub0math::f64::Length margin,
                                                 std::span<TerrainContact> out) const noexcept;

    /**
     * Bounds a casing edge's gap from below over a whole step, including the hub's rotation.
     *
     * Every point of the edge moves at most `reach = hubTravel + nodeTravel + lever * min(hubRotation, 2)` during the
     * step, where `lever` is the larger distance of the edge's ends from the hub centre: a point at hub-frame position
     * `q` moves by `(c(t) - c(0)) + R(t) (q(t) - q(0)) + (R(t) - R(0)) q(0)`, and `|R(t) - R(0)| q(0)| = 2 |q(0)|
     * |sin(dtheta / 2)| <= |q(0)| min(|dtheta|, 2)`. Interior points of the edge are convex combinations of its ends,
     * so they move no further. The gap to a fixed set is 1-Lipschitz in the edge's movement, so the gap at every
     * instant is at least the starting gap less `reach`. No trigonometry is evaluated, so the bound is the same on
     * every platform.
     *
     * @param edge The edge at the start of the step.
     * @param motion Bounds on the step's motion; non-negative.
     * @param horizon Pieces further than `reach + horizon` are not measured, so `clearance` is capped at `horizon`.
     * @return The bound, or InvalidQuery.
     */
    [[nodiscard]] ClearanceBound certifyClearance(const CasingEdge& edge, const StepMotionBound& motion,
                                                  sub0math::f64::Length horizon) const noexcept;

    /**
     * Bounds a closed ring's gap from below over a whole step: the smallest of its edges' bounds.
     *
     * @param nodes Node positions at the start of the step, in ring order; at least two.
     * @param skin Edge radius, not negative.
     * @param motion Bounds on the step's motion.
     * @param horizon As `certifyClearance`.
     * @return The bound for the whole ring, or InvalidQuery.
     */
    [[nodiscard]] ClearanceBound certifyRingClearance(std::span<const sub0math::Point2<sub0math::f64::Length>> nodes,
                                                      sub0math::f64::Length skin, const StepMotionBound& motion,
                                                      sub0math::f64::Length horizon) const noexcept;

    /**
     * Bounds a closed ring's gap from below over a step whose start and end are both known: the swept certificate of
     * the improved contact model (docs/WHEEL.md section 3.4, "M2d status").
     *
     * For each edge the gap is 1-Lipschitz in the edge's movement, and along the step's interpolant (SweptStepMotion)
     * a point is within `tau R` of its start and `(1 - tau) R` of its end, so at every instant the gap is at least
     * `max(g0 - tau R, g1 - (1 - tau) R) >= (g0 + g1 - R) / 2`, with g0 and g1 the edge's gaps at the two ends. This
     * two-sided bound certifies far more steps than the one-sided `g0 - R` (certifyRingClearance), which knows only
     * the start: a casing rolling past a corner at millimetres of clearance is certified, and a step is only refused
     * where the casing could really have touched in between (a graze, or tunnelling through a thin feature).
     * No trigonometry is evaluated.
     *
     * @param start Node positions at the start, in ring order; at least two.
     * @param end Node positions at the end, the same count.
     * @param skin Edge radius, not negative.
     * @param motion The step's motion.
     * @param horizon Gaps beyond `R + horizon` are not measured, so each bound is capped at `horizon`.
     * @param edgeBounds Empty, or one entry per edge receiving that edge's bound (capped at `horizon`), m.
     * @return The smallest bound over the edges, or InvalidQuery.
     */
    [[nodiscard]] ClearanceBound certifySweptRingClearance(std::span<const sub0math::Point2<sub0math::f64::Length>> start,
                                                           std::span<const sub0math::Point2<sub0math::f64::Length>> end,
                                                           sub0math::f64::Length skin, const SweptStepMotion& motion,
                                                           sub0math::f64::Length horizon,
                                                           std::span<double> edgeBounds = {}) const noexcept;

private:
    /** One edge's smallest gap to the pieces, measured exactly when below `cap` and reported as `cap` otherwise. */
    [[nodiscard]] double edgeGap(double ax, double ay, double bx, double by, double radius, double cap) const noexcept;
    /** One edge's clearance bound (capped at `horizon`), its reach written to `reach`; no input checks. */
    [[nodiscard]] double edgeClearance(double ax, double ay, double bx, double by, double radius,
                                       const StepMotionBound& motion, double horizon, double& reach) const noexcept;
    /** One edge's contacts appended to `out` and `set`, the smallest gap folded into `minimum`; no input checks. */
    void edgeContacts(double ax, double ay, double bx, double by, double radius, double margin,
                      std::span<TerrainContact> out, std::size_t edgeIndex, TerrainContactSet& set,
                      double& minimum) const noexcept;
    /** Index range of pieces whose x extent meets [lo, hi]. */
    void window(double lo, double hi, std::size_t& first, std::size_t& last) const noexcept;

    std::vector<double> x_, y_;             ///< Canonical vertices, m.
    std::vector<TerrainPiece> pieces_;      ///< One per non-vertical canonical segment.
    std::vector<double> pieceRight_;        ///< Right x of each piece, m (never decreasing): the window's search key.
    std::vector<double> pieceTop_;          ///< Highest point of each piece, m.
    std::vector<std::size_t> segment_;      ///< Canonical segment on top of each piece.
    double floor_ = 0.0;                    ///< Common floor, m.
    double highest_ = -std::numeric_limits<double>::infinity(); ///< Highest point of any piece, m.
};

} // namespace sub0physics
