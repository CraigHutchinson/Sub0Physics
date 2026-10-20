#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>

#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/core/id.hpp"
#include "sub0physics/multibody/description.hpp"

// Constraints (docs/MULTIBODY.md section 6, M1c): loop closures between two sites and strict unilateral limits on a
// coordinate or a distance. They are described here, compiled and validated by ModelBuilder, and honoured by the
// steppers without allocating: a velocity-level solve with the step's own factorised matrix (`J A^-1 J^T`, a few
// rows), limits as a small LCP solved by direct pivoting, steps split at limit engagement, and a projection onto
// the closure manifold after every accepted step.

namespace sub0physics::multibody
{

/** Identifies a constraint within one ModelBuilder and the Model it builds; numbered as added. */
using ConstraintId = Id<struct ConstraintIdTag, std::uint32_t>;

/** What a loop closure makes equal. */
enum class ClosureKind : std::uint8_t
{
    Point, ///< The two sites' origins coincide: three rows (world x, y, z).
    Frame  ///< The two site frames coincide, origin and orientation: six rows.
};

/**
 * A loop closure: two sites, usually on different branches of the tree, held together. A planar mechanism's
 * closure has dependent rows (the out-of-plane ones); the solver detects and drops them, so no axis mask is needed.
 */
struct LoopClosureDescription
{
    std::string_view name;                ///< Unique among constraints, non-empty; copied by the builder.
    SiteId first;                         ///< One site.
    SiteId second;                        ///< The other; on a different body (or the world).
    ClosureKind kind = ClosureKind::Point; ///< Point or frame equality.
};

/**
 * A strict unilateral limit on one joint coordinate: a travel end or a hard stop, `lower <= q <= upper`.
 * Progressive buffers stay force elements; this is the stop behind them, never overshot.
 */
struct CoordinateLimitDescription
{
    std::string_view name;            ///< Unique among constraints, non-empty; copied by the builder.
    BodyId body;                      ///< The body whose joint coordinate is limited.
    std::optional<double> lower{};    ///< Lowest coordinate (rad or m); none for no lower stop.
    std::optional<double> upper{};    ///< Highest coordinate (rad or m), above `lower`; none for no upper stop.
    double restitution = 0.0;         ///< Newton's coefficient at an impact, in [0, 1]; 0 is a plastic stop.
};

/** A strict unilateral limit on the distance between two sites: a rope (maximum) or a strut stop (minimum). */
struct DistanceLimitDescription
{
    std::string_view name;                          ///< Unique among constraints, non-empty; copied by the builder.
    SiteId first;                                   ///< One end.
    SiteId second;                                  ///< The other; on a different body (or the world).
    std::optional<sub0math::f64::Length> minimum{}; ///< Shortest distance; positive. None for no minimum.
    std::optional<sub0math::f64::Length> maximum{}; ///< Longest distance; above `minimum`. None for no maximum.
    double restitution = 0.0;                       ///< Newton's coefficient at an impact, in [0, 1].
};

/**
 * The constraint solver's levers (docs/LOD.md): targets that the solver meets, not fixed effort. Every decision
 * depends only on the state, these settings and the tick, so a replay reproduces it.
 *
 * Lengths are metres, and the same numbers are radians for angular rows (frame closures, revolute limits).
 */
struct ConstraintSettings
{
    double closureTolerance = 1e-10;  ///< Requested closure residual after projection, max |phi| [m or rad]; positive.
    std::uint32_t maxProjectionIterations = 8; ///< Guard on the projection's Gauss-Newton iterations; a step that cannot meet `closureTolerance` within it is rejected. 0 disables projection (velocity level only; the residual then drifts and is only recorded).
    double eventTolerance = 1e-9;     ///< Requested gap at a limit's engagement [m or rad]: a split step lands with the gap in [0, eventTolerance]; also the distance at which a limit counts as engaged. Positive.
    std::uint32_t maxEventIterations = 60; ///< Guard on the root find of one engagement; a step that cannot meet `eventTolerance` within it is rejected.
    std::uint32_t maxEventsPerStep = 8;    ///< Guard on the number of splits in one tick; more is a fault.
    double solverTolerance = 1e-12;   ///< Relative tolerance of the small dense solves: pivots below it times the largest are dependent rows (dropped), and the LCP's sign tests use it. Positive, below 1e-3.
    bool measureProjectionEnergy = true; ///< Whether each step measures the energy the projection changes (two energy evaluations); State::projectionEnergy stays zero when off.
};

/** A compiled loop closure. Rows `firstRow .. firstRow + rows` of the model's constraint layout are its own. */
struct LoopClosure
{
    ConstraintId id;       ///< The constraint's id.
    SiteId first;          ///< First site; the row is `first - second`.
    SiteId second;         ///< Second site.
    ClosureKind kind;      ///< Point (3 rows) or frame (6 rows).
    std::size_t firstRow;  ///< First row in the layout.
    std::size_t wrenchSlot; ///< Index of the first site's wrench in Workspace::constraintSiteWrenches; the second follows.

    /** @return The number of rows: 3 or 6. */
    [[nodiscard]] std::size_t rows() const noexcept
    {
        return kind == ClosureKind::Point ? 3U : 6U;
    }
};

/** Whether a compiled limit is on a coordinate or on a distance. */
enum class LimitKind : std::uint8_t
{
    Coordinate, ///< On a joint coordinate.
    Distance    ///< On the distance between two sites.
};

/**
 * A compiled one-sided limit, `gap = sign (x - bound) >= 0`, where x is the coordinate or the distance. A
 * description with both bounds compiles to two of these.
 */
struct UnilateralLimit
{
    ConstraintId id;          ///< The constraint's id.
    LimitKind kind;           ///< Coordinate or distance.
    std::size_t coordinate;   ///< The coordinate (LimitKind::Coordinate).
    SiteId first;             ///< First site (LimitKind::Distance).
    SiteId second;            ///< Second site (LimitKind::Distance).
    double bound;             ///< The bound [m or rad].
    double sign;              ///< +1 for a lower bound, -1 for an upper one.
    double restitution;       ///< Newton's coefficient.
    std::size_t row;          ///< Its row in the layout (after every closure row).
    std::size_t wrenchSlot;   ///< Index of the first site's wrench (LimitKind::Distance); the second follows.
};

/** A constraint's wrench on the body carrying one of its sites; plain SI in world axes, like ElementSiteWrench. */
struct ConstraintSiteWrench
{
    ConstraintId constraint; ///< The constraint.
    SiteId site;             ///< The site; the wrench acts on the site's body.
    sub0math::Vec3d force;   ///< Force on the body [N], world axes, through the site's origin.
    sub0math::Vec3d moment;  ///< Moment [N m], world axes (frame closures only).
};

/**
 * The wrench a joint transmits from the parent to its body, from computeReactions: plain SI, world axes, the
 * moment taken about the joint frame's origin. It includes the joint's own generalized force (its projection on
 * the joint's motion axis is the applied, coordinate-element and coordinate-limit force).
 */
struct JointReaction
{
    sub0math::Vec3d force;  ///< Force on the body from its parent [N], world axes.
    sub0math::Vec3d moment; ///< Moment on the body from its parent about the joint origin [N m], world axes.
};

/**
 * What the constraint solver achieved in the last step (or projection), against what ConstraintSettings
 * requested: the evidence a target-driven lever owes (docs/LOD.md).
 */
struct ConstraintReport
{
    double closureResidual = 0.0;         ///< Achieved max |phi| of the closures at the accepted state [m or rad].
    double closureVelocityResidual = 0.0; ///< Achieved max |J qd| of the closures at the accepted state [m/s or rad/s].
    double minimumGap = 0.0;              ///< Smallest limit gap at the accepted state; never negative for coordinate limits [m or rad].
    double engagementGap = 0.0;           ///< Gap at the last located engagement in the step (requested: within eventTolerance).
    std::uint32_t substeps = 0;           ///< Sub-steps the tick was split into (1 without an event).
    std::uint32_t eventIterations = 0;    ///< Root-find trials spent locating engagements in the step.
    std::uint32_t projectionIterations = 0; ///< Projection Gauss-Newton iterations in the step.
    std::uint32_t lcpPivots = 0;          ///< Active-set changes the limit LCP made in the step.
    double closureWork = 0.0;             ///< Work of the closure impulses in the step [J]; zero for an exact ideal constraint.
    double limitWork = 0.0;               ///< Work of the limit impulses in the step [J]; never positive (impact losses).
    double projectionEnergy = 0.0;        ///< Energy the projection changed in the step [J] (when measured).
};

} // namespace sub0physics::multibody
