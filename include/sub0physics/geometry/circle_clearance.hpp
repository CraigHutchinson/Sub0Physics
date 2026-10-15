#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "sub0physics/geometry/point2.hpp"
#include "sub0physics/geometry/polyline_profile.hpp"

namespace sub0physics
{

/** Whether a clearance query has an answer, and why not. */
enum class ClearanceStatus : std::uint8_t
{
    Available,     ///< Distance, gap and witnesses are valid.
    OutsideDomain, ///< The circle's footprint leaves the profile; terrain is never extrapolated.
    OnBoundary,    ///< The centre lies on the surface within tolerance; the normal is undefined.
    CentreInside   ///< The centre is inside the solid; buried contact is unsupported.
};

/** Which profile feature a witness point lies on. */
enum class FeatureKind : std::uint8_t
{
    Vertex,  ///< A polyline vertex; `index` is the vertex index.
    Segment  ///< The interior of a segment; `index` is its first vertex.
};

/** A nearest surface point and its approach direction. */
struct ClearanceWitness
{
    Point2 point;      ///< Nearest point on the surface, m.
    Point2 normal;    ///< Unit vector from that point towards the centre.
    FeatureKind kind; ///< Feature the point lies on.
    std::size_t index; ///< Vertex or segment index.
};

/** Result of a clearance query. Fields other than `status` are meaningful only when Available. */
struct CircleClearance
{
    ClearanceStatus status;                ///< Whether the remaining fields are valid.
    double distance;                       ///< Centre to nearest surface point, m.
    double gap;                            ///< distance - radius, m; negative means overlap.
    std::vector<ClearanceWitness> witnesses; ///< Every feature within the tie tolerance of the nearest.
};

/** Query of a circle against a profile. */
struct CircleQuery
{
    Point2 centre;       ///< Circle centre, m.
    double radius;       ///< Circle radius, m; positive.
    double tieTolerance; ///< Distances closer than this to the nearest count as ties, m; not a collision margin.
};

/**
 * Measures how far a circle is from a rigid profile and which features are nearest.
 *
 * Square edges produce vertex witnesses with the true approach normal, not the face normal, which is
 * what distinguishes square-edge from rounded contact. Equidistant features are all reported.
 *
 * @param profile The ground profile.
 * @param query The circle and tie tolerance.
 * @return The measured clearance, or the reason none is defined.
 */
[[nodiscard]] CircleClearance measureCircleClearance(const PolylineProfile& profile, const CircleQuery& query);

} // namespace sub0physics
