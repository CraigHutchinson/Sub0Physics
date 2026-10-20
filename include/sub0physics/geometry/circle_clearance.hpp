#pragma once

#include <cstddef>
#include <cstdint>
#include <sub0math/geometry/sphere.hpp>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>
#include <span>

#include "sub0physics/geometry/polyline_profile.hpp"
#include "sub0physics/geometry/prepared_profile.hpp"
#include "sub0physics/geometry/witness.hpp"

namespace sub0physics
{

/** Whether a clearance query has an answer, and why not. */
enum class ClearanceStatus : std::uint8_t
{
    Available,     ///< Distance, gap and witnesses are valid.
    OutsideDomain, ///< The footprint leaves the profile (never extrapolated), or a query value is non-finite, beyond
                   ///< `PolylineProfile::CoordinateLimit`, or (the tolerance) negative.
    OnBoundary,    ///< The centre lies on the surface within tolerance; the normal is undefined.
    CentreInside   ///< The centre is inside the solid; buried contact is unsupported.
};

/**
 * Result of a clearance query. Fields other than `status` are meaningful only when Available. The
 * witnesses themselves are written to the caller's buffer, so a query never allocates.
 */
struct CircleClearance
{
    ClearanceStatus status;         ///< Whether the remaining fields are valid.
    sub0math::f64::Length distance; ///< Centre to nearest surface point.
    sub0math::f64::Length gap;      ///< distance - radius; negative means overlap.
    std::size_t witnessCount;       ///< Witnesses written to the buffer: the first ties, in profile order.
    bool overflowed;                ///< More features tied than the buffer holds; the rest were dropped.
};

/** Witnesses a typical query needs: one, or a few where a corner or valley ties. */
inline constexpr std::size_t DefaultWitnessCapacity = 8;

/** Query of a circle against a profile. */
struct CircleQuery
{
    sub0math::Circle<sub0math::f64::Length> circle; ///< The circle; positive radius.
    sub0math::f64::Length tieTolerance;             ///< Distances this close to the nearest count as ties; not a collision margin.
};

/**
 * Measures how far a circle is from a rigid profile and which features are nearest.
 *
 * Square edges produce vertex witnesses with the true approach normal, not the face normal, which is
 * what distinguishes square-edge from rounded contact. Equidistant features are all reported, each once,
 * in profile order, until the buffer is full; then `overflowed` is set and the remaining ties are dropped.
 *
 * @param profile The ground profile.
 * @param query The circle and tie tolerance.
 * @param witnesses The caller's buffer for the witnesses; may be empty when only the distance is wanted
 *        (`overflowed` is then set whenever a tie exists). Its contents past `witnessCount` are unchanged.
 * @return The measured clearance, or the reason none is defined.
 */
[[nodiscard]] CircleClearance measureCircleClearance(const PolylineProfile& profile, const CircleQuery& query,
                                                    std::span<Witness2> witnesses) noexcept;

/**
 * Measures a circle's clearance from a prepared profile; the same query and result as the overload on a
 * `PolylineProfile`, faster on long profiles.
 *
 * The nearest distance is found over the horizontal window the segment under the centre bounds, with four
 * independent minima; the arithmetic of every distance is the unprepared one, so the result, including the
 * distance, gap and witnesses, is bit-identical to the other overload's.
 *
 * @param profile The prepared ground profile.
 * @param query The circle and tie tolerance.
 * @param witnesses The caller's buffer for the witnesses, as the other overload.
 * @return The measured clearance, or the reason none is defined.
 */
[[nodiscard]] CircleClearance measureCircleClearance(const PreparedProfile& profile, const CircleQuery& query,
                                                    std::span<Witness2> witnesses) noexcept;

} // namespace sub0physics
