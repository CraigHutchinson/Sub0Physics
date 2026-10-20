#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <sub0math/geometry/aabb.hpp>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/core/id.hpp"

namespace sub0physics::terrain
{

struct FeatureTag;

/** Identifies an explicit feature; unique within an asset. */
using FeatureId = Id<FeatureTag>;

/**
 * Exact geometry that a height field cannot hold: a box in the terrain frame. A face of the terrain that is
 * vertical (a riser, a block's side) is a box of zero thickness along its normal, so a height field's ramp
 * across that cell never stands in for it. The extract lists the features whose plan bounds meet it.
 */
struct BoxFeature
{
    FeatureId id;                                      ///< Identity.
    sub0math::Aabb<sub0math::f64::Length, 3> bounds;   ///< Terrain-frame bounds [m]; a face has equal bounds on one axis.
    sub0math::Vec3d faceNormal;                        ///< Outward unit normal of the face that motivated the box; zero for a solid box.
};

/** The result of `featuresMeeting`. */
struct FeatureHits
{
    std::size_t count;  ///< Indices written.
    bool overflowed;    ///< More features met the box than `out` holds.
};

/**
 * Finds the features that meet a closed plan box.
 *
 * @param features The candidates.
 * @param plan The box in the terrain frame.
 * @param out Where to write the indices into `features`.
 * @return The indices written and whether the buffer was too small.
 */
[[nodiscard]] inline FeatureHits featuresMeeting(std::span<const BoxFeature> features,
                                                 const sub0math::Aabb<sub0math::f64::Length, 2>& plan,
                                                 std::span<std::size_t> out) noexcept
{
    FeatureHits hits{0, false};
    for (std::size_t i = 0; i < features.size(); ++i)
    {
        const auto& b = features[i].bounds;
        if (b.lo.x <= plan.hi.x && b.hi.x >= plan.lo.x && b.lo.y <= plan.hi.y && b.hi.y >= plan.lo.y)
        {
            if (hits.count < out.size())
            {
                out[hits.count++] = i;
            }
            else
            {
                hits.overflowed = true;
            }
        }
    }
    return hits;
}

} // namespace sub0physics::terrain
