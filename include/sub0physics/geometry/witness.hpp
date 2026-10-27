#pragma once

#include <cstddef>
#include <sub0math/geometry/witness.hpp>
#include <sub0math/units/quantities.hpp>

namespace sub0physics
{

/** Which kind of surface feature a witness point lies on (Sub0Math's vocabulary). */
using sub0math::WitnessKind;

/**
 * The surface feature that answers a nearest-point or hit query: Sub0Math's `Witness` in metres. Queries return the
 * feature itself and never average the normals of several: at a corner the witness is the vertex with the true
 * approach direction, not a blend of its faces.
 *
 * @tparam N Dimensions of the space: 2 for a profile, 3 for a surface.
 */
template<std::size_t N>
using Witness = sub0math::Witness<sub0math::f64::Length, N>;

using Witness2 = Witness<2>; ///< A witness on a planar profile.
using Witness3 = Witness<3>; ///< A witness on a surface in space.

} // namespace sub0physics
