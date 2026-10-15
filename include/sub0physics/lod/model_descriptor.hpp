#pragma once

#include <cstdint>
#include <string_view>

namespace sub0physics
{

/**
 * Names one model implementation so that runs, records and switches can say exactly what ran.
 *
 * Fidelity is a property of a named model and its validity envelope, not a global integer that means
 * "more accurate everywhere".
 */
struct ModelDescriptor
{
    std::string_view id;  ///< Stable identifier, e.g. "axial-closed-form"; static storage.
    std::uint32_t revision; ///< Bumped whenever results can change.
    double relativeCost;  ///< Cost per step relative to the cheapest sibling; positive.
};

} // namespace sub0physics
