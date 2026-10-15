#pragma once

#include <cstdint>
#include <string_view>

namespace sub0physics
{

/** Classifies why a physical computation refused to produce a result. */
enum class FaultKind : std::uint8_t
{
    NonFiniteState, ///< A state or derivative became infinite or NaN.
    OutsideDomain,  ///< The state left the region where the model is defined.
    InvalidInput    ///< An argument was outside the contract of the call.
};

/**
 * A refused computation: its class and a fixed explanation.
 *
 * @note `detail` always refers to a string literal, so a fault can be created, copied and returned on a
 *       hot path without allocating.
 */
struct Fault
{
    FaultKind kind;          ///< What class of refusal this is.
    std::string_view detail; ///< Human-readable reason; static storage.
};

} // namespace sub0physics
