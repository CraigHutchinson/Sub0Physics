#pragma once

#include <cstdint>
#include <expected>
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

/**
 * Makes the error of a call whose argument broke its contract.
 *
 * @param detail What is wrong; must be a string literal (static storage).
 * @return An unexpected `Fault{FaultKind::InvalidInput, detail}`, ready to return from a function
 *         that yields `std::expected<T, Fault>`.
 */
[[nodiscard]] constexpr std::unexpected<Fault> invalidInput(std::string_view detail) noexcept
{
    return std::unexpected(Fault{FaultKind::InvalidInput, detail});
}

} // namespace sub0physics
