#pragma once

#include <compare>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>

namespace sub0physics
{

/**
 * A strongly typed identifier: an unsigned integer that only compares with ids of the same tag.
 *
 * Ids with different tags do not convert to each other or to the integer, and constructing one from an
 * integer is explicit. The largest value of `Rep` is the invalid sentinel, and a default-constructed id is
 * invalid, so an id that was never assigned cannot be mistaken for the first real one (which is 0).
 *
 * @tparam Tag An (incomplete) type that names the kind of thing identified.
 * @tparam Rep The unsigned integer that stores it.
 */
template<typename Tag, std::unsigned_integral Rep = std::uint32_t>
class Id
{
public:
    using Representation = Rep; ///< The stored integer type.

    /** Makes the invalid id. */
    constexpr Id() noexcept = default;

    /**
     * Makes an id from its integer.
     *
     * @param value The integer; the maximum of `Rep` is the invalid sentinel.
     */
    constexpr explicit Id(Rep value) noexcept : value_(value)
    {
    }

    /** @return The invalid id: the maximum of `Rep`. */
    [[nodiscard]] static constexpr Id invalid() noexcept
    {
        return Id{};
    }

    /** @return The stored integer. */
    [[nodiscard]] constexpr Rep value() const noexcept
    {
        return value_;
    }

    /** @return True unless this is the invalid sentinel. */
    [[nodiscard]] constexpr bool valid() const noexcept
    {
        return value_ != Sentinel;
    }

    /**
     * Finds the following id, for counters such as a generation or a revision.
     *
     * @return The id one higher. A counter that reaches the top becomes invalid rather than wrapping to a
     *         value already used, and the invalid id stays invalid.
     */
    [[nodiscard]] constexpr Id next() const noexcept
    {
        return value_ == Sentinel ? *this : Id{static_cast<Rep>(value_ + 1U)};
    }

    friend constexpr bool operator==(const Id&, const Id&) = default;
    friend constexpr auto operator<=>(const Id&, const Id&) = default;

private:
    static constexpr Rep Sentinel = std::numeric_limits<Rep>::max();

    Rep value_ = Sentinel;
};

} // namespace sub0physics

/** Hashes an id by its integer, so ids of any tag can key unordered containers. */
template<typename Tag, typename Rep>
struct std::hash<sub0physics::Id<Tag, Rep>>
{
    [[nodiscard]] std::size_t operator()(const sub0physics::Id<Tag, Rep>& id) const noexcept
    {
        return std::hash<Rep>{}(id.value());
    }
};
