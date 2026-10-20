#include <compare>
#include <concepts>
#include <cstdint>
#include <doctest.h>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

#include "sub0physics/core/id.hpp"

using namespace sub0physics;

namespace
{

using BodyLike = Id<struct BodyLikeTag>;
using SiteLike = Id<struct SiteLikeTag>;
using Wide = Id<struct WideTag, std::uint64_t>;
using Small = Id<struct SmallTag, std::uint8_t>;

template<typename A, typename B>
concept Comparable = requires(A a, B b) { a == b; };

} // namespace

// Explicit construction, and no conversions to the integer or between tags.
static_assert(!std::is_convertible_v<std::uint32_t, BodyLike>);
static_assert(!std::is_convertible_v<BodyLike, std::uint32_t>);
static_assert(std::is_constructible_v<BodyLike, std::uint32_t>);
static_assert(!std::is_convertible_v<BodyLike, SiteLike>);
static_assert(!std::is_constructible_v<SiteLike, BodyLike>);
static_assert(!Comparable<BodyLike, SiteLike>);
static_assert(Comparable<BodyLike, BodyLike>);
static_assert(std::is_trivially_copyable_v<BodyLike> && sizeof(BodyLike) == sizeof(std::uint32_t));
static_assert(sizeof(Wide) == sizeof(std::uint64_t) && sizeof(Small) == 1U);

TEST_CASE("An id starts invalid, and zero is a real id")
{
    constexpr BodyLike none;
    constexpr BodyLike first{0U};
    static_assert(!none.valid());
    static_assert(first.valid());
    static_assert(none == BodyLike::invalid());
    static_assert(none.value() == UINT32_MAX);
    CHECK(none != first);
    CHECK(Small{}.value() == 255U);
    CHECK(Wide{}.value() == UINT64_MAX);
}

TEST_CASE("Ids order by value and are three-way comparable")
{
    CHECK(BodyLike{1U} < BodyLike{2U});
    CHECK(BodyLike{2U} >= BodyLike{2U});
    CHECK(std::is_gt(BodyLike{3U} <=> BodyLike{2U}));
    CHECK(BodyLike{7U} == BodyLike{7U});
    CHECK(BodyLike{5U} < BodyLike::invalid()); // the sentinel sorts last
}

TEST_CASE("next advances a counter and never wraps onto a used value")
{
    CHECK(BodyLike{4U}.next() == BodyLike{5U});
    CHECK(BodyLike{}.next() == BodyLike::invalid()); // stays invalid
    CHECK(Small{254U}.next() == Small::invalid());   // the top counter value becomes the sentinel
    CHECK(Small{253U}.next().next() == Small::invalid());
    CHECK(Small{253U}.next().next().next() == Small::invalid()); // and stays there
    BodyLike generation{0U};
    for (int i = 0; i < 1000; ++i) generation = generation.next();
    CHECK(generation.value() == 1000U);
}

TEST_CASE("Ids hash by value and key unordered containers")
{
    std::unordered_set<BodyLike> set;
    for (std::uint32_t i = 0; i < 100U; ++i) set.insert(BodyLike{i % 10U});
    CHECK(set.size() == 10U);
    CHECK(set.contains(BodyLike{9U}));
    CHECK_FALSE(set.contains(BodyLike{10U}));
    std::unordered_map<Wide, int> map{{Wide{1ULL << 40U}, 3}};
    CHECK(map.at(Wide{1ULL << 40U}) == 3);
    CHECK(std::hash<BodyLike>{}(BodyLike{42U}) == std::hash<std::uint32_t>{}(42U));
}
