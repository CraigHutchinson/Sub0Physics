#pragma once

#include <cmath>
#include <compare>
#include <cstdint>
#include <expected>
#include <limits>
#include <sub0math/geometry/point.hpp>
#include <sub0math/role/role.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/fault.hpp"

namespace sub0physics::terrain
{

struct TerrainFrameTag;

/**
 * A position in the terrain frame T: x forward, y left, z up, metres, double precision so that a large
 * world keeps sub-micrometre positions. Distinct from a page-local point, so the two cannot be mixed.
 */
using TerrainPoint = sub0math::Role<TerrainFrameTag, sub0math::Point3<sub0math::f64::Length>>;

/** A position in a page's local frame: float metres from the page's corner and height origin. */
using LocalPoint = sub0math::Point3<sub0math::Length>;

/** Identifies a page: its index on the grid of its level, and the level (0 = finest). */
struct TileKey
{
    std::int32_t i;     ///< Index along x.
    std::int32_t j;     ///< Index along y.
    std::uint8_t level; ///< Quadtree level; spacing doubles per level.

    friend constexpr bool operator==(const TileKey&, const TileKey&) = default;
    friend constexpr auto operator<=>(const TileKey&, const TileKey&) = default;
};

/**
 * The outcome of locating a page: the key, or the fact that the point has none. `known` is false for a
 * non-finite point, a level outside the grid, or a page index beyond 32-bit range; `key` is meaningless
 * then and must not be read, as with dereferencing an empty optional.
 */
struct PageLookup
{
    TileKey key;  ///< The page, when known.
    bool known;   ///< The point lies on a page of the grid.

    /** @return `known`, so a lookup can be tested like an optional. */
    [[nodiscard]] constexpr explicit operator bool() const noexcept
    {
        return known;
    }
};

/**
 * How terrain pages tile the plane. Every page has the same number of cells; a page at level k has
 * node spacing `finestSpacing * 2^k`, so its world size doubles per level.
 *
 * Built only by `create`, which validates the levels and the spacing, so every grid in existence is
 * valid and the `1 << level` arithmetic cannot overflow.
 */
class GridSpec
{
public:
    static constexpr std::int32_t CellsPerSide = 256;   ///< Cells along each side of every page.
    static constexpr std::uint8_t MaxLevels = 16;        ///< Most levels a grid may have.
    static constexpr double MinSpacing = 1e-6;           ///< Smallest finest node spacing [m]; keeps the float reciprocal finite.
    static constexpr double MaxSpacing = 1e6;            ///< Largest finest node spacing [m]; keeps page sides far from overflow.

    GridSpec() = delete;

    /**
     * Validates and builds a grid.
     *
     * @param finestSpacing Node spacing at level 0 [m]; finite and within [MinSpacing, MaxSpacing].
     * @param levels Number of levels, 1 to MaxLevels.
     * @return The grid, or an `InvalidInput` fault.
     */
    [[nodiscard]] static constexpr std::expected<GridSpec, Fault> create(sub0math::f64::Length finestSpacing,
                                                                         std::uint8_t levels) noexcept
    {
        if (levels < 1U || levels > MaxLevels)
        {
            return invalidInput("grid needs 1 to 16 levels");
        }
        const double s = finestSpacing.si();
        if (!(s >= MinSpacing && s <= MaxSpacing)) // NaN fails both
        {
            return invalidInput("grid spacing must be finite, from 1 micrometre to 1000 km");
        }
        return GridSpec{finestSpacing, levels};
    }

    /** @return Node spacing at level 0. */
    [[nodiscard]] constexpr sub0math::f64::Length finestSpacing() const noexcept
    {
        return finestSpacing_;
    }

    /** @return The number of levels, 1 to MaxLevels. */
    [[nodiscard]] constexpr std::uint8_t levels() const noexcept
    {
        return levels_;
    }

    /**
     * Node spacing at a level.
     *
     * @param level Below `levels()`; larger values are a contract violation.
     * @return finestSpacing * 2^level.
     */
    [[nodiscard]] constexpr sub0math::f64::Length spacing(std::uint8_t level) const noexcept
    {
        return finestSpacing_ * static_cast<double>(1U << level);
    }

    /**
     * Side length of a page at a level.
     *
     * @param level Below `levels()`.
     * @return spacing(level) * CellsPerSide.
     */
    [[nodiscard]] constexpr sub0math::f64::Length pageSide(std::uint8_t level) const noexcept
    {
        return spacing(level) * static_cast<double>(CellsPerSide);
    }

    /**
     * Finds the page containing a point at a level. Points on a shared edge belong to the page on the
     * positive side.
     *
     * @param p The point; only x and y are used.
     * @param level The level to look at.
     * @return The page's key, or not known for a level outside the grid, a non-finite coordinate or a
     *         page index that does not fit 32 bits.
     */
    [[nodiscard]] PageLookup pageOf(const TerrainPoint& p, std::uint8_t level) const noexcept
    {
        if (level >= levels_)
        {
            return {TileKey{0, 0, level}, false};
        }
        const double side = pageSide(level).si();
        const double fi = std::floor(p.value().x.si() / side);
        const double fj = std::floor(p.value().y.si() / side);
        // Negated range tests so that NaN fails; the bounds are exact doubles.
        constexpr double Lowest = static_cast<double>(std::numeric_limits<std::int32_t>::min());
        constexpr double Highest = static_cast<double>(std::numeric_limits<std::int32_t>::max());
        if (!(fi >= Lowest && fi <= Highest && fj >= Lowest && fj <= Highest))
        {
            return {TileKey{0, 0, level}, false};
        }
        return {TileKey{static_cast<std::int32_t>(fi), static_cast<std::int32_t>(fj), level}, true};
    }

    /**
     * Position of a page's corner in the terrain frame (z = 0).
     *
     * @param key The page; its level below `levels()`.
     * @return Its minimum-x, minimum-y corner.
     */
    [[nodiscard]] constexpr sub0math::Point3<sub0math::f64::Length> pageCorner(const TileKey& key) const noexcept
    {
        const auto side = pageSide(key.level);
        return {side * static_cast<double>(key.i), side * static_cast<double>(key.j), sub0math::f64::Length{}};
    }

private:
    constexpr GridSpec(sub0math::f64::Length finestSpacing, std::uint8_t levels) noexcept
        : finestSpacing_(finestSpacing), levels_(levels)
    {
    }

    sub0math::f64::Length finestSpacing_; ///< Node spacing at level 0 [m].
    std::uint8_t levels_;                 ///< Number of levels.
};

/**
 * Places a page's local frame in the terrain frame: same axes, translated to the page's corner and to
 * the page's height origin. Local coordinates therefore stay small, and float keeps them exact enough
 * at any distance from the terrain origin.
 */
struct PageFrame
{
    TileKey key;                                      ///< The page.
    sub0math::Point3<sub0math::f64::Length> origin;   ///< Local (0, 0, 0) in the terrain frame.

    friend constexpr bool operator==(const PageFrame&, const PageFrame&) = default;
};

/**
 * Converts a terrain-frame point to page-local coordinates. The subtraction is done in double, then
 * narrowed, so no large magnitude is ever rounded to float.
 *
 * @param frame The page frame.
 * @param p The terrain point.
 * @return Its local coordinates.
 */
[[nodiscard]] constexpr LocalPoint toLocal(const PageFrame& frame, const TerrainPoint& p) noexcept
{
    const auto d = p.value() - frame.origin;
    return {static_cast<sub0math::Length>(d.x), static_cast<sub0math::Length>(d.y), static_cast<sub0math::Length>(d.z)};
}

/**
 * Converts page-local coordinates to the terrain frame.
 *
 * @param frame The page frame.
 * @param p The local point.
 * @return The terrain point.
 */
[[nodiscard]] constexpr TerrainPoint toTerrain(const PageFrame& frame, const LocalPoint& p) noexcept
{
    using L = sub0math::f64::Length;
    return TerrainPoint{frame.origin + sub0math::Vec3<L>{L{p.x}, L{p.y}, L{p.z}}};
}

} // namespace sub0physics::terrain
