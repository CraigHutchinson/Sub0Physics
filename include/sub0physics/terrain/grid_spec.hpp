#pragma once

#include <cmath>
#include <compare>
#include <cstdint>
#include <sub0math/geometry/point.hpp>
#include <sub0math/role/role.hpp>
#include <sub0math/units/quantities.hpp>

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
 * How terrain pages tile the plane. Every page has the same number of cells; a page at level k has
 * node spacing `finestSpacing * 2^k`, so its world size doubles per level.
 */
struct GridSpec
{
    static constexpr std::int32_t CellsPerSide = 256; ///< Cells along each side of every page.

    sub0math::f64::Length finestSpacing; ///< Node spacing at level 0; positive.
    std::uint8_t levels;                 ///< Number of levels; at least 1, at most 16.

    /**
     * Node spacing at a level.
     *
     * @param level Below `levels`.
     * @return finestSpacing * 2^level.
     */
    [[nodiscard]] constexpr sub0math::f64::Length spacing(std::uint8_t level) const noexcept
    {
        return finestSpacing * static_cast<double>(1U << level);
    }

    /**
     * Side length of a page at a level.
     *
     * @param level Below `levels`.
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
     * @param level Below `levels`.
     * @return The page's key.
     */
    [[nodiscard]] TileKey pageOf(const TerrainPoint& p, std::uint8_t level) const noexcept
    {
        const double side = pageSide(level).si();
        return {static_cast<std::int32_t>(std::floor(p.value().x.si() / side)),
                static_cast<std::int32_t>(std::floor(p.value().y.si() / side)), level};
    }

    /**
     * Position of a page's corner in the terrain frame (z = 0).
     *
     * @param key The page.
     * @return Its minimum-x, minimum-y corner.
     */
    [[nodiscard]] constexpr sub0math::Point3<sub0math::f64::Length> pageCorner(const TileKey& key) const noexcept
    {
        const auto side = pageSide(key.level);
        return {side * static_cast<double>(key.i), side * static_cast<double>(key.j), sub0math::f64::Length{}};
    }
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
