#pragma once

#include <algorithm>
#include <array>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <sub0math/geometry/aabb.hpp>
#include <sub0math/geometry/point.hpp>
#include <sub0math/geometry/regular_grid.hpp>
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
 * non-finite point, a level outside the grid, or a point outside the grid's rectangle; `key` is
 * meaningless then and must not be read, as with dereferencing an empty optional.
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
 * The outcome of locating a cell: its page, the cell's index within the page and whether there is one.
 * `known` is false as for `PageLookup`; the other fields are meaningless then.
 */
struct CellLookup
{
    TileKey page;   ///< The page that owns the cell, when known.
    std::int32_t i; ///< Cell index along x within the page, 0 to cellsPerSide - 1.
    std::int32_t j; ///< Cell index along y within the page.
    bool known;     ///< The point lies on a cell of the grid.

    /** @return `known`, so a lookup can be tested like an optional. */
    [[nodiscard]] constexpr explicit operator bool() const noexcept
    {
        return known;
    }
};

/** How samples sit on the grid; recorded when the data is imported, because it is a property of the source. */
enum class Registration : std::uint8_t
{
    NodeIsPoint, ///< A sample is the value at a grid node (pixel-is-point): cells + 1 samples per page side.
    CellCentre   ///< A sample is the value over a cell, placed at its centre (pixel-is-area): cells samples per side.
};

/**
 * The shape of a grid, for `GridSpec::create`. Designated initialisers name the fields; only the spacing
 * is required.
 */
struct GridLayout
{
    sub0math::f64::Length finestSpacing;                    ///< Cell side at level 0 [m]; within [MinSpacing, MaxSpacing].
    std::uint8_t levels = 1;                                ///< Number of levels, 1 to MaxLevels; each doubles the spacing.
    std::int32_t cellsPerSide = 256;                        ///< Cells along each side of every page, 1 to MaxCellsPerSide.
    Registration registration = Registration::NodeIsPoint;  ///< Where a sample sits within its cell.
    sub0math::Point2<sub0math::f64::Length> anchor{};       ///< Terrain-frame position of page (0, 0)'s minimum corner [m].
    std::array<std::int32_t, 2> pages{1 << 16, 1 << 16};    ///< Pages along x and y at level 0; at least 1 each.
};

/**
 * How terrain pages tile the plane. Every page has the same number of cells; a page at level k has cell
 * spacing `finestSpacing * 2^k`, so its world size doubles per level. The grid covers a finite rectangle
 * of pages that starts at the anchor and extends towards +x and +y; a point outside it has no page.
 *
 * Built on `sub0math::RegularGrid`: at each level the cells of the whole rectangle form one regular grid,
 * which gives two operations their rules.
 * - **Ownership** (`pageOf`, `cellOf`): each point belongs to exactly one cell, and so one page. Cells are
 *   half-open: a point on a shared edge belongs to the cell on the positive side, and the far edge of the
 *   grid belongs to the last cell.
 * - **Coverage** (`cellsTouching`, `pagesTouching`): the cells (pages) a closed box touches. A box that
 *   reaches an edge exactly includes the cells on both sides.
 *
 * Built only by `create`, which validates the layout, so every grid in existence is valid, the
 * `1 << level` arithmetic cannot overflow and every cell index fits 32 bits.
 */
class GridSpec
{
public:
    static constexpr std::int32_t DefaultCellsPerSide = 256;  ///< Cells per page side unless the layout says otherwise.
    static constexpr std::int32_t MaxCellsPerSide = 1 << 16;  ///< Most cells along a page side.
    static constexpr std::int32_t MaxCellsPerAxis = 1 << 30;  ///< Most cells along the whole grid at level 0.
    static constexpr std::uint8_t MaxLevels = 16;             ///< Most levels a grid may have.
    static constexpr double MinSpacing = 1e-6;                ///< Smallest finest spacing [m]; keeps the float reciprocal finite.
    static constexpr double MaxSpacing = 1e6;                 ///< Largest finest spacing [m]; keeps page sides far from overflow.
    static constexpr double MaxAnchor = 1e9;                  ///< Largest anchor coordinate magnitude [m].

    GridSpec() = delete;

    /**
     * Validates and builds a grid.
     *
     * @param layout The grid's shape: spacing, levels, cells per page, registration, anchor and extent.
     * @return The grid, or an `InvalidInput` fault saying what is out of range.
     */
    [[nodiscard]] static constexpr std::expected<GridSpec, Fault> create(const GridLayout& layout) noexcept
    {
        if (layout.levels < 1U || layout.levels > MaxLevels)
        {
            return invalidInput("grid needs 1 to 16 levels");
        }
        const double s = layout.finestSpacing.si();
        if (!(s >= MinSpacing && s <= MaxSpacing)) // NaN fails both
        {
            return invalidInput("grid spacing must be finite, from 1 micrometre to 1000 km");
        }
        if (layout.cellsPerSide < 1 || layout.cellsPerSide > MaxCellsPerSide)
        {
            return invalidInput("grid needs 1 to 65536 cells per page side");
        }
        if (layout.registration != Registration::NodeIsPoint && layout.registration != Registration::CellCentre)
        {
            return invalidInput("grid registration is not a known value");
        }
        const double ax = layout.anchor.x.si();
        const double ay = layout.anchor.y.si();
        if (!(ax >= -MaxAnchor && ax <= MaxAnchor && ay >= -MaxAnchor && ay <= MaxAnchor)) // NaN fails
        {
            return invalidInput("grid anchor must be finite and within 1e9 m");
        }
        for (const std::int32_t pages : layout.pages)
        {
            if (pages < 1 || std::int64_t{pages} * layout.cellsPerSide > MaxCellsPerAxis)
            {
                return invalidInput("grid needs at least 1 page per axis and at most 2^30 cells per axis");
            }
        }
        return GridSpec{layout};
    }

    /**
     * Builds the default-extent grid centred on the terrain origin: node-registered, 256 cells per page
     * and 65536 x 65536 pages at level 0. (At a spacing so large that the corner would pass `MaxAnchor`,
     * the anchor stops there and the grid is no longer centred.)
     *
     * @param finestSpacing Cell spacing at level 0 [m]; finite and within [MinSpacing, MaxSpacing].
     * @param levels Number of levels, 1 to MaxLevels.
     * @return The grid, or an `InvalidInput` fault.
     */
    [[nodiscard]] static constexpr std::expected<GridSpec, Fault> create(sub0math::f64::Length finestSpacing,
                                                                         std::uint8_t levels) noexcept
    {
        GridLayout layout{.finestSpacing = finestSpacing, .levels = levels};
        const double half = std::max(-MaxAnchor, -0.5 * static_cast<double>(layout.pages[0]) * static_cast<double>(layout.cellsPerSide) *
                                                      finestSpacing.si());
        layout.anchor = {sub0math::f64::Length::fromSi(half), sub0math::f64::Length::fromSi(half)};
        return create(layout);
    }

    /** @return Cell spacing at level 0. */
    [[nodiscard]] constexpr sub0math::f64::Length finestSpacing() const noexcept
    {
        return layout_.finestSpacing;
    }

    /** @return The number of levels, 1 to MaxLevels. */
    [[nodiscard]] constexpr std::uint8_t levels() const noexcept
    {
        return layout_.levels;
    }

    /** @return Cells along each side of every page. */
    [[nodiscard]] constexpr std::int32_t cellsPerSide() const noexcept
    {
        return layout_.cellsPerSide;
    }

    /** @return Where a sample sits within its cell. */
    [[nodiscard]] constexpr Registration registration() const noexcept
    {
        return layout_.registration;
    }

    /** @return Terrain-frame position of page (0, 0)'s minimum corner. */
    [[nodiscard]] constexpr sub0math::Point2<sub0math::f64::Length> anchor() const noexcept
    {
        return layout_.anchor;
    }

    /** @return Samples along each page side: cells + 1 for node registration, cells for cell-centre. */
    [[nodiscard]] constexpr std::int32_t samplesPerSide() const noexcept
    {
        return layout_.registration == Registration::NodeIsPoint ? layout_.cellsPerSide + 1 : layout_.cellsPerSide;
    }

    /**
     * Cell spacing at a level.
     *
     * @param level Below `levels()`; larger values are a contract violation.
     * @return finestSpacing * 2^level.
     */
    [[nodiscard]] constexpr sub0math::f64::Length spacing(std::uint8_t level) const noexcept
    {
        return layout_.finestSpacing * static_cast<double>(1U << level);
    }

    /**
     * Side length of a page at a level.
     *
     * @param level Below `levels()`.
     * @return spacing(level) * cellsPerSide().
     */
    [[nodiscard]] constexpr sub0math::f64::Length pageSide(std::uint8_t level) const noexcept
    {
        return spacing(level) * static_cast<double>(layout_.cellsPerSide);
    }

    /**
     * The authored rectangle in the terrain frame: from the anchor, `pages * cellsPerSide` finest cells along
     * each axis. (Coarse pages may extend past it; points outside it have no page.)
     *
     * @return The rectangle, lower corner first.
     */
    [[nodiscard]] constexpr sub0math::Aabb<sub0math::f64::Length, 2> domain() const noexcept
    {
        const auto s = layout_.finestSpacing;
        return {layout_.anchor,
                {layout_.anchor.x + s * static_cast<double>(std::int64_t{layout_.pages[0]} * layout_.cellsPerSide),
                 layout_.anchor.y + s * static_cast<double>(std::int64_t{layout_.pages[1]} * layout_.cellsPerSide)}};
    }

    /**
     * Reports whether a closed box lies inside the authored rectangle (edges included).
     *
     * @param box The box in the terrain frame.
     * @return True when it does; false for any box that leaves the rectangle or has a non-finite corner.
     */
    [[nodiscard]] constexpr bool covers(const sub0math::Aabb<sub0math::f64::Length, 2>& box) const noexcept
    {
        const auto d = domain();
        return box.lo.x >= d.lo.x && box.hi.x <= d.hi.x && box.lo.y >= d.lo.y && box.hi.y <= d.hi.y;
    }

    /**
     * Number of pages along one axis at a level: the level-0 count halved (rounding up) per level, so a
     * coarse page may extend past the authored rectangle.
     *
     * @param level Below `levels()`.
     * @param axis 0 for x, 1 for y.
     * @return The page count, at least 1.
     */
    [[nodiscard]] constexpr std::int32_t pageCount(std::uint8_t level, std::size_t axis) const noexcept
    {
        const std::int64_t pages = layout_.pages[axis];
        return static_cast<std::int32_t>((pages + (std::int64_t{1} << level) - 1) >> level);
    }

    /**
     * The regular grid of all cells at a level: the one object that ownership and coverage are defined on.
     *
     * @param level Below `levels()`.
     * @return Origin at the anchor, spacing(level), and one node more than the cells on each axis.
     */
    [[nodiscard]] constexpr sub0math::RegularGrid<sub0math::f64::Length, 2> cellGrid(std::uint8_t level) const noexcept
    {
        return {layout_.anchor, spacing(level),
                {pageCount(level, 0) * layout_.cellsPerSide + 1, pageCount(level, 1) * layout_.cellsPerSide + 1}};
    }

    /**
     * Finds the page that owns a point at a level. Points on a shared edge belong to the page on the
     * positive side.
     *
     * @param p The point; only x and y are used.
     * @param level The level to look at.
     * @return The page's key, or not known for a level outside the grid, a non-finite coordinate or a
     *         point outside the grid's rectangle.
     */
    [[nodiscard]] PageLookup pageOf(const TerrainPoint& p, std::uint8_t level) const noexcept
    {
        const CellLookup c = cellOf(p, level);
        return {c.page, c.known};
    }

    /**
     * Finds the cell that owns a point at a level, as a page and the cell's index within it.
     *
     * @param p The point; only x and y are used.
     * @param level The level to look at.
     * @return The cell, or not known as for `pageOf`.
     */
    [[nodiscard]] CellLookup cellOf(const TerrainPoint& p, std::uint8_t level) const noexcept
    {
        if (level >= layout_.levels)
        {
            return {TileKey{0, 0, level}, 0, 0, false};
        }
        const auto cell = cellGrid(level).cellOf({p.value().x, p.value().y});
        if (!cell)
        {
            return {TileKey{0, 0, level}, 0, 0, false};
        }
        const std::int32_t n = layout_.cellsPerSide;
        return {TileKey{cell.index[0] / n, cell.index[1] / n, level}, cell.index[0] % n, cell.index[1] % n, true};
    }

    /**
     * Finds the cells that a closed box touches at a level, clipped to the grid.
     *
     * @param box The box in the terrain frame; both corners finite.
     * @param level The level to look at.
     * @return Inclusive range of cell indices over the whole grid (page * cellsPerSide + cell); empty when
     *         the box misses the grid, a corner is not finite or the level is outside the grid.
     */
    [[nodiscard]] sub0math::GridRange<2> cellsTouching(const sub0math::Aabb<sub0math::f64::Length, 2>& box,
                                                       std::uint8_t level) const noexcept
    {
        if (level >= layout_.levels)
        {
            return {{0, 0}, {-1, -1}};
        }
        return cellGrid(level).cellsTouching(box);
    }

    /**
     * Finds the pages that a closed box touches at a level: the pages of the cells of `cellsTouching`.
     *
     * @param box The box in the terrain frame; both corners finite.
     * @param level The level to look at.
     * @return Inclusive range of page indices (the TileKey i and j); empty as for `cellsTouching`.
     */
    [[nodiscard]] sub0math::GridRange<2> pagesTouching(const sub0math::Aabb<sub0math::f64::Length, 2>& box,
                                                       std::uint8_t level) const noexcept
    {
        sub0math::GridRange<2> r = cellsTouching(box, level);
        if (r.empty())
        {
            return r;
        }
        for (std::size_t axis = 0; axis < 2; ++axis)
        {
            r.lo[axis] /= layout_.cellsPerSide;
            r.hi[axis] /= layout_.cellsPerSide;
        }
        return r;
    }

    /**
     * Position of a page's corner in the terrain frame (z = 0).
     *
     * @param key The page; its level below `levels()`.
     * @return Its minimum-x, minimum-y corner: the node at the page's first cell.
     */
    [[nodiscard]] constexpr sub0math::Point3<sub0math::f64::Length> pageCorner(const TileKey& key) const noexcept
    {
        const auto s = spacing(key.level);
        return {layout_.anchor.x + s * firstCell(key.i), layout_.anchor.y + s * firstCell(key.j), sub0math::f64::Length{}};
    }

    /**
     * Position of a sample in the terrain frame (z = 0): the node itself for node registration, the
     * centre of the cell for cell-centre registration.
     *
     * @param key The page.
     * @param i Sample index along x, below `samplesPerSide()`.
     * @param j Sample index along y, below `samplesPerSide()`.
     * @return Where that sample was taken.
     */
    [[nodiscard]] constexpr sub0math::Point3<sub0math::f64::Length> samplePosition(const TileKey& key, std::int32_t i,
                                                                                  std::int32_t j) const noexcept
    {
        const auto s = spacing(key.level);
        const double offset = layout_.registration == Registration::CellCentre ? 0.5 : 0.0;
        return {layout_.anchor.x + s * (firstCell(key.i) + static_cast<double>(i) + offset),
                layout_.anchor.y + s * (firstCell(key.j) + static_cast<double>(j) + offset), sub0math::f64::Length{}};
    }

private:
    constexpr explicit GridSpec(const GridLayout& layout) noexcept : layout_(layout)
    {
    }

    // Index of a page's first cell, as an exact double (the product is formed in 64 bits).
    [[nodiscard]] constexpr double firstCell(std::int32_t page) const noexcept
    {
        return static_cast<double>(std::int64_t{page} * layout_.cellsPerSide);
    }

    GridLayout layout_; ///< The validated layout.
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
