#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <sub0math/geometry/point.hpp>
#include <sub0math/geometry/ray.hpp>
#include <sub0math/geometry/regular_grid.hpp>
#include <sub0math/geometry/triangle.hpp>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/core/column_buffer.hpp"
#include "sub0physics/geometry/witness.hpp"
#include "sub0physics/terrain/corridor_builder.hpp"
#include "sub0physics/terrain/grid_spec.hpp"
#include "sub0physics/terrain/height_page.hpp"

namespace sub0physics::terrain
{

struct HeightExtractWriter;
class HeightExtract;
struct HeightRequest;
struct PageWindow;
struct HeightResult;

/** Capacity of a `HeightExtract`. */
struct HeightExtractLimits
{
    std::size_t nodes = 20000;      ///< Most node samples.
    std::size_t triangles = 32768;  ///< Most triangles (two per delivered cell); also bounds the cells.
};

/**
 * The pages a request can read, as a rectangle of page pointers in row-major order: the caller pins the pages
 * it needs in a `SourceTileCache` and lists them here, so extraction itself needs no cache and no lookup.
 */
struct PageWindow
{
    sub0math::GridRange<2> pages;              ///< The page indices the window spans (`GridSpec::pagesTouching`).
    std::uint8_t level = 0;                    ///< Level of the pages.
    std::span<const HeightPage* const> data;   ///< One entry per page, row by row in j then i; null for a page that is absent.

    /**
     * Finds a page.
     *
     * @param i Page index along x.
     * @param j Page index along y.
     * @return The page, or null outside the window or where the entry is null.
     */
    [[nodiscard]] const HeightPage* at(std::int32_t i, std::int32_t j) const noexcept
    {
        if (i < pages.lo[0] || i > pages.hi[0] || j < pages.lo[1] || j > pages.hi[1])
        {
            return nullptr;
        }
        const auto width = static_cast<std::size_t>(pages.hi[0] - pages.lo[0] + 1);
        const auto n = static_cast<std::size_t>(j - pages.lo[1]) * width + static_cast<std::size_t>(i - pages.lo[0]);
        return n < data.size() ? data[n] : nullptr;
    }
};

/** What to extract. */
struct HeightRequest
{
    std::uint8_t level = 0;                       ///< Level of the pages.
    sub0math::GridRange<2> cells;                 ///< The block, in cell indices over the whole level (`GridSpec::cellsTouching`).
    std::uint32_t stride = 1;                     ///< Fine cells per extract cell: a power of two, 1 to 256. The resolution lever.
    const CellMask* mask = nullptr;               ///< Optional: only the cells it sets are requested (over the same cell indices).
};

/** Whether an extract cell was asked for and delivered. */
enum class CellState : std::uint8_t
{
    NotRequested = 0, ///< Outside the mask: nothing was asked, nothing delivered.
    Delivered = 1,    ///< Its two triangles are in the extract.
    Unknown = 2       ///< Asked for, but a node is unauthored, the page is absent or the cell is off the grid.
};

/** How `extractHeights` ended. */
enum class HeightStatus : std::uint8_t
{
    Complete,       ///< Every requested cell was delivered.
    Partial,        ///< Some requested cells are `Unknown`; the rest is usable, and nothing was flattened or filled in.
    Budget,         ///< The block needs more nodes or triangles than the extract holds; nothing was written.
    InvalidRequest  ///< Bad level, empty block, bad stride, or a mask that does not describe the block's cells.
};

/** The result of `extractHeights`. */
struct HeightResult
{
    HeightStatus status;              ///< How it ended.
    std::size_t deliveredCells = 0;   ///< Cells with triangles.
    std::size_t unknownCells = 0;     ///< Requested cells that could not be delivered.
};

/** Where a ray meets an extract. */
struct RayResult
{
    bool hit = false;                     ///< A triangle was hit.
    sub0math::Length distance{};          ///< Distance along the ray.
    Witness3 witness{};                   ///< Face witness: patch-local point, the unit normal facing the ray's origin, and the triangle index.
};

/**
 * Terrain heights as an owned structure-of-arrays patch (the grid layout of `docs/TERRAIN.md`): node samples in
 * rows, and two triangles per delivered cell on the pages' fixed triangulation, ready for ray and
 * closest-point queries.
 *
 * **Frame.** Patch-local: node (a, b) is at x = a * spacing(), y = b * spacing() from `origin()`, and heights are
 * relative to `origin().z`. Everything is `float` metres, so a patch far from the terrain origin keeps its
 * precision; `origin()` (double) is the only way back to the terrain frame.
 *
 * **Resolution.** `spacing()` is the spacing actually achieved (finest spacing times the stride), never the
 * requested one relabelled. A cell is delivered only if every finer node under it is authored, so a coarse
 * cell never interpolates across missing data.
 *
 * Allocated once by the constructor; refilled without allocating.
 */
class HeightExtract
{
public:
    using Length = sub0math::Length;
    using TriangleF = sub0math::Triangle<Length, 3>;

    static constexpr std::uint32_t NoTriangle = 0xFFFFFFFFU; ///< `cellTriangle` entry of a cell without triangles.

    /**
     * Reserves storage for the largest extract.
     *
     * @param limits Capacities.
     */
    explicit HeightExtract(const HeightExtractLimits& limits = {});

    /** @return The capacities. */
    [[nodiscard]] HeightExtractLimits limits() const noexcept
    {
        return {nodeZ_.capacity(), triangles_.capacity()};
    }

    /** @return The terrain-frame position of node (0, 0) and the height origin (z). */
    [[nodiscard]] const sub0math::Point3<sub0math::f64::Length>& origin() const noexcept
    {
        return origin_;
    }

    /** @return The spacing achieved between neighbouring nodes. */
    [[nodiscard]] Length spacing() const noexcept
    {
        return spacing_;
    }

    /** @return Fine cells per extract cell. */
    [[nodiscard]] std::uint32_t stride() const noexcept
    {
        return stride_;
    }

    /** @return Cells along x and y. */
    [[nodiscard]] std::array<std::int32_t, 2> cells() const noexcept
    {
        return {cellsX_, cellsY_};
    }

    /** @return The level cell index of the extract's first fine cell (the block's low corner, aligned to the stride). */
    [[nodiscard]] std::array<std::int32_t, 2> firstCell() const noexcept
    {
        return {firstX_, firstY_};
    }

    /** @return Node heights relative to `origin().z`, `(cells x + 1) * (cells y + 1)` of them, row by row in y. */
    [[nodiscard]] std::span<const Length> z() const noexcept
    {
        return nodeZ_.view();
    }

    /** @return 1 where the node is authored, else 0; one per node. */
    [[nodiscard]] std::span<const std::uint8_t> nodeKnown() const noexcept
    {
        return nodeKnown_.view();
    }

    /** @return The state of each cell, row by row in y. */
    [[nodiscard]] std::span<const CellState> cellState() const noexcept
    {
        return cellState_.view();
    }

    /** @return The triangles of the delivered cells, in cell order: the lower (a, b, d) then the upper (a, d, c). */
    [[nodiscard]] std::span<const TriangleF> triangles() const noexcept
    {
        return triangles_.view();
    }

    /** @return The unit normal of each triangle, out of the solid: derived from the same triangles. */
    [[nodiscard]] std::span<const sub0math::Vec3f> triangleNormal() const noexcept
    {
        return triangleNormal_.view();
    }

    /** @return For each cell, the index of its first triangle (its second follows), or `NoTriangle` for a cell that is not delivered. */
    [[nodiscard]] std::span<const std::uint32_t> cellTriangle() const noexcept
    {
        return cellTriangle_.view();
    }

    /** @return The cell index (row-major over `cells()`) of each triangle. */
    [[nodiscard]] std::span<const std::uint32_t> triangleCell() const noexcept
    {
        return triangleCell_.view();
    }

    /** @return Words for a corridor mask, for `extractCorridor`. */
    [[nodiscard]] std::span<std::uint64_t> maskWords() noexcept
    {
        return maskWords_.mutableView();
    }

private:
    friend struct HeightExtractWriter;
    friend HeightResult extractHeights(const GridSpec&, const PageWindow&, const HeightRequest&, HeightExtract&) noexcept;

    sub0math::Point3<sub0math::f64::Length> origin_{};
    Length spacing_{};
    std::uint32_t stride_ = 1;
    std::int32_t cellsX_ = 0;
    std::int32_t cellsY_ = 0;
    std::int32_t firstX_ = 0;
    std::int32_t firstY_ = 0;
    ColumnBuffer<Length> nodeZ_;
    ColumnBuffer<std::uint8_t> nodeKnown_;
    ColumnBuffer<CellState> cellState_;
    ColumnBuffer<std::uint32_t> cellTriangle_;
    ColumnBuffer<TriangleF> triangles_;
    ColumnBuffer<sub0math::Vec3f> triangleNormal_;
    ColumnBuffer<std::uint32_t> triangleCell_;
    ColumnBuffer<std::uint64_t> maskWords_;
};

/**
 * Extracts a block of pages into node samples and triangles.
 *
 * The block is widened outward to whole extract cells (multiples of the stride on the level's cell grid), then
 * read page by page: heights as they are stored, triangles on the page's fixed diagonal, in `float`, relative to
 * the patch origin. A cell with an unauthored node is `Unknown`, never flat. Nothing is truncated: a block over
 * the capacity is refused.
 *
 * @param grid The grid the pages belong to (256 cells, node-registered).
 * @param pages The pages the request can read.
 * @param request The block, stride and optional mask.
 * @param out The result; cleared first.
 * @return The status and counts.
 */
[[nodiscard]] HeightResult extractHeights(const GridSpec& grid, const PageWindow& pages, const HeightRequest& request,
                                          HeightExtract& out) noexcept;

/**
 * Extracts the terrain under a corridor: rasterises the corridor on the level's cell grid, then extracts the
 * cells it touches (the mask's bounding block, with only the touched cells requested).
 *
 * @param grid The grid.
 * @param corridor The swept corridor.
 * @param pages The pages it can read; build the window from `GridSpec::pagesTouching(corridor.bounds(), level)`.
 * @param level Level of the pages.
 * @param stride Fine cells per extract cell, as in `HeightRequest`.
 * @param out The result; its `maskWords()` hold the corridor's mask while the call runs.
 * @return The status and counts; InvalidRequest when the corridor misses the grid or the mask words are too few.
 */
[[nodiscard]] HeightResult extractCorridor(const GridSpec& grid, const CorridorBuilder& corridor, const PageWindow& pages,
                                           std::uint8_t level, std::uint32_t stride, HeightExtract& out) noexcept;

/**
 * Casts a ray against an extract's delivered triangles, walking the cells the ray crosses.
 *
 * @param extract The patch.
 * @param ray The ray in the patch-local frame; the direction is a unit vector.
 * @param maxDistance Ignore hits farther than this.
 * @return The nearest hit with a face witness, or `hit == false`.
 */
[[nodiscard]] RayResult castRay(const HeightExtract& extract, const sub0math::Ray<sub0math::Length, 3>& ray,
                                sub0math::Length maxDistance) noexcept;

/**
 * Casts a ray against every delivered triangle, with no acceleration. The reference that `castRay` is checked
 * against, and what a tiny extract should use.
 *
 * @param extract The patch.
 * @param ray The ray in the patch-local frame.
 * @param maxDistance Ignore hits farther than this.
 * @return The nearest hit, or `hit == false`.
 */
[[nodiscard]] RayResult castRayBruteForce(const HeightExtract& extract, const sub0math::Ray<sub0math::Length, 3>& ray,
                                          sub0math::Length maxDistance) noexcept;

} // namespace sub0physics::terrain
