#pragma once

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>
#include <sub0math/geometry/regular_grid.hpp>
#include <sub0math/linear/vec.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/core/column_buffer.hpp"
#include "sub0physics/terrain/grid_spec.hpp"

namespace sub0physics::terrain
{

/** One occupied vertical interval of a column: a layer of one material. */
struct ColumnSpan
{
    sub0math::f64::Length upper; ///< Top of the span (z up) [m].
    sub0math::f64::Length lower; ///< Bottom of the span; below `upper` [m].
    std::uint8_t material;       ///< The product's material code; not interpreted here.
};

/**
 * What `extractColumns` reads from a tile of cell-centred columns (one sample per cell). A tile is
 * addressed by the cell's index within it, `j * cellsPerSide + i`.
 */
template<typename T>
concept ColumnTile = requires(const T& tile, std::size_t cell, std::size_t k) {
    { tile.known(cell) } -> std::convertible_to<bool>;                  // false: unauthored, never treated as flat
    { tile.top(cell) } -> std::same_as<sub0math::f64::Length>;          // top surface of the column
    { tile.normal(cell) } -> std::same_as<sub0math::Vec3d>;             // unit normal of the top, out of the solid
    { tile.spanCount(cell) } -> std::convertible_to<std::size_t>;       // spans of the column
    { tile.span(cell, k) } -> std::same_as<ColumnSpan>;                 // span k, top first
};

/** Capacity of a `ColumnExtract`. */
struct ColumnLimits
{
    std::size_t columns = 4096; ///< Most columns, and most unknown cells reported.
    std::size_t spans = 12288;  ///< Most spans over all columns.
};

namespace detail
{
struct ColumnWriter;
} // namespace detail

/**
 * Owned structure-of-arrays result of `extractColumns`: for every cell of the requested block, row by row,
 * its identity, centre, top surface, normal and spans. Allocated once by the constructor, never aliases a
 * tile (so eviction cannot change it), and refilled without allocating.
 *
 * Column j's spans are `[offset()[j], offset()[j + 1])` of `upper()`, `lower()` and `material()`.
 */
class ColumnExtract
{
public:
    using Length = sub0math::f64::Length;

    /**
     * Reserves storage for the largest extract.
     *
     * @param limits Capacities.
     */
    explicit ColumnExtract(const ColumnLimits& limits = {})
    {
        ids_.reserve(limits.columns);
        x_.reserve(limits.columns);
        y_.reserve(limits.columns);
        top_.reserve(limits.columns);
        normal_.reserve(limits.columns);
        offset_.reserve(limits.columns + 1U);
        upper_.reserve(limits.spans);
        lower_.reserve(limits.spans);
        material_.reserve(limits.spans);
        missing_.reserve(limits.columns);
    }

    /** @return The capacities. */
    [[nodiscard]] ColumnLimits limits() const noexcept
    {
        return {ids_.capacity(), upper_.capacity()};
    }

    /** Empties every column; the capacity stays. */
    void clear() noexcept
    {
        ids_.clear();
        x_.clear();
        y_.clear();
        top_.clear();
        normal_.clear();
        offset_.clear();
        upper_.clear();
        lower_.clear();
        material_.clear();
        missing_.clear();
    }

    /** @return Cell identity `j * cellsAcross + i` over the whole level (cells along x times j, plus i), one per column. */
    [[nodiscard]] std::span<const std::uint64_t> ids() const noexcept
    {
        return ids_.view();
    }

    /** @return Cell-centre x of each column in the terrain frame. */
    [[nodiscard]] std::span<const Length> x() const noexcept
    {
        return x_.view();
    }

    /** @return Cell-centre y of each column in the terrain frame. */
    [[nodiscard]] std::span<const Length> y() const noexcept
    {
        return y_.view();
    }

    /** @return Top surface height of each column. */
    [[nodiscard]] std::span<const Length> top() const noexcept
    {
        return top_.view();
    }

    /** @return Outward unit normal of each column's top surface. */
    [[nodiscard]] std::span<const sub0math::Vec3d> normal() const noexcept
    {
        return normal_.view();
    }

    /** @return Span offsets: columns + 1 entries. */
    [[nodiscard]] std::span<const std::uint32_t> offset() const noexcept
    {
        return offset_.view();
    }

    /** @return Top of each span. */
    [[nodiscard]] std::span<const Length> upper() const noexcept
    {
        return upper_.view();
    }

    /** @return Bottom of each span. */
    [[nodiscard]] std::span<const Length> lower() const noexcept
    {
        return lower_.view();
    }

    /** @return Material code of each span. */
    [[nodiscard]] std::span<const std::uint8_t> material() const noexcept
    {
        return material_.view();
    }

    /** @return Identities of the unauthored cells, when the result is Partial. */
    [[nodiscard]] std::span<const std::uint64_t> missing() const noexcept
    {
        return missing_.view();
    }

private:
    friend struct detail::ColumnWriter;

    ColumnBuffer<std::uint64_t> ids_;
    ColumnBuffer<Length> x_;
    ColumnBuffer<Length> y_;
    ColumnBuffer<Length> top_;
    ColumnBuffer<sub0math::Vec3d> normal_;
    ColumnBuffer<std::uint32_t> offset_;
    ColumnBuffer<Length> upper_;
    ColumnBuffer<Length> lower_;
    ColumnBuffer<std::uint8_t> material_;
    ColumnBuffer<std::uint64_t> missing_;
};

/** How `extractColumns` ended. Only Complete leaves columns in the extract. */
enum class ColumnStatus : std::uint8_t
{
    Complete,       ///< Every requested cell was delivered.
    Partial,        ///< Some cells are unauthored: `missing()` lists them and no column is delivered.
    ColumnBudget,   ///< More cells than the extract holds.
    SpanBudget,     ///< More spans than the extract holds.
    MissingTile,    ///< The lookup had no tile for a page the block needs (a caller error: pin tiles first).
    UnsupportedGrid ///< The grid is not cell-centre registered, or the level is outside it.
};

/** The result of `extractColumns`. */
struct ColumnResult
{
    ColumnStatus status;       ///< How it ended.
    std::size_t columns = 0;   ///< Columns delivered (Complete) or requested.
    std::size_t spans = 0;     ///< Spans delivered (Complete) or counted.
};

namespace detail
{

struct ColumnWriter
{
    static void prepare(ColumnExtract& out, std::size_t columns, std::size_t spans) noexcept
    {
        // Capacities were checked by the caller.
        (void)out.ids_.resize(columns);
        (void)out.x_.resize(columns);
        (void)out.y_.resize(columns);
        (void)out.top_.resize(columns);
        (void)out.normal_.resize(columns);
        (void)out.offset_.resize(columns + 1U);
        (void)out.upper_.resize(spans);
        (void)out.lower_.resize(spans);
        (void)out.material_.resize(spans);
    }

    struct Pointers
    {
        std::uint64_t* ids;
        sub0math::f64::Length* x;
        sub0math::f64::Length* y;
        sub0math::f64::Length* top;
        sub0math::Vec3d* normal;
        std::uint32_t* offset;
        sub0math::f64::Length* upper;
        sub0math::f64::Length* lower;
        std::uint8_t* material;
    };

    static Pointers pointers(ColumnExtract& out) noexcept
    {
        return {out.ids_.data(),    out.x_.data(),     out.y_.data(),     out.top_.data(),   out.normal_.data(),
                out.offset_.data(), out.upper_.data(), out.lower_.data(), out.material_.data()};
    }

    static bool addMissing(ColumnExtract& out, std::uint64_t id) noexcept
    {
        return out.missing_.pushBack(id);
    }
};

} // namespace detail

/**
 * Extracts a block of cell-centred columns into an owned structure-of-arrays result (the grid layout with
 * spans): identity, centre, top, normal and spans per cell, row by row in `j` then `i`.
 *
 * Reading runs along a row that share a page, so the page lookup and index arithmetic are paid once per
 * run, not per cell. Unauthored cells decide the outcome first: any makes the result Partial with their
 * identities and no columns, because unknown ground is never delivered as flat. Nothing is truncated: a
 * block over the column or span capacity is refused.
 *
 * @tparam Tile A `ColumnTile`.
 * @tparam Lookup Callable `const Tile*(const TileKey&)`: the tile of a page, or null.
 * @param grid A cell-centre-registered grid.
 * @param level Level of the block.
 * @param cells The block, an inclusive range of cell indices over the whole level (`GridSpec::cellsTouching`).
 * @param tileOf Finds a page's tile; must be cheap. Called at least once per run of a row.
 * @param out The result; cleared first.
 * @return The status and counts.
 */
template<ColumnTile Tile, typename Lookup>
[[nodiscard]] ColumnResult extractColumns(const GridSpec& grid, std::uint8_t level, const sub0math::GridRange<2>& cells,
                                          Lookup&& tileOf, ColumnExtract& out) noexcept
{
    out.clear();
    if (grid.registration() != Registration::CellCentre || level >= grid.levels())
    {
        return {ColumnStatus::UnsupportedGrid};
    }
    if (cells.empty())
    {
        return {ColumnStatus::Complete};
    }
    const auto n = static_cast<std::size_t>(cells.count());
    if (n > out.limits().columns)
    {
        return {ColumnStatus::ColumnBudget, n};
    }
    const std::int32_t side = grid.cellsPerSide();
    const auto across = static_cast<std::uint64_t>(grid.pageCount(level, 0)) * static_cast<std::uint64_t>(side);

    // Visits the block row by row in runs of cells that share a page: f(tile, first cell in the page, i0, i1, j).
    // PERF: page runs and indexed writes 1.6x faster than per-cell.
    bool missingTile = false;
    const auto forEachRun = [&](auto&& f) {
        for (std::int32_t j = cells.lo[1]; j <= cells.hi[1]; ++j)
        {
            for (std::int32_t i = cells.lo[0]; i <= cells.hi[0];)
            {
                const std::int32_t end = std::min(cells.hi[0], (i / side + 1) * side - 1);
                const Tile* tile = tileOf(TileKey{i / side, j / side, level});
                if (tile == nullptr)
                {
                    missingTile = true;
                    return;
                }
                f(*tile, static_cast<std::size_t>((j % side) * side + i % side), i, end, j);
                i = end + 1;
            }
        }
    };

    // One pass finds the unauthored cells and counts the spans.
    std::size_t spans = 0;
    forEachRun([&](const Tile& tile, std::size_t index, std::int32_t i0, std::int32_t i1, std::int32_t j) {
        for (std::int32_t i = i0; i <= i1; ++i, ++index)
        {
            if (!tile.known(index))
            {
                // The missing list holds as many as there are columns, and that was checked above.
                (void)detail::ColumnWriter::addMissing(out, static_cast<std::uint64_t>(j) * across + static_cast<std::uint64_t>(i));
            }
            spans += tile.spanCount(index);
        }
    });
    if (missingTile)
    {
        out.clear();
        return {ColumnStatus::MissingTile, n};
    }
    if (!out.missing().empty())
    {
        return {ColumnStatus::Partial, n, spans};
    }
    if (spans > out.limits().spans)
    {
        return {ColumnStatus::SpanBudget, n, spans};
    }
    detail::ColumnWriter::prepare(out, n, spans);

    const double h = grid.spacing(level).si();
    const double x0 = grid.anchor().x.si();
    const double y0 = grid.anchor().y.si();
    const auto p = detail::ColumnWriter::pointers(out);
    std::uint64_t* ids = p.ids;
    sub0math::f64::Length* xs = p.x;
    sub0math::f64::Length* ys = p.y;
    sub0math::f64::Length* tops = p.top;
    sub0math::Vec3d* normals = p.normal;
    std::uint32_t* offsets = p.offset;
    sub0math::f64::Length* uppers = p.upper;
    sub0math::f64::Length* lowers = p.lower;
    std::uint8_t* materials = p.material;
    std::size_t column = 0;
    std::uint32_t span = 0;
    forEachRun([&](const Tile& tile, std::size_t index, std::int32_t i0, std::int32_t i1, std::int32_t j) {
        const auto cy = sub0math::f64::Length::fromSi(y0 + h * (static_cast<double>(j) + 0.5));
        const std::uint64_t rowBase = static_cast<std::uint64_t>(j) * across;
        for (std::int32_t i = i0; i <= i1; ++i, ++index, ++column)
        {
            ids[column] = rowBase + static_cast<std::uint64_t>(i);
            xs[column] = sub0math::f64::Length::fromSi(x0 + h * (static_cast<double>(i) + 0.5));
            ys[column] = cy;
            tops[column] = tile.top(index);
            normals[column] = tile.normal(index);
            offsets[column] = span;
            const std::size_t count = tile.spanCount(index);
            for (std::size_t k = 0; k < count; ++k, ++span)
            {
                const ColumnSpan s = tile.span(index, k);
                uppers[span] = s.upper;
                lowers[span] = s.lower;
                materials[span] = s.material;
            }
        }
    });
    offsets[column] = span;
    return {ColumnStatus::Complete, n, spans};
}

} // namespace sub0physics::terrain
