#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <sub0math/config.hpp>
#include <sub0math/geometry/point.hpp>
#include <sub0math/role/role.hpp>
#include <sub0math/units/quantities.hpp>
#include <vector>

#include "sub0physics/core/id.hpp"
#include "sub0physics/fault.hpp"
#include "sub0physics/terrain/grid_spec.hpp"

namespace sub0physics::terrain
{

struct SourceRevisionTag;

/** The immutable content generation of a page; changes whenever its authored data changes. */
using SourceRevision = Id<SourceRevisionTag>;

/**
 * Conservative height bounds of a block of cells, in the page's 0.1 mm steps.
 *
 * Covers every authored node of the block. `anyUnknown` means some node has no data, so the bounds
 * cannot prove anything about that part of the block.
 */
struct HeightRange
{
    std::uint16_t lo;   ///< Lowest authored height, in steps.
    std::uint16_t hi;   ///< Highest authored height, in steps.
    bool anyKnown;      ///< At least one node is authored (otherwise lo and hi are meaningless).
    bool anyUnknown;    ///< At least one node is not authored.

    friend constexpr bool operator==(const HeightRange&, const HeightRange&) = default;
};

/**
 * A surface height, or the fact that there is none: `known` is false outside the page or where the
 * triangle under the point has an unauthored node. `height` is meaningless then and must not be read,
 * as with dereferencing an empty optional.
 */
struct HeightSample
{
    sub0math::Length height; ///< Height above the page's height origin, when known.
    bool known;              ///< The surface is defined here.

    /** @return `known`, so a sample can be tested like an optional. */
    [[nodiscard]] constexpr explicit operator bool() const noexcept
    {
        return known;
    }
};

/**
 * One page of the terrain's base surface (representation B): heights on a 257 x 257 node grid
 * (256 x 256 cells), quantised to 0.1 mm above the page's height origin, with a mask of authored nodes
 * and a min/max pyramid for clearance proofs and query pruning.
 *
 * The page size is fixed at compile time (256 cells), so the lookup has constant strides and bounds and no
 * runtime page size to load. A grid of another cell count or of cell-centre registration describes a
 * different kind of source and is refused by `buildHeightPage`; see docs/TERRAIN.md for the measurement
 * behind that choice.
 *
 * Each page holds its own edge nodes; a neighbour's matching nodes are equal by construction (stitched
 * at import), while each cell, and so its state, belongs to exactly one page.
 *
 * The surface between nodes is fixed: cell (i, j) is split along its (i, j)-(i+1, j+1) diagonal. Heights,
 * normals and contact all use that same triangulation. Unknown nodes are never interpolated across.
 */
class HeightPage
{
public:
    static constexpr std::int32_t Cells = GridSpec::DefaultCellsPerSide; ///< Cells per side: the 256 fast path.
    static constexpr std::int32_t Nodes = Cells + 1;                    ///< Nodes per side.
    static constexpr std::int32_t BlockCells = 8;                       ///< Cells per side of a pyramid level-0 block.
    static constexpr std::int32_t PyramidLevels = 6;                    ///< 32, 16, 8, 4, 2 and 1 blocks per side.
    static constexpr float StepMetres = 1e-4F;                          ///< Height quantum.

    /** @return Where the page sits in the terrain frame. */
    [[nodiscard]] const PageFrame& frame() const noexcept
    {
        return frame_;
    }

    /** @return The content generation. */
    [[nodiscard]] SourceRevision revision() const noexcept
    {
        return revision_;
    }

    /** @return Node spacing. */
    [[nodiscard]] sub0math::Length spacing() const noexcept
    {
        return spacing_;
    }

    /**
     * Reports whether a node is authored.
     *
     * @param i Node index along x, 0 to Nodes - 1.
     * @param j Node index along y, 0 to Nodes - 1.
     * @return True when the node has data.
     */
    [[nodiscard]] bool known(std::int32_t i, std::int32_t j) const noexcept
    {
        const auto n = index(i, j);
        return ((known_[n >> 6U] >> (n & 63U)) & 1U) != 0U;
    }

    /**
     * Local height of a node.
     *
     * @param i Node index along x.
     * @param j Node index along y.
     * @return Height above the page's height origin; meaningless for an unknown node.
     */
    [[nodiscard]] sub0math::Length nodeHeight(std::int32_t i, std::int32_t j) const noexcept
    {
        return toHeight(heights_[index(i, j)]);
    }

    /**
     * Height of the surface above a local position, on the page's fixed triangulation.
     *
     * @param x Local x, 0 to the page side.
     * @param y Local y, 0 to the page side.
     * @return The height, or not known outside the page or where the triangle under the point has an
     *         unknown node.
     */
    [[nodiscard]] HeightSample heightAt(sub0math::Length x, sub0math::Length y) const noexcept;

    /**
     * Heights of the surface at many local positions: the batch form of `heightAt`, structure-of-arrays in and
     * out, for extraction and for callers that sample a patch. Each sample is exactly what `heightAt` gives.
     *
     * @param xs Local x of each sample.
     * @param ys Local y of each sample; the same count as `xs`.
     * @param heights Output: the height of each sample; unspecified where the sample is not known. At least
     *        as many as `xs`.
     * @param known Output: 1 where the surface is defined, else 0 (outside the page, or under a triangle with an
     *        unknown node). At least as many as `xs`.
     * @return The number of samples that are known; 0 for mismatched or short buffers (nothing is written).
     */
    [[nodiscard]] std::size_t heightsAt(std::span<const sub0math::Length> xs, std::span<const sub0math::Length> ys,
                                        std::span<sub0math::Length> heights, std::span<std::uint8_t> known) const noexcept;

    /**
     * One row of node heights, for extraction loops that read many nodes in order.
     *
     * @param j Node row, 0 to Nodes - 1.
     * @return The `Nodes` quantised heights of the row in 0.1 mm steps (meaningless where a node is unknown).
     */
    [[nodiscard]] std::span<const std::uint16_t> nodeRow(std::int32_t j) const noexcept
    {
        return {heights_.data() + static_cast<std::size_t>(j) * static_cast<std::size_t>(Nodes), static_cast<std::size_t>(Nodes)};
    }

    /**
     * Bounds of a pyramid block.
     *
     * @param level Pyramid level, 0 (8 x 8 cells) to PyramidLevels - 1 (the whole page).
     * @param bi Block index along x, below 32 >> level.
     * @param bj Block index along y, below 32 >> level.
     * @return Its conservative height range.
     */
    [[nodiscard]] HeightRange blockRange(std::int32_t level, std::int32_t bi, std::int32_t bj) const noexcept
    {
        const std::int32_t side = 32 >> level;
        return pyramid_[static_cast<std::size_t>(level)][static_cast<std::size_t>(bj * side + bi)];
    }

    /**
     * Converts quantised steps to a local height.
     *
     * @param steps Height in 0.1 mm steps.
     * @return The local height.
     */
    [[nodiscard]] static constexpr sub0math::Length toHeight(std::uint16_t steps) noexcept
    {
        return sub0math::Length::fromSi(static_cast<float>(steps) * StepMetres);
    }

private:
    friend std::expected<HeightPage, Fault> buildHeightPage(const GridSpec& grid, const TileKey& key, SourceRevision revision,
                                                           std::span<const sub0math::f64::Length> heights,
                                                           std::span<const std::uint8_t> known);

    HeightPage(PageFrame frame, SourceRevision revision, sub0math::Length spacing)
        : frame_(frame), revision_(revision), spacing_(spacing), inverseSpacing_(1.0F / spacing.si())
    {
    }

    [[nodiscard]] static std::size_t index(std::int32_t i, std::int32_t j) noexcept
    {
        return static_cast<std::size_t>(j) * static_cast<std::size_t>(Nodes) + static_cast<std::size_t>(i);
    }

    // The lookup shared by `heightAt` and `heightsAt`: the arrays and the reciprocal spacing are arguments, so
    // a batch hoists them out of its loop.
    [[nodiscard]] static HeightSample lookup(const std::uint16_t* h, const std::uint64_t* mask, float inverseSpacing,
                                             sub0math::Length x, sub0math::Length y) noexcept;

    PageFrame frame_;
    SourceRevision revision_;
    sub0math::Length spacing_;
    float inverseSpacing_;                                        // nodes per metre
    std::vector<std::uint16_t> heights_;                          // Nodes x Nodes, row-major in j
    std::vector<std::uint64_t> known_;                            // one bit per node
    std::array<std::vector<HeightRange>, PyramidLevels> pyramid_; // level 0 = 32 x 32 blocks
};

// PERF: plain HeightSample ~1.7x faster than std::optional (GCC 16).
SUB0MATH_FORCE_INLINE HeightSample HeightPage::lookup(const std::uint16_t* h, const std::uint64_t* mask, float inverseSpacing,
                                                      sub0math::Length x, sub0math::Length y) noexcept
{
    // PERF: reciprocal spacing faster than dividing by the spacing.
    const float lx = x.si() * inverseSpacing;
    const float ly = y.si() * inverseSpacing;
    if (!(lx >= 0.0F && ly >= 0.0F && lx <= static_cast<float>(Cells) && ly <= static_cast<float>(Cells)))
    {
        return {sub0math::Length{}, false};
    }
    const std::int32_t i = std::min(static_cast<std::int32_t>(lx), Cells - 1);
    const std::int32_t j = std::min(static_cast<std::int32_t>(ly), Cells - 1);
    const float fx = lx - static_cast<float>(i);
    const float fy = ly - static_cast<float>(j);
    const auto bit = [mask](std::size_t n) { return ((mask[n >> 6U] >> (n & 63U)) & 1U) != 0U; };
    const std::size_t a = index(i, j);
    const std::size_t d = a + static_cast<std::size_t>(Nodes) + 1U;
    // Fixed diagonal (i, j)-(i+1, j+1): triangle a-b-d below it, a-d-c above.
    // PERF: two plain branches faster than a merged triangle select.
    if (fx >= fy)
    {
        const std::size_t b = a + 1U;
        if (!bit(a) || !bit(b) || !bit(d)) return {sub0math::Length{}, false};
        const float za = static_cast<float>(h[a]), zb = static_cast<float>(h[b]), zd = static_cast<float>(h[d]);
        return {sub0math::Length::fromSi((za + (zb - za) * fx + (zd - zb) * fy) * StepMetres), true};
    }
    const std::size_t c = a + static_cast<std::size_t>(Nodes);
    if (!bit(a) || !bit(d) || !bit(c)) return {sub0math::Length{}, false};
    const float za = static_cast<float>(h[a]), zd = static_cast<float>(h[d]), zc = static_cast<float>(h[c]);
    return {sub0math::Length::fromSi((za + (zd - zc) * fx + (zc - za) * fy) * StepMetres), true};
}

SUB0MATH_FORCE_INLINE HeightSample HeightPage::heightAt(sub0math::Length x, sub0math::Length y) const noexcept
{
    return lookup(heights_.data(), known_.data(), inverseSpacing_, x, y);
}

inline std::size_t HeightPage::heightsAt(std::span<const sub0math::Length> xs, std::span<const sub0math::Length> ys,
                                         std::span<sub0math::Length> heights, std::span<std::uint8_t> known) const noexcept
{
    const std::size_t n = xs.size();
    if (ys.size() != n || heights.size() < n || known.size() < n)
    {
        return 0;
    }
    const std::uint16_t* h = heights_.data();
    const std::uint64_t* mask = known_.data();
    const float inverse = inverseSpacing_;
    std::size_t count = 0;
    for (std::size_t k = 0; k < n; ++k)
    {
        const HeightSample s = lookup(h, mask, inverse, xs[k], ys[k]);
        heights[k] = s.height;
        known[k] = s.known ? 1U : 0U;
        count += s.known ? 1U : 0U;
    }
    return count;
}

/**
 * Builds a page from authored samples, validating them.
 *
 * The page's height origin is its lowest authored height, so stored heights are non-negative.
 *
 * @param grid The grid the page belongs to: 256 cells per page and node registration.
 * @param key The page.
 * @param revision Its content generation.
 * @param heights Nodes x Nodes heights in the terrain frame, row-major in j (row j holds i = 0..Nodes-1).
 * @param known Nodes x Nodes flags: nonzero where the height is authored. Unauthored heights are ignored.
 * @return The page, or an InvalidInput fault: a grid that is not 256-cell node-registered, wrong sample
 *         count, a non-finite authored height, a key outside the grid, or an authored height range over 6.5535 m (a page needs 32-bit heights).
 */
[[nodiscard]] std::expected<HeightPage, Fault> buildHeightPage(const GridSpec& grid, const TileKey& key, SourceRevision revision,
                                                              std::span<const sub0math::f64::Length> heights,
                                                              std::span<const std::uint8_t> known);

} // namespace sub0physics::terrain
