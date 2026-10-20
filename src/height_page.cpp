#include "sub0physics/terrain/height_page.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sub0physics::terrain
{

using sub0math::Length;

std::expected<HeightPage, Fault> buildHeightPage(const GridSpec& grid, const TileKey& key, SourceRevision revision,
                                                 std::span<const sub0math::f64::Length> heights,
                                                 std::span<const std::uint8_t> known)
{
    constexpr auto count = static_cast<std::size_t>(HeightPage::Nodes) * static_cast<std::size_t>(HeightPage::Nodes);
    if (heights.size() != count || known.size() != count)
    {
        return std::unexpected(Fault{FaultKind::InvalidInput, "a height page needs 257 x 257 heights and flags"});
    }
    // A GridSpec can only come from GridSpec::create, so the grid itself is valid; the key is not.
    if (key.level >= grid.levels())
    {
        return std::unexpected(Fault{FaultKind::InvalidInput, "page level is outside the grid"});
    }

    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    for (std::size_t n = 0; n < count; ++n)
    {
        if (known[n] == 0U) continue;
        const double z = heights[n].si();
        if (!std::isfinite(z))
        {
            return std::unexpected(Fault{FaultKind::InvalidInput, "authored height is not finite"});
        }
        lo = std::min(lo, z);
        hi = std::max(hi, z);
    }
    const double step = static_cast<double>(HeightPage::StepMetres);
    if (lo <= hi && (hi - lo) / step > 65535.0)
    {
        return std::unexpected(Fault{FaultKind::InvalidInput, "page height range exceeds 16-bit heights (6.5535 m)"});
    }
    const double origin = lo <= hi ? lo : 0.0;

    auto corner = grid.pageCorner(key);
    corner.z = sub0math::f64::Length::fromSi(origin);
    HeightPage page(PageFrame{key, corner}, revision, static_cast<Length>(grid.spacing(key.level)));
    page.heights_.assign(count, 0U);
    page.known_.assign((count + 63U) / 64U, 0U);
    for (std::size_t n = 0; n < count; ++n)
    {
        if (known[n] == 0U) continue;
        page.heights_[n] = static_cast<std::uint16_t>(std::lround((heights[n].si() - origin) / step));
        page.known_[n >> 6U] |= std::uint64_t{1} << (n & 63U);
    }

    // Pyramid level 0: each 8 x 8-cell block spans 9 x 9 nodes (its shared edges included), so a block's
    // range bounds every triangle inside it.
    constexpr std::int32_t blocks = HeightPage::Cells / HeightPage::BlockCells;
    auto& base = page.pyramid_[0];
    base.assign(static_cast<std::size_t>(blocks * blocks), HeightRange{0xFFFF, 0, false, false});
    for (std::int32_t bj = 0; bj < blocks; ++bj)
    {
        for (std::int32_t bi = 0; bi < blocks; ++bi)
        {
            HeightRange r{0xFFFF, 0, false, false};
            for (std::int32_t j = bj * HeightPage::BlockCells; j <= (bj + 1) * HeightPage::BlockCells; ++j)
            {
                for (std::int32_t i = bi * HeightPage::BlockCells; i <= (bi + 1) * HeightPage::BlockCells; ++i)
                {
                    if (!page.known(i, j))
                    {
                        r.anyUnknown = true;
                        continue;
                    }
                    const std::uint16_t q = page.heights_[HeightPage::index(i, j)];
                    r.lo = std::min(r.lo, q);
                    r.hi = std::max(r.hi, q);
                    r.anyKnown = true;
                }
            }
            base[static_cast<std::size_t>(bj * blocks + bi)] = r;
        }
    }
    // Coarser levels merge 2 x 2 children.
    for (std::int32_t level = 1; level < HeightPage::PyramidLevels; ++level)
    {
        const std::int32_t side = blocks >> level;
        const std::int32_t childSide = side * 2;
        const auto& children = page.pyramid_[static_cast<std::size_t>(level - 1)];
        auto& out = page.pyramid_[static_cast<std::size_t>(level)];
        out.assign(static_cast<std::size_t>(side * side), HeightRange{0xFFFF, 0, false, false});
        for (std::int32_t bj = 0; bj < side; ++bj)
        {
            for (std::int32_t bi = 0; bi < side; ++bi)
            {
                HeightRange r{0xFFFF, 0, false, false};
                for (std::int32_t cj = 0; cj < 2; ++cj)
                {
                    for (std::int32_t ci = 0; ci < 2; ++ci)
                    {
                        const auto& c = children[static_cast<std::size_t>((bj * 2 + cj) * childSide + bi * 2 + ci)];
                        r.anyUnknown = r.anyUnknown || c.anyUnknown;
                        if (!c.anyKnown) continue;
                        r.lo = std::min(r.lo, c.lo);
                        r.hi = std::max(r.hi, c.hi);
                        r.anyKnown = true;
                    }
                }
                out[static_cast<std::size_t>(bj * side + bi)] = r;
            }
        }
    }
    return page;
}

} // namespace sub0physics::terrain
