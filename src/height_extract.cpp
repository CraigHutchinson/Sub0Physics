#include "sub0physics/terrain/height_extract.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>
#include <sub0math/geometry/ray.hpp>

namespace sub0physics::terrain
{

using sub0math::Length;
using sub0math::Vec3f;

struct HeightExtractWriter
{
    static void prepare(HeightExtract& e, std::size_t nodes, std::size_t cells) noexcept
    {
        (void)e.nodeZ_.resize(nodes);
        (void)e.nodeKnown_.resize(nodes);
        (void)e.cellState_.resize(cells);
        (void)e.cellTriangle_.resize(cells);
        e.triangles_.clear();
        e.triangleNormal_.clear();
        e.triangleCell_.clear();
    }

    static void clear(HeightExtract& e) noexcept
    {
        e.nodeZ_.clear();
        e.nodeKnown_.clear();
        e.cellState_.clear();
        e.cellTriangle_.clear();
        e.triangles_.clear();
        e.triangleNormal_.clear();
        e.triangleCell_.clear();
        e.cellsX_ = e.cellsY_ = 0;
    }
};

HeightExtract::HeightExtract(const HeightExtractLimits& limits)
{
    const std::size_t cells = limits.triangles / 2U;
    nodeZ_.reserve(limits.nodes);
    nodeKnown_.reserve(limits.nodes);
    cellState_.reserve(cells);
    cellTriangle_.reserve(cells);
    triangles_.reserve(limits.triangles);
    triangleNormal_.reserve(limits.triangles);
    triangleCell_.reserve(limits.triangles);
    maskWords_.reserve(cells / 64U + 2U);
    (void)maskWords_.resize(maskWords_.capacity());
}

HeightResult extractHeights(const GridSpec& grid, const PageWindow& pages, const HeightRequest& request,
                            HeightExtract& out) noexcept
{
    HeightExtractWriter::clear(out);
    constexpr std::int32_t N = HeightPage::Cells;
    const std::uint32_t s = request.stride;
    if (grid.cellsPerSide() != N || grid.registration() != Registration::NodeIsPoint || request.level >= grid.levels() ||
        pages.level != request.level || request.cells.empty() || s == 0U || s > 256U || !std::has_single_bit(s))
    {
        return {HeightStatus::InvalidRequest};
    }
    const std::int32_t totalX = grid.pageCount(request.level, 0) * N;
    const std::int32_t totalY = grid.pageCount(request.level, 1) * N;
    const auto& r = request.cells;
    if (r.lo[0] < 0 || r.lo[1] < 0 || r.hi[0] >= totalX || r.hi[1] >= totalY)
    {
        return {HeightStatus::InvalidRequest};
    }
    const auto si = static_cast<std::int32_t>(s);
    const std::int32_t lo0 = r.lo[0] / si * si;
    const std::int32_t lo1 = r.lo[1] / si * si;
    const std::int32_t cx = (r.hi[0] / si * si + si - lo0) / si;
    const std::int32_t cy = (r.hi[1] / si * si + si - lo1) / si;
    const auto nodes = static_cast<std::size_t>(cx + 1) * static_cast<std::size_t>(cy + 1);
    const auto cells = static_cast<std::size_t>(cx) * static_cast<std::size_t>(cy);
    const HeightExtractLimits cap = out.limits();
    if (nodes > cap.nodes || cells * 2U > cap.triangles)
    {
        return {HeightStatus::Budget};
    }
    HeightExtractWriter::prepare(out, nodes, cells);

    const double fine = grid.spacing(request.level).si();
    out.stride_ = s;
    out.firstX_ = lo0;
    out.firstY_ = lo1;
    out.cellsX_ = cx;
    out.cellsY_ = cy;
    out.spacing_ = Length::fromSi(static_cast<float>(fine * static_cast<double>(s)));
    const float h = out.spacing_.si();

    // Height origin: the first page present, row by row, so it does not depend on the mask.
    double originZ = 0.0;
    for (const HeightPage* page : pages.data)
    {
        if (page != nullptr)
        {
            originZ = page->frame().origin.z.si();
            break;
        }
    }
    out.origin_ = {grid.anchor().x + grid.spacing(request.level) * static_cast<double>(lo0),
                   grid.anchor().y + grid.spacing(request.level) * static_cast<double>(lo1), sub0math::f64::Length::fromSi(originZ)};
    const auto dzOf = [originZ](const HeightPage& page) {
        return static_cast<float>(page.frame().origin.z.si() - originZ);
    };

    // Nodes. A node on a shared page edge is read from the page on the positive side, or, where that page is
    // absent, from the page below it. Each page quantises against its own height origin, so the two copies of
    // an edge node can differ by one step (0.1 mm): triangles take their corner heights from these node
    // values, so neighbouring cells meet exactly at a seam.
    Length* nodeZ = out.nodeZ_.data();
    std::uint8_t* nodeKnown = out.nodeKnown_.data();
    const auto nodeSample = [&](std::int32_t I, std::int32_t J, Length& z, std::uint8_t& known) {
        z = Length{};
        known = 0U;
        if (I > totalX || J > totalY)
        {
            return;
        }
        const std::int32_t pi = std::min(I, totalX - 1) / N;
        const std::int32_t pj = std::min(J, totalY - 1) / N;
        const HeightPage* page = pages.at(pi, pj);
        std::int32_t li = I - pi * N;
        std::int32_t lj = J - pj * N;
        if (page == nullptr)
        {
            // The page below on an edge node: local index N of the previous page.
            const std::int32_t qi = (I % N == 0 && I > 0) ? pi - 1 : pi;
            const std::int32_t qj = (J % N == 0 && J > 0) ? pj - 1 : pj;
            if (qi == pi && qj == pj)
            {
                return;
            }
            page = pages.at(qi, qj);
            if (page == nullptr)
            {
                return;
            }
            li = I - qi * N;
            lj = J - qj * N;
        }
        if (page->known(li, lj))
        {
            known = 1U;
            z = Length::fromSi(page->nodeHeight(li, lj).si() + dzOf(*page));
        }
    };
    for (std::int32_t b = 0; b <= cy; ++b)
    {
        for (std::int32_t a = 0; a <= cx; ++a)
        {
            const auto at = static_cast<std::size_t>(b) * static_cast<std::size_t>(cx + 1) + static_cast<std::size_t>(a);
            nodeSample(lo0 + a * si, lo1 + b * si, nodeZ[at], nodeKnown[at]);
        }
    }

    // Cells.
    const CellMask* mask = request.mask;
    const auto requested = [&](std::int32_t I, std::int32_t J) {
        if (mask == nullptr)
        {
            return true;
        }
        const std::int32_t i1 = std::min(mask->range.hi[0], I + si - 1);
        const std::int32_t j1 = std::min(mask->range.hi[1], J + si - 1);
        for (std::int32_t j = std::max(mask->range.lo[1], J); j <= j1; ++j)
        {
            for (std::int32_t i = std::max(mask->range.lo[0], I); i <= i1; ++i)
            {
                if (mask->contains(i, j))
                {
                    return true;
                }
            }
        }
        return false;
    };
    HeightResult result{HeightStatus::Complete};
    CellState* state = out.cellState_.data();
    std::uint32_t* cellTri = out.cellTriangle_.data();
    for (std::int32_t b = 0; b < cy; ++b)
    {
        for (std::int32_t a = 0; a < cx; ++a)
        {
            const std::int32_t I = lo0 + a * si;
            const std::int32_t J = lo1 + b * si;
            const auto cell = static_cast<std::size_t>(b) * static_cast<std::size_t>(cx) + static_cast<std::size_t>(a);
            cellTri[cell] = HeightExtract::NoTriangle;
            if (!requested(I, J))
            {
                state[cell] = CellState::NotRequested;
                continue;
            }
            const HeightPage* page = (I < totalX && J < totalY) ? pages.at(I / N, J / N) : nullptr;
            bool deliver = page != nullptr;
            const std::int32_t li = I % N;
            const std::int32_t lj = J % N;
            if (deliver)
            {
                // Every finer node under the cell must be authored: a coarse cell never interpolates across nodata.
                for (std::int32_t v = 0; v <= si && deliver; ++v)
                {
                    for (std::int32_t u = 0; u <= si && deliver; ++u)
                    {
                        deliver = page->known(li + u, lj + v);
                    }
                }
            }
            if (!deliver)
            {
                state[cell] = CellState::Unknown;
                ++result.unknownCells;
                continue;
            }
            const auto nodeAt = [&](std::int32_t u, std::int32_t v) {
                return static_cast<std::size_t>(b + v) * static_cast<std::size_t>(cx + 1) + static_cast<std::size_t>(a + u);
            };
            if (nodeKnown[nodeAt(0, 0)] == 0U || nodeKnown[nodeAt(1, 0)] == 0U || nodeKnown[nodeAt(0, 1)] == 0U ||
                nodeKnown[nodeAt(1, 1)] == 0U)
            {
                state[cell] = CellState::Unknown;
                ++result.unknownCells;
                continue;
            }
            state[cell] = CellState::Delivered;
            ++result.deliveredCells;
            const float za = nodeZ[nodeAt(0, 0)].si(), zb = nodeZ[nodeAt(1, 0)].si();
            const float zc = nodeZ[nodeAt(0, 1)].si(), zd = nodeZ[nodeAt(1, 1)].si();
            const float x0 = static_cast<float>(a) * h, x1 = static_cast<float>(a + 1) * h;
            const float y0 = static_cast<float>(b) * h, y1 = static_cast<float>(b + 1) * h;
            using P = sub0math::Point3<Length>;
            const P pa{Length::fromSi(x0), Length::fromSi(y0), Length::fromSi(za)};
            const P pb{Length::fromSi(x1), Length::fromSi(y0), Length::fromSi(zb)};
            const P pc{Length::fromSi(x0), Length::fromSi(y1), Length::fromSi(zc)};
            const P pd{Length::fromSi(x1), Length::fromSi(y1), Length::fromSi(zd)};
            cellTri[cell] = static_cast<std::uint32_t>(out.triangles_.size());
            // Capacities were checked: two triangles per cell at most.
            (void)out.triangles_.pushBack({pa, pb, pd});
            (void)out.triangles_.pushBack({pa, pd, pc});
            // Normals from the cross products of the same triangles, written out: (b-a)x(d-a), (d-a)x(c-a).
            const auto unit = [](float nx, float ny, float nz) {
                const float inv = 1.0F / std::sqrt(nx * nx + ny * ny + nz * nz);
                return Vec3f{nx * inv, ny * inv, nz * inv};
            };
            (void)out.triangleNormal_.pushBack(unit(-(zb - za), zb - zd, h));
            (void)out.triangleNormal_.pushBack(unit(zc - zd, za - zc, h));
            (void)out.triangleCell_.pushBack(static_cast<std::uint32_t>(cell));
            (void)out.triangleCell_.pushBack(static_cast<std::uint32_t>(cell));
        }
    }
    if (result.unknownCells != 0U)
    {
        result.status = HeightStatus::Partial;
    }
    return result;
}

HeightResult extractCorridor(const GridSpec& grid, const CorridorBuilder& corridor, const PageWindow& pages, std::uint8_t level,
                             std::uint32_t stride, HeightExtract& out) noexcept
{
    if (level >= grid.levels())
    {
        return {HeightStatus::InvalidRequest};
    }
    const Raster raster = corridor.rasterise(grid.cellGrid(level), out.maskWords());
    if (raster.status != RasterStatus::Complete)
    {
        HeightExtractWriter::clear(out);
        return {HeightStatus::InvalidRequest};
    }
    HeightRequest request;
    request.level = level;
    request.cells = raster.mask.range;
    request.stride = stride;
    request.mask = &raster.mask;
    return extractHeights(grid, pages, request, out);
}

namespace
{

RayResult hitOf(const HeightExtract& e, const sub0math::Ray<Length, 3>& ray, std::size_t tri, const sub0math::RayHit<Length>& hit)
{
    RayResult r;
    r.hit = true;
    r.distance = hit.distance;
    Vec3f n = e.triangleNormal()[tri];
    if (n.x * ray.direction.x + n.y * ray.direction.y + n.z * ray.direction.z > 0.0F)
    {
        n = {-n.x, -n.y, -n.z};
    }
    const auto at = ray.at(hit.distance);
    r.witness = {{sub0math::f64::Length::fromSi(static_cast<double>(at.x.si())), sub0math::f64::Length::fromSi(static_cast<double>(at.y.si())),
                  sub0math::f64::Length::fromSi(static_cast<double>(at.z.si()))},
                 {static_cast<double>(n.x), static_cast<double>(n.y), static_cast<double>(n.z)},
                 WitnessKind::Face,
                 tri};
    return r;
}

} // namespace

RayResult castRayBruteForce(const HeightExtract& e, const sub0math::Ray<Length, 3>& ray, Length maxDistance) noexcept
{
    RayResult best;
    const auto tris = e.triangles();
    for (std::size_t t = 0; t < tris.size(); ++t)
    {
        const auto hit = sub0math::intersect(ray, tris[t], maxDistance);
        if (hit && (!best.hit || hit->distance < best.distance))
        {
            best = hitOf(e, ray, t, *hit);
        }
    }
    return best;
}

RayResult castRay(const HeightExtract& e, const sub0math::Ray<Length, 3>& ray, Length maxDistance) noexcept
{
    const auto [cx, cy] = e.cells();
    const float h = e.spacing().si();
    if (cx <= 0 || cy <= 0 || !(h > 0.0F))
    {
        return {};
    }
    const float o[2] = {ray.origin.x.si(), ray.origin.y.si()};
    const float d[2] = {ray.direction.x, ray.direction.y};
    const float extent[2] = {static_cast<float>(cx) * h, static_cast<float>(cy) * h};
    const float maxD = maxDistance.si();
    float t0 = 0.0F, t1 = maxD;
    for (int axis = 0; axis < 2; ++axis)
    {
        if (std::abs(d[axis]) < 1e-12F)
        {
            if (o[axis] < 0.0F || o[axis] > extent[axis])
            {
                return {};
            }
            continue;
        }
        float ta = (0.0F - o[axis]) / d[axis];
        float tb = (extent[axis] - o[axis]) / d[axis];
        if (ta > tb)
        {
            std::swap(ta, tb);
        }
        t0 = std::max(t0, ta);
        t1 = std::min(t1, tb);
    }
    if (!(t0 <= t1))
    {
        return {};
    }
    // Amanatides-Woo walk of the cells from the entry point.
    std::int32_t c[2];
    float tMax[2], tDelta[2];
    std::int32_t step[2];
    const std::int32_t count[2] = {cx, cy};
    for (int axis = 0; axis < 2; ++axis)
    {
        const float p = o[axis] + d[axis] * t0;
        c[axis] = std::clamp(static_cast<std::int32_t>(std::floor(p / h)), 0, count[axis] - 1);
        if (std::abs(d[axis]) < 1e-12F)
        {
            step[axis] = 0;
            tMax[axis] = std::numeric_limits<float>::infinity();
            tDelta[axis] = std::numeric_limits<float>::infinity();
        }
        else
        {
            step[axis] = d[axis] > 0.0F ? 1 : -1;
            const float edge = static_cast<float>(c[axis] + (step[axis] > 0 ? 1 : 0)) * h;
            tMax[axis] = (edge - o[axis]) / d[axis];
            tDelta[axis] = h / std::abs(d[axis]);
        }
    }
    RayResult best;
    const auto tris = e.triangles();
    const auto cellTri = e.cellTriangle();
    // A hit can lie on a cell edge, so keep walking while the next cell starts no later than the best hit.
    float enter = t0;
    for (;;)
    {
        if (best.hit && enter > best.distance.si() + 1e-6F)
        {
            break;
        }
        const auto cell = static_cast<std::size_t>(c[1]) * static_cast<std::size_t>(cx) + static_cast<std::size_t>(c[0]);
        const std::uint32_t first = cellTri[cell];
        if (first != HeightExtract::NoTriangle)
        {
            for (std::uint32_t t = first; t < first + 2U; ++t)
            {
                const auto hit = sub0math::intersect(ray, tris[t], maxDistance);
                if (hit && (!best.hit || hit->distance < best.distance))
                {
                    best = hitOf(e, ray, t, *hit);
                }
            }
        }
        const int axis = tMax[0] < tMax[1] ? 0 : 1;
        if (!(tMax[axis] <= t1) || step[axis] == 0)
        {
            break;
        }
        enter = tMax[axis];
        c[axis] += step[axis];
        if (c[axis] < 0 || c[axis] >= count[axis])
        {
            break;
        }
        tMax[axis] += tDelta[axis];
    }
    return best;
}

} // namespace sub0physics::terrain
