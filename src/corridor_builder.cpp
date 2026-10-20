#include "sub0physics/terrain/corridor_builder.hpp"

#include <algorithm>
#include <cmath>

namespace sub0physics::terrain
{

namespace
{

using sub0math::f64::Length;

bool within(double v, double limit) noexcept
{
    return v >= -limit && v <= limit; // NaN fails
}

bool validPose(const PlanPose& p) noexcept
{
    return within(p.x.si(), CorridorBuilder::PoseLimit) && within(p.y.si(), CorridorBuilder::PoseLimit) &&
           within(p.yaw.si(), CorridorBuilder::PoseLimit);
}

// A vector expressed in the frame whose first axis is (c, s).
struct Frame
{
    double c, s;

    [[nodiscard]] std::array<double, 2> local(double x, double y) const noexcept
    {
        return {c * x + s * y, -s * x + c * y};
    }
};

} // namespace

std::expected<CorridorBuilder, Fault> CorridorBuilder::create(const CorridorSpec& spec) noexcept
{
    if (spec.subIntervals < 1U || spec.subIntervals > MaxSubIntervals)
    {
        return invalidInput("a corridor needs 1 to 64 sub-intervals");
    }
    if (!validPose(spec.from) || !validPose(spec.to))
    {
        return invalidInput("corridor poses must be finite and within 1e6 m");
    }
    const CorridorFootprint& f = spec.footprint;
    if (!(f.halfLength.si() > 0.0 && f.halfLength.si() <= PoseLimit) || !(f.halfWidth.si() > 0.0 && f.halfWidth.si() <= PoseLimit) ||
        !(f.behind.si() >= 0.0 && f.behind.si() <= PoseLimit))
    {
        return invalidInput("corridor footprint needs positive half sizes and a lookbehind of at least zero");
    }

    const double ax = spec.from.x.si(), ay = spec.from.y.si(), ayaw = spec.from.yaw.si();
    const double dx = spec.to.x.si() - ax, dy = spec.to.y.si() - ay;
    const double dyaw = spec.to.yaw.si() - ayaw;
    const double k = static_cast<double>(spec.subIntervals);
    // Lookbehind opposite the displacement; along the start heading when stationary.
    const double distance = std::hypot(dx, dy);
    const double ux = distance != 0.0 ? dx / distance : std::cos(ayaw);
    const double uy = distance != 0.0 ? dy / distance : std::sin(ayaw);
    const double wx = -f.behind.si() * ux, wy = -f.behind.si() * uy;
    const double hl = f.halfLength.si(), hw = f.halfWidth.si();
    const double delta = std::abs(dyaw) / (2.0 * k); // largest heading departure from a sub-interval's middle
    const double rotX = hl * delta * delta / 2.0 + hw * delta;
    const double rotY = hw * delta * delta / 2.0 + hl * delta;

    CorridorBuilder out;
    out.count_ = spec.subIntervals;
    for (std::size_t m = 0; m < out.count_; ++m)
    {
        const double s = (static_cast<double>(m) + 0.5) / k;
        const double yaw = ayaw + dyaw * s;
        const Frame frame{std::cos(yaw), std::sin(yaw)};
        const auto w = frame.local(wx, wy);
        const auto d = frame.local(dx, dy);
        const double slack = RoundingSlack.si();
        const double halfX = hl + rotX + std::abs(d[0]) / (2.0 * k) + std::abs(w[0]) / 2.0 + slack;
        const double halfY = hw + rotY + std::abs(d[1]) / (2.0 * k) + std::abs(w[1]) / 2.0 + slack;
        // Centre: the pose at the middle of the sub-interval, moved half-way to the lookbehind copy.
        const double lx = w[0] / 2.0, ly = w[1] / 2.0;
        CorridorBox& box = out.boxes_[m];
        box.centre = {Length::fromSi(ax + dx * s + frame.c * lx - frame.s * ly), Length::fromSi(ay + dy * s + frame.s * lx + frame.c * ly)};
        box.halfExtents = {Length::fromSi(halfX), Length::fromSi(halfY)};
        box.axes.columns[0] = {frame.c, frame.s};
        box.axes.columns[1] = {-frame.s, frame.c};
        const auto bb = boundingBox(box);
        if (m == 0)
        {
            out.bounds_ = bb;
        }
        else
        {
            out.bounds_.lo = {std::min(out.bounds_.lo.x, bb.lo.x), std::min(out.bounds_.lo.y, bb.lo.y)};
            out.bounds_.hi = {std::max(out.bounds_.hi.x, bb.hi.x), std::max(out.bounds_.hi.y, bb.hi.y)};
        }
    }
    return out;
}

Raster CorridorBuilder::rasterise(const sub0math::RegularGrid<Length, 2>& grid, std::span<std::uint64_t> words) const noexcept
{
    const sub0math::GridRange<2> range = grid.cellsTouching(bounds_);
    if (range.empty())
    {
        return {RasterStatus::OutsideGrid, {range, {}}, 0U};
    }
    const std::size_t needed = wordsFor(range);
    if (words.size() < needed)
    {
        return {RasterStatus::BufferTooSmall, {range, {}}, 0U};
    }
    std::fill_n(words.begin(), static_cast<std::ptrdiff_t>(needed), std::uint64_t{0});
    const auto width = static_cast<std::size_t>(range.hi[0] - range.lo[0] + 1);
    const Length half = grid.spacing / 2.0;
    std::size_t cells = 0;
    for (const CorridorBox& box : boxes())
    {
        const sub0math::AabbOverlapTester<Length, 2> overlaps(box);
        const sub0math::GridRange<2> candidates = grid.cellsTouching(boundingBox(box));
        for (std::int32_t j = std::max(candidates.lo[1], range.lo[1]); j <= std::min(candidates.hi[1], range.hi[1]); ++j)
        {
            for (std::int32_t i = std::max(candidates.lo[0], range.lo[0]); i <= std::min(candidates.hi[0], range.hi[0]); ++i)
            {
                const auto corner = grid.nodeAt({i, j});
                if (!overlaps.overlapsCentred({corner.x + half, corner.y + half}, {half, half}))
                {
                    continue;
                }
                const auto n = static_cast<std::size_t>(j - range.lo[1]) * width + static_cast<std::size_t>(i - range.lo[0]);
                const std::uint64_t bit = std::uint64_t{1} << (n & 63U);
                std::uint64_t& word = words[n >> 6U];
                cells += (word & bit) == 0U ? 1U : 0U;
                word |= bit;
            }
        }
    }
    return {RasterStatus::Complete, {range, words.first(needed)}, cells};
}

std::expected<TightCorridor, Fault> tightCorridor(CorridorSpec spec, const sub0math::RegularGrid<sub0math::f64::Length, 2>& grid,
                                                 std::span<std::uint64_t> words, double relativeGain) noexcept
{
    if (!(relativeGain >= 0.0 && relativeGain <= 1.0))
    {
        return invalidInput("the corridor gain target must be from 0 to 1");
    }
    const auto cellsFor = [&](std::uint32_t k) -> std::expected<std::size_t, Fault> {
        spec.subIntervals = k;
        const auto corridor = CorridorBuilder::create(spec);
        if (!corridor)
        {
            return std::unexpected(corridor.error());
        }
        const Raster raster = corridor->rasterise(grid, words);
        if (raster.status == RasterStatus::BufferTooSmall)
        {
            return invalidInput("the corridor mask scratch is too small");
        }
        return raster.cells; // zero when the corridor misses the grid
    };
    std::uint32_t k = 1;
    auto cells = cellsFor(k);
    if (!cells)
    {
        return std::unexpected(cells.error());
    }
    while (k < CorridorBuilder::MaxSubIntervals)
    {
        const auto finer = cellsFor(k * 2U);
        if (!finer)
        {
            return std::unexpected(finer.error());
        }
        if (static_cast<double>(*finer) > static_cast<double>(*cells) * (1.0 - relativeGain))
        {
            break;
        }
        k *= 2U;
        cells = finer;
    }
    return TightCorridor{k, *cells};
}

} // namespace sub0physics::terrain
