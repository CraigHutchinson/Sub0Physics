#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <sub0math/geometry/aabb.hpp>
#include <sub0math/geometry/oriented_box.hpp>
#include <sub0math/geometry/regular_grid.hpp>
#include <sub0math/units/quantities.hpp>

#include "sub0physics/fault.hpp"

namespace sub0physics::terrain
{

/** A wheel's pose on the ground plane in the terrain frame: where the axle projects, and the heading. */
struct PlanPose
{
    sub0math::f64::Length x;  ///< Forward [m].
    sub0math::f64::Length y;  ///< Left [m].
    sub0math::f64::Angle yaw; ///< Heading, counter-clockwise from +x [rad].
};

/**
 * The ground footprint a wheel needs terrain under: a rectangle around the pose, plus a lookbehind copy of
 * it shifted against the direction of travel.
 */
struct CorridorFootprint
{
    sub0math::f64::Length halfLength; ///< Half size along the heading; positive.
    sub0math::f64::Length halfWidth;  ///< Half size across the heading; positive.
    sub0math::f64::Length behind;     ///< Lookbehind against the displacement (along the start heading when stationary); not negative.
};

/** What `CorridorBuilder::create` needs: a pose pair, a footprint and how finely to split the motion. */
struct CorridorSpec
{
    PlanPose from;                    ///< Pose at the start of the interval.
    PlanPose to;                      ///< Pose at its end; position and heading interpolate linearly from `from`.
    CorridorFootprint footprint;      ///< What is swept.
    std::uint32_t subIntervals = 8;   ///< Boxes to split the motion into, 1 to CorridorBuilder::MaxSubIntervals.
};

/** The oriented box of one sub-interval: a rectangle in the terrain plane with its own heading. */
using CorridorBox = sub0math::OrientedBox<sub0math::f64::Length, 2>;

/** Whether `CorridorBuilder::rasterise` produced a mask, and if not why. */
enum class RasterStatus : std::uint8_t
{
    Complete,       ///< The mask holds every touched cell.
    OutsideGrid,    ///< The corridor misses the grid, or its bounds are not finite; the mask is empty.
    BufferTooSmall  ///< The caller's words cannot hold the range; nothing was written.
};

/**
 * A rectangle of cells with one bit each, in the caller's storage: bit `(j - lo.j) * width + (i - lo.i)` of the
 * words, least significant bit first.
 */
struct CellMask
{
    sub0math::GridRange<2> range;            ///< Cells the bits describe; empty when there are none.
    std::span<const std::uint64_t> words;    ///< The bits.

    /**
     * Reads one cell.
     *
     * @param i Cell index along x.
     * @param j Cell index along y.
     * @return True when the cell is set; false for a cell outside the range.
     */
    [[nodiscard]] bool contains(std::int32_t i, std::int32_t j) const noexcept
    {
        if (i < range.lo[0] || i > range.hi[0] || j < range.lo[1] || j > range.hi[1])
        {
            return false;
        }
        const auto width = static_cast<std::size_t>(range.hi[0] - range.lo[0] + 1);
        const auto n = static_cast<std::size_t>(j - range.lo[1]) * width + static_cast<std::size_t>(i - range.lo[0]);
        return ((words[n >> 6U] >> (n & 63U)) & 1U) != 0U;
    }
};

/** Result of `CorridorBuilder::rasterise`. */
struct Raster
{
    RasterStatus status; ///< Whether `mask` is valid.
    CellMask mask;       ///< The covered cells, in the caller's storage.
    std::size_t cells;   ///< Cells set in the mask.
};

/**
 * Builds the terrain a wheel needs for one interval of motion, as cells of a regular grid.
 *
 * The motion from a pose pair (position and heading interpolated linearly, as a solver step does) is split
 * into sub-intervals. For each one `create` makes one oriented box that encloses the footprint at every pose
 * in the sub-interval, so the corridor is never an axis-aligned rectangle around the whole path: a turning or
 * diagonal wheel gets a thin corridor, not its bounding box. `rasterise` then marks every cell a box touches,
 * by the separating-axis test against the cell.
 *
 * **Why the boxes contain the motion.** For a sub-interval, the box is aligned with the heading at its middle.
 * The footprint (with its lookbehind copy) has a centre offset and half extents in that frame. Over the
 * sub-interval the heading differs from the middle by at most `d = |dyaw| / (2k)`: a point at local
 * coordinates (qx, qy) then moves by at most `|qx| d^2/2 + |qy| d` along x and `|qy| d^2/2 + |qx| d` along y
 * (from `1 - cos e <= e^2/2` and `|sin e| <= |e|`). The position moves by at most half the sub-interval's
 * displacement, taken per axis in the box frame. The half extents grow by those amounts, and by
 * `RoundingSlack` for floating-point rounding. The bound holds for any rotation, however large.
 *
 * Coverage is closed: a cell the footprint only touches is included, as in `sub0math::RegularGrid`.
 */
class CorridorBuilder
{
public:
    static constexpr std::uint32_t MaxSubIntervals = 64;              ///< Most boxes one corridor may have.
    static constexpr double PoseLimit = 1e6;                          ///< Largest pose or footprint magnitude [m].
    static constexpr sub0math::f64::Length RoundingSlack = sub0math::f64::Length::fromSi(1e-9); ///< Added to every half extent.

    /**
     * Validates a spec and builds the sub-interval boxes.
     *
     * @param spec The pose pair, footprint and number of sub-intervals.
     * @return The corridor, or an InvalidInput fault: a non-finite or huge pose, a non-positive half size,
     *         a negative lookbehind, or a sub-interval count outside 1 to MaxSubIntervals.
     */
    [[nodiscard]] static std::expected<CorridorBuilder, Fault> create(const CorridorSpec& spec) noexcept;

    /** @return The boxes in time order; at least one. */
    [[nodiscard]] std::span<const CorridorBox> boxes() const noexcept
    {
        return {boxes_.data(), count_};
    }

    /** @return The smallest axis-aligned box around all the boxes. */
    [[nodiscard]] sub0math::Aabb<sub0math::f64::Length, 2> bounds() const noexcept
    {
        return bounds_;
    }

    /**
     * Marks the cells of a grid that the corridor touches.
     *
     * @param grid The cells to mark, such as `GridSpec::cellGrid(level)`.
     * @param words Storage for the mask: at least `wordsFor` the range of `bounds()` on the grid.
     * @return The mask over the grid's cells that the bounds touch (cleared first, then set), or why none.
     */
    [[nodiscard]] Raster rasterise(const sub0math::RegularGrid<sub0math::f64::Length, 2>& grid,
                                   std::span<std::uint64_t> words) const noexcept;

    /**
     * Finds how many words a mask over a range needs.
     *
     * @param range The cells; may be empty.
     * @return Number of 64-bit words, at least 1 when the range is not empty.
     */
    [[nodiscard]] static constexpr std::size_t wordsFor(const sub0math::GridRange<2>& range) noexcept
    {
        return static_cast<std::size_t>((range.count() + 63) / 64);
    }

private:
    CorridorBuilder() noexcept = default;

    std::array<CorridorBox, MaxSubIntervals> boxes_{};
    std::size_t count_ = 0;
    sub0math::Aabb<sub0math::f64::Length, 2> bounds_{};
};

} // namespace sub0physics::terrain
