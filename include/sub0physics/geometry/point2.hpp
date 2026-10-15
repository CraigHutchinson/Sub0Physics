#pragma once

namespace sub0physics
{

/** A point or vector in a vertical plane: `x` along the ground, `y` up. Metres unless stated. */
struct Point2
{
    double x; ///< Horizontal coordinate.
    double y; ///< Vertical coordinate.

    /** Compares coordinates exactly. */
    friend constexpr bool operator==(const Point2&, const Point2&) noexcept = default;
};

} // namespace sub0physics
