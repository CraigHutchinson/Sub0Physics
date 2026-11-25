#pragma once

#include <cmath>
#include <limits>

namespace sub0physics::detail
{

// Some standard libraries do not yet support std::isfinite in constant evaluation.
[[nodiscard]] constexpr bool isFiniteParameter(double value) noexcept
{
    if consteval
    {
        constexpr double Limit = std::numeric_limits<double>::max();
        return value >= -Limit && value <= Limit;
    }
    else
    {
        return std::isfinite(value);
    }
}

} // namespace sub0physics::detail
