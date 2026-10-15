#pragma once

#include <initializer_list>
#include <sub0math/units.hpp>
#include <utility>
#include <vector>

#include "sub0physics/geometry/polyline_profile.hpp"

// Test shorthand: profile points written as (x, y) pairs in metres.

namespace test_profiles
{

inline sub0physics::ProfilePoint at(double x, double y)
{
    return {sub0math::metres(x), sub0math::metres(y)};
}

inline std::vector<sub0physics::ProfilePoint> metresPoints(std::initializer_list<std::pair<double, double>> points)
{
    std::vector<sub0physics::ProfilePoint> out;
    out.reserve(points.size());
    for (const auto& [x, y] : points)
    {
        out.push_back(at(x, y));
    }
    return out;
}

} // namespace test_profiles
