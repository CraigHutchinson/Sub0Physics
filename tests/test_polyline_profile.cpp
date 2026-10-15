#include <cmath>
#include <doctest.h>
#include <limits>
#include <vector>

#include "sub0physics/geometry/polyline_profile.hpp"
#include "test_profiles.hpp"

using namespace sub0physics;
using sub0math::metres;
using test_profiles::metresPoints;

namespace
{

std::string_view rejection(std::vector<ProfilePoint> points)
{
    const auto profile = PolylineProfile::create(std::move(points));
    REQUIRE_FALSE(profile);
    CHECK(profile.error().kind == FaultKind::InvalidInput);
    return profile.error().detail;
}

} // namespace

TEST_CASE("Profile validation rejects malformed polylines")
{
    CHECK(rejection(metresPoints({{0.0, 0.0}})).find("two points") != std::string_view::npos);
    CHECK(rejection(metresPoints({{0.0, 0.0}, {std::numeric_limits<double>::quiet_NaN(), 0.0}})).find("finite") != std::string_view::npos);
    CHECK(rejection(metresPoints({{0.0, 0.0}, {1e101, 0.0}})).find("limit") != std::string_view::npos);
    CHECK(rejection(metresPoints({{0.0, 0.0}, {1.0, 0.0}, {0.5, 0.0}})).find("overhang") != std::string_view::npos);
    CHECK(rejection(metresPoints({{0.0, 0.0}, {1e-7, 0.0}})).find("micrometre") != std::string_view::npos);
    CHECK(rejection(metresPoints({{0.0, 0.0}, {0.0, 1.0}, {0.0, 0.5}, {1.0, 0.5}})).find("reverse") != std::string_view::npos);
    CHECK(rejection(metresPoints({{0.0, 0.0}, {0.0, 1.0}})).find("horizontal") != std::string_view::npos);
}

TEST_CASE("Surface height interpolates and takes the top of vertical faces")
{
    auto profile = PolylineProfile::create(metresPoints({{-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.1}, {0.2, 0.1}, {0.2, 0.0}, {1.0, 0.0}}));
    REQUIRE(profile);
    CHECK(profile->surfaceHeight(metres(-0.5)).si() == 0.0);
    CHECK(profile->surfaceHeight(metres(0.0)).si() == doctest::Approx(0.1));
    CHECK(profile->surfaceHeight(metres(0.1)).si() == doctest::Approx(0.1));
    CHECK(profile->surfaceHeight(metres(0.2)).si() == doctest::Approx(0.1));
    CHECK(profile->surfaceHeight(metres(0.6)).si() == 0.0);
    CHECK(profile->surfaceHeight(metres(-1.0)).si() == 0.0);
    CHECK(profile->surfaceHeight(metres(1.0)).si() == 0.0);
    CHECK(std::isinf(profile->surfaceHeight(metres(1.5)).si()));

    auto ramp = PolylineProfile::create(metresPoints({{0.0, 0.0}, {2.0, 1.0}}));
    REQUIRE(ramp);
    CHECK(ramp->surfaceHeight(metres(1.0)).si() == doctest::Approx(0.5));
}
