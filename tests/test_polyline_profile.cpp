#include <cmath>
#include <doctest.h>
#include <limits>
#include <vector>

#include "sub0physics/geometry/polyline_profile.hpp"

using namespace sub0physics;

namespace
{

std::string_view rejection(std::vector<Point2> points)
{
    const auto profile = PolylineProfile::create(std::move(points));
    REQUIRE_FALSE(profile);
    CHECK(profile.error().kind == FaultKind::InvalidInput);
    return profile.error().detail;
}

} // namespace

TEST_CASE("Profile validation rejects malformed polylines")
{
    CHECK(rejection({{0.0, 0.0}}).find("two points") != std::string_view::npos);
    CHECK(rejection({{0.0, 0.0}, {std::numeric_limits<double>::quiet_NaN(), 0.0}}).find("finite") != std::string_view::npos);
    CHECK(rejection({{0.0, 0.0}, {1e101, 0.0}}).find("limit") != std::string_view::npos);
    CHECK(rejection({{0.0, 0.0}, {1.0, 0.0}, {0.5, 0.0}}).find("overhang") != std::string_view::npos);
    CHECK(rejection({{0.0, 0.0}, {1e-7, 0.0}}).find("micrometre") != std::string_view::npos);
    CHECK(rejection({{0.0, 0.0}, {0.0, 1.0}, {0.0, 0.5}, {1.0, 0.5}}).find("reverse") != std::string_view::npos);
    CHECK(rejection({{0.0, 0.0}, {0.0, 1.0}}).find("horizontal") != std::string_view::npos);
}

TEST_CASE("Surface height interpolates and takes the top of vertical faces")
{
    auto profile = PolylineProfile::create({{-1.0, 0.0}, {0.0, 0.0}, {0.0, 0.1}, {0.2, 0.1}, {0.2, 0.0}, {1.0, 0.0}});
    REQUIRE(profile);
    CHECK(profile->surfaceHeight(-0.5) == 0.0);
    CHECK(profile->surfaceHeight(0.0) == doctest::Approx(0.1));
    CHECK(profile->surfaceHeight(0.1) == doctest::Approx(0.1));
    CHECK(profile->surfaceHeight(0.2) == doctest::Approx(0.1));
    CHECK(profile->surfaceHeight(0.6) == 0.0);
    CHECK(profile->surfaceHeight(-1.0) == 0.0);
    CHECK(profile->surfaceHeight(1.0) == 0.0);
    CHECK(std::isinf(profile->surfaceHeight(1.5)));

    auto ramp = PolylineProfile::create({{0.0, 0.0}, {2.0, 1.0}});
    REQUIRE(ramp);
    CHECK(ramp->surfaceHeight(1.0) == doctest::Approx(0.5));
}
