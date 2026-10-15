#include <cmath>
#include <doctest.h>
#include <limits>

#include "sub0physics/verification/error_statistics.hpp"

using namespace sub0physics;

TEST_CASE("Error statistics report peak and RMS against a limit")
{
    ErrorStatistics stats(0.5);
    CHECK_FALSE(stats.passes());
    stats.add(0.3);
    stats.add(-0.4);
    CHECK(stats.maxAbsolute() == doctest::Approx(0.4));
    CHECK(stats.rms() == doctest::Approx(std::sqrt((0.09 + 0.16) / 2.0)));
    CHECK(stats.count() == 2U);
    CHECK(stats.passes());
    stats.add(0.6);
    CHECK_FALSE(stats.passes());
}

TEST_CASE("A non-finite difference can never pass")
{
    ErrorStatistics stats(1.0);
    stats.add(std::numeric_limits<double>::quiet_NaN());
    stats.add(0.0);
    CHECK_FALSE(stats.isFinite());
    CHECK_FALSE(stats.passes());
}
