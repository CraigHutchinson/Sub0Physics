#include <doctest.h>

#include "sub0physics/verification/channel_table.hpp"

using namespace sub0physics;

TEST_CASE("Channel table stores contiguous columns up to capacity")
{
    ChannelTable<3> table(2);
    CHECK(table.empty());
    CHECK(table.append({1.0, 2.0, 3.0}));
    CHECK(table.append({4.0, 5.0, 6.0}));
    CHECK_FALSE(table.append({7.0, 8.0, 9.0}));
    CHECK(table.size() == 2U);

    const auto second = table.column(1);
    REQUIRE(second.size() == 2U);
    CHECK(second[0] == 2.0);
    CHECK(second[1] == 5.0);
    CHECK(table.row(1) == ChannelTable<3>::Row{4.0, 5.0, 6.0});
}
