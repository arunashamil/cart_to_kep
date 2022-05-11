
#include "gtest/gtest.h"
#include "../cart_to_kep.hpp"

TEST(TestCCK, Test)
{
    std::array<double,6> cartesian_data{-104162698.3, 269165908, 25673581.8, -491.300, 515.530, -131.348};
    auto res = cart_to_kep(cartesian_data);
    for (int i = 0; i < 6; ++i) {
        std::cout << res[i] << " ";
    }
    std::array<double, 6> data = {178995.7, 0.93314, 35.1 * M_PI / 180., M_PI - 61.6 * M_PI / 180., 0.7  * M_PI / 180, 170.4 * M_PI / 180};
    ASSERT_EQ(res, data);
}