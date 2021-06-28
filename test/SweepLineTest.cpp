#include <gtest/gtest.h>

#include "Generators.hpp"
#include "SweepLine.hpp"
#include "NaiveYao.hpp"

TEST(SweepLineTest, Ordering) { // 12/2/2020 -> 737761
    constexpr tBox BOUNDS{{0, 0},
                          {1, 1}};
    constexpr tIndex nPoints = 10;
    constexpr tDim K = 6;


    Uniform uni;
    auto points = uni.generate(nPoints, BOUNDS);

    SweepLine<K> sl;
    auto is = sl(points, BOUNDS);

    NaiveYao<K> nav;
    auto exp = nav(points);

    ASSERT_TRUE(checkGraph(is, exp));

}