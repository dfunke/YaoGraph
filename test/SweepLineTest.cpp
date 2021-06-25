#include <gtest/gtest.h>

#include "SweepLine.hpp"
#include "Generators.hpp"

TEST(SweepLineTest, Ordering) { // 12/2/2020 -> 737761
    constexpr tBox BOUNDS{{0, 0},
                          {1, 1}};
    constexpr tIndex nPoints = 10;
    constexpr tDim K = 6;


    Uniform uni;
    auto points = uni.generate(nPoints, BOUNDS);

    SweepLine<K> sl;
    sl(points, BOUNDS);
}