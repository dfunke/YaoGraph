//
// Created by Daniel Funke on 20.05.22.
//

#include <gtest/gtest.h>

#include "Utils/ASSERT.hpp"
#include "Generators.hpp"

constexpr tBox BOUNDS{{0, 0},
                      {1, 1}};
constexpr tIndex nPoints = 1000;
constexpr tIndex GenSeed = SEED;

bool checkWithinBounds(const tPoints &points, const tBox &bounds) {
    for (auto const &p : points) {
        if (p[X] < bounds.low[X] || p[X] > bounds.high[X]) {
            return false;
        }
        if (p[Y] < bounds.low[Y] || p[Y] > bounds.high[Y]) {
            return false;
        }
    }

    return true;
}

TEST(GeneratorTest, Uniform) {

    Uniform gen(GenSeed);
    auto points = gen.generate(nPoints, BOUNDS);

    EXPECT_EQ(points.size(), nPoints);
    EXPECT_TRUE(checkWithinBounds(points, BOUNDS));
}

TEST(GeneratorTest, Gaussian) {

    Gaussian gen(GenSeed);
    auto points = gen.generate(nPoints, BOUNDS);

    EXPECT_EQ(points.size(), nPoints);
    EXPECT_TRUE(checkWithinBounds(points, BOUNDS));
}

TEST(GeneratorTest, Grid) {

    Grid gen(GenSeed);
    auto points = gen.generate(nPoints, BOUNDS);

    EXPECT_LE(points.size(), nPoints);
    EXPECT_TRUE(checkWithinBounds(points, BOUNDS));
}

TEST(GeneratorTest, Road) {

    Road gen(GenSeed);
    auto points = gen.generate(nPoints, BOUNDS);

    EXPECT_EQ(points.size(), nPoints);
    EXPECT_TRUE(checkWithinBounds(points, BOUNDS));
}