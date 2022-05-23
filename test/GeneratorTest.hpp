//
// Created by Daniel Funke on 20.05.22.
//

#include <gtest/gtest.h>

#include "Utils/ASSERT.hpp"
#include "Generators.hpp"

constexpr tBox GeneratorTestBOUNDS{{0, 0},
                      {1, 1}};
constexpr tIndex GeneratorTestN = 1000;
constexpr tIndex GeneratorTestGenSeed = SEED;

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

template<typename Distribution, bool EQ = true>
void performTest() {
    Road gen(GeneratorTestGenSeed);
    auto points = gen.generate(GeneratorTestN, GeneratorTestBOUNDS);

    if constexpr(EQ)
        EXPECT_EQ(points.size(), GeneratorTestN);
    else
        EXPECT_LE(points.size(), GeneratorTestN);
    EXPECT_TRUE(checkWithinBounds(points, GeneratorTestBOUNDS));
}

TEST(GeneratorTest, Uniform) {
    performTest<Uniform>();
}

TEST(GeneratorTest, Gaussian) {
    performTest<Gaussian>();
}

TEST(GeneratorTest, Grid) {
    performTest<Grid, false>();
}

TEST(GeneratorTest, Road) {
    performTest<Road>();
}