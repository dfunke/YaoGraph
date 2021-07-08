#include <gtest/gtest.h>

#include "Generators.hpp"
#include "NaiveYao.hpp"
#include "SweepLine.hpp"

TEST(SweepLineTest, Ordering) {
    constexpr tBox BOUNDS{{0, 0},
                          {1, 1}};
    constexpr tIndex nPoints = 1e2;
    constexpr tDim K = 6;


    Uniform uni;
    auto points = uni.generate(nPoints, BOUNDS);

    SweepLine<K> sl;
    auto is = sl(points, BOUNDS);

    NaiveYao<K> nav;
    auto exp = nav(points);

    auto [valid, invalidVertices] = checkGraph(is, exp);

#ifdef WITH_CAIRO
    if (!valid) {
        Painter painter(BOUNDS, 1000);
        painter.draw(points);
        painter.draw(exp, points);
        painter.setColor(1, 0, 0);
        for (auto idx : invalidVertices) {
            painter.draw(idx, is[idx], points);
            painter.drawCones(points[idx], K);
        }
        painter.save("invalidVertices");
    }
#endif

    EXPECT_TRUE(valid);
}