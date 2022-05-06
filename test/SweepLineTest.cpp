#include <gtest/gtest.h>

#include "Generators.hpp"
#include "NaiveYao.hpp"
#include "SweepLine.hpp"

#include "utils/CGALKernel.hpp"
#include "utils/InexactKernel.hpp"

TEST(SweepLineTest, Ordering) {
    constexpr tBox BOUNDS{{0, 0},
                          {1, 1}};
    constexpr tIndex nPoints = 1000;
    constexpr tDim K = 6;

    Uniform uni;
    auto points = uni.generate(nPoints, BOUNDS);

    std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1);

    SweepLine<K, CGALKernel<ExactPredicatesInexactConstructions>> sl;
    auto is = sl(points, BOUNDS);

#ifndef VTUNE

    NaiveYao<K, CGALKernel<ExactPredicatesInexactConstructions>> nav;
    auto exp = nav(points);

    auto [valid, invalidVertices] = checkGraph(is, exp);

#ifdef WITH_CAIRO
    if (!valid) {
        Painter painter(BOUNDS, 1000);
        painter.draw(points, false);
        painter.draw(exp, points);
        for (auto idx : invalidVertices) {

            for (tDim k = 0; k < K; ++k) {
                if (exp[idx].neighbor[k] != is[idx].neighbor[k]) {

                    painter.setColor(0, 1, 0);
                    painter.draw(points[exp[idx].neighbor[k]], exp[idx].neighbor[k], true);

                    painter.setColor(1, 0, 0);
                    painter.draw(points[is[idx].neighbor[k]], is[idx].neighbor[k], true);
                }
            }

            painter.setColor(1, 0, 0);
            painter.draw(idx, is[idx], points);
            painter.drawCones(points[idx], K);
        }
        painter.save("invalidVertices");
    }
#endif// ifdef WITH_CAIRO
#else // ifndef VTUNE
    bool valid = true;
#endif// ifndef VTUNE

    EXPECT_TRUE(valid);
}