#include <gtest/gtest.h>

//#define LOG_DEBUG
//#undef LOG_DEBUG

#include "utils/Logging.hpp"

#include "Generators.hpp"
#include "NaiveYao.hpp"
#include "SweepLine.hpp"

#include "utils/CGALKernel.hpp"
#include "utils/InexactKernel.hpp"

TEST(SweepLineTest, Test) {
    constexpr tBox BOUNDS{{0, 0},
                          {1, 1}};
    constexpr tIndex nPoints = 20000;
    constexpr tDim K = 6;

    Uniform uni(SEED);
    auto points = uni.generate(nPoints, BOUNDS);

    SweepLine<K, InexactKernel> sl;
    auto is = sl(points, BOUNDS);

#ifndef VTUNE

    NaiveYao<K, CGALKernel<ExactPredicatesInexactConstructions>> nav;
    auto exp = nav(points, BOUNDS);

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

TEST(NaiveYaoTest, Test) {
    constexpr tBox BOUNDS{{0, 0},
                          {1, 1}};
    constexpr tIndex nPoints = 1000;
    constexpr tDim K = 6;

    Uniform uni(SEED);
    auto points = uni.generate(nPoints, BOUNDS);

    NaiveYao<K, InexactKernel> naive_ie;
    auto ie = naive_ie(points, BOUNDS);

    NaiveYao<K, CGALKernel<ExactPredicatesInexactConstructions>> naive_cgalIE;
    auto cgalIE = naive_cgalIE(points, BOUNDS);

    NaiveYao<K, CGALKernel<ExactPredicatesInexactConstructions>> naive_cgalE;
    auto cgalE = naive_cgalE(points, BOUNDS);

    auto [valid_ie_cgalE, invalidVertices_ie_cgalE] = checkGraph(ie, cgalE);
    auto [valid_cgalIE_cgalE, invalidVertices_cgalIE_cgalE] = checkGraph(cgalIE, cgalE);

#ifdef WITH_CAIRO
    if (!valid_ie_cgalE) {

        auto &is = ie;
        auto &exp = cgalE;
        auto &invalidVertices = invalidVertices_ie_cgalE;

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
        painter.save("invalidVertices_ie_cgalE");
    }

    if (!valid_cgalIE_cgalE) {

        auto &is = cgalIE;
        auto &exp = cgalE;
        auto &invalidVertices = invalidVertices_cgalIE_cgalE;

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
        painter.save("invalidVertices_cgalIE_cgalE");
    }
#endif// ifdef WITH_CAIRO

    EXPECT_TRUE(valid_ie_cgalE && valid_cgalIE_cgalE);
}