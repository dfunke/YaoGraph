#include <gtest/gtest.h>

#include "Utils/ASSERT.hpp"

//#define LOG_DEBUG
//#undef LOG_DEBUG
#include "Utils/Logging.hpp"

#include "Algorithms/GridYao.hpp"
#include "Algorithms/NaiveYao.hpp"
#include "Algorithms/SweepLine.hpp"
#include "Generators.hpp"

#include "Kernels/CGALKernel.hpp"
#include "Kernels/InexactKernel.hpp"

constexpr tDim K = 6;
constexpr tBox YaoTestBOUNDS{{0, 0},
                      {1, 1}};
constexpr tIndex YaoTestN = 1000;
constexpr tIndex YaoTestGenSeed = SEED;

using Dist = Grid;

template<typename IsAlgorithm, typename ExpAlgorithm, typename Distribution>
bool performTest() {

    Distribution gen(YaoTestGenSeed);
    auto points = gen.generate(YaoTestN, YaoTestBOUNDS);

    IsAlgorithm isAlg;
    auto is = isAlg(points, YaoTestBOUNDS);

#ifndef VTUNE

    ExpAlgorithm expAlg;
    auto exp = expAlg(points, YaoTestBOUNDS);

    auto [valid, invalidVertices] = checkGraph(is, exp);

#ifdef WITH_CAIRO
    if (!valid) {
        Painter painter(YaoTestBOUNDS, 1000);
        painter.draw(points, false);
        painter.draw(exp, points);
        for (auto idx : invalidVertices) {

            for (tDim k = 0; k < K; ++k) {
                if (exp[idx].neighbor[k] != is[idx].neighbor[k]) {

                    if (exp[idx].neighbor[k] != INF_IDX) {
                        painter.setColor(0, 1, 0);
                        painter.draw(points[exp[idx].neighbor[k]], exp[idx].neighbor[k], true);
                    }

                    if (is[idx].neighbor[k] != INF_IDX) {
                        painter.setColor(1, 0, 0);
                        painter.draw(points[is[idx].neighbor[k]], is[idx].neighbor[k], true);
                    }
                }
            }

            painter.setColor(1, 0, 0);
            painter.draw(idx, is[idx], points);
            painter.drawCones(points[idx], K);
        }
        painter.save("invalidVertices_" + IsAlgorithm::name());
    }
#endif// ifdef WITH_CAIRO
#else // ifndef VTUNE
    bool valid = true;
#endif// ifndef VTUNE

    return valid;
}

TEST(YaoGraphTest, SweeplineInexact) {
    ASSERT_TRUE((performTest<SweepLine<K, InexactKernel>, NaiveYao<K, CGALKernel<ExactPredicatesExactConstructions>>, Dist>()));
}

TEST(YaoGraphTest, SweeplineCGALInexact) {
    ASSERT_TRUE((performTest<SweepLine<K, CGALKernel<ExactPredicatesInexactConstructions>>, NaiveYao<K, CGALKernel<ExactPredicatesExactConstructions>>, Dist>()));
}

TEST(YaoGraphTest, SweeplineCGALExact) {
    ASSERT_TRUE((performTest<SweepLine<K, CGALKernel<ExactPredicatesExactConstructions>>, NaiveYao<K, CGALKernel<ExactPredicatesExactConstructions>>, Dist>()));
}

TEST(YaoGraphTest, NaiveInexact) {
    ASSERT_TRUE((performTest<NaiveYao<K, InexactKernel>, NaiveYao<K, CGALKernel<ExactPredicatesExactConstructions>>, Dist>()));
}

TEST(YaoGraphTest, NaiveCGALInexact) {
    ASSERT_TRUE((performTest<NaiveYao<K, CGALKernel<ExactPredicatesInexactConstructions>>, NaiveYao<K, CGALKernel<ExactPredicatesExactConstructions>>, Dist>()));
}

TEST(YaoGraphTest, GridYaoInexact) {
    ASSERT_TRUE((performTest<GridYao<K, InexactKernel, 100>, NaiveYao<K, CGALKernel<ExactPredicatesExactConstructions>>, Dist>()));
}

TEST(YaoGraphTest, GridYaoCGALInexact) {
    ASSERT_TRUE((performTest<GridYao<K, CGALKernel<ExactPredicatesInexactConstructions>, 100>, NaiveYao<K, CGALKernel<ExactPredicatesExactConstructions>>, Dist>()));
}

TEST(YaoGraphTest, GridYaoCGALExact) {
    ASSERT_TRUE((performTest<GridYao<K, CGALKernel<ExactPredicatesExactConstructions>, 100>, NaiveYao<K, CGALKernel<ExactPredicatesExactConstructions>>, Dist>()));
}
