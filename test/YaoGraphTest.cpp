#include <gtest/gtest.h>

#include "Utils/ASSERT.hpp"

//#define LOG_DEBUG
//#undef LOG_DEBUG
#include "Utils/Logging.hpp"

#include "Algorithms/GridYao.hpp"
#include "Algorithms/NaiveYao.hpp"
#include "Algorithms/SweepLine.hpp"
#include "Generators/Generators.hpp"

#include "Kernels/CGALKernel.hpp"
#include "Kernels/InexactKernel.hpp"

constexpr tDim K = 6;
constexpr tBox YaoTestBOUNDS{{0, 0},
                             {1, 1}};
constexpr tIndex YaoTestN = 1e5;
constexpr tIndex YaoTestGenSeed = SEED;

using Dist = Uniform;

template<typename IsAlgorithm, typename ExpAlgorithm, typename Distribution, typename... Args>
bool performTest(Args... args) {

    Distribution gen(YaoTestGenSeed);
    auto [points, bounds] = gen.generate(YaoTestN, YaoTestBOUNDS);

    IsAlgorithm isAlg;
    auto is = isAlg(K, points, bounds, args...);

#ifndef VTUNE

    ExpAlgorithm expAlg;
    auto exp = expAlg(K, points, bounds);

    auto [valid, invalidVertices] = checkGraph(is, exp);

#ifdef WITH_CAIRO
    if (!valid) {
        Painter basePainter(YaoTestBOUNDS, 1000);
        basePainter.draw(points, true);
        basePainter.draw(exp, points);

        for (auto idx : invalidVertices) {

            Painter painter(basePainter);

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
            painter.setColor(0, 0, 1);
            painter.drawCones(points[idx], K);

            painter.save("invalidVertices_" + IsAlgorithm::name() + "_" + std::to_string(idx));
        }
    }
#endif// ifdef WITH_CAIRO
#else // ifndef VTUNE
    bool valid = true;
#endif// ifndef VTUNE

    return valid;
}

TEST(YaoGraphTest, SweeplineInexact) {
    ASSERT_TRUE((performTest<SweepLine<InexactKernel>, NaiveYao<CGALKernel<ExactPredicatesExactConstructions>>, Dist>()));
}

TEST(YaoGraphTest, SweeplineCGALInexact) {
    ASSERT_TRUE((performTest<SweepLine<CGALKernel<ExactPredicatesInexactConstructions>>, NaiveYao<CGALKernel<ExactPredicatesExactConstructions>>, Dist>()));
}

TEST(YaoGraphTest, SweeplineCGALExact) {
    ASSERT_TRUE((performTest<SweepLine<CGALKernel<ExactPredicatesExactConstructions>>, NaiveYao<CGALKernel<ExactPredicatesExactConstructions>>, Dist>()));
}

TEST(YaoGraphTest, NaiveInexact) {
    ASSERT_TRUE((performTest<NaiveYao<InexactKernel>, NaiveYao<CGALKernel<ExactPredicatesExactConstructions>>, Dist>()));
}

TEST(YaoGraphTest, NaiveCGALInexact) {
    ASSERT_TRUE((performTest<NaiveYao<CGALKernel<ExactPredicatesInexactConstructions>>, NaiveYao<CGALKernel<ExactPredicatesExactConstructions>>, Dist>()));
}

TEST(YaoGraphTest, GridYaoInexact) {
    ASSERT_TRUE((performTest<GridYao<InexactKernel>, NaiveYao<CGALKernel<ExactPredicatesExactConstructions>>, Dist>(100)));
}

TEST(YaoGraphTest, GridYaoCGALInexact) {
    ASSERT_TRUE((performTest<GridYao<CGALKernel<ExactPredicatesInexactConstructions>>, NaiveYao<CGALKernel<ExactPredicatesExactConstructions>>, Dist>(100)));
}

TEST(YaoGraphTest, GridYaoCGALExact) {
    ASSERT_TRUE((performTest<GridYao<CGALKernel<ExactPredicatesExactConstructions>>, NaiveYao<CGALKernel<ExactPredicatesExactConstructions>>, Dist>(100)));
}
