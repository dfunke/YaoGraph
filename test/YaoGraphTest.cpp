#include <gtest/gtest.h>

//#define LOG_DEBUG
//#undef LOG_DEBUG

#include "utils/Logging.hpp"

#include "Generators.hpp"
#include "NaiveYao.hpp"
#include "SweepLine.hpp"

#include "utils/CGALKernel.hpp"
#include "utils/InexactKernel.hpp"

constexpr tDim K = 6;
constexpr tBox BOUNDS{{0, 0},
                      {1, 1}};
constexpr tIndex nPoints = 2000;

template<typename IsAlgorithm, typename ExpAlgorithm, typename Distribution>
bool performTest() {

    Distribution gen(SEED);
    auto points = gen.generate(nPoints, BOUNDS);

    IsAlgorithm isAlg;
    auto is = isAlg(points, BOUNDS);

#ifndef VTUNE

    ExpAlgorithm expAlg;
    auto exp = expAlg(points, BOUNDS);

    auto [valid, invalidVertices] = checkGraph(is, exp);

#ifdef WITH_CAIRO
    if (!valid) {
        Painter painter(BOUNDS, 1000);
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

            break;
        }
        painter.save("invalidVertices_" + IsAlgorithm::name() + ".png");
    }
#endif// ifdef WITH_CAIRO
#else // ifndef VTUNE
    bool valid = true;
#endif// ifndef VTUNE

    return valid;
}

TEST(YaoGraphTest, SweeplineInexact) {
    ASSERT_TRUE((performTest<SweepLine<K, InexactKernel>, NaiveYao<K, CGALKernel<ExactPredicatesExactConstructions>>, Uniform>()));
}

TEST(YaoGraphTest, SweeplineCGALInexact) {
    ASSERT_TRUE((performTest<SweepLine<K, CGALKernel<ExactPredicatesInexactConstructions>>, NaiveYao<K, CGALKernel<ExactPredicatesExactConstructions>>, Uniform>()));
}

TEST(YaoGraphTest, SweeplineCGALExact) {
    ASSERT_TRUE((performTest<SweepLine<K, CGALKernel<ExactPredicatesExactConstructions>>, NaiveYao<K, CGALKernel<ExactPredicatesExactConstructions>>, Uniform>()));
}

TEST(YaoGraphTest, NaiveInexact) {
    ASSERT_TRUE((performTest<NaiveYao<K, InexactKernel>, NaiveYao<K, CGALKernel<ExactPredicatesExactConstructions>>, Uniform>()));
}

TEST(YaoGraphTest, NaiveCGALInexact) {
    ASSERT_TRUE((performTest<NaiveYao<K, CGALKernel<ExactPredicatesInexactConstructions>>, NaiveYao<K, CGALKernel<ExactPredicatesExactConstructions>>, Uniform>()));
}
