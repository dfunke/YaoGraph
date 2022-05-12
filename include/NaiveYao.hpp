//
// Created by Daniel Funke on 18.03.21.
//

#pragma once

#include "Predicates.hpp"
#include "Types.hpp"
#include "utils/CGALKernel.hpp"
#include "utils/InexactKernel.hpp"

template<tDim C, typename Kernel>
class NaiveYaoBase {

public:
    using tVertex = tYaoVertex<C>;
    using tGraph = tYaoGraph<tVertex>;

    using tPoint = typename Kernel::Point;

    static std::string name() {
        return "NaiveYao_" + Kernel::name();
    }
};

template<tDim C, typename Kernel>
class NaiveYao;

template<tDim C>
class NaiveYao<C, InexactKernel> : public NaiveYaoBase<C, InexactKernel> {

private:
    using Kernel = InexactKernel;
    using Base = NaiveYaoBase<C, Kernel>;

public:
    auto operator()(const tPoints &points, [[maybe_unused]] const tBox &bounds) const {
        typename Base::tGraph g(points.size());

        for (tIndex i = 0; i < points.size(); ++i) {
            for (tIndex j = 0; j < points.size(); ++j) {

                if (i == j) continue;

                auto d = Kernel::distance2(Kernel::mkPoint(points[i]), Kernel::mkPoint(points[j]));
                tDim sec = std::floor(atan2P(points[i], points[j]) / (2 * M_PI / C));

                if (d < g[i].distance[sec]) {
                    g[i].neighbor[sec] = j;
                    g[i].distance[sec] = Kernel::to_float(d);
                }
            }
        }

        return g;
    }
};

template<tDim C>
class NaiveYao<C, CGALKernel<ExactPredicatesInexactConstructions>> : public NaiveYaoBase<C, CGALKernel<ExactPredicatesInexactConstructions>> {

private:
    using Kernel = CGALKernel<ExactPredicatesInexactConstructions>;
    using Base = NaiveYaoBase<C, Kernel>;

public:
    auto operator()(const tPoints &points, [[maybe_unused]] const tBox &bounds) const {
        typename Base::tGraph g(points.size());

        auto rays = Kernel::computeCones(C);
        std::vector<typename Kernel::Point> kPoints;
        kPoints.reserve(points.size());

        for (tIndex i = 0; i < points.size(); ++i) {
            kPoints.emplace_back(Kernel::mkPoint(points[i]));
        }

        for (tIndex i = 0; i < points.size(); ++i) {

            std::vector<typename Kernel::Line> cLines;
            cLines.reserve(C);
            for (tDim k = 0; k < C; ++k) {
                cLines.emplace_back(kPoints[i], rays[k]);
            }

            for (tIndex j = 0; j < points.size(); ++j) {

                if (i == j) continue;

                tDim sec, secGuess;
                sec = secGuess = std::floor(atan2P(points[i], points[j]) / (2 * M_PI / C));
                [[maybe_unused]] bool found = false;
                for (; sec != (secGuess - 1) % C; ++sec) {
                    auto osL = cLines[sec].oriented_side(kPoints[j]);
                    auto osR = cLines[(sec + 1) % C].oriented_side(kPoints[j]);

                    if ((osL == CGAL::ON_POSITIVE_SIDE || osL == CGAL::ON_ORIENTED_BOUNDARY) && osR == CGAL::ON_NEGATIVE_SIDE) {
                        found = true;
                        break;
                    }
                }
                assert(found);

                if (g[i].neighbor[sec] == INF_IDX || Kernel::compareDistance(kPoints[i], kPoints[j], kPoints[g[i].neighbor[sec]])) {
                    g[i].neighbor[sec] = j;
                    g[i].distance[sec] = Kernel::to_float(Kernel::distance2(kPoints[i], kPoints[j]));//TODO: exact solution
                }
            }
        }

        return g;
    }
};

template<tDim C>
class NaiveYao<C, CGALKernel<ExactPredicatesExactConstructions>> : public NaiveYaoBase<C, CGALKernel<ExactPredicatesExactConstructions>> {

private:
    using Kernel = CGALKernel<ExactPredicatesExactConstructions>;
    using Base = NaiveYaoBase<C, Kernel>;

public:
    auto operator()(const tPoints &points, [[maybe_unused]] const tBox &bounds) const {
        typename Base::tGraph g(points.size());

        auto rays = Kernel::computeCones(C);
        std::vector<typename Kernel::Point> kPoints;
        kPoints.reserve(points.size());

        for (tIndex i = 0; i < points.size(); ++i) {
            kPoints.emplace_back(Kernel::mkPoint(points[i]));
        }

        for (tIndex i = 0; i < points.size(); ++i) {

            std::vector<typename Kernel::Line> cLines;
            cLines.reserve(C);

            std::array<Kernel::Float, C> dists;

            for (tDim k = 0; k < C; ++k) {
                cLines.emplace_back(kPoints[i], rays[k]);
                dists[k] = std::numeric_limits<tIFloat>::max();
            }

            for (tIndex j = 0; j < points.size(); ++j) {

                if (i == j) continue;

                tDim sec, secGuess;
                sec = secGuess = std::floor(atan2P(points[i], points[j]) / (2 * M_PI / C));
                [[maybe_unused]] bool found = false;
                for (; sec != (secGuess - 1) % C; ++sec) {
                    auto osL = cLines[sec].oriented_side(kPoints[j]);
                    auto osR = cLines[(sec + 1) % C].oriented_side(kPoints[j]);

                    if ((osL == CGAL::ON_POSITIVE_SIDE || osL == CGAL::ON_ORIENTED_BOUNDARY) && osR == CGAL::ON_NEGATIVE_SIDE) {
                        found = true;
                        break;
                    }
                }
                assert(found);

                if (g[i].neighbor[sec] == INF_IDX || Kernel::compareDistance(kPoints[i], kPoints[j], dists[sec])) {
                    g[i].neighbor[sec] = j;
                    dists[sec] = Kernel::distance2(kPoints[i], kPoints[j]);
                    g[i].distance[sec] = Kernel::to_float(dists[sec]);//TODO: exact solution
                }
            }
        }

        return g;
    }
};