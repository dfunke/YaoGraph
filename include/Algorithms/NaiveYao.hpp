//
// Created by Daniel Funke on 18.03.21.
//

#pragma once

#include "Kernels/CGALKernel.hpp"
#include "Kernels/InexactKernel.hpp"
#include "Predicates.hpp"
#include "Types.hpp"

template<tDim C, typename Kernel>
class NaiveYao {

public:
    using tVertex = tYaoVertex<C, typename Kernel::Float>;
    using tGraph = tYaoGraph<tVertex>;

    using tPoint = typename Kernel::Point;

    static std::string name() {
        return "NaiveYao_" + Kernel::name();
    }

    auto operator()(const tPoints &iPoints, [[maybe_unused]] const tBox &bounds) const {
        tGraph g(iPoints.size());

        auto rays = Kernel::computeCones(C);
        auto kPoints = Kernel::mkPoints(iPoints);

        for (tIndex i = 0; i < kPoints.size(); ++i) {

            auto cLines = Kernel::computePointCones(kPoints[i], rays);

            for (tIndex j = 0; j < kPoints.size(); ++j) {

                if (i == j) continue;

                auto sec = Kernel::getCone(kPoints[i], kPoints[j], cLines);

                if (g[i].neighbor[sec] == INF_IDX || Kernel::compareDistance(kPoints[i], kPoints[j], kPoints[g[i].neighbor[sec]], g[i].distance[sec])) {
                    g[i].neighbor[sec] = j;
                    g[i].distance[sec] = Kernel::distance2(kPoints[i], kPoints[j]);
                }
            }
        }

        return g;
    }
};