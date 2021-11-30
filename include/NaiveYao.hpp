//
// Created by Daniel Funke on 18.03.21.
//

#pragma once

#include "Predicates.hpp"
#include "Types.hpp"
#include "utils/InexactKernel.hpp"

template<tDim C, typename Kernel>
class NaiveYao {

public:
    using tVertex = tYaoVertex<C>;
    using tGraph = tYaoGraph<tVertex>;

    using tPoint = typename Kernel::Point;

    tGraph operator()(const tPoints &points) const {
        tGraph g(points.size());

        for (tIndex i = 0; i < points.size(); ++i) {
            for (tIndex j = 0; j < points.size(); ++j) {

                if (i == j) continue;

                //auto d = InexactKernel::distance2(points[i], points[j]);
                auto d = Kernel::distance2(Kernel::mkPoint(points[i]), Kernel::mkPoint(points[j]));
                tDim sec = std::floor(atan2P(points[i], points[j]) / (2 * M_PI / C));

                if (d < g[i].distance[sec]) {
                    g[i].neighbor[sec] = j;
                    g[i].distance[sec] = d;
                }
            }
        }

        return g;
    }
};