//
// Created by Daniel Funke on 18.03.21.
//

#pragma once

#include "Types.hpp"
#include "Predicates.hpp"

template <tDim K>
class NaiveYao {

public:
    using tVertex = tYaoVertex<K>;
    using tGraph = tYaoGraph<tVertex>;

    tGraph operator()(const tPoints& points) const {
        tGraph g(points.size());

        for(tIndex i = 0; i < points.size(); ++i){
            for(tIndex j = 0; j < points.size(); ++j){

                if(i == j) continue;

                auto d = distance2(points[i], points[j]);
                tDim sec = std::floor(atan2P(points[i], points[j]) / (2 * M_PI / K));

                if(d < g[i].distance[sec]){
                    g[i].neighbor[sec] = j;
                    g[i].distance[sec] = d;
                }

            }
        }

        return g;
    }

};