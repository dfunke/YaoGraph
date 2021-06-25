//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include <CGAL/Construct_theta_graph_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "Types.hpp"

template<tDim Cones>
class CGAL_Theta2D {

public:
    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point = K::Point_2;
    using Direction = K::Direction_2;
    using Graph = boost::adjacency_list<boost::listS,
                                        boost::vecS,
                                        boost::undirectedS,
                                        Point>;
    using Theta = CGAL::Construct_theta_graph_2<K, Graph>;

public:
    auto operator()(const tPoints &points) {

        Theta theta(Cones);
        Graph g;

        auto transform = [](const auto &p) -> auto {
            return Point(p[0], p[1]);
        };

        theta(boost::make_transform_iterator(points.begin(), transform),
              boost::make_transform_iterator(points.end(), transform),
              g);

        return g;
    }
};