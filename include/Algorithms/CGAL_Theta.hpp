//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include <CGAL/Construct_theta_graph_2.h>

#include "Kernels/CGALKernel.hpp"
#include "Types.hpp"

template<typename Kernel>
class CGAL_Theta2D {

public:
    using K = Kernel;
    using Point = typename K::Point_2;
    using Direction = typename K::Direction_2;
    using Graph = boost::adjacency_list<boost::listS,
                                        boost::vecS,
                                        boost::undirectedS,
                                        Point>;
    using Theta = CGAL::Construct_theta_graph_2<K, Graph>;

    static std::string name() {
        return "CGALTheta_CGAL" + KernelName<K>::name();
    }

public:
    auto operator()(const tDim &K, const tPoints &points, [[maybe_unused]] const tBox &bounds) {

        Theta theta(K);
        Graph g;

        auto transform = [](const auto &p) -> auto{
            return Point(p[0], p[1]);
        };

        theta(boost::make_transform_iterator(points.begin(), transform),
              boost::make_transform_iterator(points.end(), transform),
              g);

        return g;
    }
};