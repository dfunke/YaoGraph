//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include <CGAL/Construct_yao_graph_2.h>

#include "Kernels/CGALKernel.hpp"
#include "Types.hpp"

template<tDim Cones, typename Kernel>
class CGAL_Yao2D {

public:
    using K = Kernel;
    using Point = typename K::Point_2;
    using Direction = typename K::Direction_2;
    using Graph = boost::adjacency_list<boost::listS,
                                        boost::vecS,
                                        boost::undirectedS,
                                        Point>;
    using Yao = CGAL::Construct_yao_graph_2<K, Graph>;

    static std::string name() {
        return "CGALYao_CGAL" + KernelName<K>::name();
    }

public:
    auto operator()(const tPoints &points, [[maybe_unused]] const tBox &bounds) {

        Yao yao(Cones);
        Graph g;

        auto transform = [](const auto &p) -> auto {
            return Point(p[0], p[1]);
        };

        yao(boost::make_transform_iterator(points.begin(), transform),
            boost::make_transform_iterator(points.end(), transform),
            g);

        return g;
    }
};