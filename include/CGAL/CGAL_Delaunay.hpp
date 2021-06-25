//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "Types.hpp"

class CGAL_Delaunay2D {

public:
    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using DT = CGAL::Delaunay_triangulation_2<K>;
    using Point = DT::Point;

public:
    auto operator()(const tPoints &points) {

        DT dt;

        auto transform = [](const auto &p) -> auto {
            return Point(p[0], p[1]);
        };

        dt.insert(boost::make_transform_iterator(points.begin(), transform),
                  boost::make_transform_iterator(points.end(), transform));

        return dt;
    }
};