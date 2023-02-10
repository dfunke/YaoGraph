//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include "GeneratorBase.hpp"

class Grid : public Generator<Grid> {

    friend class Generator<Grid>;

public:
    Grid(const tIndex &seed) : Generator(seed) {}

    std::string name() const override {
        return "grid";
    }

    std::tuple<tPoints, tBox> _generate(const tIndex n, const tBox &bounds) override {
        tPoints points;
        points.reserve(n);

        tIndex pointsPerDim = std::floor(std::sqrt(n));

        tIFloatVector pointSpacing;
        for (tDim d = 0; d < bounds.high.size(); ++d) {
            pointSpacing[d] = (bounds.high[d] - bounds.low[d]) / pointsPerDim;
        }

        for (tIndex i = 0; i < pointsPerDim; ++i) {
            for (tIndex j = 0; j < pointsPerDim; ++j) {
                tIFloatVector p;
                p[X] = bounds.low[X] + (pointSpacing[X] / 2) + pointSpacing[X] * i;
                p[Y] = bounds.low[Y] + (pointSpacing[Y] / 2) + pointSpacing[Y] * j;

                ASSERT(bounds.contains(p));
                points.emplace_back(p);
            }
        }

        return std::make_tuple(points, bounds);
    }
};
