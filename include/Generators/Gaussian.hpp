//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include "GeneratorBase.hpp"

class Gaussian : public Generator<Gaussian> {

    friend class Generator<Gaussian>;

public:
    Gaussian(const tIndex &seed) : Generator(seed),
                                   dist(0, .25) {}

    std::string name() const override {
        return "gaussian";
    }

    std::tuple<tPoints, tBox> _generate(const tIndex n, const tBox &bounds) override {
        tPoints points;
        points.reserve(n);

        auto halfLength = .5 * (bounds.high - bounds.low);
        auto midpoint = bounds.low + halfLength;

        for (tIndex i = 0; i < n; ++i) {

            tIFloatVector p;
            do {
                for (tDim d = 0; d < D; ++d) {
                    p[d] = midpoint[d] + halfLength[d] * rand();
                }
            } while (!bounds.contains(p));

            ASSERT(bounds.contains(p));
            points.emplace_back(p);
        }

        return std::make_tuple(points, bounds);
    }

private:
    std::normal_distribution<tIFloat> dist;
};