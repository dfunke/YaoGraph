//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include "GeneratorBase.hpp"

class Uniform : public Generator<Uniform> {

    friend class Generator<Uniform>;

public:
    Uniform(const tIndex &seed) : Generator(seed),
                                  dist(0, std::nextafter(1, std::numeric_limits<tIFloat>::max())) {}

    std::string name() const override {
        return "uni";
    }

    std::tuple<tPoints, tBox> _generate(const tIndex n, const tBox &bounds) override {
        tPoints points;
        points.reserve(n);

        for (tIndex i = 0; i < n; ++i) {

            tIFloatVector p;
            for (tDim d = 0; d < D; ++d) {
                p[d] = bounds.low[d] + (bounds.high[d] - bounds.low[d]) * rand();
            }

            ASSERT(bounds.contains(p));
            points.emplace_back(p);
        }

        return std::make_tuple(points, bounds);
    }

private:
    std::uniform_real_distribution<tIFloat> dist;
};