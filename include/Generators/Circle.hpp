//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include "GeneratorBase.hpp"

class Circle : public Generator<Circle> {

    friend class Generator<Circle>;

public:
    Circle(const tIndex &seed) : Generator(seed),
                                 dist(0, std::nextafter(1, std::numeric_limits<tIFloat>::max())) {}

    std::string name() const override {
        return "circle";
    }

    std::tuple<tPoints, tBox> _generate(const tIndex n, const tBox &bounds) override {
        tPoints points;
        points.reserve(n);

        tIFloatVector radius;
        auto center = bounds.low;

        for (uint d = 0; d < D; ++d) {
            center[d] += (bounds.high[d] - bounds.low[d]) / 2;
            radius[d] = .45 * (bounds.high[d] - bounds.low[d]);//give it some space at the edge
        }

        // we generate normal distributed points and map them to points on a unit circle/sphere
        // these points are uniform distributed on the circle/sphere
        // we then map the circle/sphere points to the ellipse/ellipsoid
        // the resulting points are not uniform distributed, but close enough

        for (tIndex i = 1; i <= n; ++i) {
            tIFloatVector p;

            //normal distributed point around origin
            tIFloat r = 0;//radius
            for (uint d = 0; d < D; ++d) {
                p[d] = rand();
                r += p[d] * p[d];
            }

            r = std::sqrt(r);

            // now map them to the ellipse/ellipsoid
            for (uint d = 0; d < D; ++d) {
                //            origin      distortion  point on circle/sphere
                p[d] = center[d] + radius[d] * (p[d] / r);
            }


            points.push_back(p);
        }

        return std::make_tuple(points, bounds);
    }

private:
    std::normal_distribution<tIFloat> dist;
};