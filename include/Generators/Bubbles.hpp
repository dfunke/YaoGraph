//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include "GeneratorBase.hpp"

class Bubbles : public Generator<Bubbles> {

    friend class Generator<Bubbles>;

public:
    Bubbles(const tIndex &seed) : Generator(seed),
                                  dist(0, std::nextafter(1, std::numeric_limits<tIFloat>::max())) {}

    std::string name() const override {
        return "bubbles";
    }

    std::tuple<tPoints, tBox> _generate(const tIndex n, const tBox &bounds) override {
        tPoints points;
        points.reserve(n);

        auto gMidpoint = bounds.low;

        for (uint d = 0; d < D; ++d) {
            gMidpoint[d] += (bounds.high[d] - bounds.low[d]) / 2;
        }

        std::vector<tBox> bubbles;
        bubbles.resize(pow(2, D));

        for (uint i = 0; i < pow(2, D); ++i) {
            for (uint d = 0; d < D; ++d) {
                bubbles[i].low[d] =
                        i & (1 << d) ? gMidpoint[d] : bounds.low[d];
                bubbles[i].high[d] =
                        i & (1 << d) ? bounds.high[d] : gMidpoint[d];
            }
        }

        tIndex ppB = n / bubbles.size();// points per bubble
        for (uint b = 0; b < bubbles.size(); ++b) {
            const auto &bubble = bubbles[b];

            auto midpoint = bubble.low;
            auto span = bubble.high;

            for (uint d = 0; d < D; ++d) {
                midpoint[d] += (bubble.high[d] - bubble.low[d]) / 2;
                span[d] -= bubble.low[d];
            }

            for (tIndex i = b * ppB + 1; i <= (b + 1) * ppB; ++i) {
                tIFloatVector p;

                do {
                    p = midpoint;

                    for (uint d = 0; d < D; ++d) {
                        p[d] += span[d] * rand();
                    }
                } while (!bounds.contains(p));


                points.push_back(p);
            }
        }

        return std::make_tuple(points, bounds);
    }

private:
    std::uniform_real_distribution<tIFloat> dist;
};