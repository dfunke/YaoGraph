//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include <random>

#include "Predicates.hpp"
#include "Types.hpp"

template<typename Dist>
class Generator {

public:
    Generator(const tIndex &seed) : gen(seed) {}

protected:
    auto rand() {
        return static_cast<Dist &>(*this).dist(gen);
    }

private:
    std::mt19937 gen;
};

class Uniform : public Generator<Uniform> {

    friend class Generator<Uniform>;

public:
    Uniform(const tIndex &seed) : Generator(seed),
                                  dist(0, std::nextafter(1, std::numeric_limits<tIFloat>::max())) {}

    static std::string name() {
        return "uni";
    }

    tPoints generate(const tIndex n, const tBox &bounds) {
        tPoints points;
        points.reserve(n);

        for (tIndex i = 0; i < n; ++i) {

            tIFloatVector p;
            for (tDim d = 0; d < D; ++d) {
                p[d] = bounds.low[d] + (bounds.high[d] - bounds.low[d]) * rand();
            }

            points.emplace_back(p);
        }

        return points;
    }

private:
    std::uniform_real_distribution<tIFloat> dist;
};

class Gaussian : public Generator<Gaussian> {

    friend class Generator<Gaussian>;

public:
    Gaussian(const tIndex &seed) : Generator(seed),
                                   dist(0, .25) {}

    static std::string name() {
        return "gaussian";
    }

    tPoints generate(const tIndex n, const tBox &bounds) {
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

            points.emplace_back(p);
        }

        return points;
    }

private:
    std::normal_distribution<tIFloat> dist;
};

class Grid : public Generator<Grid> {

    friend class Generator<Grid>;

public:
    Grid(const tIndex &seed) : Generator(seed) {}

    static std::string name() {
        return "grid";
    }

    tPoints generate(const tIndex n, const tBox &bounds) {
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
                p[X] = bounds.low[X] + pointSpacing[X] * i;
                p[Y] = bounds.low[Y] + pointSpacing[Y] * j;

                points.emplace_back(p);
            }
        }

        return points;
    }
};
