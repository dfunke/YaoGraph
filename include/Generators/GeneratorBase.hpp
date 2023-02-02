//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include <algorithm>
#include <fstream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>

#include "Predicates.hpp"
#include "Types.hpp"
#include "Utils/ASSERT.hpp"

class GeneratorBase {
public:
    std::tuple<tPoints, tBox> generate(const tIndex n, const tBox &bounds) {
        auto [p, b] = _generate(n, bounds);
        p.setDistName(name());
        p.setSeed(seed());

        return std::make_tuple(p, b);
    }

    virtual ~GeneratorBase() = default;

    virtual std::string name() const = 0;
    virtual tIndex seed() const = 0;
    virtual void setSeed(const tIndex &seed) = 0;

protected:
    virtual std::tuple<tPoints, tBox> _generate(const tIndex n, const tBox &bounds) = 0;

public:
    static tPoints rescalePoints(const tPoints &inPoints, const tBox &oldBounds, const tBox &newBounds) {
        tPoints outPoints;
        outPoints.reserve(inPoints.size());

        for (const auto &i : inPoints) {

            tIFloatVector p;
            for (tDim d = 0; d < D; ++d) {
                p[d] = newBounds.low[d] + (newBounds.high[d] - newBounds.low[d]) * ((i[d] - oldBounds.low[d]) / static_cast<tIFloat>((oldBounds.high[d] - oldBounds.low[d])));
            }

            ASSERT(newBounds.contains(p));
            outPoints.emplace_back(p);
        }

        return outPoints;
    }
};

template<typename Dist>
class Generator : public GeneratorBase {

public:
    Generator(const tIndex &seed) : gen(seed), _seed(seed) {}

protected:
    auto rand() {
        return static_cast<Dist &>(*this).dist(gen);
    }

    // See https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
    std::vector<tIndex> FisherYatesShuffle(const tIndex &sample_size, const tIndex &pop_size) {
        if (sample_size > pop_size) {
            std::vector<tIndex> sample(pop_size);
            std::iota(sample.begin(), sample.end(), 0);

            return sample;
        }

        std::vector<tIndex> sample(sample_size);

        for (tIndex i = 0; i != pop_size; ++i) {
            auto lDist = std::uniform_int_distribution<tIndex>(0, i);
            tIndex j = lDist(gen);
            if (j < sample.size()) {
                if (i < sample.size()) {
                    sample[i] = sample[j];
                }
                sample[j] = i;
            }
        }
        return sample;
    }

    template<typename T, typename Compare>
    std::vector<tIndex> sort_indices(const std::vector<T> &v, Compare comp) {

        // initialize original index locations
        std::vector<tIndex> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);

        // sort indexes based on comparing values in v
        // using std::stable_sort instead of std::sort
        // to avoid unnecessary index re-orderings
        // when v contains elements of equal values
        std::stable_sort(idx.begin(), idx.end(),
                         [&v, &comp](tIndex i1, tIndex i2) { return comp(v[i1], v[i2]); });

        return idx;
    }

    template<typename T, tIndex D2>
    std::tuple<tPoints, tBox> extract_points(const tIndex &n, const std::vector<std::array<T, D2>> &inPoints,
                                             std::array<std::vector<tIndex>, D2> coordSort) {
        ASSERT(inPoints.size() > n);

        // sample random point to start with
        std::array<tIndex, D2> seedPoint;
        for (tDim d = 0; d < D2; ++d) {
            seedPoint[d] = static_cast<tIndex>(std::floor((.25 + .5 * rand()) * coordSort[d].size()));
        }

        tIFloat windowSize = .0;
        tIFloat windowSizeInc = .01;
        const tIFloat maxWindowSize = .5;

        std::vector<tIndex> pointIdx;

        while (pointIdx.size() < n && windowSize <= maxWindowSize) {
            pointIdx.clear();//TODO: this can be done smarter

            windowSize += windowSizeInc;

            std::array<std::vector<tIndex>, D2> windowIdx;

            for (tDim d = 0; d < D2; ++d) {

                auto begin = coordSort[d].begin();
                std::advance(begin,
                             static_cast<tIndex>(std::floor(std::max(0.0,
                                                                     seedPoint[X] - (coordSort[d].size() * windowSize)))));

                auto end = coordSort[d].begin();
                std::advance(end,
                             static_cast<tIndex>(std::floor(std::min(static_cast<double>(coordSort[d].size()),
                                                                     seedPoint[X] + (coordSort[d].size() * windowSize)))));

                std::copy(begin, end, std::back_inserter(windowIdx[d]));
                std::sort(windowIdx[d].begin(), windowIdx[d].end());

                if (windowIdx[d].size() < n) {
                    continue;
                }
            }

            std::array<std::vector<tIndex>, D2> windowIdxIS;
            windowIdxIS[0] = windowIdx[0];

            for (tDim d = 1; d < D2; ++d) {
                std::set_intersection(windowIdxIS[d - 1].begin(), windowIdxIS[d - 1].end(), windowIdx[d].begin(), windowIdx[d].end(), std::back_inserter(windowIdxIS[d]));
            }

            pointIdx = windowIdxIS[D2 - 1];

            if (pointIdx.size() > 1.05 * n) {
                windowSize -= windowSizeInc;
                windowSizeInc /= 2;
                pointIdx.clear();
            }
        }

        tPoints points;
        points.reserve(pointIdx.size());

        tIFloatVector minPoint = {std::numeric_limits<int>::max(), std::numeric_limits<int>::max()};
        tIFloatVector maxPoint = {std::numeric_limits<int>::min(), std::numeric_limits<int>::min()};

        for (const auto &i : pointIdx) {
            tIFloatVector p;
            for (tDim d = 0; d < D; ++d) {
                p[d] = inPoints[i][d];
            }

            if (p[X] < minPoint[X]) minPoint[X] = p[X];
            if (p[Y] < minPoint[Y]) minPoint[Y] = p[Y];
            if (p[X] > maxPoint[X]) maxPoint[X] = p[X];
            if (p[Y] > maxPoint[Y]) maxPoint[Y] = p[Y];

            points.emplace_back(p);
        }

        return std::make_tuple(points, tBox{minPoint, maxPoint});
    }

public:
    tIndex seed() const override {
        return _seed;
    }

    void setSeed(const tIndex &seed) override {
        _seed = seed;
        gen.seed(seed);
    }


private:
    std::mt19937 gen;
    tIndex _seed;
};
