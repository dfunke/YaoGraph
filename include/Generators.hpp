//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include <algorithm>
#include <fstream>
#include <numeric>
#include <random>
#include <sstream>

#include "Predicates.hpp"
#include "Types.hpp"

template<typename Dist>
class Generator {

public:
    Generator(const tIndex &seed) : gen(seed) {}
    virtual tPoints generate(const tIndex n, const tBox &bounds) = 0;

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

    tPoints generate(const tIndex n, const tBox &bounds) override {
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

    tPoints generate(const tIndex n, const tBox &bounds) override {
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

    tPoints generate(const tIndex n, const tBox &bounds) override {
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

        return points;
    }
};

class Road : public Generator<Road> {

    friend class Generator<Road>;

public:
    Road(const tIndex &seed) : Generator(seed) {}

    static std::string name() {
        return "road";
    }

    tPoints generate(const tIndex n, const tBox &bounds) override {

        std::string fileName = ROAD_DATA "/USA-road-d.NY.co";
        std::ifstream file(fileName, std::ios::in);
        std::string line;

        using tIntPoint = std::array<int, 2>;
        std::vector<tIntPoint> inPoints;
        tIntPoint minPoint = {std::numeric_limits<int>::max(), std::numeric_limits<int>::max()};
        tIntPoint maxPoint = {std::numeric_limits<int>::min(), std::numeric_limits<int>::min()};

        while (std::getline(file, line)) {
            if (line.starts_with('c')) {
                // comment line
                continue;
            }
            if (line.starts_with('p')) {
                // problem line, occurs before any input points by definition
                // p aux co N
                auto lastSpace = line.find_last_of(' ');
                tIndex pointsInFile = std::stoul(line.substr(lastSpace + 1));
                inPoints.resize(pointsInFile);
            }
            if (line.starts_with('v')) {
                // vertex line
                // v IDX X Y

                std::istringstream iss(line);
                char c;
                int idx, x, y;
                if (!(iss >> c >> idx >> x >> y)) { continue; }// error in file

                inPoints[idx] = {x, y};

                if (x < minPoint[X]) minPoint[X] = x;
                if (y < minPoint[Y]) minPoint[Y] = y;
                if (x > maxPoint[X]) maxPoint[X] = x;
                if (y > maxPoint[Y]) maxPoint[Y] = y;
            }
        }

        tPoints points;
        points.reserve(n);

        auto samples = FisherYatesShuffle(n, inPoints.size());

        for (const auto &i : samples) {

            tIFloatVector p;
            for (tDim d = 0; d < D; ++d) {
                p[d] = bounds.low[d] + (bounds.high[d] - bounds.low[d]) * ((inPoints[i][d] - minPoint[d]) / static_cast<tIFloat>((maxPoint[d] - minPoint[d])));
            }

            ASSERT(bounds.contains(p));
            points.emplace_back(p);
        }

        return points;
    }

private:
    std::uniform_int_distribution<tIndex> dist;
};

class Stars : public Generator<Stars> {

    friend class Generator<Stars>;

public:
    Stars(const tIndex &seed) : Generator(seed) {}

    static std::string name() {
        return "stars";
    }

    tPoints generate(const tIndex n, const tBox &bounds) override {

        std::string fileName = STAR_DATA "/gaia_1.points";
        std::ifstream file(fileName, std::ios::in);
        std::string line;

        using tInPoint = std::array<double, 3>;
        std::vector<tInPoint> inPoints;

        while (std::getline(file, line)) {
            // star line
            // x y z

            std::istringstream iss(line);
            double x, y, z;
            if (!(iss >> x >> y >> z)) { continue; }// error in file

            inPoints.push_back({x, y, z});
        }

        auto xSort = sort_indices(inPoints, [](const tInPoint &a, const tInPoint &b) { return a[X] < b[X]; });
        auto ySort = sort_indices(inPoints, [](const tInPoint &a, const tInPoint &b) { return a[Y] < b[Y]; });
        auto zSort = sort_indices(inPoints, [](const tInPoint &a, const tInPoint &b) { return a[Z] < b[Z]; });

        // use 25 - 75 quantiles for X and Y to reject outliers
        double xMin = inPoints[xSort[xSort.size() * .25]][X];
        double xMax = inPoints[xSort[xSort.size() * .75]][X];

        double yMin = inPoints[ySort[xSort.size() * .25]][Y];
        double yMax = inPoints[ySort[xSort.size() * .75]][Y];

        std::vector<tInPoint> slicePoints;
        slicePoints.reserve(n);

        tInPoint minPoint = {std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max()};
        tInPoint maxPoint = {std::numeric_limits<int>::min(), std::numeric_limits<int>::min(), std::numeric_limits<int>::min()};

        // start with 45 - 55 quantile for zslice and grow if necessary
        double zSliceSize = .05;

        while (slicePoints.size() < n && zSliceSize <= .25) {
            slicePoints.clear();//TODO: this can be done smarter

            for (tIndex i = zSort[zSort.size() * (.5 - zSliceSize)]; i < zSort[zSort.size() * (.5 + zSliceSize)]; ++i) {
                const auto &p = inPoints[i];
                if (xMin <= p[X] && p[X] <= xMax && yMin <= p[Y] && p[Y] <= yMax) {
                    slicePoints.push_back(p);

                    if (p[X] < minPoint[X]) minPoint[X] = p[X];
                    if (p[Y] < minPoint[Y]) minPoint[Y] = p[Y];
                    if (p[X] > maxPoint[X]) maxPoint[X] = p[X];
                    if (p[Y] > maxPoint[Y]) maxPoint[Y] = p[Y];
                }
            }

            if (slicePoints.size() < n) {
                zSliceSize += .05;
            }
        }

        tPoints points;
        points.reserve(n);

        auto samples = FisherYatesShuffle(n, slicePoints.size());

        for (const auto &i : samples) {
            tIFloatVector p;
            for (tDim d = 0; d < D; ++d) {
                p[d] = bounds.low[d] + (bounds.high[d] - bounds.low[d]) * ((slicePoints[i][d] - minPoint[d]) / static_cast<tIFloat>((maxPoint[d] - minPoint[d])));
            }

            ASSERT(bounds.contains(p));
            points.emplace_back(p);
        }

        return points;
    }

private:
    std::uniform_int_distribution<tIndex> dist;
};