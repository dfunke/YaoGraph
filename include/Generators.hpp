//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include <fstream>
#include <random>
#include <sstream>
#include <algorithm>

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
                p[X] = bounds.low[X] + pointSpacing[X] * i;
                p[Y] = bounds.low[Y] + pointSpacing[Y] * j;

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

        std::string fileName = "/Users/funke/Development/GeoGraph/data/USA-road-d.NY.co";
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

        for (const auto & i : samples) {

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
    // See https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
    std::vector<tIndex> FisherYatesShuffle(const tIndex &sample_size, const tIndex &pop_size) {
        if(sample_size > pop_size){
            std::vector<tIndex> sample(pop_size);
            std::iota(sample.begin(), sample.end(), 0);

            return sample;
        }

        std::vector<tIndex> sample(sample_size);

        for (tIndex i = 0; i != pop_size; ++i) {
            dist = std::uniform_int_distribution<tIndex>(0, i);
            tIndex j = rand();
            if (j < sample.size()) {
                if (i < sample.size()) {
                    sample[i] = sample[j];
                }
                sample[j] = i;
            }
        }
        return sample;
    }

private:
    std::uniform_int_distribution<tIndex> dist;
};