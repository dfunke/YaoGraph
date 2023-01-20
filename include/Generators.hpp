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

    template<typename T, tIndex D2, typename Trans>
    std::tuple<tPoints, tBox> extract_points(const tIndex &n, const std::vector<std::array<T, D2>> &inPoints, Trans &trans) {
        ASSERT(inPoints.size() > n);

        std::array<std::vector<tIndex>, D2> coordSort;

        for (tDim d = 0; d < D2; ++d) {
            coordSort[d] = sort_indices(inPoints, [d](const std::array<T, D2> &a, const std::array<T, D2> &b) { return a[d] < b[d]; });
        }

        // sample random point to start with
        std::array<tIndex, D2> seedPoint;
        for (tDim d = 0; d < D2; ++d) {
            seedPoint[d] = static_cast<tIndex>(std::floor((.25 + .5 * rand()) * coordSort[d].size()));
        }

        tIFloat windowSize = .0;
        const tIFloat windowSizeInc = .01;
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
            if (pointIdx.size() < n) {
                continue;
            }
        }

        tPoints points;
        points.reserve(pointIdx.size());

        tIFloatVector minPoint = {std::numeric_limits<int>::max(), std::numeric_limits<int>::max()};
        tIFloatVector maxPoint = {std::numeric_limits<int>::min(), std::numeric_limits<int>::min()};

        for (const auto &i : pointIdx) {
            tIFloatVector p;
            for (tDim d = 0; d < D; ++d) {
                p[d] = trans(inPoints[i][d]);
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

private:
    std::mt19937 gen;
    tIndex _seed;
};

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

class Road : public Generator<Road> {

    friend class Generator<Road>;

public:
    Road(const tIndex &seed) : Generator(seed) {}

    std::string name() const override {
        return "road";
    }

    std::tuple<tPoints, tBox> _generate(const tIndex n, [[maybe_unused]] const tBox &bounds) override {

        std::string fileName = ROAD_DATA "/USA-road-d.USA.co";
        //std::string fileName = ROAD_DATA "/USA-road-d.NY.co";
        std::ifstream file(fileName, std::ios::in);
        std::string line;

        using tIntPoint = std::array<int, 2>;
        std::vector<tIntPoint> inPoints;

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

                inPoints[idx - 1] = {x, y};
            }
        }

        auto trans = [](const int i) {
            return i / 1e6;
        };

        return extract_points(n, inPoints, trans);
    }

private:
    std::uniform_int_distribution<tIndex> dist;
};

class Stars : public Generator<Stars> {

    friend class Generator<Stars>;

public:
    Stars(const tIndex &seed) : Generator(seed) {}

    std::string name() const override {
        return "stars";
    }

    std::tuple<tPoints, tBox> _generate(const tIndex n, [[maybe_unused]] const tBox &bounds) override {

        std::string fileName = STAR_DATA "/gaia_1331909727.points";
        //std::string fileName = STAR_DATA "/gaia_1.points";
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

        auto trans = [](const double i) {
            return i;
        };

        return extract_points(n, inPoints, trans);
    }

private:
    std::uniform_real_distribution<tIFloat> dist;
};

std::unique_ptr<GeneratorBase> getGen(const char &dist, const tIndex &seed) {
    // [_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar]

    switch (dist) {
        case 'u':
            return std::make_unique<Uniform>(seed);
        case 'g':
            return std::make_unique<Gaussian>(seed);
        case 'd':
            return std::make_unique<Grid>(seed);
        case 'r':
            return std::make_unique<Road>(seed);
        case 's':
            return std::make_unique<Stars>(seed);

        default:
            return nullptr;
    }
}

std::string getGenName(const char &dist) {
    // [_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar]

    switch (dist) {
        case 'u':
            return "uni";
        case 'g':
            return "gaussian";
        case 'd':
            return "grid";
        case 'r':
            return "road";
        case 's':
            return "stars";

        default:
            return nullptr;
    }
}