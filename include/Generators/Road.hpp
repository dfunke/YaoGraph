//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include "GeneratorBase.hpp"

class Road : public Generator<Road> {

    friend class Generator<Road>;

public:
    Road(const tIndex &seed) : Generator(seed) {

        //std::string fileName = ROAD_DATA "/USA-road-d.USA.co";
        std::string fileName = ROAD_DATA "/USA-road-d.NY.co";
        std::ifstream file(fileName, std::ios::in);
        std::string line;

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
                filePoints.resize(pointsInFile);
            }
            if (line.starts_with('v')) {
                // vertex line
                // v IDX X Y

                std::istringstream iss(line);
                char c;
                int idx, x, y;
                if (!(iss >> c >> idx >> x >> y)) { continue; }// error in file

                filePoints[idx - 1] = {x / 1e6, y / 1e6};
            }
        }

        for (tDim d = 0; d < 2; ++d) {
            coordSort[d] = sort_indices(filePoints, [d](const tFilePoint &a, const tFilePoint &b) { return a[d] < b[d]; });
        }
    }

    std::string name() const override {
        return "road";
    }

    std::tuple<tPoints, tBox> _generate(const tIndex n, [[maybe_unused]] const tBox &bounds) override {

        return extract_points(n, filePoints, coordSort);
    }

private:
    std::uniform_int_distribution<tIndex> dist;

    using tFilePoint = std::array<double, 2>;
    std::vector<tFilePoint> filePoints;
    std::array<std::vector<tIndex>, 2> coordSort;
};
