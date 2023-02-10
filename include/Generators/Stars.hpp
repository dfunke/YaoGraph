//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include "GeneratorBase.hpp"

class Stars : public Generator<Stars> {

    friend class Generator<Stars>;

public:
    Stars(const tIndex &seed) : Generator(seed) {

        //std::string fileName = STAR_DATA "/gaia_1331909727.points";
        std::string fileName = STAR_DATA "/gaia_1.points";
        std::ifstream file(fileName, std::ios::in);
        std::string line;

        while (std::getline(file, line)) {
            // star line
            // x y z

            std::istringstream iss(line);
            double x, y, z;
            if (!(iss >> x >> y >> z)) { continue; }// error in file

            filePoints.push_back({x, y, z});
        }

        for (tDim d = 0; d < 3; ++d) {
            coordSort[d] = sort_indices(filePoints, [d](const tFilePoint &a, const tFilePoint &b) { return a[d] < b[d]; });
        }
    }

    std::string name() const override {
        return "stars";
    }

    std::tuple<tPoints, tBox> _generate(const tIndex n, [[maybe_unused]] const tBox &bounds) override {

        return extract_points(n, filePoints, coordSort);
    }

private:
    std::uniform_real_distribution<tIFloat> dist;

    using tFilePoint = std::array<double, 3>;
    std::vector<tFilePoint> filePoints;
    std::array<std::vector<tIndex>, 3> coordSort;
};