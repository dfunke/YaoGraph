#include <filesystem>
#include <fstream>
#include <iostream>

#include "Utils/ASSERT.hpp"

#include "Generators.hpp"
#include "Types.hpp"

// constants
constexpr tBox BOUNDS{{0, 0},
                      {1, 1}};
constexpr tIndex minN = 1e3;
constexpr tIndex maxN = 1e7;
constexpr tDim Cones = 6;
constexpr tIndex cellOcc = 1e2;
constexpr tDim RepsPerI = 3;
constexpr tDim RepsPerN = 3;

const tIndex Seeds[] = {8158, 14030, 18545, 20099, 24065, 35700, 37197, 38132, 59135, 60315};
const char Dists[] = {'u', 'g', 'd', 'r', 's'};

void generate(const char &dist, const tIndex &n, const tIndex &seed, const tBox &iBounds) {

    auto gen = getGen(dist, seed);

    std::cout << "Generating " << gen->name() << " with " << n << " points" << std::endl;

    auto [points, oBounds] = gen->generate(n, iBounds);
    ASSERT(points.size() >= n);

    std::stringstream dir;
    dir << DATA_DIR << "/" << gen->name();

    std::filesystem::path p(dir.str());
    std::filesystem::create_directories(p);

    std::ofstream file(p.append("points_" + gen->name() + "_" + std::to_string(n) + "_" + std::to_string(seed) + ".csv"), std::ios::out | std::ios::trunc);

    // header
    file << "# n " << points.size() << std::endl;
    file << "# b " << oBounds.low[X] << " " << oBounds.low[Y] << " " << oBounds.high[X] << " " << oBounds.high[Y] << std::endl;

    for (auto &p : points) {
        file << p[X] << " " << p[Y] << std::endl;
    }
}

int main(int argc, char **argv) {
    for (auto g : Dists) {
        for (auto s : Seeds) {
            for (tIndex nPoints = minN; nPoints <= maxN; nPoints += 3 * pow(10, floor(log10(nPoints)))) {
                generate(g, nPoints, s, BOUNDS);
            }
        }
    }
}
