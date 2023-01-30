#include <filesystem>
#include <fstream>
#include <iostream>

#include "contrib/popl.hpp"

#include "Utils/ASSERT.hpp"

#include "Generators/Generators.hpp"
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
const char Dists[] = {'u', 'g', 'd', 'r', 's', 'c', 'b'};

void generate(GeneratorBase &gen, const tIndex &n, const tBox &iBounds) {

    std::cout << "Generating " << gen.name() << " with " << n << " points" << std::endl;

    auto [points, oBounds] = gen.generate(n, iBounds);
    ASSERT(points.size() >= n);

    std::stringstream dir;
    dir << DATA_DIR << "/" << gen.name();

    std::filesystem::path p(dir.str());
    std::filesystem::create_directories(p);

    std::ofstream file(p.append("points_" + gen.name() + "_" + std::to_string(n) + "_" + std::to_string(gen.seed()) + ".csv"), std::ios::out | std::ios::trunc);

    // header
    file << "# n " << points.size() << std::endl;
    file << "# b " << oBounds.low[X] << " " << oBounds.low[Y] << " " << oBounds.high[X] << " " << oBounds.high[Y] << std::endl;

    for (auto &p : points) {
        file << p[X] << " " << p[Y] << std::endl;
    }
}

int main(int argc, char **argv) {

    popl::OptionParser op("Point Generator");
    auto sHelp = op.add<popl::Switch>("h", "help", "produce help message");

    // generate points
    auto sAllDists = op.add<popl::Switch>("-a", "--all", "generate all available distributions");
    auto oN = op.add<popl::Value<tIndex>>("n", "n", "number of points to generate");
    auto oMinN = op.add<popl::Value<tIndex>>("", "minN", "minimum number of points to generate", minN);
    auto oMaxN = op.add<popl::Value<tIndex>>("", "maxN", "maxium number of points to generate", maxN);
    auto oSeed = op.add<popl::Value<tIndex>>("s", "seed", "seed for RNG", Seeds[0]);
    auto oInst = op.add<popl::Value<tIndex>>("i", "inst", "number of instances");

    op.parse(argc, argv);

    if (sHelp->is_set()) {
        std::cout << op << "\n";
        std::cout << "list of distributions to generate [_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar, _c_ircle, _b_ubbles]" << std::endl;

        return 0;
    }

    std::vector<tIndex> seeds;
    if (oInst->is_set()) {
        // if oInst is set, oSeed is ignored
        for (uint i = 0; i < oInst->value(); ++i) {
            seeds.push_back(Seeds[i]);
        }
    } else {
        seeds.push_back(oSeed->value());
    }

    std::vector<char> dists;
    if (sAllDists->is_set()) {
        for (auto d : Dists) {
            dists.push_back(d);
        }
    } else if (op.non_option_args().size()) {
        for (uint i = 0; i < op.non_option_args().size(); ++i) {
            dists.push_back(op.non_option_args()[i][0]);
        }
    } else {
        std::cout << "Distribution must be specified" << std::endl;
    }

    if (oN->is_set()) {
        oMinN->set_value(oN->value());
        oMaxN->set_value(oN->value());
    }

    for (auto g : dists) {
        auto gen = getGen(g, 0);

        if (!gen) {
            std::cout << "Invalid distribution specified: " << g << std::endl;
            continue;
        }

        for (auto s : seeds) {
            for (tIndex nPoints = oMinN->value(); nPoints <= oMaxN->value(); nPoints += 3 * pow(10, floor(log10(nPoints)))) {
                gen->setSeed(s);
                generate(*gen, nPoints, BOUNDS);
            }
        }
    }
}
