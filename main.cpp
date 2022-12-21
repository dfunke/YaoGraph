#include <fstream>
#include <iostream>

#include "contrib/popl.hpp"

#include "Utils/ASSERT.hpp"

#include "Generators.hpp"
#include "Types.hpp"
#include "Utils/Timer.hpp"

#include "Algorithms/GridYao.hpp"
#include "Algorithms/NaiveYao.hpp"
#include "Algorithms/SweepLine.hpp"

#ifdef WITH_CAIRO
#include "Painter.hpp"
#endif

#ifdef WITH_CGAL
#include "Algorithms/CGAL_Yao.hpp"
#include "Kernels/CGALKernel.hpp"
#endif

// constants
constexpr tBox BOUNDS{{0, 0},
                      {1, 1}};
constexpr tIndex minN = 1e3;
constexpr tIndex maxN = 1e5;
constexpr tDim Cones = 6;
constexpr tIndex cellOcc = 1e2;
constexpr tDim RepsPerI = 3;
constexpr tDim RepsPerN = 3;

const tIndex Seeds[] = {8158, 14030, 18545, 20099, 24065, 35700, 37197, 38132, 59135, 60315};

#define DISTS Gaussian, Uniform, Grid, Road, Stars

std::unique_ptr<GeneratorBase> getGen(const char &dist, const tIndex& seed) {
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

std::tuple<tPoints, tBox> readPointsFile(const std::string &fileName) {
    std::ifstream file(fileName, std::ios::in);
    std::string line;

    tPoints points;
    tIFloatVector minPoint = {std::numeric_limits<tIFloat>::max(), std::numeric_limits<tIFloat>::max()};
    tIFloatVector maxPoint = {std::numeric_limits<tIFloat>::min(), std::numeric_limits<tIFloat>::min()};

    while (std::getline(file, line)) {
        // vertex line
        // X Y

        std::istringstream iss(line);
        tIFloat x, y;
        if (!(iss >> x >> y)) { continue; }// error in file

        points.push_back({x, y});

        if (x < minPoint[X]) minPoint[X] = x;
        if (y < minPoint[Y]) minPoint[Y] = y;
        if (x > maxPoint[X]) maxPoint[X] = x;
        if (y > maxPoint[Y]) maxPoint[Y] = y;
    }

    return std::make_tuple(points, tBox{minPoint, maxPoint});
}

template<typename Algorithm, typename... Args>
auto runAlg(const tDim &K, const tPoints &points, const tBox &bounds, Args... args) {
    std::cout << "Generating Yao graph of " << points.size() << " points with " << Algorithm::name() << std::endl;
    auto result = Timer<Algorithm>::time(K, points, bounds, args...);
    std::cout << "Time: " << std::get<0>(result) << "ms" << std::endl;

    return std::get<1>(result);
}

template<typename Algorithm, typename Distribution, typename... Args>
void benchmarkImpl(Args... args) {
    std::ofstream file("benchmark_" + Algorithm::name() + ".csv", std::ios::out | std::ios::app);

    std::cout << "Benchmarking " << Algorithm::name() << " with " << Distribution::name() << " distribution" << std::endl;

    // header
    file << "# dist n seed rep t" << std::endl;
    std::cout << "dist n seed rep t" << std::endl;

    for (tIndex nPoints = minN; nPoints <= maxN; nPoints += 3 * pow(10, floor(log10(nPoints)))) {
        for (tDim rpn = 0; rpn < RepsPerN; ++rpn) {

            Distribution gen(Seeds[rpn]);
            auto points = gen.generate(nPoints, BOUNDS);

            for (tDim rpi = 0; rpi < RepsPerI; ++rpi) {
                auto result = Timer<Algorithm>::time(Cones, points, BOUNDS, args...);
                file << Distribution::name() << " " << nPoints << " " << Seeds[rpn] << " " << rpi << " " << std::get<0>(result) << std::endl;
                std::cout << Distribution::name() << " " << nPoints << " " << Seeds[rpn] << " " << rpi << " " << std::get<0>(result) << std::endl;
            }
        }
    }
}

template<typename Algorithm, typename Distribution, typename... FDists, typename... Args>
void benchmark(Args... args) {
    benchmarkImpl<Algorithm, Distribution>(args...);
    if constexpr (sizeof...(FDists) > 0) {
        benchmark<Algorithm, FDists...>(args...);
    }
}

void benchmark() {

#ifdef WITH_CGAL
    benchmark<CGAL_Yao2D<ExactPredicatesInexactConstructions>, DISTS>();
    benchmark<CGAL_Yao2D<ExactPredicatesExactConstructions>, DISTS>();
#endif

    benchmark<NaiveYao<InexactKernel>, DISTS>();
#ifdef WITH_CGAL
    benchmark<NaiveYao<CGALKernel<ExactPredicatesInexactConstructions>>, DISTS>();
    benchmark<NaiveYao<CGALKernel<ExactPredicatesExactConstructions>>, DISTS>();
#endif

    benchmark<GridYao<InexactKernel>, DISTS>(cellOcc);
#ifdef WITH_CGAL
    benchmark<GridYao<CGALKernel<ExactPredicatesInexactConstructions>>, DISTS>(cellOcc);
    benchmark<GridYao<CGALKernel<ExactPredicatesExactConstructions>>, DISTS>(cellOcc);
#endif

    benchmark<SweepLine<InexactKernel>, DISTS>();
#ifdef WITH_CGAL
    benchmark<SweepLine<CGALKernel<ExactPredicatesInexactConstructions>>, DISTS>();
    benchmark<SweepLine<CGALKernel<ExactPredicatesExactConstructions>>, DISTS>();
#endif
}

int main(int argc, char **argv) {
    popl::OptionParser op("Yao Graph generator");
    auto sHelp = op.add<popl::Switch>("h", "help", "produce help message");

    // benchmark
    auto sBenchmark = op.add<popl::Switch>("b", "benchmark", "run benchmark suite");

#ifndef NDEBUG
    auto sDebug = op.add<popl::Switch>("x", "debug", "find offending instance");
#endif

    // generate points
    auto oDist = op.add<popl::Value<char>>("d", "dist", "point distribution [_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar]", 'u');
    auto oN = op.add<popl::Value<tIndex>>("n", "n", "number of points to generate");
    auto oSeed = op.add<popl::Value<tIndex>>("s", "seed", "seed for RNG", Seeds[0]);

    // points file
    auto oPointsFile = op.add<popl::Value<std::string>>("f", "infile", "file with points");

    // algorithm to use
    auto oAlg = op.add<popl::Value<char>>("a", "alg", "algorithm to use [_s_weepline, _g_rid, _n_aive]", 's');
    auto oKern = op.add<popl::Value<char>>("c", "kernel", "kernel to use [_i_nexact, CGALExact_p_redicatesInexactConstructions, CGALExactpredicatesInexact_c_onstructions]", 'i');
    auto oK = op.add<popl::Value<tDim>>("k", "cones", "number of cones, default 6", 6);
    auto oCellOcc = op.add<popl::Value<tIndex>>("g", "cellOcc", "number of points per cell (grid algorithm)", cellOcc);

    // output
    auto oOutFile = op.add<popl::Value<std::string>>("o", "outfile", "file for graph output");
    auto sStdOut = op.add<popl::Switch>("p", "stdout", "write graph to stdout");

    op.parse(argc, argv);

    if (sHelp->is_set()) {
        std::cout << op << "\n";

        return 0;
    }

    if (sBenchmark->is_set()) {
        benchmark();

        return 0;
    }

#ifndef NDEBUG
    if (sDebug->is_set()) {
        ASSERT(oN->is_set());

        tIndex seed = 1;
        while (seed != 0) {
            Uniform gen(seed);
            auto points = gen.generate(oN->value(), BOUNDS);

            std::cout << "seed: " << seed << " n: " << points.size() << std::endl;

            SweepLine<CGALKernel<ExactPredicatesInexactConstructions>> alg;
            auto res = alg(Cones, points, BOUNDS);

            seed++;
        }

        return 0;
    }
#endif

    // get points
    tPoints points;
    tBox bounds;
    if (oN->is_set()) {
        auto gen = getGen(oDist->value(), oSeed->value());
        points = gen->generate(oN->value(), BOUNDS);
        bounds = BOUNDS;
    } else if (oPointsFile->is_set()) {
        std::tie(points, bounds) = readPointsFile(oPointsFile->value());
    } else {
        std::cout << "Either specify input point file or point generation option" << std::endl;
        return 0;
    }

    // run algorithm
    tYaoGraph graph(points.size(), oK->value());
    switch (oAlg->value()) {
        case 's':
            switch (oKern->value()) {
                case 'i':
                    graph = runAlg<SweepLine<InexactKernel>>(oK->value(), points, bounds);
                    break;
                case 'p':
                    graph = runAlg<SweepLine<CGALKernel<ExactPredicatesInexactConstructions>>>(oK->value(), points, bounds);
                    break;
                case 'c':
                    graph = runAlg<SweepLine<CGALKernel<ExactPredicatesExactConstructions>>>(oK->value(), points, bounds);
                    break;
            }

            break;
        case 'g':
            switch (oKern->value()) {
                case 'i':
                    graph = runAlg<GridYao<InexactKernel>>(oK->value(), points, bounds, oCellOcc->value());
                    break;
                case 'p':
                    graph = runAlg<GridYao<CGALKernel<ExactPredicatesInexactConstructions>>>(oK->value(), points, bounds, oCellOcc->value());
                    break;
                case 'c':
                    graph = runAlg<GridYao<CGALKernel<ExactPredicatesExactConstructions>>>(oK->value(), points, bounds, oCellOcc->value());
                    break;
            }

            break;
        case 'n':
            switch (oKern->value()) {
                case 'i':
                    graph = runAlg<NaiveYao<InexactKernel>>(oK->value(), points, bounds);
                    break;
                case 'p':
                    graph = runAlg<NaiveYao<CGALKernel<ExactPredicatesInexactConstructions>>>(oK->value(), points, bounds);
                    break;
                case 'c':
                    graph = runAlg<NaiveYao<CGALKernel<ExactPredicatesExactConstructions>>>(oK->value(), points, bounds);
                    break;
            }
    }

    if (oOutFile->is_set()) {
        std::ofstream file(oOutFile->value(), std::ios::out | std::ios::trunc);
        file << graph << std::endl;
    }

    if (sStdOut->is_set()) {
        std::cout << graph << std::endl;
    }
}
