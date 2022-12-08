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

std::unique_ptr<GeneratorBase> getGen(const char &dist) {
    // [_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar]

    switch (dist) {
        case 'u':
            return std::make_unique<Uniform>(Seeds[0]);
        case 'g':
            return std::make_unique<Gaussian>(Seeds[0]);
        case 'd':
            return std::make_unique<Grid>(Seeds[0]);
        case 'r':
            return std::make_unique<Road>(Seeds[0]);
        case 's':
            return std::make_unique<Stars>(Seeds[0]);

        default:
            return nullptr;
    }
};

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

template<typename Algorithm>
auto runAlg(const tPoints &points, const tBox &bounds) {
    std::cout << "Generating Yao graph of " << points.size() << " points with " << Algorithm::name() << std::endl;
    auto result = Timer<Algorithm>::time(points, bounds);
    std::cout << "Time: " << std::get<0>(result) << "ms" << std::endl;

    return std::get<1>(result);
}

template<typename Algorithm, typename Distribution>
void benchmarkImpl() {
    std::ofstream file("benchmark_" + Algorithm::name() + ".csv", std::ios::out | std::ios::app);

    std::cout << "Benchmarking " << Algorithm::name() << std::endl;
    
    // header
    file << "# dist n seed rep t" << std::endl;
    std::cout << "dist n seed rep t" << std::endl;

    for (tIndex nPoints = minN; nPoints <= maxN; nPoints += 3 * pow(10, floor(log10(nPoints)))) {
        for (tDim rpn = 0; rpn < RepsPerN; ++rpn) {

            Distribution gen(Seeds[rpn]);
            auto points = gen.generate(nPoints, BOUNDS);

            for (tDim rpi = 0; rpi < RepsPerI; ++rpi) {
                auto result = Timer<Algorithm>::time(points, BOUNDS);
                file << Distribution::name() << " " << nPoints << " " << Seeds[rpn] << " " << rpi << " " << std::get<0>(result) << std::endl;
                std::cout << Distribution::name() << " " << nPoints << " " << Seeds[rpn] << " " << rpi << " " << std::get<0>(result) << std::endl;
            }
        }
    }
}

template<typename Algorithm, typename Distribution, typename... FDists>
void benchmark() {
    benchmarkImpl<Algorithm, Distribution>();
    if constexpr (sizeof...(FDists) > 0) {
        benchmark<Algorithm, FDists...>();
    }
}

void benchmark() {

#ifdef WITH_CGAL
    benchmark<CGAL_Yao2D<Cones, ExactPredicatesInexactConstructions>, DISTS>();
    benchmark<CGAL_Yao2D<Cones, ExactPredicatesExactConstructions>, DISTS>();
#endif

    benchmark<NaiveYao<Cones, InexactKernel>, DISTS>();
#ifdef WITH_CGAL
    benchmark<NaiveYao<Cones, CGALKernel<ExactPredicatesInexactConstructions>>, DISTS>();
    benchmark<NaiveYao<Cones, CGALKernel<ExactPredicatesExactConstructions>>, DISTS>();
#endif

    benchmark<GridYao<Cones, InexactKernel, cellOcc>, DISTS>();
#ifdef WITH_CGAL
    benchmark<GridYao<Cones, CGALKernel<ExactPredicatesInexactConstructions>, cellOcc>, DISTS>();
    benchmark<GridYao<Cones, CGALKernel<ExactPredicatesExactConstructions>, cellOcc>, DISTS>();
#endif

    benchmark<SweepLine<Cones, InexactKernel>, DISTS>();
#ifdef WITH_CGAL
    benchmark<SweepLine<Cones, CGALKernel<ExactPredicatesInexactConstructions>>, DISTS>();
    benchmark<SweepLine<Cones, CGALKernel<ExactPredicatesExactConstructions>>, DISTS>();
#endif
}

int main(int argc, char **argv) {
    popl::OptionParser op("Yao Graph generator");
    auto sHelp = op.add<popl::Switch>("h", "help", "produce help message");

    // benchmark
    auto sBenchmark = op.add<popl::Switch>("b", "benchmark", "run benchmark suite");

    // generate points
    auto oDist = op.add<popl::Value<char>>("d", "dist", "point distribution [_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar]");
    auto oN = op.add<popl::Value<tIndex>>("n", "n", "number of points to generate");

    // points file
    auto oPointsFile = op.add<popl::Value<std::string>>("f", "file", "file with points");

    // algorithm to use
    auto oAlg = op.add<popl::Value<char>>("a", "alg", "algorithm to use [_s_weepline, _g_rid, _n_aive]", 's');
    auto oKern = op.add<popl::Value<char>>("k", "kernel", "kernel to use [_i_nexact, CGALExact_p_redicatesInexactConstructions, CGALExactpredicatesInexact_c_onstructions]", 'i');

    op.parse(argc, argv);

    if (sHelp->is_set()) {
        std::cout << op << "\n";

        return 0;
    }

    if (sBenchmark->is_set()) {
        benchmark();

        return 0;
    }

    // get points
    tPoints points;
    tBox bounds;
    if (oDist->is_set() && oN->is_set()) {
        auto gen = getGen(oDist->value());
        points = gen->generate(oN->value(), BOUNDS);
        bounds = BOUNDS;
    } else if (oPointsFile->is_set()) {
        std::tie(points, bounds) = readPointsFile(oPointsFile->value());
    }

    // run algorithm
    switch (oAlg->value()) {
        case 's':
            switch (oKern->value()) {
                case 'i':
                    runAlg<SweepLine<Cones, InexactKernel>>(points, bounds);
                    break;
                case 'p':
                    runAlg<SweepLine<Cones, CGALKernel<ExactPredicatesInexactConstructions>>>(points, bounds);
                    break;
                case 'c':
                    runAlg<SweepLine<Cones, CGALKernel<ExactPredicatesExactConstructions>>>(points, bounds);
                    break;
            }

            break;
        case 'g':
            switch (oKern->value()) {
                case 'i':
                    runAlg<GridYao<Cones, InexactKernel, cellOcc>>(points, bounds);
                    break;
                case 'p':
                    runAlg<GridYao<Cones, CGALKernel<ExactPredicatesInexactConstructions>, cellOcc>>(points, bounds);
                    break;
                case 'c':
                    runAlg<GridYao<Cones, CGALKernel<ExactPredicatesExactConstructions>, cellOcc>>(points, bounds);
                    break;
            }

            break;
        case 'n':
            switch (oKern->value()) {
                case 'i':
                    runAlg<NaiveYao<Cones, InexactKernel>>(points, bounds);
                    break;
                case 'p':
                    runAlg<NaiveYao<Cones, CGALKernel<ExactPredicatesInexactConstructions>>>(points, bounds);
                    break;
                case 'c':
                    runAlg<NaiveYao<Cones, CGALKernel<ExactPredicatesExactConstructions>>>(points, bounds);
                    break;
            }
    }
}
