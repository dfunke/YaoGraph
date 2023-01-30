#include <fstream>
#include <iostream>

#include "contrib/popl.hpp"

#include "Utils/ASSERT.hpp"

#include "Generators/Generators.hpp"
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
constexpr tIndex cellOcc = 1e2;
constexpr tDim Reps = 1;

const tIndex Seeds[] = {8158, 14030, 18545, 20099, 24065, 35700, 37197, 38132, 59135, 60315};

std::tuple<tPoints, tBox> readPointsFile(const std::string &fileName) {
    std::ifstream file(fileName, std::ios::in);
    std::string line;

    tPoints points;
    tIFloatVector minPoint = {std::numeric_limits<tIFloat>::max(), std::numeric_limits<tIFloat>::max()};
    tIFloatVector maxPoint = {std::numeric_limits<tIFloat>::min(), std::numeric_limits<tIFloat>::min()};

    while (std::getline(file, line)) {

        if (line.starts_with("# n")) {
            // number of points line
            // # n XXX
            auto lastSpace = line.find_last_of(' ');
            tIndex pointsInFile = std::stoul(line.substr(lastSpace + 1));
            points.reserve(pointsInFile);

            continue;
        }

        if (line.starts_with("# b")) {
            // bounding box line
            // # b minX minY maxX maxY
            tIFloat minX, minY, maxX, maxY;
            std::stringstream iss(line.substr(3));

            if (!(iss >> minX >> minY >> maxX >> maxY)) { continue; }// error in file

            minPoint[X] = minX;
            minPoint[Y] = minY;
            maxPoint[X] = maxX;
            maxPoint[Y] = maxY;

            continue;
        }

        //  vertex line
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

bool gBenchmarkMode = false;

template<typename Algorithm, typename... Args>
auto runAlg(const tDim &K, const tPoints &points, const tBox &bounds, Args... args) {

    std::ofstream file;
    if (gBenchmarkMode) {
        file = std::ofstream("benchmark_" + Algorithm::name() + ".csv", std::ios::out | std::ios::app);

        std::cout << "Benchmarking " << Algorithm::name() << " with " << points.getDistName() << " distribution" << std::endl;

        // header
        file << "# k dist n seed rep t" << std::endl;
        std::cout << "k dist n seed rep t" << std::endl;
    }

    tYaoGraph graph(points.size(), K);

    tDim lReps = gBenchmarkMode ? Reps : 1;
    for (tDim rpi = 0; rpi < lReps; ++rpi) {

        if (!gBenchmarkMode) {
            std::cout << "Generating Yao graph of " << points.size() << " points with " << Algorithm::name() << std::endl;
        }

        auto result = Timer<Algorithm>::time(K, points, bounds, args...);

        if (gBenchmarkMode) {
            file << K << " " << points.getDistName() << " " << points.size() << " " << points.getSeed() << " " << rpi << " " << std::get<0>(result) << std::endl;
            std::cout << K << " " << points.getDistName() << " " << points.size() << " " << points.getSeed() << " " << rpi << " " << std::get<0>(result) << std::endl;
        } else {
            std::cout << "Time: " << std::get<0>(result) << "ms" << std::endl;
        }

        graph = std::get<1>(result);
    }

    return graph;
}

int main(int argc, char **argv) {
    popl::OptionParser op("Yao Graph generator");
    auto sHelp = op.add<popl::Switch>("h", "help", "produce help message");

    // benchmark
    auto sBenchmark = op.add<popl::Switch>("b", "benchmark", "run benchmark suite");

    // generate points
    auto oDist = op.add<popl::Value<char>>("d", "dist", "point distribution [_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar, _c_ircle, _b_ubbles]", 'u');
    auto oN = op.add<popl::Value<tIndex>>("n", "n", "number of points to generate");
    auto oSeed = op.add<popl::Value<tIndex>>("s", "seed", "seed for RNG", Seeds[0]);

    // points file
    auto oPointsFile = op.add<popl::Value<std::string>>("f", "infile", "file with points");

    // algorithm to use
    auto oAlg = op.add<popl::Value<char>>("a", "alg", "algorithm to use [_s_weepline, _g_rid, _n_aive, _c_gal]", 's');
    auto oKern = op.add<popl::Value<char>>("c", "kernel", "kernel to use [_i_nexact, CGALExact_p_redicatesInexactConstructions, CGALExactpredicatesInexact_c_onstructions]", 'i');
    auto oK = op.add<popl::Value<tDim>>("k", "cones", "number of cones, default 6", 6);
    auto oCellOcc = op.add<popl::Value<tIndex>>("g", "cellOcc", "number of points per cell (grid algorithm)", cellOcc);

    // output
    auto oOutFile = op.add<popl::Value<std::string>>("o", "outfile", "file for graph output");
    auto sStdOut = op.add<popl::Switch>("p", "stdout", "write graph to stdout");
#ifdef WITH_CAIRO
    auto oImageOut = op.add<popl::Value<std::string>>("i", "image", "write image of Yao graph to specified file");
#endif

    //verify
    auto sVerify = op.add<popl::Switch>("v", "verify", "verify graph (-v grid-based, -v -v naive yao (slow))");

    op.parse(argc, argv);

    if (sHelp->is_set()) {
        std::cout << op << "\n";

        return 0;
    }

    if (sBenchmark->is_set()) {
        gBenchmarkMode = true;
    }

    // get points
    tPoints points;
    tBox bounds;
    if (oN->is_set()) {
        auto gen = getGen(oDist->value(), oSeed->value());
        std::tie(points, bounds) = gen->generate(oN->value(), BOUNDS);
    } else if (oPointsFile->is_set()) {
        std::tie(points, bounds) = readPointsFile(oPointsFile->value());

        if (oDist->is_set()) {
            points.setDistName(getGenName(oDist->value()));
        }

        if (oSeed->is_set()) {
            points.setSeed(oSeed->value());
        }
    } else {
        std::cout << "Either specify input point file or point generation option" << std::endl;
        return 0;
    }

#ifdef WITH_CAIRO
    if (oImageOut->is_set()) {
        points = GeneratorBase::rescalePoints(points, bounds, BOUNDS);
        bounds = BOUNDS;
    }
#endif

    // run algorithm
    tYaoGraph graph(points.size(), oK->value());
    switch (oAlg->value()) {
        case 's':
            switch (oKern->value()) {
                case 'i':
                    graph = runAlg<SweepLine<InexactKernel>>(oK->value(), points, bounds);
                    break;
#ifdef WITH_CGAL
                case 'p':
                    graph = runAlg<SweepLine<CGALKernel<ExactPredicatesInexactConstructions>>>(oK->value(), points, bounds);
                    break;
                case 'c':
                    graph = runAlg<SweepLine<CGALKernel<ExactPredicatesExactConstructions>>>(oK->value(), points, bounds);
                    break;
#endif
            }

            break;
        case 'g':
            switch (oKern->value()) {
                case 'i':
                    graph = runAlg<GridYao<InexactKernel>>(oK->value(), points, bounds, oCellOcc->value());
                    break;
#ifdef WITH_CGAL
                case 'p':
                    graph = runAlg<GridYao<CGALKernel<ExactPredicatesInexactConstructions>>>(oK->value(), points, bounds, oCellOcc->value());
                    break;
                case 'c':
                    graph = runAlg<GridYao<CGALKernel<ExactPredicatesExactConstructions>>>(oK->value(), points, bounds, oCellOcc->value());
                    break;
#endif
            }

            break;
        case 'n':
            switch (oKern->value()) {
                case 'i':
                    graph = runAlg<NaiveYao<InexactKernel>>(oK->value(), points, bounds);
                    break;
#ifdef WITH_CGAL
                case 'p':
                    graph = runAlg<NaiveYao<CGALKernel<ExactPredicatesInexactConstructions>>>(oK->value(), points, bounds);
                    break;
                case 'c':
                    graph = runAlg<NaiveYao<CGALKernel<ExactPredicatesExactConstructions>>>(oK->value(), points, bounds);
                    break;
#endif
            }
#ifdef WITH_CGAL
            break;
        case 'c':
            switch (oKern->value()) {
                case 'p':
                    runAlg<CGAL_Yao2D<ExactPredicatesInexactConstructions>>(oK->value(), points, bounds);
                    break;
                case 'c':
                    runAlg<CGAL_Yao2D<ExactPredicatesExactConstructions>>(oK->value(), points, bounds);
                    break;
            }
#endif
    }

    if (oOutFile->is_set()) {
        std::ofstream file(oOutFile->value(), std::ios::out | std::ios::trunc);
        file << graph << std::endl;
    }

    if (sStdOut->is_set()) {
        std::cout << graph << std::endl;
    }

#ifdef WITH_CAIRO
    if (oImageOut->is_set()) {
        Painter painter(bounds, 1000);
        painter.draw(points, false);
        painter.save(oImageOut->value() + ".points");
        painter.draw(graph, points);
        painter.save(oImageOut->value() + ".graph");
    }
#endif

    if (sVerify->is_set()) {
        tYaoGraph exp(points.size(), oK->value());

        if (sVerify->count() == 1) {
            switch (oKern->value()) {
                case 'i':
                    exp = runAlg<GridYao<InexactKernel>>(oK->value(), points, bounds, oCellOcc->value());
                    break;
                case 'p':
                    exp = runAlg<GridYao<CGALKernel<ExactPredicatesInexactConstructions>>>(oK->value(), points, bounds, oCellOcc->value());
                    break;
                case 'c':
                    exp = runAlg<GridYao<CGALKernel<ExactPredicatesExactConstructions>>>(oK->value(), points, bounds, oCellOcc->value());
                    break;
            }
        } else if (sVerify->count() == 2) {
            switch (oKern->value()) {
                case 'i':
                    exp = runAlg<NaiveYao<InexactKernel>>(oK->value(), points, bounds);
                    break;
                case 'p':
                    exp = runAlg<NaiveYao<CGALKernel<ExactPredicatesInexactConstructions>>>(oK->value(), points, bounds);
                    break;
                case 'c':
                    exp = runAlg<NaiveYao<CGALKernel<ExactPredicatesExactConstructions>>>(oK->value(), points, bounds);
                    break;
            }
        }

        auto [valid, invalidVertices] = checkGraph(graph, exp);

#ifdef WITH_CAIRO
        if (!valid) {
            Painter basePainter(bounds, 1000);
            basePainter.draw(points, true);
            basePainter.draw(exp, points);

            for (auto idx : invalidVertices) {

                Painter painter(basePainter);

                for (tDim k = 0; k < oK->value(); ++k) {
                    if (exp[idx].neighbor[k] != graph[idx].neighbor[k]) {

                        if (exp[idx].neighbor[k] != INF_IDX) {
                            painter.setColor(0, 1, 0);
                            painter.draw(points[exp[idx].neighbor[k]], exp[idx].neighbor[k], true);
                        }

                        if (graph[idx].neighbor[k] != INF_IDX) {
                            painter.setColor(1, 0, 0);
                            painter.draw(points[graph[idx].neighbor[k]], graph[idx].neighbor[k], true);
                        }
                    }
                }

                painter.setColor(1, 0, 0);
                painter.draw(idx, graph[idx], points);
                painter.setColor(0, 0, 1);
                painter.drawCones(points[idx], oK->value());

                painter.save("invalidVertices_" + std::to_string(idx));
            }
        }
#endif// ifdef WITH_CAIRO

        return !valid;
    }
}
