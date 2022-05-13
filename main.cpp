#include <fstream>
#include <iostream>

#include "Generators.hpp"
#include "Timer.hpp"
#include "Types.hpp"

#include "GridYao.hpp"
#include "NaiveYao.hpp"
#include "SweepLine.hpp"

#ifdef WITH_CAIRO
#include "Painter.hpp"
#endif

#ifdef WITH_CGAL
#include "CGAL/CGAL_Delaunay.hpp"
#include "CGAL/CGAL_Theta.hpp"
#include "CGAL/CGAL_Yao.hpp"
#include "utils/CGALKernel.hpp"
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

template<typename Algorithm, typename Distribution>
void benchmark() {
    std::ofstream file("benchmark_" + Algorithm::name() + ".csv", std::ofstream::out | std::ofstream::app);

    // header
    file << "dist n seed rep t" << std::endl;
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

int main() {

#ifdef WITH_CGAL
    benchmark<CGAL_Yao2D<Cones, ExactPredicatesInexactConstructions>, Gaussian>();
    benchmark<CGAL_Yao2D<Cones, ExactPredicatesExactConstructions>, Gaussian>();
#endif

    benchmark<NaiveYao<Cones, InexactKernel>, Gaussian>();
#ifdef WITH_CGAL
    benchmark<NaiveYao<Cones, CGALKernel<ExactPredicatesInexactConstructions>>, Gaussian>();
    benchmark<NaiveYao<Cones, CGALKernel<ExactPredicatesExactConstructions>>, Gaussian>();
#endif

    benchmark<GridYao<Cones, InexactKernel, cellOcc>, Gaussian>();
#ifdef WITH_CGAL
    benchmark<GridYao<Cones, CGALKernel<ExactPredicatesInexactConstructions>, cellOcc>, Gaussian>();
    benchmark<GridYao<Cones, CGALKernel<ExactPredicatesExactConstructions>, cellOcc>, Gaussian>();
#endif

    benchmark<SweepLine<Cones, InexactKernel>, Gaussian>();
#ifdef WITH_CGAL
    benchmark<SweepLine<Cones, CGALKernel<ExactPredicatesInexactConstructions>>, Gaussian>();
    benchmark<SweepLine<Cones, CGALKernel<ExactPredicatesExactConstructions>>, Gaussian>();
#endif

    return 0;
}
