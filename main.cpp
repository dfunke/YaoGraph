#include <fstream>
#include <iostream>

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

template<typename Algorithm, typename Distribution>
void benchmarkImpl() {
    std::ofstream file("benchmark_" + Algorithm::name() + ".csv", std::ios::out | std::ios::app);

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

int main() {

#ifdef WITH_CGAL
    benchmark<CGAL_Yao2D<Cones, ExactPredicatesInexactConstructions>, Gaussian, Uniform, Grid>();
    benchmark<CGAL_Yao2D<Cones, ExactPredicatesExactConstructions>, Gaussian, Uniform, Grid>();
#endif

    benchmark<NaiveYao<Cones, InexactKernel>, Gaussian, Uniform, Grid>();
#ifdef WITH_CGAL
    benchmark<NaiveYao<Cones, CGALKernel<ExactPredicatesInexactConstructions>>, Gaussian, Uniform, Grid>();
    benchmark<NaiveYao<Cones, CGALKernel<ExactPredicatesExactConstructions>>, Gaussian, Uniform, Grid>();
#endif

    benchmark<GridYao<Cones, InexactKernel, cellOcc>, Gaussian, Uniform, Grid>();
#ifdef WITH_CGAL
    benchmark<GridYao<Cones, CGALKernel<ExactPredicatesInexactConstructions>, cellOcc>, Gaussian, Uniform, Grid>();
    benchmark<GridYao<Cones, CGALKernel<ExactPredicatesExactConstructions>, cellOcc>, Gaussian, Uniform, Grid>();
#endif

    benchmark<SweepLine<Cones, InexactKernel>, Gaussian, Uniform, Grid>();
#ifdef WITH_CGAL
    benchmark<SweepLine<Cones, CGALKernel<ExactPredicatesInexactConstructions>>, Gaussian, Uniform, Grid>();
    benchmark<SweepLine<Cones, CGALKernel<ExactPredicatesExactConstructions>>, Gaussian, Uniform, Grid>();
#endif

    return 0;
}
