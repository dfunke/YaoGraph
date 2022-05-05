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
#endif

int main() {

    constexpr tBox BOUNDS{{0, 0},
                          {1, 1}};
    constexpr tIndex minN = 1e3;
    constexpr tIndex maxN = 1e5;
    constexpr tDim Cones = 6;
    constexpr tIndex cellOcc = 1e3;

    std::cout << "n naive grid sl";
#ifdef WITH_CGAL
    std::cout << " cgal_del cgal_yao cgal_theta";
#endif
    std::cout << std::endl;

    for (tIndex nPoints = minN; nPoints <= maxN; nPoints += 3 * pow(10, floor(log10(nPoints)))) {
        Uniform uni;
        auto points = uni.generate(nPoints, BOUNDS);

        auto naive = Timer<NaiveYao<Cones, CGALKernel>>::time(points);

        auto grid = Timer<GridYao<Cones>>::time(points, BOUNDS, cellOcc);
        checkGraph(std::get<1>(grid), std::get<1>(naive));

        auto sl = Timer<SweepLine<Cones, CGALKernel>>::time(points, BOUNDS);
        checkGraph(std::get<1>(sl), std::get<1>(naive));

#ifdef WITH_CAIRO
        Painter basePainter(BOUNDS, 1000);
        basePainter.draw(points);
        basePainter.save("01_points");

        Painter naivePainter(basePainter);
        naivePainter.draw(std::get<1>(naive), points);
        naivePainter.save("02_naive");

        Painter gridPainter(basePainter);
        gridPainter.draw(std::get<1>(grid), points);
        gridPainter.save("03_grid");

        Painter slPainter(basePainter);
        slPainter.draw(std::get<1>(sl), points);
        slPainter.save("04_sl");
#endif

#ifdef WITH_CGAL
        auto del = Timer<CGAL_Delaunay2D>::time(points);
        auto yao = Timer<CGAL_Yao2D<Cones>>::time(points);
        auto theta = Timer<CGAL_Theta2D<Cones>>::time(points);
#endif

        std::cout << nPoints
                  << " " << std::get<0>(naive)
                  << " " << std::get<0>(grid)
                  << " " << std::get<0>(sl)
#ifdef WITH_CGAL
                  << " " << std::get<0>(del)
                  << " " << std::get<0>(yao)
                  << " " << std::get<0>(theta)
#endif
                  << std::endl;
    }

    return 0;
}
