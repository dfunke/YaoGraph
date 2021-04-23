#include <iostream>

#include "Types.hpp"
#include "Generators.hpp"
#include "Timer.hpp"

#include "NaiveYao.hpp"
#include "GridYao.hpp"
#include "SweepLine.hpp"

#ifdef WITH_CAIRO
#include "Painter.hpp"
#endif

#ifdef WITH_CGAL
#include "CGAL/CGAL_Delaunay.hpp"
#include "CGAL/CGAL_Yao.hpp"
#include "CGAL/CGAL_Theta.hpp"
#endif

int main() {

    constexpr tBox BOUNDS{{0, 0},
                          {1, 1}};
    constexpr tIndex minN = 1e2;
    constexpr tIndex maxN = 1e2;
    constexpr tDim Cones = 6;
    constexpr tIndex cellOcc = 1e3;

    std::cout << "n naive grid";
#ifdef WITH_CGAL
    std::cout << " cgal_del cgal_yao cgal_theta";
#endif
    std::cout << std::endl;

    for (tIndex nPoints = minN; nPoints <= maxN; nPoints += 3 * pow(10, floor(log10(nPoints)))) {
        Uniform uni;
        auto points = uni.generate(nPoints, BOUNDS);

        auto naive = Timer<NaiveYao<Cones>>::time(points);
        auto grid = Timer<GridYao<Cones>>::time(points, BOUNDS, cellOcc);
        auto sl = Timer<SweepLine<Cones>>::time(points);

        checkGraph(std::get<1>(grid), std::get<1>(naive));

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

        tIndex i = 71;
        auto exp = std::get<1>(naive)[i];
        auto is = std::get<1>(grid)[i];
        Painter detailPainter(basePainter);
        detailPainter.drawCones(points[i], Cones);

        detailPainter.setColor(.5, .5, .5, .5);
        detailPainter.drawGrid(25 * 25);

        detailPainter.setColor(1, 0, 0);
        detailPainter.draw(i, exp, points);

        detailPainter.setColor(0, 0, 0);
        detailPainter.draw(i, is, points);
        detailPainter.save("04_detail");

        std::cout << "is:  " << is << std::endl;
        std::cout << "exp: " << exp << std::endl;
        for (tDim k = 0; k < Cones; ++k) {
            auto nis = is.neighbor[k];
            auto nexp = exp.neighbor[k];
            std::cout << "\t" << i << "-" << nis;
            if (nis != tIndex(-1)) {
                std::cout << ": " << atan2P(points[i], points[nis]) << " rad, sector "
                          << std::floor(atan2P(points[i], points[nis]) / (2 * M_PI / Cones))
                          << " dist " << is.distance[k];
            }
            if (nexp != nis && nexp != tIndex(-1)) {
                std::cout << " EXP : " << atan2P(points[i], points[nexp]) << " rad, sector "
                          << std::floor(atan2P(points[i], points[nexp]) / (2 * M_PI / Cones))
                          << " dist " << exp.distance[k];
            }
            std::cout << std::endl;
        }
#endif

#ifdef WITH_CGAL
        auto del = Timer<CGAL_Delaunay2D>::time(points);
        auto yao = Timer<CGAL_Yao2D<Cones>>::time(points);
        auto theta = Timer<CGAL_Theta2D<Cones>>::time(points);
#endif

        std::cout << nPoints
                  << " " << std::get<0>(naive)
                  << " " << std::get<0>(grid)
                  #ifdef WITH_CGAL
                  << " " << std::get<0>(del)
                  << " " << std::get<0>(yao)
                  << " " << std::get<0>(theta)
                  #endif
                  << std::endl;

    }

    return 0;
}
