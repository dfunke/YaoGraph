#include <iostream>

#include "Types.hpp"
#include "Generators.hpp"
#include "Timer.hpp"

#include "NaiveYao.hpp"

#ifdef BUILD_CGAL
#include "CGAL/CGAL_Delaunay.hpp"
#include "CGAL/CGAL_Yao.hpp"
#include "CGAL/CGAL_Theta.hpp"
#endif

int main() {

    constexpr tBox BOUNDS{{0, 0},
                          {1, 1}};
    constexpr tIndex minN = 1e3;
    constexpr tIndex maxN = 1e4;
    constexpr tDim Cones = 6;

    std::cout << "n naive cgal_del cgal_yao cgal_theta" << std::endl;
    for (tIndex nPoints = minN; nPoints <= maxN; nPoints += 3 * pow(10, floor(log10(nPoints)))) {
        Uniform uni;
        auto points = uni.generate(nPoints, BOUNDS);

        auto naive = Timer<NaiveYao<Cones>>::time(points);

#ifdef BUILD_CGAL
        auto del = Timer<CGAL_Delaunay2D>::time(points);
        auto yao = Timer<CGAL_Yao2D<Cones>>::time(points);
        auto theta = Timer<CGAL_Theta2D<Cones>>::time(points);
#else
        int del = 0;
        int yao = 0;
        int theta = 0;
#endif

        std::cout << nPoints << " " << naive << " " << del << " " << yao << " " << theta << std::endl;

    }

    return 0;
}
