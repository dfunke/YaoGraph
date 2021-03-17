#include "Types.hpp"
#include "Generators.hpp"

#include "Delaunay.hpp"
#include "Yao.hpp"
#include "Theta.hpp"

#include "Timer.hpp"

int main() {

    constexpr tBox BOUNDS{{0, 0},
                          {1, 1}};
    constexpr tIndex minN = 1e3;
    constexpr tIndex maxN = 1e6;

    std::cout << "n del yao theta" << std::endl;
    for (tIndex nPoints = minN; nPoints <= maxN; nPoints += 3*pow(10, floor(log10(nPoints)))) {
        Uniform uni;
        auto points = uni.generate(nPoints, BOUNDS);

        auto del = Timer<Delaunay2D>::time(points);
        auto yao = Timer<Yao2D>::time(points);
        auto theta = Timer<Theta2D>::time(points);

        std::cout << nPoints << " " << del << " " << yao << " " << theta << std::endl;
    }

    return 0;
}
