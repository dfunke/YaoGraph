#pragma once

#include <array>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using tCoord = float;
using tDim = short;
using tIndex = std::size_t;

constexpr tDim D = 2;
constexpr int SEED = 1986;

using tPoint = std::array<tCoord, D>;
using tPoints = std::vector<tPoint>;
using tPointSet = std::unordered_set<tIndex>;

tCoord distance(const tPoint &a, const tPoint &b) {
    tCoord dist = 0;

    for (tDim d = 0; d < D; ++d) {
        dist += (b[d] - a[d]) * (b[d] - a[d]);
    }

    return std::sqrt(dist);
}

struct tBox {
    tPoint low;
    tPoint high;
};