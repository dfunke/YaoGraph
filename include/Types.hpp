#pragma once

#include <array>
#include <unordered_set>
#include <vector>
#include <limits>

using tCoord = float;
using tDim = unsigned short;
using tIndex = std::size_t;

constexpr tDim D = 2;
constexpr int SEED = 1986;

using tPoint = std::array<tCoord, D>;
using tPoints = std::vector<tPoint>;
using tPointSet = std::unordered_set<tIndex>;

struct tBox {
    tPoint low;
    tPoint high;
};

template<tDim K>
struct tYaoVertex {

    tYaoVertex() {
        for (tDim i = 0; i < neighbor.size(); ++i) {
            neighbor[i] = tIndex(-1);
            distance[i] = std::numeric_limits<tCoord>::max();
        }
    }

    std::array<tIndex, K> neighbor;
    std::array<tCoord, K> distance;
};

template<typename V>
using tYaoGraph = std::vector<V>;
