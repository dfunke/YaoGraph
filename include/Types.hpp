#pragma once

#include <array>
#include <unordered_set>
#include <vector>
#include <limits>
#include <type_traits>

using tFloat = float;
using tDim = unsigned short;
using tIndex = std::size_t;

constexpr tDim D = 2;
constexpr int SEED = 1986;

using tFloatVector = std::array<tFloat, D>;
using tIndexVector = std::array<tIndex, D>;
using tPoints = std::vector<tFloatVector>;
using tPointSet = std::unordered_set<tIndex>;

struct tBox {
    tFloatVector low;
    tFloatVector high;
};

template<tDim K>
struct tYaoVertex {

    tYaoVertex() {
        for (tDim i = 0; i < neighbor.size(); ++i) {
            neighbor[i] = tIndex(-1);
            distance[i] = std::numeric_limits<tFloat>::max();
        }
    }

    std::array<tIndex, K> neighbor;
    std::array<tFloat, K> distance;
};

template<typename V>
using tYaoGraph = std::vector<V>;
