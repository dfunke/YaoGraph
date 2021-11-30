#pragma once

#include <array>
#include <cassert>
#include <limits>
#include <type_traits>
#include <unordered_set>
#include <vector>

using tFloat = double;
using tDim = unsigned short;
using tIndex = std::size_t;

constexpr tDim X = 0;
constexpr tDim Y = 1;

constexpr tDim D = 2;
constexpr int SEED = 1986;

using tFloatVector = std::array<tFloat, D>;
using tIndexVector = std::array<tIndex, D>;
using tPoints = std::vector<tFloatVector>;
using tIndexSet = std::unordered_set<tIndex>;

template<typename F>
struct MaxError;

template<>
struct MaxError<float> {
    static constexpr float value = 1e-10;
};

template<>
struct MaxError<double> {
    static constexpr double value = 1e-22;
};

struct tBox {
    tFloatVector low;
    tFloatVector high;

    bool contains(const tFloatVector &p) const {
        assert(p.size() == low.size());
        for (tDim d = 0; d < p.size(); ++d) {
            if (p[d] < low[d] || high[d] < p[d]) {
                return false;
            }
        }

        return true;
    }
};

template<tDim C>
struct tYaoVertex {

    static const tDim K = C;

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
