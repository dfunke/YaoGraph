#pragma once

#include <array>
#include <limits>
#include <type_traits>
#include <unordered_set>
#include <vector>

using tIFloat = double;// Inexact Float
using tDim = unsigned short;
using tIndex = std::size_t;

constexpr tDim X = 0;
constexpr tDim Y = 1;
constexpr tDim Z = 2;
constexpr tIndex INF_IDX = tIndex(-1);

constexpr tDim D = 2;
constexpr int SEED = 1986;

using tIFloatVector = std::array<tIFloat, D>;
using tIndexVector = std::array<tIndex, D>;
using tPoints = std::vector<tIFloatVector>;
using tIndexSet = std::unordered_set<tIndex>;

enum tOrientedSide {
    RIGHT = -1,
    LINE,
    LEFT
};

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
    tIFloatVector low;
    tIFloatVector high;

    bool contains(const tIFloatVector &p) const {
        ASSERT(p.size() == low.size());
        for (tDim d = 0; d < p.size(); ++d) {
            if (p[d] < low[d] || high[d] < p[d]) {
                return false;
            }
        }

        return true;
    }
};

struct tYaoVertex {

    explicit tYaoVertex(const tDim &k) : neighbor(k, INF_IDX), distance(k, std::numeric_limits<tIFloat>::max()) {}

    std::vector<tIndex> neighbor;
    std::vector<tIFloat> distance;
};

struct tYaoGraph : public std::vector<tYaoVertex> {
    tYaoGraph(const tIndex &n, const tDim &k) : std::vector<tYaoVertex>(n, tYaoVertex(k)) {}
};
