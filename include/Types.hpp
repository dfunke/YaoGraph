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

struct tBox {
    tPoint low;
    tPoint high;
};
