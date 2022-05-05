#pragma once

#include <cmath>

#include "Types.hpp"

tFloat dot(const tFloatVector &a, const tFloatVector &b) {
    tFloat dot = 0;
    for (tDim d = 0; d < a.size(); ++d) {
        dot += a[d] * b[d];
    }
    return dot;
}

tFloat wrapAngle(const tFloat &a) {
    const auto TwoPi = 2 * M_PI;
    return a - TwoPi * std::floor(a / TwoPi);
}

tFloat atan2P(const tFloat &y, const tFloat &x) {
    auto rad = std::atan2(y, x);
    return rad >= 0 ? rad : 2 * M_PI + rad;
}

tFloat atan2P(const tFloatVector &a, const tFloatVector &b) {
    return atan2P(b[1] - a[1], b[0] - a[0]);
}

template<typename TA, typename TB, std::size_t D,
         typename = typename std::enable_if<std::is_arithmetic<TA>::value, TA>::type,
         typename = typename std::enable_if<std::is_arithmetic<TB>::value, TB>::type>
auto operator+(const std::array<TA, D> &a, const std::array<TB, D> &b) {
    auto result = a;
    for (tDim d = 0; d < a.size(); ++d) {
        result[d] += b[d];
    }
    return result;
}

template<typename TA, typename TB, std::size_t D,
         typename = typename std::enable_if<std::is_arithmetic<TA>::value, TA>::type,
         typename = typename std::enable_if<std::is_arithmetic<TB>::value, TB>::type>
auto operator*(const std::array<TA, D> &a, const std::array<TB, D> &b) {
    auto result = a;
    for (tDim d = 0; d < a.size(); ++d) {
        result[d] *= b[d];
    }
    return result;
}

template<typename TA, typename TB, std::size_t D,
         typename = typename std::enable_if<std::is_arithmetic<TA>::value, TA>::type,
         typename = typename std::enable_if<std::is_arithmetic<TB>::value, TB>::type>
auto operator-(const std::array<TA, D> &a, const std::array<TB, D> &b) {
    auto result = a;
    for (tDim d = 0; d < a.size(); ++d) {
        result[d] -= b[d];
    }
    return result;
}

template<typename TA, typename TB, std::size_t D,
         typename = typename std::enable_if<std::is_arithmetic<TA>::value, TA>::type,
         typename = typename std::enable_if<std::is_arithmetic<TB>::value, TB>::type>
auto operator*(const TA &a, const std::array<TB, D> &b) {
    auto result = b;
    for (tDim d = 0; d < b.size(); ++d) {
        result[d] *= a;
    }
    return result;
}

template<tDim K>
bool operator==(const tYaoVertex<K> &a, const tYaoVertex<K> &b) {
    for (tDim d = 0; d < K; ++d) {
        if (a.neighbor[d] != b.neighbor[d]) {
            return false;
        }
    }

    return true;
}

template<tDim K>
std::ostream &operator<<(std::ostream &os, const tYaoVertex<K> &v) {

    char sep = 0;

    for (tDim d = 0; d < K; ++d) {
        os << sep << (v.neighbor[d] != tIndex(-1) ? std::to_string(v.neighbor[d]) : "I");
        sep = ' ';
    }

    return os;
}

std::ostream &operator<<(std::ostream &os, const tFloatVector &v) {
    //    os << "(" << v[X] << ", " << v[Y] << ")";
    os << v[X] << " " << v[Y];
    return os;
}

template<typename YaoGraph>
std::tuple<bool, tIndexSet> checkGraph(const YaoGraph &is, const YaoGraph &exp) {

    bool valid = true;
    tIndexSet invalidVertices;

    if (is.size() != exp.size()) {
        valid = false;
        std::cerr << "size error: is " << is.size() << " exp " << exp.size() << std::endl;
    }

    for (tIndex i = 0; i < is.size(); ++i) {
        if (is[i] != exp[i]) {
            valid = false;
            invalidVertices.insert(i);

            std::cerr << "vertex error " << i << ":\n\tis:  " << is[i] << "\n\texp: " << exp[i] << std::endl;

            for (tIndex k = 0; k < YaoGraph::value_type::K; ++k) {
                if (is[i].neighbor[k] != exp[i].neighbor[k]) {
                    std::cerr << "\tcone " << k << ": is: " << is[i].neighbor[k] << " (" << is[i].distance[k]
                              << ") exp: " << exp[i].neighbor[k] << " (" << exp[i].distance[k] << ")" << std::endl;
                }
            }
        }
    }

    if (!valid) {
        std::cerr << "graph not valid";
        if (!invalidVertices.empty()) {
            std::cerr << " - " << invalidVertices.size() << "/" << is.size() << " vertices invalid";
        }
        std::cerr << std::endl;
    }

    return {valid, invalidVertices};
}