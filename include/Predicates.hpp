#pragma once

#include <cmath>

#include "Types.hpp"

tCoord distance2(const tPoint &a, const tPoint &b) {
    tCoord dist2 = 0;
    for (tDim d = 0; d < a.size(); ++d) {
        dist2 += (b[d] - a[d]) * (b[d] - a[d]);
    }
    return dist2;
}

tCoord distance(const tPoint &a, const tPoint &b) {
    return std::sqrt(distance2(a, b));
}

tCoord atan2P(const tCoord &y, const tCoord &x) {
    auto rad = std::atan2(y, x);
    return rad >= 0 ? rad : 2 * M_PI + rad;
}

tCoord atan2P(const tPoint &a, const tPoint &b) {
    return atan2P(b[1] - a[1], b[0] - a[0]);
}