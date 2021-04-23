//
// Created by Daniel Funke on 18.03.21.
//

#pragma once

#include <queue>

#include "Predicates.hpp"
#include "Types.hpp"

template <tDim K> class SweepLine {

public:
  using tVertex = tYaoVertex<K>;
  using tGraph = tYaoGraph<tVertex>;

  tGraph operator()(const tPoints &points) const {
    tGraph g(points.size());

    for (tDim k = 0; k < K; ++k) {
      sweepline(points, k * (2 * M_PI / K), (k + 1) * (2 * M_PI / K), g);
    }

    return g;
  }

private:
  void sweepline(const tPoints &points, tFloat lTheta, tFloat uTheta,
                 tGraph &graph) const {

    enum pqType { P, IS };

    struct pqItem {
      tFloat key;
      pqType type;
      tFloatVector coords;
    };

    auto pqCmp = [](const pqItem &a, const pqItem &b) { return a.key < b.key; };

    tFloat slAngle = 0.5 * (M_PI + lTheta + uTheta);
    tFloat slDir = M_PI + .5 * (lTheta + uTheta);
    tFloatVector slVec = {std::cos(slDir), std::sin(slDir)};

    auto prj = [&slVec](const tFloatVector &p) {
      auto projection = slVec - ((dot(p, slVec) / dot(p, p)) * p);
      return dot(projection, projection);
    };

    std::priority_queue<pqItem, std::vector<pqItem>, decltype(pqCmp)> pq(pqCmp);

    for (const auto &p : points) {
      pq.emplace(prj(p), pqType::P, p);
    }

    while (!pq.empty()) {

      auto i = pq.top();
      pq.pop();

      switch (i.type) {
      case pqType::P:
        break;
      case pqType::IS:
        break;
      };
    }
  }
};