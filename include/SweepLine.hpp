//
// Created by Daniel Funke on 18.03.21.
//

#pragma once

#include <queue>
#include <list>

#include "Predicates.hpp"
#include "Types.hpp"

#ifdef WITH_CAIRO

#include "Painter.hpp"

#endif

template<tDim K>
class SweepLine {

public:
    using tVertex = tYaoVertex<K>;
    using tGraph = tYaoGraph<tVertex>;

    tGraph operator()(const tPoints &points, const tBox &bounds) const {
        tGraph g(points.size());

        for (tDim k = 0; k < K; ++k) {
            sweepline(points, k, g, bounds);
        }

        return g;
    }

private:

    struct Ray {
        tFloatVector p;
        tFloat angle;

        tIndex leftRegion;
        tIndex rightRegion;

        Ray(const tFloatVector &p_, const tFloat &angle_) : p(p_),
                                                            angle(angle_),
                                                            leftRegion(tIndex(-1)),
                                                            rightRegion(tIndex(-1)) {}

        Ray(const tFloatVector &p_, const tFloat &angle_, const tIndex &lr, const tIndex &rr) : p(p_),
                                                                                                angle(angle_),
                                                                                                leftRegion(lr),
                                                                                                rightRegion(rr) {}

        bool leftOf(const tFloatVector &x) const {
            return std::signbit(
                    (std::cos(angle) - p[X] * (x[Y] - p[Y]) - (std::sin(angle) - p[Y]) * (x[X] - p[X])));
        }

        tFloatVector intersection(const Ray &b, const tBox &bounds) const {
            if (std::fmod(std::fmod(angle - b.angle, M_PI) + M_PI, M_PI) == 0)
                throw std::domain_error("parallel lines");

            if (std::fmod(std::fmod(angle, M_PI) + M_PI, M_PI) == M_PI_2) {
                // vertical line this's x
                return {p[X], std::tan(b.angle) * (p[X] - b.p[X]) + b.p[Y]};
            } else if (std::fmod(std::fmod(b.angle, M_PI) + M_PI, M_PI) == M_PI_2) {
                // vertical line at b's x
                return {b.p[X], std::tan(angle) * (b.p[X] - p[X]) + p[Y]};
            }

            tFloat m0 = std::tan(angle); // Line 0: y = m0 (x - x0) + y0
            tFloat m1 = std::tan(b.angle); // Line 1: y = m1 (x - x1) + y1

            tFloat x = ((m0 * p[X] - m1 * b.p[X]) - (p[Y] - b.p[Y])) / (m0 - m1);
            tFloat y = m0 * (x - p[X]) + p[Y];

            if (!bounds.contains({x, y})) {
                throw std::domain_error("outside bounds");
            }

            return {x, y};
        }

    };

    struct RayUnion {
        Ray upperRay;
        Ray lowerRay;

        RayUnion(const Ray &upper, const Ray &lower) : upperRay(upper), lowerRay(lower) {
            // the starting point of the lower ray must be on the upper ray
            assert(lowerRay.p[Y] == std::tan(upperRay.angle) * (lowerRay.p[X] - upperRay.p[X]) + upperRay.p[Y]);
        }

        bool leftOf(const tFloatVector &x) const {
            // union is only active as long as SL is above lowerRay, only upperRay needs to be checked
            return upperRay.leftOf(x);
        }

        tFloatVector intersection(const Ray &b, const tBox &bounds) const {
            try {
                auto is = upperRay.intersection(b, bounds);
                // check whehter IS is before starting point of lowerRay
                if (distance2(upperRay.p, is) < distance2(upperRay.p, lowerRay.p)) {
                    return is;
                }
            } catch (std::domain_error &e) {}

            // upper ray does not intersect b OR intersection is below starting point of lower ray
            return lowerRay.intersection(b, bounds);
        }

        tFloatVector intersection(const RayUnion &b, const tBox &bounds) const {
            try {
                auto is = upperRay.intersection(b.upperRay, bounds);
                // check whehter IS is before starting point of lowerRay for both rays
                if (distance2(upperRay.p, is) < distance2(upperRay.p, lowerRay.p)
                    && distance2(b.upperRay.p, is) < distance2(b.upperRay.p, b.lowerRay.p)) {
                    return is;
                }
            } catch (std::domain_error &e) {}

            try {
                auto is = upperRay.intersection(b.lowerRay, bounds);
                // check whehter IS is before starting point of lowerRay for our ray
                if (distance2(upperRay.p, is) < distance2(upperRay.p, lowerRay.p)) {
                    return is;
                }
            } catch (std::domain_error &e) {}

            try {
                auto is = lowerRay.intersection(b.upperRay, bounds);
                // check whehter IS is before starting point of b's lowerRay
                if (distance2(b.upperRay.p, is) < distance2(b.upperRay.p, b.lowerRay.p)) {
                    return is;
                }
            } catch (std::domain_error &e) {}

            // upper ray does not intersect b OR intersection is below starting point of lower ray
            return lowerRay.intersection(b.lowerRay, bounds);
        }
    };

    struct Event {

        enum Type {
            Input, Intersection, Deletion
        };

        using tHandle = typename std::list<Ray>::const_iterator;

        Event(const tFloatVector &pos_) : type(Type::Deletion), pos(pos_), idx(tIndex(-1)) {}

        Event(const tFloatVector &pos_, const tIndex &idx_) : type(Type::Input), pos(pos_), idx(idx_) {}

        Event(const tFloatVector &pos_, const tHandle &left, const tHandle &right) : type(Type::Intersection),
                                                                                     pos(pos_),
                                                                                     idx(tIndex(-1)),
                                                                                     leftRay(left), rightRay(right) {}

        Type type;
        tFloatVector pos;

        tIndex idx;

        tHandle leftRay;
        tHandle rightRay;

    };

    struct Sweepline {

        tFloat slDirection;
        tFloatVector slDirVector;

        Sweepline(const tFloat slDir_) : slDirection(slDir_),
                                         slDirVector({std::cos(slDirection), std::sin(slDirection)}) {}

        tFloat prj(const tFloatVector &p) {
            return -(p[X] * slDirVector[X] + p[Y] * slDirVector[Y]) / dot(slDirVector, slDirVector);
        }

        std::list<Ray> rays;

        // return iterator to ray immediatly right of point p
        auto find(const tFloatVector &p) const {
            auto it = rays.begin();
            while (it != rays.end() && it->leftOf(p))
                ++it;

            return it;
        }

#ifdef WITH_CAIRO

        void draw(const tFloatVector &pos, Painter &painter) const {

            //draw sweepline
            painter.drawLine(pos, {pos[X] + std::cos(slDirection + tFloat(M_PI_2)),
                                   pos[Y] + std::sin(slDirection + tFloat(M_PI_2))});
            painter.drawLine(pos, {pos[X] + std::cos(slDirection - tFloat(M_PI_2)),
                                   pos[Y] + std::sin(slDirection - tFloat(M_PI_2))});

            //draw direction
            painter.setDash();
            painter.drawLine(pos, {pos[X] + slDirVector[X], pos[Y] + slDirVector[Y]});
            painter.unsetDash();

            //draw rays
            for (const auto &r : rays) {
                painter.drawLine(r.p, {r.p[X] + std::cos(r.angle), r.p[Y] + std::sin(r.angle)});
            }

        }

#endif

    };

    void sweepline(const tPoints &points, tDim k, tGraph &graph, const tBox &bounds) const {

        using pqItem = std::pair<tFloat, Event>;
        auto pqCmp = [](const pqItem &a, const pqItem &b) { return a.first < b.first; };

        tFloat lTheta = k * (2 * M_PI / K);
        tFloat uTheta = (k + 1) * (2 * M_PI / K);

        Sweepline sl(M_PI + .5 * (lTheta + uTheta));
        std::priority_queue<pqItem, std::vector<pqItem>, decltype(pqCmp)> pq(pqCmp);
        std::priority_queue<pqItem, std::vector<pqItem>, decltype(pqCmp)> delPQ(pqCmp);

        for (tIndex i = 0; i < points.size(); ++i) {
            //std::cout << i << ": " << points[i][X] << ", " << points[i][Y]  << " prj: " << sl.prj(points[i]) << std::endl;
            pq.push({sl.prj(points[i]), Event({points[i], i})});
        }

#ifdef WITH_CAIRO
        Painter basePainter(bounds, 1000);
#endif
        tIndex idx = 0;
        while (!pq.empty()) {

            auto p = pq.top().second;
            std::cout << idx << ": (" << p.pos[X] << ", " << p.pos[Y] << ") key: " << pq.top().first << std::endl;

            basePainter.draw(p.pos, idx);
            sl.draw(p.pos, basePainter);

            switch (p.type) {
                case Event::Type::Input: {

                    auto itBr = sl.find(p.pos); // right ray
                    auto itBl = (itBr == sl.rays.begin() ? sl.rays.end() : std::prev(itBr)); // left ray

                    // add graph edge
                    if (itBr != sl.rays.end() && itBr->leftRegion != tIndex(-1)) {
                        assert(itBl->rightRegion == itBr->leftRegion);

                        graph[p.idx].neighbor[k] = itBr->leftRegion;
                        graph[p.idx].distance[k] = distance2(p.pos, points[itBr->leftRegion]);
                    }

                    // check if Bl and Br intersect
                    if (itBl != sl.rays.end() && itBr != sl.rays.end()) {
                        try {
                            auto is = itBl->intersection(*itBr, bounds);
                            delPQ.push({sl.prj(is), Event({is})});
                        } catch (std::domain_error &e) {}
                    }

                    // create new rays and insert them into SL
                    auto itBln = sl.rays.insert(itBr, Ray({p.pos, lTheta + M_PI,
                                                           itBl != sl.rays.end() ? itBl->rightRegion : tIndex(-1),
                                                           p.idx}));
                    auto itBrn = sl.rays.insert(itBr, Ray({p.pos, uTheta + M_PI, p.idx,
                                                           itBr != sl.rays.end() ? itBr->leftRegion : tIndex(-1)}));

                    // insert intersection points into PQ
                    if (itBl != sl.rays.end()) {
                        try {
                            auto is = itBl->intersection(*itBln, bounds);
                            pq.push({sl.prj(is), Event({is}, itBl, itBln)});
                        } catch (std::domain_error &e) {}
                    }

                    if (itBr != sl.rays.end()) {
                        try {
                            auto is = itBr->intersection(*itBrn, bounds);
                            pq.push({sl.prj(is), Event({is}, itBrn, itBr)});
                        } catch (std::domain_error &e) {}
                    }

                    break;
                }

                case Event::Type::Intersection: {

                    // check whehter IS is deleted
                    if (delPQ.top().second.pos == p.pos) {
                        delPQ.pop();
                        break;
                    }

                    auto itBl = p.leftRay;
                    auto itBr = p.rightRay;

                    assert(std::next(itBl) == itBr);
                    assert(itBl == std::prev(itBr));

                    // delete intersection points from PQ
                    if (itBl != sl.rays.begin()) {
                        try {
                            auto itBll = std::prev(itBl);
                            auto is = itBl->intersection(*itBll, bounds);
                            delPQ.push({sl.prj(is), Event({is})});
                        } catch (std::domain_error &e) {}
                    }

                    if (itBr != sl.rays.end()) {
                        try {
                            auto itBrr = std::next(itBr);
                            if (itBrr != sl.rays.end()) {
                                auto is = itBr->intersection(*itBrr, bounds);
                                delPQ.push({sl.prj(is), Event({is})});
                            }
                        } catch (std::domain_error &e) {}
                    }

                    Ray rL({p.pos, lTheta + M_PI, itBl->leftRegion, itBr->rightRegion});
                    Ray rR({p.pos, uTheta + M_PI, itBl->leftRegion, itBr->rightRegion});

                    // bisector line between A and B
                    auto pL = points[itBl->leftRegion];
                    auto pR = points[itBr->rightRegion];
                    auto pMid = 0.5 * (pL + pR);
                    tFloat aBs = std::atan2(pR[Y] - pL[Y], pR[X] - pL[X]) + M_PI_2; //TODO check angle orientation
                    Ray Bs({pMid, aBs, itBl->leftRegion, itBr->rightRegion});

                    // check for intersections of Bln, Blr and BS
                    bool BsL = false;
                    tFloatVector pBsL;
                    bool BsR = false;
                    tFloatVector pBsR;
                    try {
                        pBsL = Bs.intersection(rL, bounds);
                        BsL = true;
                    } catch (std::domain_error &e) {}

                    try {
                        pBsR = Bs.intersection(rR, bounds);
                        BsR = true;
                    } catch (std::domain_error &e) {}

                    //remove old Rays from list
                    sl.rays.erase(itBl);
                    auto itInsertPos = sl.rays.erase(itBr);
                    auto itBn = sl.rays.end();

                    if (!BsL && !BsR) {
                        //bisector intersects no ray from P
                        if (sl.prj(pR) < sl.prj(pL)) {
                            //right point is further from v -> right ray becomes border
                            itBn = sl.rays.insert(itInsertPos, rR);
                        } else {
                            //left point is further from v -> left ray becomes border
                            itBn = sl.rays.insert(itInsertPos, rL);
                        }
                    } else {
                        if (BsL && BsR) {
                            //bisector intersects both rays - > must be in same point v
                            assert(pBsL == pBsR);
                            assert(pBsL == p.pos);
                            assert(pBsR == p.pos);
                            // boundary originates at v with bisector angle
                            itBn = sl.rays.insert(itInsertPos, Ray({p.pos, aBs, itBl->leftRegion, itBr->rightRegion}));
                        } else {
                            if (BsL) {
                                //bisector intersects left ray
                                RayUnion ru(rL, Ray({pBsL, aBs, itBl->leftRegion, itBr->rightRegion}));
                            } else {
                                //bisector intersects right ray
                                RayUnion ru(rR, Ray({pBsR, aBs, itBl->leftRegion, itBr->rightRegion}));
                            }
                        }
                    }
                    assert(itBn != sl.rays.end());

                    break;
                }

                case Event::Type::Deletion: {
                    break;
                }
            }

            pq.pop();
            ++idx;
        }

        basePainter.save("01_points_" + std::to_string(k));
    }

};