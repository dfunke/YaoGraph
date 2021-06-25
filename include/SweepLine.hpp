//
// Created by Daniel Funke on 18.03.21.
//

#pragma once

#include <list>
#include <queue>

#include "Predicates.hpp"
#include "Types.hpp"

#ifdef WITH_CAIRO

#include "Painter.hpp"

#endif

struct Ray {
    tFloatVector p;
    tFloat angle;

    tIndex leftRegion;
    tIndex rightRegion;

    // extension of ray with another ray
    std::unique_ptr<Ray> ext;// initialized with nullptr

    Ray(const tFloatVector &p_, const tFloat &angle_)
        : p(p_), angle(angle_), leftRegion(tIndex(-1)), rightRegion(tIndex(-1)) {}

    Ray(const tFloatVector &p_, const tFloat &angle_, const tIndex &lr,
        const tIndex &rr)
        : p(p_), angle(angle_), leftRegion(lr), rightRegion(rr) {}

    Ray(const tFloatVector &p_, const tFloat &angle_, const tIndex &lr,
        const tIndex &rr, const Ray &ext_)
        : p(p_), angle(angle_), leftRegion(lr), rightRegion(rr),
          ext(std::make_unique<Ray>(ext_)) {

        // the starting point of the lower ray must be on the upper ray
        assert(ext->p[Y] == std::tan(angle) * (ext->p[X] - p[X]) + p[Y]);
        assert(leftRegion == ext->leftRegion);
        assert(rightRegion == ext->rightRegion);
    }

    Ray(const Ray &o)
        : p(o.p), angle(o.angle), leftRegion(o.leftRegion),
          rightRegion(o.rightRegion) {
        if (o.ext) {
            ext = std::make_unique<Ray>(*o.ext);
        }
    }

    Ray &operator=(const Ray &o) {
        p = o.p;
        angle = o.angle;
        leftRegion = o.leftRegion;
        rightRegion = o.rightRegion;

        if (o.ext) {
            ext = std::make_unique<Ray>(*o.ext);
        }

        return *this;
    }

    // convenience constructor for ray unions
    Ray(const Ray &upper, const Ray &lower) : Ray(upper.p, upper.angle, upper.leftRegion, upper.rightRegion, lower) {}

    bool leftOf(const tFloatVector &x) const {
        // we only consider main ray, when the starting point of lower ray is
        // swept, this ray will be replaced by it

        return std::signbit((std::cos(angle) - p[X] * (x[Y] - p[Y]) -
                             (std::sin(angle) - p[Y]) * (x[X] - p[X])));
    }

    tFloatVector intersection(const Ray &b, const tBox &bounds) const {
        return intersection(*this, b, bounds);
    }

public:
    static tFloatVector intersection(const Ray &a, const Ray &b, const tBox &bounds) {
        if (!a.ext && !b.ext) {
            return isRR(a, b, bounds);
        } else if (a.ext && b.ext) {
            return isUU(a, b, bounds);
        } else if (a.ext) {
            return isUR(a, b, bounds);
        } else {// b.ext
            return isUR(b, a, bounds);
        }
    }

private:
    static tFloatVector isRR(const Ray &a, const Ray &b, const tBox &bounds) {
        if (std::fmod(std::fmod(a.angle - b.angle, M_PI) + M_PI, M_PI) == 0)
            throw std::domain_error("parallel lines");

        if (std::fmod(std::fmod(a.angle, M_PI) + M_PI, M_PI) == M_PI_2) {
            // vertical line this's x
            return {a.p[X], std::tan(b.angle) * (a.p[X] - b.p[X]) + b.p[Y]};
        } else if (std::fmod(std::fmod(b.angle, M_PI) + M_PI, M_PI) == M_PI_2) {
            // vertical line at b's x
            return {b.p[X], std::tan(a.angle) * (b.p[X] - a.p[X]) + a.p[Y]};
        }

        tFloat m0 = std::tan(a.angle);// Line 0: y = m0 (x - x0) + y0
        tFloat m1 = std::tan(b.angle);// Line 1: y = m1 (x - x1) + y1

        tFloat x = ((m0 * a.p[X] - m1 * b.p[X]) - (a.p[Y] - b.p[Y])) / (m0 - m1);
        tFloat y = m0 * (x - a.p[X]) + a.p[Y];

        if (!bounds.contains({x, y})) {
            throw std::domain_error("outside bounds");
        }

        return {x, y};
    }

    static tFloatVector isUR(const Ray &ua, const Ray &b, const tBox &bounds) {
        assert(ua.ext);

        try {
            // check for intersection of main ray of ua and b
            auto is = isRR(ua, b, bounds);
            // check whether IS is before starting point of extension ray of ur
            if (distance2(ua.p, is) < distance2(ua.p, ua.ext->p)) {
                return is;
            }
        } catch (std::domain_error &e) {
        }

        // upper ray of ur does not intersect r OR intersection is below starting point of extension ray
        return isRR(*ua.ext, b, bounds);
    }

    static tFloatVector isUU(const Ray &ua, const Ray &ub, const tBox &bounds) {
        assert(ua.ext && ub.ext);

        try {
            // check for intersection of main ray of ua and ub
            auto is = isRR(ua, ub, bounds);
            // check whether IS is before starting point of lowerRay for both rays
            if (distance2(ua.p, is) < distance2(ua.p, ua.ext->p) &&
                distance2(ub.p, is) < distance2(ub.p, ub.ext->p)) {
                return is;
            }
        } catch (std::domain_error &e) {
        }

        try {
            // check for intersection of main ray of ua and lower ray of ub
            auto is = isRR(ua, *ub.ext, bounds);
            // check whether IS is before starting point of lowerRay of ua
            if (distance2(ua.p, is) < distance2(ua.p, ua.ext->p)) {
                return is;
            }
        } catch (std::domain_error &e) {
        }

        try {
            // check for intersection of lower ray of ua and main ray of ub
            auto is = isRR(*ua.ext, ub, bounds);
            // check whether IS is before starting point of ub's lowerRay
            if (distance2(ub.p, is) < distance2(ub.p, ub.ext->p)) {
                return is;
            }
        } catch (std::domain_error &e) {
        }

        // check for intersection of lower rays of ua and ub
        return isRR(*ua.ext, *ub.ext, bounds);
    }
};

struct SweeplineDS {

    using tRays = std::list<Ray>;
    using tConstRayHandle = typename tRays::const_iterator;
    using tRayHandle = typename tRays::iterator;

    tFloat slDirection;
    tFloatVector slDirVector;
    tRays slRays;

    SweeplineDS(const tFloat slDir_)
        : slDirection(slDir_),
          slDirVector({std::cos(slDirection), std::sin(slDirection)}) {}

    auto begin() const { return tConstRayHandle(slRays.begin()); }

    auto end() const { return tConstRayHandle(slRays.end()); }

    auto begin() { return tRayHandle(slRays.begin()); }

    auto end() { return tRayHandle(slRays.end()); }

    auto insert(const tRayHandle &pos, const Ray &ray) {
        return slRays.insert(pos, ray);
    }

    auto erase(const tRayHandle &pos) {
        return slRays.erase(pos);
    }

    tFloat prj(const tFloatVector &p) {
        return -(p[X] * slDirVector[X] + p[Y] * slDirVector[Y]) /
               dot(slDirVector, slDirVector);
    }

    // return iterator to ray immediatly right of point p
    auto find(const tFloatVector &p) {// TODO: make const?
        auto it = begin();
        while (it != end() && it->leftOf(p))
            ++it;

        return it;
    }

#ifdef WITH_CAIRO

    void draw(const tFloatVector &pos, Painter &painter) const {

        // draw sweepline
        painter.drawLine(pos, {pos[X] + std::cos(slDirection + tFloat(M_PI_2)),
                               pos[Y] + std::sin(slDirection + tFloat(M_PI_2))});
        painter.drawLine(pos, {pos[X] + std::cos(slDirection - tFloat(M_PI_2)),
                               pos[Y] + std::sin(slDirection - tFloat(M_PI_2))});

        // draw direction
        painter.setDash();
        painter.drawLine(pos, {pos[X] + slDirVector[X], pos[Y] + slDirVector[Y]});
        painter.unsetDash();

        // draw rays
        for (const auto &r : slRays) {
            if (!r.ext) {
                painter.drawLine(
                        r.p, {r.p[X] + std::cos(r.angle), r.p[Y] + std::sin(r.angle)});
            } else {
                painter.drawLine(r.p, r.ext->p);
                painter.drawLine(r.ext->p, {r.ext->p[X] + std::cos(r.ext->angle), r.ext->p[Y] + std::sin(r.ext->angle)});
            }
        }
    }

#endif
};

struct Event {

    enum Type { Input,
                Intersection,
                Deletion };

    using tRayHandle = typename SweeplineDS::tRayHandle;

    Event(const tFloatVector &pos_)
        : type(Type::Deletion), pos(pos_), idx(tIndex(-1)) {}

    Event(const tFloatVector &pos_, const tIndex &idx_)
        : type(Type::Input), pos(pos_), idx(idx_) {}

    Event(const tFloatVector &pos_, const tRayHandle &left)
        : type(Type::Deletion), pos(pos_), leftRay(left) {}

    Event(const tFloatVector &pos_, const tRayHandle &left,
          const tRayHandle &right)
        : type(Type::Intersection), pos(pos_), idx(tIndex(-1)), leftRay(left),
          rightRay(right) {}

    Type type;
    tFloatVector pos;

    tIndex idx;

    tRayHandle leftRay;
    tRayHandle rightRay;
};

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
    void sweepline(const tPoints &points, tDim k, tGraph &graph,
                   const tBox &bounds) const {

        using pqItem = std::pair<tFloat, Event>;
        auto pqCmp = [](const pqItem &a, const pqItem &b) {
            return a.first < b.first;
        };

        tFloat lTheta = k * (2 * M_PI / K);
        tFloat uTheta = (k + 1) * (2 * M_PI / K);

        SweeplineDS sl(M_PI + .5 * (lTheta + uTheta));
        std::priority_queue<pqItem, std::vector<pqItem>, decltype(pqCmp)> pq(pqCmp);
        std::priority_queue<pqItem, std::vector<pqItem>, decltype(pqCmp)> delPQ(
                pqCmp);

        for (tIndex i = 0; i < points.size(); ++i) {
            // std::cout << i << ": " << points[i][X] << ", " << points[i][Y]  << "
            // prj: " << sl.prj(points[i]) << std::endl;
            pq.push({sl.prj(points[i]), Event({points[i], i})});
        }

#ifdef WITH_CAIRO
        Painter basePainter(bounds, 1000);
#endif
        tIndex idx = 0;
        while (!pq.empty()) {

            auto p = pq.top().second;
            std::cout << idx << ": (" << p.pos[X] << ", " << p.pos[Y]
                      << ") key: " << pq.top().first << std::endl;

#ifdef WITH_CAIRO
            basePainter.draw(p.pos, idx);
            sl.draw(p.pos, basePainter);
#endif

            switch (p.type) {
                case Event::Type::Input: {

                    auto itBr = sl.find(p.pos);// right ray
                    auto itBl =
                            (itBr == sl.begin() ? sl.end() : std::prev(itBr));// left ray

                    // add graph edge
                    if (itBr != sl.end() && itBr->leftRegion != tIndex(-1)) {
                        assert(itBl->rightRegion == itBr->leftRegion);

                        graph[p.idx].neighbor[k] = itBr->leftRegion;
                        graph[p.idx].distance[k] = distance2(p.pos, points[itBr->leftRegion]);
                    }

                    // check if Bl and Br intersect
                    if (itBl != sl.end() && itBr != sl.end()) {
                        try {
                            auto is = itBl->intersection(*itBr, bounds);
                            delPQ.push({sl.prj(is), Event({is})});
                        } catch (std::domain_error &e) {
                        }
                    }

                    // create new rays and insert them into SL
                    auto itBln = sl.insert(
                            itBr,
                            Ray({p.pos, lTheta + tFloat(M_PI),
                                 itBl != sl.end() ? itBl->rightRegion : tIndex(-1), p.idx}));
                    auto itBrn = sl.insert(
                            itBr, Ray({p.pos, uTheta + tFloat(M_PI), p.idx,
                                       itBr != sl.end() ? itBr->leftRegion : tIndex(-1)}));

                    // insert intersection points into PQ
                    if (itBl != sl.end()) {
                        try {
                            auto is = itBl->intersection(*itBln, bounds);
                            pq.push({sl.prj(is), Event({is, itBl, itBln})});
                        } catch (std::domain_error &e) {
                        }
                    }

                    if (itBr != sl.end()) {
                        try {
                            auto is = itBr->intersection(*itBrn, bounds);
                            pq.push({sl.prj(is), Event({is, itBrn, itBr})});
                        } catch (std::domain_error &e) {
                        }
                    }

                    break;
                }

                case Event::Type::Intersection: {

                    // check whehter IS is deleted
                    if (!delPQ.empty() && delPQ.top().second.pos == p.pos) {
                        delPQ.pop();
                        break;
                    }

                    auto itBl = p.leftRay;
                    auto itBr = p.rightRay;

                    assert(std::next(itBl) == itBr);
                    assert(itBl == std::prev(itBr));

                    // delete intersection points from PQ
                    if (itBl != sl.begin()) {
                        try {
                            auto itBll = std::prev(itBl);
                            auto is = itBl->intersection(*itBll, bounds);
                            delPQ.push({sl.prj(is), Event({is})});
                        } catch (std::domain_error &e) {
                        }
                    }

                    if (itBr != sl.end()) {
                        try {
                            auto itBrr = std::next(itBr);
                            if (itBrr != sl.end()) {
                                auto is = itBr->intersection(*itBrr, bounds);
                                delPQ.push({sl.prj(is), Event({is})});
                            }
                        } catch (std::domain_error &e) {
                        }
                    }

                    Ray rL({p.pos, lTheta + tFloat(M_PI), itBl->leftRegion,
                            itBr->rightRegion});
                    Ray rR({p.pos, uTheta + tFloat(M_PI), itBl->leftRegion,
                            itBr->rightRegion});

                    // bisector line between A and B
                    auto pL = points[itBl->leftRegion];
                    auto pR = points[itBr->rightRegion];
                    auto pMid = 0.5 * (pL + pR);
                    tFloat aBs = std::atan2(pR[Y] - pL[Y], pR[X] - pL[X]) +
                                 M_PI_2;// TODO check angle orientation
                    Ray Bs({pMid, aBs, itBl->leftRegion, itBr->rightRegion});

                    // check for intersections of Bln, Blr and BS
                    bool BsL = false;
                    tFloatVector pBsL;
                    bool BsR = false;
                    tFloatVector pBsR;
                    try {
                        pBsL = Bs.intersection(rL, bounds);
                        BsL = true;
                    } catch (std::domain_error &e) {
                    }

                    try {
                        pBsR = Bs.intersection(rR, bounds);
                        BsR = true;
                    } catch (std::domain_error &e) {
                    }

                    // remove old Rays from list
                    sl.erase(itBl);
                    auto itInsertPos = sl.erase(itBr);
                    auto itBn = sl.end();

                    if (!BsL && !BsR) {
                        // bisector intersects no ray from P
                        if (sl.prj(pR) < sl.prj(pL)) {
                            // right point is further from v -> right ray becomes border
                            itBn = sl.insert(itInsertPos, rR);
                        } else {
                            // left point is further from v -> left ray becomes border
                            itBn = sl.insert(itInsertPos, rL);
                        }
                    } else {
                        if (BsL && BsR) {
                            // bisector intersects both rays - > must be in same point v
                            assert(pBsL == pBsR);
                            assert(pBsL == p.pos);
                            assert(pBsR == p.pos);
                            // boundary originates at v with bisector angle
                            itBn = sl.insert(itInsertPos, Ray({p.pos, aBs, itBl->leftRegion,
                                                               itBr->rightRegion}));
                        } else {
                            if (BsL) {
                                // bisector intersects left ray
                                itBn = sl.insert(itInsertPos,
                                                 Ray({rL, Ray({pBsL, aBs, itBl->leftRegion,
                                                               itBr->rightRegion})}));
                                pq.push({sl.prj(pBsL), Event({pBsL, itBn})});
                            } else {
                                // bisector intersects right ray
                                itBn = sl.insert(itInsertPos,
                                                 Ray({rR, Ray({pBsR, aBs, itBl->leftRegion,
                                                               itBr->rightRegion})}));
                                pq.push({sl.prj(pBsR), Event({pBsR, itBn})});
                            }
                        }
                    }
                    assert(itBn != sl.end());// some boundary was inserted into SL

                    // insert intersection points into PQ
                    if (itBn != sl.begin()) {
                        try {
                            auto itL = std::prev(itBn);
                            auto is = itBn->intersection(*itL, bounds);
                            pq.push({sl.prj(is), Event({is, itL, itBn})});
                        } catch (std::domain_error &e) {
                        }
                    }

                    auto itR = std::next(itBn);
                    if (itR != sl.end()) {
                        try {
                            auto is = itBn->intersection(*itR, bounds);
                            pq.push({sl.prj(is), Event({is, itBn, itR})});
                        } catch (std::domain_error &e) {
                        }
                    }

                    break;
                }

                case Event::Type::Deletion: {

                    auto itB = p.leftRay;// we store the ray to be deleted as left ray
                    *itB = *(itB->ext);
                    // we replace the RayUnion with its lower ray, so all intersections pointers should still be valid
                    // all intersections processed after this point will be with lower ray

                    break;
                }
            }

            pq.pop();
            ++idx;
        }

#ifdef WITH_CAIRO
        basePainter.save("01_points_" + std::to_string(k));
#endif
    }
};