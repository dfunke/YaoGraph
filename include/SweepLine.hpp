//
// Created by Daniel Funke on 18.03.21.
//

#pragma once

#include <list>
#include <queue>
#include <memory>
#include <sstream>

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

    Ray(const tFloatVector &p_, const tFloat &angle_, const tIndex &lr,
        const tIndex &rr)
        : p(p_), angle(angle_), leftRegion(lr), rightRegion(rr) {}

    Ray(const tFloatVector &p_, const tFloat &angle_, const tIndex &lr,
        const tIndex &rr, const Ray &ext_)
        : p(p_), angle(angle_), leftRegion(lr), rightRegion(rr),
          ext(std::make_unique<Ray>(ext_)) {

        // the starting point of the lower ray must be on the upper rayl;
        std::cout << "exp: " << ext->p[Y] << "\nis:  " << std::tan(angle) * (ext->p[X] - p[X]) + p[Y] << std::endl;
        assert(std::abs(ext->p[Y] - (std::tan(angle) * (ext->p[X] - p[X]) + p[Y])) < 0.001);//TODO make more meaningful
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

        return (((p[X] + std::cos(angle)) - p[X]) * (x[Y] - p[Y]) - ((p[Y] + std::sin(angle)) - p[Y]) * (x[X] - p[X])) > 0;
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

#ifdef WITH_CAIRO

    void draw(Painter &painter) const {
        if (!ext) {
            painter.drawLine(p, {p[X] + std::cos(angle), p[Y] + std::sin(angle)});
        } else {
            painter.drawLine(p, ext->p);
            painter.drawLine(ext->p, {ext->p[X] + std::cos(ext->angle), ext->p[Y] + std::sin(ext->angle)});
        }
    }

#endif

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
        } catch (std::domain_error &e) {}

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
        } catch (std::domain_error &e) {}

        try {
            // check for intersection of main ray of ua and lower ray of ub
            auto is = isRR(ua, *ub.ext, bounds);
            // check whether IS is before starting point of lowerRay of ua
            if (distance2(ua.p, is) < distance2(ua.p, ua.ext->p)) {
                return is;
            }
        } catch (std::domain_error &e) {}

        try {
            // check for intersection of lower ray of ua and main ray of ub
            auto is = isRR(*ua.ext, ub, bounds);
            // check whether IS is before starting point of ub's lowerRay
            if (distance2(ub.p, is) < distance2(ub.p, ub.ext->p)) {
                return is;
            }
        } catch (std::domain_error &e) {}

        // check for intersection of lower rays of ua and ub
        return isRR(*ua.ext, *ub.ext, bounds);
    }
};

std::ostream &operator<<(std::ostream &os, const Ray &r) {
    os << "p: " << r.p << " angle: " << r.angle << " boundary: "
       << (r.leftRegion != tIndex(-1) ? std::to_string(r.leftRegion) : "INF") << "/"
       << (r.rightRegion != tIndex(-1) ? std::to_string(r.rightRegion) : "INF");
    if (r.ext) {
        os << " extension: " << *r.ext;
    }

    return os;
}

std::string to_string(const Ray &r) {
    std::stringstream os;
    os << "p: " << r.p << " angle: " << r.angle << " boundary: "
       << (r.leftRegion != tIndex(-1) ? std::to_string(r.leftRegion) : "INF") << "/"
       << (r.rightRegion != tIndex(-1) ? std::to_string(r.rightRegion) : "INF");
    if (r.ext) {
        os << " extension: " << *r.ext;
    }

    return os.str();
}

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
        return (p[X] * slDirVector[X] + p[Y] * slDirVector[Y]) /
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
            r.draw(painter);
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
            return a.first > b.first;
        };

        tFloat lTheta = k * (2 * M_PI / K);
        tFloat uTheta = (k + 1) * (2 * M_PI / K);

        SweeplineDS sl(M_PI + .5 * (lTheta + uTheta));
        std::priority_queue<pqItem, std::vector<pqItem>, decltype(pqCmp)> pq(pqCmp);
        std::priority_queue<pqItem, std::vector<pqItem>, decltype(pqCmp)> delPQ(pqCmp);

        for (tIndex i = 0; i < points.size(); ++i) {
            pq.push({sl.prj(points[i]), Event({points[i], i})});
        }

#ifdef WITH_CAIRO

        constexpr std::array<float, 3> BASE = {0, 0, 0};
        constexpr std::array<float, 3> SL = {0, 0, 1};
        constexpr std::array<float, 3> POINT = {1, 0, 0};
        constexpr std::array<float, 3> IS = {0, 1, 0};
        constexpr std::array<float, 3> BOUND = {0, 1, 1};

        Painter basePainter(bounds, 1000);
        basePainter.setColor(BASE);
        basePainter.draw(points, true);
#endif
        tIndex idx = 0;
        while (!pq.empty()) {

            auto cKey = pq.top().first;
            auto cPoint = pq.top().second;
            pq.pop();

            std::cout << idx << ": " << cPoint.idx << " (" << cPoint.pos[X] << ", " << cPoint.pos[Y]
                      << "): key: " << cKey << std::endl;

#ifdef WITH_CAIRO
            //basePainter.draw(cPoint.pos, idx);

            Painter stepPainter(basePainter);

            stepPainter.setColor(POINT);
            stepPainter.draw(cPoint.pos, idx, false);

            stepPainter.setColor(SL);
            sl.draw(cPoint.pos, stepPainter);
            stepPainter.save("img_k" + std::to_string(k) + "_s" + std::to_string(idx));
#endif

            switch (cPoint.type) {
                case Event::Type::Input: {

                    std::cout << idx << " type: input point" << std::endl;

                    auto itBr = sl.find(cPoint.pos);// right ray
                    auto itBl =
                            (itBr == sl.begin() ? sl.end() : std::prev(itBr));// left ray
                    assert(itBl == sl.end() || itBl->leftOf(cPoint.pos));

                    std::cout << idx << " left ray: " << (itBl != sl.end() ? to_string(*itBl) : "NULL") << std::endl;
                    std::cout << idx << " right ray: " << (itBr != sl.end() ? to_string(*itBr) : "NULL") << std::endl;

                    // add graph edge
                    if (itBr != sl.end() && itBr->leftRegion != tIndex(-1)) {
                        assert(itBl->rightRegion == itBr->leftRegion);

                        graph[cPoint.idx].neighbor[k] = itBr->leftRegion;
                        graph[cPoint.idx].distance[k] = distance2(cPoint.pos, points[itBr->leftRegion]);

                        std::cout << idx << " edge added: (" << cPoint.idx << ", " << itBr->leftRegion
                                  << ") w: " << distance2(cPoint.pos, points[itBr->leftRegion]) << std::endl;
                    }

                    // check if Bl and Br intersect, check only required if they don't originate from same point
                    if (itBl != sl.end() && itBr != sl.end() && itBl->p != itBr->p) {
                        try {
                            auto is = itBl->intersection(*itBr, bounds);
                            delPQ.push({sl.prj(is), Event({is})});
                            std::cout << idx << " deleted intersection point at " << is << " key: " << sl.prj(is) << std::endl;
                        } catch (std::domain_error &e) {}
                    }

                    // create new rays and insert them into SL
                    auto itBln = sl.insert(
                            itBr,
                            Ray({cPoint.pos, lTheta + tFloat(M_PI),
                                 itBl != sl.end() ? itBl->rightRegion : tIndex(-1), cPoint.idx}));
                    auto itBrn = sl.insert(
                            itBr, Ray({cPoint.pos, uTheta + tFloat(M_PI), cPoint.idx,
                                       itBr != sl.end() ? itBr->leftRegion : tIndex(-1)}));

                    std::cout << idx << " left ray " << *itBln << std::endl;
                    std::cout << idx << " right ray " << *itBrn << std::endl;

                    // insert intersection points into PQ
                    if (itBl != sl.end()) {
                        try {
                            auto is = itBl->intersection(*itBln, bounds);
                            if (sl.prj(is) > cKey) {
                                pq.push({sl.prj(is), Event({is, itBl, itBln})});
                                std::cout << idx << " added left intersection point at " << is << " key: " << sl.prj(is) << std::endl;
                            }
                        } catch (std::domain_error &e) {}
                    }

                    if (itBr != sl.end()) {
                        try {
                            auto is = itBrn->intersection(*itBr, bounds);
                            if (sl.prj(is) > cKey) {
                                pq.push({sl.prj(is), Event({is, itBrn, itBr})});
                                std::cout << idx << " added right intersection point at " << is << " key: " << sl.prj(is) << std::endl;
                            }
                        } catch (std::domain_error &e) {}
                    }

                    break;
                }

                case Event::Type::Intersection: {

                    std::cout << idx << " type: intersection point" << std::endl;

                    // check for point that are not processed on delPQ
                    while (!delPQ.empty() && delPQ.top().first < cKey) {
                        delPQ.pop();
                    }

                    if (!delPQ.empty() && delPQ.top().second.pos == cPoint.pos) {
                        delPQ.pop();
                        std::cout << idx << " previously deleted" << std::endl;
                        break;
                    }

                    auto itBl = cPoint.leftRay;
                    auto itBr = cPoint.rightRay;

#ifdef WITH_CAIRO
                    stepPainter.setColor(BOUND);
                    itBl->draw(stepPainter);
                    itBr->draw(stepPainter);
                    stepPainter.save("img_k" + std::to_string(k) + "_s" + std::to_string(idx));
#endif

                    assert(std::next(itBl) == itBr);
                    assert(itBl == std::prev(itBr));

                    std::cout << idx << " left ray: " << (itBl != sl.end() ? to_string(*itBl) : "NULL") << std::endl;
                    std::cout << idx << " right ray: " << (itBr != sl.end() ? to_string(*itBr) : "NULL") << std::endl;

                    // delete intersection points from PQ
                    if (itBl != sl.begin()) {
                        try {
                            auto itBll = std::prev(itBl);

                            if (itBl->p != itBll->p) {
                                auto is = itBll->intersection(*itBl, bounds);
                                delPQ.push({sl.prj(is), Event({is})});
                                std::cout << idx << " deleted left intersection point at " << is << " key: " << sl.prj(is) << std::endl;
                            }
                        } catch (std::domain_error &e) {}
                    }

                    if (itBr != sl.end()) {
                        try {
                            auto itBrr = std::next(itBr);

                            if (itBrr != sl.end() && itBr->p != itBrr->p) {
                                auto is = itBr->intersection(*itBrr, bounds);
                                delPQ.push({sl.prj(is), Event({is})});
                                std::cout << idx << " deleted right intersection point at " << is << " key: " << sl.prj(is) << std::endl;
                            }
                        } catch (std::domain_error &e) {}
                    }

                    Ray rL({cPoint.pos, lTheta + tFloat(M_PI), itBl->leftRegion,
                            itBr->rightRegion});
                    Ray rR({cPoint.pos, uTheta + tFloat(M_PI), itBl->leftRegion,
                            itBr->rightRegion});

                    // bisector line between A and B
                    auto pL = points[itBl->leftRegion];
                    auto pR = points[itBr->rightRegion];
                    auto pMid = 0.5 * (pL + pR);
                    tFloat aBs = std::atan2(pR[Y] - pL[Y], pR[X] - pL[X]);
                    aBs += (((sl.slDirection - aBs) < M_PI) ? 1 : -1) * M_PI_2;// TODO check angle orientation
                    Ray Bs({pMid, aBs, itBl->leftRegion, itBr->rightRegion});

#ifdef WITH_CAIRO
                    stepPainter.setColor(IS);
                    rL.draw(stepPainter);
                    rR.draw(stepPainter);
                    stepPainter.drawLine(pL, pR);
                    stepPainter.drawSquare(pMid);
                    Bs.draw(stepPainter);
#endif

                    std::cout << "left ray: " << rL << std::endl;
                    std::cout << "right ray: " << rR << std::endl;
                    std::cout << "bisector: " << Bs << std::endl;

                    // check for intersections of Bln, Blr and BS
                    bool BsL = false;
                    tFloatVector pBsL;
                    bool BsR = false;
                    tFloatVector pBsR;
                    try {
                        pBsL = Bs.intersection(rL, bounds);
                        if (sl.prj(pBsL) >= cKey) {// check whether IS is before SL
                            BsL = true;
                        }
                    } catch (std::domain_error &e) {}

                    try {
                        pBsR = Bs.intersection(rR, bounds);
                        if (sl.prj(pBsR) >= cKey) {// check whether IS is before SL
                            BsR = true;
                        }
                    } catch (std::domain_error &e) {}

#ifdef WITH_CAIRO
                    if (BsL) {
                        stepPainter.setColor(IS);
                        stepPainter.draw(pBsL, 0);
                    }

                    if (BsR) {
                        stepPainter.setColor(IS);
                        stepPainter.draw(pBsR, 1);
                    }
#endif

                    // check whether rays have extension before deletion: delete delete event
                    if (itBl->ext) {
                        delPQ.push({sl.prj(itBl->ext->p), Event({itBl->ext->p})});
                        std::cout << idx << " deleted Deletion event " << itBl->ext->p << " key: " << sl.prj(itBl->ext->p) << std::endl;
                    }
                    if (itBr->ext) {
                        delPQ.push({sl.prj(itBr->ext->p), Event({itBr->ext->p})});
                        std::cout << idx << " deleted Deletion event " << itBr->ext->p << " key: " << sl.prj(itBr->ext->p) << std::endl;
                    }

                    // remove old Rays from list
                    sl.erase(itBl);
                    auto itInsertPos = sl.erase(itBr);
                    auto itBn = sl.end();

#ifdef WITH_CAIRO
                    basePainter.drawSquare(cPoint.pos);
                    basePainter.drawLine(itBl->p, cPoint.pos);
                    basePainter.drawLine(itBr->p, cPoint.pos);
#endif

                    if (!BsL && !BsR) {
                        std::cout << idx << " case a) no intersection" << std::endl;
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
                            assert(pBsL == cPoint.pos);
                            assert(pBsR == cPoint.pos);
                            std::cout << idx << " case c) both intersect" << std::endl;
                            // boundary originates at v with bisector angle
                            itBn = sl.insert(itInsertPos, Ray({cPoint.pos, aBs, itBl->leftRegion,
                                                               itBr->rightRegion}));
                        } else {
                            std::cout << idx << " case b) one intersection" << std::endl;
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

                    std::cout << "found boundary: " << *itBn << std::endl;

                    // insert intersection points into PQ
                    if (itBn != sl.begin()) {
                        try {
                            auto itL = std::prev(itBn);
                            auto is = itL->intersection(*itBn, bounds);
                            if (sl.prj(is) > cKey) {// only consider point if not yet swept
                                pq.push({sl.prj(is), Event({is, itL, itBn})});
                                std::cout << idx << " added left intersection point at " << is << " key: " << sl.prj(is) << std::endl;
                            }
                        } catch (std::domain_error &e) {}
                    }

                    auto itR = std::next(itBn);
                    if (itR != sl.end()) {
                        try {
                            auto is = itBn->intersection(*itR, bounds);
                            if (sl.prj(is) > cKey) {// only consider point if not yet swept
                                pq.push({sl.prj(is), Event({is, itBn, itR})});
                                std::cout << idx << " added right intersection point at " << is << " key: " << sl.prj(is) << std::endl;
                            }
                        } catch (std::domain_error &e) {}
                    }

#ifdef WITH_CAIRO
                    stepPainter.save("img_k" + std::to_string(k) + "_s" + std::to_string(idx));
#endif

                    break;
                }

                case Event::Type::Deletion: {

                    std::cout << idx << " type: deletion point" << std::endl;

                    // check for point that are not processed on delPQ
                    while (!delPQ.empty() && delPQ.top().first < cKey) {
                        delPQ.pop();
                    }

                    if (!delPQ.empty() && delPQ.top().second.pos == cPoint.pos) {
                        delPQ.pop();
                        std::cout << idx << " previously deleted" << std::endl;
                        break;
                    }

                    auto itB = cPoint.leftRay;// we store the ray to be deleted as left ray
                    assert(itB != sl.end());
                    assert(itB->ext);

#ifdef WITH_CAIRO
                    basePainter.drawLine(itB->p, itB->ext->p);
#endif

                    *itB = *(itB->ext);
                    // we replace the RayUnion with its lower ray, so all intersections pointers should still be valid
                    // all intersections processed after this point will be with lower ray

                    break;
                }
            }

            ++idx;

            std::cout << std::endl;
        }
    }
};
