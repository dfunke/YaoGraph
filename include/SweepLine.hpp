//
// Created by Daniel Funke on 18.03.21.
//

#pragma once

#include <list>
#include <memory>
#include <queue>
#include <sstream>

#include "Predicates.hpp"
#include "Types.hpp"
#include "utils/ListIndexTree.hpp"

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

        // the starting points of upper and lower ray should be different
        assert(!approxEQ(p, ext->p));

        // the starting point of the lower ray must be on the upper rayl;
        assert(approxEQ(ext->p, {ext->p[X], std::tan(angle) * (ext->p[X] - p[X]) + p[Y]}));
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

    struct tIntersectionRetVal {
        bool valid;
        tFloatVector pos;
    };

    tIntersectionRetVal intersection(const Ray &b, const tBox &bounds) const {
        return intersection(*this, b, bounds);
    }

public:
    static tIntersectionRetVal intersection(const Ray &a, const Ray &b, const tBox &bounds) {
        if (!a.ext && !b.ext) {
            return isRR(a, b, bounds);
        } else if (a.ext && b.ext) {
            return isUU(a, b, bounds);
        } else if (a.ext) {
            return isUR(a, b, bounds);
        } else {// b.ext
            return isRU(a, b, bounds);
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
    static tIntersectionRetVal isRR(const Ray &a, const Ray &b, const tBox &bounds) {
        if (std::fmod(std::fmod(a.angle - b.angle, M_PI) + M_PI, M_PI) == 0) {
            return {false, {}};
        }

        if (std::fmod(std::fmod(a.angle, M_PI) + M_PI, M_PI) == M_PI_2) {
            // vertical line this's x
            return {true, {a.p[X], std::tan(b.angle) * (a.p[X] - b.p[X]) + b.p[Y]}};
        } else if (std::fmod(std::fmod(b.angle, M_PI) + M_PI, M_PI) == M_PI_2) {
            // vertical line at b's x
            return {true, {b.p[X], std::tan(a.angle) * (b.p[X] - a.p[X]) + a.p[Y]}};
        }

        tFloat m0 = std::tan(a.angle);// Line 0: y = m0 (x - x0) + y0
        tFloat m1 = std::tan(b.angle);// Line 1: y = m1 (x - x1) + y1

        tFloat x = ((m0 * a.p[X] - m1 * b.p[X]) - (a.p[Y] - b.p[Y])) / (m0 - m1);
        tFloat y = m0 * (x - a.p[X]) + a.p[Y];

        if (!bounds.contains({x, y})) {
            return {false, {}};
        }

        return {true, {x, y}};
    }

    static tIntersectionRetVal isUR(const Ray &ua, const Ray &b, const tBox &bounds) {
        assert(ua.ext);

        // check for intersection of main ray of ua and b
        auto is = isRR(ua, b, bounds);
        // check whether IS is before starting point of extension ray of ur
        if (is.valid && distance2(ua.p, is.pos) < distance2(ua.p, ua.ext->p)) {
            return is;
        }

        // upper ray of ur does not intersect r OR intersection is below starting point of extension ray
        is = isRR(*ua.ext, b, bounds);
        // check whether IS is after starting point of extension ray of ur
        if (is.valid && distance2(ua.p, is.pos) >= distance2(ua.p, ua.ext->p)) {
            return is;
        }

        return {false, {}};
    }

    static tIntersectionRetVal isRU(const Ray &a, const Ray &ub, const tBox &bounds) {
        assert(ub.ext);

        // check for intersection of a and main ray of ub
        auto is = isRR(a, ub, bounds);
        // check whether IS is before starting point of extension ray of ur
        if (is.valid && distance2(ub.p, is.pos) < distance2(ub.p, ub.ext->p)) {
            return is;
        }

        // upper ray of ur does not intersect r OR intersection is below starting point of extension ray
        is = isRR(a, *ub.ext, bounds);
        // check whether IS is after starting point of extension ray of ur
        if (is.valid && distance2(ub.p, is.pos) >= distance2(ub.p, ub.ext->p)) {
            return is;
        }

        return {false, {}};
    }

    static tIntersectionRetVal isUU(const Ray &ua, const Ray &ub, const tBox &bounds) {
        assert(ua.ext && ub.ext);

        // check for intersection of main ray of ua and ub
        auto is = isRR(ua, ub, bounds);
        // check whether IS is before starting point of lowerRay of ua
        // check whether IS is before starting point of lowerRay of ub
        if (is.valid &&
            distance2(ua.p, is.pos) < distance2(ua.p, ua.ext->p) &&
            distance2(ub.p, is.pos) < distance2(ub.p, ub.ext->p)) {

            return is;
        }


        // check for intersection of main ray of ua and lower ray of ub
        is = isRR(ua, *ub.ext, bounds);
        // check whether IS is before starting point of lowerRay of ua
        // check whether IS is after starting point of lowerRay of ub
        if (is.valid &&
            distance2(ua.p, is.pos) < distance2(ua.p, ua.ext->p) &&
            distance2(ub.p, is.pos) >= distance2(ub.p, ub.ext->p)) {
            return is;
        }


        // check for intersection of lower ray of ua and main ray of ub
        is = isRR(*ua.ext, ub, bounds);
        // check whether IS is after starting point of ua's lowerRay
        // check whether IS is before starting point of ub's lowerRay
        if (is.valid &&
            distance2(ua.p, is.pos) >= distance2(ua.p, ua.ext->p) &&
            distance2(ub.p, is.pos) < distance2(ub.p, ub.ext->p)) {
            return is;
        }


        // check for intersection of lower rays of ua and ub
        is = isRR(*ua.ext, *ub.ext, bounds);
        // check whether IS is after starting point of ua's lowerRay
        // check whether IS is after starting point of ub's lowerRay
        if (is.valid &&
            distance2(ua.p, is.pos) >= distance2(ua.p, ua.ext->p) &&
            distance2(ub.p, is.pos) >= distance2(ub.p, ub.ext->p)) {
            return is;
        }

        return {false, {}};
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

    using tRays = SearchTree<Ray>;
    using tRayHandle = typename tRays::Iterator;

    tFloat slDirection;
    tFloatVector slDirVector;
    tRays slRays;

    SweeplineDS(const tFloat slDir_)
        : slDirection(slDir_),
          slDirVector({std::cos(slDirection), std::sin(slDirection)}) {}


    tRayHandle begin() { return slRays.begin(); }
    tRayHandle end() { return slRays.end(); }

    auto insert(tRayHandle &pos, const Ray &ray) {
        return slRays.insert(pos, ray);
    }

    auto erase(tRayHandle &pos) {
        return slRays.erase(pos);
    }

    tFloat prj(const tFloatVector &p) {
        return (p[X] * slDirVector[X] + p[Y] * slDirVector[Y]) /
               dot(slDirVector, slDirVector);
    }

    tRayHandle linearFind(const tFloatVector &p) {// TODO: make const?
        auto exp = begin();
        while (exp != end() && exp->leftOf(p)) {
            ++exp;
        }

        return exp;
    }

    // return iterator to ray immediatly right of point p
    tRayHandle find(const tFloatVector &p) {// TODO: make const?

        auto cmp = [](const tFloatVector &p, const Ray &r) {
            return !r.leftOf(p);
        };

        auto it = slRays.find(p, cmp);

        assert(it == end() || !it->leftOf(p));
        assert(std::prev(it) == end() || std::prev(it)->leftOf(p));
        assert(it == linearFind(p));

        return it;
    }

#ifdef WITH_CAIRO

    void draw(const tFloatVector &pos, Painter &painter) {

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
        for (auto &r : slRays) {
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
        [[maybe_unused]] constexpr std::array<float, 3> SL = {0, 0, 1};
        [[maybe_unused]] constexpr std::array<float, 3> POINT = {1, 0, 0};
        [[maybe_unused]] constexpr std::array<float, 3> IS = {0, 1, 0};
        [[maybe_unused]] constexpr std::array<float, 3> BOUND = {0, 1, 1};

        Painter basePainter(bounds, 1000);
        basePainter.setColor(BASE);
        basePainter.draw(points, true);
#endif
        tIndex idx = 0;
        while (!pq.empty()) {

            auto cKey = pq.top().first;
            auto cPoint = pq.top().second;
            pq.pop();

            // std::cout << idx << ": " << cPoint.idx << " (" << cPoint.pos[X] << ", " << cPoint.pos[Y] << "): key: " << cKey << std::endl;

#if defined(WITH_CAIRO) && defined(PAINT_STEPS)
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

                    // std::cout << idx << " type: input point" << std::endl;

                    auto itBr = sl.find(cPoint.pos);                              // right ray
                    auto itBl = (itBr == sl.begin() ? sl.end() : std::prev(itBr));// left ray
                    assert(itBl == sl.end() || itBl->leftOf(cPoint.pos));

                    // std::cout << idx << " left ray: " << (itBl != sl.end() ? to_string(*itBl) : "NULL") << std::endl;
                    // std::cout << idx << " right ray: " << (itBr != sl.end() ? to_string(*itBr) : "NULL") << std::endl;

                    // add graph edge
                    if (itBr != sl.end() && itBr->leftRegion != tIndex(-1)) {
                        assert(itBl->rightRegion == itBr->leftRegion);

                        graph[cPoint.idx].neighbor[k] = itBr->leftRegion;
                        graph[cPoint.idx].distance[k] = distance2(cPoint.pos, points[itBr->leftRegion]);

                        // std::cout << idx << " edge added: (" << cPoint.idx << ", " << itBr->leftRegion << ") w: " << distance2(cPoint.pos, points[itBr->leftRegion]) << std::endl;
                    }

                    // check if Bl and Br intersect, check only required if they don't originate from same point
                    if (itBl != sl.end() && itBr != sl.end() && itBl->p != itBr->p) {
                        auto is = itBl->intersection(*itBr, bounds);
                        if (is.valid) {
                            delPQ.push({sl.prj(is.pos), Event({is.pos})});
                            // std::cout << idx << " deleted intersection point at " << is.pos << " key: " << sl.prj(is.pos) << std::endl;
                        }
                    }

                    // create new rays and insert them into SL
                    auto itBln = sl.insert(
                            itBr,
                            Ray({cPoint.pos, lTheta + tFloat(M_PI),
                                 itBl != sl.end() ? itBl->rightRegion : tIndex(-1), cPoint.idx}));
                    auto itBrn = sl.insert(
                            itBr, Ray({cPoint.pos, uTheta + tFloat(M_PI), cPoint.idx,
                                       itBr != sl.end() ? itBr->leftRegion : tIndex(-1)}));

                    // std::cout << idx << " left ray " << *itBln << std::endl;
                    // std::cout << idx << " right ray " << *itBrn << std::endl;

                    // insert intersection points into PQ
                    if (itBl != sl.end()) {
                        auto is = itBl->intersection(*itBln, bounds);
                        if (is.valid && sl.prj(is.pos) > cKey) {
                            pq.push({sl.prj(is.pos), Event({is.pos, itBl, itBln})});
                            // std::cout << idx << " added left intersection point at " << is.pos << " key: " << sl.prj(is.pos) << std::endl;
                        }
                    }

                    if (itBr != sl.end()) {
                        auto is = itBrn->intersection(*itBr, bounds);
                        if (is.valid && sl.prj(is.pos) > cKey) {
                            pq.push({sl.prj(is.pos), Event({is.pos, itBrn, itBr})});
                            // std::cout << idx << " added right intersection point at " << is.pos << " key: " << sl.prj(is.pos) << std::endl;
                        }
                    }

                    break;
                }

                case Event::Type::Intersection: {

                    // std::cout << idx << " type: intersection point" << std::endl;

                    // check for point that are not processed on delPQ
                    while (!delPQ.empty() && delPQ.top().first < cKey) {
                        delPQ.pop();
                    }

                    if (!delPQ.empty() && delPQ.top().second.pos == cPoint.pos) {
                        delPQ.pop();
                        // std::cout << idx << " previously deleted" << std::endl;
                        break;
                    }

                    auto itBl = cPoint.leftRay;
                    auto itBr = cPoint.rightRay;

#if defined(WITH_CAIRO) && defined(PAINT_STEPS)
                    stepPainter.setColor(BOUND);
                    itBl->draw(stepPainter);
                    itBr->draw(stepPainter);
                    stepPainter.save("img_k" + std::to_string(k) + "_s" + std::to_string(idx));
#endif

                    assert(std::next(itBl) == itBr);
                    assert(itBl == std::prev(itBr));

                    // std::cout << idx << " left ray: " << (itBl != sl.end() ? to_string(*itBl) : "NULL") << std::endl;
                    // std::cout << idx << " right ray: " << (itBr != sl.end() ? to_string(*itBr) : "NULL") << std::endl;

                    // delete intersection points from PQ
                    if (itBl != sl.begin()) {
                        auto itBll = std::prev(itBl);

                        if (itBl->p != itBll->p) {
                            auto is = itBll->intersection(*itBl, bounds);
                            if (is.valid) {
                                delPQ.push({sl.prj(is.pos), Event({is.pos})});
                                // std::cout << idx << " deleted left intersection point at " << is.pos << " key: " << sl.prj(is.pos) << std::endl;
                            }
                        }
                    }

                    if (itBr != sl.end()) {
                        auto itBrr = std::next(itBr);

                        if (itBrr != sl.end() && itBr->p != itBrr->p) {
                            auto is = itBr->intersection(*itBrr, bounds);
                            if (is.valid) {
                                delPQ.push({sl.prj(is.pos), Event({is.pos})});
                                // std::cout << idx << " deleted right intersection point at " << is.pos << " key: " << sl.prj(is.pos) << std::endl;
                            }
                        }
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
                    aBs += (((sl.slDirection - aBs) < M_PI) ? 1 : -1) * M_PI_2;// TODO: check angle orientation
                    Ray Bs({pMid, aBs, itBl->leftRegion, itBr->rightRegion});

#if defined(WITH_CAIRO) && defined(PAINT_STEPS)
                    stepPainter.setColor(IS);
                    rL.draw(stepPainter);
                    rR.draw(stepPainter);
                    stepPainter.drawLine(pL, pR);
                    stepPainter.drawSquare(pMid);
                    Bs.draw(stepPainter);
#endif

                    // std::cout << idx << " left ray: " << rL << std::endl;
                    // std::cout << idx << " right ray: " << rR << std::endl;
                    // std::cout << idx << " bisector: " << Bs << std::endl;

                    // check for intersections of Bln, Blr and BS
                    auto BsL = Bs.intersection(rL, bounds);
                    auto BsR = Bs.intersection(rR, bounds);

                    if (BsL.valid && BsR.valid && approxEQ(BsL.pos, BsR.pos) && approxEQ(BsL.pos, cPoint.pos)) { //TODO better way to test?
                        assert(approxEQ(BsR.pos, cPoint.pos));
                        // we don't check for sweepline projection to avoid floating point problems
                    } else {
                        if (BsL.valid && approxLT(sl.prj(BsL.pos), cKey)) {// check whether IS is before SL
                            BsL.valid = false;
                        }
                        if (BsR.valid && approxLT(sl.prj(BsR.pos), cKey)) {// check whether IS is before SL
                            BsR.valid = false;
                        }
                    }

#if defined(WITH_CAIRO) && defined(PAINT_STEPS)
                    if (BsL.valid) {
                        stepPainter.setColor(IS);
                        stepPainter.draw(BsL.pos, 0);
                    }

                    if (BsR.valid) {
                        stepPainter.setColor(IS);
                        stepPainter.draw(BsR.pos, 1);
                    }
#endif

                    // check whether rays have extension before deletion: delete delete event
                    if (itBl->ext) {
                        delPQ.push({sl.prj(itBl->ext->p), Event({itBl->ext->p})});
                        // std::cout << idx << " deleted Deletion event " << itBl->ext->p << " key: " << sl.prj(itBl->ext->p) << std::endl;
                    }
                    if (itBr->ext) {
                        delPQ.push({sl.prj(itBr->ext->p), Event({itBr->ext->p})});
                        // std::cout << idx << " deleted Deletion event " << itBr->ext->p << " key: " << sl.prj(itBr->ext->p) << std::endl;
                    }

#ifdef WITH_CAIRO
                    basePainter.drawSquare(cPoint.pos);
                    basePainter.drawLine(itBl->p, cPoint.pos);
                    basePainter.drawLine(itBr->p, cPoint.pos);
#endif

                    // remove old Rays from list
                    sl.erase(itBl);
                    auto itInsertPos = sl.erase(itBr);
                    auto itBn = sl.end();

                    if (!BsL.valid && !BsR.valid) {
                        // std::cout << idx << " case a) no intersection" << std::endl;
                        // bisector intersects no ray from P
                        if (sl.prj(pR) < sl.prj(pL)) {
                            // right point is further from v -> right ray becomes border
                            itBn = sl.insert(itInsertPos, rR);
                        } else {
                            // left point is further from v -> left ray becomes border
                            itBn = sl.insert(itInsertPos, rL);
                        }
                    } else {
                        if (BsL.valid && BsR.valid) {
                            // bisector intersects both rays - > must be in same point v
                            assert(approxEQ(BsL.pos, BsR.pos));
                            assert(approxEQ(BsL.pos, cPoint.pos));
                            assert(approxEQ(BsR.pos, cPoint.pos));
                            // std::cout << idx << " case c) both intersect" << std::endl;
                            // boundary originates at v with bisector angle
                            Bs.p = cPoint.pos;
                            itBn = sl.insert(itInsertPos, Bs);
                        } else {
                            // std::cout << idx << " case b) one intersection" << std::endl;
                            if (BsL.valid) {
                                // bisector intersects left ray
                                // boundary to intersection point, then ray with bisector angle
                                Bs.p = BsL.pos;
                                if (!approxEQ(rL.p, Bs.p)) {
                                    itBn = sl.insert(itInsertPos, Ray({rL, Bs}));
                                    pq.push({sl.prj(Bs.p), Event({Bs.p, itBn})});
                                } else {// if they are almost equal immediately use bisector ray
                                    itBn = sl.insert(itInsertPos, Bs);
                                }
                            } else {
                                // bisector intersects right ray
                                // boundary to intersection point, then ray with bisector angle
                                Bs.p = BsR.pos;
                                if (!approxEQ(rR.p, Bs.p)) {
                                    itBn = sl.insert(itInsertPos, Ray({rR, Bs}));
                                    pq.push({sl.prj(Bs.p), Event({Bs.p, itBn})});
                                } else {// if they are almost equal immediately use bisector ray
                                    itBn = sl.insert(itInsertPos, Bs);
                                }
                            }
                        }
                    }
                    assert(itBn != sl.end());// some boundary was inserted into SL

                    // std::cout << "found boundary: " << *itBn << std::endl;

                    // insert intersection points into PQ
                    if (itBn != sl.begin()) {
                        auto itL = std::prev(itBn);
                        auto is = itL->intersection(*itBn, bounds);
                        if (is.valid && sl.prj(is.pos) > cKey) {// only consider point if not yet swept
                            pq.push({sl.prj(is.pos), Event({is.pos, itL, itBn})});
                            // std::cout << idx << " added left intersection point at " << is.pos << " key: " << sl.prj(is.pos) << std::endl;
                        }
                    }

                    auto itR = std::next(itBn);
                    if (itR != sl.end()) {
                        auto is = itBn->intersection(*itR, bounds);
                        if (is.valid && sl.prj(is.pos) > cKey) {// only consider point if not yet swept
                            pq.push({sl.prj(is.pos), Event({is.pos, itBn, itR})});
                            // std::cout << idx << " added right intersection point at " << is.pos << " key: " << sl.prj(is.pos) << std::endl;
                        }
                    }

#if defined(WITH_CAIRO) && defined(PAINT_STEPS)
                    stepPainter.save("img_k" + std::to_string(k) + "_s" + std::to_string(idx));
#endif

                    break;
                }

                case Event::Type::Deletion: {

                    // std::cout << idx << " type: deletion point" << std::endl;

                    // check for point that are not processed on delPQ
                    while (!delPQ.empty() && delPQ.top().first < cKey) {
                        delPQ.pop();
                    }

                    if (!delPQ.empty() && delPQ.top().second.pos == cPoint.pos) {
                        delPQ.pop();
                        // std::cout << idx << " previously deleted" << std::endl;
                        break;
                    }

                    auto itB = cPoint.leftRay;// we store the ray to be deleted as left ray
                    assert(itB != sl.end());
                    assert(itB->ext);
                    assert(itB->leftRegion == itB->ext->leftRegion);
                    assert(itB->rightRegion == itB->ext->rightRegion);
                    assert(!itB->ext->ext);

#ifdef WITH_CAIRO
                    basePainter.drawLine(itB->p, itB->ext->p);
#endif

//                    std::cout << idx << " old ray: " << *itB << std::endl;

                    // we replace the RayUnion with its lower ray, so all intersections pointers should still be valid
                    // all intersections processed after this point will be with lower ray
                    itB->p = itB->ext->p;
                    itB->angle = itB->ext->angle;
                    // and delete the lower ray
                    itB->ext.reset();

//                    std::cout << idx << " new ray: " << *itB << std::endl;

                    assert(!itB->ext);

                    break;
                }
            }

            ++idx;

            // std::cout << std::endl;
        }

#ifdef WITH_CAIRO
        basePainter.save("img_k" + std::to_string(k));
#endif
    }
};
