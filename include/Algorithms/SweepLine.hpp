//
// Created by Daniel Funke on 18.03.21.
//

#pragma once

#include <functional>
#include <list>
#include <memory>
#include <queue>
#include <sstream>
#include <unordered_map>

#include "Predicates.hpp"
#include "Types.hpp"
#include "Utils/ListIndexTree.hpp"
#include "Utils/PriorityQueue.hpp"

#include "Utils/Logging.hpp"

#include "Kernels/InexactKernel.hpp"
#ifdef WITH_CGAL
#include "Kernels/CGALKernel.hpp"
#endif

#include "contrib/robin_hood.h"

#ifdef WITH_CAIRO
#include "Painter.hpp"
//#define PAINT_STEPS
#endif

template<typename Kernel>
class SweepLine {

    using tEFloat = typename Kernel::Float;// possibly exact float
    using tDirection = typename Kernel::Direction;
    using tPoint = typename Kernel::Point;

    using tRay = typename Kernel::Ray;

    struct SweeplineDS {

        using tRays = SearchTree<tRay>;
        using tRayHandle = typename tRays::Iterator;

        tDirection slDirection;
        tRays slRays;

        SweeplineDS(const tDirection slDir_)
            : slDirection(slDir_) {}


        tRayHandle begin() { return slRays.begin(); }
        tRayHandle end() { return slRays.end(); }

        auto insert(tRayHandle &pos, const tRay &ray) {
            return slRays.insert(pos, ray);
        }

        auto insert_pair(tRayHandle &pos, const tRay &left, const tRay &right) {
            return slRays.insert_pair(pos, left, right);
        }

        auto erase(tRayHandle &pos) {
            return slRays.erase(pos);
        }

        auto replace(tRayHandle &pos, const tRay &ray) {
            return slRays.replace(pos, ray);
        }

        tEFloat prj(const tPoint &p) {
            return slDirection.prj(p);
        }

        tRayHandle linearFind(const tPoint &p) {// TODO: make const?
            auto exp = begin();
            while (exp != end() && exp->orientedSide(p) == tOrientedSide::LEFT) {
                ++exp;
            }

            return exp;
        }

        // return iterator to ray immediatly right of point p
        tRayHandle find(const tPoint &p) {// TODO: make const?

            auto cmp = [](const tPoint &p, const tRay &r) {
                return r.orientedSide(p) != tOrientedSide::LEFT;
            };

            auto it = slRays.find(p, cmp);

#ifndef NDEBUG
            //TODO debug code
            auto itLin = linearFind(p);
            if (it != itLin) {
                std::cout << (it != end() ? to_string(*it) : "NULL") << std::endl;
                std::cout << (itLin != end() ? to_string(*itLin) : "NULL") << std::endl;
            }
#endif

            ASSERT(it == end() || it->orientedSide(p) != tOrientedSide::LEFT);
            ASSERT(std::prev(it) == end() || std::prev(it)->orientedSide(p) != tOrientedSide::RIGHT);
            ASSERT(it == linearFind(p));

            return it;
        }

#ifdef WITH_CAIRO

        void draw(const tPoint &pos, Painter &painter) {

            // draw sweepline
            painter.setDash();
            painter.drawLine(pos, {pos[X] + std::cos(slDirection.angle() + tIFloat(M_PI_2)),
                                   pos[Y] + std::sin(slDirection.angle() + tIFloat(M_PI_2))});
            painter.drawLine(pos, {pos[X] + std::cos(slDirection.angle() - tIFloat(M_PI_2)),
                                   pos[Y] + std::sin(slDirection.angle() - tIFloat(M_PI_2))});
            painter.unsetDash();

            // draw direction
            painter.setDash();
            //painter.drawLine(pos, {pos[X] + slDirVector[X], pos[Y] + slDirVector[Y]});
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

        Event(const tPoint &pos_)
            : type(Type::Deletion), pos(pos_), idx(INF_IDX) {}

        Event(const tPoint &pos_, const tIndex &idx_)
            : type(Type::Input), pos(pos_), idx(idx_) {}

        Event(const tPoint &pos_, const tRayHandle &left)
            : type(Type::Deletion), pos(pos_), idx(INF_IDX), leftRay(left) {}

        Event(const tPoint &pos_, const tRayHandle &left,
              const tRayHandle &right)
            : type(Type::Intersection), pos(pos_), idx(INF_IDX), leftRay(left),
              rightRay(right) {}

        Type type;
        tPoint pos;

        tIndex idx;

        tRayHandle leftRay;
        tRayHandle rightRay;
    };

public:
    static std::string name() {
        return "Sweepline_" + Kernel::name();
    }

    auto operator()(const tDim &K, const tPoints &points, const tBox &bounds) const {
        tYaoGraph g(points.size(), K);

        auto rays = Kernel::computeCones(K);

        for (tDim k = 0; k < K; ++k) {
            sweepline(points, k, K, rays, g, bounds);
        }

        return g;
    }

private:
    // Only for pairs of std::hash-able types for simplicity.
    // You can of course template this struct to allow other hash functions
    struct pair_hash {
        template<class IT>
        std::size_t operator()(const std::pair<IT, IT> &p) const {
            auto h1 = std::hash<typename IT::const_it_pointer>{}(p.first.addr());
            auto h2 = std::hash<typename IT::const_it_pointer>{}(p.second.addr());

            // from Boost HashCombine
            return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h2 >> 2));
        }
    };

    struct it_hash {
        template<class IT>
        std::size_t operator()(const IT &p) const {
            return std::hash<typename IT::const_it_pointer>{}(p.addr());
        }
    };

    void sweepline(const tPoints &iPoints, tDim k, const tDim &K, const auto &rays, tYaoGraph &graph,
                   const tBox &iBounds) const {

        using pqItem = std::pair<tEFloat, Event>;
        auto pqCmp = [](const pqItem &a, const pqItem &b) {
            return a.first > b.first;
        };

        using PQ = PriQueueAdapter<pqItem, decltype(pqCmp)>;
        using isKey = std::pair<typename Event::tRayHandle, typename Event::tRayHandle>;
        robin_hood::unordered_map<isKey, typename PQ::handle, pair_hash> isMap;
        robin_hood::unordered_map<typename Event::tRayHandle, typename PQ::handle, it_hash> extMap;

        auto bounds = Kernel::mkBBox(iBounds);

        tDirection lRay(-rays[k]);          // range [0..2PI]
        tDirection rRay(-rays[(k + 1) % K]);// range [0..2PI]
        //ASSERT(lRay.angle() <= rRay.angle());

        tDirection slDir(Kernel::Bisector(lRay, rRay));// range [0..2PI]
        //ASSERT(lRay.angle() <= slDir.angle() && slDir.angle() <= rRay.angle());

        SweeplineDS sl(slDir);// range [0..2PI]

        typename PQ::arrType pqArr;
        pqArr.reserve(iPoints.size());
        for (tIndex i = 0; i < iPoints.size(); ++i) {
            auto p = Kernel::mkPoint(iPoints[i]);
            pqArr.push_back({sl.prj(p), Event(p, i)});
        }

        PQ pq(std::move(pqArr), std::sqrt(iPoints.size()));

#ifdef WITH_CAIRO

        [[maybe_unused]] constexpr std::array<float, 3> BLACK = {0, 0, 0};
        [[maybe_unused]] constexpr std::array<float, 3> BLUE = {0, 0, 1};
        [[maybe_unused]] constexpr std::array<float, 3> RED = {1, 0, 0};
        [[maybe_unused]] constexpr std::array<float, 3> GREEN = {0, 1, 0};
        [[maybe_unused]] constexpr std::array<float, 3> CYAN = {0, 1, 1};
        [[maybe_unused]] constexpr std::array<float, 3> YELLOW = {1, 1, 0};
        [[maybe_unused]] constexpr std::array<float, 3> ORANGE = {.8, .4, 0};
        [[maybe_unused]] constexpr std::array<float, 3> PINK = {.7, .3, .7};

        Painter basePainter(iBounds, 1000);
        basePainter.setColor(BLACK);
        basePainter.draw(iPoints, true);
#endif

        LOG("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl);
        LOG("Processing cone " << k << " (" << lRay.angle() << ", " << rRay.angle() << ") - sweepline " << slDir.angle() << std::endl);

        tIndex idx = 0;
        while (!pq.empty()) {

            auto cKey = pq.top().first;
            auto cPoint = pq.top().second;
            //            pq.pop(); // we need the pq for drawing

            LOG(idx << ": "
                    << cPoint.idx << " (" << cPoint.pos << "): key: " << cKey << std::endl);

#if defined(WITH_CAIRO) && defined(PAINT_STEPS)
            //basePainter.draw(cPoint.pos, idx);

            Painter stepPainter(basePainter);

            stepPainter.setColor(RED);
            stepPainter.draw(cPoint.pos, idx, false);

            stepPainter.setColor(BLUE);
            sl.draw(cPoint.pos, stepPainter);

            stepPainter.setColor(ORANGE);
            for (const auto &is : isMap) {
                stepPainter.draw(pq.get_key(is.second).second.pos, 0, false);
            }

            stepPainter.setColor(PINK);
            for (const auto &ext : extMap) {
                stepPainter.draw(pq.get_key(ext.second).second.pos, 0, false);
            }

            std::ostringstream stepFilename;
            stepFilename << "img_k" << k << "_s" << std::setfill('0')
                         << std::setw(static_cast<int>(std::ceil(std::log10(iPoints.size()) + 2))) << idx;
            stepPainter.save(stepFilename.str());
#endif
            // now pop PQ after drawing
            pq.pop();

            auto handleRayIntersection = [&isMap, &extMap, &pq, &bounds, &sl, &idx, &cKey](auto &leftRay, auto &rightRay, bool left) {
                auto is = leftRay->intersection(*rightRay, bounds);
                if (is) {
                    const tPoint *p = std::get_if<tPoint>(&*is);
                    if (p && sl.prj(*p) > cKey) {
                        isMap[std::make_pair(leftRay, rightRay)] = pq.push({sl.prj(*p), Event({*p, leftRay, rightRay})});
                        LOG(idx << ": "
                                << " added " << (left ? "left" : "right")
                                << " intersection point at (" << *p << ") key: " << sl.prj(*p) << std::endl);
                    }
                    //TODO handle ray
                    tRay *r = std::get_if<tRay>(&*is);
                    if (r) {
                        if (r->direction() != leftRay->direction()) {
                            r->invert();
                        }

                        LOG(idx << ": "
                                << "RI left: " << *leftRay << std::endl);
                        LOG(idx << ": "
                                << "RI right: " << *rightRay << std::endl);
                        LOG(idx << ": "
                                << "RI IS Ray: " << *r << std::endl);

                        ASSERT(leftRay->rightRegion == rightRay->leftRegion);
                        ASSERT(r->direction() == leftRay->direction() && r->direction() == rightRay->direction());

			//TODO handle extensions correctly

                        if (left) {
                            if (rightRay->isExtended()) {
                                auto pqR = extMap.find(rightRay);
                                ASSERT(pqR != extMap.end());
                                pq.remove(pqR->second);
                                LOG(idx << ": "
                                        << " RI deleted Deletion event " << rightRay->extOrigin() << " key: " << sl.prj(rightRay->extOrigin()) << std::endl);
                                extMap.erase(pqR);
                            }

                            // remove rightRay, as leftRay replaces it
                            // point iterator to right ray to left ray

                            leftRay->leftRegion = leftRay->leftRegion;
                            leftRay->rightRegion = rightRay->rightRegion;
                            sl.erase(rightRay);
                            rightRay = leftRay;
                        } else {
                            if (leftRay->isExtended()) {
                                auto pqL = extMap.find(leftRay);
                                ASSERT(pqL != extMap.end());
                                pq.remove(pqL->second);
                                LOG(idx << ": "
                                        << " RI deleted Deletion event " << leftRay->extOrigin() << " key: " << sl.prj(leftRay->extOrigin()) << std::endl);
                                extMap.erase(pqL);
                            }

                            // remove leftRay, as rightRay replaces it
                            // point iterator to left ray to right ray

                            rightRay->leftRegion = leftRay->leftRegion;
                            rightRay->rightRegion = rightRay->rightRegion;
                            sl.erase(leftRay);
                            leftRay = rightRay;
                        }

                        LOG(idx << ": "
                                << "new Ray: " << (left ? *leftRay : *rightRay) << std::endl);
                    }
                }
            };

            switch (cPoint.type) {
                case Event::Type::Input: {

                    LOG(idx << ": "
                            << " type: input point" << std::endl);

                    auto itBr = sl.find(cPoint.pos);// right ray

                    if (itBr != sl.end() && itBr->orientedSide(cPoint.pos) == tOrientedSide::LINE) {// check whether point is ON right boundary
                        if (itBr->direction() == lRay) {
                            LOG(idx << ": point ON right boundary that is a left ray" << std::endl);
                            LOG(idx << ": old Ray " << *itBr << std::endl);

                            itBr = std::next(itBr);

                            LOG(idx << ": new Ray " << *itBr << std::endl);
                        }
                    }

                    auto itBl = (itBr == sl.begin() ? sl.end() : std::prev(itBr));// left ray
                    ASSERT(itBl == sl.end() || itBl->orientedSide(cPoint.pos) != tOrientedSide::RIGHT);

                    LOG(idx << ": "
                            << " left ray: " << (itBl != sl.end() ? to_string(*itBl) : "NULL") << std::endl);
                    LOG(idx << ": "
                            << " right ray: " << (itBr != sl.end() ? to_string(*itBr) : "NULL") << std::endl);

#if defined(WITH_CAIRO) && defined(PAINT_STEPS)
                    stepPainter.setColor(CYAN);
                    if (itBl != sl.end()) itBl->draw(stepPainter);
                    stepPainter.setColor(GREEN);
                    if (itBr != sl.end()) itBr->draw(stepPainter);
                    stepPainter.save(stepFilename.str());
#endif

                    // add graph edge
                    // if p lies on right boundary, it does not belong to this cone
                    if (itBr != sl.end() && itBr->leftRegion != INF_IDX) {
                        ASSERT(itBl->rightRegion == itBr->leftRegion);

                        auto r = itBr->orientedSide(cPoint.pos) != tOrientedSide::LINE ? itBr->leftRegion : itBr->rightRegion;
                        if (r != INF_IDX) {
                            graph[cPoint.idx].neighbor[k] = r;
                            graph[cPoint.idx].distance[k] = Kernel::to_float(Kernel::distance2(cPoint.pos, Kernel::mkPoint(iPoints[r])));
                        }

                        LOG(idx << ": "
                                << " edge added: (" << cPoint.idx << ", " << itBr->leftRegion << ") w: " << Kernel::distance2(cPoint.pos, Kernel::mkPoint(iPoints[itBr->leftRegion])) << std::endl);

#ifdef WITH_CAIRO
                        basePainter.drawLine(itBr->origin(), cPoint.pos);
#endif
                    }

                    // check if Bl and Br intersect, check only required if they don't originate from same point
                    if (itBl != sl.end() && itBr != sl.end()) {
                        auto itIs = isMap.find(std::make_pair(itBl, itBr));
                        if (itIs != isMap.end()) {
                            LOG(idx << ": "
                                    << " deleted intersection point between " << *itIs->first.first << " and " << *itIs->first.second << std::endl);
                            pq.remove(itIs->second);
                            isMap.erase(itIs);
                        }
                    }

                    // create new rays and insert them into SL
                    auto itBln = sl.insert_pair(
                            itBr,
                            tRay({cPoint.pos, lRay,
                                  itBl != sl.end() ? itBl->rightRegion : INF_IDX, cPoint.idx}),
                            tRay({cPoint.pos, rRay, cPoint.idx,
                                  itBr != sl.end() ? itBr->leftRegion : INF_IDX}));
                    auto itBrn = std::next(itBln);

                    ASSERT(itBln == std::prev(itBrn));
                    ASSERT(itBrn == std::next(itBln));

                    LOG(idx << ": "
                            << " left ray " << *itBln << std::endl);
                    LOG(idx << ": "
                            << " right ray " << *itBrn << std::endl);

                    // insert intersection points into PQ
                    if (itBl != sl.end()) {
                        handleRayIntersection(itBl, itBln, true);
                    }

                    if (itBr != sl.end()) {
                        handleRayIntersection(itBrn, itBr, false);
                    }

                    break;
                }

                case Event::Type::Intersection: {

                    LOG(idx << ": "
                            << " type: intersection point" << std::endl);

                    auto itBl = cPoint.leftRay;
                    auto itBr = cPoint.rightRay;

#if defined(WITH_CAIRO) && defined(PAINT_STEPS)
                    stepPainter.setColor(CYAN);
                    itBl->draw(stepPainter);
                    itBr->draw(stepPainter);
                    stepPainter.save(stepFilename.str());
#endif

                    ASSERT(std::next(itBl) == itBr);
                    ASSERT(itBl == std::prev(itBr));

                    LOG(idx << ": "
                            << " left ray: " << (itBl != sl.end() ? to_string(*itBl) : "NULL") << std::endl);
                    LOG(idx << ": "
                            << " right ray: " << (itBr != sl.end() ? to_string(*itBr) : "NULL") << std::endl);

                    // delete intersection point from hash map
                    [[maybe_unused]] auto chk = isMap.erase(std::make_pair(itBl, itBr));
                    ASSERT(chk);

                    // delete intersection points from PQ
                    if (itBl != sl.begin()) {
                        auto itBll = std::prev(itBl);
                        ASSERT(itBll != sl.end());

                        auto itIs = isMap.find(std::make_pair(itBll, itBl));
                        if (itIs != isMap.end()) {
                            LOG(idx << ": "
                                    << " deleted intersection point between " << *itIs->first.first << " and " << *itIs->first.second << std::endl);
                            pq.remove(itIs->second);
                            isMap.erase(itIs);
                        }
                    }

                    if (itBr != sl.end()) {
                        auto itBrr = std::next(itBr);

                        if (itBrr != sl.end()) {
                            auto itIs = isMap.find(std::make_pair(itBr, itBrr));
                            if (itIs != isMap.end()) {
                                LOG(idx << ": "
                                        << " deleted intersection point between " << *itIs->first.first << " and " << *itIs->first.second << std::endl);
                                pq.remove(itIs->second);
                                isMap.erase(itIs);
                            }
                        }
                    }

                    tRay rL({cPoint.pos, lRay, itBl->leftRegion, itBr->rightRegion});
                    tRay rR({cPoint.pos, rRay, itBl->leftRegion, itBr->rightRegion});

                    // bisector line between A and B
                    auto pL = Kernel::mkPoint(iPoints[itBl->leftRegion]);
                    auto pR = Kernel::mkPoint(iPoints[itBr->rightRegion]);

                    auto pMid = Kernel::Midpoint(pL, pR);
                    auto aBs = Kernel::Bisector(pL, pR, sl.slDirection);

                    // ASSERT(Kernel::approxEQ(pMid, Kernel::Midpoint(pR, pL)));
                    // ASSERT(Kernel::approxEQ(aBs.angle(), Kernel::Bisector(pR, pL, sl.slDirection).angle()));

                    tRay Bs(pMid, aBs, itBl->leftRegion, itBr->rightRegion);

#if defined(WITH_CAIRO) && defined(PAINT_STEPS)
                    stepPainter.setColor(GREEN);
                    rL.draw(stepPainter);
                    rR.draw(stepPainter);

                    stepPainter.setColor(YELLOW);
                    stepPainter.drawLine(pL, pR);
                    stepPainter.drawSquare(pMid);
                    Bs.draw(stepPainter);
#endif

                    LOG(idx << ": "
                            << " left ray: " << rL << std::endl);
                    LOG(idx << ": "
                            << " right ray: " << rR << std::endl);
                    LOG(idx << ": "
                            << " bisector: " << Bs << std::endl);

                    // check for intersections of Bln, Blr and BS
                    auto BsL = Bs.intersection(rL, bounds);
                    const tPoint *pBsL = nullptr;
                    if (BsL) {
                        ASSERT(std::holds_alternative<tPoint>(*BsL));
                        pBsL = std::get_if<tPoint>(&*BsL);
                    }

                    auto BsR = Bs.intersection(rR, bounds);
                    const tPoint *pBsR = nullptr;
                    if (BsR) {
                        ASSERT(std::holds_alternative<tPoint>(*BsR));
                        pBsR = std::get_if<tPoint>(&*BsR);
                    }

                    //TODO can these be rays as well?

                    if (pBsL && pBsR && Kernel::approxEQ(*pBsL, *pBsR) && Kernel::approxEQ(*pBsL, cPoint.pos)) {//TODO better way to test?
                        ASSERT(Kernel::approxEQ(*pBsR, cPoint.pos));
                        // we don't check for sweepline projection to avoid floating point problems
                    } else {
                        if (pBsL && Kernel::approxLT(sl.prj(*pBsL), cKey)) {// check whether IS is before SL
                            pBsL = nullptr;
                        }
                        if (pBsR && Kernel::approxLT(sl.prj(*pBsR), cKey)) {// check whether IS is before SL
                            pBsR = nullptr;
                        }
                    }

#if defined(WITH_CAIRO) && defined(PAINT_STEPS)
                    if (pBsL) {
                        stepPainter.setColor(RED);
                        stepPainter.draw(*pBsL, 0);
                    }

                    if (pBsR) {
                        stepPainter.setColor(RED);
                        stepPainter.draw(*pBsR, 1);
                    }
#endif

                    // check whether rays have extension before deletion: delete delete event
                    if (itBl->isExtended()) {
                        auto pqBl = extMap.find(itBl);
                        ASSERT(pqBl != extMap.end());
                        pq.remove(pqBl->second);
                        LOG(idx << ": "
                                << " deleted Deletion event " << itBl->extOrigin() << " key: " << sl.prj(itBl->extOrigin()) << std::endl);
                        extMap.erase(pqBl);
                    }
                    if (itBr->isExtended()) {
                        auto pqBr = extMap.find(itBr);
                        ASSERT(pqBr != extMap.end());
                        pq.remove(pqBr->second);
                        LOG(idx << ": "
                                << " deleted Deletion event " << itBr->extOrigin() << " key: " << sl.prj(itBr->extOrigin()) << std::endl);
                        extMap.erase(pqBr);
                    }

#ifdef WITH_CAIRO
                    basePainter.drawSquare(cPoint.pos);
                    basePainter.drawLine(itBl->origin(), cPoint.pos);
                    basePainter.drawLine(itBr->origin(), cPoint.pos);
#endif

                    // remove old Rays from list
                    sl.erase(itBl);
                    //auto itInsertPos = sl.erase(itBr);
                    auto itBn = sl.end();

                    if (!pBsL && !pBsR) {
                        LOG(idx << ": "
                                << " case a) no intersection" << std::endl);
                        // bisector intersects no ray from P
                        if (sl.prj(pR) < sl.prj(pL)) {
                            // right point is further from v -> right ray becomes border
                            itBn = sl.replace(itBr, rR);
                        } else {
                            // left point is further from v -> left ray becomes border
                            itBn = sl.replace(itBr, rL);
                        }
                    } else {
                        if (pBsL && pBsR) {
                            // bisector intersects both rays - > must be in same point v
                            ASSERT(Kernel::approxEQ(*pBsL, *pBsR));
                            ASSERT(Kernel::approxEQ(*pBsL, cPoint.pos));
                            ASSERT(Kernel::approxEQ(*pBsR, cPoint.pos));
                            LOG(idx << ": "
                                    << " case c) both intersect" << std::endl);
                            // boundary originates at v with bisector angle
                            Bs.setOrigin(cPoint.pos);
                            itBn = sl.replace(itBr, Bs);
                        } else {
                            LOG(idx << ": "
                                    << " case b) one intersection" << std::endl);
                            if (pBsL) {
                                // bisector intersects left ray
                                // boundary to intersection point, then ray with bisector angle
                                Bs.setOrigin(*pBsL);
                                if (!Kernel::approxEQ(rL.origin(), Bs.origin())) {
                                    itBn = sl.replace(itBr, tRay({rL, Bs}));
                                    extMap[itBn] = pq.push({sl.prj(Bs.origin()), Event(Bs.origin(), itBn)});
                                    LOG(idx << ": "
                                            << " added deletion point at (" << Bs.origin() << ") key: " << sl.prj(Bs.origin())
                                            << " for " << *itBn << std::endl);
                                } else {// if they are almost equal immediately use bisector ray
                                    itBn = sl.replace(itBr, Bs);
                                }
                            } else {
                                // bisector intersects right ray
                                // boundary to intersection point, then ray with bisector angle
                                Bs.setOrigin(*pBsR);
                                if (!Kernel::approxEQ(rR.origin(), Bs.origin())) {
                                    itBn = sl.replace(itBr, tRay({rR, Bs}));
                                    extMap[itBn] = pq.push({sl.prj(Bs.origin()), Event(Bs.origin(), itBn)});
                                    LOG(idx << ": "
                                            << " added deletion point at (" << Bs.origin() << ") key: " << sl.prj(Bs.origin())
                                            << " for " << *itBn << std::endl);
                                } else {// if they are almost equal immediately use bisector ray
                                    itBn = sl.replace(itBr, Bs);
                                }
                            }
                        }
                    }
                    ASSERT(itBn != sl.end());// some boundary was inserted into SL

                    LOG(idx << ": "
                            << "found boundary: " << *itBn << std::endl);

                    // insert intersection points into PQ

                    auto saveItBn = itBn;
                    auto saveItL = itBn;
                    if (itBn != sl.begin()) {
                        auto itL = std::prev(itBn);
                        saveItL = itL;
                        handleRayIntersection(itL, itBn, true);
                    }

                    auto itR = std::next(itBn);
                    if (itR != sl.end()) {
                        handleRayIntersection(itBn, itR, false);
                    }

                    //check if itBn was modified and re-route intersection points
                    if (itBn != saveItBn && saveItL != saveItBn) {
                        // saveItL is different from saveItBn so there was a previous IT

                        // the left neighbor should still be the same with the new itBn
                        ASSERT(saveItL == std::prev(itBn));

                        auto itIs = isMap.find(std::make_pair(saveItL, saveItBn));
                        if (itIs != isMap.end()) {
                            auto oldEvent = pq.get_key(itIs->second);
                            pq.remove(itIs->second);
                            isMap.erase(itIs);

                            isMap[std::make_pair(saveItL, itBn)] = pq.push({sl.prj(oldEvent.second.pos),
                                                                            Event({oldEvent.second.pos, saveItL, itBn})});

                            LOG(idx << ": "
                                    << " updated intersection point between " << *saveItL << " and " << *itBn << std::endl);
                        }
                    }

#if defined(WITH_CAIRO) && defined(PAINT_STEPS)
                    stepPainter.save(stepFilename.str());
#endif

                    break;
                }

                case Event::Type::Deletion: {

                    LOG(idx << ": "
                            << " type: deletion point" << std::endl);

                    auto itB = cPoint.leftRay;// we store the ray to be deleted as left ray
                    ASSERT(itB != sl.end());
                    ASSERT(itB->isExtended());
                    //ASSERT(itB->leftRegion == itB->ext->leftRegion);
                    //ASSERT(itB->rightRegion == itB->ext->rightRegion);
                    //ASSERT(!itB->ext->ext);

                    // delete extension point from hash map
                    [[maybe_unused]] auto chk = extMap.erase(itB);
                    ASSERT(chk);

#ifdef WITH_CAIRO
                    basePainter.drawLine(itB->origin(), itB->extOrigin());
#endif

                    LOG(idx << ": "
                            << " old ray: " << *itB << std::endl);

                    // we replace the RayUnion with its lower ray, so all intersections pointers should still be valid
                    // all intersections processed after this point will be with lower ray
                    itB->foldExtension();

                    LOG(idx << " new ray: " << *itB << std::endl);

                    ASSERT(!itB->isExtended());

                    break;
                }
            }

            ++idx;

            LOG(std::endl);
        }

#ifdef WITH_CAIRO
        basePainter.save("img_k" + std::to_string(k));
#endif
    }
};
