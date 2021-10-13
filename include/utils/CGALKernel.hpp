//
// Created by funke on 10/1/21.
//

#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Ray_2.h>

#include "Predicates.hpp"
#include "Types.hpp"

#ifdef WITH_CAIRO
#include "Painter.hpp"
#endif

struct CGAL_Ray {

    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using CGAL_Point = K::Point_2;
    using CGAL_Angle = CGAL::Direction_2<K>;
    using CGAL_Segment = CGAL::Segment_2<K>;
    using CGAL_Box = CGAL::Bbox_2;

    CGAL::Ray_2<K> iRay;

    tIndex leftRegion;
    tIndex rightRegion;

    // extension of ray with a segment before ray
    std::unique_ptr<CGAL_Segment> ext;// initialized with nullptr

    CGAL_Ray(const CGAL_Point &p_, const CGAL_Angle &angle_,
             const tIndex &lr, const tIndex &rr)
        : iRay(p_, angle_), leftRegion(lr), rightRegion(rr) {}

    CGAL_Ray(const CGAL_Ray &o)
        : iRay(o.iRay), leftRegion(o.leftRegion),
          rightRegion(o.rightRegion) {
        if (o.ext) {
            ext = std::make_unique<CGAL_Segment>(*o.ext);
        }
    }

    CGAL_Ray &operator=(const CGAL_Ray &o) {
        iRay = o.iRay;
        leftRegion = o.leftRegion;
        rightRegion = o.rightRegion;

        if (o.ext) {
            ext = std::make_unique<CGAL_Segment>(*o.ext);
        }

        return *this;
    }

    bool leftOf(const CGAL_Point &x) const {
        // we only consider main ray, when the starting point of lower ray is
        // swept, this ray will be replaced by it

        return iRay.supporting_line().oriented_side(x) == CGAL::ON_NEGATIVE_SIDE;
    }

    struct tIntersectionRetVal {
        bool valid;
        tFloatVector pos;
    };

    tIntersectionRetVal intersection(const CGAL_Ray &b, const CGAL_Box &bounds) const {
        return intersection(*this, b, bounds);
    }

public:
    static tIntersectionRetVal intersection(const CGAL_Ray &a, const CGAL_Ray &b, const CGAL_Box &bounds) {
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
                    painter.drawLine(iRay.start(), iRay.point(1));
                } else {
                    painter.drawLine(ext->start(), ext->end());
                    painter.drawLine(iRay.start(), iRay.point(1));
                }
    }

#endif

private:
    static tIntersectionRetVal isRR(const CGAL_Ray &a, const CGAL_Ray &b, const CGAL_Box &bounds) {

        auto result = CGAL::intersection(a.iRay, b.iRay);

        if (result) {
            const CGAL_Point *p = boost::get<CGAL_Point>(&*result);
            if (p && CGAL::do_intersect(*p, bounds)) {
                return {true, {p->x(), p->y()}};
            } else {
                // could also be a Segment_2 or a Ray_2 TODO: handle correctly
                return {false, {}};
            }
        } else {
            return {false, {}};
        }
    }

    static tIntersectionRetVal isUR(const CGAL_Ray &ua, const CGAL_Ray &b, const CGAL_Box &bounds) {
        assert(ua.ext);

        // check for intersection of pre-segment of ua and b
        auto result = CGAL::intersection(*ua.ext, b.iRay);

        if (result) {
            if (const CGAL_Point *p = boost::get<CGAL_Point>(&*result)) {
                assert(CGAL::do_intersect(*p, bounds));
                return {true, {p->x(), p->y()}};
            } else {
                // could also be a Segment_2 or a Ray_2 TODO: handle correctly
                return {false, {}};
            }
        } else {
            // check intersecton of main ray of ua and b
            return isRR(ua, b, bounds);
        }
    }

    static tIntersectionRetVal isRU(const CGAL_Ray &a, const CGAL_Ray &ub, const CGAL_Box &bounds) {
        assert(ub.ext);

        // check for intersection of pre-segment of ua and b
        auto result = CGAL::intersection(a.iRay, *ub.ext);

        if (result) {
            if (const CGAL_Point *p = boost::get<CGAL_Point>(&*result)) {
                assert(CGAL::do_intersect(*p, bounds));
                return {true, {p->x(), p->y()}};
            } else {
                // could also be a Segment_2 or a Ray_2 TODO: handle correctly
                return {false, {}};
            }
        } else {
            // check intersecton of a and main ray of ub
            return isRR(a, ub, bounds);
        }
    }

    static tIntersectionRetVal isUU(const CGAL_Ray &ua, const CGAL_Ray &ub, const CGAL_Box &bounds) {
        assert(ua.ext && ub.ext);

        // check for intersection of pre-segment of ua and ub
        auto result = CGAL::intersection(*ua.ext, *ub.ext);
        if (result) {
            if (const CGAL_Point *p = boost::get<CGAL_Point>(&*result)) {
                assert(CGAL::do_intersect(*p, bounds));
                return {true, {p->x(), p->y()}};
            }
        }


        // check for intersection of main ray of ua and pre-segment of ub
        result = CGAL::intersection(ua.iRay, *ub.ext);
        if (result) {
            if (const CGAL_Point *p = boost::get<CGAL_Point>(&*result)) {
                assert(CGAL::do_intersect(*p, bounds));
                return {true, {p->x(), p->y()}};
            }
        }

        // check for intersection of pre-segment of ua and main ray of ub
        result = CGAL::intersection(*ua.ext, ub.iRay);
        if (result) {
            if (const CGAL_Point *p = boost::get<CGAL_Point>(&*result)) {
                assert(CGAL::do_intersect(*p, bounds));
                return {true, {p->x(), p->y()}};
            }
        }

        // check for intersection of main rays of ua and ub
        return isRR(ua, ub, bounds);
    }
};

std::ostream &operator<<(std::ostream &os, const CGAL_Ray &r) {
    os << "boundary: "
       << (r.leftRegion != tIndex(-1) ? std::to_string(r.leftRegion) : "INF") << "/"
       << (r.rightRegion != tIndex(-1) ? std::to_string(r.rightRegion) : "INF");
    if (r.ext) {
        os << " extension: " << *r.ext;
    }

    return os;
}

std::string to_string(const CGAL_Ray &r) {
    std::stringstream os;
    os << "boundary: "
       << (r.leftRegion != tIndex(-1) ? std::to_string(r.leftRegion) : "INF") << "/"
       << (r.rightRegion != tIndex(-1) ? std::to_string(r.rightRegion) : "INF");
    if (r.ext) {
        os << " extension: " << *r.ext;
    }

    return os.str();
}