//
// Created by funke on 10/1/21.
//

#pragma once

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
    using CGAL_Box = CGAL::Bbox_2<K>;

    CGAL::Ray_2<K> iRay;

    tIndex leftRegion;
    tIndex rightRegion;

    // extension of ray with another ray
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

    bool leftOf(const tFloatVector &x) const {
        // we only consider main ray, when the starting point of lower ray is
        // swept, this ray will be replaced by it

        return (((p[X] + std::cos(angle)) - p[X]) * (x[Y] - p[Y]) - ((p[Y] + std::sin(angle)) - p[Y]) * (x[X] - p[X])) > 0;
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
            painter.drawLine(p, {p[X] + std::cos(angle), p[Y] + std::sin(angle)});
        } else {
            painter.drawLine(p, ext->p);
            painter.drawLine(ext->p, {ext->p[X] + std::cos(ext->angle), ext->p[Y] + std::sin(ext->angle)});
        }
    }

#endif

private:
    static tIntersectionRetVal isRR(const CGAL_Ray &a, const CGAL_Ray &b, const CGAL_Box &bounds) {
        if (std::fmod(std::fmod(a.angle - b.angle, M_PI) + M_PI, M_PI) == 0) {
            return {false, {}};
        }

        if (std::fmod(std::fmod(a.angle, M_PI) + M_PI, M_PI) == M_PI_2) {
            // vertical line this's x
            return {true, {a.p[X], b.tanAngle * (a.p[X] - b.p[X]) + b.p[Y]}};
        } else if (std::fmod(std::fmod(b.angle, M_PI) + M_PI, M_PI) == M_PI_2) {
            // vertical line at b's x
            return {true, {b.p[X], a.tanAngle * (b.p[X] - a.p[X]) + a.p[Y]}};
        }

        tFloat m0 = a.tanAngle;// Line 0: y = m0 (x - x0) + y0
        tFloat m1 = b.tanAngle;// Line 1: y = m1 (x - x1) + y1

        tFloat x = ((m0 * a.p[X] - m1 * b.p[X]) - (a.p[Y] - b.p[Y])) / (m0 - m1);
        tFloat y = m0 * (x - a.p[X]) + a.p[Y];

        if (!bounds.contains({x, y})) {
            return {false, {}};
        }

        return {true, {x, y}};
    }

    static tIntersectionRetVal isUR(const CGAL_Ray &ua, const CGAL_Ray &b, const CGAL_Box &bounds) {
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

    static tIntersectionRetVal isRU(const CGAL_Ray &a, const CGAL_Ray &ub, const CGAL_Box &bounds) {
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

    static tIntersectionRetVal isUU(const CGAL_Ray &ua, const CGAL_Ray &ub, const CGAL_Box &bounds) {
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

std::ostream &operator<<(std::ostream &os, const CGAL_Ray &r) {
    os << "p: " << r.p << " angle: " << r.angle << " boundary: "
       << (r.leftRegion != tIndex(-1) ? std::to_string(r.leftRegion) : "INF") << "/"
       << (r.rightRegion != tIndex(-1) ? std::to_string(r.rightRegion) : "INF");
    if (r.ext) {
        os << " extension: " << *r.ext;
    }

    return os;
}

std::string to_string(const CGAL_Ray &r) {
    std::stringstream os;
    os << "p: " << r.p << " angle: " << r.angle << " boundary: "
       << (r.leftRegion != tIndex(-1) ? std::to_string(r.leftRegion) : "INF") << "/"
       << (r.rightRegion != tIndex(-1) ? std::to_string(r.rightRegion) : "INF");
    if (r.ext) {
        os << " extension: " << *r.ext;
    }

    return os.str();
}