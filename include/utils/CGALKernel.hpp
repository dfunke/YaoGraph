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

class CGALKernel {

public:
    using K = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Float = tFloat;
    using Vector = CGAL::Vector_2<K>;
    using Point = K::Point_2;
    using Segment = CGAL::Segment_2<K>;
    using Box = CGAL::Bbox_2;

    class Direction {
    public:
        Direction(const tFloat _dir) : dir(Vector(std::cos(_dir), std::sin(_dir))) {}

        Float prj(const Point &p) {
            return 0;
        }

        auto &base() const {
            return dir;
        }

    private:
        CGAL::Direction_2<K> dir;
    };

    struct Ray {
    private:
        CGAL::Ray_2<K> iRay;


    public:
        tIndex leftRegion;
        tIndex rightRegion;

    private:
        // extension of ray with a segment before ray
        std::unique_ptr<Segment> ext;// initialized with nullptr

    public:
        Ray(const Point &p_, const Direction &angle_,
            const tIndex &lr, const tIndex &rr)
            : iRay(p_, angle_.base()), leftRegion(lr), rightRegion(rr) {}

        Ray(const Ray &o)
            : iRay(o.iRay), leftRegion(o.leftRegion),
              rightRegion(o.rightRegion) {
            if (o.ext) {
                ext = std::make_unique<Segment>(*o.ext);
            }
        }

        std::string str() const {
            std::stringstream os;

            os << (leftRegion != tIndex(-1) ? std::to_string(leftRegion) : "INF") << "/"
               << (rightRegion != tIndex(-1) ? std::to_string(rightRegion) : "INF");
            if (ext) {
                os << " EXT";
            }

            return os.str();
        }

        Ray &operator=(const Ray &o) {
            iRay = o.iRay;
            leftRegion = o.leftRegion;
            rightRegion = o.rightRegion;

            if (o.ext) {
                ext = std::make_unique<Segment>(*o.ext);
            }

            return *this;
        }

        bool leftOf(const Point &x) const {
            // extension is before ray, as long as its active we check against it
            if (ext) {
                return ext->supporting_line().oriented_side(x) == CGAL::ON_NEGATIVE_SIDE;
            } else {
                return iRay.supporting_line().oriented_side(x) == CGAL::ON_NEGATIVE_SIDE;
            }
        }

        bool isExtended() const {
            return (bool) ext;
        }

        Point origin() const {
            if (ext) {
                return ext->start();
            } else {
                return iRay.start();
            }
        }

        Point extOrigin() const {
            assert(ext);
            return iRay.start();
        }

        void setOrigin(const Point &p_) {
            iRay = CGAL::Ray_2<K>(p_, iRay.direction());
        }

        void foldExtension() {
            assert(ext);

            // delete the upper ray
            ext.reset();

            assert(!ext);
        }

        struct tIntersectionRetVal {
            bool valid;
            tFloatVector pos;
        };

        tIntersectionRetVal intersection(const Ray &b, const Box &bounds) const {
            return intersection(*this, b, bounds);
        }

    public:
        static tIntersectionRetVal intersection(const Ray &a, const Ray &b, const Box &bounds) {

            if (a.iRay.start() == b.iRay.start()) {
                // if both rays originate at the same point, we don't consider it an intersection
                return {false, {}};
            } else {
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
        static tIntersectionRetVal isRR(const Ray &a, const Ray &b, const Box &bounds) {

            auto result = CGAL::intersection(a.iRay, b.iRay);

            if (result) {
                const Point *p = boost::get<Point>(&*result);
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

        static tIntersectionRetVal isUR(const Ray &ua, const Ray &b, const Box &bounds) {
            assert(ua.ext);

            // check for intersection of pre-segment of ua and b
            auto result = CGAL::intersection(*ua.ext, b.iRay);

            if (result) {
                if (const Point *p = boost::get<Point>(&*result)) {
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

        static tIntersectionRetVal isRU(const Ray &a, const Ray &ub, const Box &bounds) {
            assert(ub.ext);

            // check for intersection of pre-segment of ua and b
            auto result = CGAL::intersection(a.iRay, *ub.ext);

            if (result) {
                if (const Point *p = boost::get<Point>(&*result)) {
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

        static tIntersectionRetVal isUU(const Ray &ua, const Ray &ub, const Box &bounds) {
            assert(ua.ext && ub.ext);

            // check for intersection of pre-segment of ua and ub
            auto result = CGAL::intersection(*ua.ext, *ub.ext);
            if (result) {
                if (const Point *p = boost::get<Point>(&*result)) {
                    assert(CGAL::do_intersect(*p, bounds));
                    return {true, {p->x(), p->y()}};
                }
            }


            // check for intersection of main ray of ua and pre-segment of ub
            result = CGAL::intersection(ua.iRay, *ub.ext);
            if (result) {
                if (const Point *p = boost::get<Point>(&*result)) {
                    assert(CGAL::do_intersect(*p, bounds));
                    return {true, {p->x(), p->y()}};
                }
            }

            // check for intersection of pre-segment of ua and main ray of ub
            result = CGAL::intersection(*ua.ext, ub.iRay);
            if (result) {
                if (const Point *p = boost::get<Point>(&*result)) {
                    assert(CGAL::do_intersect(*p, bounds));
                    return {true, {p->x(), p->y()}};
                }
            }

            // check for intersection of main rays of ua and ub
            return isRR(ua, ub, bounds);
        }
    };
};

std::ostream &operator<<(std::ostream &os, const CGALKernel::Ray &r) {
    os << r.str();
    return os;
}
