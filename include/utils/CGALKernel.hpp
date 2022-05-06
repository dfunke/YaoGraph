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

using ExactPredicatesInexactConstructions = CGAL::Exact_predicates_inexact_constructions_kernel;
using ExactPredicatesExactConstructions = CGAL::Exact_predicates_exact_constructions_kernel;

template<typename K>
class CGALKernel {

public:
    using Float = typename K::FT;
    using Vector = CGAL::Vector_2<K>;
    using Point = typename K::Point_2;
    using Segment = CGAL::Segment_2<K>;
    using Box = CGAL::Bbox_2;

    class Direction {
    public:
        Direction(const tFloat _dir) : dir(Vector(std::cos(_dir), std::sin(_dir))) {}
        Direction(const Vector &v) : dir(v) {}
        Direction(const CGAL::Direction_2<K> &d) : dir(d) {}

        Float prj(const Point &p) {
            return CGAL::scalar_product(Vector(p[X], p[Y]), dir.vector()) / CGAL::scalar_product(dir.vector(), dir.vector());
        }

        tFloat angle() const {
            return atan2P(CGAL::to_double(dir.dy()), CGAL::to_double(dir.dx()));
        }

        auto &base() const {
            return dir;
        }

        Direction perp() const {
            return {dir.perpendicular(CGAL::CLOCKWISE)};// TODO: check angle orientation
        }

        Direction perp(const Direction &ref) const {
            auto p = perp().dir.vector();
            if (CGAL::angle(p, ref.dir.vector()) == CGAL::OBTUSE) {
                p *= -1;
            }
            return Direction(p);
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

            // use fabs(atan2P) to eliminate -0 output

            if (ext) {
                os << " p: " << ext->start() << " a: " << std::fabs(atan2P(CGAL::to_double(ext->direction().dy()), CGAL::to_double(ext->direction().dx())));
                os << " EXT: ";
            }
            os << " p: " << iRay.start() << " a: " << std::fabs(atan2P(CGAL::to_double(iRay.direction().dy()), CGAL::to_double(iRay.direction().dx())));

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

        // convenience constructor for ray unions
        Ray(const Ray &upper, const Ray &lower) {
            *this = lower;
            this->ext = std::make_unique<Segment>(upper.origin(), lower.origin());
        }


        bool leftOf(const Point &x) const {
            // extension is before ray, as long as its active we check against it
            if (ext) {
                return ext->supporting_line().oriented_side(x) == CGAL::ON_POSITIVE_SIDE;
            } else {
                return iRay.supporting_line().oriented_side(x) == CGAL::ON_POSITIVE_SIDE;
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
            Point pos;
        };

        tIntersectionRetVal intersection(const Ray &b, const Box &bounds) const {
            return intersection(*this, b, bounds);
        }

    public:
        static tIntersectionRetVal intersection(const Ray &a, const Ray &b, const Box &bounds) {

            if (a.origin() == b.origin()) {
                // we ignore rays originating at he same source
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
                if (p) {
                    if (CGAL::do_intersect(*p, bounds)) {
                        return {true, {p->x(), p->y()}};
                    } else {
                        return {false, {}};
                    }
                } else {
                    // could also be a Segment_2 or a Ray_2 TODO: handle correctly
                    assert(false);
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

    static Point mkPoint(const tFloatVector &p) {
        return Point(p[X], p[Y]);
    }

    static Box mkBBox(const tBox &box) {
        return Box(box.low[X], box.low[Y], box.high[X], box.high[Y]);
    }

    static Point Midpoint(const Point &a, const Point &b) {
        return CGAL::midpoint(a, b);
    }

    static Direction Bisector(const Point &l, const Point &r) {
        return Direction(CGAL::bisector(l, r).direction());
    }

    static Direction Bisector(const Point &l, const Point &r, const Direction &ref) {
        auto b = CGAL::bisector(l, r).direction().vector();
        if (CGAL::angle(b, ref.base().vector()) == CGAL::OBTUSE) {
            b *= -1;
        }
        return Direction(b);
    }

    static Float distance2(const Point &a, const Point &b) {
        return CGAL::squared_distance(a, b);
    }

    //    static Float distance(const Point &a, const Point &b) {
    //        return CGAL::sqrt(distance2(a, b));
    //    }

    static bool approxEQ(const Point &a, const Point &b) {
        return distance2(a, b) < MaxError<tFloat>::value;//TODO: use more meaningful test
    }

    static bool approxEQ(const Float &a, const Float &b) {
        return std::fabs(CGAL::to_double(a) - CGAL::to_double(b)) < MaxError<tFloat>::value;//TODO: use more meaningful test
    }

    static bool approxLT(const Float &a, const Float &b) {
        return a - b < 0;//TODO: use more meaningful test
    }

    static bool approxGT(const Float &a, const Float &b) {
        return a - b > 0;//TODO: use more meaningful test
    }

    static tFloat to_float(const Float &x) {
        return CGAL::to_double(x);
    }
};

//TODO: there has to be a better way for template deduction
std::ostream &operator<<(std::ostream &os, const typename CGALKernel<ExactPredicatesInexactConstructions>::Ray &r) {
    os << r.str();
    return os;
}
std::string to_string(const typename CGALKernel<ExactPredicatesInexactConstructions>::Ray &r) {
    return r.str();
}

std::ostream &operator<<(std::ostream &os, const typename CGALKernel<ExactPredicatesExactConstructions>::Ray &r) {
    os << r.str();
    return os;
}
std::string to_string(const typename CGALKernel<ExactPredicatesExactConstructions>::Ray &r) {
    return r.str();
}
