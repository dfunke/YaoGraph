//
// Created by funke on 10/1/21.
//

#pragma once

#include <CGAL/Compute_cone_boundaries_2.h>
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
struct KernelName;

template<>
struct KernelName<ExactPredicatesExactConstructions> {
    static std::string name() {
        return "ExactPredExactCon";
    }
};

template<>
struct KernelName<ExactPredicatesInexactConstructions> {
    static std::string name() {
        return "ExactPredInexactCon";
    }
};

template<typename K>
class CGALKernel {

public:
    using Float = typename K::FT;
    using Vector = CGAL::Vector_2<K>;
    using Point = typename K::Point_2;
    using Segment = CGAL::Segment_2<K>;
    using Box = CGAL::Bbox_2;
    using Line = CGAL::Line_2<K>;

    static std::string name() {
        return "CGAL" + KernelName<K>::name();
    }

    class Direction {
    public:
        Direction(const tIFloat _dir) : dir(Vector(std::cos(_dir), std::sin(_dir))) {}
        Direction(const Vector &v) : dir(v) {}
        Direction(const CGAL::Direction_2<K> &d) : dir(d) {}

        Float prj(const Point &p) {
            return CGAL::scalar_product(Vector(p[X], p[Y]), dir.vector()) / CGAL::scalar_product(dir.vector(), dir.vector());
        }

        tIFloat angle() const {
            return atan2P(CGAL::to_double(dir.dy()), CGAL::to_double(dir.dx()));
        }

        auto &base() const {
            return dir;
        }

        Direction operator-() const {
            return Direction(-dir);
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

        bool operator==(const Direction &o) const {
            return dir == o.dir;
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

            os << (leftRegion != INF_IDX ? std::to_string(leftRegion) : "INF_IDX") << "/"
               << (rightRegion != INF_IDX ? std::to_string(rightRegion) : "INF_IDX");

            // use fabs(atan2P) to eliminate -0 output

            if (ext) {
                os << " p: " << ext->start() << " a: " << std::fabs(atan2P(CGAL::to_double(ext->direction().dy()), CGAL::to_double(ext->direction().dx())));
                os << " EXT: ";
            }
            os << " p: " << iRay.start() << " a: " << std::fabs(atan2P(CGAL::to_double(iRay.direction().dy()), CGAL::to_double(iRay.direction().dx())));

            return os.str();
        }

        std::string sstr() const {
            std::stringstream os;

            // use fabs(angle()) to eliminate -0 output

            os << (leftRegion != INF_IDX ? std::to_string(leftRegion) : "INF_IDX") << "/"
               << (rightRegion != INF_IDX ? std::to_string(rightRegion) : "INF_IDX");

            return os.str();
        }

        Ray &operator=(const Ray &o) {
            iRay = o.iRay;
            leftRegion = o.leftRegion;
            rightRegion = o.rightRegion;

            ext.reset();
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


        tOrientedSide orientedSide(const Point &x) const {
            // extension is before ray, as long as its active we check against it
            if (ext) {
                return static_cast<tOrientedSide>(ext->supporting_line().oriented_side(x));
            } else {
                return static_cast<tOrientedSide>(iRay.supporting_line().oriented_side(x));
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

        Direction direction() const {
            if (ext) {
                return ext->direction();
            } else {
                return iRay.direction();
            }
        }

        void invert() {
            if (ext) {
                *ext = ext->opposite();
            }
            iRay = iRay.opposite();
        }


        Point extOrigin() const {
            ASSERT(ext);
            return iRay.start();
        }

        void setOrigin(const Point &p_) {
            iRay = CGAL::Ray_2<K>(p_, iRay.direction());
        }

        void foldExtension() {
            ASSERT(ext);

            // delete the upper ray
            ext.reset();

            ASSERT(!ext);
        }

        using tIntersectionRetVal = std::optional<std::variant<Point, Ray>>;

        tIntersectionRetVal intersection(const Ray &b, const Box &bounds) const {
            return intersection(*this, b, bounds);
        }

    public:
        static tIntersectionRetVal intersection(const Ray &a, const Ray &b, const Box &bounds) {

            if (a.origin() == b.origin()) {
                // we ignore rays originating at he same source
                return std::nullopt;
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
                painter.draw(iRay.start(), 0, false);
                painter.drawLine(iRay.start(), iRay.point(2));
            } else {
                painter.draw(ext->start(), 0, false);
                painter.drawLine(ext->start(), ext->end());
                painter.drawSquare(ext->end());
                painter.drawLine(iRay.start(), iRay.point(2));
            }
        }

#endif

    private:
        template<typename A, typename B>
        static tIntersectionRetVal CGALintersection(const A &a, const B &b, const Box &bounds) {

            auto result = CGAL::intersection(a, b);

            if (result) {
                const Point *p = boost::get<Point>(&*result);
                if (p) {
                    // intersection is a single point, check whether its within our bounds
                    //TODO bound check required at all?
                    if (CGAL::do_intersect(*p, bounds)) {
                        return *p;
                    } else {
                        return std::nullopt;
                    }
                } else {
                    if constexpr (std::is_same_v<A, CGAL::Ray_2<K>> && std::is_same_v<B, CGAL::Ray_2<K>>) {
                        const CGAL::Ray_2<K> *r = boost::get<CGAL::Ray_2<K>>(&*result);
                        if (r) {
                            // one roy lies on other ray
                            // direction must be the same
                            ASSERT(a.direction() == b.direction());
                            ASSERT(r->direction() == a.direction());
                            return Ray(r->source(), r->direction(), INF_IDX, INF_IDX);
                        }
                    }

                    const CGAL::Segment_2<K> *s = boost::get<CGAL::Segment_2<K>>(&*result);
                    if (s) {
                        // a) one ray lies on other ray with opposite direction, segment between origin's
                        // both inputs are Rays
                        // => this should not happen
                        // b) extension segment lies on ray
                        // one input must be Segment
                        ASSERT((std::is_same_v<A, CGAL::Segment_2<K>> || std::is_same_v<B, CGAL::Segment_2<K>>) );
                        return Ray(s->source(), s->direction(), INF_IDX, INF_IDX);
                    }
                    return std::nullopt;
                }
            } else {
                return std::nullopt;
            }
        }

        static tIntersectionRetVal isRR(const Ray &a, const Ray &b, const Box &bounds) {
            return CGALintersection(a.iRay, b.iRay, bounds);
        }

        static tIntersectionRetVal isUR(const Ray &ua, const Ray &b, const Box &bounds) {
            ASSERT(ua.ext);

            // check for intersection of pre-segment of ua and b
            auto result = CGALintersection(*ua.ext, b.iRay, bounds);

            if (result) {
                return result;
            } else {
                // check intersecton of main ray of ua and b
                return isRR(ua, b, bounds);
            }
        }

        static tIntersectionRetVal isRU(const Ray &a, const Ray &ub, const Box &bounds) {
            ASSERT(ub.ext);

            // check for intersection of pre-segment of ua and b
            auto result = CGALintersection(a.iRay, *ub.ext, bounds);

            if (result) {
                return result;
            } else {
                // check intersecton of a and main ray of ub
                return isRR(a, ub, bounds);
            }
        }

        static tIntersectionRetVal isUU(const Ray &ua, const Ray &ub, const Box &bounds) {
            ASSERT(ua.ext && ub.ext);

            // check for intersection of pre-segment of ua and ub
            auto result = CGALintersection(*ua.ext, *ub.ext, bounds);
            if (result) {
                return result;
            }

            // check for intersection of main ray of ua and pre-segment of ub
            result = CGALintersection(ua.iRay, *ub.ext, bounds);
            if (result) {
                return result;
            }

            // check for intersection of pre-segment of ua and main ray of ub
            result = CGALintersection(*ua.ext, ub.iRay, bounds);
            if (result) {
                return result;
            }

            // check for intersection of main rays of ua and ub
            return isRR(ua, ub, bounds);
        }
    };

    static Point mkPoint(const tIFloatVector &p) {
        return Point(p[X], p[Y]);
    }

    static std::vector<Point> mkPoints(const std::vector<tIFloatVector> &p) {
        std::vector<Point> points;
        points.reserve(p.size());
        for (const auto &i : p) {
            points.push_back(mkPoint(i));
        }

        return points;
    }

    static Box mkBBox(const tBox &box) {
        return Box(box.low[X], box.low[Y], box.high[X], box.high[Y]);
    }

    static Point Midpoint(const Point &a, const Point &b) {
        return CGAL::midpoint(a, b);
    }

    static Direction Bisector(const Direction &l, const Direction &r) {
        return Direction(0.5 * (l.base().vector() + r.base().vector()));
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
        return distance2(a, b) < MaxError<tIFloat>::value;//TODO: use more meaningful test
    }

    static bool approxEQ(const Float &a, const Float &b) {
        return std::fabs(CGAL::to_double(a) - CGAL::to_double(b)) < MaxError<tIFloat>::value;//TODO: use more meaningful test
    }

    static bool approxLT(const Float &a, const Float &b) {
        return a - b < 0;//TODO: use more meaningful test
    }

    static bool approxGT(const Float &a, const Float &b) {
        return a - b > 0;//TODO: use more meaningful test
    }

    static tIFloat to_float(const Float &x) {
        return CGAL::to_double(x);
    }

    static tIFloat to_float_exact(const Float &x) {
        if constexpr (std::is_same_v<K, ExactPredicatesExactConstructions>) {
            x.exact();
        }
        return CGAL::to_double(x);
    }

    static bool compareDistance(const Point &origin, const Point &newPoint, [[maybe_unused]] const Point &oldPoint, const Float &oldDist) {
        if constexpr (std::is_same_v<K, ExactPredicatesInexactConstructions>) {
            return compareDistance(origin, newPoint, oldPoint);
        } else {
            // we have exact constructions
            return compareDistance(origin, newPoint, oldDist);
        }
    }

    static bool compareDistance(const Point &origin, const Point &newPoint, const Point &oldPoint) {
        return CGAL::compare_distance_to_point(origin, newPoint, oldPoint) == CGAL::SMALLER;
    }

    static bool compareDistance(const Point &origin, const Point &newPoint, const Float &oldDist) {
        return CGAL::compare_squared_distance(origin, newPoint, oldDist) == CGAL::SMALLER;
    }

    static auto computeCones(const tDim &cones) {

        CGAL::Compute_cone_boundaries_2<K> functor;
        std::vector<CGAL::Direction_2<K>> rays(cones);
        functor(cones, CGAL::Direction_2<K>(1, 0), rays.begin());

        return rays;
    }

    static auto computePointCones(const Point &p, const std::vector<CGAL::Direction_2<K>> &rays) {
        std::vector<Line> cLines;
        cLines.reserve(rays.size());
        for (tDim k = 0; k < rays.size(); ++k) {
            cLines.emplace_back(p, rays[k]);
        }

        return cLines;
    }

    static tDim getCone(const Point &origin, const Point &newPoint, const std::vector<Line> &cones) {

        tDim sec, secGuess;
        sec = secGuess = std::floor(atan2P(to_float(newPoint[1] - origin[1]), to_float(newPoint[0] - origin[0])) / (2 * M_PI / cones.size()));

#ifndef NDEBUG
        bool found = false;
#endif

        for (; sec < cones.size() || sec % cones.size() < secGuess; ++sec) {
            auto osL = cones[sec % cones.size()].oriented_side(newPoint);
            auto osR = cones[(sec + 1) % cones.size()].oriented_side(newPoint);

            if ((osL == CGAL::ON_POSITIVE_SIDE || osL == CGAL::ON_ORIENTED_BOUNDARY) && osR == CGAL::ON_NEGATIVE_SIDE) {
#ifndef NDEBUG
                found = true;
#endif
                break;
            }
        }
        ASSERT(found);

        return sec % cones.size();
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
