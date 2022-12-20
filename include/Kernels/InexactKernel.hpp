//
// Created by funke on 10/1/21.
//

#pragma once

#include <memory>
#include <optional>
#include <sstream>
#include <variant>

#include "Predicates.hpp"
#include "Types.hpp"

#ifdef WITH_CAIRO
#include "Painter.hpp"
#endif

class InexactKernel {

public:
    using Float = tIFloat;
    using Vector = tIFloatVector;
    using Point = tIFloatVector;

    static std::string name() {
        return "InexactKernel";
    }

    class Direction {
    public:
        Direction(const tIFloat _dir) : dir(_dir), tanDir(std::tan(dir)), vec({std::cos(dir), std::sin(dir)}) {}
        Direction(const Vector _vec) : Direction(atan2P(_vec[Y], _vec[X])) {}

        Float prj(const Point &p) {
            return (p[X] * vec[X] + p[Y] * vec[Y]) / dot(vec, vec);
        }

        Float tan() const {
            return tanDir;
        }

        Float cos() const {
            return std::cos(dir);
        }

        Float sin() const {
            return std::sin(dir);
        }

        tIFloat angle() const {
            return dir;
        }

        Direction operator-() const {
            return Direction(wrapAngle(dir + M_PI));
        }

        Direction perp() const {
            return {wrapAngle(dir + M_PI_2)};// TODO: check angle orientation
        }

        Direction perp(const Direction &ref) const {
            auto p = perp();
            return angleBetween(p.dir, ref.dir) <= M_PI_2 ? p : Direction(wrapAngle(p.dir + M_PI));
        }

        bool operator==(const Direction &o) const {
            return dir == o.dir;
        }

    private:
        Float dir;
        Float tanDir;
        Vector vec;
    };

    class Line {
    public:
        Line(const Point &_a, const Point &_b) : a(_a), b(_b) {}

    private:
        Point a;
        Point b;
    };

    struct Ray {
    private:
        Point p;
        Direction dir;

    public:
        tIndex leftRegion;
        tIndex rightRegion;

    private:
        // extension of ray with another ray
        std::unique_ptr<Ray> ext;// initialized with nullptr

    public:
        Ray(const Point &p_, const Direction &angle_, const tIndex &lr, const tIndex &rr)
            : p(p_), dir(angle_), leftRegion(lr), rightRegion(rr) {}

        Ray(const Point &p_, const Direction &angle_, const tIndex &lr,
            const tIndex &rr, const Ray &ext_)
            : p(p_), dir(angle_), leftRegion(lr), rightRegion(rr),
              ext(std::make_unique<Ray>(ext_)) {

            // the starting points of upper and lower ray should be different
            ASSERT(!approxEQ(p, ext->p));

            // the starting point of the lower ray must be on the upper rayl;
            ASSERT(approxEQ(ext->p, {ext->p[X], dir.tan() * (ext->p[X] - p[X]) + p[Y]}));
            ASSERT(leftRegion == ext->leftRegion);
            ASSERT(rightRegion == ext->rightRegion);
        }

        Ray(const Ray &o)
            : p(o.p), dir(o.dir), leftRegion(o.leftRegion),
              rightRegion(o.rightRegion) {
            if (o.ext) {
                ext = std::make_unique<Ray>(*o.ext);
            }
        }

        std::string str() const {
            std::stringstream os;

            // use fabs(angle()) to eliminate -0 output

            os << (leftRegion != INF_IDX ? std::to_string(leftRegion) : "INF_IDX") << "/"
               << (rightRegion != INF_IDX ? std::to_string(rightRegion) : "INF_IDX")
               << " p: " << p << " a: " << std::fabs(dir.angle());
            if (ext) {
                os << " EXT: "
                   << " p: " << ext->p << " a: " << std::fabs(ext->dir.angle());
                ;
            }

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
            p = o.p;
            dir = o.dir;
            leftRegion = o.leftRegion;
            rightRegion = o.rightRegion;

            if (o.ext) {
                ext = std::make_unique<Ray>(*o.ext);
            }

            return *this;
        }

        // convenience constructor for ray unions
        Ray(const Ray &upper, const Ray &lower) : Ray(upper.p, upper.dir, upper.leftRegion, upper.rightRegion, lower) {}

        bool isExtended() const {
            return (bool) ext;
        }

        Point origin() const {
            return p;
        }

        Direction direction() const {
            return dir;
        }

        Point extOrigin() const {
            ASSERT(ext);
            return ext->p;
        }

        void setOrigin(const Point &p_) {
            p = p_;
        }

        void foldExtension() {
            ASSERT(ext);

            // we replace the RayUnion with its lower ray, so all intersections pointers should still be valid
            // all intersections processed after this point will be with lower ray
            p = ext->p;
            dir = ext->dir;
            // and delete the lower ray
            ext.reset();

            ASSERT(!ext);
        }

        tOrientedSide orientedSide(const Point &x) const {
            // we only consider main ray, when the starting point of lower ray is
            // swept, this ray will be replaced by it
            Float res = (((p[X] + dir.cos()) - p[X]) * (x[Y] - p[Y]) - ((p[Y] + dir.sin()) - p[Y]) * (x[X] - p[X]));
            return res > 0 ? tOrientedSide::LEFT : res < 0 ? tOrientedSide::RIGHT
                                                           : tOrientedSide::LINE;
        }

        using tIntersectionRetVal = std::optional<std::variant<Point, Ray>>;

        tIntersectionRetVal intersection(const Ray &b, const tBox &bounds) const {
            return intersection(*this, b, bounds);
        }

    public:
        static tIntersectionRetVal intersection(const Ray &a, const Ray &b, const tBox &bounds) {

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
                painter.drawLine(p, {p[X] + dir.cos(), p[Y] + dir.sin()});
            } else {
                painter.drawLine(p, ext->p);
                painter.drawLine(ext->p, {ext->p[X] + ext->dir.cos(), ext->p[Y] + ext->dir.sin()});
            }
        }

#endif

    private:
        static tIntersectionRetVal isRR(const Ray &a, const Ray &b, const tBox &bounds) {
            if (std::fmod(std::fmod(a.dir.angle() - b.dir.angle(), M_PI) + M_PI, M_PI) == 0) {
                // parallel lines
                //TODO compute ray or segment
                return std::nullopt;
            }

            if (std::fmod(std::fmod(a.dir.angle(), M_PI) + M_PI, M_PI) == M_PI_2) {
                // vertical line this's x
                return Point{a.p[X], b.dir.tan() * (a.p[X] - b.p[X]) + b.p[Y]};
            } else if (std::fmod(std::fmod(b.dir.angle(), M_PI) + M_PI, M_PI) == M_PI_2) {
                // vertical line at b's x
                return Point{b.p[X], a.dir.tan() * (b.p[X] - a.p[X]) + a.p[Y]};
            }

            Float m0 = a.dir.tan();// Line 0: y = m0 (x - x0) + y0
            Float m1 = b.dir.tan();// Line 1: y = m1 (x - x1) + y1

            Float x = ((m0 * a.p[X] - m1 * b.p[X]) - (a.p[Y] - b.p[Y])) / (m0 - m1);
            Float y = m0 * (x - a.p[X]) + a.p[Y];

            if (!bounds.contains({x, y})) {
                return std::nullopt;
            }

            return Point{x, y};
        }

        static tIntersectionRetVal isUR(const Ray &ua, const Ray &b, const tBox &bounds) {
            ASSERT(ua.ext);

            // check for intersection of main ray of ua and b
            auto is = isRR(ua, b, bounds);
            // check whether IS is before starting point of extension ray of ur
            if (is) {
                const Point *p = std::get_if<Point>(&*is);
                if (p && distance2(ua.p, *p) < distance2(ua.p, ua.ext->p)) {
                    return is;
                }
                //TODO handle ray return
            }

            // upper ray of ur does not intersect r OR intersection is below starting point of extension ray
            is = isRR(*ua.ext, b, bounds);
            // check whether IS is after starting point of extension ray of ur
            if (is) {
                const Point *p = std::get_if<Point>(&*is);
                if (p && distance2(ua.p, *p) >= distance2(ua.p, ua.ext->p)) {
                    return is;
                }
            }

            return std::nullopt;
        }

        static tIntersectionRetVal isRU(const Ray &a, const Ray &ub, const tBox &bounds) {
            ASSERT(ub.ext);

            // check for intersection of a and main ray of ub
            auto is = isRR(a, ub, bounds);
            // check whether IS is before starting point of extension ray of ur
            if (is) {
                const Point *p = std::get_if<Point>(&*is);
                if (p && distance2(ub.p, *p) < distance2(ub.p, ub.ext->p)) {
                    return is;
                }
            }

            // upper ray of ur does not intersect r OR intersection is below starting point of extension ray
            is = isRR(a, *ub.ext, bounds);
            // check whether IS is after starting point of extension ray of ur
            if (is) {
                const Point *p = std::get_if<Point>(&*is);
                if (p && distance2(ub.p, *p) >= distance2(ub.p, ub.ext->p)) {
                    return is;
                }
            }

            return std::nullopt;
        }

        static tIntersectionRetVal isUU(const Ray &ua, const Ray &ub, const tBox &bounds) {
            ASSERT(ua.ext && ub.ext);

            // check for intersection of main ray of ua and ub
            auto is = isRR(ua, ub, bounds);
            // check whether IS is before starting point of lowerRay of ua
            // check whether IS is before starting point of lowerRay of ub
            if (is) {
                const Point *p = std::get_if<Point>(&*is);
                if (p &&
                    distance2(ua.p, *p) < distance2(ua.p, ua.ext->p) &&
                    distance2(ub.p, *p) < distance2(ub.p, ub.ext->p)) {

                    return is;
                }
            }

            // check for intersection of main ray of ua and lower ray of ub
            is = isRR(ua, *ub.ext, bounds);
            // check whether IS is before starting point of lowerRay of ua
            // check whether IS is after starting point of lowerRay of ub
            if (is) {
                const Point *p = std::get_if<Point>(&*is);
                if (p &&
                    distance2(ua.p, *p) < distance2(ua.p, ua.ext->p) &&
                    distance2(ub.p, *p) >= distance2(ub.p, ub.ext->p)) {
                    return is;
                }
            }

            // check for intersection of lower ray of ua and main ray of ub
            is = isRR(*ua.ext, ub, bounds);
            // check whether IS is after starting point of ua's lowerRay
            // check whether IS is before starting point of ub's lowerRay
            if (is) {
                const Point *p = std::get_if<Point>(&*is);
                if (p &&
                    distance2(ua.p, *p) >= distance2(ua.p, ua.ext->p) &&
                    distance2(ub.p, *p) < distance2(ub.p, ub.ext->p)) {
                    return is;
                }
            }

            // check for intersection of lower rays of ua and ub
            is = isRR(*ua.ext, *ub.ext, bounds);
            // check whether IS is after starting point of ua's lowerRay
            // check whether IS is after starting point of ub's lowerRay
            if (is) {
                const Point *p = std::get_if<Point>(&*is);
                if (p &&
                    distance2(ua.p, *p) >= distance2(ua.p, ua.ext->p) &&
                    distance2(ub.p, *p) >= distance2(ub.p, ub.ext->p)) {
                    return is;
                }
            }

            return std::nullopt;
        }
    };

    static Point mkPoint(const tIFloatVector &p) {
        return p;
    }

    static std::vector<Point> mkPoints(const std::vector<tIFloatVector> &p) {
        return p;
    }

    static tBox mkBBox(const tBox &box) {
        return box;
    }

    static Point Midpoint(const Point &a, const Point &b) {
        return 0.5 * (a + b);
    }

    static Direction Bisector(const Direction &l, const Direction &r) {
        auto dir = Direction(wrapAngle(0.5 * (l.angle() + r.angle())));
        return l.angle() < dir.angle() && dir.angle() < r.angle() ? dir : -dir;
    }

    static Direction Bisector(const Point &l, const Point &r) {
        return Direction(r - l).perp();
    }

    static Direction Bisector(const Point &l, const Point &r, const Direction &ref) {
        return Direction(r - l).perp(ref);
    }

    static Float distance2(const Point &a, const Point &b) {
        Float dist2 = 0;
        for (tDim d = 0; d < a.size(); ++d) {
            dist2 += (b[d] - a[d]) * (b[d] - a[d]);
        }
        return dist2;
    }

    //    static Float distance(const Point &a, const Point &b) {
    //        return std::sqrt(distance2(a, b));
    //    }

    static bool approxEQ(const Point &a, const Point &b) {
        return distance2(a, b) < MaxError<Float>::value;//TODO: use more meaningful test
    }

    static bool approxEQ(const Float &a, const Float &b) {
        return std::fabs(a - b) < MaxError<Float>::value;//TODO: use more meaningful test
    }

    static bool approxLT(const Float &a, const Float &b) {
        return a - b < MaxError<Float>::value;//TODO: use more meaningful test
    }

    static bool approxGT(const Float &a, const Float &b) {
        return a - b > MaxError<Float>::value;//TODO: use more meaningful test
    }

    static tIFloat to_float(const Float &x) {
        return x;
    }

    static tIFloat to_float_exact(const Float &x) {
        return x;
    }

    static bool compareDistance(const Point &origin, const Point &newPoint, [[maybe_unused]] const Point &oldPoint, const Float &oldDist) {
        return compareDistance(origin, newPoint, oldDist);
    }

    static bool compareDistance(const Point &origin, const Point &newPoint, const Point &oldPoint) {
        return distance2(origin, newPoint) < distance2(origin, oldPoint);
    }

    static bool compareDistance(const Point &origin, const Point &newPoint, const Float &oldDist) {
        return distance2(origin, newPoint) < oldDist;
    }

    static auto computeCones(const tDim &cones) {
        std::vector<Direction> rays;

        for (tDim k = 0; k < cones; ++k) {
            rays.emplace_back(k * (2 * M_PI / cones));
        }

        return rays;
    }

    static auto computePointCones([[maybe_unused]] const Point &p, const std::vector<Direction> &rays) {
        return rays;
    }

    static tDim getCone(const Point &origin, const Point &newPoint, const std::vector<Direction> &cones) {
        return std::floor(atan2P(to_float(newPoint[1] - origin[1]), to_float(newPoint[0] - origin[0])) / (2 * M_PI / cones.size()));
    }
};

std::ostream &operator<<(std::ostream &os, const InexactKernel::Ray &r) {
    os << r.str();
    return os;
}

std::string to_string(const InexactKernel::Ray &r) {
    return r.str();
}
