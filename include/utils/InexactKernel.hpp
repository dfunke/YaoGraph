//
// Created by funke on 10/1/21.
//

#pragma once

#include "Predicates.hpp"
#include "Types.hpp"

#ifdef WITH_CAIRO
#include "Painter.hpp"
#endif

class InexactKernel {

public:
    using Float = tFloat;
    using Vector = tFloatVector;
    using Point = tFloatVector;

    class Direction {
    public:
        Direction(const Float _dir) : dir(_dir), tanDir(std::tan(dir)), vec({std::cos(dir), std::sin(dir)}) {}
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

        Float angle() const {
            return dir;
        }

        Direction perp() const {
            return {wrapAngle(dir + M_PI_2)};// TODO: check angle orientation
        }

        Direction perp(const Direction &ref) const {
            auto p = perp();
            return angleBetween(p.dir, ref.dir) <= M_PI_2 ? p : Direction(wrapAngle(p.dir + M_PI));
        }

    private:
        Float dir;
        Float tanDir;
        Vector vec;
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
            assert(!approxEQ(p, ext->p));

            // the starting point of the lower ray must be on the upper rayl;
            assert(approxEQ(ext->p, {ext->p[X], dir.tan() * (ext->p[X] - p[X]) + p[Y]}));
            assert(leftRegion == ext->leftRegion);
            assert(rightRegion == ext->rightRegion);
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

            os << (leftRegion != tIndex(-1) ? std::to_string(leftRegion) : "INF") << "/"
               << (rightRegion != tIndex(-1) ? std::to_string(rightRegion) : "INF")
               << " p: " << p << " a: " << dir.angle();
            if (ext) {
                os << " EXT: "
                   << " p: " << ext->p << " a: " << ext->dir.angle();
                ;
            }

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

        Point extOrigin() const {
            assert(ext);
            return ext->p;
        }

        void setOrigin(const Point &p_) {
            p = p_;
        }

        void foldExtension() {
            assert(ext);

            // we replace the RayUnion with its lower ray, so all intersections pointers should still be valid
            // all intersections processed after this point will be with lower ray
            p = ext->p;
            dir = ext->dir;
            // and delete the lower ray
            ext.reset();

            assert(!ext);
        }

        bool leftOf(const tFloatVector &x) const {
            // we only consider main ray, when the starting point of lower ray is
            // swept, this ray will be replaced by it
            return (((p[X] + dir.cos()) - p[X]) * (x[Y] - p[Y]) - ((p[Y] + dir.sin()) - p[Y]) * (x[X] - p[X])) > 0;
        }

        struct tIntersectionRetVal {
            bool valid;
            Point pos;
        };

        tIntersectionRetVal intersection(const Ray &b, const tBox &bounds) const {
            return intersection(*this, b, bounds);
        }

    public:
        static tIntersectionRetVal intersection(const Ray &a, const Ray &b, const tBox &bounds) {

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
                return {false, {}};
            }

            if (std::fmod(std::fmod(a.dir.angle(), M_PI) + M_PI, M_PI) == M_PI_2) {
                // vertical line this's x
                return {true, {a.p[X], b.dir.tan() * (a.p[X] - b.p[X]) + b.p[Y]}};
            } else if (std::fmod(std::fmod(b.dir.angle(), M_PI) + M_PI, M_PI) == M_PI_2) {
                // vertical line at b's x
                return {true, {b.p[X], a.dir.tan() * (b.p[X] - a.p[X]) + a.p[Y]}};
            }

            tFloat m0 = a.dir.tan();// Line 0: y = m0 (x - x0) + y0
            tFloat m1 = b.dir.tan();// Line 1: y = m1 (x - x1) + y1

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

    static Point mkPoint(const tFloatVector &p) {
        return p;
    }

    static tBox mkBBox(const tBox &box) {
        return box;
    }

    static Point Midpoint(const Point &a, const Point &b) {
        return 0.5 * (a + b);
    }

    static Direction Bisector(const Point &l, const Point &r) {
        return Direction(r - l).perp();
    }

    static Direction Bisector(const Point &l, const Point &r, const Direction &ref) {
        return Direction(r - l).perp(ref);
    }

    static Float distance2(const Point &a, const Point &b) {
        tFloat dist2 = 0;
        for (tDim d = 0; d < a.size(); ++d) {
            dist2 += (b[d] - a[d]) * (b[d] - a[d]);
        }
        return dist2;
    }

    static Float distance(const Point &a, const Point &b) {
        return std::sqrt(distance2(a, b));
    }

    static bool approxEQ(const Point &a, const Point &b) {
        return distance2(a, b) < MaxError<tFloat>::value;//TODO: use more meaningful test
    }

    static bool approxLT(const Float &a, const Float &b) {
        return a - b < MaxError<tFloat>::value;//TODO: use more meaningful test
    }

    static bool approxGT(const Float &a, const Float &b) {
        return a - b > MaxError<tFloat>::value;//TODO: use more meaningful test
    }
};

std::ostream &operator<<(std::ostream &os, const InexactKernel::Ray &r) {
    os << r.str();
    return os;
}

std::string to_string(const InexactKernel::Ray &r) {
    return r.str();
}
