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
        Direction(const tFloat _dir) : dir(_dir), vec({std::cos(dir), std::sin(dir)}) {}

        Float prj(const Point &p) {
            return (p[X] * vec[X] + p[Y] * vec[Y]) / dot(vec, vec);
        }

    private:
        Float dir;
        Vector vec;
    };

    struct Ray {
    private:
        tFloatVector p;
        tFloat angle;
        tFloat tanAngle;

    public:
        tIndex leftRegion;
        tIndex rightRegion;

    private:
        // extension of ray with another ray
        std::unique_ptr<Ray> ext;// initialized with nullptr

    public:
        Ray(const tFloatVector &p_, const tFloat &angle_, const tFloat &tanAngle_, const tIndex &lr,
            const tIndex &rr)
            : p(p_), angle(angle_), tanAngle(tanAngle_), leftRegion(lr), rightRegion(rr) {}

        Ray(const tFloatVector &p_, const tFloat &angle_, const tIndex &lr,
            const tIndex &rr)
            : Ray(p_, angle_, std::tan(angle_), lr, rr) {}

        Ray(const tFloatVector &p_, const tFloat &angle_, const tFloat &tanAngle_, const tIndex &lr,
            const tIndex &rr, const Ray &ext_)
            : p(p_), angle(angle_), tanAngle(tanAngle_), leftRegion(lr), rightRegion(rr),
              ext(std::make_unique<Ray>(ext_)) {

            // the starting points of upper and lower ray should be different
            assert(!approxEQ(p, ext->p));

            // the starting point of the lower ray must be on the upper rayl;
            assert(approxEQ(ext->p, {ext->p[X], tanAngle * (ext->p[X] - p[X]) + p[Y]}));
            assert(leftRegion == ext->leftRegion);
            assert(rightRegion == ext->rightRegion);
        }

        Ray(const tFloatVector &p_, const tFloat &angle_, const tIndex &lr,
            const tIndex &rr, const Ray &ext_)
            : Ray(p_, angle_, std::tan(angle_), lr, rr, ext_) {}

        Ray(const Ray &o)
            : p(o.p), angle(o.angle), tanAngle(o.tanAngle), leftRegion(o.leftRegion),
              rightRegion(o.rightRegion) {
            if (o.ext) {
                ext = std::make_unique<Ray>(*o.ext);
            }
        }

        std::string str() const {
            std::stringstream os;

            os << "p: " << p << " angle: " << angle << " boundary: "
               << (leftRegion != tIndex(-1) ? std::to_string(leftRegion) : "INF") << "/"
               << (rightRegion != tIndex(-1) ? std::to_string(rightRegion) : "INF");
            if (ext) {
                os << " extension: " << ext->str();
            }

            return os.str();
        }

        Ray &operator=(const Ray &o) {
            p = o.p;
            angle = o.angle;
            tanAngle = o.tanAngle;
            leftRegion = o.leftRegion;
            rightRegion = o.rightRegion;

            if (o.ext) {
                ext = std::make_unique<Ray>(*o.ext);
            }

            return *this;
        }

        // convenience constructor for ray unions
        Ray(const Ray &upper, const Ray &lower) : Ray(upper.p, upper.angle, upper.tanAngle, upper.leftRegion, upper.rightRegion, lower) {}

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
};

std::ostream &operator<<(std::ostream &os, const InexactKernel::Ray &r) {
    os << r.str();
    return os;
}
