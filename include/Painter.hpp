#pragma once

#include <cairomm/cairomm.h>

#include "Predicates.hpp"
#include "Types.hpp"

#if (CAIROMM_MAJOR_VERSION == 1 and CAIROMM_MINOR_VERSION >= 16)
#define MY_CAIRO_FORMAT Cairo::ImageSurface::Format::ARGB32
#define MY_CAIRO_SLANT Cairo::ToyFontFace::Slant::NORMAL
#define MY_CAIRO_WEIGHT Cairo::ToyFontFace::Weight::NORMAL
#else
#define MY_CAIRO_FORMAT Cairo::FORMAT_ARGB32
#define MY_CAIRO_SLANT Cairo::FONT_SLANT_NORMAL
#define MY_CAIRO_WEIGHT Cairo::FONT_WEIGHT_NORMAL
#endif

#ifdef WITH_CGAL
#include <CGAL/Point_2.h>
#include <CGAL/number_utils.h>
#endif

class Painter {

public:
    Painter(const tBox &_bounds, uint _resolution) {
        bounds = _bounds;

        for (uint d = 0; d < 2; ++d) {
            img.low[d] = 0;// image starts at 0
            img.high[d] = img.low[d] + (bounds.high[d] - bounds.low[d]) * _resolution;
        }

        cs = Cairo::ImageSurface::create(MY_CAIRO_FORMAT, imgDim(0), imgDim(1));
        cr = Cairo::Context::create(cs);

        // draw background white
        cr->set_source_rgb(1, 1, 1);
        cr->paint();

        cr->set_line_width(1.0);
        cr->set_source_rgb(0, 0, 0);

        // set font options
        cr->select_font_face("serif", MY_CAIRO_SLANT,
                             MY_CAIRO_WEIGHT);
        cr->set_font_size(FONT_SIZE);

        // draw bounds
        cr->rectangle(translatePoint(bounds.low[0], 0),
                      translatePoint(bounds.low[1], 1),
                      translateLength(bounds.high[0] - bounds.low[0], 0),
                      translateLength(bounds.high[1] - bounds.low[1], 1));
        cr->stroke();
    }

    Painter(const Painter &a) {
        bounds = a.bounds;
        img = a.img;

        cs = Cairo::ImageSurface::create(MY_CAIRO_FORMAT, imgDim(0), imgDim(1));
        cr = Cairo::Context::create(cs);

        // set font options
        cr->select_font_face("serif", MY_CAIRO_SLANT,
                             MY_CAIRO_WEIGHT);
        cr->set_font_size(FONT_SIZE);

        // copy a's surface
        cr->save();

        cr->set_source(a.cs, 0, 0);
        cr->paint();

        cr->restore();
    }

    void save(const std::string &file) const {
        cs->write_to_png(file + ".png");
    }

    void setColor(const std::array<float, 3> &c) {
        setColor(c[0], c[1], c[2]);
    }

    void setColor(float r, float g, float b,
                  float alpha = 1) {
        cr->set_source_rgba(r, g, b, alpha);
    }

    void draw(const tFloatVector &point) {
        cr->arc(translatePoint(point[0], 0),
                translatePoint(point[1], 1), 5, 0, 2 * M_PI);
        cr->fill();
    }

    void drawSquare(const tFloatVector &point) {
        cr->rectangle(translatePoint(point[0], 0) - 5,
                      translatePoint(point[1], 1) - 5, 10, 10);
        cr->fill();
    }

#ifdef WITH_CGAL
    template<typename K>
    void drawSquare(const CGAL::Point_2<K> &point) {
        cr->rectangle(translatePoint(CGAL::to_double(point[0]), 0) - 5,
                      translatePoint(CGAL::to_double(point[1]), 1) - 5, 10, 10);
        cr->fill();
    }
#endif

    void draw(const tFloatVector &point, const tIndex &idx, bool label = true) {
        cr->arc(translatePoint(point[0], 0),
                translatePoint(point[1], 1), 5, 0, 2 * M_PI);
        cr->fill();

        if (label) {
            cr->move_to(translatePoint(point[0], 0) + 7,
                        translatePoint(point[1], 1) + 7);
            cr->set_font_size(10);
            cr->show_text(std::to_string(idx));
            cr->set_font_size(FONT_SIZE);
        }
    }

    void draw(const tPoints &points, bool label = true) {
        for (tIndex i = 0; i < points.size(); ++i) {
            draw(points[i], i, label);
        }
    }

    void draw(const tIndexSet &pointSet, const tPoints &points) {
        for (const auto &p : pointSet) {
            draw(points[p], p);
        }
    }

    void drawCircle(const tFloatVector &center, const tFloat &radius) {

        // draw center
        cr->arc(translatePoint(center[0], 0),
                translatePoint(center[1], 1), 5, 0, 2 * M_PI);
        cr->fill();

        cr->arc(translatePoint(center[0], 0),
                translatePoint(center[1], 1),
                translateLength(radius, 0), 0, 2 * M_PI);
        cr->stroke();
    }

    void drawLine(const tFloatVector &a, const tFloatVector &b) {
        cr->move_to(translatePoint(a[0], 0), translatePoint(a[1], 1));
        cr->line_to(translatePoint(b[0], 0), translatePoint(b[1], 1));
        cr->stroke();
    }

#ifdef WITH_CGAL
    template<typename K>
    void drawLine(const CGAL::Point_2<K> &a, const CGAL::Point_2<K> &b) {
        cr->move_to(translatePoint(CGAL::to_double(a[0]), 0), translatePoint(CGAL::to_double(a[1]), 1));
        cr->line_to(translatePoint(CGAL::to_double(b[0]), 0), translatePoint(CGAL::to_double(b[1]), 1));
        cr->stroke();
    }
#endif

    void drawCone(const tFloatVector &apex, const tDim &k, const tDim &K, const tFloat length = 1) {
        tFloat lowerPhi = k * 2 * M_PI / K;
        tFloat upperPhi = (k + 1) * 2 * M_PI / K;

        drawCone(apex, lowerPhi, upperPhi, length);
    }

    void drawCone(const tFloatVector &apex, const tFloat &lowerPhi, const tFloat &upperPhi, const tFloat length = 1) {
        drawLine(apex, apex + tFloatVector({length * std::cos(lowerPhi),
                                            length * std::sin(lowerPhi)}));
        drawLine(apex, apex + tFloatVector({length * std::cos(upperPhi),
                                            length * std::sin(upperPhi)}));
    }

    void drawGrid(const tIndex _nCells) {
        auto numCellsPerDim = std::ceil(std::sqrt(_nCells));

        tFloatVector cellSize;
        for (tDim d = 0; d < bounds.high.size(); ++d) {
            cellSize[d] = (bounds.high[d] - bounds.low[d]) / numCellsPerDim;
        }

        for (tIndex x = 0; x < numCellsPerDim; ++x) {
            for (tIndex y = 0; y < numCellsPerDim; ++y) {
                cr->rectangle(translatePoint(x * cellSize[0], 0),
                              translatePoint(y * cellSize[1], 1),
                              translateLength(cellSize[0], 0),
                              translateLength(cellSize[1], 1));
            }
        }

        cr->stroke();
    }

    void drawCones(const tFloatVector &apex, const tDim &K, const tFloat length = 1) {
        for (tDim k = 0; k < K; ++k) {
            drawCone(apex, k, K, length);
        }
    }

    template<tDim K>
    void draw(const tIndex &i, const tYaoVertex<K> &v, const tPoints &points) {
        for (tDim k = 0; k < K; ++k) {
            if (v.neighbor[k] != tIndex(-1)) {
                drawLine(points[i], points[v.neighbor[k]]);
            }
        }
    }

    template<typename V>
    void draw(const tYaoGraph<V> &g, const tPoints &points) {
        for (tIndex i = 0; i < g.size(); ++i) {
            draw(i, g[i], points);
        }
    }

    void setDash() {
        cr->set_dash(std::vector<double>({2, 2}), 0);
    }

    void unsetDash() {
        cr->unset_dash();
    }

private:
    float imgDim(const uint dim) {
        return img.high[dim] - img.low[dim];
    }

    float translatePoint(float in, uint dim) {
        return ((in - bounds.low[dim]) /
                ((bounds.high[dim] - bounds.low[dim]))) *
               imgDim(dim);
    }

    float translateLength(float in, uint dim) {
        return (in / ((bounds.high[dim] - bounds.low[dim]))) * imgDim(dim);
    }

private:
    Cairo::RefPtr<Cairo::ImageSurface> cs;
    Cairo::RefPtr<Cairo::Context> cr;

    tBox bounds;
    tBox img;

    static constexpr uint FONT_SIZE = 15;
};