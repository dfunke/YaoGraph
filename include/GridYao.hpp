//
// Created by Daniel Funke on 18.03.21.
//

#pragma once

#include "Types.hpp"
#include "Predicates.hpp"

template<tDim K>
class GridYao {

private:

    class Grid {

    private:
        friend class GridYao<K>;

        using tGridCell = std::vector<tIndex>;
        using tGrid = std::vector<tGridCell>;

    public:

        Grid(const tBox &_bounds, const tIndex _nCells) : bounds(_bounds) {
            numCellsPerDim = std::ceil(std::sqrt(_nCells));
            numCells = numCellsPerDim * numCellsPerDim;

            for (tDim d = 0; d < bounds.high.size(); ++d) {
                cellSize[d] = (bounds.high[d] - bounds.low[d]) / numCellsPerDim;
            }

            cells.resize(numCells);
        }

        void buildGrid(const tPoints &points) {
            for (tIndex i = 0; i < points.size(); ++i) {
                tIndex idx = getIndex(points[i]);
                cells[idx].push_back(i);
            }
        }

        tIndexVector getIndexVector(const tFloatVector &p) const {
            tIndexVector idx;
            for (tDim d = 0; d < idx.size(); ++d) {
                idx[d] = std::floor(p[d] / cellSize[d]);
            }
            return idx;
        }

        tIndex getIndex(const tIndexVector &idx) const {
            tIndex i = 0;
            for (tDim d = idx.size() - 1; d < idx.size(); --d) {
                i += std::pow(numCellsPerDim, d) * idx[d];
            }
            return i;
        }

        tIndex getIndex(const tFloatVector &p) const {
            return getIndex(getIndexVector(p));
        }

        const tGridCell &operator[](const tIndex &i) const {
            return cells[i];
        }

        tGridCell &operator[](const tIndex &i) {
            return cells[i];
        }

        const tGridCell &operator[](const tIndexVector &i) const {
            return cells[getIndex(i)];
        }

        tGridCell &operator[](const tIndexVector &i) {
            return cells[getIndex(i)];
        }

    private:
        tIndex numCellsPerDim;
        tIndex numCells;
        tFloatVector cellSize;
        tBox bounds;
        tGrid cells;
    };

public:
    using tVertex = tYaoVertex<K>;
    using tGraph = tYaoGraph<tVertex>;

    tGraph operator()(const tPoints &points, const tBox &bounds, const tIndex &cellOcc) const {
        tGraph graph(points.size());
        Grid grid(bounds, std::ceil(points.size() / static_cast<tFloat>(cellOcc)));
        grid.buildGrid(points);

        for (tIndex i = 0; i < points.size(); ++i) {

            for (tIndex radius = 0;
                 !isFinalized(graph[i], radius, grid.cellSize) && radius < grid.numCellsPerDim; ++radius) {
                searchRadius(i, radius, graph, grid, points);
            }
        }

        return graph;
    }

private:

    bool isFinalized(const tVertex &v, const tIndex &radius, const tFloatVector &cellSize) const {

        auto minLength = *std::min_element(cellSize.begin(), cellSize.end());

        for (tDim k = 0; k < v.neighbor.size(); ++k) {
            if (v.neighbor[k] == tIndex(-1)) {
                // no neighbor found so far
                return false;
            } else {
                // we found a neighbor for given search radius, can we find a closer one still
                if (v.distance[k] > std::pow((radius - 1) * minLength, 2)) {
                    // found neighbor is further away than any point for next cell radius could be
                    return false;
                }
            }
        }

        return true;
    }

    void searchRadius(const tIndex &point, const int &radius,
                      tGraph &graph, const Grid &grid, const tPoints &points) const {

        auto idx = grid.getIndexVector(points[point]);

        for (int dx = -radius; dx <= radius; ++dx) {
            for (int dy = -radius; dy <= radius; ++dy) {

                if (std::max(std::abs(dx), std::abs(dy)) < radius) {
                    // skip inner cells, already visited in prior calls;
                    continue;
                }

                auto cell = idx + std::array<int, idx.size()>({dx, dy});
                if (std::all_of(cell.begin(), cell.end(),
                                [&grid](const auto &x) { return x < grid.numCellsPerDim; })) {
                    visitGridCell(point, cell, graph, grid, points);
                }
            }
        }
    }

    void visitGridCell(const tIndex &point, const tIndexVector &cell,
                       tGraph &graph, const Grid &grid, const tPoints &points) const {

        auto &vertex = graph[point];

        for (const auto &p : grid[cell]) {

            if (p == point) continue;

            auto d = distance2(points[point], points[p]);
            tDim sec = std::floor(atan2P(points[point], points[p]) / (2 * M_PI / K));

            if (d < vertex.distance[sec]) {
                vertex.neighbor[sec] = p;
                vertex.distance[sec] = d;
            }

        }

    }

};