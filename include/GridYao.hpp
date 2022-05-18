//
// Created by Daniel Funke on 18.03.21.
//

#pragma once

#include "Predicates.hpp"
#include "Types.hpp"
#include "utils/InexactKernel.hpp"

template<tDim C, typename Kernel, tIndex CellOcc>
class GridYao {

private:
    using tEFloat = typename Kernel::Float;
    using tEFloatVector = std::array<tEFloat, D>;

    using tKPoint = typename Kernel::Point;
    using tKPoints = std::vector<tKPoint>;

    class Grid {

    private:
        friend class GridYao<C, Kernel, CellOcc>;

        using tGridCell = std::vector<tIndex>;
        using tGrid = std::vector<tGridCell>;

    public:
        Grid(const tBox &_bounds, const tIndex _nCells) : bounds(_bounds) {
            numCellsPerDim = std::ceil(std::sqrt(_nCells));
            numCells = numCellsPerDim * numCellsPerDim;

            for (tDim d = 0; d < bounds.high.size(); ++d) {
                cellSize[d] = (bounds.high[d] - bounds.low[d]) / numCellsPerDim;
            }
            minCellSize = *std::min_element(cellSize.begin(), cellSize.end());

            cells.resize(numCells);
        }

        void buildGrid(const tKPoints &points) {
            for (tIndex i = 0; i < points.size(); ++i) {
                tIndex idx = getIndex(points[i]);
                cells[idx].push_back(i);
            }
        }

        tIndexVector getIndexVector(const tKPoint &p) const {
            tIndexVector idx;
            for (tDim d = 0; d < idx.size(); ++d) {
                idx[d] = std::floor(Kernel::to_float_exact(p[d] / cellSize[d]));//TODO exact?
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

        tIndex getIndex(const tKPoint &p) const {
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
        tEFloatVector cellSize;
        tEFloat minCellSize;
        tBox bounds;
        tGrid cells;
    };

public:
    using tVertex = tYaoVertex<C, tEFloat>;
    using tGraph = tYaoGraph<tVertex>;

    static std::string name() {
        return "GridYao" + std::to_string(CellOcc) + "_" + Kernel::name();
    }

    tGraph operator()(const tPoints &iPoints, const tBox &bounds) const {
        tGraph graph(iPoints.size());

        auto rays = Kernel::computeCones(C);
        auto kPoints = Kernel::mkPoints(iPoints);

        Grid grid(bounds, std::ceil(kPoints.size() / static_cast<tIFloat>(CellOcc)));
        grid.buildGrid(kPoints);

        for (tIndex i = 0; i < kPoints.size(); ++i) {

            auto cLines = Kernel::computePointCones(kPoints[i], rays);

            for (tIndex radius = 0;
                 !isFinalized(graph[i], radius, grid.minCellSize) && radius < grid.numCellsPerDim; ++radius) {
                searchRadius(i, radius, graph, grid, kPoints, cLines);
            }
        }

        return graph;
    }

private:
    bool isFinalized(const tVertex &v, const tIndex &radius, const tEFloat &minCellSize) const {

        auto d = static_cast<int>(radius - 1) * minCellSize;//TODO better type handling, exact computation;
        auto d2 = d * d;

        for (tDim k = 0; k < v.neighbor.size(); ++k) {
            if (v.neighbor[k] == INF_IDX) {
                // no neighbor found so far
                return false;
            } else {
                // we found a neighbor for given search radius, can we find a closer one still
                if (v.distance[k] > d2) {
                    // found neighbor is further away than any point for next cell radius could be
                    return false;
                }
            }
        }

        return true;
    }

    void searchRadius(const tIndex &point, const int &radius,
                      tGraph &graph, const Grid &grid, const tKPoints &points, const auto &cLines) const {

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
                    visitGridCell(point, cell, graph, grid, points, cLines);
                }
            }
        }
    }

    void visitGridCell(const tIndex &point, const tIndexVector &cell,
                       tGraph &graph, const Grid &grid, const tKPoints &points, const auto &cLines) const {

        auto &vertex = graph[point];

        for (const auto &p : grid[cell]) {

            if (p == point) continue;

            auto sec = Kernel::getCone(points[point], points[p], cLines);

            if (vertex.neighbor[sec] == INF_IDX || Kernel::compareDistance(points[point], points[p], points[vertex.neighbor[sec]], vertex.distance[sec])) {
                vertex.neighbor[sec] = p;
                vertex.distance[sec] = Kernel::distance2(points[point], points[p]);
            }
        }
    }
};