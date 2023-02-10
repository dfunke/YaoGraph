//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include <memory>

#include "Bubbles.hpp"
#include "Circle.hpp"
#include "Gaussian.hpp"
#include "Grid.hpp"
#include "Road.hpp"
#include "Stars.hpp"
#include "Uniform.hpp"

std::unique_ptr<GeneratorBase> getGen(const char &dist, const tIndex &seed) {
    // [_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar]

    switch (dist) {
        case 'u':
            return std::make_unique<Uniform>(seed);
        case 'g':
            return std::make_unique<Gaussian>(seed);
        case 'd':
            return std::make_unique<Grid>(seed);
        case 'r':
            return std::make_unique<Road>(seed);
        case 's':
            return std::make_unique<Stars>(seed);
        case 'b':
            return std::make_unique<Bubbles>(seed);
        case 'c':
            return std::make_unique<Circle>(seed);

        default:
            return nullptr;
    }
}

std::string getGenName(const char &dist) {
    // [_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar]

    switch (dist) {
        case 'u':
            return "uni";
        case 'g':
            return "gaussian";
        case 'd':
            return "grid";
        case 'r':
            return "road";
        case 's':
            return "stars";
        case 'c':
            return "circle";
        case 'b':
            return "bubbles";

        default:
            return nullptr;
    }
}