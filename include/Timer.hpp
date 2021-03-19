//
// Created by Daniel Funke on 17.03.21.
//

#pragma once

#include <chrono>

template<typename T>
class Timer {

public:

    template<typename... Args>
    static auto time(Args... args) {
        T obj;

        auto t1 = std::chrono::high_resolution_clock::now();
        auto res = obj(args...);
        auto t2 = std::chrono::high_resolution_clock::now();

        return std::make_tuple(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count(), res);
    }

};