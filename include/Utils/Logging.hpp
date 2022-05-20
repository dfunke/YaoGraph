#pragma once

#ifdef LOG_DEBUG
#define LOG(msg) \
    std::cout << msg
#else
#define LOG(msg) \
    ((void)(0))
#endif