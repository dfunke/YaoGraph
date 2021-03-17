cmake_minimum_required(VERSION 3.17)
project(GeoGraph)

# default to Release building
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message("Defaulting CMAKE_BUILD_TYPE to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
include_directories("${PROJECT_SOURCE_DIR}/include/")

# set default for MARCH
if (NOT MARCH)
    set(MARCH native)
endif ()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSTIONS ON)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSTIONS ON)

# setup modern c++ flags
string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=${MARCH} -pedantic -Wall -Wextra")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=${MARCH} -pedantic -Wall -Wextra")

find_package(CGAL)

include_directories(SYSTEM ${CGAL_INCLUDE_DIRS} ${INCLUDE_DIRECTORIES})

add_executable(GeoGraph main.cpp)
target_link_libraries(GeoGraph ${CGAL_LIBRARIES})
