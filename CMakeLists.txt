cmake_minimum_required(VERSION 3.12)
project(GeoGraph)

# default to Release building
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message("Defaulting CMAKE_BUILD_TYPE to Release")
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type")
endif ()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
include_directories("${PROJECT_SOURCE_DIR}/include/")

# set default for MARCH
if (NOT MARCH)
    set(MARCH native CACHE STRING "ARCH to use")
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

add_executable(GeoGraph GeoGraph.cpp)
add_executable(Generator Generator.cpp)

option(WITH_CGAL "Build CGAL algorithms for comparision" OFF)
if (WITH_CGAL)
    add_compile_definitions(WITH_CGAL)
    find_package(CGAL)
    include_directories(SYSTEM ${CGAL_INCLUDE_DIRS} ${INCLUDE_DIRECTORIES})
    target_link_libraries(GeoGraph ${CGAL_LIBRARIES})
endif ()

option(WITH_CAIRO "Build with Cairo for drawing" OFF)
if (WITH_CAIRO)
    add_compile_definitions(WITH_CAIRO)
    find_package(Cairomm REQUIRED)
    include_directories(SYSTEM "/usr/local/include" ${INCLUDE_DIRECTORIES})
    include_directories(SYSTEM ${Cairomm_INCLUDE_DIRS} ${INCLUDE_DIRECTORIES})
    target_link_libraries(GeoGraph ${Cairomm_LIBRARIES})
endif ()

if (CMAKE_BUILD_TYPE MATCHES Debug)
    message("debug mode - enable logging")
    add_compile_definitions(LOG_DEBUG)
endif ()

option(WITH_VTUNE "Build for Profiling" OFF)
if (WITH_VTUNE)
    add_compile_definitions(VTUNE)
endif ()

find_package(GTest)

if (GTEST_FOUND)
    add_subdirectory(test)
endif ()

set(ROAD_DIR "${PROJECT_SOURCE_DIR}/data" CACHE PATH "Path to road networks")
add_compile_definitions(ROAD_DATA="${ROAD_DIR}")

set(DATA_DIR "${PROJECT_SOURCE_DIR}/data" CACHE PATH "Path to graph data")
add_compile_definitions(DATA_DIR="${DATA_DIR}")

set(STAR_DIR "/global_data/funke/Gaia/" CACHE PATH "Path to gaia data")
add_compile_definitions(STAR_DATA="${STAR_DIR}")