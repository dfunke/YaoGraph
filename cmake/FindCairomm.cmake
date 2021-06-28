# - Try to find Cairomm 1.16
# Once done, this will define
#
#  Cairomm_FOUND - system has Cairomm
#  Cairomm_INCLUDE_DIRS - the Cairomm include directories
#  Cairomm_LIBRARIES - link these to use Cairomm

include(LibFindMacros)

# Dependencies
libfind_package(Cairomm Cairo)
libfind_package(Cairomm SigC++)

foreach (dir ${CMAKE_SYSTEM_PREFIX_PATH})
    #    message(INFO ${dir})
    file(GLOB Cairomm_VERSION_PATH "${dir}/include/cairomm*")
    #    message(INFO ${Cairomm_VERSION_PATH})
    if (NOT ${Cairomm_VERSION_PATH} STREQUAL "")
        break()
    endif ()
endforeach ()

#message(INFO ${Cairomm_VERSION_PATH})
get_filename_component(Cairomm_VERSION_DIR ${Cairomm_VERSION_PATH} NAME)
#message(INFO ${Cairomm_VERSION_DIR})
string(REPLACE "cairomm-" "" Cairomm_VERSION ${Cairomm_VERSION_DIR})
#message(INFO ${Cairomm_VERSION})


# Use pkg-config to get hints about paths
libfind_pkg_check_modules(Cairomm_PKGCONF ${Cairomm_VERSION_DIR})

# Main include dir
find_path(Cairomm_INCLUDE_DIR
        NAMES cairomm/cairomm.h
        PATHS ${Cairomm_PKGCONF_INCLUDE_DIRS}
        PATH_SUFFIXES ${Cairomm_VERSION_DIR}
        )

find_path(Cairommconfig_INCLUDE_DIR
        NAMES cairommconfig.h
        PATHS ${Cairomm_PKGCONF_INCLUDE_DIRS}
        PATH_SUFFIXES ${Cairomm_VERSION_DIR}
        )

libfind_library(Cairomm cairomm ${Cairomm_VERSION})

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(Cairomm_PROCESS_INCLUDES Cairomm_INCLUDE_DIR Cairommconfig_INCLUDE_DIR SigC++_INCLUDE_DIRS Cairo_INCLUDE_DIRS)
set(Cairomm_PROCESS_LIBS Cairomm_LIBRARY Cairo_LIBRARIES SigC++_LIBRARIES)
libfind_process(Cairomm)

