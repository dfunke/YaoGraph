# - Try to find SigC++-3.0
# Once done, this will define
#
#  SigC++_FOUND - system has SigC++
#  SigC++_INCLUDE_DIRS - the SigC++ include directories
#  SigC++_LIBRARIES - link these to use SigC++

include(LibFindMacros)

foreach (dir ${CMAKE_SYSTEM_PREFIX_PATH})
    #    message(INFO ${dir})
    file(GLOB SigC++_VERSION_PATH "${dir}/include/sigc++*")
    #    message(INFO ${SigC++_VERSION_PATH})
    if (NOT ${SigC++_VERSION_PATH} STREQUAL "")
        break()
    endif ()
endforeach ()

#message(INFO ${SigC++_VERSION_PATH})
get_filename_component(SigC++_VERSION_DIR ${SigC++_VERSION_PATH} NAME)
#message(INFO ${SigC++_VERSION_DIR})
string(REPLACE "sigc++-" "" SigC++_VERSION ${SigC++_VERSION_DIR})
#message(INFO ${SigC++_VERSION})

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(SigC++_PKGCONF ${SigC++_VERSION_DIR})

# Main include dir
find_path(SigC++_INCLUDE_DIR
        NAMES sigc++/sigc++.h
        PATHS ${SigC++_PKGCONF_INCLUDE_DIRS}
        PATH_SUFFIXES ${SigC++_VERSION_DIR}
        )

# Glib-related libraries also use a separate config header, which is in lib dir
find_path(SigC++Config_INCLUDE_DIR
        NAMES sigc++config.h
        PATHS ${SigC++_PKGCONF_INCLUDE_DIRS} /usr
        PATH_SUFFIXES lib/${SigC++_VERSION_DIR}/include
        )

libfind_library(SigC++ sigc ${SigC++_VERSION})

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(SigC++_PROCESS_INCLUDES SigC++_INCLUDE_DIR SigC++Config_INCLUDE_DIR)
set(SigC++_PROCESS_LIBS SigC++_LIBRARY)
libfind_process(SigC++)

