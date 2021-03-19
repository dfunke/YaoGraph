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

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(Cairomm_PKGCONF cairomm-1.16)

# Main include dir
find_path(Cairomm_INCLUDE_DIR
  NAMES cairomm/cairomm.h
  PATHS ${Cairomm_PKGCONF_INCLUDE_DIRS}
  PATH_SUFFIXES cairomm-1.16
)

find_path(Cairommconfig_INCLUDE_DIR
  NAMES cairommconfig.h
  PATHS ${Cairomm_PKGCONF_INCLUDE_DIRS}
  PATH_SUFFIXES cairomm-1.16
)

libfind_library(Cairomm cairomm 1.16)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(Cairomm_PROCESS_INCLUDES Cairomm_INCLUDE_DIR Cairommconfig_INCLUDE_DIR SigC++_INCLUDE_DIRS Cairo_INCLUDE_DIRS)
set(Cairomm_PROCESS_LIBS Cairomm_LIBRARY Cairo_LIBRARIES SigC++_LIBRARIES)
libfind_process(Cairomm)

