# - Try to find KaHIP
# Once done, this will define
#
#  KaHIP_FOUND - system has KaHIP
#  --KaHIP_INCLUDE_DIRS - the KaHIP include directories
#  EJDB_LIBRARIES - link these to use KaHIP

include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(KaHIP_PKGCONF libkahip)

# Include dir
find_path(KaHIP_INCLUDE_DIR
  NAMES kaHIP_interface.h
  PATHS ${KaHIP_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(KaHIP_LIBRARY
  NAMES kahip
  PATHS ${KaHIP_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(KaHIP_PROCESS_INCLUDES KaHIP_INCLUDE_DIR)
set(KaHIP_PROCESS_LIBS KaHIP_LIBRARY)
libfind_process(KaHIP)

