# - Try to find EJDB
# Once done, this will define
#
#  EJDB_FOUND - system has SigC++
#  --EJDB_INCLUDE_DIRS - the SigC++ include directories
#  EJDB_LIBRARIES - link these to use SigC++

include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(EJDB_PKGCONF libejdb)

# Include dir
find_path(EJDB_INCLUDE_DIR
  NAMES ejdb.h
  PATHS ${EJDB_PKGCONF_INCLUDE_DIRS}
  PATH_SUFFIXES ejdb
)

# Finally the library itself
find_library(EJDB_LIBRARY
  NAMES ejdb
  PATHS ${EJDB_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(EJDB_PROCESS_INCLUDES EJDB_INCLUDE_DIR)
set(EJDB_PROCESS_LIBS EJDB_LIBRARY)
libfind_process(EJDB)

