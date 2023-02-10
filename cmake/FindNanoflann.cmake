# - Try to find nanoflann
# Once done, this will define
#
#  Nanoflann_FOUND - system has Nanoflann
#  Nanoflann_INCLUDE_DIRS - the Nanoflann include directories

# Include dir
find_path(Nanoflann_INCLUDE_DIR
  NAMES nanoflann.hpp
)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Nanoflann  DEFAULT_MSG
                                  Nanoflann_INCLUDE_DIR)

set(Nanoflann_INCLUDE_DIRS Nanoflann_INCLUDE_DIR)
