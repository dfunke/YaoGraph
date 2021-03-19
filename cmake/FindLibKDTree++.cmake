# - Try to find libkdtree++
# Once done, this will define
#
#  LibKDTree_FOUND - system has KaHIP
#  --KaHIP_INCLUDE_DIRS - the KaHIP include directories
#  EJDB_LIBRARIES - link these to use KaHIP

# Include dir
find_path(LibKDTree_INCLUDE_DIR
  NAMES kdtree++/kdtree.hpp
)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LibKDTree  DEFAULT_MSG
                                  LibKDTree_INCLUDE_DIR)

set(LibKDTree_INCLUDE_DIRS LibKDTree_INCLUDE_DIR)
