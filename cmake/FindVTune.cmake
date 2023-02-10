
# Copyright (c) 2013-2014 Stefan.Eilemann@epfl.ch
# finds the ittnotify API

#  VTUNE_FOUND - system has VTUNE
#  VTUNE_INCLUDE_DIRS - the VTUNE include directories
#  VTUNE_LIBRARIES - link these to use VTUNE

include(LibFindMacros)

find_program(VTUNE_EXECUTABLE amplxe-cl
  HINTS $ENV{VTUNE_AMPLIFIER_XE_2015_DIR}
  PATH_SUFFIXES bin64
)

if(NOT VTUNE_EXECUTABLE)
  return()
endif()

get_filename_component(VTUNE_DIR ${VTUNE_EXECUTABLE} PATH)
set(VTUNE_DIR "${VTUNE_DIR}/..")

# Include dir
find_path(VTUNE_INCLUDE_DIR
  NAMES ittnotify.h
  PATHS ${VTUNE_DIR}
  PATH_SUFFIXES include
)

# Finally the library itself
find_library(VTUNE_LIBRARY
  NAMES ittnotify
  PATHS ${VTUNE_DIR}
  PATH_SUFFIXES lib64
)

if(UNIX)
  list(APPEND VTUNE_LIBRARY dl)
endif()

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(VTUNE_PROCESS_INCLUDES VTUNE_INCLUDE_DIR)
set(VTUNE_PROCESS_LIBS VTUNE_LIBRARY)
libfind_process(VTUNE)