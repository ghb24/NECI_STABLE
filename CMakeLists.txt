cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

project(neci Fortran CXX C)

# Get access to custom cmake modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include(neci_system)

# TODO:
#   i) Add source code in subdirs
#  ii) Add templating tools
# iii) Add toolchain files
#  iv) Add executables
#   v) Add unit tests (are they ctest compatible?)

# Declare options that this project uses


# Include source directory

add_subdirectory(src)

# And we are done

message("")
message("  +-------------------------------------------------+")
message("  | NECI configuration now complete                 |")
message("  | VERSION ${NECI_VERSION}                         |")
message("  | SHAID ${NECI_SHAID}                             |")
message("  +-------------------------------------------------+")
message("")
message("  You can now do 'make' to compile the software.")
