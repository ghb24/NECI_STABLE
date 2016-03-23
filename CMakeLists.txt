cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

project(neci Fortran CXX C)

# Get access to custom cmake modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include(neci_system)

# TODO:
#    i) Add source code in subdirs
#   ii) Add templating tools
#  iii) Add toolchain files
#   iv) Add executables
#    v) Add unit tests (are they ctest compatible?)
#   vi) Add searching for packages within options (e.g. MPI...)
#  vii) Pick up VERSION.cmake
# viii) Add configuration name to _CONFIG define
#   ix) Add option to cmake gui
#    x) Attack GIT to get version info
#   xi) List all configured (a) libraries, (b) executables at the end
#  xii) Summary of options, etc, at the end.

# Declare options that this project uses

neci_add_option(
    FEATURE MPI
    DEFAULT ON
    DESCRIPTION "Use MPI for parallel compilation")

neci_add_option(
    FEATURE FFTW
    DEFAULT ON
    DESCRIPTION "Enable functionality requiring FFTW"
    REQUIRED_PACKAGES FFTW)


# Global (always on) defines
list(APPEND NECI_GLOBAL_DEFINES
        
    __Linux
    HAVE_SSE2
    MAXMEM=99999
    DSFMT_MEXP=19937
    _CONFIG="cmake-TODO"
)



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
