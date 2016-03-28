cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

project(NECI Fortran CXX C)

# Get access to custom cmake modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include(neci_system)

# TODO:
#    i) Correct 32-bit support
#   ii) Add build hostname to the output
#  iii) Add toolchain files
#   iv) We shouldn't have to rebuild everything for just changing version info --> restrict to header
#          Use -- set_property( SOURCE <> PROPERTY COMPILE_DEFINITIONS ... )
#    v) Add unit tests (are they ctest compatible?)
#   vi) Add warnings option
#   ix) Add option to cmake gui
#  xii) Summary of options, etc, at the end. (see feature_summary, ecbuild_feature).
#  xiv) hdf5
#   xv) Which targets should be included in ALL (see EXCLUDE_FROM_ALL)

# Declare options that this project uses
# <Package>_NECI has a custom finder that does a bit more work than normal...

find_package( MPI_NECI )

neci_add_option(
    FEATURE MPI
    DEFAULT ON
    DESCRIPTION "Use MPI for parallel compilation"
    REQUIRED_PACKAGES MPI_NECI )

neci_add_option(
    FEATURE FFTW
    DEFAULT ON
    DESCRIPTION "Enable functionality requiring FFTW"
    REQUIRED_PACKAGES FFTW )

neci_add_option(
    FEATURE SHARED_MEMORY
    DEFAULT ON
    DESCRIPTION "Use shared memory to store integrals"
    REQUIRED_PACKAGES LibRT )

neci_add_option(
    FEATURE HDF5
    DEFAULT ON
    DESCRIPTION "Use shared memory to store integrals"
    REQUIRED_PACKAGES HDF5_NECI ) # HDF5_NECI adds the Fortran COMPONENT to the default (and the ability to build)

# REQUIRED packages
# <Package>_NECI has a custom finder that does a bit more work than normal...

find_package( LAPACK_NECI REQUIRED )

# Build-type configuration

set(${PROJECT_NAME}_CONFIGURATIONS neci kneci dneci mneci)
set(${PROJECT_NAME}_neci_DEFINITIONS HElement_t=real)
set(${PROJECT_NAME}_kneci_DEFINITIONS HElement_t=complex __CMPLX)
set(${PROJECT_NAME}_dneci_DEFINITIONS HElement_t=real __DOUBLERUN)
set(${PROJECT_NAME}_mneci_DEFINITIONS HElement_t=real __PROG_NUMRUNS)

# Global (always on) defines

list(APPEND NECI_GLOBAL_DEFINES
    __Linux
    HAVE_SSE2
    MAXMEM=99999
    DSFMT_MEXP=19937
    _CONFIG="cmake-${CMAKE_Fortran_COMPILER_ID}-$<$<BOOL:${HAVE_MPI}>:par>-${CMAKE_BUILD_TYPE}"
    _VCS_VER="${${PROJECT_NAME}_GIT_SHA1}"

    # Global defines that depend on options
    $<$<BOOL:${HAVE_MPI}>:PARALLEL>
    $<$<BOOL:${HAVE_SHARED_MEMORY}>:__SHARED_MEM>
    $<$<NOT:$<BOOL:${HAVE_FFTW}>>:DISABLE_FFTW>
    $<$<BOOL:${HAVE_HDF5}>:__USE_HDF>
)


# Include source directory

add_subdirectory(src)

# And we are done

neci_print_summary()
