cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

project(NECI Fortran CXX C)

# Get access to custom cmake modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include(neci_system)

# TODO:
#    i) We shouldn't have to rebuild everything for just changing version info --> restrict to header
#          Use -- set_property( SOURCE <> PROPERTY COMPILE_DEFINITIONS ... )
#   ii) Which targets should be included in ALL (see EXCLUDE_FROM_ALL)

# TEST WITH:
# toolchain on archer
# PGI

# Declare options that this project uses
# <Package>_NECI has a custom finder that does a bit more work than normal...


  if (CMAKE_VERSION VERSION_LESS "3.1")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  else ()
    set (CMAKE_CXX_STANDARD 11)
  endif ()
  
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
    DEFAULT off
    DESCRIPTION "Use shared memory to store integrals"
    REQUIRED_PACKAGES HDF5_NECI ) # HDF5_NECI adds the Fortran COMPONENT to the default (and the ability to build)

neci_add_option(
    FEATURE MOLCAS
    DEFAULT OFF
    DESCRIPTION "Build NECI with molcas functionality" )

neci_add_option(
    FEATURE SHARED
    DEFAULT OFF
    DESCRIPTION "Build using shared libraries, rather than static linking" )

# Special option overrides

if( HAVE_SHARED_MEMORY AND CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" )
    message(FATAL_ERROR "Shared memory not supported in PGI builds. Run cmake again with -DENABLE_SHARED_MEMORY=OFF")
endif()

# REQUIRED packages
# <Package>_NECI has a custom finder that does a bit more work than normal...

find_package( LAPACK_NECI REQUIRED )

# Build-type configuration

set(${PROJECT_NAME}_CONFIGURATIONS neci kneci dneci mneci kdneci kmneci)
set(${PROJECT_NAME}_neci_DEFINITIONS HElement_t=real)
set(${PROJECT_NAME}_kneci_DEFINITIONS HElement_t=complex __CMPLX)
set(${PROJECT_NAME}_dneci_DEFINITIONS HElement_t=real __DOUBLERUN)
set(${PROJECT_NAME}_mneci_DEFINITIONS HElement_t=real __PROG_NUMRUNS)
set(${PROJECT_NAME}_kdneci_DEFINITIONS HElement_t=complex __CMPLX __DOUBLERUN)
set(${PROJECT_NAME}_kmneci_DEFINITIONS HElement_t=complex __CMPLX __PROG_NUMRUNS)


# Global (always on) defines

list(APPEND NECI_GLOBAL_DEFINES
    __Linux
    MAXMEM=99999
    DSFMT_MEXP=19937
    _CONFIG="cmake-${CMAKE_Fortran_COMPILER_ID}-$<$<BOOL:${HAVE_MPI}>:par>-${CMAKE_BUILD_TYPE}"
    _VCS_VER="${${PROJECT_NAME}_GIT_SHA1}"

    # Only use SSE2 if not disabled
    $<$<NOT:$<BOOL:${${PROJECT_NAME}_DISABLE_SSE2}>>:HAVE_SSE2>

    # Is this a 32/64bit build
    $<$<BOOL:${${PROJECT_NAME}_64BIT_BUILD}>:__INT64>
    $<$<BOOL:${${PROJECT_NAME}_64BIT_BUILD}>:POINTER8>

    # Is this a debug build
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},DEBUG>:__DEBUG>

    # Global defines that depend on options
    $<$<BOOL:${HAVE_MPI}>:PARALLEL>
    $<$<BOOL:${HAVE_SHARED_MEMORY}>:__SHARED_MEM>
    $<$<NOT:$<BOOL:${HAVE_FFTW}>>:DISABLE_FFTW>
    $<$<BOOL:${HAVE_HDF5}>:__USE_HDF>
    $<$<BOOL:${HAVE_MOLCAS}>:_MOLCAS_>

    # Identify the compiler to the code
    $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},Intel>:__IFORT>
    $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},GNU>:__GFORTRAN__>
    $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},PGI>:__PGI>
)

# Add the source directory as an include path for everything
include_directories( ${CMAKE_SOURCE_DIR}/src )

# Include source directory

add_subdirectory(src)
add_subdirectory(unit_tests)

# And we are done

neci_print_summary()

