cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

if (POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif ()

if (POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
endif ()

project(NECI Fortran CXX C)

# Get access to custom cmake modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})


# install git hooks
execute_process( COMMAND "${PROJECT_SOURCE_DIR}/tools/install_hooks.sh" ${PROJECT_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(Python3_FIND_UNVERSIONED_NAMES FIRST)

# We require at least ifort17, previous versions do not support all required options
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    if("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS 17)
        message(FATAL_ERROR "ifort 17 or newer required")
    endif()
endif()

include(neci_system)

# TODO:
#    i) We shouldn't have to rebuild everything for just changing version info --> restrict to header
#          Use -- set_property( SOURCE <> PROPERTY COMPILE_DEFINITIONS ... )
#   ii) Which targets should be included in ALL (see EXCLUDE_FROM_ALL)

# TEST WITH:
# toolchain on archer
# PGI

find_package( MPI REQUIRED )
list(APPEND NECI_GLOBAL_DEFINES USE_MPI)

# Declare options that this project uses
# <Package>_NECI has a custom finder that does a bit more work than normal...

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
    DEFAULT OFF
    DESCRIPTION "Use HDF5 file format for storage"
    REQUIRED_PACKAGES HDF5_NECI ) # HDF5_NECI adds the Fortran COMPONENT to the default (and the ability to build)

neci_add_option(
    FEATURE SHARED
    DEFAULT OFF
    DESCRIPTION "Build using shared libraries, rather than static linking" )

neci_add_option(
    FEATURE DOC
    DEFAULT OFF
    DESCRIPTION "Enable the documentation. (Has additional requirements like pandoc or latexmk)")

# Special option overrides

if( HAVE_SHARED_MEMORY AND CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" )
    message(FATAL_ERROR "Shared memory not supported in PGI builds. Run cmake again with -DENABLE_SHARED_MEMORY=OFF")
endif()

# REQUIRED packages
# <Package>_NECI has a custom finder that does a bit more work than normal...

find_package( BLAS REQUIRED )
find_package( LAPACK REQUIRED )

# Build-type configuration

set(${PROJECT_NAME}_CONFIGURATIONS neci kneci dneci mneci kdneci kmneci)
add_custom_target(programs)
add_dependencies(programs ${${PROJECT_NAME}_CONFIGURATIONS})
set(${PROJECT_NAME}_neci_DEFINITIONS HElement_t=real)

set(${PROJECT_NAME}_kneci_DEFINITIONS HElement_t=complex CMPLX_)
set(${PROJECT_NAME}_dneci_DEFINITIONS HElement_t=real DOUBLERUN_)
set(${PROJECT_NAME}_mneci_DEFINITIONS HElement_t=real PROG_NUMRUNS_)
set(${PROJECT_NAME}_kdneci_DEFINITIONS HElement_t=complex CMPLX_ DOUBLERUN_)
set(${PROJECT_NAME}_kmneci_DEFINITIONS HElement_t=complex CMPLX_ PROG_NUMRUNS_)

set(${PROJECT_NAME}_neci_WARNERR TRUE)
set(${PROJECT_NAME}_dneci_WARNERR TRUE)
set(${PROJECT_NAME}_mneci_WARNERR TRUE)

# Global (always on) defines

list(APPEND NECI_GLOBAL_DEFINES
    Linux_
    MAXMEM=99999
    DSFMT_MEXP=19937
    _CONFIG="cmake-${CMAKE_Fortran_COMPILER_ID}-$<$<BOOL:${HAVE_MPI}>:par>-${CMAKE_BUILD_TYPE}"
    _VCS_VER="${${PROJECT_NAME}_GIT_SHA1}"

    # Only use SSE2 if not disabled
    $<$<NOT:$<BOOL:${${PROJECT_NAME}_DISABLE_SSE2}>>:HAVE_SSE2>

    # Is this a 32/64bit build
    $<$<BOOL:${${PROJECT_NAME}_64BIT_BUILD}>:INT64_>
    $<$<BOOL:${${PROJECT_NAME}_64BIT_BUILD}>:POINTER8>

    # Is this a debug build
    $<$<STREQUAL:${CMAKE_BUILD_TYPE},DEBUG>:DEBUG_>

    # Global defines that depend on options
    $<$<BOOL:${HAVE_SHARED_MEMORY}>:SHARED_MEM_>
    $<$<NOT:$<BOOL:${HAVE_FFTW}>>:DISABLE_FFTW>
    $<$<BOOL:${HAVE_HDF5}>:USE_HDF_>

    # Identify the compiler to the code
    $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},Intel>:IFORT_>
    $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},IntelLLVM>:INTELLLVM_>
    $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},GNU>:GFORTRAN_>
    $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},PGI>:PGI_>
)

# Add the source directory as an include path for everything
include_directories( ${PROJECT_SOURCE_DIR}/src )

# Include source directory

add_subdirectory(src)

# Don't do NECI tests if NECI is added to other project.
# BUILD_TESTING option comes from `include(CTest)`
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    add_subdirectory(unit_tests)
endif()

if (HAVE_DOC)
    add_subdirectory(docs)
endif()

add_subdirectory(python)

# And we are done
neci_print_summary()

