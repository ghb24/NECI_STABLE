# disallow in-source build
# ========================

# Check for failed attempts to build in the source tree
if (EXISTS ${CMAKE_SOURCE_DIR}/CMakeCache.txt )
    message( FATAL_ERROR "Project ${PROJECT_NAME} contains a CMakeCache.txt inside source tree [${CMAKE_SOURCE_DIR}/CMakeCache.txt].\n Please remove it and
    make sure that source tree is prestine and clean of unintended files, before retrying." )
endif()

get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

if(${srcdir} STREQUAL ${bindir})
    message("")
    message("######################################################")
    message("You are attempting to build in your source directory (${srcdir}).")
    message("You must run cmake from a different build directory.")
    message("######################################################")
    message( FATAL_ERROR "${PROJECT_NAME} requires an out of source build.\n Please create a separate build directory and run 'cmake path/to/project [options]' from there.")
endif()

# Hard ensure that cmake is new enough (if project forgot...)
# ===========================================================

set( NECI_CMAKE_MINIMUM "2.8.4" )
if( ${CMAKE_VERSION} VERSION_LESS ${NECI_CMAKE_MINIMUM} )
    message(FATAL_ERROR "${PROJECT_NAME} requires at least CMake ${NECI_CMAKE_MINIMUM} -- you are using ${CMAKE_COMMAND} [${CMAKE_VERSION}]\n Please, get a newer version of CMake @ www.cmake.org" )
endif()

# Include our further macros (only if this is the top project)
# ============================================================
if( PROJECT_NAME STREQUAL CMAKE_PROJECT_NAME )

    # CMake provided macro sets

    include(CMakeParseArguments)

    # Our local macros

    include( neci_add_option )
    include( neci_add_library )

endif()
