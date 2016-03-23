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

    # Add extra macros from the contrib/ directory (externally sourced finders)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/contrib" )

    # CMake provided macro sets

    include(CMakeParseArguments)

    # Our local macros

    include( neci_add_option )
    include( neci_add_library )
    include( neci_print_summary )

    include( ${CMAKE_CURRENT_LIST_DIR}/contrib/GetGitRevisionDescription.cmake )

endif()

# ========================================================================================================
# Initialise the project
# ========================================================================================================

# Reset the list of target libraries and executables
set ( ${PROJECT_NAME}_ALL_LIBS "" CACHE INTERNAL "")
set ( ${PROJECT_NAME}_ALL_EXES "" CACHE INTERNAL "")

# This should be a git project, but it is possible that people will copy it, so check for that

if( EXISTS ${PROJECT_SOURCE_DIR}/.git )
    get_git_head_revision( GIT_REFSPEC ${PROJECT_NAME}_GIT_SHA1 )
    if ( NOT ${PROJECT_NAME}_GIT_SHA1 )
        message( WARNING "WARNING: Unable to get git SHA1 ID. Using placeholder")
        set( ${PROJECT_NAME}_GIT_SHA1 "GIT-UNKNOWN" )
    endif()
else()
    message( WARNING "WARNING: Source code not under version control")
endif()

# Read and parse the version file
set( ${PROJECT_NAME}_VERSION_STR "0.0.0" )
if ( EXISTS ${PROJECT_SOURCE_DIR}/VERSION.cmake )
    include( ${PROJECT_SOURCE_DIR}/VERSION.cmake )
else()
    message( WARNING "WARNING: Version number not set. Using default ${${PROJECT_NAME}_VERSION_STR}" )
endif()
