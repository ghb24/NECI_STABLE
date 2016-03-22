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


