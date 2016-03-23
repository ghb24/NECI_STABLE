# Inspired by ecbuild_add_executable in ecbuild by the ECMWF

# neci_add_executable
# ===================
#
# Add an executable with a given list of source files. ::
#
#   neci_add_executable(   TARGET <name>
#                          SOURCES <source1> [<source2> ...]
#                        [ LIBS <lib1> [<lib2> ...] ]
#
# Options
# -------
#
# TARGET : required
#   target name
#
# SOURCES : required
#   list of source files
#
# LIBS : optional
#   list of libraries to link against (CMake targets, or external libraries)
#
##############################################################################


macro( neci_add_executable )

    set( options  )
    set( single_value_args TARGET )
    set( multi_value_args SOURCES LIBS )

    cmake_parse_arguments( _p "${options}" "${single_value_args}" "${multi_value_args}" ${_FIRST_ARG} ${ARGN} )

    # Check that the target option (required) is set appropriately

    if( NOT _p_TARGET  )
      message( FATAL_ERROR "The call to neci_add_executable() doesn't specify the TARGET." )
    endif()

    if( NOT _p_SOURCES )
      message( FATAL_ERROR "The call to neci_add_executable() doesn't specify the SOURCES." )
    endif()

    # Actually create the executable

    add_executable( ${_p_TARGET} ${_p_SOURCES} )

    # Add the link libraries

    if ( _p_LIBS )
      list( REMOVE_DUPLICATES _p_LIBS )
      foreach( lib ${_p_LIBS} )
        target_link_libraries( ${_p_TARGET} ${lib} )
      endforeach()
    endif()

    # Where do we put the Fortran modules?

    set_property( TARGET ${_p_TARGET} PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/${_p_TARGET} )

    # Where do the files get built to

    set_property( TARGET ${_p_TARGET} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

    # Add to the global list of libraries
    set( ${PROJECT_NAME}_ALL_EXES ${${PROJECT_NAME}_ALL_EXES} ${_p_TARGET} CACHE INTERNAL "" )

endmacro()
