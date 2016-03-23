# Inspired by ecbuild_add_library in ecbuild by the ECMWF

# neci_add_library
# ===================
#
# Add a library with a given list of source files. ::
#
#   neci_add_library( TARGET <name>
#                        SOURCES <source1> [<source2> ...]
#                        [ TYPE SHARED|STATIC|MODULE ]
# Options
# -------
#
# TARGET : required
#   target name
#
# SOURCES : required
#   list of source files
#
# TYPE : optional
#   library type, one of:
#
# DEFINITIONS : optional
#   list of definitions to add to preprocessor defines
#
#   :SHARED: libraries are linked dynamically and loaded at runtime
#   :STATIC: archives of object files for use when linking other targets.
#
##############################################################################


macro( neci_add_library )

    set( options  )
    set( single_value_args TARGET TYPE )
    set( multi_value_args SOURCES DEFINITIONS )

    cmake_parse_arguments( _p "${options}" "${single_value_args}" "${multi_value_args}" ${_FIRST_ARG} ${ARGN} )

    # Check that the target option (required) is set appropriately

    if( NOT _p_TARGET  )
      message( FATAL_ERROR "The call to neci_add_library() doesn't specify the TARGET." )
    endif()

    if( NOT _p_SOURCES )
      message( FATAL_ERROR "The call to neci_add_library() doesn't specify the SOURCES." )
    endif()

    # Check that the type of the library is correctly specified

    if ( _p_TYPE )
      if ( NOT _p_TYPE MATCHES "STATIC" AND NOT _p_TYPE MATCHES "SHARED" )
        message( FATAL_ERROR "Library ${_p_TARGET}: type must be one of [ STATIC | SHARED ]" )
      endif()
      message( STATUS "Library ${_p_TARGET}: library type is ${_p_TYPE}" )
    else()
      message( FATAL_ERROR "Library ${_p_TARGET}: Type not specified. Must be one of [ STATIC | SHARED ]" )
    endif()

    # Actually add the library to the cmake build

    add_library( ${_p_TARGET} ${_p_TYPE} ${_p_SOURCES} )

    # Add definitions to the compliation

    if (DEFINED _p_DEFINITIONS )
        get_property( _target_defs TARGET ${_p_TARGET} PROPERTY COMPILE_DEFINITIONS )
        list( APPEND _target_defs ${_p_DEFINITIONS} )
        message( STATUS "Library ${_p_TARGET} using definitions: ${_target_defs}" )
        set_property( TARGET ${_p_TARGET} PROPERTY COMPILE_DEFINITIONS ${_target_defs} )
    endif()

    # Add to the global list of libraries
    #list( APPEND ${PROJECT_NAME}_ALL_LIBS ${_p_TARGET} )
    set( ${PROJECT_NAME}_ALL_LIBS ${${PROJECT_NAME}_ALL_LIBS} ${_p_TARGET} CACHE INTERNAL "" )

endmacro()
