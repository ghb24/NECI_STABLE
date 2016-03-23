# Inspired by ecbuild_add_library in ecbuild by the ECMWF

# neci_add_library
# ===================
#
# Add a library with a given list of source files. ::
#
#   neci_add_library( TARGET <name>
#                        SOURCES <source1> [<source2> ...]
#                        [ TEMPLATED_SOURCES <source1> [<source2> ...] ]
#                        [ TYPE SHARED|STATIC|MODULE ]
#                        [ PRIVATE_INCLUDEs <dir1> [<dir2> ...] ]
#                        [ DEFINITIONS <define1> [<define3> ...] ]
#                        [ OUTPUT_NAME <name> ]
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
# TEMPLATED_SOURCES : optional
#   list of source files to be passed through f90_template.py (note that this ensures that they are
#   not redefined multiple times).
#
# TYPE : optional
#   library type, one of:
#
# PRIVATE_INCLUDES : optional
#   list of paths to add to include directories that will NOT be exported to other projects. Currently
#   the PUBLIC_INCLUDES functionality is not implemented in this project.
#
# DEFINITIONS : optional
#   list of definitions to add to preprocessor defines
#
# OUTPUT_NAME : optional
#   set the name of the output file (so it may differ from the target name)
#
#   :SHARED: libraries are linked dynamically and loaded at runtime
#   :STATIC: archives of object files for use when linking other targets.
#
##############################################################################


macro( neci_add_library )

    set( options  )
    set( single_value_args TARGET TYPE OUTPUT_NAME )
    set( multi_value_args SOURCES DEFINITIONS TEMPLATED_SOURCES PRIVATE_INCLUDES )

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

    # Add .F90.template files if supplied

    set( ${_p_TARGET}_TEMPLATED_SOURCES )
    if ( _p_TEMPLATED_SOURCES )

      # Ensure that the templates get put somewhere unique for each target
      set( _template_dir ${CMAKE_BINARY_DIR}/templated/${_p_TARGET} )
      file( MAKE_DIRECTORY ${_template_dir} )

      get_filename_component( _templater_tool ${CMAKE_SOURCE_DIR}/tools/f90_template.py ABSOLUTE )

      foreach(_templated_file ${_p_TEMPLATED_SOURCES})
        get_filename_component( _templated_file_base ${_templated_file} NAME_WE )
        get_filename_component( _templated_file_absolute ${_templated_file} ABSOLUTE) 
        set( _templated_target_file ${_template_dir}/${_templated_file_base}.F90 )
        list( APPEND ${_p_TARGET}_TEMPLATED_SOURCES ${_templated_target_file} )
        add_custom_command(
            COMMAND ${_templater_tool} ${_templated_file_absolute} ${_templated_target_file}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            OUTPUT ${_templated_target_file}
            DEPENDS ${_templated_file_absolute} )
      endforeach()
    endif()

    # Actually add the library to the cmake build

    add_library( ${_p_TARGET} ${_p_TYPE} ${_p_SOURCES} ${${_p_TARGET}_TEMPLATED_SOURCES} )

    # Add definitions to the compliation

    if (DEFINED _p_DEFINITIONS )
        get_property( _target_defs TARGET ${_p_TARGET} PROPERTY COMPILE_DEFINITIONS )
        list( APPEND _target_defs ${_p_DEFINITIONS} )
        message( STATUS "Library ${_p_TARGET} using definitions: ${_target_defs}" )
        set_property( TARGET ${_p_TARGET} PROPERTY COMPILE_DEFINITIONS ${_target_defs} )
    endif()

    # Add (private) includes

    if( DEFINED _p_PRIVATE_INCLUDES )
      list( REMOVE_DUPLICATES _p_PRIVATE_INCLUDES )
      foreach( include_path ${_p_PRIVATE_INCLUDES} )
        if( "${CMAKE_VERSION}" VERSION_LESS "2.8.11" ) # PRIVATE functionality doesn't exist before 2.8.11
          target_include_directories( ${_p_TARGET} PUBLIC ${include_path} )
        else()
          target_include_directories( ${_p_TARGET} PRIVATE ${include_path} )
        endif()
      endforeach()
    endif()

    # Have we manually set an output name?

    if ( _p_OUTPUT_NAME )
      message( STATUS "Library ${_p_TARGET}: Output name is ${_p_OUTPUT_NAME}" )
      set_property( TARGET ${_p_TARGET} PROPERTY OUTPUT_NAME ${_p_OUTPUT_NAME} )
    endif()

    # Add to the global list of libraries
    #list( APPEND ${PROJECT_NAME}_ALL_LIBS ${_p_TARGET} )
    set( ${PROJECT_NAME}_ALL_LIBS ${${PROJECT_NAME}_ALL_LIBS} ${_p_TARGET} CACHE INTERNAL "" )

endmacro()
