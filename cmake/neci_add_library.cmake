# Inspired by ecbuild_add_library in ecbuild by the ECMWF

# neci_add_library
# ===================
#
# Add a library with a given list of source files. ::
#
#   neci_add_library(   TARGET <name>
#                       SOURCES <source1> [<source2> ...]
#                       LINKER_LANGUAGE <(C|CXX|Fortran)>
#                     [ TEMPLATED_SOURCES <source1> [<source2> ...] ]
#                     [ FYPP_SOURCES <source1> [<source2> ...] ]
#                     [ TYPE SHARED|STATIC|MODULE ]
#                     [ PRIVATE_INCLUDEs <dir1> [<dir2> ...] ]
#                     [ LIBS <lib1> [<lib2> ...] ]
#                     [ DEFINITIONS <define1> [<define3> ...] ]
#                     [ OUTPUT_NAME <name> ]
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
# LINKER_LANGUAGE : required
#   what language should be used to link the executable. If not specified, cmake will determine this automatically
#   Values: C, CXX, Fortran
#
#   required so that <project>_<lang>_linker_flags can be made to work
#
# TEMPLATED_SOURCES : optional
#   list of source files to be passed through f90_template.py (note that this ensures that they are
#   not redefined multiple times).
#
# FYPP_SOURCES : optional
#   list of source files to be passed through fypp (note that this ensures that they are
#   not redefined multiple times).
#
# RELAX_WARNINGS : optional
#   list of old source files which require relaxed warning policies.
#
# TYPE : optional
#   library type, one of:
#
# PRIVATE_INCLUDES : optional
#   list of paths to add to include directories that will NOT be exported to other projects. Currently
#   the PUBLIC_INCLUDES functionality is not implemented in this project.
#
# LIBS : optional
#   list of libraries to link against (CMake targets, or external libraries)
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

    set( options )
    set( single_value_args TARGET TYPE OUTPUT_NAME LINKER_LANGUAGE WARNERR)
    set( multi_value_args SOURCES DEFINITIONS TEMPLATED_SOURCES FYPP_SOURCES RELAX_WARNINGS PRIVATE_INCLUDES LIBS )

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
      find_package(Python3 REQUIRED)

      # Ensure that the templates get put somewhere unique for each target
      set( _template_dir ${CMAKE_BINARY_DIR}/templated/${_p_TARGET} )
      file( MAKE_DIRECTORY ${_template_dir} )

      get_filename_component( _templater_tool ${PROJECT_SOURCE_DIR}/tools/f90_template.py ABSOLUTE )

      foreach(_templated_file ${_p_TEMPLATED_SOURCES})
        get_filename_component( _templated_file_base ${_templated_file} NAME_WE )
        get_filename_component( _templated_file_absolute ${_templated_file} ABSOLUTE)
        set( _templated_target_file ${_template_dir}/${_templated_file_base}.F90 )
        list( APPEND ${_p_TARGET}_TEMPLATED_SOURCES ${_templated_target_file} )
        add_custom_command(
            COMMAND ${Python3_EXECUTABLE} ${_templater_tool} ${_templated_file_absolute} ${_templated_target_file}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            OUTPUT ${_templated_target_file}
            DEPENDS ${_templated_file_absolute} )
        set_property(SOURCE ${_templated_target_file}  PROPERTY COMPILE_FLAGS ${${PROJECT_NAME}_Fortran_relaxed_WARNING_FLAGS})
      endforeach()
    endif()


    # Add fypp templating.
    # Check if it exists in $PATH, if not init the fypp git submodule
    # Convert all templated files to f90 files.
    set( ${_p_TARGET}_FYPP_SOURCES )
    if (_p_FYPP_SOURCES )
        find_package(Python3 REQUIRED)
        # This will also search in default locations e.g. $PATH, env ...
        find_program (_fypp fypp ${PROJECT_SOURCE_DIR}/External/fypp/bin/)
        if (NOT _fypp)
            message(STATUS "Submodule update")
            find_package(Git QUIET)
            if(GIT_FOUND)
                execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                RESULT_VARIABLE GIT_SUBMOD_RESULT)
                if(NOT GIT_SUBMOD_RESULT EQUAL "0")
                    message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules.")
                endif()
                find_program (_fypp fypp ${PROJECT_SOURCE_DIR}/External/fypp/bin/)
            else()
                message(FATAL_ERROR "Git not found")
            endif()
        endif()
        set( _fypp_dir ${CMAKE_BINARY_DIR}/fypp/${_p_TARGET} )
        file( MAKE_DIRECTORY ${_fypp_dir} )

        set(_fypp_options)
        list( APPEND _fypp_options -m itertools -m functools)


        foreach(_fypp_file ${_p_FYPP_SOURCES})
            get_filename_component( _fypp_file_base ${_fypp_file} NAME_WE )
            get_filename_component( _fypp_file_absolute ${_fypp_file} ABSOLUTE)
            set( _fypp_target_file ${_fypp_dir}/${_fypp_file_base}.F90 )
            list( APPEND ${_p_TARGET}_FYPP_SOURCES ${_fypp_target_file} )
            add_custom_command(
                COMMAND ${Python3_EXECUTABLE} ${_fypp} ${_fypp_options} ${_fypp_file_absolute} ${_fypp_target_file}
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                OUTPUT ${_fypp_target_file}
                DEPENDS ${_fypp_file_absolute})
        endforeach()
    endif()


    # Actually add the library to the cmake build
    add_library( ${_p_TARGET} ${_p_TYPE} ${_p_SOURCES} ${${_p_TARGET}_FYPP_SOURCES} ${${_p_TARGET}_TEMPLATED_SOURCES})


    # Add definitions to the compliation
    if (DEFINED _p_DEFINITIONS )
        get_property( _target_defs TARGET ${_p_TARGET} PROPERTY COMPILE_DEFINITIONS )
        list( APPEND _target_defs ${_p_DEFINITIONS} )
        message( STATUS "Library ${_p_TARGET} using definitions: ${_target_defs}" )
        set_property( TARGET ${_p_TARGET} PROPERTY COMPILE_DEFINITIONS ${_target_defs} )
    endif()


    # set the warn-error flags
    if (DEFINED _p_WARNERR )
        if( HAVE_WARNINGS )
        foreach( _lang C CXX Fortran )
            if( CMAKE_${_lang}_COMPILER_LOADED AND DEFINED ${PROJECT_NAME}_${_lang}_WARN_ERROR_FLAG )
                target_compile_options(${_p_TARGET} PRIVATE $<$<COMPILE_LANGUAGE:${_lang}>:${${PROJECT_NAME}_${_lang}_WARN_ERROR_FLAG}>)
            endif()
        endforeach()
        endif()
    endif()

    if( HAVE_WARNINGS )
#         set( ${_p_TARGET}_RELAX_WARNINGS )
        if (_p_RELAX_WARNINGS )
            set_property(SOURCE ${_p_RELAX_WARNINGS} PROPERTY COMPILE_FLAGS ${${PROJECT_NAME}_Fortran_relaxed_WARNING_FLAGS})
        endif()
    endif()


    # Add (private) includes

    if( DEFINED _p_PRIVATE_INCLUDES )
        list( REMOVE_DUPLICATES _p_PRIVATE_INCLUDES )
        foreach( include_path ${_p_PRIVATE_INCLUDES} )
            target_include_directories( ${_p_TARGET} PRIVATE ${include_path} )
        endforeach()
    endif()

    # If we have special compilation flags for F77 files, then add them here

    if( ${PROJECT_NAME}_F77_FLAGS )
        foreach( _file ${_p_SOURCES} )
            get_filename_component( _extension ${_file} EXT )
            if( _extension STREQUAL ".F" )
                set_source_files_properties( ${_file} PROPERTIES COMPILE_FLAGS ${${PROJECT_NAME}_F77_FLAGS} )
            endif()
        endforeach()
    endif()

    # We need c++11 for the molpro plugin

    set_property(TARGET ${_p_TARGET} PROPERTY CXX_STANDARD 11)

    # Add the link libraries

    if ( _p_LIBS )
	    list( REMOVE_DUPLICATES _p_LIBS )
        foreach( _lib ${_p_LIBS} )
            target_link_libraries( ${_p_TARGET} ${_lib} )
        endforeach()
    endif()

    # Have we manually set an output name?

    if ( _p_OUTPUT_NAME )
        message( STATUS "Library ${_p_TARGET}: Output name is ${_p_OUTPUT_NAME}" )
        set_property( TARGET ${_p_TARGET} PROPERTY OUTPUT_NAME ${_p_OUTPUT_NAME} )
    endif()

    # Where do we put the Fortran modules?

    set_property( TARGET ${_p_TARGET} PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/${_p_TARGET} )
    target_include_directories( ${_p_TARGET} PRIVATE ${CMAKE_BINARY_DIR}/modules/${_p_TARGET} )

    # Where do the files get built to

    set_property( TARGET ${_p_TARGET} PROPERTY LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
    set_property( TARGET ${_p_TARGET} PROPERTY ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )

    # Global dependencies

    if ( DEFINED ${PROJECT_NAME}_GLOBAL_DEPENDENCIES AND NOT ${PROJECT_NAME}_GLOBAL_DEPENDENCIES STREQUAL "" )
        add_dependencies( ${_p_TARGET} ${${PROJECT_NAME}_GLOBAL_DEPENDENCIES} )
    endif()

    # Specify the linker language manually

    if( NOT DEFINED _p_LINKER_LANGUAGE OR NOT _p_LINKER_LANGUAGE MATCHES "(C|CXX|Fortran)" )
        message( FATAL_ERROR "LINKER_LANGUAGE not set for library: ${_p_TARGET}" )
    endif()

    set_property( TARGET ${_p_TARGET} PROPERTY LINKER_LANGUAGE ${_p_LINKER_LANGUAGE} )
    message(STATUS "Library ${_p_TARGET}: Setting linker language to ${_p_LINKER_LANGUAGE}" )
    if( DEFINED NECI_${_p_LINKER_LANGUAGE}_${_p_TYPE}_LINK_LIBRARIES )
        target_link_libraries( ${_p_TARGET} ${NECI_${_p_LINKER_LANGUAGE}_${_p_TYPE}_LINK_LIBRARIES} )
        message(STATUS "Library ${_p_TARGET}: Adding link libraries ${NECI_${_p_LINKER_LANGUAGE}_${_p_TYPE}_LINK_LIBRARIES}" )
    endif()
    if( DEFINED NECI_${_p_LINKER_LANGUAGE}_LINK_LIBRARIES )
        target_link_libraries( ${_p_TARGET} ${NECI_${_p_LINKER_LANGUAGE}_LINK_LIBRARIES} )
        message(STATUS "Library ${_p_TARGET}: Adding link libraries ${NECI_${_p_LINKER_LANGUAGE}_LINK_LIBRARIES}" )
    endif()
    if( DEFINED NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS )
        target_link_libraries( ${_p_TARGET} ${NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS} )
        message(STATUS "Library ${_p_TARGET}: Adding linker flags ${NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS}" )
    endif()
    if( DEFINED NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS_${CMAKE_BUILD_TYPE} )
        target_link_libraries( ${_p_TARGET} ${NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS_${CMAKE_BUILD_TYPE}} )
        message(STATUS "Library ${_p_TARGET}: Adding linker flags ${NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS_${CMAKE_BUILD_TYPE}}" )
    endif()

    # Add to the global list of libraries
    #list( APPEND ${PROJECT_NAME}_ALL_LIBS ${_p_TARGET} )
    set( ${PROJECT_NAME}_ALL_LIBS ${${PROJECT_NAME}_ALL_LIBS} ${_p_TARGET} CACHE INTERNAL "" )

endmacro()
