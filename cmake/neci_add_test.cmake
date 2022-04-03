# Inspired by ecbuild_add_test in ecbuild by the ECMWF

# neci_add_test
# ===================
#
# Add a test with a given list of source files. ::
#
#   neci_add_test(   TARGET <name>
#                    SOURCES <source1> [<source2> ...]
#                    LINKER_LANGUAGE <(C|CXX|Fortran)>
#                  [ MPI <nprocessors> ]
#                  [ LIBS <lib1> [<lib2> ...]
#                  [ META_TARGET <target> ]
#		   [ DEFINITIONS <define1> [<define2> ... ] ]
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
# WARNERR : optional
#   Treat warnings as errrors. Defaults to off.
#
# MPI : optional
#   number of MPI tasks to use.
#
#   If greater than 1, and MPI is not available, the test is disabled
#
# LIBS : optional
#   list of libraries to link against (CMake targets, or external libraries)
#
# META_TARGET : optional
#   add the test to the specified meta target (and create that target if required)
#
# DEFINITIONS : optional
#   list of definitions to add to preprocessor defines
#	same as in neci_add_library, since unit-tests may also depend on build-target!
#       but i cannot get it running yet with the preprocessor running across the file first..
#
# PRIVATE_INCLUDES : optional
#   list of paths to add to include directories that will NOT be exported to other projects. Currently
#   the PUBLIC_INCLUDES functionality is not implemented in this project.
#
##############################################################################

macro( neci_add_test )

    set( options )
    set( single_value_args TARGET LINKER_LANGUAGE MPI META_TARGET WARNERR)
    set( multi_value_args SOURCES LIBS DEFINITIONS PRIVATE_INCLUDES)

    cmake_parse_arguments( _p "${options}" "${single_value_args}" "${multi_value_args}" ${_FIRST_ARG} ${ARGN} )

    if( _p_UNPARSED_ARGUMENTS )
        message( FATAL_ERROR "Unknown keywords given to neci_add_test(): \"${_p_UNPARSED_ARGUMENTS}\"" )
    endif()

    set( _test_dir ${CMAKE_CURRENT_BINARY_DIR}/${_p_TARGET} )

    # TODO: Check for (and use) MPI. Modify _test_command to make uso of MPI if desired
    set( _test_command ${_p_TARGET} )

    set( _test_arguments "" )

    # Build the executable

    add_executable( ${_p_TARGET} ${_p_SOURCES} )
    # make each target reside in its own work directory
    set_target_properties( ${_p_TARGET} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${_test_dir})


    # Add definitions to the compliation
    if (DEFINED _p_DEFINITIONS )
        get_property( _target_defs TARGET ${_p_TARGET} PROPERTY COMPILE_DEFINITIONS )
        list( APPEND _target_defs ${_p_DEFINITIONS} )
        message( DEBUG "Library ${_p_TARGET} using definitions: ${_target_defs}" )
        set_property( TARGET ${_p_TARGET} PROPERTY COMPILE_DEFINITIONS ${_target_defs} )
    endif()

    # Add the link libraries
    if( _p_LIBS )
        list( REMOVE_DUPLICATES _p_LIBS )
        foreach( _lib ${_p_LIBS} )
            target_link_libraries( ${_p_TARGET} ${_lib} )

            # If this library is a Fortran library, we probably need access to its module files...
            get_property( _lib_ll TARGET ${_lib} PROPERTY LINKER_LANGUAGE )
            if( _lib_ll STREQUAL "Fortran" )
                get_property( _lib_modules TARGET ${_lib} PROPERTY Fortran_MODULE_DIRECTORY )
                target_include_directories( ${_p_TARGET} PRIVATE ${_lib_modules} )
            endif()
        endforeach()
    endif()

    # add MPI include dir
    target_include_directories(${_p_TARGET} PRIVATE ${MPI_Fortran_INCLUDE_PATH})

    # Specify the linker language (and additional associated properties)

    if( NOT DEFINED _p_LINKER_LANGUAGE OR NOT _p_LINKER_LANGUAGE MATCHES "(C|CXX|Fortran)" )
      message( FATAL_ERROR "LINKER_LANGUAGE not set for executable: ${_p_TARGET}" )
    endif()

    set_property( TARGET ${_p_TARGET} PROPERTY LINKER_LANGUAGE ${_p_LINKER_LANGUAGE} )
    if( DEFINED NECI_${_p_LINKER_LANGUAGE}_EXE_LINK_LIBRARIES )
      target_link_libraries( ${_p_TARGET} ${NECI_${_p_LINKER_LANGUAGE}_EXE_LINK_LIBRARIES} )
      message(DEBUG "Executable ${_p_TARGET}: Adding link libraries ${NECI_${_p_LINKER_LANGUAGE}_LINK_LIBRARIES}" )
    endif()
    if( DEFINED NECI_${_p_LINKER_LANGUAGE}_LINK_LIBRARIES )
      target_link_libraries( ${_p_TARGET} ${NECI_${_p_LINKER_LANGUAGE}_LINK_LIBRARIES} )
      message(DEBUG "Executable ${_p_TARGET}: Adding link libraries ${NECI_${_p_LINKER_LANGUAGE}_LINK_LIBRARIES}" )
    endif()
    if( DEFINED NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS )
      target_link_libraries( ${_p_TARGET} ${NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS} )
      message(DEBUG "Executable ${_p_TARGET}: Adding linker flags ${NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS}" )
    endif()
    if( DEFINED NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS_${CMAKE_BUILD_TYPE} )
      target_link_libraries( ${_p_TARGET} ${NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS_${CMAKE_BUILD_TYPE}} )
      message(DEBUG "Library ${_p_TARGET}: Adding linker flags ${NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS_${CMAKE_BUILD_TYPE}}" )
    endif()

    # Add (private) includes
    if( DEFINED _p_PRIVATE_INCLUDES )
        list( REMOVE_DUPLICATES _p_PRIVATE_INCLUDES )
        foreach( include_path ${_p_PRIVATE_INCLUDES} )
            target_include_directories( ${_p_TARGET} PRIVATE ${include_path} )
        endforeach()
    endif()

    # set the warn-error flags
    if (DEFINED _p_WARNERR )
        if(HAVE_WARNINGS )
        foreach( _lang C CXX Fortran )
            if( CMAKE_${_lang}_COMPILER_LOADED AND DEFINED ${PROJECT_NAME}_${_lang}_WARN_ERROR_FLAG )
                target_compile_options(${_p_TARGET} PRIVATE $<$<COMPILE_LANGUAGE:${_lang}>:${${PROJECT_NAME}_${_lang}_WARN_ERROR_FLAG}>)
            endif()
        endforeach()
        endif()
    endif()

    # Add the testing infrastructure to the executable

    add_test( NAME ${_p_TARGET} COMMAND ${_test_command} "${_test_arguments}" ${_test_dir} )

    set_tests_properties( ${_p_TARGET} PROPERTIES WORKING_DIRECTORY ${_test_dir} )


    message(DEBUG "Added unit test: ${_p_TARGET}" )

    # Add to the global list of tests

    set( ${PROJECT_NAME}_ALL_TESTS ${${PROJECT_NAME}_ALL_TESTS} ${_p_TARGET} CACHE INTERNAL "" )

    # Add to the metatarget if necessary

    if( _p_META_TARGET  )
        list(FIND ${PROJECT_NAME}_ALL_META_TARGETS ${_p_META_TARGET} _idx)
        if( NOT DEFINED _idx OR "${_idx}" EQUAL -1 )
            add_custom_target(${_p_META_TARGET})
            set( ${PROJECT_NAME}_ALL_META_TARGETS ${${PROJECT_NAME}_ALL_META_TARGETS} ${_p_META_TARGET} CACHE INTERNAL "")
        endif()
        add_dependencies( ${_p_META_TARGET} ${_p_TARGET} )
        message(DEBUG "Added test ${_p_TARGET} to meta target ${_p_META_TARGET}")
    endif()

endmacro()
