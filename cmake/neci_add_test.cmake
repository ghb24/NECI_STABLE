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
##############################################################################

macro( neci_add_test )

    set( options )
    set( single_value_args TARGET LINKER_LANGUAGE MPI META_TARGET )
    set( multi_value_args SOURCES LIBS )

    cmake_parse_arguments( _p "${options}" "${single_value_args}" "${multi_value_args}" ${_FIRST_ARG} ${ARGN} )

    if( _p_UNPARSED_ARGUMENTS )
        message( FATAL_ERROR "Unknown keywords given to neci_add_test(): \"${_p_UNPARSED_ARGUMENTS}\"" )
    endif()

    set( _test_dir ${CMAKE_CURRENT_BINARY_DIR} )

    # TODO: Check for (and use) MPI. Modify _test_command to make uso of MPI if desired
    set( _test_command ${_p_TARGET} )

    set( _test_arguments "" )

    # Build the executable

    add_executable( ${_p_TARGET} ${_p_SOURCES} )

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

    # Specify the linker language (and additional associated properties)

    if( NOT DEFINED _p_LINKER_LANGUAGE OR NOT _p_LINKER_LANGUAGE MATCHES "(C|CXX|Fortran)" )
      message( FATAL_ERROR "LINKER_LANGUAGE not set for executable: ${_p_TARGET}" )
    endif()

    set_property( TARGET ${_p_TARGET} PROPERTY LINKER_LANGUAGE ${_p_LINKER_LANGUAGE} )
    if( DEFINED NECI_${_p_LINKER_LANGUAGE}_EXE_LINK_LIBRARIES )
      target_link_libraries( ${_p_TARGET} ${NECI_${_p_LINKER_LANGUAGE}_EXE_LINK_LIBRARIES} )
      message(STATUS "Executable ${_p_TARGET}: Adding link libraries ${NECI_${_p_LINKER_LANGUAGE}_LINK_LIBRARIES}" )
    endif()
    if( DEFINED NECI_${_p_LINKER_LANGUAGE}_LINK_LIBRARIES )
      target_link_libraries( ${_p_TARGET} ${NECI_${_p_LINKER_LANGUAGE}_LINK_LIBRARIES} )
      message(STATUS "Executable ${_p_TARGET}: Adding link libraries ${NECI_${_p_LINKER_LANGUAGE}_LINK_LIBRARIES}" )
    endif()
    if( DEFINED NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS )
      target_link_libraries( ${_p_TARGET} ${NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS} )
      message(STATUS "Executable ${_p_TARGET}: Adding linker flags ${NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS}" )
    endif()
    if( DEFINED NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS_${CMAKE_BUILD_TYPE} )
      target_link_libraries( ${_p_TARGET} ${NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS_${CMAKE_BUILD_TYPE}} )
      message(STATUS "Library ${_p_TARGET}: Adding linker flags ${NECI_${_p_LINKER_LANGUAGE}_LINKER_FLAGS_${CMAKE_BUILD_TYPE}}" )
    endif()

    # Add the testing infrastructure to the executable

    add_test( ${_p_TARGET} ${_test_command} "${_test_arguments}" ${_test_dir} )

    set_tests_properties( ${_p_TARGET} PROPERTIES WORKING_DIRECTORY ${_test_dir} )

    message(STATUS "Added unit test: ${_p_TARGET}" )

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
        message(STATUS "Added test ${_p_TARGET} to meta target ${_p_META_TARGET}")
    endif()

endmacro()
