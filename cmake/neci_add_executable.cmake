# Inspired by ecbuild_add_executable in ecbuild by the ECMWF

# neci_add_executable
# ===================
#
# Add an executable with a given list of source files. ::
#
#   neci_add_executable(   TARGET <name>
#                          SOURCES <source1> [<source2> ...]
#                          LINKER_LANGUAGE <(C|CXX|Fortran)>
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
# LINKER_LANGUAGE : required
#   what language should be used to link the executable. If not specified, cmake will determine this automatically
#   Values: C, CXX, Fortran
#
#   required so that <project>_<lang>_linker_flags can be made to work
#
# LIBS : optional
#   list of libraries to link against (CMake targets, or external libraries)
#
##############################################################################


macro( neci_add_executable )

    set( options  )
    set( single_value_args TARGET LINKER_LANGUAGE )
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

    message(DEBUG "Adding executable: ${_p_TARGET}")
    add_executable( ${_p_TARGET} ${_p_SOURCES} )

    # Add the link libraries

    if ( _p_LIBS )
	  list( REMOVE_DUPLICATES _p_LIBS )
      foreach( _lib ${_p_LIBS} )
		target_link_libraries( ${_p_TARGET} ${_lib} )
      endforeach()
    endif()

    # Where do we put the Fortran modules?

    set_property( TARGET ${_p_TARGET} PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/${_p_TARGET} )

    # Where do the files get built to

    set_property( TARGET ${_p_TARGET} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )

    # Molpro plugin requires c++11

    set_property( TARGET ${_p_TARGET} PROPERTY CXX_STANDARD 11)

    # Global dependencies

    if ( DEFINED ${PROJECT_NAME}_GLOBAL_DEPENDENCIES AND NOT ${PROJECT_NAME}_GLOBAL_DEPENDENCIES STREQUAL "" )
        add_dependencies( ${_p_TARGET} ${${PROJECT_NAME}_GLOBAL_DEPENDENCIES} )
    endif()

    # Specify the linker language manually

    if( NOT DEFINED _p_LINKER_LANGUAGE OR NOT _p_LINKER_LANGUAGE MATCHES "(C|CXX|Fortran)" )
      message( FATAL_ERROR "LINKER_LANGUAGE not set for executable: ${_p_TARGET}" )
    endif()

    set_property( TARGET ${_p_TARGET} PROPERTY LINKER_LANGUAGE ${_p_LINKER_LANGUAGE} )
    message(DEBUG "Executable ${_p_TARGET}: Setting linker language to ${_p_LINKER_LANGUAGE}" )
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

    # Add to the global list of libraries
    set( ${PROJECT_NAME}_ALL_EXES ${${PROJECT_NAME}_ALL_EXES} ${_p_TARGET} CACHE INTERNAL "" )

endmacro()
