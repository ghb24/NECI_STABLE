#
# Automatic detection of various clusters
# --> Make automagical modifications to the build configuration, set global variables, and/or enforce that
#     we are using certain toolchain files

# Cluster detection ultimately depends on hostname checking ... imprecise, but hey.
site_name(_hostname)

# Special settings for archer

if( ${_hostname} MATCHES "eslogin[0-9]" )

	message("***************************************************************************************")
	message("Compilation on Archer detected")
	message("-- Archer is a somewhat weird system, with wrapper scripts that defeat autodetection")
	message("--")
	message("-- Ensure that the correct modules are loaded, then:")
	message("-- Use the module cray-hdf5-parallel for HDF5 support")
	message("-- Run cmake with -DCMAKE_TOOLCHAIN_FILE=<neci_dir>/toolchains/archer.cmake")
	message("-- Normal MPI detection is overridden. If the compiler is not specified, this will fail")
	message("***************************************************************************************")

    set( ${PROJECT_NAME}_KNOWN_CLUSTER "Archer" )

    if ( DEFINED CMAKE_TOOLCHAIN_FILE )
        get_filename_compontent( _toolchain_file ${CMAKE_TOOLCHAIN_FILE} NAME )
    endif()

    if ( NOT DEFINED CMAKE_TOOLCHAIN_FILE OR
         NOT _toolchain_file STREQUAL "archer.cmake" )
        message(FATAL_ERROR "Compiling on Archer requires the archer.cmake toolchain file")
    endif()

endif()


if( ${_hostname} MATCHES "hydra[0-9]" )
    
    set( ${PROJECT_NAME}_KNOWN_CLUSTER "Hydra" )

    if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
        message( STATUS "Changing RELEASE optimisation flags (Fortran) as Hydra is an inhomogeneous cluster" )
        string( REPLACE "-xHost" "" _fortran_tmp "${CMAKE_Fortran_FLAGS_RELEASE}" )
        set( CMAKE_Fortran_FLAGS_RELEASE "${_fortran_tmp} -xavx" )
	endif()
    if( CMAKE_C_COMPILER_ID STREQUAL "Intel" )
        message( STATUS "Changing RELEASE optimisation flags (C) as Hydra is an inhomogeneous cluster" )
        string( REPLACE "-xHost" "" _c_tmp "${CMAKE_C_FLAGS_RELEASE}" )
        set( CMAKE_Fortran_C_RELEASE "${_c_tmp} -xavx" )
    endif()
    if( CMAKE_CXX_COMPILER_ID STREQUAL "Intel" )
        message( STATUS "Changing RELEASE optimisation flags (CXX) as Hydra is an inhomogeneous cluster" )
        string( REPLACE "-xHost" "" _cxx_tmp "${CMAKE_CXX_FLAGS_RELEASE}" )
        set( CMAKE_Fortran_CXX_RELEASE "${_cxx_tmp} -xavx" )
    endif()

endif()

set( ${PROJECT_NAME}_BUILD_HOSTNAME ${_hostname} )
mark_as_advanced( _hostname ${PROJECT_NAME}_BUILD_HOSTNAME )

