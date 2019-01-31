# Instructs CMake to find the various different MPI packages that are available, including the logic that
# we have inherited from the older NECI cmake files.
#
# Defines:
# MPI_NECI_FOUND - True if MPI is found
# MPI_NECI_LIBRARIES - List of libraries
# MPI_NECI_Fortran_LIBRARIES - libraries for linking with Fortran
# MPI_NECI_INCLUDE_PATH - Where to look for headers

if ( NOT MPI_NECI_FOUND )

    find_package( MPI )

    # Map the output of the find to the MPI_NECI finder
    set( MPI_NECI_FOUND ${MPI_FOUND})

    # If we are using ifort, and we have not found the mpiifort wrapper, then it is normally a good
    # idea to force things. Otherwise FindMPI can end up finding the GNU stuff, especially when using
    # mpi.ibm
    #
    # This HACK can be easily disabled in a toolchain file, or on the commandline, by setting the variable
    # ${PROJECT_NAME}_SIMPLE_MPI_SEARCH
    if ( MPI_FOUND AND CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
        if ( NOT DEFINED ENABLE_INTEL_MPI_OVERRIDE OR ENABLE_INTEL_MPI_OVERRIDE )
            get_filename_component( _wrapper_nm ${MPI_Fortran_COMPILER} NAME )
            if( NOT _wrapper_nm STREQUAL "mpiifort" )
                message(STATUS "=======================================================================")
                message(STATUS "")
                message(STATUS "  Intel MPI wrapper override: FORCING the use of mpiifort, mpiicc, mpiicpc")
                message(STATUS "")
                message(STATUS "  If this is not what you want disable by running CMake with -DENABLE_INTEL_MPI_OVERRIDE=OFF ")
                message(STATUS "  Alternatively set the ENABLE_INTEL_MPI_OVERRIDE variable to OFF in a toolchain file.")
                message(STATUS "")
                message(STATUS "=======================================================================")
                set( MPI_Fortran_COMPILER mpiifort )
                set( MPI_C_COMPILER mpiicc )
                set( MPI_CXX_COMPILER mpiicpc )
                find_package( MPI )
                if ( MPI_FOUND )
					set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} CACHE STRING "" FORCE)
					set( CMAKE_C_COMPILER ${MPI_C_COMPILER} CACHE STRING "" FORCE)
					set( CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE STRING "" FORCE)
                else()
                    message( WARNING "No appropriate MPI wrapper found for ifort" )
                endif()
            endif()
        endif()
    endif()

    if( MPI_FOUND)
        set(MPI_NECI_LIBRARIES ${MPI_LIBRARIES})
        set(MPI_NECI_Fortran_LIBRARIES ${MPI_Fortran_LIBRARIES})
        set(MPI_NECI_INCLUDE_PATH ${MPI_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
	set(HAVE_MPI ON)
        if (NOT MPI_NECI_FIND_QUIETLY )
            message(STATUS "MPI found")
        endif()
    else()
        if ( MPI_NECI_FIND_REQUIRED )
            message(FATAL_ERROR "Package MPI required, but not found")
        else()
            if (NOT MPI_NECI_FIND_QUIETLY )
                message(STATUS "Package MPI not found")
            endif()
        endif()
    endif()

endif()
