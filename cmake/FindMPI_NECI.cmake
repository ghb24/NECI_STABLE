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
    if ( MPI_FOUND )
        set(MPI_NECI_LIBRARIES ${MPI_LIBRARIES})
        set(MPI_NECI_Fortran_LIBRARIES ${MPI_Fortran_LIBRARIES})
        set(MPI_NECI_INCLUDE_PATH ${MPI_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
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
