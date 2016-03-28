# Instructs CMake to find the various different packages that support BLAS/LAPACK. A lot of this work is
# done by the default finders, but insufficient to be truly robust across systems.
#
# Defines:
# LAPACK_NECI_FOUND - True if BLAS/LAPACK provider is found
# LAPACK_NECI_LIBRARIES - List of libraries
# LAPACK_NECI_INCLUDE_PATH - Where to look for headers

if ( NOT LAPACK_NECI_FOUND )

    find_package( MKL )
    if( MKL_FOUND )
        set( LAPACK_NECI_FOUND ${MKL_FOUND} )
        set( LAPACK_NECI_LIBRARIES ${MKL_LIBRARIES} )
        set( LAPACK_NECI_INCLUDE_PATH ${MKL_INCLUDE_DIR} )
    else()
        find_package( LAPACK )
        set( LAPACK_NECI_FOUND ${LAPACK_FOUND} )
        if( LAPACK_FOUND )
            set( LAPACK_NECI_LIBRARIES ${LAPACK_LIBRARIES} )
        endif()
    endif()


    # Map the output of the find to the LAPACK_NECI finder
    if ( LAPACK_NECI_FOUND )
        if (NOT LAPACK_NECI_FIND_QUIETLY )
            message(STATUS "BLAS/LAPACK provider found")
        endif()
    else()
        if ( LAPACK_NECI_FIND_REQUIRED )
            message(FATAL_ERROR "BLAS/LAPACK provider required, but not found")
        else()
            if (NOT LAPACK_NECI_FIND_QUIETLY )
                message(STATUS "BLAS/LAPACK provider not found")
            endif()
        endif()
    endif()

endif()

