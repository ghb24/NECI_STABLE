# - Find FFTW
# Find the native FFTW includes and library
#
# If FFTW_PATH is specified, this will be used as the search location
#
# FFTW_INCLUDES - where to find fftw3.h
# FFTW_LIBRARIES - List of libraries when using FFTW.
# FFTW_FOUND - True if FFTW found.

if ( NOT FFTW_FOUND )

    if( DEFINED FFTW_PATH )
        find_path (FFTW_INCLUDES fftw3.h PATHS ${FFTW_PATH} ${FFTW_PATH}/include NO_DEFAULT_PATH)
        find_library (FFTW_LIBRARIES NAMES fftw3 PATHS ${FFTW_PATH} ${FFTW_PATH}/lib NO_DEFAULT_PATH)
    endif()

    find_path (FFTW_INCLUDES fftw3.h)

    find_library (FFTW_LIBRARIES NAMES fftw3)

    # handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
    # all listed variables are TRUE
    include( FindPackageHandleStandardArgs )
    find_package_handle_standard_args( FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES )

    mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)

endif()
