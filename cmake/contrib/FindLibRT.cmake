#
# Try to find librt library
#
# Once done this will define
# LIBRT_FOUND
# LIBRT_LIBRARIES
# Unecessary duplicates
# LibRT_FOUND
# LibRT_LIBRARIES
#

#if ( NOT LIBRT_FOUND )
#
#    find_library(LIBRT_LIBRARY rt)
#
#    set(LIBRT_LIBRARIES ${LIBRT_LIBRARY})
#
#    include(FindPackageHandleStandardArgs)
#    find_package_handle_standard_args(LibRT DEFAULT_MSG LIBRT_LIBRARY)
#
#    mark_as_advanced( LIBRT_LIBRARY LIBRT_LIBRARIES )
#
#endif()

# Instead of using a normal finder, rely on the fact that this is part of the POSIX standard. As a result:
#
#  i) It _will_ b there
# ii) The finders may behave oddly

if( NOT LibRT_FOUND )

    set( LibRT_FOUND ON )

    if( UNIX AND NOT APPLE )
        set( LIBRT_LIBRARIES rt )
    else()
        set( LIBRT_LIBRARIES "" )
    endif()
    set( LibRT_LIBRARIES ${LIBRT_LIBRARIES} )

    mark_as_advanced( LIBRT_LIBRARIES LibRT_LIBRARIES )

endif()
set(LIBRT_FOUND ${LibRT_FOUND})
