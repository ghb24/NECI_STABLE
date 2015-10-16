#
# Try to find librt library
#
# Once done this will define
#  LIBRT_FOUND
#  LIBRT_LIBRARIES
#

find_library(LIBRT_LIBRARY rt)

set(LIBRT_LIBRARIES ${LIBRT_LIBRARY})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(LibRT DEFAULT_MSG LIBRT_LIBRARY)
