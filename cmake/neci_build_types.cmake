#
# Ensure that all of the correct build types have been defined, and that their flags get combined correctly.
#

# i) Test debug mode in here (so that there is a well defined debug definitions test)
# ii) If debug mode is on, we should automagically add warnings.

set( CMAKE_CONFIGURATION_TYPES DEBUG RELEASE CLUSTER CACHE TYPE INTERNAL FORCE )

# We want to have a default build type!
if( NOT DEFINED CMAKE_BUILD_TYPE OR NOT CMAKE_BUILD_TYPE )
    message( STATUS "Build type not specified. Using RELEASE as default" )
    set( CMAKE_BUILD_TYPE "RELEASE" CACHE STRING "Type of build, options are ${CMAKE_CONFIGURATION_TYPES}" FORCE )
endif()

# Enforce consistent naming standards (all upper case)
string( TOUPPER ${CMAKE_BUILD_TYPE} _build_type )

# No easy way to test if value in array!?! (until cmake 3.0.3)
set( _valid_type OFF )
foreach( _ty ${CMAKE_CONFIGURATION_TYPES} )
    if ( _build_type STREQUAL ${_ty} )
        set( _valid_type ON)
    endif()
endforeach()

if( NOT _valid_type )
    message( FATAL_ERROR "Specified build type not valid. Choices are: ${CMAKE_CONFIGURATION_TYPES}" )
endif()
set( CMAKE_BUILD_TYPE ${_build_type} )

