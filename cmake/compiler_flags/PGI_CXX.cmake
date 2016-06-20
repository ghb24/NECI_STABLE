
# For some reason, the CMake auto-detection often picks up pgcpp as the PGI c++ compiler, but that is
# really just the preprocessor. Replace it with the actual compiler
get_filename_component( _cxx_nm ${CMAKE_CXX_COMPILER} NAME )
if ( _cxx_nm STREQUAL "pgcpp" )
    message( STATUS "Forcing the use of pgc++ not pgcc" )
    set( CMAKE_CXX_COMPILER pgc++ )
endif()

# Special defines for the PGI C++ compiler

set( ${PROJECT_NAME}_CXX_FLAGS "-mcmodel=medium" )
set( ${PROJECT_NAME}_CXX_FLAGS_CLUSTER "-Mipa=fast" )

# Warning flags

set( ${PROJECT_NAME}_CXX_WARNING_FLAGS "-Minform=warn" )

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_CXX_FLAGS "-pc=32" )

set( ${PROJECT_NAME}_64BIT_CXX_FLAGS "-pc=64" )
