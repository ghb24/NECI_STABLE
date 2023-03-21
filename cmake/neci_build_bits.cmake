# neci_build_bits
# ===============
#
# Behaviour may, or may not, depend on the size of a pointer (32 or 64 bits), and the subsequent
# behaviour based on that

if ( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    set( _bits 64 )
else()
    set( _bits 32 )
endif()
message(STATUS "Creating a ${_bits}-bit build")

# Create globally available macros for accessing status

set( ${PROJECT_NAME}_32BIT_BUILD OFF )
set( ${PROJECT_NAME}_64BIT_BUILD OFF )
set( ${PROJECT_NAME}_${_bits}BIT_BUILD ON )

# If we have defined 32/64 bit compiler definitions, add them to the defines. Otherwise use defaults.

# Update bit-specific compile flags

foreach( _type C CXX Fortran )

    # Provide some defaults that work with almost all compilers
    if ( NOT DEFINED ${PROJECT_NAME}_${_bits}BIT_${_type}_FLAGS )
        set( ${PROJECT_NAME}_${_bits}BIT_${_type}_FLAGS "-m${_bits}" )
    endif()

    if ( DEFINED ${PROJECT_NAME}_${_bits}BIT_FLAGS AND NOT ${PROJECT_NAME}_${_bits}BIT_FLAGS STREQUAL "" )
        message(STATUS "Adding ${_bits}-bit flags for ${_type}: ${${PROJECT_NAME}_${_bits}BIT_FLAGS}" )
        set( CMAKE_${_type}_FLAGS "${CMAKE_${_type}_FLAGS} ${${PROJECT_NAME}_${_bits}BIT_FLAGS}" )
    endif()

    if ( DEFINED ${PROJECT_NAME}_${_bits}BIT_${_type}_FLAGS )
        message(STATUS "Adding ${_bits}-bit flags for ${_type}: ${${PROJECT_NAME}_${_bits}BIT_${_type}_FLAGS}" )
        set( CMAKE_${_type}_FLAGS "${CMAKE_${_type}_FLAGS} ${${PROJECT_NAME}_${_bits}BIT_${_type}_FLAGS}" )
    endif()
endforeach()

# And the same for linker flags

# Provide some defaults that work with almost all compilers
#if ( NOT DEFINED ${PROJECT_NAME}_${_bits}BIT_LINKER_FLAGS )
#    set( ${PROJECT_NAME}_${_bits}BIT_LINKER_FLAGS "-m${_bits}" )
#endif()

foreach( _type C EXE SHARED STATIC MODULE)

    if ( DEFINED ${PROJECT_NAME}_${_bits}BIT_LINKER_FLAGS AND
            NOT ${PROJECT_NAME}_${_bits}BIT_LINKER_FLAGS STREQUAL "" )
        message(STATUS "Adding ${_bits}-bit linker flags for ${_type}: ${${PROJECT_NAME}_${_bits}BIT_LINKER_FLAGS}" )
        set( CMAKE_${_type}_LINKER_FLAGS "${CMAKE_${_type}_LINKER_FLAGS} ${${PROJECT_NAME}_${_bits}BIT_LINKER_FLAGS}" )
    endif()

    if ( DEFINED ${PROJECT_NAME}_${_bits}BIT_${_type}_LINKER_FLAGS )
        message(STATUS "Adding ${_bits}-bit flags for ${_type}: ${${PROJECT_NAME}_${_bits}BIT_${_type}_LINKER_FLAGS}" )
        set( CMAKE_${_type}_LINKER_FLAGS "${CMAKE_${_type}_LINKER_FLAGS} ${${PROJECT_NAME}_${_bits}BIT_${_type}_LINKER_FLAGS}" )
    endif()
endforeach()
