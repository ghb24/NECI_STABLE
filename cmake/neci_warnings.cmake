#
# Add an option to switch on warnings. This should be enabled by default.
# =======================================================================

neci_add_option(
    FEATURE WARNINGS
    DEFAULT OFF
    DESCRIPTION "Enable compilation warnings (to the maximal degree)" )

if ( HAVE_WARNINGS )

    message( STATUS "Enabling compilation warnings" )

    # Allow workarounds for false-positive warnings.
    list(APPEND NECI_GLOBAL_DEFINES _WARNING_WORKAROUND_)

    # Add compiler warnings for each language that they have been defined for.
    foreach( _lang C CXX Fortran )
      if( CMAKE_${_lang}_COMPILER_LOADED )
        if( DEFINED ${PROJECT_NAME}_${_lang}_WARNING_FLAGS AND
            NOT ${PROJECT_NAME}_${_lang}_WARNING_FLAGS STREQUAL "" )
          set( CMAKE_${_lang}_FLAGS "${CMAKE_${_lang}_FLAGS} ${${PROJECT_NAME}_${_lang}_WARNING_FLAGS}" )
        endif()
      endif()
    endforeach()
endif()


neci_add_option(
    FEATURE WARN_ERROR
    DEFAULT OFF
    DESCRIPTION "Treat warnings as error.")

if ( HAVE_WARN_ERROR)
    message( STATUS "Treat Warnings as errors." )
    if( NOT HAVE_WARNINGS )
        message(FATAL_ERROR "ENABLE_WARN_ERROR=ON requires ENABLE_WARNINGS=ON.")
    endif()
    foreach( _lang C CXX Fortran )
        if( CMAKE_${_lang}_COMPILER_LOADED AND DEFINED ${PROJECT_NAME}_${_lang}_WARN_ERROR_FLAG )
            set( CMAKE_${_lang}_FLAGS "${CMAKE_${_lang}_FLAGS} ${${PROJECT_NAME}_${_lang}_WARN_ERROR_FLAG}" )
        endif()
    endforeach()

endif()
