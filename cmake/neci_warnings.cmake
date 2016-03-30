#
# Add an option to switch on warnings. This should be enabled by default.
# =======================================================================

neci_add_option(
    FEATURE WARNINGS
    DEFAULT ON
    DESCRIPTION "Enable compilation warnings (to the maximal degree)" )

if( HAVE_WARNINGS )

    message( STATUS "Enabling compilation warnings" )

    # Add compiler warnings for each language that they have been defined for.

    foreach( _lang C CXX Fortran )
      if( CMAKE_${_lang}_COMPILER_LOADED )
        if( DEFINED ${PROJECT_NAME}_${_lang}_WARNING_FLAGS AND
            NOT ${PROJECT_NAME}_${_lang}_WARNING_FLAGS STREQUAL "" )
          message("SETTING ${_lang} .. ${CMAKE_${_lang}_FLAGS} .. ${${PROJECT_NAME}_${_lang}_WARNING_FLAGS} " )
          set( CMAKE_${_lang}_FLAGS "${CMAKE_${_lang}_FLAGS} ${${PROJECT_NAME}_${_lang}_WARNING_FLAGS}" )
        endif()
      endif()
    endforeach()

endif()
