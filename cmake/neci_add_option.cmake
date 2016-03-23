# Inspired by ecbuild_add_option in ecbuild by the ECMWF

# neci_add_option
# ===============
#
# Add a CMake configuration option
#
#   neci_add_option(   FEATURE <name>
#                    [ DEFAULT ON|OFF ]
#                    [ DESCRIPTION <description> ])
#                    [ REQUIRED_PACKAGES <package1> [<package2> ...] ]
#
# Options
# -------
#
# FEATURE : required
#   name of the feature / option
#
# DEFAULT : optional, defaults to ON
#  if set to ON, the feature is enabled even if not explicitly requested
#
# DESCRIPTION : optional
#  string describing the feature (shown in summary and stored in the cache)
#
# REQUIRED_PACKAGES : optional
#  list of packages required to be found for this feature to be enabled.
#
# Usage
# -----
#
# Features should be enabled with ``-DENABLE_<feature>=[ON|OFF]``. If they are enabled with REQUIRE
# then an error will occur if the package is not found.
#
# Results in the cmake variable ``HAVE_<feature>`` being set to ON
#
# Unlike eckit, we don't yet have the ability to automagically turn a feature off if a package
# is not available.
# TODO: Turn off options when a package is not available

macro( neci_add_option )

    # set( options ADVANCED ) 
    set( single_value_args FEATURE DEFAULT DESCRIPTION )
    set( multi_value_args  REQUIRED_PACKAGES ) # CONDITION

    cmake_parse_arguments( _p "" "${single_value_args}" "${multi_value_args}" ${_FIRST_ARG} ${ARGN} )

    if( _p_UNPARSED_ARGUMENTS )
      message(FATAL_ERROR "Unknown arguments passed to neci_add_option(): \"${_p_UNPARSED_ARGUMENTS}\"")
    endif()
  
    # Check (rquired) feature parameter
    if( NOT _p_FEATURE )
      message(FATAL_ERROR "The call to neci_add_option() doesn't specify the FEATURE.")
    endif()

    # Check DEFAULT parameter

    if( NOT DEFINED _p_DEFAULT )
      set( _p_DEFAULT ON )
    else()
      if( NOT _p_DEFAULT MATCHES "[Oo][Nn]" AND NOT _p_DEFAULT MATCHES "[Oo][Ff][Ff]" )
        message(FATAL_ERROR "In macro neci_add_option(), DEFAULT is either ON or OFF: \"${_p_DEFAULT}\"")
      endif()
    endif()

    message( STATUS "Option ENABLE_${_p_FEATURE} defaults to ${_p_DEFAULT}" )

    # check if user provided value

    get_property( _in_cache CACHE ENABLE_${_p_FEATURE} PROPERTY VALUE )

    if ( ENABLE_${_p_FEATURE} MATCHES "REQUIRE" )
      set( ENABLE_${_p_FEATURE} ON CACHE BOOL "" FORCE )
      message( STATUS "Option ENABLE_${_p_FEATURE} was required" )
      set( ${_p_FEATURE}_user_provided_input 1 CACHE BOOL "" FORCE )
    elseif ( NOT ENABLE_${_p_FEATURE} STREQUAL "" AND _in_cache )
      message( STATUS "Option ENABLE_${_p_FEATURE} was found in the cache" )
      set( ${_p_FEATURE}_user_provided_input 1 CACHE BOOL "" )
    else()
      message( STATUS "Option ENABLE_${_p_FEATURE} was not found in the cache" )
      set( ${_p_FEATURE}_user_provided_input 0 CACHE BOOL "" )
    endif()

    mark_as_advanced( ${_p_FEATURE}_user_provided_input )

    # And set a parameter to determine if we have it. Note that we change the name to HAVE_<feature>
    # as ENABLE_<feature> is for the user provided option, and we want to reserve the right to have
    # some logic between these values

    if ( ${_p_FEATURE}_user_provided_input )
      set ( HAVE_${_p_FEATURE} ${ENABLE_${_p_FEATURE}} )
    else()
      set ( HAVE_${_p_FEATURE} ${_p_DEFAULT} )
    endif()

    # If we want to enable a package, then check that its required packages exist. If they do not,
    # then we need to disable the package (unless REQUIRED is set, in which case we return an error.

    if ( HAVE_${_p_FEATURE} )
        
      set( ${_p_FEATURE}_packages_found ON )
      set( ${_p_FEATURE}_failed_list "" )

      foreach( pkg ${_p_REQUIRED_PACKAGES} )
        message( STATUS "Option ENABLE_${_p_FEATURE}: Searching for dependent package: ${pkg}" )

        find_package( ${pkg} )
        if( NOT ${pkg}_FOUND )
          set( ${_p_FEATURE}_packages_found OFF )
          list(APPEND ${_p_FEATURE}_failed_list ${pkg} )
        endif()
      endforeach()

      mark_as_advanced( ${_p_FEATURE}_packages_found )
      if ( NOT ${_p_FEATURE}_packages_found )
        message(STATUS "Option ENABLE_${_p_FEATURE}: Dependent packages missing: ${${_p_FEATURE}_failed_list}" )
        if ( ${_p_FEATURE}_user_provided_input )
          message( FATAL_ERROR "Option ENABLE_${_p_FEATURE} has been required." )
        else()
          set( HAVE_${_p_FEATURE} OFF )
        endif()
      endif()

    endif()
    
    # And finally some pretty output.

    if ( HAVE_${_p_FEATURE} )
      message( STATUS "Feature ${_p_FEATURE} enabled." )
    else()
      message( STATUS "Feature ${_p_FEATURE} disabled." )
    endif()

endmacro()
