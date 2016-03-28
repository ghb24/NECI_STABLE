# neci_compiler_flags
# ===================
#
# Ensure that customised compile flags can be used for each of the relevant compilers
#
# i) We can load flags for the specified compilers from complier_flags/Compiler_Language.cmake
#
#     - Flags _OVERRIDE_ those set by default
#     - e.g. NECI_CXX_FLAGS_DEBUG --> overrides CMAKE_CXX_FLAGS_DEBUG
#
# ii) We can specify override flags in any toolchain files
#
# iii) We can also override linker flags both by language (NECI_Fortran_LINK_FLAGS) and by object type
#      (NECI_EXE_LINKER_FLAGS_DEBUG)
#
# n.b. It is also possible to use the NECI_<lang>_LINK_LIBRARIES command to force additional
#      linker libraries to be added in the add_library and add_executable stages.

macro( neci_compiler_flags _lang )

    if( CMAKE_${_lang}_COMPILER_LOADED )

        set( _lang_compiler_defs_file "${CMAKE_CURRENT_LIST_DIR}/compiler_flags/${CMAKE_${_lang}_COMPILER_ID}_${_lang}.cmake" )

        message(STATUS "${_lang} complier: Attempting to load flags from ${_lang_compiler_defs_file}")

        include( ${_lang_compiler_defs_file} OPTIONAL )

    endif()

    # Set the CMAKE flags from these
    foreach( _btype DEBUG RELEASE RELWITHDEBINFO CLUSTER )
      if ( DEFINED NECI_${_lang}_FLAGS_${_btype} )
        message(STATUS "${_lang} compiler: Overriding ${_btype} flags")
        set( CMAKE_${_lang}_FLAGS_${_btype} ${NECI_${_lang}_FLAGS_${_btype}} )
      endif()
      mark_as_advanced( CMAKE_${_lang}_FLAGS_${_btype} )
    endforeach()

    if ( DEFINED NECI_${_lang}_FLAGS )
      message(STATUS "${_lang} compiler: Overriding global flags")
      set( CMAKE_${_lang}_FLAGS ${NECI_${_lang}_FLAGS} )
    endif()
    mark_as_advanced( CMAKE_${_lang}_FLAGS )

    if ( DEFINED NECI_${_lang}_LINK_FLAGS )
      message(STATUS "${_lang} compiler: Overriding global linker flags")
      set( CMAKE_${_lang}_LINK_FLAGS ${NECI_${_lang}_LINK_FLAGS} )
    endif()
    mark_as_advanced( CMAKE_${_lang}_FLAGS )

endmacro()


### OVERRIDE compiler flags. This must be an override because cmake forcibly sets them.
foreach( _lang C CXX Fortran )
  neci_compiler_flags( ${_lang} )
endforeach()

foreach( _btype DEBUG RELEASE RELWITHDEBINFO CLUSTER )
    foreach( _obj EXE SHARED MODULE )
      if ( NECI_${_obj}_LINKER_FLAGS_${_btype} )
        set( CMAKE_${_obj}_LINKER_FLAGS_${_btype} ${NECI_${_obj}_LINKER_FLAGS_${_btype} )
      endif()
    endforeach()
endforeach()
