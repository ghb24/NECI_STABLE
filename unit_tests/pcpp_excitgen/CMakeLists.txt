
foreach( _config ${${PROJECT_NAME}_CONFIGURATIONS} )

    neci_add_test(
        TARGET test_${_config}_aliasTables
        SOURCES test_aliasTables.F90
        LIBS lib${_config} fruit
        LINKER_LANGUAGE Fortran
        META_TARGET test_${_config}
	DEFINITIONS
	    "${${PROJECT_NAME}_${_config}_DEFINITIONS}"
    )

endforeach()
