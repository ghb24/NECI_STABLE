
foreach( _config ${${PROJECT_NAME}_CONFIGURATIONS} )

    # this is a bit cumbersome.. the whole unit test setup can i optimize that 
    # or template it? i guess and try that! 
    neci_add_test(
        TARGET test_${_config}_k_space_hubbard
        SOURCES test_k_space_hubbard.F90
        LIBS lib${_config} fruit
        LINKER_LANGUAGE Fortran
        META_TARGET test_${_config}
	DEFINITIONS
	    "${${PROJECT_NAME}_${_config}_DEFINITIONS}"
    )




endforeach()
