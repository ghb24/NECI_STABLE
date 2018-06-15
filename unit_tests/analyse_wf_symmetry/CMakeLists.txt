
foreach( _config ${${PROJECT_NAME}_CONFIGURATIONS} )

    # this is a bit cumbersome.. the whole unit test setup can i optimize that 
    # or template it? i guess and try that! 
    neci_add_test(
    	TARGET test_${_config}_analyse_wf_symmetry.F90
        SOURCES test_analyse_wf_symmetry.F90
        LIBS lib${_config} fruit
        LINKER_LANGUAGE Fortran
        META_TARGET test_${_config}
	DEFINITIONS
	    "${${PROJECT_NAME}_${_config}_DEFINITIONS}"
    )




endforeach()
