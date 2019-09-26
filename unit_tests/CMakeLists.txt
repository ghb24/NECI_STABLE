macro (add_unit_test test_name)

          neci_add_test(
    	    TARGET test_neci_${test_name}
	    SOURCES test_${test_name}.F90
	    LIBS libneci fruit
	    LINKER_LANGUAGE Fortran
	    META_TARGET test_neci
	    DEFINITIONS
	    	    "${${PROJECT_NAME}_neci_DEFINITIONS}"
	)

endmacro()

# Make some stuff available to all the tests
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# Build FRUIT for unit testing

add_subdirectory( fruit_src )

# And start adding tests
add_subdirectory( det_bit_ops )

# finally get unit testing integrated in software projects 
add_subdirectory( molecular_tc)

# finally do some unit-testing also with my newly added features 
add_subdirectory( back_spawn )

add_subdirectory( back_spawn_excit_gen ) 

add_subdirectory( ueg_excit_gens )

#add_subdirectory( tJ_model )

#add_subdirectory( k_space_hubbard )

add_subdirectory( lattice_models_utils )

add_subdirectory( cc_amplitudes )

add_subdirectory( cepa_shifts )

add_subdirectory( analyse_wf_symmetry )

add_subdirectory( pcpp_excitgen)

