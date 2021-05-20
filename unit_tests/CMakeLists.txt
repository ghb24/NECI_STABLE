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
add_subdirectory( real_space_hubbard )

add_subdirectory( lattice_mod )

add_subdirectory( molecular_tc)

add_subdirectory( back_spawn )

add_subdirectory( back_spawn_excit_gen )

add_subdirectory( ueg_excit_gens )

add_subdirectory( tJ_model )

add_subdirectory( gasci )

add_subdirectory( k_space_hubbard )

add_subdirectory( lattice_models_utils )

add_subdirectory( cc_amplitudes )

add_subdirectory( cepa_shifts )

add_subdirectory( guga )

add_subdirectory( impurity_excit_gen )

add_subdirectory( pcpp_excitgen)

add_subdirectory( gdata_io )

add_subdirectory( test_loop )

add_subdirectory( guga_pchb )

add_subdirectory( orb_idx_and_exc_types )

add_subdirectory( sets_mod )

add_subdirectory( shared_array )

add_subdirectory( lib )
