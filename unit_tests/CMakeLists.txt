macro (add_unit_test test_name)
      foreach( _config ${${PROJECT_NAME}_CONFIGURATIONS} )

          neci_add_test(
    	    TARGET test_${_config}_${test_name}
	    SOURCES test_${test_name}.F90
	    LIBS lib${_config} fruit
	    LINKER_LANGUAGE Fortran
	    META_TARGET test_${_config}
	    DEFINITIONS
	    	    "${${PROJECT_NAME}_${_config}_DEFINITIONS}"
	)

	endforeach()
endmacro()

# Make some stuff available to all the tests
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# Build FRUIT for unit testing

add_subdirectory( fruit_src )

# And start adding tests
add_subdirectory( det_bit_ops )

# finally get unit testing integrated in software projects 
#add_subdirectory( real_space_hubbard )

#add_subdirectory( lattice_mod )


add_subdirectory( impurity_excit_gen )

add_subdirectory( molecular_tc)

# finally do some unit-testing also with my newly added features 
add_subdirectory( back_spawn )

add_subdirectory( back_spawn_excit_gen ) 

add_subdirectory( umathash )

add_subdirectory( ueg_excit_gens )

#add_subdirectory( tJ_model )

#add_subdirectory( k_space_hubbard )

add_subdirectory( lattice_models_utils )

add_subdirectory( cc_amplitudes )

add_subdirectory( cepa_shifts )

add_subdirectory( analyse_wf_symmetry )

add_subdirectory( pcpp_excitgen)

