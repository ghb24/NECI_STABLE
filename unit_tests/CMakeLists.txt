macro (add_unit_test _config test_name)
        neci_add_test(
    	    TARGET test_${_config}_${test_name}
	    SOURCES
            test_${test_name}.F90
	    LIBS
            lib${_config} fruit
	    LINKER_LANGUAGE
            Fortran
	    META_TARGET
            test_${_config}
        PRIVATE_INCLUDES
            "${MPI_Fortran_INCLUDE_PATH}"
            "${HDF5_NECI_INCLUDE_PATH}"
	    DEFINITIONS
            "${NECI_GLOBAL_DEFINES}"
	    	"${${PROJECT_NAME}_${_config}_DEFINITIONS}"
        WARNERR
            ${${PROJECT_NAME}_${_config}_WARNERR}
	)
endmacro()

macro (add_unit_test_all_configs test_name)
    foreach( _config ${${PROJECT_NAME}_CONFIGURATIONS} )
        add_unit_test( ${_config} ${test_name} )
    endforeach()
endmacro()

# Make some stuff available to all the tests
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# Build FRUIT for unit testing

add_subdirectory( fruit_src )

# And start adding tests
# finally get unit testing integrated in software projects

add_subdirectory( det_bit_ops )

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

add_subdirectory( SD_spin_purification )

add_subdirectory( lib )

add_subdirectory( input_parser )

add_subdirectory( util_mod )
