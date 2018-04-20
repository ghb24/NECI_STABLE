
# Make some stuff available to all the tests
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# Build FRUIT for unit testing

add_subdirectory( fruit_src )

# And start adding tests
add_subdirectory( det_bit_ops )

# finally get unit testing integrated in software projects 
add_subdirectory( real_space_hubbard )

add_subdirectory( lattice_mod )


add_subdirectory( impurity_excit_gen )

# finally do some unit-testing also with my newly added features 
add_subdirectory( back_spawn )

add_subdirectory( back_spawn_excit_gen ) 

add_subdirectory( ueg_excit_gens )

add_subdirectory( tJ_model )

add_subdirectory( k_space_hubbard )

add_subdirectory( lattice_models_utils )

add_subdirectory( cc_amplitudes )

add_subdirectory( cepa_shifts )
