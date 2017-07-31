
# Make some stuff available to all the tests
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# Build FRUIT for unit testing

add_subdirectory( fruit_src )

# And start adding tests
add_subdirectory( det_bit_ops )

# finally do some unit-testing also with my newly added features 
add_subdirectory( back_spawn )

add_subdirectory( back_spawn_excit_gen ) 

add_subdirectory( ueg_excit_gens )
