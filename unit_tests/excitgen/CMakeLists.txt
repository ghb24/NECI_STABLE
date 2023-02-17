add_library(
    pchb_excitgen_test_helper
    pchb_excitgen_test_helper.F90
)

set(FRUIT_LIB fruit)
set(_config neci)
target_link_libraries(pchb_excitgen_test_helper lib${_config} ${FRUIT_LIB})
get_property(fruit_mods TARGET ${FRUIT_LIB} PROPERTY Fortran_MODULE_DIRECTORY)
target_include_directories(pchb_excitgen_test_helper PRIVATE ${fruit_mods})
get_property(lib_mods TARGET lib${_config} PROPERTY Fortran_MODULE_DIRECTORY)
target_include_directories(pchb_excitgen_test_helper PRIVATE ${lib_mods})


add_unit_test(neci pcpp_excitgen)
add_unit_test(neci pchb_excitgen_rhf_hermitian pchb_excitgen_test_helper)
add_unit_test(neci pchb_excitgen_rhf_nonhermitian pchb_excitgen_test_helper)
add_unit_test(neci pchb_excitgen_uhf_hermitian pchb_excitgen_test_helper)
add_unit_test(neci pchb_excitgen_uhf_nonhermitian pchb_excitgen_test_helper)
