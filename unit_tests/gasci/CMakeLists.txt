add_library(
    gasci_pchb_test_helper
    gasci_pchb_test_helper.F90
)


set(FRUIT_LIB fruit)
set(_config neci)
target_link_libraries(gasci_pchb_test_helper lib${_config} ${FRUIT_LIB})
get_property(fruit_mods TARGET ${FRUIT_LIB} PROPERTY Fortran_MODULE_DIRECTORY)
target_include_directories(gasci_pchb_test_helper PRIVATE ${fruit_mods})
get_property(lib_mods TARGET lib${_config} PROPERTY Fortran_MODULE_DIRECTORY)
target_include_directories(gasci_pchb_test_helper PRIVATE ${lib_mods})

add_unit_test(neci gasci)
add_unit_test(neci gasci_util)
add_unit_test(neci gasci_on_the_fly_heat_bath)
add_unit_test(neci gasci_disconnected)
add_unit_test(neci gasci_discarding)
add_unit_test(neci gasci_supergroup_index)
add_unit_test(neci gasci_pchb_rhf_hermitian gasci_pchb_test_helper)
add_unit_test(neci gasci_pchb_rhf_nonhermitian gasci_pchb_test_helper)
add_unit_test(neci gasci_pchb_uhf_hermitian gasci_pchb_test_helper)
add_unit_test(neci gasci_pchb_uhf_nonhermitian gasci_pchb_test_helper)
