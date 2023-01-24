neci_add_library(
    TARGET gasci_pchb_test_helper
    SOURCES gasci_pchb_test_helper.F90
    LINKER_LANGUAGE Fortran
    TYPE STATIC
    LIBS lib${_config} fruit
    PRIVATE_INCLUDES
        "${MPI_Fortran_INCLUDE_PATH}"
        "${HDF5_NECI_INCLUDE_PATH}"
    DEFINITIONS
        "${NECI_GLOBAL_DEFINES}"
        "${${PROJECT_NAME}_${_config}_DEFINITIONS}"
    )

target_link_libraries(gasci_pchb_test_helper fruit)

add_unit_test(neci gasci)
add_unit_test(neci gasci_util)
add_unit_test(neci gasci_on_the_fly_heat_bath)
add_unit_test(neci gasci_disconnected)
add_unit_test(neci gasci_discarding)
add_unit_test(neci gasci_supergroup_index)
add_unit_test(neci gasci_pchb)
