foreach( _config mneci)
    neci_add_test(
        TARGET
            test_${_config}_loop
        SOURCES
            test_loop.F90 test_loop_testcases.F90 test_loop_helpers.F90
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
endforeach()
