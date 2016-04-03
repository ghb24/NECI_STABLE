
foreach( _config ${${PROJECT_NAME}_CONFIGURATIONS} )

    neci_add_test(
        TARGET test_${_config}_countbits
        SOURCES test_countbits.F90
        LIBS lib${_config} fruit
        LINKER_LANGUAGE Fortran
        META_TARGET test_${_config}
    )

endforeach()
