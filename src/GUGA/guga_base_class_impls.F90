#include "macros.h"
submodule (guga_base_class) guga_base_class_impls
    use gasci, only: GAS_specification
    use guga_pchb_class, only: PCHB_t
    use guga_gas_pchb_main, only: GAS_PCHB_ExcGenerator_t

    better_implicit_none
contains

    module subroutine init_guga_class()
        if (allocated(GUGA_selected_exc_generator)) then
            if (GUGA_selected_exc_generator == GUGA_exc_generator_vals%FCI_PCHB) then
                allocate(PCHB_t :: GUGA_exc_generator)
                select type(GUGA_exc_generator)
                type is (PCHB_t)
                block
                    use guga_pchb_class, only: options
                    call GUGA_exc_generator%init(options)
                end block
                end select
            else if (GUGA_selected_exc_generator == GUGA_exc_generator_vals%GAS_PCHB) then
                allocate(GAS_PCHB_ExcGenerator_t :: GUGA_exc_generator)
                select type(GUGA_exc_generator)
                type is (GAS_PCHB_ExcGenerator_t)
                block
                    use guga_gas_pchb_main, only: options
                    call GUGA_exc_generator%init(GAS_specification, options)
                end block
                end select
            end if
        end if
    end subroutine
end submodule
