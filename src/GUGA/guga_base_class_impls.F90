#include "macros.h"
submodule (guga_base_class) guga_base_class_impls
    use guga_pchb_class, only: PCHB_t, options

    better_implicit_none
contains

    module subroutine init_guga_class()

        if (GUGA_selected_exc_generator == GUGA_exc_generator_vals%FCI_PCHB) then
            allocate(PCHB_t :: guga_exc_generator)
            select type(guga_exc_generator)
            type is (PCHB_t)
                call guga_exc_generator%init(options)
            end select
        end if
    end subroutine
end submodule
