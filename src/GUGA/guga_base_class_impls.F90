#include "macros.h"
submodule (guga_base_class) guga_base_class_impls
    use SystemData, only: t_guga_pchb
    use guga_pchb_class, only: PCHB_t

    better_implicit_none
contains

    module subroutine init_guga_class()
        if (t_guga_pchb) then
            allocate(PCHB_t :: guga_exc_generator)
            select type(guga_exc_generator)
            type is (PCHB_t)
                call guga_exc_generator%init()
            end select
        end if
    end subroutine
end submodule
