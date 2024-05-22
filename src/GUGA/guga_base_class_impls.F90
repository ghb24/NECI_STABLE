#include "macros.h"
submodule (guga_base_class) guga_base_class_impls
    use generalized_supergroup, only: indexer, AlsoGUGA_SuperGroupIndexer_t
    use guga_pchb_class, only: PCHB_t
    use guga_gas_pchb_main, only: GAS_PCHB_ExcGenerator_t, options

    better_implicit_none
contains

    module subroutine init_guga_class()
        routine_name("guga_pchb_class::init_guga_class")
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
                    select type(indexer)
                    class is (AlsoGUGA_SuperGroupIndexer_t)
                        call GUGA_exc_generator%init(indexer, options)
                    class default
                        call stop_all(this_routine, "GUGA PCHB requires an indexer compatible to it.")
                    end select
                end block
                end select
            end if
        end if
    end subroutine
end submodule
