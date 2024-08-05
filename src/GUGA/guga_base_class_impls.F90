#include "macros.h"
submodule (guga_base_class) guga_base_class_impls
    use property_vector_index, only: indexer, AlsoGUGA_PropertyIndexer_t
    use SystemData, only: nBasis
    use guga_pchb_class, only: PCHB_t
    use guga_prop_vec_pchb_main, only: PropVec_PCHB_ExcGenerator_t, PropVec_PCHB_user_input, decide_on_PCHB_options

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
                    use guga_pchb_class, only: decide_on_PCHB_options, FCI_PCHB_user_input
                    call GUGA_exc_generator%init(decide_on_PCHB_options(FCI_PCHB_user_input, nBasis, nEl))
                end block
                end select
            else if (GUGA_selected_exc_generator == GUGA_exc_generator_vals%PropVec_PCHB) then
                allocate(PropVec_PCHB_ExcGenerator_t :: GUGA_exc_generator)
                select type(GUGA_exc_generator)
                type is (PropVec_PCHB_ExcGenerator_t)
                block
                    select type(indexer)
                    class is (AlsoGUGA_PropertyIndexer_t)
                        call GUGA_exc_generator%init(indexer, decide_on_PCHB_options(PropVec_PCHB_user_input, nBasis, nEl))
                    class default
                        call stop_all(this_routine, "GUGA PCHB requires an indexer compatible to it.")
                    end select
                end block
                end select
            end if
        end if
    end subroutine
end submodule
