module test_pq_property_index_mod

    use fruit, only: assert_equals, assert_true, assert_false, run_test_case
    use constants, only: dp, int64, n_int
    use util_mod, only: operator(.div.), operator(.isclose.), near_zero
    use util_mod, only: factrl, swap, cumsum, custom_findloc

    use gasci, only: LocalGASSpec_t, CumulGASSpec_t, FlexibleGASSpec_t
    use gasci_supergroup_index, only: GAS_SuperGroupIndexer_t
    use composition_utils, only: composition_idx, get_compositions, composition_from_idx
    use pq_property_index, only: PQ_indexer_t

    use sets_mod, only: subset
    implicit none
    private
    public :: test_pq_driver

contains

    subroutine test_pq_driver()
        call run_test_case(test_is_subset, "test_is_subset")
        call run_test_case(test_contained_in_P, "test_contained_in_P")
        call run_test_case(test_allowed_excitations, "test_allowed_excitations")
    end subroutine

    subroutine test_is_subset()
        block
            type(LocalGASSpec_t) :: GAS_spec, P_spec
            type(PQ_indexer_t) :: indexer

            GAS_spec = LocalGASSpec_t(n_min=[1, 1], n_max=[3, 3], spat_GAS_orbs=[1, 1, 2, 2], nEl=4)
            P_spec = LocalGASSpec_t(n_min=[2, 2], n_max=[2, 2], spat_GAS_orbs=[1, 1, 2, 2], nEl=4)

            indexer = PQ_indexer_t(GAS_spec=GAS_spec, P_spec=P_spec)

            !P space is subset of total_space
            call assert_true(subset(indexer%P_space%allowed_composition_indices, &
                                    indexer%total_space%allowed_composition_indices) &
                             )

            !But not the other way around
            call assert_false(subset(indexer%total_space%allowed_composition_indices, &
                                     indexer%P_space%allowed_composition_indices) &
                              )

        end block
    end subroutine

    subroutine test_contained_in_P()
        block

            type(LocalGASSpec_t) :: GAS_spec, P_spec
            type(PQ_indexer_t) :: indexer

            GAS_spec = LocalGASSpec_t(n_min=[1, 1], n_max=[3, 3], spat_GAS_orbs=[1, 1, 2, 2], nEl=4)
            P_spec = LocalGASSpec_t(n_min=[2, 2], n_max=[2, 2], spat_GAS_orbs=[1, 1, 2, 2], nEl=4)

            indexer = PQ_indexer_t(GAS_spec=GAS_spec, P_spec=P_spec)

            call assert_true(size(indexer%contained_in_P) == 3)
            call assert_true(size(indexer%contained_in_P) == size(indexer%total_space%allowed_composition_indices))

            call assert_true(count(indexer%contained_in_P) == 1)
            call assert_true(count(indexer%contained_in_P) == size(indexer%P_space%allowed_composition_indices))
        end block
    end subroutine

    subroutine test_allowed_excitations()
        block
            integer :: i, j, k
            logical :: allowed_excitation(9)
            type(LocalGASSpec_t) :: GAS_spec, P_spec
            type(PQ_indexer_t) :: indexer

            GAS_spec = LocalGASSpec_t(n_min=[1, 1], n_max=[3, 3], spat_GAS_orbs=[1, 1, 2, 2], nEl=4)
            P_spec = LocalGASSpec_t(n_min=[2, 2], n_max=[2, 2], spat_GAS_orbs=[1, 1, 2, 2], nEl=4)

            indexer = PQ_indexer_t(GAS_spec=GAS_spec, P_spec=P_spec)

            k = 1
            do i = 1, 3
                do j = 1, 3
                    allowed_excitation(k) = indexer%contained_in_P(i) .or. indexer%contained_in_P(j)
                    k = k + 1
                end do
            end do

            !There should be 5 allowed excitations (or rather transitions between prop_vecs) 
            !P space to anywhere, and Q->P twice
            call assert_true(count(allowed_excitation) == 5)

        end block
    end subroutine

end module test_pq_property_index_mod

program test_pq_program

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
                     get_failed_count
    use util_mod, only: stop_all
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_pq_property_index_mod, only: test_pq_driver

    implicit none
    integer :: failed_count
    block

        call MPIInit(.false.)

        call init_fruit()

        call test_pq_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_pq_program', 'failed_tests')

        call MPIEnd(.false.)
    end block

end program test_pq_program
