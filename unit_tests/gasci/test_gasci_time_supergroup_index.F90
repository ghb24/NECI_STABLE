module test_gasci_supergroup_index_mod
    use fruit
    use constants, only: dp, int64, n_int, iout
    use util_mod, only: operator(.div.), operator(.isclose.), near_zero, choose
    use util_mod, only: factrl, intswap, cumsum

    use gasci, only: GASSpec_t
    use gasci_general_pchb
    use gasci_supergroup_index

    implicit none
    private
    public :: test_time_supergroup_indexer_class

contains

    subroutine test_time_supergroup_indexer_class()
        use timing_neci, only: timer, get_total_time, set_timer, halt_timer

        integer, parameter :: n_iter = 1000000
        integer, allocatable :: determinant(:)
        integer, allocatable :: sg_indices(:), comp_indices(:)
        integer :: n_exc

        determinant =  [1, 2, 3, 4, 5, 6, &
                        13, 14, 15, 16, 17, 18, &
                        25, 26, 27, 28, 29, 30, &
                        37, 38, 39, 40, 41, 42, &
                        49, 50, 51, 52, 53, 54]
        allocate(sg_indices(n_iter), comp_indices(n_iter))


        do n_exc = 0, 4
            block
                type(SuperGroupIndexer_t) :: indexer
                type(GASSpec_t) :: GAS_spec

                GAS_spec = get_GAS_spec(n_exc)
                indexer = SuperGroupIndexer_t(GAS_spec)

                block
                    type(timer) :: SG_timer
                    integer :: i


                    call set_timer(SG_timer)
                    do i = 1, n_iter
                        sg_indices(i) = indexer%idx_nI(determinant)
                    end do
                    call halt_timer(SG_timer)

                    write(iout, *) sum(sg_indices)
                    write(iout, *)
                    write(iout, *) 'n exc', n_exc
                    write(iout, *) 'Total time', get_total_time(SG_timer)
                    write(iout, *) 'time per supergroup idx in mikro seconds', get_total_time(SG_timer) / real(n_iter, dp) * 1e6_dp
                    write(iout, *)
                end block

                block
                    type(timer) :: SG_timer
                    integer :: i
                    integer, allocatable :: supergroup(:)

                    supergroup = GAS_spec%count_per_GAS(determinant)

                    call set_timer(SG_timer)
                    do i = 1, n_iter
                        sg_indices(i) = composition_idx(supergroup)
                    end do
                    call halt_timer(SG_timer)

                    write(iout, *) sum(sg_indices)
                    write(iout, *)
                    write(iout, *) 'n exc', n_exc
                    write(iout, *) 'Total time', get_total_time(SG_timer)
                    write(iout, *) 'time per composition idx in mikro seconds', get_total_time(SG_timer) / real(n_iter, dp) * 1e6_dp
                    write(iout, *)
                end block
            end block
        end do

        contains

        pure function get_GAS_spec(n_interspace_excitations) result(GAS_spec)
            integer, intent(in) :: n_interspace_excitations
            type(GASSpec_t) :: GAS_spec

            integer :: i, j

            associate(n => n_interspace_excitations)
                GAS_spec = GASSpec_t(&
                    n_min=[6 - n, 12 - n, 18 - n, 24 - n, 30], &
                    n_max=[6 + n, 12 + n, 18 + n, 24 + n, 30], &
                    spat_GAS_orbs = [([(j, i = 1, 6)], j = 1, 5)])
            end associate
        end function

    end subroutine

end module


program test_gasci_program

    use mpi
    use fruit
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_gasci_supergroup_index_mod, only: test_time_supergroup_indexer_class


    implicit none
    integer :: failed_count, err

    integer :: n
    block

        call MPIInit(.false.)

        call init_fruit()

        call test_gasci_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_gasci_program', 'failed_tests')

        call MPIEnd(.false.)
    end block

contains

    subroutine test_gasci_driver()
        call run_test_case(test_time_supergroup_indexer_class, "test_time_supergroup_indexer_class")
    end subroutine
end program test_gasci_program
