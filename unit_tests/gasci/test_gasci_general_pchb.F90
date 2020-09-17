module test_gasci_general_pchb
    use fruit
    use constants, only: dp, int64, n_int
    use util_mod, only: operator(.div.), operator(.isclose.), near_zero, choose
    use util_mod, only: factrl, intswap, cumsum
    use orb_idx_mod, only: calc_spin_raw, sum, SpinOrbIdx_t
    use excitation_types, only: Excitation_t

    use gasci, only: GASSpec_t
    use gasci_general_pchb
    use gasci_general, only: gen_all_excits

    use sltcnd_mod, only: dyn_sltcnd_excit
    use unit_test_helper_excitgen, only: test_excitation_generator, &
        init_excitgen_test, finalize_excitgen_test, generate_random_integrals, &
        FciDumpWriter_t
    use unit_test_helpers, only: run_excit_gen_tester
    implicit none
    private
    public :: test_pgen, test_partitioning, test_supergroup_offsets



contains


    subroutine test_pgen()
        use gasci, only: global_GAS_spec => GAS_specification
        use SystemData, only: tGASSpinRecoupling
        use FciMCData, only: pSingles, pDoubles, pParallel
        type(GASSpec_t) :: GAS_spec
        integer, parameter :: det_I(6) = [1, 2, 3, 7, 8, 10]

        logical :: successful
        integer, parameter :: n_iters=5 * 10**6

        pParallel = 0.05_dp
        pSingles = 0.3_dp
        pDoubles = 1.0_dp - pSingles

        call assert_true(tGASSpinRecoupling)

        GAS_spec = GASSpec_t(n_min=[2, size(det_I)], n_max=[4, size(det_I)], &
                             spat_GAS_orbs=[1, 1, 1, 2, 2, 2])
        call assert_true(GAS_spec%is_valid())
        call assert_true(GAS_spec%contains(det_I))
        global_GAS_spec = GAS_spec

        call init_excitgen_test(size(det_I), FciDumpWriter_t(random_fcidump, 'FCIDUMP'))
        call run_excit_gen_tester( &
            gen_general_GASCI_pchb, 'discarding GASCI implementation, random fcidump', &
            opt_nI=det_I, opt_n_iters=n_iters, &
            gen_all_excits=gen_all_excits, &
            problem_filter=is_problematic,&
            successful=successful)
        call assert_true(successful)
        call finalize_excitgen_test()

    contains

        subroutine random_fcidump(iunit)
            integer, intent(in) :: iunit
            integer :: n_spat_orb, iGAS

            n_spat_orb = sum([(GAS_spec%GAS_size(iGAS), iGAS = 1, GAS_spec%nGAS())]) .div. 2

            call generate_random_integrals(&
                iunit, n_el=size(det_I), n_spat_orb=n_spat_orb, &
                sparse=1.0_dp, sparseT=1.0_dp, total_ms=sum(calc_spin_raw(det_I)))
        end subroutine

        logical function is_problematic(det_I, exc, pgen_diagnostic)
            type(SpinOrbIdx_t), intent(in) :: det_I
            class(Excitation_t), intent(in) :: exc
            real(dp), intent(in) :: pgen_diagnostic
            is_problematic = (abs(1.0_dp - pgen_diagnostic) >= 0.15_dp) &
                              .and. .not. near_zero(dyn_sltcnd_excit(det_I%idx, exc, .true.))
        end function

    end subroutine test_pgen

    subroutine test_partitioning()
        integer(int64), allocatable :: partitions(:, :), cn_min(:), cn_max(:)

        integer :: i, n, k
        logical :: correct

        do n = 1, 10
            do k = 1, 10
                partitions = get_partitions(k, n)
                correct = .true.
                do i = 1, size(partitions, 2)
                    if (i /= get_partition_index(partitions(:, i))) then
                        correct = .false.
                        exit
                    end if
                end do
                call assert_true(correct)
            end do
        end do

!         block
!             logical :: problem
!
!             write(*, *)
!             partitions = get_partitions(3, 3)
!             cn_min = [0, 1, 3]
!             cn_max = [2, 2, 3]
!             do i = 1, size(partitions, 2)
!                 problem = .not. (all(cn_min <= cumsum(partitions(:, i))) .and. all(cumsum(partitions(:, i)) <= cn_max))
!                 write(*, *) '*', partitions(:, i), merge('!!!', '   ', problem)
!                 write(*, *) '>', cumsum(partitions(:, i)), merge('!!!', '   ', problem)
!             end do
!         end block

        block
            write(*, *)
            partitions = get_partitions(2, 2)

            do i = 1, size(partitions, 2)
                write(*, *) partitions(:, i)
            end do
        end block


        block
            logical :: correct

            write(*, *)
            partitions = get_partitions(3, 3)
            cn_min = [0, 1, 3]
            cn_max = [2, 2, 3]

            write(*, *) new_get_n_partitions(3, 3, int(cn_min), int(cn_max))

            do i = 1, size(partitions, 2)
                correct = all(cn_min <= cumsum(partitions(:, i))) .and. all(cumsum(partitions(:, i)) <= cn_max)

                if (correct) then
                    write(*, *) partitions(:, i)
!                     write(*, *) new_get_partition_index(int(partitions(:, i)), int(cn_min), int(cn_max)), partitions(:, i)
!                     write(*, *) '>', cumsum(partitions(:, i))
                end if
            end do
        end block


!         block
!             logical :: correct
!
!             write(*, *)
!             partitions = get_partitions(5, 30)
!             cn_min = [5, 11, 17, 23, 30]
!             cn_max = [7, 13, 19, 25, 30]
!
!             do i = 1, size(partitions, 2)
!                 correct = all(cn_min <= cumsum(partitions(:, i))) .and. all(cumsum(partitions(:, i)) <= cn_max)
!
!                 if (correct) then
!                     write(*, *) '*', partitions(:, i)
!                     write(*, *) '>', cumsum(partitions(:, i))
!                 end if
!             end do
!         end block
    end subroutine


    subroutine test_supergroup_offsets()
        type(GASSpec_t) :: GAS_spec
        integer, parameter :: nEl = 6

        GAS_spec = GASSpec_t(n_min=[2, nEl], n_max=[4, nEl], &
                             spat_GAS_orbs=[1, 1, 1, 2, 2, 2])
        call assert_true(GAS_spec%is_valid())



    end subroutine
end module test_gasci_general_pchb

program test_gasci_program

    use mpi
    use fruit
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_gasci_general_pchb, only: test_pgen, test_partitioning, test_supergroup_offsets


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
!         call run_test_case(test_pgen, "test_pgen")
        call run_test_case(test_partitioning, "test_partition_vectors")
        call run_test_case(test_supergroup_offsets, "test_supergroup_offsets")
    end subroutine
end program test_gasci_program
