module time_gasci_general_pchb
    use fruit
    use constants, only: dp, int64, n_int, iout
    use util_mod, only: operator(.div.), operator(.isclose.), near_zero, choose
    use util_mod, only: factrl, intswap, cumsum
    use orb_idx_mod, only: calc_spin_raw, sum, SpinOrbIdx_t
    use excitation_types, only: Excitation_t

    use gasci, only: GASSpec_t
    use gasci_general_pchb, only: GAS_PCHB_ExcGenerator_t
    use gasci_supergroup_index, only: SuperGroupIndexer_t

    use sltcnd_mod, only: dyn_sltcnd_excit
    use unit_test_helper_excitgen, only: test_excitation_generator, &
        init_excitgen_test, finalize_excitgen_test, generate_random_integrals, &
        FciDumpWriter_t
    use unit_test_helpers, only: run_excit_gen_tester

    use SystemData, only: nEl

    use timing_neci, only: timer, get_total_time, set_timer, halt_timer

    implicit none
    integer, parameter :: global_n_iter = 1 * 10**7, global_sg_iter = 1 * 10**7

contains

    subroutine time_FePor()
        use SystemData, only: tGASSpinRecoupling
        use FciMCData, only: pSingles, pDoubles, pParallel
        use LoggingData, only: nPrintTimer
        type(GASSpec_t) :: GAS_spec
        type(GAS_PCHB_ExcGenerator_t) :: exc_generator
        integer, parameter :: &
            det_I(32) = &
                [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, &
                 33, 34, 35, 36, 37, 39, 43, 44, 45, 46, 47, 48, 49, 50]

        logical :: successful
        integer :: sg_sum
        integer, parameter :: n_iter = global_n_iter, n_sg_iter = global_sg_iter
        integer :: n_exc

        type(FciDumpWriter_t) :: fcidump_writers(1)
        real(dp), parameter :: p_singles(2) = [1.00_dp, 1.00_dp]
        logical, parameter :: discarding_singles(2) = [.false., .true.]
        character(20), parameter :: names(2) = [character(20) :: 'PC singles', 'discarding singles']
!         real(dp), parameter :: p_singles(1) = [1.0_dp]
!         logical, parameter :: discarding_singles(1) = [.false.]

        integer :: i_single_exc_gen, i_fc_wr, i_p_singles

        fcidump_writers = [FciDumpWriter_t(FePor_fcidump, 'FCIDUMP')]
        pParallel = 0.40_dp

        call assert_true(tGASSpinRecoupling)

        do i_p_singles = 1, size(p_singles)
            pSingles = p_singles(i_p_singles)
            pDoubles = 1.0_dp - pSingles
            do i_fc_wr = 1, size(fcidump_writers)
                do n_exc = 0, 0


                    GAS_spec = get_GAS_spec(n_exc)
                    call assert_true(GAS_spec%is_valid())
                    call assert_true(GAS_spec%contains_det(det_I))

                    call init_excitgen_test(det_I, fcidump_writers(i_fc_wr))

                    call exc_generator%init(GAS_spec, use_lookup=.false., create_lookup=.false., recoupling=.true., discarding_singles=discarding_singles(i_p_singles))

                    call run_excit_gen_tester( &
                        exc_generator, names(i_p_singles), &
                        opt_nI=det_I, opt_n_iters=n_iter, &
                        problem_filter=is_problematic,&
                        successful=successful)

                    call exc_generator%finalize()
                    call assert_true(successful)


                    block
                        type(timer) :: SG_timer
                        type(SuperGroupIndexer_t) :: indexer
                        integer :: i
                        indexer = SuperGroupIndexer_t(GAS_spec)

                        sg_sum = 0
                        call set_timer(SG_timer)
                        do i = 1, n_sg_iter
                            sg_sum = indexer%idx_nI(det_I) + sg_sum
                        end do
                        call halt_timer(SG_timer)

                        write(iout, *) sg_sum
                        write(iout, *)
                        write(iout, *) 'n exc', n_exc
                        write(iout, *) 'Total time', get_total_time(SG_timer)
                        write(iout, *) 'time per supergroup idx in mikro seconds', get_total_time(SG_timer) / real(n_sg_iter, dp) * 1e6_dp
                        write(iout, *)
                    end block

!                     block
!                         use gasci, only: global_GAS_spec => GAS_specification
!                         use gasci_discarding, only: GAS_DiscardingGenerator_t
!                         type(GAS_DiscardingGenerator_t) :: exc_generator
!
!                         call exc_generator%init(GAS_spec)
!                         call run_excit_gen_tester( &
!                             exc_generator, 'discarding GASCI implementation, random fcidump', &
!                             opt_nI=det_I, opt_n_iters=n_iter, &
!                             problem_filter=is_problematic,&
!                             successful=successful)
!                         call assert_true(successful)
!                         call exc_generator%finalize()
!                     end block
                    call finalize_excitgen_test()
                end do
            end do
        end do
    contains

        subroutine FePor_fcidump(iunit)
            integer, intent(in) :: iunit

            integer :: in_file, stat
            character(len=500) :: line

            open(file='./tripl.FciDmp', newunit=in_file, action='read', status='old')
                read(in_file, '(A)', iostat=stat) line
                do while (stat == 0)
                    write(iunit, '(A)') line
                    read(in_file, '(A)', iostat=stat) line
                end do
            close(in_file)
        end subroutine


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
!             is_problematic = abs(1._dp - pgen_diagnostic) >= 0.15_dp &
!                               .and. .not. near_zero(dyn_sltcnd_excit(det_I%idx, exc, .true.))
            is_problematic = .false.
        end function

        pure function get_GAS_spec(n_interspace_excitations) result(GAS_spec)
            integer, intent(in) :: n_interspace_excitations
            type(GASSpec_t) :: GAS_spec

            integer :: i, j

            associate(n => n_interspace_excitations)
                GAS_spec = GASSpec_t(&
                    n_min=[18 - n, 32], &
                    n_max=[18 + n, 32], &
                    spat_GAS_orbs = [(1, i = 1, 16), (2, i = 1, 18)])
            end associate
        end function
    end subroutine

end module


program time_gasci_program

    use mpi
    use fruit
    use Parallel_neci, only: MPIInit, MPIEnd
    use time_gasci_general_pchb, only: time_FePor


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
    write(*, *) 'I switched off the tests for the automatic unit test pipeline'
    write(*, *) '   because it takes too long.'
    write(*, *) 'You have to uncomment the following line in the code.'
!         call run_test_case(time_FePor, "test_pgen")
    end subroutine
end program
