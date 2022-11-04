module test_time_pchb
    use fruit
    use constants, only: dp, int64, n_int
    use util_mod, only: operator(.div.), operator(.isclose.), near_zero, choose
    use util_mod, only: factrl, swap, cumsum
    use orb_idx_mod, only: calc_spin_raw, sum, SpinOrbIdx_t
    use excitation_types, only: Excitation_t

    use gasci, only: LocalGASSpec_t
    use gasci_pchb, only: GAS_PCHB_ExcGenerator_t, possible_GAS_singles
    use gasci_discarding, only: GAS_DiscardingGenerator_t
    use excitation_generators, only: ExcitationGenerator_t

    use sltcnd_mod, only: dyn_sltcnd_excit
    use unit_test_helper_excitgen, only: test_excitation_generator, &
        init_excitgen_test, finalize_excitgen_test, generate_random_integrals, &
        FciDumpWriter_t
    use unit_test_helpers, only: run_excit_gen_tester

    use SystemData, only: nEl
    implicit none
    private
    public :: test_pgen, test_FePor

contains

    subroutine test_pgen()
        use FciMCData, only: pSingles, pDoubles, pParallel
        type(GAS_PCHB_ExcGenerator_t) :: GAS_PCHB
        type(GAS_DiscardingGenerator_t) :: GAS_discarding
        type(LocalGASSpec_t) :: GAS_spec
        integer, parameter :: det_I(6) = [1, 2, 3, 11, 12, 14]

        logical :: successful
        integer, parameter :: n_interspace_exc = 0
        integer, parameter :: n_iters=10**4
        real(dp), parameter :: p_singles(2) = [1.0_dp, 0.0_dp]

        integer :: i_singles

        pParallel = 0.10_dp

        GAS_spec = LocalGASSpec_t(n_min=[3 - n_interspace_exc, size(det_I)], n_max=[3 + n_interspace_exc, size(det_I)], &
                             spat_GAS_orbs=[1, 1, 1, 1, 1, 2, 2, 2, 2, 2])

        call assert_true(GAS_spec%is_valid())
        call assert_true(GAS_spec%contains_conf(det_I))

!         do i_singles = 1, size(p_singles)
!             pSingles = p_singles(i_singles)
!             pDoubles = 1.0_dp - pSingles
!
!             call init_excitgen_test(det_I, FciDumpWriter_t(Li2_FCIDUMP, 'FCIDUMP'))
!             call GAS_discarding%init(GAS_spec)
!             call run_excit_gen_tester( &
!                 GAS_discarding, 'general implementation, Li2 like system', &
!                 opt_nI=det_I, &
!                 opt_n_dets=n_iters, &
!                 problem_filter=is_problematic,&
!                 successful=successful)
!             call GAS_discarding%finalize()
!             call assert_true(successful)
!             call finalize_excitgen_test()
!         end do

        do i_singles = 1, size(p_singles)
            pSingles = p_singles(i_singles)
            pDoubles = 1.0_dp - pSingles

            call init_excitgen_test(det_I, FciDumpWriter_t(Li2_FCIDUMP, 'FCIDUMP'))
            call GAS_PCHB%init(GAS_spec, use_lookup=.false., create_lookup=.false., recoupling=.true., &
                                    used_singles_generator=possible_GAS_singles%BITMASK_UNIFORM)
            call run_excit_gen_tester( &
                GAS_PCHB, 'general implementation, Li2 like system', &
                opt_nI=det_I, &
                opt_n_dets=n_iters, &
                problem_filter=is_problematic,&
                successful=successful)
            call GAS_PCHB%finalize()
            call assert_true(successful)
            call finalize_excitgen_test()
        end do


    contains

        subroutine Li2_FCIDUMP(iunit)
            integer, intent(in) :: iunit
            integer :: n_spat_orb, iGAS
            call copy_file('Li2.FciDmp', iunit)
        end subroutine

        logical function is_problematic(det_I, exc, pgen_diagnostic)
            type(SpinOrbIdx_t), intent(in) :: det_I
            class(Excitation_t), intent(in) :: exc
            real(dp), intent(in) :: pgen_diagnostic
            is_problematic = .false.
!             is_problematic = (abs(1.0_dp - pgen_diagnostic) >= 0.15_dp) &
!                               .and. .not. near_zero(dyn_sltcnd_excit(det_I%idx, exc, .true.))
        end function
    end subroutine test_pgen

    subroutine test_FePor()
        use FciMCData, only: pSingles, pDoubles, pParallel
        type(GAS_PCHB_ExcGenerator_t) :: GAS_PCHB
        type(GAS_DiscardingGenerator_t) :: GAS_discarding
        type(LocalGASSpec_t) :: GAS_spec
        integer, parameter :: det_I(32) = [&
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, &
            33, 34, 35, 37, 38, 39, &
            43, 44, 45, 46, 47, 48, 49, 50]

        logical :: successful
        integer, parameter :: n_iters=10**7
        real(dp), parameter :: p_singles(2) = [1.0_dp, 0.0_dp]

        integer :: i_singles, i

        pParallel = 0.40_dp

        GAS_spec = LocalGASSpec_t(n_min=[18, size(det_I)], n_max=[18, size(det_I)], &
                             spat_GAS_orbs=[(1, i = 1, 16), (2, i = 1, 18)])

        call assert_true(GAS_spec%is_valid())
        call assert_true(GAS_spec%contains_conf(det_I))

!         do i_singles = 1, size(p_singles)
!             pSingles = p_singles(i_singles)
!             pDoubles = 1.0_dp - pSingles
!
!             call init_excitgen_test(det_I, FciDumpWriter_t(FePor_FCIDUMP, 'FCIDUMP'))
!             call GAS_discarding%init(GAS_spec)
!             call run_excit_gen_tester( &
!                 GAS_discarding, 'general implementation, Li2 like system', &
!                 opt_nI=det_I, &
!                 opt_n_dets=n_iters, &
!                 problem_filter=is_problematic,&
!                 successful=successful)
!             call GAS_discarding%finalize()
!             call assert_true(successful)
!             call finalize_excitgen_test()
!         end do

        do i_singles = 1, size(p_singles)
            pSingles = p_singles(i_singles)
            pDoubles = 1.0_dp - pSingles

            call init_excitgen_test(det_I, FciDumpWriter_t(FePor_FCIDUMP, 'FCIDUMP'))
            call GAS_PCHB%init(GAS_spec, use_lookup=.false., create_lookup=.false., recoupling=.true., &
                                    used_singles_generator=possible_GAS_singles%BITMASK_UNIFORM)
            call run_excit_gen_tester( &
                GAS_PCHB, 'general implementation, Li2 like system', &
                opt_nI=det_I, &
                opt_n_dets=n_iters, &
                problem_filter=is_problematic,&
                successful=successful)
            call GAS_PCHB%finalize()
            call assert_true(successful)
            call finalize_excitgen_test()
        end do


    contains

        subroutine FePor_FCIDUMP(iunit)
            integer, intent(in) :: iunit
            integer :: n_spat_orb, iGAS
            call copy_file('tripl.FciDmp', iunit)
        end subroutine

        logical function is_problematic(det_I, exc, pgen_diagnostic)
            type(SpinOrbIdx_t), intent(in) :: det_I
            class(Excitation_t), intent(in) :: exc
            real(dp), intent(in) :: pgen_diagnostic
            is_problematic = .false.
!             is_problematic = (abs(1.0_dp - pgen_diagnostic) >= 0.15_dp) &
!                               .and. .not. near_zero(dyn_sltcnd_excit(det_I%idx, exc, .true.))
        end function
    end subroutine test_FePor

    subroutine copy_file(path, iunit)
        character(*), intent(in) :: path
        integer, intent(in) :: iunit

        integer :: file_id
        character(:), allocatable :: line

        open(file=path, newunit=file_id, status='old', form='formatted', action='read')

        do while (read_line(file_id, line))
            write(iunit, '(A)') line
        end do

        close(file_id)

        contains


        logical function read_line(file_id, line)
            integer, intent(in) :: file_id
            character(:), allocatable, intent(out) :: line
            character(*), parameter :: this_routine = 'read_line'
            integer :: iread
            character(1024) :: tmp_line
            read_line = .false.
            read(file_id, "(A)", iostat=iread) tmp_line
            if (iread > 0) then
                line = ''
                call stop_all(this_routine, 'Error in read_next')
            else if (is_iostat_end(iread)) then
                line = ''
            else
                line = trim(tmp_line)
                read_line = .true.
            end if
        end function

    end subroutine


end module


program test_time_pchb_prog

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case
    use util_mod, only: stop_all
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_time_pchb, only: test_pgen, test_FePor


    implicit none
    integer :: failed_count

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
! We want to make sure it compiles, but don't run in the default test-suite
!         call run_test_case(test_pgen, "test_pgen")
!         call run_test_case(test_FePor, "test_pgen")
    end subroutine
end program
