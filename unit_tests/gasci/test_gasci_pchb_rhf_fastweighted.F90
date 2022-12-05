module test_gasci_pchb_rhf_fastweighted
    use fruit, only: assert_equals, assert_true, assert_false
    use constants, only: dp, int64, n_int, maxExcit
    use util_mod, only: operator(.div.), operator(.isclose.), near_zero
    use util_mod, only: factrl, swap, cumsum
    use orb_idx_mod, only: calc_spin_raw, sum, SpinOrbIdx_t
    use SystemData, only: nEl
    use excitation_types, only: Excitation_t

    use gasci, only: LocalGASSpec_t
    use gasci_pchb_main, only: GAS_PCHB_ExcGenerator_t, GAS_PCHB_options_t, &
        GAS_PCHB_options_vals
    use gasci_pchb_doubles_main, only: PCHB_DoublesOptions_t
    use gasci_singles_main, only: GAS_PCHB_SinglesOptions_t

    use excitation_generators, only: ExcitationGenerator_t

    use sltcnd_mod, only: dyn_sltcnd_excit_old
    use unit_test_helper_excitgen, only: test_excitation_generator, &
        init_excitgen_test, finalize_excitgen_test, generate_random_integrals, &
        FciDumpWriter_t
    use unit_test_helpers, only: run_excit_gen_tester

    use SystemData, only: nEl
    implicit none
    private
    public :: test_pgen_RHF_hermitian
            ! test_pgen_UHF_hermitian, &
            !   test_pgen_RHF_nonhermitian, test_pgen_UHF_nonhermitian

    type :: random_fcidump_writer_t
        logical :: UHF
        logical :: hermitian
    contains
        private
        procedure, public :: random_fcidump_member
    end type random_fcidump_writer_t


contains

    subroutine test_pgen_RHF_hermitian()
        call test_pgen_general(.false., .true.)
    end subroutine test_pgen_RHF_hermitian

    ! subroutine test_pgen_UHF_hermitian()
    !     call test_pgen_general(.true., .true.)
    ! end subroutine test_pgen_UHF_hermitian
    !
    ! subroutine test_pgen_RHF_nonhermitian()
    !     call test_pgen_general(.false., .false.)
    ! end subroutine test_pgen_RHF_nonhermitian
    !
    ! subroutine test_pgen_UHF_nonhermitian()
    !     call test_pgen_general(.true., .false.)
    ! end subroutine test_pgen_UHF_nonhermitian

    subroutine test_pgen_general(UHF, hermitian)
        use FciMCData, only: pSingles, pDoubles, pParallel
        logical, intent(in) :: UHF, hermitian
        type(GAS_PCHB_ExcGenerator_t) :: exc_generator
        type(LocalGASSpec_t) :: GAS_spec
        type(random_fcidump_writer_t) :: dumpwriter
        integer, parameter :: det_I(6) = [1, 2, 3, 7, 8, 10]

        logical :: successful
        integer :: n_interspace_exc
        integer, parameter :: n_iters=10**7

        character(len=128) :: message
        write(message, *) "Failed for UHF=", UHF, ", hermitian=", hermitian

        dumpwriter = random_fcidump_writer_t(UHF=UHF, hermitian=hermitian)

        pParallel = 0.05_dp
        pSingles = 0.2_dp
        pDoubles = 1.0_dp - pSingles


        do n_interspace_exc = 0, 1
            GAS_spec = LocalGASSpec_t(n_min=[3, 3] - n_interspace_exc, n_max=[3, 3] + n_interspace_exc, &
                                 spat_GAS_orbs=[1, 1, 1, 2, 2, 2])
            call assert_true(GAS_spec%is_valid())
            call assert_true(GAS_spec%contains_conf(det_I))

            call init_excitgen_test(det_I, FciDumpWriter_t(random_fcidump, 'FCIDUMP'))
            call exc_generator%init(&
                GAS_spec, &
                options=GAS_PCHB_options_t(&
                    GAS_PCHB_SinglesOptions_t(&
                        GAS_PCHB_options_vals%singles%algorithm%BITMASK_UNIFORM &
                    ), &
                    PCHB_DoublesOptions_t( &
                        GAS_PCHB_options_vals%doubles%particle_selection%UNIFORM, &
                        GAS_PCHB_options_vals%doubles%hole_selection%RHF_FAST_WEIGHTED &
                    ), &
                    UHF=UHF, &
                    use_lookup=.false. &
                ) &
            )



            call run_excit_gen_tester( &
                exc_generator, 'general implementation, Li2 like system', &
                opt_nI=det_I, &
                opt_n_dets=n_iters, &
                problem_filter=is_problematic,&
                successful=successful)
            call exc_generator%finalize()
            call assert_true(successful, trim(message))
            call finalize_excitgen_test()
        end do


        do n_interspace_exc = 0, 1
            GAS_spec = LocalGASSpec_t(n_min=[3, 3] - n_interspace_exc, n_max=[3, 3] + n_interspace_exc, &
                                 spat_GAS_orbs=[1, 1, 1, 2, 2, 2])
            call assert_true(GAS_spec%is_valid())
            call assert_true(GAS_spec%contains_conf(det_I))

            call init_excitgen_test(det_I, FciDumpWriter_t(random_fcidump, 'FCIDUMP'))
            call exc_generator%init(&
                GAS_spec, &
                options=GAS_PCHB_options_t(&
                    GAS_PCHB_SinglesOptions_t(&
                        GAS_PCHB_options_vals%singles%algorithm%BITMASK_UNIFORM &
                    ), &
                    PCHB_DoublesOptions_t( &
                        GAS_PCHB_options_vals%doubles%particle_selection%PC_WEIGHTED, &
                        GAS_PCHB_options_vals%doubles%hole_selection%RHF_FAST_WEIGHTED &
                    ), &
                    UHF=UHF, &
                    use_lookup=.false. &
                ) &
            )
            call run_excit_gen_tester( &
                exc_generator, 'general implementation, Li2 like system', &
                opt_nI=det_I, &
                opt_n_dets=n_iters, &
                problem_filter=is_problematic,&
                successful=successful)
            call exc_generator%finalize()
            call assert_true(successful, trim(message))
            call finalize_excitgen_test()
        end do

        do n_interspace_exc = 0, 1
            GAS_spec = LocalGASSpec_t(n_min=[3, 3] - n_interspace_exc, n_max=[3, 3] + n_interspace_exc, &
                                 spat_GAS_orbs=[1, 1, 1, 2, 2, 2])
            call assert_true(GAS_spec%is_valid())
            call assert_true(GAS_spec%contains_conf(det_I))

            call init_excitgen_test(det_I, FciDumpWriter_t(random_fcidump, 'FCIDUMP'))
            call exc_generator%init(&
                GAS_spec, &
                options=GAS_PCHB_options_t(&
                    GAS_PCHB_SinglesOptions_t(&
                        GAS_PCHB_options_vals%singles%algorithm%BITMASK_UNIFORM &
                    ), &
                    PCHB_DoublesOptions_t( &
                        GAS_PCHB_options_vals%doubles%particle_selection%PC_WEIGHTED_APPROX, &
                        GAS_PCHB_options_vals%doubles%hole_selection%RHF_FAST_WEIGHTED &
                    ), &
                    UHF=UHF, &
                    use_lookup=.false. &
                ) &
            )
            call run_excit_gen_tester( &
                exc_generator, 'general implementation, Li2 like system', &
                opt_nI=det_I, &
                opt_n_dets=n_iters, &
                problem_filter=is_problematic,&
                successful=successful)
            call exc_generator%finalize()
            call assert_true(successful, trim(message))
            call finalize_excitgen_test()
        end do

    contains

        subroutine random_fcidump(iunit)
            integer, intent(in) :: iunit
            ! dumpwriter comes from the host subroutine
            call dumpwriter%random_fcidump_member(iunit, GAS_spec, det_I)
        end subroutine random_fcidump

        logical function is_problematic(nI, exc, ic, pgen_diagnostic)
            integer, intent(in) :: nI(nEl), exc(2, maxExcit), ic
            real(dp), intent(in) :: pgen_diagnostic
            is_problematic = &
                (abs(1._dp - pgen_diagnostic) >= 0.15_dp) &
                .and. .not. near_zero(dyn_sltcnd_excit_old(nI, ic, exc, .true.))
        end function
    end subroutine test_pgen_general

    subroutine random_fcidump_member(this, iunit, GAS_spec, det_I)
        ! UHF, hermitian
        class(random_fcidump_writer_t), intent(inout) :: this
        type(LocalGASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: det_I(:)
        integer, intent(in) :: iunit
        integer :: n_spat_orb, iGAS

        n_spat_orb = sum([(GAS_spec%GAS_size(iGAS), iGAS = 1, GAS_spec%nGAS())]) .div. 2

        call generate_random_integrals(&
            iunit, n_el=size(det_I), n_spat_orb=n_spat_orb, &
            sparse=0.5_dp, sparseT=0.5_dp, total_ms=sum(calc_spin_raw(det_I)), &
            UHF=this%UHF, hermitian=this%hermitian)
    end subroutine

end module


program test_gasci_program

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case
    use util_mod, only: stop_all
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_gasci_pchb_rhf_fastweighted, only: test_pgen_RHF_hermitian


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
        call run_test_case(test_pgen_RHF_hermitian, "test_pgen_RHF_hermitian")
        ! TODO(@jph): Good luck with that ;-)
        ! call run_test_case(test_pgen_UHF_hermitian, "test_pgen_UHF_hermitian")
        ! call run_test_case(test_pgen_RHF_nonhermitian, "test_pgen_RHF_nonhermitian")
        ! call run_test_case(test_pgen_UHF_nonhermitian, "test_pgen_UHF_nonhermitian")
    end subroutine test_gasci_driver
end program test_gasci_program
