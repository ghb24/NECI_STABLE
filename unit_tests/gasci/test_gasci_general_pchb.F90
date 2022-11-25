module test_gasci_general_pchb
    use fruit, only: assert_equals, assert_true, assert_false
    use constants, only: dp, int64, n_int, maxExcit
    use util_mod, only: operator(.div.), operator(.isclose.), near_zero
    use util_mod, only: factrl, swap, cumsum
    use orb_idx_mod, only: calc_spin_raw, sum, SpinOrbIdx_t
    use SystemData, only: nEl
    use excitation_types, only: Excitation_t

    use gasci, only: LocalGASSpec_t
    use gasci_pchb_general, only: GAS_PCHB_ExcGenerator_t, possible_GAS_singles
    use gasci_pc_select_particles, only: PCHB_particle_selections
    use excitation_generators, only: ExcitationGenerator_t

    use sltcnd_mod, only: dyn_sltcnd_excit_old
    use unit_test_helper_excitgen, only: test_excitation_generator, &
        init_excitgen_test, finalize_excitgen_test, generate_random_integrals, &
        FciDumpWriter_t
    use unit_test_helpers, only: run_excit_gen_tester

    use SystemData, only: nEl
    implicit none
    private
    public :: test_pgen_rhf_hermitian

    type :: random_fcidump_writer_t
        logical :: uhf
        logical :: hermitian
    contains
        private
        procedure, public :: random_fcidump_member
    end type random_fcidump_writer_t


contains

    ! @jph more tests here, at least (these can be done through one function):
    ! - [x] rhf hermitian
    ! - [ ] uhf hermitian
    ! - [ ] rhf nonhermitian
    ! - [ ] uhf nonhermitian

    subroutine test_pgen_general(uhf, hermitian)
        logical :: uhf, hermitian
        ! dumpwriter = random_fcidump_writer_t(uhf=uhf, hermitian=hermitian)
        ! @jph

    end subroutine test_pgen_general

    subroutine test_pgen_rhf_hermitian()
        use FciMCData, only: pSingles, pDoubles, pParallel
        type(GAS_PCHB_ExcGenerator_t) :: exc_generator
        type(LocalGASSpec_t) :: GAS_spec
        type(random_fcidump_writer_t) :: dumpwriter
        integer, parameter :: det_I(6) = [1, 2, 3, 7, 8, 10]

        logical :: successful
        integer :: n_interspace_exc
        integer, parameter :: n_iters=10**7
        logical :: is_uhf = .false.

        dumpwriter = random_fcidump_writer_t(uhf=is_uhf, hermitian=.true.)

        pParallel = 0.05_dp
        pSingles = 0.2_dp
        pDoubles = 1.0_dp - pSingles


        do n_interspace_exc = 0, 1
            GAS_spec = LocalGASSpec_t(n_min=[3, 3] - n_interspace_exc, n_max=[3, 3] + n_interspace_exc, &
                                 spat_GAS_orbs=[1, 1, 1, 2, 2, 2])
            call assert_true(GAS_spec%is_valid())
            call assert_true(GAS_spec%contains_conf(det_I))

            call init_excitgen_test(det_I, FciDumpWriter_t(random_fcidump, 'FCIDUMP'))
            call exc_generator%init(GAS_spec, use_lookup=.false., create_lookup=.false., &
                                    used_singles_generator=possible_GAS_singles%PC_UNIFORM, &
                                    PCHB_particle_selection=PCHB_particle_selections%UNIFORM, &
                                    is_uhf=is_uhf)
            call run_excit_gen_tester( &
                exc_generator, 'general implementation, Li2 like system', &
                opt_nI=det_I, &
                opt_n_dets=n_iters, &
                problem_filter=is_problematic,&
                successful=successful)
            call exc_generator%finalize()
            call assert_true(successful)
            call finalize_excitgen_test()
        end do


        do n_interspace_exc = 0, 1
            GAS_spec = LocalGASSpec_t(n_min=[3, 3] - n_interspace_exc, n_max=[3, 3] + n_interspace_exc, &
                                 spat_GAS_orbs=[1, 1, 1, 2, 2, 2])
            call assert_true(GAS_spec%is_valid())
            call assert_true(GAS_spec%contains_conf(det_I))

            call init_excitgen_test(det_I, FciDumpWriter_t(random_fcidump, 'FCIDUMP'))
            call exc_generator%init(GAS_spec, use_lookup=.false., create_lookup=.false., &
                                    used_singles_generator=possible_GAS_singles%PC_UNIFORM, &
                                    PCHB_particle_selection=PCHB_particle_selections%PC_WEIGHTED, &
                                    is_uhf=is_uhf)
            call run_excit_gen_tester( &
                exc_generator, 'general implementation, Li2 like system', &
                opt_nI=det_I, &
                opt_n_dets=n_iters, &
                problem_filter=is_problematic,&
                successful=successful)
            call exc_generator%finalize()
            call assert_true(successful)
            call finalize_excitgen_test()
        end do

        do n_interspace_exc = 0, 1
            GAS_spec = LocalGASSpec_t(n_min=[3, 3] - n_interspace_exc, n_max=[3, 3] + n_interspace_exc, &
                                 spat_GAS_orbs=[1, 1, 1, 2, 2, 2])
            call assert_true(GAS_spec%is_valid())
            call assert_true(GAS_spec%contains_conf(det_I))

            call init_excitgen_test(det_I, FciDumpWriter_t(random_fcidump, 'FCIDUMP'))
            call exc_generator%init(GAS_spec, use_lookup=.false., create_lookup=.false., &
                                    used_singles_generator=possible_GAS_singles%PC_UNIFORM, &
                                    PCHB_particle_selection=PCHB_particle_selections%PC_WEIGHTED_APPROX, &
                                    is_uhf=is_uhf)
            call run_excit_gen_tester( &
                exc_generator, 'general implementation, Li2 like system', &
                opt_nI=det_I, &
                opt_n_dets=n_iters, &
                problem_filter=is_problematic,&
                successful=successful)
            call exc_generator%finalize()
            call assert_true(successful)
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
    end subroutine test_pgen_rhf_hermitian

    subroutine random_fcidump_member(this, iunit, GAS_spec, det_I)
        ! uhf, hermitian
        ! @jph TODO the uhf, hermitian cannot stay here if it is to work with the way the tests work right now
        class(random_fcidump_writer_t), intent(inout) :: this
        type(LocalGASSpec_t), intent(in) :: GAS_spec
        integer, intent(in) :: det_I(:)
        integer, intent(in) :: iunit
        integer :: n_spat_orb, iGAS

        ! @ jph *temporary*
        logical :: uhf = .false., hermitian = .true.

        n_spat_orb = sum([(GAS_spec%GAS_size(iGAS), iGAS = 1, GAS_spec%nGAS())]) .div. 2

        call generate_random_integrals(&
            iunit, n_el=size(det_I), n_spat_orb=n_spat_orb, &
            sparse=0.5_dp, sparseT=0.5_dp, total_ms=sum(calc_spin_raw(det_I)), &
            uhf=this%uhf, hermitian=this%hermitian)
    end subroutine

end module


program test_gasci_program

    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case
    use util_mod, only: stop_all
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_gasci_general_pchb, only: test_pgen_rhf_hermitian


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
        call run_test_case(test_pgen_rhf_hermitian, "test_pgen")
    end subroutine
end program test_gasci_program
