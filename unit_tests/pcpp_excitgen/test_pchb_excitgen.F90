program test_pchb_excitgen
    use constants, only: dp, n_int, maxExcit
    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case, assert_true, &
        assert_equals
    use util_mod, only: stop_all
    use unit_test_helper_excitgen, only: set_ref, free_ref, calc_pgen, &
        generate_random_integrals, FciDumpWriter_t, test_excitation_generator, &
        init_excitgen_test
    use Parallel_neci, only: MPIInit, MPIEnd
    use pchb_excitgen, only: PCHB_FCI_excit_generator_t
    use orb_idx_mod, only: beta
    use procedure_pointers, only: generate_excitation

    implicit none

    type(PCHB_FCI_excit_generator_t) :: PCHB_FCI

    call MPIInit(.false.)
    call init_fruit()
    call pchb_test_driver()
    call fruit_summary()
    call fruit_finalize()
    block
        integer :: failed_count
        call get_failed_count(failed_count)
        if (failed_count /= 0) call stop_all('test_pchb_excitgen', 'failed_tests')
    end block
    call MPIEnd(.false.)

contains

    subroutine pchb_test_driver()
        implicit none
        real(dp) :: pTot, pNull
        integer :: numEx, nFound, i
        ! There can be some excitations with really low matrix elements -> we need a lot
        ! of samples to hit all
        integer, parameter :: nSamples = 1000000

        ! set the excitation generator to pchb
        generate_excitation => gen_rand_excit_pchb
        calc_pgen => calc_pgen_pchb

        ! prepare an excitation generator test
        call init_excitgen_test(ref_det=[(i, i=1, 5)], fcidump_writer=FciDumpWriter_t(random_fcidump, 'FCIDUMP'))

        ! prepare the pchb excitgen: set the weights/map-table
        call set_ref()
        call PCHB_FCI%init()

        ! test the excitation generator
        call test_excitation_generator(nSamples, pTot, pNull, numEx, nFound, .true.)
        ! make sure all excits have been found
        call assert_equals(numEx, nFound)
        ! and the total prob is 1.0
        call assert_true(abs(1.0 - pTot) < 0.05)

        ! free memory
        call free_ref()
        call PCHB_FCI%finalize()
    end subroutine pchb_test_driver

    function calc_pgen_pchb(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
        use SymExcitDataMod, only: scratchSize
        use bit_rep_data, only: NIfTot
        use SystemData, only: nel
        implicit none
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, 2), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)

        real(dp) :: pgen
        integer :: help_ex(2, maxExcit)
        help_ex(:, :2) = ex

        pgen = PCHB_FCI%get_pgen(nI, ilutI, help_ex, ic, ClassCount2, ClassCountUnocc2)
    end function

    subroutine gen_rand_excit_pchb(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                   ex, tParity, pGen, hel, store, part_type)

        use SystemData, only: nel
        use bit_rep_data, only: NIfTot
        use FciMCData, only: excit_gen_store_type
        implicit none

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        call PCHB_FCI%gen_exc(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                              ex, tParity, pGen, hel, store, part_type)
    end subroutine

    subroutine random_fcidump(iunit)
        integer, intent(in) :: iunit
        call generate_random_integrals( &
            iunit, n_el=5, n_spat_orb=12, sparse=0.9_dp, sparseT=0.1_dp, total_ms=beta)
    end subroutine
end program test_pchb_excitgen
