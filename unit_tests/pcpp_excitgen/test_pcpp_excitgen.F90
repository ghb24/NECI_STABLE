#include "macros.h"
program test_pcpp_excitgen
    use constants
    use Parallel_neci, only: MPIInit, MPIEnd
    use fruit
    use pcpp_excitgen
    use unit_test_helper_excitgen
    use orb_idx_mod, only: beta
    use procedure_pointers, only: generate_excitation
    use FciMCData, only: projEDet
    implicit none

    call MPIInit(.false.)
    call init_fruit()
    call pcpp_test_driver()
    call fruit_summary()
    call fruit_finalize()
    call MPIEnd(.false.)

contains

    !------------------------------------------------------------------------------------------!

    subroutine pcpp_test_driver()
        implicit none
        real(dp) :: pTot, pNull
        integer :: numEx, nFound, i
        integer, parameter :: nSamples = 200000

        ! set the excitation we want to test
        generate_excitation => gen_rand_excit_pcpp
        ! prepare everything for testing the excitgen
        call init_excitgen_test(ref_det=[(i, i = 1, 5)], fcidump_writer=FciDumpWriter_t(random_fcidump, 'FCIDUMP'))

        ! prepare the pcpp excitation generator: get the precomputed weights

        call set_ref()
        call init_pcpp_excitgen()

        ! test the mapping of electrons
        call test_elec_mapping()
        ! run the excitgen test: do nSamples excitations, and compare them with all possible excits
        calc_pgen => calc_pgen_local
        call test_excitation_generator(nSamples,pTot,pNull,numEx,nFound,.true.)
        ! make sure all excitations are found
        call assert_equals(numEx,nFound)
        ! make sure the probability is normalized - we don't require pNull to match exactly (that would mean generating all possible null excitations)
        ! so the threshold is rather loose
        ! for in-depth insight, study the output and consider the ratio of pGen and the number of occurences per excitation (should fluctuate around 1)
        call assert_true(abs(1.0-pTot)<0.01)

        call free_ref()
        call finalize_excitgen_test()
    end subroutine pcpp_test_driver

    subroutine test_elec_mapping()
        ! test that the mapping of electrons does what its supposed to do
        implicit none
        ! map to a determinant that has 2 electrons excited
        integer :: nI(nel)
        integer :: elec_map(nel)
        integer(n_int) :: ilut(0:NIfTot)
        ! create the excited determinant
        nI = (/1,4,5,8,9/)
        ! create the map between this det and the ref
        call encodeBitDet(nI,ilut)
        elec_map = create_elec_map(ilut)
        ! now map eletrons: start with the first electron
        ! account for spin-conservation: is nel+1 odd or even?
        call assert_equals(1,elec_map(1))
        ! the second
        call assert_equals(4,elec_map(2))
        ! the third
        call assert_equals(5,elec_map(3))
        ! the fourth
        call assert_equals(8,elec_map(4))
        ! the fifth
        call assert_equals(9,elec_map(5))

    end subroutine test_elec_mapping

    subroutine random_fcidump(iunit)
        integer, intent(in) :: iunit
        call generate_random_integrals(&
            iunit, n_el=5, n_spat_orb=12, sparse=0.9_dp, sparseT=0.1_dp, total_ms=beta)
    end subroutine random_fcidump

    function calc_pgen_local(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
            use constants
            use SymExcitDataMod, only: scratchSize
            use bit_rep_data, only: NIfTot
            use SystemData, only: nel
            implicit none
            integer, intent(in) :: nI(nel)
            integer(n_int), intent(in) :: ilutI(0:NIfTot)
            integer, intent(in) :: ex(2, 2), ic
            integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)

            real(dp) :: pgen

            unused_var(nI)
            unused_var(ClassCount2)
            unused_var(ClassCountUnocc2)

            pgen = calc_pgen_pcpp(ilutI, ex, ic)

        end function calc_pgen_local

end program test_pcpp_excitgen
