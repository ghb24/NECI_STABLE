program test_pcpp_excitgen
  use constants
  use fruit
  use Parallel_neci, only: MPIInit, MPIEnd
  use pchb_excitgen, only: gen_rand_excit_pchb, PCHB_FCI
  use unit_test_helper_excitgen
  use orb_idx_mod, only: beta
  use procedure_pointers, only: generate_excitation
  implicit none

  call MPIInit(.false.)
  call init_fruit()
  call pchb_test_driver()
  call fruit_summary()
  call fruit_finalize()
  call MPIEnd(.false.)

contains

  subroutine pchb_test_driver()
    implicit none
    real(dp) :: pTot, pNull
    integer :: numEx, nFound
    ! There can be some excitations with really low matrix elements -> we need a lot
    ! of samples to hit all
    integer, parameter :: nSamples = 1000000

    ! set the excitation generator to pchb
    generate_excitation => gen_rand_excit_pchb
    calc_pgen => calc_pgen_pchb

    ! prepare an excitation generator test
    call init_excitgen_test(n_el=5, fcidump_writer=FciDumpWriter_t(random_fcidump, 'FCIDUMP'))

    ! prepare the pchb excitgen: set the weights/map-table
    call set_ref()
    call PCHB_FCI%init()

    ! test the excitation generator
    call test_excitation_generator(nSamples,pTot,pNull,numEx,nFound,.true.)
    ! make sure all excits have been found
    call assert_equals(numEx,nFound)
    ! and the total prob is 1.0
    call assert_true(abs(1.0-pTot) < 0.05)

    ! free memory
    call free_ref()
    call PCHB_FCI%finalize()
  end subroutine pchb_test_driver

  function calc_pgen_pchb(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
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

      pgen = PCHB_FCI%calc_pgen(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
  end function

  subroutine random_fcidump(iunit)
    integer, intent(in) :: iunit
    call generate_random_integrals(&
      iunit, n_el=5, n_spat_orb=12, sparse=0.9_dp, sparseT=0.1_dp, total_ms=beta)
  end subroutine
end program test_pcpp_excitgen
