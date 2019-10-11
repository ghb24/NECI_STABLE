program test_pcpp_excitgen
  use constants
  use fruit
  use pchb_excitgen
  use unit_test_helper_excitgen
  use procedure_pointers, only: generate_excitation
  implicit none

  call init_fruit()
  call pchb_test_driver()
  call fruit_summary()
  call fruit_finalize()

contains

  subroutine pchb_test_driver()
    implicit none
    real(dp) :: pTot, pNull
    integer :: numEx, nFound
    ! There can be some excitations with really low matrix elements -> we need a lot
    ! of samples to hit all
    integer, parameter :: nSamples = 10000000

    ! set the excitation generator to pchb
    generate_excitation => gen_rand_excit_pchb
    calc_pgen => calc_pgen_pchb
    ! prepare an excitation generator test
    call init_excitgen_test()

    ! prepare the pchb excitgen: set the weights/map-table
    call set_ref()
    call init_pchb_excitgen()

    ! test the excitation generator
    call test_excitation_generator(nSamples,pTot,pNull,numEx,nFound,.true.)
    ! make sure all excits have been found
    call assert_equals(numEx,nFound)
    ! and the total prob is 1.0
    call assert_true(abs(1.0-pTot) < 0.01)

    ! free memory
    call free_ref()
    call finalize_excitgen_test()
  end subroutine pchb_test_driver

end program test_pcpp_excitgen
