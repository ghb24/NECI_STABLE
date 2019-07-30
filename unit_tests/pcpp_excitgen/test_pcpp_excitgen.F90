program test_pcpp_excitgen
  use constants
  use fruit
  use pcpp_excitgen
  use unit_test_helper_excitgen
  use procedure_pointers, only: generate_excitation
  implicit none

  call init_fruit()
  call pcpp_test_driver()
  call fruit_summary()
  call fruit_finalize()
  
contains

  !------------------------------------------------------------------------------------------!

  subroutine pcpp_test_driver()
    real(dp) :: pTot, pNull
    integer, parameter :: nSamples = 1000

    ! set the excitation we want to test
    generate_excitation => gen_rand_excit_pcpp
    ! prepare everything for testing the excitgen
    call init_excitgen_test()

    ! prepare the pcpp excitation generator: get the precomputed weights
    call set_ref()    
    call init_pcpp_excitgen()

    ! run the test: do nSamples excitations, and compare them with all possible excits
    call test_excitation_generator(nSamples,pTot,pNull)
    call free_ref()
    call finalize_excitgen_test()    
  end subroutine pcpp_test_driver
  
end program test_pcpp_excitgen
