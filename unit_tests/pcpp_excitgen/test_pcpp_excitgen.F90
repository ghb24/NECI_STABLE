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
    implicit none
    real(dp) :: pTot, pNull
    integer, parameter :: nSamples = 500000

    ! set the excitation we want to test
    generate_excitation => gen_rand_excit_pcpp
    ! prepare everything for testing the excitgen
    call init_excitgen_test()

    ! prepare the pcpp excitation generator: get the precomputed weights
    call set_ref()    
    call init_pcpp_excitgen()

    ! test the mapping of electrons
    call test_elec_mapping()
    ! run the excitgen test: do nSamples excitations, and compare them with all possible excits
    call test_excitation_generator(nSamples,pTot,pNull)
    call free_ref()
    call finalize_excitgen_test()    
  end subroutine pcpp_test_driver

  subroutine test_elec_mapping()
    ! test that the mapping of electrons does what its supposed to do
    implicit none
    ! map to a determinant that has 2 electrons excited
    integer :: nI(nel)
    integer(n_int) :: ilut(0:NIfTot)
    integer :: ex(2,2)
    logical :: tPar
    ! create the excited determinant
    nI = projEDet(:,1)
    ex(1,1) = 1
    ex(1,2) = 2
    ex(2,1) = nel + 1
    ex(2,2) = nel + 2    
    call FindExcitDet(ex,nI,2,tPar)
    call EncodeBitDet(nI,ilut)
    ! now map eletrons: start with the first electron
    ! account for spin-conservation: is nel+1 odd or even?
    call assert_equals(nel+1+mod(nel,2),map_elec_from_ref(ilut,1))
    ! the second
    call assert_equals(nel+2-mod(nel,2),map_elec_from_ref(ilut,2))
    ! the third
    call assert_equals(projEDet(3,1),map_elec_from_ref(ilut,3))
    
  end subroutine test_elec_mapping
  
end program test_pcpp_excitgen
