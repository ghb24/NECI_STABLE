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
    integer :: numEx, nFound
    integer, parameter :: nSamples = 100000

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
    call test_excitation_generator(nSamples,pTot,pNull,numEx,nFound,.false.)
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

end program test_pcpp_excitgen
