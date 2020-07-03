#include "macros.h"

program test_impurity_excit_gen
    use impurity_models
    use fruit
    use unit_test_helper_excitgen

    implicit none

    integer, parameter :: nBasis_ = 6
    integer, parameter :: nel_ = 3

    call init_fruit()
    ! this is where the tests go
    call impurity_excit_gen_test_driver
    call fruit_summary
    call fruit_finalize

contains

    subroutine impurity_excit_gen_test_driver
        implicit none
        integer, parameter :: nSamples = 1000
        real(dp) :: pTot, pNull
        integer :: nFound, numEx
        call init_excitgen_test(nBasis_/2, nel_)
        call init_impurity_tests
        call test_excitation_generator(nSamples, pTot, pNull, numEx, nFound, .false.)
        call assert_equals(numEx,nFound)
        call assert_true(abs(1.0-pTot) < 0.001)
    end subroutine impurity_excit_gen_test_driver

    !------------------------------------------------------------------------------------------!

    subroutine init_impurity_tests
      use SystemData, only: nBasis, nel, t_complex_ints
      use OneEInts, only: tmat2d, tCPMDSymTMat, tOneElecDiag
      use procedure_pointers, only: get_umat_el
      use Integrals_neci, only: get_umat_el_tumat2d
      use bit_rep_data, only: niftot, nifdbo
      use dSFMT_interface, only: dSFMT_init
      implicit none
      t_complex_ints = .false.
      nel = 3      
      nBasis = 6
      niftot = 8
      nifdbo = 0
      tCPMDSymTMat = .false.
      tOneElecDiag = .false.
      allocate(tmat2d(6,6))
      tmat2d = 0.0
      tmat2d(1,3) = 0.1
      tmat2d(1,5) = 0.1
      tmat2d(2,4) = 0.1
      tmat2d(2,6) = 0.1
      tmat2d(3,1) = 0.1
      tmat2d(5,1) = 0.1
      tmat2d(4,2) = 0.1
      tmat2d(6,2) = 0.1
      tmat2d(1,1) = -0.5
      tmat2d(2,2) = -0.5
      tmat2d(3,3) = 0.3
      tmat2d(4,4) = 0.3
      tmat2d(5,5) = -0.3
      tmat2d(5,5) = -0.3

      get_umat_el => get_umat_el_test
      generate_excitation => gen_excit_impurity_model      
      call dSFMT_init(3)
      call setupImpurityExcitgen()
      
    end subroutine init_impurity_tests

    subroutine generate_test_ilut(ilut,nI)
        use bit_rep_data, only: niftot
        use SystemData, only: nel

        integer(n_int), intent(out) :: ilut(0:niftot)
        integer,intent(out) :: nI(nel)
        integer :: i

        nI = (/ 1,3,6/)
        ilut = 0
        do i=1,nel
            ilut(0) = ibset(ilut(0),nI(i)-1)
        enddo
    end subroutine generate_test_ilut


end program
