#include "macros.h"

program test_impurity_excit_gen
  use impurityModels
  use fruit
  use fruit_util

  implicit none

  call init_fruit()
  ! this is where the tests go
  call impurity_excit_gen_test_driver
  call fruit_summary
  call fruit_finalize

  contains

    subroutine impurity_excit_gen_test_driver
      implicit none
      call init_impurity_tests
      call test_orbital_lists
      call test_bath_excitation
    end subroutine impurity_excit_gen_test_driver

!------------------------------------------------------------------------------------------!

    subroutine init_impurity_tests
      use SystemData, only: nBasis, nel, tReltvy, tExch, G1, basisFN
      use OneEInts, only: tmat2d, tCPMDSymTMat, tOneElecDiag
      use procedure_pointers, only: get_umat_el
      use Integrals_neci, only: get_umat_el_tumat2d
      use bit_rep_data, only: niftot, tuseflags, nifdbo
      use dSFMT_interface, only: dSFMT_init
      implicit none
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

      tReltvy = .false.
      tExch = .false.
      allocate(G1(nBasis))
      G1(1)%Ms=1
      G1(2)%Ms=-1
      G1(3)%Ms=1
      G1(4)%Ms=-1
      G1(5)%Ms=1
      G1(6)%Ms=-1

      tuseflags = .false.

      get_umat_el => get_umat_el_test
      call dSFMT_init(3)

    end subroutine init_impurity_tests

!------------------------------------------------------------------------------------------!

    subroutine test_orbital_lists
      use FciMCData, only: nBath, nImp
      logical :: isBath(nBasis)

      call constructBath(isBath)
      call assert_true(.not. isBath(1))
      call assert_true(.not. isBath(2))
      call assert_true(isBath(3))
      call assert_true(isBath(4))
      call assert_true(isBath(5))
      call assert_true(isBath(6))

      call generateOrbitalLists(isBath)

      call assert_equals(nBath,4)
      call assert_equals(nImp,2)
    end subroutine test_orbital_lists

    subroutine test_bath_excitation
      use constants, only: dp, lenof_sign
      use bit_rep_data, only: niftot
      implicit none
      integer(n_int) :: ilut(0:niftot)
      integer :: dest
      real(dp) ::  pGen
      call generate_test_ilut(ilut)
      call hamiltonian_weighted_pick_single_bath(1,dest,pGen,ilut)
      call assert_true(dest .ge. 3)
    end subroutine test_bath_excitation

!------------------------------------------------------------------------------------------!

    function get_umat_el_test(a,b,c,d) result(uel)
      use constants, only: dp
      implicit none
      integer, intent(in) :: a,b,c,d
      HElement_t(dp) :: uel
      
      uel = 0.0
      if(((a .eq. 1) .and. (a .eq. c)) .and. (b .eq. d) .and. (b .eq. 2)) uel = 1
      if(((a .eq. 1) .and. (a .eq. d)) .and. (b .eq. c) .and. (b .eq. 2)) uel = -1
      if(((a .eq. 2) .and. (a .eq. c)) .and. (b .eq. d) .and. (b .eq. 1)) uel = 1
      if(((a .eq. 2) .and. (a .eq. d)) .and. (b .eq. c) .and. (b .eq. 1)) uel = -1
    end function get_umat_el_test

    subroutine generate_test_ilut(ilut)
      use bit_rep_data, only: niftot
      integer(n_int), intent(out) :: ilut(0:niftot)
      integer :: nI(nel), i

      nI = (/ 1,3,6/)
      ilut = 0
      do i=1,nel
         ilut(0) = ibset(ilut(0),nI(i)-1)
      enddo
    end subroutine generate_test_ilut
end program test_impurity_excit_gen
