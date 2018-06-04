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
      call test_imp_excitation
      call test_single_excitation
      call test_gen_excit_impurity_model
    end subroutine impurity_excit_gen_test_driver

!------------------------------------------------------------------------------------------!

    subroutine init_impurity_tests
      use SystemData, only: nBasis, nel
      use OneEInts, only: tmat2d, tCPMDSymTMat, tOneElecDiag
      use procedure_pointers, only: get_umat_el
      use Integrals_neci, only: get_umat_el_tumat2d
      use bit_rep_data, only: niftot, nifdbo
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

      get_umat_el => get_umat_el_test
      call dSFMT_init(3)

    end subroutine init_impurity_tests

!------------------------------------------------------------------------------------------!

    subroutine test_imp_excitation
      use constants, only: dp
      use bit_rep_data, only: niftot
      use SystemData, only: nel
      implicit none
      integer(n_int) :: ilut(0:niftot)
      integer :: dest, nI(nel)
      real(dp) :: pGen
      
      dest = 0
      call generate_test_ilut(ilut,nI)
      call hamiltonian_weighted_pick_single(1,dest,pGen,ilut)
      call assert_true((dest .eq. 1) .and. (dest > 0))
    end subroutine test_imp_excitation

!------------------------------------------------------------------------------------------!

    subroutine test_single_excitation
      use constants, only: dp
      use bit_rep_data, only: niftot
      use DetBitOps, only: FindBitExcitLevel
      use SystemData, only: nel

      implicit none
      integer(n_int) :: ilut(0:niftot), ilutnJ(0:niftot)
      integer :: nI(nel), nJ(nel), ExcitMat(2,2)
      logical :: tParity
      real(dp) :: pGen
      
      call generate_test_ilut(ilut,nI)
      call generate_imp_single_excitation(nI,ilut,nJ,ilutnJ,ExcitMat,tParity,pGen)
      call assert_true(FindBitExcitLevel(ilut,ilutnJ)==1)
      call assert_true(ExcitMat(2,1) <= nImp)
    end subroutine test_single_excitation

!------------------------------------------------------------------------------------------!

    subroutine test_double_excitation
      use constants, only: dp
      use bit_rep_data, only: niftot
      use SystemData, only: nel
      use DetBitOps, only: FindBitExcitLevel
      
      implicit none
      integer(n_int) :: ilut(0:niftot), ilutnJ(0:niftot)
      integer :: nI(nel), nJ(nel), ex(2,2)
      logical :: tParity
      real(dp) :: pGen

      call generate_test_ilut(ilut,nI)
      call generate_imp_double_excitation(nI,ilut,nJ,ilutnJ,ex,tParity,pGen)
      call assert_true(FindBitExcitLevel(ilut,ilutnJ)==2)
      call assert_true(ex(2,1) <= nImp)
      call assert_true(ex(2,2) <= nImp)
    end subroutine test_double_excitation

!------------------------------------------------------------------------------------------!

    subroutine test_gen_excit_impurity_model
      use constants, only: dp
      use bit_rep_data, only: niftot
      use SystemData, only: nel
      use FciMCData, only: pSingles, excit_gen_store_type
      use DetBitOps, only: FindBitExcitLevel
      implicit none

      integer(n_int) :: ilut(0:niftot), ilutnJ(0:niftot)
      integer :: nI(nel), nJ(nel), ex(2,2), IC, exFlag
      logical :: tParity
      real(dp) :: pGen
      HElement_t(dp) :: HElGen
      type(excit_gen_store_type), target :: store

      call generate_test_ilut(ilut,nI)
      pSingles = 0.0
      call gen_excit_impurity_model(nI,ilut,nJ,ilutnJ,exFlag,IC,ex,tParity,pGen,&
           HElGen, store)
      ! our test model does not support double excitations, the impurity is too small
      call assert_true(FindBitExcitLevel(ilut,ilutnJ)==1)
    end subroutine test_gen_excit_impurity_model

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

!------------------------------------------------------------------------------------------!

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


end program test_impurity_excit_gen
