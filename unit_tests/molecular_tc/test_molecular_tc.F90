! unit tests for molecular tc
#include "macros.h"
program test_molecular_tc
  use constants
  use SystemData, only: nBasis, tStoreSpinOrbs, nel, G1, nullBasisFn
  use bit_rep_data, only: NIfTot
  use LMat_mod
  use fruit
  
  implicit none

  call init_fruit()
  
  call molecular_tc_test_driver()
  
  call fruit_summary()
  call fruit_finalize()

  contains

    subroutine molecular_tc_test_driver()
      implicit none
      call setup_tests()

      call run_lmat_test()
      call run_sltcnd_test()
      call run_excitgen_test()

      call clear_resources()
    end subroutine molecular_tc_test_driver
    
    subroutine setup_tests()
      use procedure_pointers, only: get_umat_el
      use sltcnd_mod, only: nullUMat, initSltCndPtr
      use SystemData, only: nOccAlpha, nOccBeta, AA_elec_pairs, AB_elec_pairs, &
           BB_elec_pairs, par_elec_pairs, tNoSymGenRandExcits, t_mol_3_body
      use dSFMT_interface, only: dSFMT_init
      use FciMCData, only: pParallel
      use tc_three_body_data, only: pTriples
      implicit none

      integer :: i

      nBasis = 14
      tStoreSpinOrbs = .false.
      nel = 3

      allocate(G1(nBasis))
      G1 = nullBasisFn

      do i = 1, nBasis
         if(is_beta(i)) then
            G1(i)%ms = -1
         else
            G1(i)%ms = 1
         end if
      end do
      
      get_umat_el => nullUMat

      t_mol_3_body = .true.
      call initSltCndPtr()

      call init_dummy()

      nOccAlpha = 2
      nOccBeta = 1
      pParallel = 0.5
      AA_elec_pairs = nOccAlpha**2
      BB_elec_pairs = nOccBeta**2
      AB_elec_pairs = nOccAlpha*nOccBeta
      par_elec_pairs = AA_elec_pairs + BB_elec_pairs

      tNoSymGenRandExcits = .true.

      NIfTot = 2

      call dSFMT_init(4)
      
      pTriples = 1.0
    end subroutine setup_tests

    subroutine run_lmat_test()
      implicit none

      call readLMat()
    end subroutine run_lmat_test
   
    subroutine run_sltcnd_test()
      use sltcnd_mod
      use, INTRINSIC :: IEEE_ARITHMETIC
      implicit none
      
      integer :: nI(nel), nJ(nel)
      HElement_t(dp) :: matel
      
      nI = (/2,4,6/)

      nJ = (/8,10,14/)
      matel = sltcnd_compat(nI,nJ,3)
      print *, "Triple matrix element", matel
      
      nJ = (/4,10,14/)
      print *, "Double matrix element", sltcnd_compat(nI,nJ,2)

      nJ = (/2,4,14/)
      print *, "Single matrix element", sltcnd_compat(nI,nJ,1)

      print *, "Diagonal matrix element", sltcnd_compat(nI,nI,0)
    end subroutine run_sltcnd_test

    subroutine run_excitgen_test()
      use tc_three_body_excitgen
      use tc_three_body_data, only: pgen3B, pgen2B, pgen1B, pgen0B, pTriples
      use DeterminantData, only: write_det
      implicit none

      integer :: nI(nel), nJ(nel), exFlag, ic, ExcitMat(2,3)
      logical :: tParity
      real(dp) :: pgen
      HElement_t(dp) :: helgen
      type(excit_gen_store_type), target :: store
      integer(n_int) :: ilut(0:NIfTot), ilutJ(0:NIfTot)
     
      nI = (/1,2,4/)
      call setup_mol_tc_excitgen(nI)
      pTriples = 1.0

      call assert_equals(pgen0B,0.0_dp)
      call assert_equals(pgen1B,1.0/60.0_dp)
      call assert_equals(pgen2B,0.0_dp)
      call assert_equals(pgen3B,0.0_dp)

      ilut = 0_n_int
      ilutJ = 0_n_int
      call gen_excit_mol_tc(nI, ilut, nJ, ilutJ, exFlag, ic, ExcitMat, &
           tParity, pgen, helgen, store)

      print *, "Generated: ", nJ
      print *, "Prob: ", pgen

      call assert_equals(pgen,calc_pgen_triple(nI,ExcitMat))
    end subroutine run_excitgen_test

    subroutine init_dummy()
      use OneEInts, only: TMat2D, tOneElecDiag
      use UMatCache, only: UMat2D, UMat2DExch
      use SymData
      use SymExcitDataMod, only: SpinOrbSymLabel
      implicit none
      tOneElecDiag = .false.

      allocate(TMat2D(nBasis,nBasis))
      TMat2D = 0.0

      allocate(Symclasses(nBasis))
      Symclasses = 1
      allocate(symlabelintscum(1))
      symlabelintscum = 0
      allocate(StateSymMap(nBasis))
      StateSymMap = 1

      allocate(UMat2D(nBasis,nBasis))
      allocate(UMat2DExch(nBasis,nBasis))
      UMat2D = 0.0
      UMat2DExch = 0.0

      allocate(spinorbsymlabel(nBasis), source = 0)
      
    end subroutine init_dummy

    subroutine clear_resources()
      use OneEInts, only: TMat2D
      use UMatCache, only: UMat2D, UMat2DExch
      use SymData
      implicit none
      
      deallocate(UMat2DExch)
      deallocate(UMat2D)
      deallocate(StateSymMap)
      deallocate(symlabelintscum)
      deallocate(Symclasses)
      deallocate(TMat2D)
    end subroutine clear_resources

    
end program test_molecular_tc
