#include "macros.h"
module tc_precomputed_excitgen
  use constants
  use tc_three_body_data
  use aliasSampling
  use bit_reps, only: niftot
  use SystemData, only: nel, nBasis, G1, BRR, symmax, Symmetry
  use sym_mod, only: symprod, symconj
  use FciMCData, only: excit_gen_store_type, pDoubles, pSingles, projEDet, ilutRef
  use dSFMT_interface , only : genrand_real2_dSFMT
  use sltcnd_mod, only: sltcnd_excit
  use util_mod, only: binary_search_first_ge
  implicit none

  ! these are the samplers used for generating single excitations
  ! we have one hole sampler for each orbital
  type(aliasSampler_t), allocatable :: single_hole_sampler(:)
  ! we have a single electron sampler, the electron is always chosen first
  type(aliasSampler_t) :: single_elec_sampler

  ! these are the samplers used for generating double excitations
  type(aliasSampler_t) :: double_elec_one_sampler
  type(aliasSampler_t), allocatable :: double_elec_two_sampler(:)

  type(aliasSampler_t), allocatable :: double_hole_one_sampler(:)
  type(aliasSampler_t), allocatable :: double_hole_two_sampler(:,:)

  integer, allocatable :: refDet(:)

contains

  !------------------------------------------------------------------------------------------!
  ! Main routines generating the excitation
  !------------------------------------------------------------------------------------------!

  ! this is the interface routine that calls the corresponding generator to create
  ! a single/double or potentially triple excitation
  subroutine gen_ran_excit_pcpp(nI, ilut, nJ, ilutnJ, exFlag, ic, ExcitMat, tParity, pGen, &
       HElGen, store, part_type)
    implicit none
    ! The interface is common to all excitation generators, see proc_ptrs.F90
    integer, intent(in) :: nI(nel), exFlag
    integer(n_int), intent(in) :: iLut(0:niftot)
    integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(out) :: pGen
    type(excit_gen_store_type), intent(inout), target :: store

    ! Not used
    integer(n_int), intent(out) :: ilutnJ(0:niftot)
    HElement_t(dp), intent(out) :: HElGen

    integer, intent(in), optional :: part_type

    real(dp) :: r

    ! decide whether to generate a single or double excitation
    r = genrand_real2_dSFMT()
    if(r < pDoubles) then
       call generate_double_pcpp(nI, ilut, nJ, excitMat, tParity, pGen)
       IC = 2
       pGen = pGen * pDoubles
    else
       call generate_single_pcpp(nI, ilut, nJ, excitMat, tParity, pGen)
       IC = 1
       pGen = pGen * pSingles
    endif
    
  end subroutine gen_ran_excit_pcpp

  !------------------------------------------------------------------------------------------!

  subroutine generate_double_pcpp(nI, ilut, nJ, excitMat, tParity, pGen)
    implicit none
    ! given the initial determinant (both as nI and ilut), create a random single excitation
    ! given by nJ/ilutnJ/excitMat with probability pGen. tParity indicates the fermi sign
    ! picked up by applying the excitation operator
    integer, intent(in) :: nI(nel)
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    integer, intent(out) :: nJ(nel)
    integer, intent(out) :: ExcitMat(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(out) :: pGen
    real(dp) :: pSGen1, pSGen2, pTGen1, pTGen2
    integer :: src1, src2
    integer :: tgt1, tgt2
    integer :: elec1, elec2
    type(Symmetry) :: tgtSym
    ! mapping electrons in the reference to orbitals in the current det is not 1-to-1
    ! so we need to remove an orbital after it has been mapped to, this is what ilutEx is for
    integer(n_int) :: ilutEx(0:NIfTot)

    call double_elec_one_sampler%sample(src1,pSGen1)
    src1 = map_elec_from_ref(ilut, src1)
    call double_elec_two_sampler(src1)%sample(src2,pSGen2)

    ! we use ilutEx for the second mapping: remove the already-mapped-to orbital, so we
    ! cannot map to it a second time
    ilutEx = ilut
    clr_orb(ilutEx,src1)
    src2 = map_elec_from_ref(ilutEx, src2)

    if(src2 < src1) call intSwap(src1,src2)

    call double_hole_one_sampler(src1)%sample(tgt1,pTGen1)
    ! generation probability  so far to ensure it has a valid value on return in any case
    pGen = pSGen1 * pSGen2 * pTGen1
    if(abort_excit(tgt1)) return
    ! we need a specific symmetry now
    tgtSym = symprod(G1(src1)%Sym,G1(src2)%Sym)
    tgtSym = symprod(tgtSym,symconj(G1(tgt1)%Sym))
    call double_hole_two_sampler(src2,tgtSym%s)%sample(tgt2,pTGen2)

    ! Update the generation probability
    pGen = pGen * pTGen2
    if(abort_excit(tgt2)) return

    elec1 = binary_search_first_ge(nI,src1)
    elec2 = binary_search_first_ge(nI,src2)
    call make_double(nI, nJ, elec1, elec2, tgt1, tgt2, ExcitMat, tParity)
       
  contains

    function abort_excit(tgt) result(abort)
      ! check if the target orbital is valid
      ! Input: tgt - orbital we want to know about: Is an excitation to this possible
      ! Output: abort - true if there is no allowed excitation
      implicit none
      integer, intent(in) :: tgt
      logical :: abort

      abort = IsOcc(ilut,tgt)      
      if(abort) then
         nJ = 0
         ExcitMat = 0
         tParity = .false.
      endif
      
    end function abort_excit
  end subroutine generate_double_pcpp

  !------------------------------------------------------------------------------------------!

  subroutine generate_single_pcpp(nI, ilut, nJ, excitMat, tParity, pGen)
    implicit none
    ! given the initial determinant (both as nI and ilut), create a random double excitation
    ! given by nJ/ilutnJ/excitMat with probability pGen. tParity indicates the fermi sign
    ! picked up by applying the excitation operator
    integer, intent(in) :: nI(nel)
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    integer, intent(out) :: nJ(nel)
    integer, intent(out) :: ExcitMat(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(out) :: pGen

    integer :: src, elec, tgt
    real(dp) :: pHole

    ! get a random electron
    call single_elec_sampler%sample(src,pGen)
    
    ! map the electron to the current determinant    
    src = map_elec_from_ref(ilut, src)
    
    ! get a random associated orbital    
    call single_hole_sampler(src)%sample(tgt,pHole)

    if(IsOcc(ilut,tgt)) then
       ! invalidate the excitation, we hit an occupied orbital
       nJ = 0
       excitMat = 0
       tParity = .false.
    else
       elec = binary_search_first_ge(nI,src)
       call make_single(nI, nJ, elec, tgt, excitMat, tParity)
    endif
    ! add the probability to find this hole from this electron
    pGen = pGen * pHole
    
  end subroutine generate_single_pcpp

  !------------------------------------------------------------------------------------------!
  ! Functions that map orbital and electron indices between reference and current determinant
  !------------------------------------------------------------------------------------------!  

  pure function map_elec_from_ref(ilut, iElec) result(src)
    ! Transfer an electron index within the reference to an orbital index
    ! within the current determinant nI
    ! Input: nI - current determinant
    !        iElec - index of the electron to be mapped within the reference
    ! Output: src - orbital of the mapped electron
    implicit none
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    integer, intent(in) :: iElec
    integer :: src
    integer(n_int) :: ilutReference(0:NIfTot)
    integer :: refOrb
    integer :: iOrb

    ilutReference = ilutRef(:,1)
    refOrb = refDet(iElec)
    ! if the corresponding orbital is occupied in ilut, no mapping is required
    if(IsOcc(ilut,refOrb)) then
       src = refOrb
       return
    endif

    ! count through the orbitals in energy order
    ! - BRR from SystemData.F90 is a list of the orbitals in energy order
    do iOrb = 1, nBasis
       ! we pick the first orbital that is occupied in ilut, unoccupied in ref and
       ! has the same spin as the unmapped orbital
       if(IsOcc(ilut,BRR(iOrb)) .and. .not.IsOcc(ilutReference,BRR(iOrb)) .and. &
            same_spin(BRR(iOrb),refOrb)) then
          src = BRR(iOrb)
          return
       endif
    end do

  end function map_elec_from_ref

  pure function map_orb_from_ref(nI, iOrb) result(tgt)
    ! UNUSED
    implicit none
    integer, intent(in) :: nI(nel)
    integer, intent(in) :: iOrb
    integer :: tgt
  end function map_orb_from_ref
  
  !------------------------------------------------------------------------------------------!
  ! Initialization routines for the pcpp excitation generator
  !------------------------------------------------------------------------------------------!

  subroutine init_pcpp_excitgen()
    implicit none
    
    allocate(refDet(nel))
    refDet = projEDet(:,1)
    
    call init_pcpp_doubles_excitgen()
    call init_pcpp_singles_excitgen()
   
  end subroutine init_pcpp_excitgen
  
  !------------------------------------------------------------------------------------------!

  subroutine init_pcpp_doubles_excitgen()
    implicit none

    call setup_elec_one_sampler()
    call setup_elec_two_sampler()

    call setup_hole_one_sampler()
    call setup_hole_two_sampler()
  contains

    subroutine setup_elec_one_sampler()
      integer :: i
      integer :: a,b,j
      real(dp) :: w(nel)
      integer :: ex(2,2)
      logical :: tPar
      integer :: iEl, jEl

      do iEl = 1, nel
         w(iEl) = 0
         do jEl = 1, nel
            if(i.ne.j) then
               i = refDet(iEl)
               j = refDet(jEl)
               do a = 1, nBasis
                  if(.not.any(a.eq.(/i,j/))) then
                     do b = 1, nBasis
                        if(.not.any(b.eq.(/a,i,j/))) then
                           call set_ex(ex,i,j,a,b)
                           w(i) = w(i) + abs(sltcnd_excit(refDet,2,ex,tPar))
                        endif
                     end do
                  end if
               end do
            end if
         end do
      end do

      call apply_lower_bound(w)
      call double_elec_one_sampler%setupSampler(w)
    end subroutine setup_elec_one_sampler

    !------------------------------------------------------------------------------------------!

    subroutine setup_elec_two_sampler()
      implicit none
      real(dp) :: w(nel)
      integer :: ex(2,2)
      logical :: tPar
      integer :: aerr
      integer :: i,j,a,b
      integer :: jEl
      
      allocate(double_elec_two_sampler(nBasis))
      do i = 1, nBasis
         w = 0.0_dp
         do jEl = 1, nel
            j = refDet(jEl)
            if(i.ne.j) then
               do a = 1, nBasis
                  if(.not.any(a.eq.(/i,j/))) then
                     do b = 1, nBasis
                        if(.not.any(b.eq.(/i,j,a/))) then
                           call set_ex(ex,i,j,a,b)
                           w(jEl) = w(jEl) + abs(sltcnd_excit(refDet,2,ex,tPar))
                        end if
                     end do
                  end if
               end do
            end if
         end do

         ! to prevent bias, a lower bound for the probabilities is set
         call apply_lower_bound(w)

         ! except for picking the same electron twice, that is actually forbidden
         do jEl = 1, nEl
            if(refDet(jEl).eq.i) w(jEl) = 0.0_dp
         end do
         call double_elec_two_sampler(i)%setupSampler(w)
      end do
    end subroutine setup_elec_two_sampler
    !------------------------------------------------------------------------------------------!

    subroutine setup_hole_one_sampler()
      ! generate precomputed probabilities for picking a hole given a selected electron
      ! this is for picking the first hole where no symmetry restrictions apply
      implicit none
      real(dp) :: w(nBasis)
      integer :: i,a
      logical :: tPar
      integer :: aerr

      allocate(double_hole_one_sampler(nBasis), stat = aerr)
      do i = 1, nBasis
         w = 0.0_dp
         do a = 1, nBasis
            ! only same-spin excitations from i -> a
            if(a.ne.i.and.same_spin(a,i)) &
                 w(a) = pp_weight_function(i,a)
         end do
         call double_hole_one_sampler(i)%setupSampler(w)
      end do
    end subroutine setup_hole_one_sampler

    !------------------------------------------------------------------------------------------!

    subroutine setup_hole_two_sampler()
      ! generate precomputed probabilities for picking hole number 2 given a selected electron
      ! this is for picking the second hole where symmetry restrictions apply
      implicit none
      real(dp) :: w(nBasis,symmax)
      integer :: j,b,iSym
      logical :: tPar
      integer :: aerr        

      ! there is one table for each symmetry and each starting orbital
      allocate(double_hole_two_sampler(nBasis,symmax), stat = aerr)
      do j = 1, nBasis
         w = 0.0_dp
         do b = 1, nBasis
            ! only same-spin and symmetry-allowed excitations from j -> b
            if(b.ne.j.and.same_spin(b,j)) &
                 w(b,G1(b)%Sym%s) = pp_weight_function(j,b)
         end do

         do iSym = 1, symmax
            call double_hole_two_sampler(j,iSym)%setupSampler(w(:,iSym))
         end do
      end do

    end subroutine setup_hole_two_sampler

    !------------------------------------------------------------------------------------------!

    pure subroutine set_ex(ex,i,j,a,b)
      implicit none
      integer, intent(out) :: ex(2,2)
      integer, intent(in) :: i,j,a,b

      ex(1,1) = refDet(i)
      ex(1,2) = refDet(j)
      ex(2,1) = a
      ex(2,2) = b
    end subroutine set_ex

  end subroutine init_pcpp_doubles_excitgen

  !------------------------------------------------------------------------------------------!

  subroutine init_pcpp_singles_excitgen()
    ! create the probability distributions for drawing single excitations
    ! The normalization 1/n_{jb} used by Neufeld and Thom seems to be just a constant, we
    ! omit it for now
    implicit none

    call setup_elecs_sampler()
    call setup_holes_sampler()
    
  contains

  !------------------------------------------------------------------------------------------!

    subroutine setup_elecs_sampler()
      ! the probability distribution for selection of electrons
      ! creates a sampler for electron indices within the reference determinant
      ! these later have to be transferred to the current determinant
      
      ! even though strictly speaking, these are sums of the hole probabilities,
      ! expressing them in terms of the latter would be an unwanted dependency,
      ! since this relation is merely a matter of choice of the algorithm and should not
      ! be reflected in the design of the implementation
      implicit none
      real(dp) :: w(nel)
      integer :: i, a
      integer :: aerr
      integer :: refOrb
      
      do i = 1, nel
         w(i) = 0
         refOrb = refDet(i)
         do a = 1, nBasis
            if(G1(a)%Sym%s.eq.G1(refOrb)%Sym%s) then
               w(i) = w(i) + acc_doub_matel(refOrb,a)
            end if
         end do
      end do
      ! load the probabilites for electron selection into the alias table
      call single_elec_sampler%setupSampler(w)
      
    end subroutine setup_elecs_sampler

    !------------------------------------------------------------------------------------------!

    subroutine setup_holes_sampler()
      ! the probability distributions for selection of holes, given the electron orbital
      implicit none
      real(dp) :: w(nBasis)
      integer :: i,a
      integer :: aerr

      ! there is one table for given source orbital
      allocate(single_hole_sampler(nBasis), stat = aerr)

      ! each table has probabilities for all given virtual orbitals a (some of them might
      ! be 0, this has no additional cost)    
      do i = 1, nBasis
         w = 0.0_dp
         do a = 1, nBasis
            ! store the accumulated matrix elements (= un-normalized probability) with
            ! the corresponding symmetry (if spins of a/i are different, w is 0)
            w(a) = acc_doub_matel(i,a)
         end do

         call single_hole_sampler(i)%setupSampler(w)
      end do

    end subroutine setup_holes_sampler

    !------------------------------------------------------------------------------------------!
    
    function acc_doub_matel(src, tgt) result(prob)
      ! Accumulate all single excitation matrix elements connecting
      ! D_(j,src)^(b,tgt) for all (j,b), where D is the reference
      ! Input: src - orbital we want to excite from (any orbital is allowed)
      !        tgt - orbital we want to excite to (any orbital is allowed)
      ! Output: prob - Accumulated matrix elements of singles
      !                If src/tgt have different spin, the result is 0, if they have
      !                different symmetry, it is not necessarily
      !
      ! IMPORTANT: The matrix elements are calculated as
      ! <D_j^b|H|D_(j,src)^(b,tgt)> := h_src^tgt + sum_(k occ in D_j^b) <tgt k|src k> - <tgt k|k src>
      ! That is, the left side is to be understood symbolic, there is no actual excitation
      ! from src to tgt
      implicit none
      integer, intent(in) :: src, tgt
      real(dp) :: prob
      integer :: b, j
      integer :: ex(2,2), nI(nel)
      logical :: tPar

      prob = 0

      if(same_spin(src,tgt)) then
         do b = 1, nBasis
            ! loop over all non-occupied orbitals
            if(.not.any(b.eq.refDet(:))) then
               do j = 1, nel
                  ! get the excited determinant D_j^b used for the matrix element
                  call make_single(refDet(:), nI, j, b, ex, tPar)

                  ! this is a symbolic excitation, we do NOT require src to be occupied
                  ! we just use the formula for excitation matrix elements
                  ex(1,1) = src
                  ex(2,1) = tgt
                  prob = prob + abs(sltcnd_excit(nI,1,ex,tPar))
               end do
            endif
         end do
      endif
      
    end function acc_doub_matel
  !------------------------------------------------------------------------------------------!
  end subroutine init_pcpp_singles_excitgen

  !------------------------------------------------------------------------------------------!
  ! Finalization routines
  !------------------------------------------------------------------------------------------!
  
  subroutine finalize_pcpp_excitgen()
    implicit none
    integer :: j
    deallocate(refDet)

    call single_elec_sampler%samplerDestructor()
    
    call clear_sampler_array(single_hole_sampler)
    call double_elec_one_sampler%samplerDestructor()
    call clear_sampler_array(double_elec_two_sampler)
    call clear_sampler_array(double_hole_one_sampler)
    do j = 1, size(double_hole_two_sampler)
       call clear_sampler_array(double_hole_two_sampler(:,j))
    end do
  contains

    subroutine clear_sampler_array(arr)
      ! call the destructor on all elements of an array, then deallocate it
      type(aliasSampler_t), allocatable :: arr(:)

      integer :: i
      
      do i = 1, size(arr)
         call arr(i)%samplerDestructor()
      end do
      deallocate(arr)
    end subroutine clear_sampler_array
    
  end subroutine finalize_pcpp_excitgen
  

  !------------------------------------------------------------------------------------------!
  ! Auxiliary functions
  !------------------------------------------------------------------------------------------!

  pure elemental function gtSpin(spinOrb) result(spinInd)
    ! get a spin index encoding the spin of the input orbital
    ! Input: spinOrb - spin orbital which spin index to get
    ! Output: spinInd - 1 if MS is +1
    !                   2 if MS is -1
    implicit none
    integer, intent(in) :: spinOrb
    integer :: spinInd

    spinInd = mod(spinOrb,2)
  end function gtSpin

  !------------------------------------------------------------------------------------------!

  pure subroutine apply_lower_bound(w)
    ! even if some excitation is not possible in the reference frame, it might be in the
    ! current determinant, so we enforce a minimum value on the probabilities
    ! Input/Output: w - on input, list of probabilities
    !                 - on output, same list with a minimal value enforced
    implicit none
    real(dp), intent(inout) :: w(:)
    real(dp) :: mVal
    integer :: i

    mVal = 0.001*minVal(w,w>eps)
    do i = 1, size(w)
       w(i) = max(w(i),mVal)
    end do
    
  end subroutine apply_lower_bound

  !------------------------------------------------------------------------------------------!
  
  function pp_weight_function(i,a) result(w)
    ! Given an excitation, return the power-pitzer weights
    ! Can be tweaked to handle 3-body excitations
    ! Input: i - selected electron
    !        a - possible orbital to excite to
    ! Output: w - approximate weight of this excitation
    implicit none
    
    integer, intent(in) :: i,a
    real(dp) :: w
    integer :: ex(2,2)

    ex(1,1) = i
    ex(1,2) = a
    ex(2,1) = i
    ex(2,2) = a
    w = sqrt(abs(sltcnd_excit(refDet,2,ex,.false.)))
  end function pp_weight_function
  
!------------------------------------------------------------------------------------------!

    pure subroutine intswap(a,b)
      integer, intent(inout) :: a,b
      integer :: tmp
      
      tmp = a
      a = b 
      b = tmp
    end subroutine intswap
  
end module tc_precomputed_excitgen
