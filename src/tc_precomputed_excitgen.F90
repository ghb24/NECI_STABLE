#include "macros.h"
module pcpp_excitgen
  use constants
  use tc_three_body_data
  use aliasSampling
  use bit_reps, only: niftot
  use SystemData, only: nel, nBasis, G1, BRR, symmax, Symmetry
  use sym_mod, only: symprod, symconj
  use DetBitOps, only: EncodeBitDet
  use FciMCData, only: excit_gen_store_type, pDoubles, pSingles, projEDet, ilutRef
  use dSFMT_interface , only : genrand_real2_dSFMT
  use Integrals_neci, only: get_umat_el
  use UMatCache, only: gtID
  use sltcnd_mod, only: sltcnd_excit
  use util_mod, only: binary_search_first_ge
  use get_excit, only: make_double, make_single
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
  subroutine gen_rand_excit_pcpp(nI, ilut, nJ, ilutnJ, exFlag, ic, ExcitMat, tParity, pGen, &
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
    integer :: elec_map(nel)

    ! create the map for the electrons
    if(.not.store%tFilled) then
       store%elec_map = create_elec_map(ilut)
       store%tFilled = .true.
    endif
    elec_map = store%elec_map
    ! decide whether to generate a single or double excitation
    r = genrand_real2_dSFMT()
    if(r < pDoubles) then
       call generate_double_pcpp(nI, elec_map, ilut, nJ, excitMat, tParity, pGen)
       IC = 2
       pGen = pGen * pDoubles
    else
       call generate_single_pcpp(nI, elec_map, ilut, nJ, excitMat, tParity, pGen)
       IC = 1
       pGen = pGen * pSingles
    endif

    ! assign ilutnJ
    call EncodeBitDet(nJ,ilutnJ)

  end subroutine gen_rand_excit_pcpp

  !------------------------------------------------------------------------------------------!

  subroutine generate_double_pcpp(nI, elec_map, ilut, nJ, excitMat, tParity, pGen)
    implicit none
    ! given the initial determinant (both as nI and ilut), create a random single excitation
    ! given by nJ/ilutnJ/excitMat with probability pGen. tParity indicates the fermi sign
    ! picked up by applying the excitation operator
    integer, intent(in) :: nI(nel)
    integer, intent(in) :: elec_map(nel)
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    integer, intent(out) :: nJ(nel)
    integer, intent(out) :: ExcitMat(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(out) :: pGen
    ! temporary storage for the probabilities in each step
    real(dp) :: pSGen1, pSGen2, pTGen1, pTGen2
    ! temporary storage for the probabilities in the swapped steps
    real(dp) :: pSSwap1, pSSwap2, pTSwap1, pTSwap2
    ! temporary storage for the unmapped electrons
    integer :: umElec1, umElec2, swapOrb1
    ! chosen source/target orbitals
    integer :: src1, src2
    integer :: tgt1, tgt2
    integer :: elec1, elec2
    ! symmetry to enforce for the last orbital
    type(Symmetry) :: tgtSym
    ! map between reference and current der for electrons
    character(*), parameter :: t_r = "generate_double_pcpp"

    call double_elec_one_sampler%sample(umElec1,pSGen1)
    src1 = elec_map(umElec1)
    
    ! in very rare cases, no mapping is possible in the first place
    ! then, abort
    if(invalid_mapping(src1)) then
       ! the pgen here is just for bookkeeping purpose, it is never
       ! used
       pGen = pSGen1
       return
    endif
    
    call double_elec_two_sampler(src1)%sample(umElec2,pSGen2)
    src2 = elec_map(umElec2)
    
    ! it is possible to not be able to map the second electron if
    ! the first mapping occupied the only available slot
    if(invalid_mapping(src2,src1)) then
       pGen = pSGen1 * pSGen2
       return
    endif

    if(src2 < src1) call intSwap(src1,src2)

    ! we could also have drawn the electrons the other way around
    pSSwap1 = double_elec_one_sampler%getProb(umElec2)
    swapOrb1 = elec_map(umElec2)
    pSSwap2 = double_elec_two_sampler(swapOrb1)%getProb(umElec1)
    ! we now have industuingishible src1/src2, add the probabilites
    ! for drawing them either way
    pGen = pSGen1 * pSGen2 + pSSwap1 * pSSwap2

    call double_hole_one_sampler(src1)%sample(tgt1,pTGen1)
    ! generation probability  so far to ensure it has a valid value on return in any case
    if(abort_excit(tgt1)) then
       pGen = pGen * pTGen1
       return
    endif
    ! we need a specific symmetry now
    tgtSym = getTgtSym(tgt1)
    call double_hole_two_sampler(src2,tgtSym%s)%sample(tgt2,pTGen2)

    if(abort_excit(tgt2,tgt1)) then
       pGen = pGen * pTGen1 * pTGen2
       return
    endif
    ! Update the generation probability
    ! We could have drawn the target orbitals the other way around
    ! -> adapt pGen
    pTSwap1 = double_hole_one_sampler(src1)%getProb(tgt2)
    tgtSym = getTgtSym(tgt2)
    pTSwap2 = double_hole_two_sampler(src2,tgtSym%s)%getProb(tgt1)
    pGen = pGen * (pTGen1 * pTGen2 + pTSwap1 * pTSwap2)

    ! generate the output determinant
    elec1 = binary_search_first_ge(nI,src1)
    elec2 = binary_search_first_ge(nI,src2)
    call make_double(nI, nJ, elec1, elec2, tgt1, tgt2, ExcitMat, tParity)

  contains

    function getTgtSym(tgt) result(sym)
      ! return the symmetry of the last target orbital given a first
      ! target orbital tgt
      ! Input: tgt - first target orbital
      ! Output: sym - symmetry of the missing orbital
      implicit none
      integer, intent(in) :: tgt
      type(symmetry) :: sym

      sym= symprod(G1(src1)%Sym,G1(src2)%Sym)
      sym = symprod(sym,symconj(G1(tgt)%Sym))
      
    end function getTgtSym

    function invalid_mapping(src,src2) result(abort)
      ! check if the mapping was successful
      ! Input: src - electron we want to know about: did the mapping succeed?
      !        src2 - other chosen electron
      ! Output: abort - true if the mapping failed
      implicit none

      integer, intent(in) :: src
      integer, optional, intent(in) :: src2
      logical :: abort

      abort = src.eq.0
      if(present(src2)) abort = abort .or. (src.eq.src2)
      if(abort) then
         nJ = 0
         tParity = .false.
         ExcitMat = 0
         ExcitMat(1,1) = src
         if(present(src2)) ExcitMat(1,2) = src2
      endif
      
    end function invalid_mapping

    function abort_excit(tgt,tgt2) result(abort)
      ! check if the target orbital is valid
      ! Input: tgt - orbital we want to know about: Is an excitation to this possible
      !        tgt2 - second target orbital if already obtained
      ! Output: abort - true if there is no allowed excitation
      implicit none
      integer, intent(in) :: tgt
      integer, optional, intent(in) :: tgt2
      logical :: abort
      
      abort = IsOcc(ilut,tgt) .or. (tgt.eq.0)
      if(present(tgt2)) abort = abort .or. tgt==tgt2
      if(abort) then
         nJ = 0
         ExcitMat(1,1) = src1
         ExcitMat(1,2) = src2
         ExcitMat(2,1) = tgt         
         if(present(tgt2)) then
            ExcitMat(2,2) = tgt2
         else
            ExcitMat(2,2) = 0
         endif
         tParity = .false.
      endif

    end function abort_excit
  end subroutine generate_double_pcpp

  !------------------------------------------------------------------------------------------!

  subroutine generate_single_pcpp(nI, elec_map, ilut, nJ, excitMat, tParity, pGen)
    implicit none
    ! given the initial determinant (both as nI and ilut), create a random double excitation
    ! given by nJ/ilutnJ/excitMat with probability pGen. tParity indicates the fermi sign
    ! picked up by applying the excitation operator
    integer, intent(in) :: nI(nel)
    integer, intent(in) :: elec_map(nel)
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
    src = elec_map(src)

    ! get a random associated orbital    
    call single_hole_sampler(src)%sample(tgt,pHole)

    if(IsOcc(ilut,tgt)) then
       ! invalidate the excitation, we hit an occupied orbital
       nJ = 0
       excitMat = 0
       excitMat(1,1) = src
       excitMat(2,1) = tgt
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

  pure function create_elec_map(ilut) result(map)
    ! Create a map to transfer orbitals between the current det (nI)
    ! and the reference determinant
    ! Input: nI - current determinant
    ! Output: map - list of orbitals where the n-th electron goes to
    implicit none
    integer(n_int), intent(in) :: ilut(0:NifTot)
    integer :: map(nel)
    integer :: i, j
    integer :: ms
    integer(n_int) :: excitedOrbs(0:NifTot)

    ! an ilut of orbitals present in ilut and not in ref
    excitedOrbs = iand(ilut,not(ilutRef(:,1)))

    do i = 1, nel
       ! occupied orbitals get mapped to themselves
       if(IsOcc(ilut,refDet(i))) then
          map(i) = refDet(i)
       else
          ! unoccupied orbitals get mapped to the next occupied orbital
          ! only look at orbitals with the right spin
          ms = (3 + G1(refDet(i))%MS)/2
          do j = ms, nBasis, 2
             ! we check the orbitals in energetical order
             ! utilize that BRR(j) has the same spin as j
             if(IsOcc(excitedOrbs,BRR(j))) then
                map(i) = BRR(j)
                ! remove the selected orb from the candidates
                clr_orb(excitedOrbs,BRR(j))
                exit
             endif
          end do
       endif
    end do
    
  end function create_elec_map

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
            if(iEl.ne.jEl) then
               i = refDet(iEl)
               j = refDet(jEl)
               do a = 1, nBasis
                  if(.not.any(a.eq.(/i,j/))) then
                     do b = 1, nBasis
                        if(.not.any(b.eq.(/a,i,j/))) then
                           call set_ex(ex,i,j,a,b)
                           w(iEl) = w(iEl) + abs(sltcnd_excit(refDet,2,ex,tPar))
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
      real(dp) :: w(nBasis,0:symmax-1)
      integer :: j,b,iSym
      logical :: tPar
      integer :: aerr        

      ! there is one table for each symmetry and each starting orbital
      allocate(double_hole_two_sampler(nBasis,0:symmax-1), stat = aerr)
      do j = 1, nBasis
         w = 0.0_dp
         do b = 1, nBasis
            ! only same-spin and symmetry-allowed excitations from j -> b
            if(b.ne.j.and.same_spin(b,j)) &
                 w(b,G1(b)%Sym%s) = pp_weight_function(j,b)
         end do

         do iSym = 0, symmax-1
            call double_hole_two_sampler(j,iSym)%setupSampler(w(:,iSym))
         end do
      end do

    end subroutine setup_hole_two_sampler

    !------------------------------------------------------------------------------------------!

    pure subroutine set_ex(ex,i,j,a,b)
      implicit none
      integer, intent(out) :: ex(2,2)
      integer, intent(in) :: i,j,a,b

      ex(1,1) = i
      ex(1,2) = j
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
      integer :: iEl

      do iEl = 1, nel
         w(iEl) = 0
         i = refDet(iEl)
         do a = 1, nBasis
            w(iEl) = w(iEl) + acc_doub_matel(i,a)
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
      ! symmetry has to be preserved
      implicit none
      integer, intent(in) :: src, tgt
      real(dp) :: prob
      integer :: b, j
      integer :: ex(2,2), nI(nel)
      logical :: tPar

      prob = 0

      if(symAllowed(src,tgt)) then
         do b = 1, nBasis
            ! loop over all non-occupied orbitals
            if(.not.any(b.eq.refDet(:))) then
               do j = 1, nel
                  ! get the excited determinant D_j^b used for the matrix element
                  if(symAllowed(refDet(j),b)) then
                     call make_single(refDet(:), nI, j, b, ex, tPar)
                     ! this is a symbolic excitation, we do NOT require src to be occupied
                     ! we just use the formula for excitation matrix elements
                     ex(1,1) = src
                     ex(2,1) = tgt
                     prob = prob + abs(sltcnd_excit(nI,1,ex,tPar))
                  endif
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
    integer :: j,k
    deallocate(refDet)

    call single_elec_sampler%samplerDestructor()

    call clear_sampler_array(single_hole_sampler)
    call double_elec_one_sampler%samplerDestructor()
    call clear_sampler_array(double_elec_two_sampler)
    call clear_sampler_array(double_hole_one_sampler)
    do j = 1, size(double_hole_two_sampler,1)
       do k = 1, size(double_hole_two_sampler,2)
          call double_hole_two_sampler(j,k)%samplerDestructor()
       end do
    end do
    deallocate(double_hole_two_sampler)
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
!    if(G1(a)%MS.eq.G1(i)%MS) then
!       w = sqrt(abs(get_umat_el(gtID(i),gtID(a),gtID(a),gtID(i))))
!    else
!       w = 0.0_dp
!    endif


    ex(1,1) = i
    ex(1,2) = a
    ex(2,1) = i
    ex(2,2) = a
    w = sqrt(abs(sltcnd_excit(refDet,2,ex,.false.)))
  end function pp_weight_function

  !------------------------------------------------------------------------------------------!

  function symAllowed(a,b) result(allowed)
    implicit none
    integer, intent(in) :: a,b
    logical :: allowed
    
    allowed = same_spin(a,b) .and. (G1(a)%Sym%s.eq.G1(b)%Sym%s)
  end function symAllowed

    !------------------------------------------------------------------------------------------!


  pure subroutine intswap(a,b)
    integer, intent(inout) :: a,b
    integer :: tmp

    tmp = a
    a = b 
    b = tmp
  end subroutine intswap

end module pcpp_excitgen
