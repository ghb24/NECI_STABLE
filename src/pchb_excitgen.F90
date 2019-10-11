#include "macros.h"
module pchb_excitgen
  use constants
  use SystemData, only: nel, nBasis, G1, tStoreSpinOrbs, AB_elec_pairs, par_elec_pairs
  use bit_rep_data, only: NIfTot
  use dSFMT_interface, only: genrand_real2_dSFMT
  use get_excit, only: make_double, exciteIlut
  use excit_gens_int_weighted, only: pick_biased_elecs
  use FciMCData, only: pSingles, excit_gen_store_type, nInvalidExcits, nValidExcits, &
       projEDet, pParallel, pDoubles
  use sltcnd_mod, only: sltcnd_excit
  use UMatCache, only: gtID
  use aliasSampling, only: aliasSamplerArray_t
  use util_mod, only: fuseIndex, linearIndex
  use GenRandSymExcitNUMod, only: construct_class_counts, createSingleExcit, &
       calc_pgen_symrandexcit2
  use SymExcitDataMod, only: pDoubNew, scratchSize
  implicit none

  type(aliasSamplerArray_t) :: pchb_sampler
  integer, allocatable :: tgtOrbs(:,:)

  contains

    ! this is the interface routine: for singles, use the uniform excitgen
    ! for doubles, the precomputed heat-bath weights
    subroutine gen_rand_excit_pchb(nI, ilutI, nJ, ilutJ, exFlag, ic, ex, tpar, &
         pgen, helgen, store, part_type)
      implicit none
    ! The interface is common to all excitation generators, see proc_ptrs.F90
      integer, intent(in) :: nI(nel), exFlag
      integer(n_int), intent(in) :: ilutI(0:NIfTot)
      integer, intent(out) :: nJ(nel), ic, ex(2,2)
      integer(n_int), intent(out) :: ilutJ(0:NIfTot)
      logical, intent(out) :: tpar
      real(dp), intent(out) :: pGen
      HElement_t(dp), intent(out) :: HElGen
      type(excit_gen_store_type), intent(inout), target :: store
      integer, intent(in), optional :: part_type

      helgen = 0.0_dp

      if(genrand_real2_dSFMT() < pSingles) then
         ic = 1
         ! default to uniform singles
         if(.not.store%tFilled) then
            call construct_class_counts(nI, store%ClassCountOcc, &
                 store%ClassCountUnocc)
            store%tFilled = .true.
         endif
         pDoubNew = 1.0-pSingles
         call createSingleExcit(nI,nJ,store%ClassCountOcc,store%classCountUnocc,ilutI,&
              ex,tPar,pGen)

         ! set the output ilut
         ilutJ = ilutI
         clr_orb(ilutJ, ex(1,1))
         set_orb(ilutJ, ex(2,1))
      else
         ic = 2
         ! use precomputed weights to generate doubles
         call generate_double_pchb(nI,ilutI,nJ,ilutJ,ex,tpar,pgen)
         pgen = pgen * (1.0-pSingles)

         if(IsNullDet(nJ)) then
            nInvalidExcits = nInvalidExcits + 1
         else
            nValidExcits = nValidExcits + 1
         endif
      end if

    end subroutine gen_rand_excit_pchb

    !------------------------------------------------------------------------------------------!

    subroutine generate_double_pchb(nI,ilutI,nJ,ilutJ,ex,tpar,pgen)
      ! given the initial determinant (both as nI and ilut), create a random double
      ! excitation using the hamiltonian matrix elements as weights
      ! Input: nI - determinant to excite from
      !        elec_map - map to translate electron picks to orbitals
      !        ilut - determinant to excite from in ilut format
      !        nJ - on return, excited determinant
      !        excitMat - on return, excitation matrix nI -> nJ
      !        tParity - on return, the parity of the excitation nI -> nJ
      !        pGen - on return, the probability of generating the excitation nI -> nJ

      integer, intent(in) :: nI(nel)
      integer(n_int), intent(in) :: ilutI(0:NIfTot)
      integer, intent(out) :: nJ(nel)
      integer(n_int), intent(out) :: ilutJ(0:NIfTot)
      integer, intent(out) :: ex(2,2)
      real(dp), intent(out) :: pgen
      logical, intent(out) :: tpar

      integer :: elecs(2), src(2), sym_prod, ispn, sum_ml, ij
      integer :: orbs(2), srcID(2), ab
      real(dp) :: pGenHoles
      logical :: invalid

      ! first, pick two random elecs
      call pick_biased_elecs(nI,elecs,src,sym_prod,ispn,sum_ml,pgen)

      ! convert to spatial orbitals if required
      srcID = gtID(src)

      invalid = .false.
      ! use the sampler for this electron pair -> order of src electrons does not matter
      ij = fuseIndex(src(1),src(2))
      ! get a pair of orbitals using the precomputed weights
      call pchb_sampler%aSample(ij,ab,pGenHoles)
      ! split the index ab (using a table containing mapping ab -> (a,b))
      orbs = tgtOrbs(:,ab)

      ! check if the picked orbs are a valid choice - if they are the same, match one
      ! occupied orbital or are zero (maybe because there are no allowed picks for
      ! the given source) abort
      invalid = (any(orbs==0) .or. any(orbs(1) == nI) &
           .or. any(orbs(2) == nI)) .or. (orbs(1) == orbs(2))

      pGen = pGen * pGenHoles
      if(invalid) then
         ! if 0 is returned, there are no excitations for the chosen elecs
         ! -> return nulldet
         nJ = 0
         ilutJ = 0_n_int
         ex(2,:) = orbs
         ex(1,:) = src
      else
         ! else, construct the det from the chosen orbs/elecs

         call make_double(nI,nJ,elecs(1),elecs(2),orbs(1),orbs(2),ex,tpar)

         ilutJ = exciteIlut(ilutI,src,orbs)
      endif

    end subroutine generate_double_pchb

  !------------------------------------------------------------------------------------------!

    function calc_pgen_pchb(nI, ex, ic, ClassCount2, ClassCountUnocc2) result(pgen)
      implicit none
      integer, intent(in) :: nI(nel)
      integer, intent(in) :: ex(2,2), ic
      integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)

      real(dp) :: pgen

      if(ic==1) then
         ! single excitations are the job of the uniform excitgen
         call calc_pgen_symrandexcit2(nI,ex,ic,ClassCount2, ClassCountUnocc2, pDoubles, pGen)
      else if(ic==2) then
         pGen = pDoubles * calc_double_pgen_pchb(ex)
      else
         pgen = 0.0_dp
      endif

    end function calc_pgen_pchb

  !------------------------------------------------------------------------------------------!

    function calc_double_pgen_pchb(ex) result(pgen)
      implicit none
      integer, intent(in) :: ex(2,2)
      real(dp) :: pgen
      integer :: ab, ij

      ! the probability of picking the two electrons: they are chosen uniformly
      if (is_alpha(ex(1,1)) .eqv. is_alpha(ex(1,2))) then
         pgen = pParallel / par_elec_pairs
      else
         pgen = (1.0_dp - pParallel) / AB_elec_pairs
      end if

      ! look up the probability for this excitation in the sampler
      ij = fuseIndex(ex(1,1),ex(1,2))
      ab = fuseIndex(ex(2,1),ex(2,2))
      pgen = pgen * pchb_sampler%aGetProb(ij,ab)

    end function calc_double_pgen_pchb

  !------------------------------------------------------------------------------------------!

    subroutine init_pchb_excitgen()
      ! initialize the pchb excitation generator
      ! this does two things:
      ! 1. setup the lookup table for the mapping ab -> (a,b)
      ! 2. setup the alias table for picking ab given ij with probability ~<ij|H|ab>
      implicit none
      integer :: ab, a, b, abMax
      integer :: aerr
      integer(int64) :: memCost

      write(iout,*) "Allocating PCHB excitation generator objects"
      ! total memory cost
      memCost = 0_int64
      ! initialize the mapping ab -> (a,b)
      abMax = fuseIndex(nBasis,nBasis)
      allocate(tgtOrbs(2,0:abMax), stat = aerr)
      do a = 1, nBasis
         do b = 1, a-1
            ab = fuseIndex(a,b)
            tgtOrbs(1,ab) = b
            tgtOrbs(2,ab) = a
         end do
      end do

      ! enable catching exceptions
      tgtOrbs(:,0) = 0

      ! setup the alias table
      call setup_pchb_sampler()

      write(iout,*) "Finished excitation generator initialization"
      write(iout,*) "Excitation generator requires", real(memCost,dp)/2.0_dp**30, "GB of memory"
      ! this is some bias used internally by CreateSingleExcit - not used here
      pDoubNew = 0.0
    contains

      subroutine setup_pchb_sampler()
        implicit none
        integer :: i,j
        integer :: ij, ijMax
        integer :: ex(2,2)
        real(dp), allocatable :: w(:)
        ! number of possible source orbital pairs
        ijMax = fuseIndex(nBasis,nBasis)
        call pchb_sampler%setupSamplerArray(int(ijMax,int64),int(abMax,int64))
        memCost = memCost + abMax*ijMax*24
        write(iout,*) "Generating samplers for PCHB excitation generator"
        ! weights per pair
        allocate(w(abMax), stat = aerr)
        do i = 1, nBasis
           ex(1,1) = i
           ! as we order a,b, we can assume j < i (j==i is not possible)
           do j = 1, i-1
              w = 0.0_dp
              ex(1,2) = j
              ! for each (i,j), get all matrix elements <ij|H|ab> and use them as
              ! weights to prepare the sampler
              do a = 1, nBasis
                 ex(2,2) = a
                 do b = 1, a-1
                    ab = fuseIndex(a,b)
                    ! ex(2,:) is in ascending order
                    ex(2,1) = b
                    ! use the actual matrix elements as weights
                    w(ab) = abs(sltcnd_excit(projEDet(:,1),2,ex,.false.))
                 end do
              end do
              ij = fuseIndex(i,j)
              call pchb_sampler%setupEntry(ij,w)
           end do
        end do

        deallocate(w)
      end subroutine setup_pchb_sampler
    end subroutine init_pchb_excitgen

  !------------------------------------------------------------------------------------------!

    subroutine finalize_pchb_excitgen()
      ! deallocate the sampler and the mapping ab -> (a,b)
      implicit none

      call pchb_sampler%samplerArrayDestructor()
      deallocate(tgtOrbs)

    end subroutine finalize_pchb_excitgen

  end module pchb_excitgen
