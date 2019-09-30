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
  use shared_memory_mpi
  use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
  use UMatCache, only: gtID, numBasisIndices
  use aliasSampling, only: aliasSamplerArray_t
  use util_mod, only: fuseIndex, linearIndex
  use GenRandSymExcitNUMod, only: construct_class_counts, createSingleExcit, &
       calc_pgen_symrandexcit2
  use SymExcitDataMod, only: pDoubNew, scratchSize
  use ParallelHelper, only: iProcIndex_intra
  implicit none

  type(aliasSamplerArray_t) :: pchb_sampler
  integer(TagIntType) :: tagTgtOrbs, tagPCHBSampler, tagAllowed
  ! the mappings sampling index -> (a,b) and ab -> sampling index
  integer(int32), pointer :: allowedOrbs(:,:), tgtOrbs(:,:,:)
  integer(MPIArg) :: win_tgt, win_allowed

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
      ! catch exception of no possible excitation
      if(ab == 0) then
         orbs = 0
      else
         orbs = tgtOrbs(:,ij,ab)
      endif

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
      integer :: ab, ij, tgt

      ! the probability of picking the two electrons: they are chosen uniformly
      if (is_alpha(ex(1,1)) .eqv. is_alpha(ex(1,2))) then
         pgen = pParallel / par_elec_pairs
      else
         pgen = (1.0_dp - pParallel) / AB_elec_pairs
      end if

      ! look up the probability for this excitation in the sampler
      ij = fuseIndex(ex(1,1),ex(1,2))
      ab = fuseIndex(ex(2,1),ex(2,2))
      ! the sampling index
      tgt = allowedOrbs(ij,ab)
      if(tgt == 0) then
         ! sampling index 0 means ab is not in the sampled set, hence pgen is 0
         pgen = 0.0_dp
      else
         pgen = pgen * pchb_sampler%aGetProb(ij,tgt)
      endif

    end function calc_double_pgen_pchb

  !------------------------------------------------------------------------------------------!

    subroutine init_pchb_excitgen()
      ! initialize the pchb excitation generator
      ! this does two things:
      ! 1. setup the lookup table for the mapping ab -> (a,b)
      ! 2. setup the alias table for picking ab given ij with probability ~<ij|H|ab>
      implicit none
      integer :: ab, a, b
      integer :: abMax, ijMax
      integer :: nBI
      integer :: aerr
      integer(int64) :: memCost
      character(*), parameter :: t_r = "init_pchb_excitgen"

      ! total memory cost
      memCost = 0_int64
      write(iout,*) "Initializing pchb excitation generator"
      ! number of spatial orbitals
      nBI= numBasisIndices(nBasis)
      ! number of possible source/target orbital pairs
      ijMax = fuseIndex(nBasis,nBasis)
      abMax = fuseIndex(nBasis,nBasis)
      ! setup the alias table
      call setup_pchb_sampler()

      ! this is some bias used internally by CreateSingleExcit - not used here
      pDoubNew = 0.0

      ! tell about the memory cost
      write(iout,*) "Finished PCHB excitation generator initialization"
      write(iout,*) "PCHB excitation generator requires", memCost/2.0_dp**30, "GB of memory"
    contains

      subroutine setup_pchb_sampler()
        ! This first determines the required memory for each pair of electrons,
        ! and lists the possible target orbitals. Then, the samplers are allocated
        ! accordingly, and seeded with these weights. Only nonzero weights are considered,
        ! this both reduces memory cost and fixes an extremely rare bug where a zero-weight excitation
        ! can be generated due to rounding errors.
        implicit none
        integer :: i,j
        integer :: ij
        integer :: ex(2,2)
        real(dp), allocatable :: w(:)
        integer(int64), allocatable :: abAllowed(:)
        integer(int64) :: tgtOrbsCost, allowedOrbsCost

        ! initialize the mapping ab -> (a,b)
        call shared_allocate_mpi(win_tgt, tgtOrbs, (/2_int64,int(ijMax,int64),int(abMax,int64)/))
        tgtOrbsCost = 2*ijMax*abMax
        call LogMemAlloc("tgtOrbs",int(tgtOrbsCost),4,t_r,tagTgtOrbs)
        memCost = memCost + tgtOrbsCost*4

        ! and the inverse mapping: (a,b) -> ab (contiguously sampled index)
        call shared_allocate_mpi(win_allowed, allowedOrbs, (/int(ijMax,int64),int(abMax,int64)/))
        allowedOrbsCost = ijMax*abMax
        call LogMemAlloc("abAllowed",int(allowedOrbsCost),4,t_r,tagAllowed)
        memCost = memCost + allowedOrbsCost*4

        ! enable catching exceptions
        tgtOrbs = 0
        allowedOrbs = 0

        write(iout,*) "Determining memory requirements"

        if(iProcIndex_intra == 0) then
           ! the number of possible (i.e. allowed) excitations per ij
           allocate(abAllowed(ijMax))
           abAllowed= 0_int64
           ! probe the required size of
           do i = 1, nBasis
              ex(1,1) = i
              ! as we order a,b, we can assume j < i (j==i is not possible)
              do j = 1, i-1
                 ex(1,2) = j
                 ! for each (i,j), get all matrix elements <ij|H|ab> and use them as
                 ! weights to prepare the sampler
                 ij = fuseIndex(i,j)
                 do a = 1, nBasis
                    ex(2,2) = a
                    do b = 1, a-1
                       ab = fuseIndex(a,b)
                       ! ex(2,:) is in ascending order
                       ex(2,1) = b
                       ! use the actual matrix elements as weights - only store nonzero weights
                       if(get_weight(ex) > eps) then
                          abAllowed(ij) = abAllowed(ij) + 1
                          ! memorize the target orbitals
                          tgtOrbs(1,ij,abAllowed(ij)) = b
                          tgtOrbs(2,ij,abAllowed(ij)) = a
                          ! memorize the index of ab in the sampling range
                          allowedOrbs(ij,ab) = abAllowed(ij)
                       endif
                    end do
                 end do
              end do
           end do

           ! allocate the sampler and set the internal pointers
           call pchb_sampler%setupSamplerArray(int(ijMax,int64),abAllowed)
           ! Log the allocation
           call LogMemAlloc("pchb sampler",int(sum(abAllowed)),8,t_r,tagPCHBSampler)
           ! Inform about memory cost
           memCost = memCost + sum(abAllowed)*24.0_dp

           write(iout,*) "Setting up sampler"

           ! now, create the probability tables
           do i = 1, nBasis
              ! the loop over electrons works the same as in probing
              ex(1,1) = i
              do j = 1, i-1
                 ex(1,2) = j
                 ij = fuseIndex(i,j)
                 ! weights per pair
                 allocate(w(abAllowed(ij)), stat = aerr)

                 ! loop over all pairs for which pgen is nonzero
                 do ab = 1, abAllowed(ij)
                    ex(2,:) = tgtOrbs(:,ij,ab)
                    ! and get the (nonzero) weight
                    w(ab) = get_weight(ex)
                 end do

                 call pchb_sampler%setupEntry(ij,w)
                 deallocate(w)
              end do
           end do

           deallocate(abAllowed)
        end if

      end subroutine setup_pchb_sampler

      function get_weight(ex) result(weight)
        ! Defines how excitations are weighted
        ! Input: ex - excitation matrix
        ! Output weight - the unnormalized probability of drawing this excitation with given electrons
        implicit none
        integer, intent(in) :: ex(2,2)
        real(dp) :: weight
        ! weight the excitations with the matrix elements
        weight = abs(sltcnd_excit(projEDet(:,1),2,ex,.false.))
      end function get_weight
    end subroutine init_pchb_excitgen

  !------------------------------------------------------------------------------------------!

    subroutine finalize_pchb_excitgen()
      ! deallocate the sampler and the mapping ab -> (a,b)
      implicit none
      character(*), parameter :: t_r = "finalize_pchb_sampler"

      call pchb_sampler%samplerArrayDestructor()
      call LogMemDealloc(t_r, tagPCHBSampler)

      if(associated(tgtOrbs)) call shared_deallocate_mpi(win_tgt, tgtOrbs)
      call LogMemDealloc(t_r, tagTgtOrbs)

      if(associated(allowedOrbs)) call shared_deallocate_mpi(win_allowed, allowedOrbs)
      call LogMemDealloc(t_r, tagAllowed)

    end subroutine finalize_pchb_excitgen

  end module pchb_excitgen
