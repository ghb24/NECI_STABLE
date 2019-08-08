#include "macros.h"
module tc_three_body_excitgen
  use constants
  use SystemData, only: nel, nOccAlpha, nOccBeta, nBasis, G1, t_exclude_3_body_excits
  use lattice_mod, only: sort_unique
  use bit_rep_data, only: NIfTot
  use k_space_hubbard, only: make_triple
  use tc_three_body_data
  use FciMCData, only: excit_gen_store_type, pDoubles, pSingles
  use dSFMT_interface, only: genrand_real2_dSFMT
  use lattice_models_utils, only: make_ilutJ
  use util_mod, only: choose
  use excit_gens_int_weighted, only: pick_biased_elecs
  use GenRandSymExcitNUMod, only: calc_pgen_symrandexcit2, ScratchSize, &
       createSingleExcit, createDoubExcit, construct_class_counts
  use SymExcitDataMod, only: pDoubNew
  use procedure_pointers, only: generate_two_body_excitation
  contains

    subroutine gen_excit_mol_tc(nI, ilut, nJ, ilutJ, exFlag, ic, ExcitMat, &
         tParity, pGen, HelGen, store, part_type)
      implicit none
      integer, intent(in) :: nI(nel), exFlag
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      integer, intent(out) :: nJ(nel), IC, ExcitMat(2,maxExcit)
      logical, intent(out) :: tParity
      real(dp), intent(out) :: pGen
      HElement_t(dp), intent(out) :: HElGen
      type(excit_gen_store_type), intent(inout), target :: store
      integer, intent(in), optional :: part_type
      integer(n_int), intent(out) :: ilutJ(0:NIfTot)
      character(*), parameter :: this_routine = 'gen_excit_mol_tc'

      real(dp) :: r

      r = genrand_real2_dSFMT()
      ! select if a triple shall be generate
      if(r < pTriples) then
         call generate_triple_excit(nI, ilut, nJ, ilutJ, ExcitMat, tParity, pGen, &
              HelGen, store)
         pGen = pGen * pTriples
         IC = 3
      else
         call generate_two_body_excitation(nI, ilut, nJ, ilutJ, exFlag, ic, ExcitMat, &
              tParity, pGen, HelGen, store, part_type)
         pGen = pGen * (1.0 - pTriples)
      endif
    end subroutine gen_excit_mol_tc

!------------------------------------------------------------------------------------------!

    function calc_pgen_mol_tc(nI, ilutI, ex, ic, ClassCount, &
         ClassCountUnocc, pDoub) result(pgen)
      implicit none
      integer, intent(in) :: nI(nel), ex(2,ic), ic
      integer, intent(in) :: ClassCount(ScratchSize)
      integer, intent(in) :: ClassCountUnocc(ScratchSize)
      integer(n_int), intent(in) :: ilutI(0:NIfTot)
      real(dp), intent(in) :: pDoub
      real(dp) :: pgen

      character(*), parameter :: t_r = "calc_pgen_mol_tc"

      pgen = 0.0_dp
      if(ic < 3) then
         ! for singles/doubles, delegate the calculation to the existing routines
         call calc_pgen_symrandexcit2(nI, ex, ic, ClassCount, ClassCountUnocc, pDoub, pgen)
         ! and take into account the bias for triples
      else if(ic == 3) then
         ! else, use the local routine
         pgen = calc_pgen_triple(nI, ex)
      endif
      
    end function calc_pgen_mol_tc

!------------------------------------------------------------------------------------------!

    function calc_pgen_triple(nI, ex) result(pgen)
      ! get the probability to get excitation `ex` from a determinant `nI`
      implicit none
      integer, intent(in) :: nI(nel), ex(2,3)
      real(dp) :: pgen
      integer :: ms, i
      character(*), parameter :: t_r = "calc_pgen_triple"

      ! get the spin 
      ms = 0
      ! sum up the spin of the single orbitals
      do i = 1, 3
         ms = ms + G1(ex(2,i))%ms
      end do
      
      ! start with pTriples
      pgen = pTriples
      ! then, add the spin bias
      if(ms == -3) then
         pgen = pgen * p0A * pgen3B
      else if(ms == -1) then
         pgen = pgen * p2B * pgen2B
      else if(ms == 1) then
         pgen = pgen * pgen1B * (1 - p2B - p0A - p0B)
      else if(ms == 3) then
         pgen = pgen * p0B * pgen0B
      else
         call stop_all(t_r, "Invalid spin")
      endif

    end function calc_pgen_triple

!------------------------------------------------------------------------------------------!

    subroutine setup_mol_tc_excitgen(HF)
      implicit none
      integer, intent(in) :: HF(nel)

      ! initialize the biases and auxiliary variables for the molecular
      ! transcorrelated 3-body excitation generator

      call init_mol_tc_biases(HF)
      call precompute_pgen()
    end subroutine setup_mol_tc_excitgen

!------------------------------------------------------------------------------------------!

    subroutine precompute_pgen()
      implicit none
      
      ! set the number of unoccupied alpha/beta
      nUnoccAlpha = nBasis/2 - nOccAlpha
      nUnoccBeta = nBasis/2 - nOccBeta

      ! the number of valid triple excitations is just given by the binomial coefficients
      pgen3B = nOccBeta * (nOccBeta - 1) * (nOccBeta - 2) * nUnoccBeta * (nUnoccBeta - 1) &
           * (nUnoccBeta - 2)
      pgen3B = scaleInvert(36.0_dp, pgen3B)

      pgen2B = nOccBeta * (nOccBeta - 1) * nOccAlpha * nUnoccBeta * (nUnoccBeta - 1) * nUnoccAlpha
      pgen2B = scaleInvert(4.0_dp, pgen2B)

      pgen1B = nOccBeta * nOccAlpha * (nOccAlpha - 1) * nUnoccBeta * nUnoccAlpha * (nUnoccAlpha - 1)
      pgen1B = scaleInvert(4.0_dp, pgen1B)

      pgen0B = nOccAlpha * (nOccAlpha - 1) * (nOccAlpha - 2) * nUnoccAlpha * (nUnoccAlpha - 1) &
           * (nUnoccAlpha - 2)
      pgen0B = scaleInvert(36.0_dp, pgen0B)

      contains
        pure function scaleInvert(scl, p) result(sp)
          real(dp), intent(in) :: scl, p
          real(dp) :: sp

          if(p > eps) then
             sp = scl / p
          else
             sp = 0.0_dp
          end if
        end function scaleInvert
    end subroutine precompute_pgen

!------------------------------------------------------------------------------------------!

    subroutine init_mol_tc_biases(HF)
      implicit none
      ! reference determinant for initializing the biases
      integer, intent(in) :: HF(nel)
      real(dp) :: normalization

      ! if we read in a value, use that one
      if(abs(pTriples) < eps) then
         pTriples = 0.1
      endif
      ! for 2 electrons, there are obviously no 
      ! triple excitations
      if(nel.eq.2) t_exclude_3_body_excits = .true.
      if(t_exclude_3_body_excits) then
         pTriples = 0.0_dp
         return
      endif
      write(iout,*) "pTriples set to ", pTriples
      ! pSingles and pDoubles add to 1, and pTriples is an additional bias not to do
      ! a two-body excitation
      
      ! scale the probabilities with the number of possible picks
      normalization = choose(nel,3)
      p0A = choose(nOccBeta,3)/normalization
      p0B = choose(noccAlpha,3)/normalization
      p2B = choose(nOccBeta,2)*nOccAlpha/normalization
      p1B = 1.0_dp - p0A - p0B - p2B
    end subroutine init_mol_tc_biases

!------------------------------------------------------------------------------------------!

    subroutine generate_triple_excit(nI, ilutI, nJ, ilutJ, ExcitMat, tParity, pGen, &
         HelGen, store)
      implicit none
      integer, intent(in) :: nI(nel)
      integer(n_int), intent(in) :: ilutI(0:NIfTot)
      integer, intent(out) :: nJ(nel), ExcitMat(2,maxExcit)
      logical, intent(out) :: tParity
      real(dp), intent(out) :: pGen
      HElement_t(dp), intent(out) :: HElGen
      type(excit_gen_store_type), intent(inout), target :: store
      integer(n_int), intent(out) :: ilutJ(0:NIfTot)

      integer :: sym_prod, src(3), tgt(3), elecs(3)
      integer :: ms

      HElGen = 0.0_dp

      ! first, pick three electrons at random
      call pick_three_elecs(nI, elecs, src, sym_prod, pgen, ms)

      ! if three electrons can be picked
      if(src(3) .ne. 0) then
         ! get three unoccupied orbitals with the same ms
         call pick_three_orbs(nI, tgt, pgen, ms)
         
         ! and create a triple excitation
         call make_triple(nI, nJ, elecs, tgt, ExcitMat, tParity)
      
         ilutJ = make_ilutJ(ilutI, ExcitMat, 3)
      else
         ! else, the excitation is invalid
         nJ = 0
         ilutJ = 0_n_int
         ExcitMat = 0
         tParity = 0
      endif
      
    end subroutine generate_triple_excit

!------------------------------------------------------------------------------------------!

    subroutine pick_three_elecs(nI, elecs, src, sym_prod, pgen, ms)
      ! picks three random electrons from nI, biased towards 0, 1, 2 or 3 beta electrons
      implicit none
      integer, intent(in) :: nI(nel)
      integer, intent(out) :: elecs(3), src(3), sym_prod, ms
      real(dp), intent(out) :: pgen
      integer :: ispn, sum_ml
      integer :: nOcc
      real(dp) :: r
      logical :: pickAlpha

      call pick_biased_elecs(nI, elecs(1:2), src(1:2), sym_prod, ispn, sum_ml, &
           pgen,p0B+p0A,p0B)
      if(ispn .eq. 3 .or. ispn .eq. 1) then
         ! all elecs have the same spin
         pickAlpha = ispn .eq. 3
         if(pickAlpha) then
            nOcc = nOccAlpha
            ms = 3
         else
            nOcc = nOccBeta
            ms = -3
         endif
         call get_missing_elec(nI,elecs,nOcc,2,pickAlpha,pgen)
         pgen = pgen * 3.0_dp
      else
         ! first picked one alpha and the beta
         r = genrand_real2_dSFMT()
         if(r < p2B) then
            ! we have 2 beta elecs + 1 alpha
            ! then pick the second beta
            call get_missing_elec(nI,elecs,nOccBeta,1,.false.,pgen)
            ms = -1
            pgen = pgen * p2B/(1.0_dp - (p0B+p0A))
         else
            ! we have 2 alpha elecs + 1 beta
            ! then pick the second alpha
            call get_missing_elec(nI,elecs,nOccAlpha,1,.true.,pgen)
            ms = 1
            pgen = pgen * p1B/(1.0_dp - (p0B+p0A))
         endif
      endif

      ! sort the generated electrons
      elecs = sort_unique(elecs)
      ! check for invalid excitation
      ! after sorting, invalid electrons are definitly at position 1
      if(elecs(1) == 0) then
         src = 0
      else
         src = nI(elecs)
      end if

    end subroutine pick_three_elecs

!------------------------------------------------------------------------------------------!

    subroutine get_missing_elec(nI, elecs, nOcc, nPicked, tAlpha, pgen)
      implicit none
      ! after picking two electrons using the symrandexcit routines,
      ! get a third one with the missing spin
      integer, intent(in) :: nI(nel)
      integer, intent(inout) :: elecs(3) ! picked electrons, elecs(1:2) has to be assigned on entry
      integer, intent(in) :: nOcc    ! number of occupied orbs of the type of the missing one
      integer, intent(in) :: nPicked ! number of already picked elecs
      logical, intent(in) :: tAlpha  ! if the missing electron is alpha
      real(dp), intent(inout) :: pgen ! probability of generation
      real(dp) :: r
      integer :: index, i, count
      
      elecs(3) = 0
      if(nOcc <= nPicked) then
         ! there are not enought elecs available - invalid excitation
         return
      endif

      r = genrand_real2_dSFMT()
      ! pick a random index of an electron
      index = int((nOcc-nPicked)*r) + 1

      ! get the corresponding electron
      count = 0
      do i = 1, nel
         if(tAlpha .neqv. is_beta(nI(i))) then
            count = count + 1
            if(count == index) elecs(3) = i
         end if
      end do
      ! the picked electrons are not counted, so skip their indices
      
      ! if we need, skip the first index
      call skipElec(1)
      ! if we need, skip the second index
      call skipElec(2)
      ! note that elecs(2) > elecs(1), so we cannot end up on elecs(1) at the end of
      ! skipping

      ! uniformly chosen
      pgen = pgen / (nOcc - nPicked)

      contains 
        
        subroutine skipElec(ind)
          implicit none
          integer, intent(in) :: ind

          ! if we need to skip an index
          if( (tAlpha .neqv. is_beta(nI(elecs(ind))))) then
             ! if we are above the index, we need to add 1 more, because we did not
             ! take elecs(ind) into account when picking elecs(3)
             if(elecs(3) >= elecs(ind)) then
                ! jump to the next electron with the right spin
                elecs(3) = elecs(3) + 1
                do while(.not. (tAlpha .neqv. is_beta(nI(elecs(3))) ) )
                   elecs(3) = elecs(3) + 1
                end do
             end if
          endif
        end subroutine skipElec
    end subroutine get_missing_elec

!------------------------------------------------------------------------------------------!

    subroutine pick_three_orbs(nI, tgt, pgen, ms)
      implicit none
      ! picks three random unoccupied orbitals, given the occupied orbitals
      integer, intent(in) :: nI(nel), ms
      integer, intent(out) :: tgt(3)
      real(dp), intent(inout) :: pgen

      integer :: i, msCur, msOrb

      msCur = ms
      do i = 0, 2
         ! get the ms of this orb
         ! we take ms = 1 until the total leftover ms is negative
         if(msCur > 0) then
            msOrb = 1
         else
         ! then we start taking ms = -1
            msOrb = -1
         endif
         call get_rand_orb(nI,tgt,msOrb,i, pgen)
         ! the remaining ms
         msCur = msCur - msOrb
      end do

      ! adjust the probability by taking permutations into account
      pgen = pgen * 4 * abs(ms)
      
    end subroutine pick_three_orbs

!------------------------------------------------------------------------------------------!

    subroutine get_rand_orb(nI, tgt, ms, nPicked, pgen)
      implicit none
      integer, intent(inout) :: tgt(3)
      integer, intent(in) :: ms, nPicked, nI(nel)
      real(dp), intent(inout) :: pgen

      integer :: pool ! available orbitals
      integer :: iOrb, i
      real(dp) :: r

      ! get the number of possible orbitals
      if(ms > 0) then
         pool = nUnoccAlpha
      else
         pool = nUnoccBeta
      end if
      ! we need to see how many same spin orbs have been picked so far
      do i = 1, nPicked
         if((ms > 0) .neqv. is_beta(tgt(i))) pool = pool - 1
      end do
      
      ! pick a random index
      r = genrand_real2_dSFMT()
      iOrb = spinOrb(int(r*pool) + 1)

      ! check if the orb is already targeted
      call skipPicked()
      ! assign the orbital
      tgt(nPicked+1) = iOrb

      ! adjust the probability
      pgen = pgen / pool
      
      ! we need to sort tgt (see above)
      tgt(1:(nPicked+1) ) = sort_unique(tgt(1:(nPicked+1) ))

      contains
        pure function spinOrb(orb) result(sorb)
          implicit none
          integer, intent(in) :: orb
          integer :: sorb

          sorb = 2*orb + (ms-1)/2
        end function spinOrb

        subroutine skipPicked()
          implicit none
          integer :: i
          integer :: invalidOrbs(nel+nPicked)

          invalidOrbs(1:nel) = nI
          invalidOrbs((nel+1):(nel+nPicked)) = tgt(1:nPicked)
          
          invalidOrbs = sort_unique(invalidOrbs)

          do i = 1, (nPicked+nel)
             ! check if the orb is already targeted
             ! assumes tgt is sorted
             if(invalidOrbs(i) <= iOrb .and. G1(invalidOrbs(i))%ms == ms) iOrb = iOrb + 2
          end do
        end subroutine skipPicked
    end subroutine get_rand_orb

!------------------------------------------------------------------------------------------!

end module tc_three_body_excitgen
