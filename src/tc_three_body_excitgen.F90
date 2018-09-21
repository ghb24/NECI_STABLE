module tc_three_body_excitgen
  use constants
  use SystemData, only: nel
  use lattice_mod, only: sort_unique
  use bit_rep_data, only: NIfTot
  use k_space_hubbard, only: make_triple
  use tc_three_body_data
  contains

    subroutine gen_excit_mol_tc(nI, ilut, nJ, ilutJ, exFlag, ic, ExcitMat, &
         ic, ExcitMat, tParity, pGen, HelGen, store, part_type)
      implicit none
      integer, intent(in) :: nI(nel), exFlag
      integer(n_int), intent(in) :: ilutI(0:NIfTot)
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
         call generate_triple_excit(nI, ilutI, nJ, ilutJ, ExcitMat, tParity, pGen, &
              HelGen, store, part_type)
         pGen = pGen * pTriples
         IC = 3
      else
         ! if no triple is generated, fall back to the normal excitation generation
         ! we need a uniform excitgen as evaluating matrix elements is too expensive
         call gen_rand_excit(nI, ilutI, nJ, ilutJ, exFlag, ic, ExcitMat, tParity, &
              pGen, HelGen, store, part_type)
         pGen = pGen * (1 - pTriples)
      endif

    end subroutine gen_excit_mol_tc

!------------------------------------------------------------------------------------------!

    subroutine generate_triple_excit(nI, ilutI, nJ, ilutJ, ExcitMat, tParity, pGen, &
         HelGen, store, part_type)
      implicit none
      integer, intent(in) :: nI(nel)
      integer(n_int), intent(in) :: ilutI(0:NIfTot)
      integer, intent(out) :: nJ(nel), IC, ExcitMat(2,3)
      logical, intent(out) :: tParity
      real(dp), intent(out) :: pGen
      HElement_t(dp), intent(out) :: HElGen
      type(excit_gen_store_type), intent(inout), target :: store
      integer, intent(in), optional :: part_type
      integer(n_int), intent(out) :: ilutJ(0:NIfTot)

      integer :: sym_prod, src(3), tgt(3)
      integer :: ms

      call pick_three_elecs(nI, elecs, src, sym_prod, pgen, ms)

      call pick_three_orbs(nI, ilutI, tgt, pgen, ms)

      call make_triple(nI, nJ, elecs, tgt, ExcitMat, tParity)
      
      ilutJ = make_ilutJ(ilutI, ExcitMat, 3)
      
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

      r = genrand_real2_dSFMT()
      if(r < p0B .or. r > (p0B+p1B+p2B)) then
         ! all elecs have the same spin
         pickAlpha = r < p0B ! if we pick 3 alphas or 3 betas
         ! first, pick 2 elecs with the normal routines
         call pick_biased_elecs_tpar(nI, elecs(1:2, src(1:2), sym_prod, ispn, sum_ml, &
              pgen,.true.,pickAlpha)
         if(pickAlpha) then
            nOcc = nOccAlpha
            ms = 3
         else
            nOcc = nOccBeta
            ms = -3
         endif
         call get_missing_elec(elecs,nOcc,2,pickAlpha,pgen)        
      else
         ! first pick one alpha and the beta
         call pick_biased_elecs_tpar(nI, elecs(1:2), src(1:2), sym_prod, ispn, sum_ml, &
              pgen, .false.,.false.)
         if(r > (p0B + p1B)) then
            ! we have 2 beta elecs + 1 alpha
            ! then pick the second beta
            call get_missing_elec(elecs,nOccBeta,1,.false.,pgen)
            ms = -1
         else
            ! we have 2 alpha elecs + 1 beta
            ! then pick the second alpha
            call get_missing_elecs(elecs,nOccAlpha,1,.true.,pgen)
            ms = 1
         endif
      endif

      ! sort the generated electrons
      elecs = sort_unique(elecs)
      if(elecs(3) == 0) then
         src = 0
      else
         src(3) = nI(elecs)
      end if

    end subroutine pick_three_elecs

!------------------------------------------------------------------------------------------!

    subroutine get_missing_elec(elecs, nOcc, nPicked, tAlpha, pgen)
      implicit none
      ! after picking two electrons using the symrandexcit routines,
      ! get a third one with the missing spin
      integer, intent(inout) :: elecs(3) ! picked electrons, elecs(1:2) has to be assigned on entry
      integer, intent(in) :: nOcc    ! number of occupied orbs of the type of the missing one
      integer, intent(in) :: nPicked ! number of already picked elecs
      logical, intent(in) :: tAlpha  ! if the missing electron is alpha
      real(dp), intent(inout) :: pgen ! probability of generation
      real(dp) :: r
      
      elecs(3) = 0
      if(nOcc <= nPicked) then
         ! there are not enought elecs available - invalid excitation
         return
      endif

      r = genrand_real2_dSFMT()
      ! pick a random index of an electron
      index = int((nOcc-nPicked)*r) + 1
      ! the picked electrons are not counted, so skip their indices
      if((tAlpha .xor. is_beta(elecs(1))) .or. (nPicked .eq. 2)) &
         if(index >= elecs(1)) index = index + 1
      if((tAlpha .xor. is_alpha(elecs(1))) .or. (nPicked .eq. 2)) &
         if(index >= elecs(2)) index = index + 1
      ! note that elecs(2) > elecs(1), so we cannot end up on elecs(1) at the end of
      ! skipping

      count = 0
      do i = 1, nel
         if(tAlpha .xor. is_beta(nI(i))) then
            count = count + 1
            if(count == index) elec = i
         end if
      end do

      ! uniformly chosen
      pgen = pgen / (nOcc - nPicked)
    end subroutine get_missing_elec

!------------------------------------------------------------------------------------------!

    subroutine pick_three_orbs(nI, ilutI, tgt, pgen, ms)
      implicit none
      ! picks three random unoccupied orbitals, given the occupied orbitals
      integer, intent(in) :: nI(nel), ms
      integer(n_int), intent(in) :: ilutI(0:NIfTot)
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
         call get_rand_orb(tgt,msOrb,i, pgen)
         ! the remaining ms
         msCur = msCur - msOrb
      end do

      ! adjust the probability by taking permutations into account
      pgen = pgen * 2 * abs(ms)
      
    end subroutine pick_three_orbs

!------------------------------------------------------------------------------------------!

    subroutine get_rand_orb(tgt, ms, nPicked, pgen)
      implicit none
      integer, intent(inout) :: tgt(3)
      integer, intent(in) :: ms, nPicked
      real(dp) :: intent(inout) :: pgen

      integer :: pool ! available orbitals
      integer :: iOrb, i
      real(dp) :: iOrb

      ! get the number of possible orbitals
      if(ms > 0) then
         pool = nUnoccAlpha
      else
         pool = nUnoccBeta
      end if
      ! we need to see how many same spin orbs have been picked so far
      do i = 1, nPicked
         if((ms > 0) .xor. is_beta(tgt(i))) pool = pool - 1
      end do
      
      ! pick a random index
      r = genrand_real2_dSFMT()
      iOrb = int(r*pool) + 1

      do i = 1, nPicked
         ! check if the orb is already occuped
         ! assumes tgt is sorted
         if(tgt(i) .eq. (2*iOrb + (ms - 1)/2 ) ) iOrb = iOrb + 1
      end do
      ! assign the orbital
      tgt(nPicked+1) = iOrb

      ! adjust the probability
      pgen = pgen / pool
      
      ! we need to sort tgt (see above)
      tgt = sort_unqiue(tgt)
    end subroutine get_rand_orb

!------------------------------------------------------------------------------------------!

end module tc_three_body_excitgen
