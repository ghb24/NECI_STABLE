#include "macros.h"
module gasci
  use SystemData, only: tNConservingGAS, tSpinConservingGAS, nBasis, nel
  use constants
  use util_mod, only: get_free_unit, binary_search_first_ge
  use bit_rep_data, only: NIfTot, NIfD
  use dSFMT_interface, only: genrand_real2_dSFMT
  use FciMCData, only: pDoubles
  use get_excit, only: make_double, make_single
  use Determinants, only: get_helement
  use excit_gens_int_weighted, only: pick_biased_elecs, pgen_select_orb
  use sltcnd_mod, only: sltcnd_excit
  implicit none

  integer :: nGAS
  ! bitmasks containing the active spaces (stored in the same format as an ilut)
  ! also have spin-resolved bitmasks
  integer(n_int), allocatable :: gasOrbs(:,:), gasSpinOrbs(:,:,:)
  ! integer list of orbitals for each active space plus respective sizes
  integer, allocatable :: gasSpinOrbList(:,:,:), gasSize(:)
  ! lookup table containing the active space for each orbital
  integer, allocatable :: gasTable(:)
  
  contains

    subroutine loadGAS()
      integer :: gas_unit, iOrb, nOrbs
      integer :: iGAS(1000)
      ! integer with all odd (ms=2) bits set to 1
#ifdef __INT64
      integer(n_int), parameter :: oddBits = 6148914691236517205_n_int
#else
      integer(n_int), parameter :: oddBits = 1431655765_n_int
#endif

      gas_unit = get_free_unit()
      open(gas_unit, file="GASOrbs",status='old')
      nOrbs = nBasis/2
      read(gas_unit,*) iGAS(1:nOrbs)
      nGAS = maxval(iGAS)
      allocate(gasOrbs(0:NIfD,nGAS))
      do iOrb = 1, nOrbs
         ! now write the orbitals read in to the current GAS
         ! set both the alpha- and the beta-spin orbital
         call setOrb(gasOrbs(:,iGAS(iOrb)),2*iOrb)
         call setOrb(gasOrbs(:,iGAS(iOrb)),2*iOrb-1)
      end do

      ! now set up the auxiliary gas lookup tables
      allocate(gasSpinOrbs(0:NIfD,nGAS,2))
      ! gasSpinOrbs is the same as gasOrbs, but spin resolved
      ! convention: odd numbers have ms=1, even have ms=2
      gasSpinOrbs(:,:,1) = iand(gasOrbs,ishft(oddBits,1))
      gasSpinOrbs(:,:,2) = iand(gasOrbs,oddBits)

      allocate(gasTable(nBasis))
      ! gasTable contains the active space index for each spin orbital
      ! it is the same as iGAS for spin orbitals
      do iOrb = 1, nOrbs
         gasTable(2*iOrb) = iGAS(iOrb)
         gasTable(2*iOrb-1) = iGAS(iOrb)
      end do

      allocate(gasSpinOrbList(nBasis,nGAS,2))
      allocate(gasSize(nGAS))
      gasSize = 0
      ! an integer list version of gasSpinOrbs (instead of binary the format)
      do iOrb = 1, nOrbs
         gasSize(iGAS(iOrb)) = gasSize(iGAS(iOrb)) + 1
         gasSpinOrbList(gasSize(iGAS(iOrb)),iGAS(iOrb),1) = 2*iOrb
         gasSpinOrbList(gasSize(iGAS(iOrb)),iGAS(iOrb),2) = 2*iOrb-1
      end do

      do iOrb = 1, nGAS
         print *, "Number of orbs in GAS", iOrb, "is", sum(popCnt(gasOrbs(:,iOrb)))
      end do

      contains 

        subroutine splitLine(line,vals,n)
          implicit none
          character(*), intent(in) :: line
          integer, intent(out) :: vals(:)
          integer, intent(out) :: n

          integer :: status, buffer(1000)
          
          n = 1
          do 
             ! use a buffer to catch exceptions
             read(line,'(A)',iostat=status) buffer(1:n)
             if(status.ne.0) exit
             vals(1:n) = buffer(1:n)
             n=n+1
          end do
          ! we overcounted by 1 because we started at 1
          n = n - 1
        end subroutine splitLine

       subroutine setOrb(ilut,orb)
          implicit none
          integer(n_int), intent(inout) :: ilut(0:NIfD)
          integer, intent(in) :: orb

          integer :: pos
          ! get the index of the integer containing this orbital
          pos = (orb-1) / bits_n_int
          ilut(pos) = ibset(ilut(pos),mod(orb-1,bits_n_int))
        end subroutine setOrb
    end subroutine loadGAS

!------------------------------------------------------------------------------------------!

    subroutine clearGAS()
      implicit none
      
      deallocate(gasSize)
      deallocate(gasSpinOrbList)
      deallocate(gasTable)
      deallocate(gasSpinOrbs)
      deallocate(gasOrbs)
    end subroutine clearGAS

!------------------------------------------------------------------------------------------!
    
    function isValidExcit(ilutI,ilutJ) result(valid)
      ! check if the excitation from ilutI to ilutJ is valid within the GAS
      implicit none
      integer(n_int), intent(in) :: ilutI(0:NIfTot), ilutJ(0:NIfTot)
      integer(n_int) :: gasI(0:NIfD), gasJ(0:NIfD)
      logical :: valid

      integer :: i, x(0:NIfD)
      ! integers with all even bits set
#ifdef __INT64
      integer(n_int), parameter :: oddBits = 6148914691236517205_n_int
#else
      integer(n_int), parameter :: oddBits = 1431655765_n_int
#endif

      valid = .true.
      ! safety check: do the gasOrbs exist
      if(.not. allocated(gasOrbs)) return
      do i = 1, nGAS
         gasI = gasComponent(ilutI,i)
         gasJ = gasComponent(ilutJ,i)
         x = popCnt(gasI)
         x = popCnt(gasJ)
         ! check if ilutI and ilutJ have the same number of electrons in GAS-i
         if(sum(popCnt(gasI)) .ne. sum(popCnt(gasJ))) valid = .false.
         if(tSpinConservingGAS) then
            ! check if the number of even bits set (=number of beta electrons) is the
            ! same in both GAS
            if(sum(popCnt(iand(gasI,oddBits))) .ne. sum(popCnt(iand(gasJ,oddBits)))) &
                 valid = .false.
         endif
      end do

      contains 
        
        function gasComponent(ilut, i) result(gasIlut)
          integer(n_int), intent(in) :: ilut(0:NIfTot)
          integer, intent(in) :: i
          
          integer(n_int) :: gasIlut(0:NIfD)
          
          gasIlut = iand(ilut(0:NIfD),gasOrbs(0:NIfD,i))
          
        end function gasComponent
        
    end function isValidExcit

!------------------------------------------------------------------------------------------!

    subroutine generate_nGAS_excitation(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                          ex, tParity, pGen, hel, store, part_type)
      ! particle number conserving GAS excitation generator:
      ! we create only excitations, that preserver the number of electrons within
      ! each active space
      use SystemData, only: nel
      use FciMCData, only: excit_gen_store_type
      use constants
      implicit none

      integer, intent(in) :: nI(nel), exFlag
      integer(n_int), intent(in) :: ilutI(0:NIfTot)
      integer, intent(out) :: nJ(nel), ic, ex(2,2)
      integer(n_int), intent(out) :: ilutJ(0:NifTot)
      real(dp), intent(out) :: pGen
      logical, intent(out) :: tParity
      HElement_t(dp), intent(out) :: hel
      type(excit_gen_store_type), intent(inout), target :: store
      integer, intent(in), optional :: part_type

      real(dp) :: r
      integer :: exFlag_unused
      integer :: pt_unused

      ! parameters of the function pointer which are not needed
      exFlag_unused = exFlag
      pt_unused = part_type
      hel = 0.0_dp

      ! single or double excitation?
      r = genrand_real2_dSFMT()
      if(r < pDoubles) then
         call generate_nGAS_double(nI,ilutI,nJ,ilutJ,ex,tParity,pgen)
         pgen = pgen * pDoubles
         ic = 2
      else
         call generate_nGAS_single(nI,ilutI,nJ,ilutJ,ex,tParity,pgen)
         pgen = pgen * (1.0_dp - pDoubles)
         ic = 1
      endif

    end subroutine generate_nGAS_excitation

!------------------------------------------------------------------------------------------!

    subroutine generate_nGAS_single(nI,ilutI,nJ,ilutJ,ex,par,pgen)
      implicit none
      integer, intent(in) :: nI(nel)
      integer(n_int), intent(in) :: ilutI(0:NIfTot)
      integer, intent(out) :: nJ(nel), ex(2,2)
      integer(n_int), intent(out) :: ilutJ(0:NIfTot)
      logical, intent(out) :: par
      real(dp), intent(out) :: pgen

      integer :: elec, tgt, src, ms
      real(dp) :: r

      ! we assume that each active space has a possible single excitation
      ! (this is very mild because active spaces without such are trivial)
      ! -> pick any random electron
      r = genrand_real2_dSFMT()
      elec = int(r*nel)+1
      src = nI(elec)
      ! adjust pgen
      pgen = 1.0_dp / nel
      ! from the same active space, get a hole
      ms = get_spin(src)
      tgt = pick_weighted_hole(nI, src, 0, 0, 1, ms, gasTable(src), pgen)

      if(tgt == 0) then
         nJ(1) = 0
         ilutJ = 0_n_int
         return
      endif

      call make_single(nI, nJ, elec, tgt, ex, par)
      ilutJ = ilutI
      clr_orb (ilutJ, src)
      set_orb (ilutJ, tgt)

    end subroutine generate_nGAS_single

!------------------------------------------------------------------------------------------!

    subroutine generate_nGAS_double(nI,ilutI,nJ,ilutJ,ex,par,pgen)
      implicit none
      integer, intent(in) :: nI(nel)
      integer(n_int), intent(in) :: ilutI(0:NIfTot)
      integer, intent(out) :: nJ(nel), ex(2,2)
      integer(n_int), intent(out) :: ilutJ(0:NIfTot)
      logical, intent(out) :: par
      real(dp), intent(out) :: pgen

      integer :: elecs(2), src(2), sym_product, ispn, sum_ml, tgt(2)
      integer :: srcGAS(2)
      integer :: ms, nJBase(nel)
      real(dp) :: r, pgen_pick1, pgen_pick2
      logical :: tExchange
      ! assuming that there are possible excitations within each active space,
      ! pick two random electrons (we would not include a full/empty space, so
      ! the assumption is very mild)
      call pick_biased_elecs(nI, elecs, src, sym_product, ispn, sum_ml, pgen)
      
      ! active spaces we consider
      srcGAS = gasTable(src)

      tExchange = (ispn == 2) .and. (.not. tSpinConservingGAS .or. srcGAS(1) == srcGAS(2))

      ! pick an empty orb from each active space chosen
      ! number of empty orbs in this active space
      r = genrand_real2_dSFMT()
      ! get the spin of the target orbital
      if(tExchange) then
         if(r>0.5_dp) then 
            ms = 1
            r = 2*(r-0.5_dp)
         else
            ms = 2
            r = 2*(0.5_dp-r)
         endif
         pgen = pgen * 0.5_dp
      else
         ms = get_spin(src(1))
      end if
      tgt = 0
      ! the first hole is chosen randomly from the first active space
      tgt(1) = pick_hole_from_active_space(ilutI, nI, srcGAS(1), ms, r, pgen_pick1)

      if(tExchange) then
         ! we picked the first spin randomly, now the second elec has the opposite spin
         ms = 3 - ms
      else
         ms = get_spin(src(2))
      endif
      ! the second hole is chosen in a weighted fashion
      tgt(2) = pick_weighted_hole(nI, src(1), src(2), tgt(1), 2, ms, srcGAS(2), pgen_pick1)
      if(any(tgt==0) .or. tgt(1) == tgt(2)) then
         call zeroResult()
         return
      end if

      ! if both excits are in the same GAS, we could also have picked them the other
      ! way around
      if(srcGAS(1)==srcGAS(2)) then
         nJBase = nI
         nJBase(elecs(1)) = tgt(2)
         pgen_pick2 = get_pgen_pick_weighted_hole(nI, src(1), src(2), tgt(2), tgt(1)) * &
              get_pgen_pick_hole_from_active_space(ilutI,srcGAS(2),get_spin(tgt(2)))
         pgen = pgen *(pgen_pick1 + pgen_pick2)
      else
         pgen = pgen * pgen_pick1
      endif

      call make_double(nI,nJ,elecs(1),elecs(2),tgt(1),tgt(2),ex,par)
      ilutJ = ilutI
      clr_orb (ilutJ, src(1))
      clr_orb (ilutJ, src(2))
      set_orb (ilutJ, tgt(1))
      set_orb (ilutJ, tgt(2))

      contains
        
        subroutine zeroResult()
          implicit none
          pgen = pgen * pgen_pick1
          ex(1,:) = sort_unique(src)
          ex(2,:) = tgt
          nJ(1) = 0
          ilutJ = 0_n_int
        end subroutine zeroResult

    end subroutine generate_nGAS_double

!------------------------------------------------------------------------------------------!

    function get_pgen_pick_weighted_hole(nI, src1, src2, tgt1, tgt2) result(pgenVal)
      implicit none
      integer, intent(in) :: nI(nel)
      integer, intent(in) :: src1, src2, tgt1, tgt2
      
      real(dp) :: pgenVal
      real(dp) :: cSum(gasSize(gasTable(tgt2)))
      integer :: gasList(gasSize(gasTable(tgt2))), nOrbs
      integer :: gasInd2

      nOrbs = gasSize(gasTable(tgt2))
      gasList = gasSpinOrbList(1:nOrbs,gasTable(tgt2),get_spin(tgt2))

      cSum = get_cumulative_list(gasList, nI, src1, src2, tgt1, 2)
      ! we know gasList contains tgt2, so we can look up its index with binary search
      gasInd2 = binary_search_first_ge(gasList,tgt2)
      if(gasInd2==1) then
         pgenVal = cSum(1)/cSum(nOrbs)
      else
         pgenVal = (cSum(gasInd2) - cSum(gasInd2-1))/cSum(nOrbs)
      endif

    end function get_pgen_pick_weighted_hole

!------------------------------------------------------------------------------------------!

    function get_pgen_pick_hole_from_active_space(ilut, srcGASInd, ms) result(pgenVal)
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfD)
      integer, intent(in) :: srcGASInd, ms
      
      real(dp) :: pgenVal
      integer :: nEmpty
      nEmpty = gasSize(srcGASInd) - sum(popCnt(iand(ilut(0:NIfD),gasSpinOrbs(0:NIfD,srcGASInd,ms))))
      ! adjust pgen
      pgenVal = 1.0_dp / real(nEmpty,dp)
    end function get_pgen_pick_hole_from_active_space

!------------------------------------------------------------------------------------------!

    function get_cumulative_list(gasList, nI, src1, src2, tgt1, ic) result(cSum)
      implicit none
      integer, intent(in) :: gasList(:)
      integer, intent(in) :: nI(nel), ic, src1, src2, tgt1
      real(dp) :: cSum(size(gasList))

      integer :: ex(2,2), i, nOrbs

      nOrbs = size(gasList)
      ex(1,1) = src1
      if(ic == 2) then
         ex(1,2) = src2
         ex(2,1) = tgt1
      else
         ex(1,2) = 0
         ex(2,1) = 0
      endif
      ! build the cumulative list of matrix elements <src|H|tgt>
      call addToCumulative(1,0.0_dp)
      do i = 2, nOrbs
         call addToCumulative(i,cSum(i-1))
      end do

      contains

        subroutine addToCumulative(i, base)
          implicit none
          integer, intent(in) :: i
          real(dp), intent(in) :: base

          if(.not. ex(2,1)==gasList(i) .and. .not. any(nI==gasList(i))) then
             ex(2,ic) = gasList(i)
             cSum(i) = abs(sltcnd_excit(nI,ic,ex,.false.)) + base
          else
             cSum(i) = base
          endif
        end subroutine addToCumulative
    end function get_cumulative_list

!------------------------------------------------------------------------------------------!

    function pick_weighted_hole(nI, src1, src2, tgt1, ic, ms, srcGASInd, pgen) result(tgt)
      implicit none
      ! pick a hole of nI with spin ms from the active space with index srcGASInd
      ! the random number is to be supplied as r
      ! nI is the source determinant, nJBase the one from which we obtain the ket of
      ! the matrix element by single excitation 
      integer, intent(in) :: nI(nel)
      integer, intent(in) :: src1, src2, tgt1, ic, ms, srcGASInd
      real(dp), intent(inout) :: pgen
      integer :: tgt
      integer :: nOrbs
      real(dp) :: r
      real(dp) :: cSum(gasSize(srcGASInd))
      integer :: gasList(gasSize(srcGASInd))
      
      ! initialize auxiliary variables
      nOrbs = gasSize(srcGASInd)
      gasList = gasSpinOrbList(1:nOrbs,srcGASInd,ms)
      ! build the cumulative list of matrix elements <src|H|tgt>
      cSum = get_cumulative_list(gasList, nI, src1, src2, tgt1, ic)
      
      ! now, pick with the weight from the cumulative list
      r = genrand_real2_dSFMT() * cSum(nOrbs)

      ! there might not be such an excitation
      if(cSum(nOrbs) > 0) then
         ! find the index of the target orbital in the gasList
         tgt = binary_search_first_ge(cSum,r) 

         ! adjust pgen with the probability for picking tgt from the cumulative list
         if(tgt==1) then
            pgen = pgen * cSum(1)/cSum(nOrbs)
         else
            pgen = pgen * (cSum(tgt) - cSum(tgt-1))/cSum(nOrbs)
         endif

         ! convert to global orbital index
         tgt = gasList(tgt)
      else
         tgt = 0
      endif
      
    end function pick_weighted_hole

!------------------------------------------------------------------------------------------!

    function pick_hole_from_active_space(ilutI, nI, srcGASInd, ms, r, pgen) result(tgt)
      implicit none
      ! pick a hole of ilutI with spin ms from the active space with index srcGASInd
      ! the random number is to be supplied as r
      integer(n_int), intent(in) :: ilutI(0:NIfD)
      integer, intent(in) :: nI(nel)
      integer, intent(in) :: ms, srcGASInd
      real(dp), intent(in) :: r
      real(dp), intent(out) :: pgen
      integer :: tgt
      integer :: nEmpty, nOrb

      nEmpty = gasSize(srcGASInd) - sum(popCnt(iand(ilutI(0:NIfD),gasSpinOrbs(0:NIfD,srcGASInd,ms))))
      ! adjust pgen
      pgen = 1.0_dp / real(nEmpty,dp)
      ! TODO: check if there are enough empty orbs
      ! Re: We can safely assume there is always an empty orb in each active space

      ! index of the target orb
      nOrb = int(r*nEmpty)+1
      call skipOrb(nOrb,gasSpinOrbs(0:NIfD,srcGASInd,ms))
      nOrb = gasSpinOrbList(nOrb,srcGASInd,ms)
      tgt = nOrb

    contains 

      subroutine skipOrb(nOrb,gasIlut)
        ! convert an unoccupied active space orbital index to an
        ! active space orbital index (i.e. correct for occ. orbs)
        implicit none
        integer, intent(inout) :: nOrb
        integer(n_int), intent(in) :: gasIlut(0:NIfD)
        integer :: j, globalOrbIndex

        do j = 1, nel
           globalOrbIndex = gasSpinOrbList(nOrb,srcGASInd,ms)
           if(nI(j) > globalOrbIndex) return
           ! check if an occ. orb is in the target active space           
           if(validTarget(nI(j),gasIlut)) &
                ! if yes, we skip this orbital, increase nOrb by 1
                nOrb = nOrb + 1
        end do
      end subroutine skipOrb

      function validTarget(orb,gasIlut) result(valid)
        ! check if an orbital is a valid target for an excitation within an active space
        implicit none
        integer, intent(in) :: orb
        integer(n_int), intent(in) :: gasIlut(0:NIfD)
        
        logical :: valid

        ! is the orbital in the acitve space?
        valid = btest(gasIlut((orb-1)/bits_n_int),mod(orb-1,bits_n_int))
      end function validTarget
    end function pick_hole_from_active_space

!------------------------------------------------------------------------------------------!

    function sort_unique(list) result(output)
      ! sorts an array of unique integers, 
      ! written by Werner Dobrautz
        integer, intent(in) :: list(:)
        integer, allocatable :: output(:)

        integer :: i, min_val,  max_val, unique(size(list))

        unique = 0
        i = 0 
        min_val = minval(list) - 1
        max_val = maxval(list) 

        do while (min_val < max_val) 
            i = i + 1 
            min_val = minval(list, mask=list>min_val) 
            unique(i) = min_val 
        end do
        allocate(output(i), source=unique(1:i))

    end function sort_unique

end module gasci
