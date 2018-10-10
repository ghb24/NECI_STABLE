module cachedExcitgen
  use constants
  use SystemData, only: nel, nBasis, G1
  use bit_rep_data, only: NIfTot
  use dSFMT_interface, only: genrand_real2_dSFMT
  use UMatHash
  use get_excit, only: make_double, exciteIlut
  use excit_gens_int_weighted, only: gen_single_4ind_ex, pick_biased_elecs
  use FciMCData, only: pSingles, excit_gen_store_type
  use UMatCache, only: gtID
  use DetBitOps, only: DetBitEq, EncodeBitDet
  implicit none

  contains 

    subroutine gen_excit_hel_cached(nI, ilutI, nJ, ilutJ, exFlag, ic, ex, tpar, &
         pgen, helgen, store, part_type)
      implicit none
      
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
         call gen_single_4ind_ex(nI,ilutI,nJ,ilutJ,ex,tpar,pgen)
         pgen = pgen * pSingles
      else
         ic = 2
         call generate_double_cached(nI,ilutI,nJ,ilutJ,ex,pgen,tpar)
         pgen = pgen * (1-pSingles)
      end if

    end subroutine gen_excit_hel_cached

    !------------------------------------------------------------------------------------------!

    subroutine generate_double_cached(nI,ilutI,nJ,ilutJ,ex,pgen,tpar)
      integer, intent(in) :: nI(nel)
      integer(n_int), intent(in) :: ilutI(0:NIfTot)
      integer, intent(out) :: nJ(nel)
      integer(n_int), intent(out) :: ilutJ(0:NIfTot)
      integer, intent(out) :: ex(2,2)
      real(dp), intent(out) :: pgen
      logical, intent(out) :: tpar

      integer :: elecs(2), src(2), sym_prod, ispn, sum_ml
      integer :: orbs(2),ms(2), srcID(2)
      integer(n_int) :: ilutK(0:NIfTot)
      logical :: invalid

      ! first, pick two random elecs
      call pick_biased_elecs(nI,elecs,src,sym_prod,ispn,sum_ml,pgen)
      ! TODO bias towards large/small
      
      ! convert to spatial orbitals if required
      srcID = gtID(src)

      invalid = .false.
      ! check if we have a spin-parallel or spin-opposite excitation
      if(ispn == 2) then
         ! use the spin-opposite cache
         call selectRSFromCSUM(srcID(1),srcID(2),orbs,pgen,CumSparseUMat,&
              nPQ, rsPQ, posPQ)

      else
         call selectRSFromCSUM(srcID(1),srcID(2),orbs,pgen,CumSparseUMatPar,&
              nPQPar, rsPQPar, posPQPar)  
      endif

      ! convert the spatial orbitals picked from the cached CSUM to spin orbs
      orbs(1) = getSpinOrb(orbs(1),G1(src(1))%ms)
      orbs(2) = getSpinOrb(orbs(2),G1(src(2))%ms)

      invalid = (any(orbs==0) .or. any(orbs(1) == nI) &
           .or. any(orbs(2) == nI)) .or. (orbs(1) == orbs(2))

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

      contains 

        subroutine checkSpinVsOcc(ind)
          integer, intent(in) :: ind

          if(any(orbs(ind)==srcID)) then
             if(any(nI == 2*orbs(ind))) then
                if(any(nI == 2*orbs(ind) - 1) ) then
                   invalid = .true.
                else
                   ms(ind) = -1               
                endif
             else
                ms(ind) = 1
             endif
          endif

        end subroutine checkSpinVsOcc

    end subroutine generate_double_cached

    function getSpinOrb(spatialOrb, ms) result(spinOrb)
      implicit none
      integer, intent(in) :: spatialOrb, ms
      integer :: spinOrb

      spinOrb = 2*spatialOrb + (ms-1)/2
      
    end function getSpinOrb

end module cachedExcitgen
