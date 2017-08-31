#include "macros.h"

module impurityModels
use SystemData, only: nBasis, nel
use bit_rep_data, only: niftot
use procedure_pointers, only: get_umat_el
use OneEInts, only: getTMatEl
use FciMCData, only: pBath, nBath, nImp, pSingles
use dSFMT_interface, only: genrand_real2_dSFMT
use constants, only: dp, n_int, eps, bits_n_int
use util_mod, only: abs_l1
use util_mod_numerical, only: binary_search_first_ge
! This module contains utility for treating impurity models, i.e. systems 
! which are seperable into a (large) noninteracting bath and a (small) interacting
! impurity

! In isBath, isBath(I)==.true. means that the orbital I is a bath orbital

implicit none 
integer, allocatable :: bathOrbitals(:), impurityOrbitals(:)

interface
   function sumFunc(i,j,k,l,nI)
     use constants, only: dp
     use SystemData, only: nel
     HElement_t(dp) :: sumFunc
     integer, intent(in) :: i,j,k,l,nI(nel)
   end function sumFunc
end interface
contains 
!------------------------------------------------------------------------------------------!

  subroutine constructBath(isBath)
    implicit none    
    ! This subroutine identifies the bath-sites for a given UMAT/TMAT
    logical, intent(out) :: isBath(nBasis)
    integer :: i,j,k,l

    isBath = .true.
    do i=1,nBasis
       ! for each orbital, check if it is in the bath, that is, check if there are 
       ! UMAT entries for that orbital
       do j=1,nBasis
          do k=1,nBasis
             do l=1,nBasis
                if(abs(get_umat_el(i,j,k,l))>eps) then
                   isBath(i)=.false.
                endif
             enddo
          enddo
       enddo
    enddo
  end subroutine constructBath

!------------------------------------------------------------------------------------------!

  subroutine generateOrbitalLists(isBath)
    implicit none    
    logical, intent(in) :: isBath(nBasis)
    integer :: i, iBath, iImp

    iImp = 0
    iBath = 0
    allocate(bathOrbitals(nBasis))
    allocate(impurityOrbitals(nBasis))
    do i = 1,nBasis
       if(isBath(i)) then
          iBath = iBath + 1
          bathOrbitals(iBath) = i
       else
          iImp = iImp + 1
          impurityOrbitals(iImp) = i
       endif
    enddo
    nBath = iBath
    nImp = iImp
  end subroutine generateOrbitalLists

!------------------------------------------------------------------------------------------!

  subroutine verifyBath(isBath)
    ! Here, we check if there are one-electron integrals between the bathsites. If there are,
    ! we move the corresponding bathsites to the impurity
    implicit none    
    logical, intent(inout) :: isBath(nBasis)
    integer :: i,j
    integer :: connectionCount(nBasis)
    integer :: maxConnection, posMax
    character(25), parameter :: this_routine = 'verifyBath'

    do
       ! We keep moving the most connected orbital to the impurity until the bath is a true bath
       connectionCount = 0
       do i=1,nBasis
          ! for each bath orbital, check if it really does not have 1-el integrals with other bathsites
          if(isBath(i)) then
             do j=1,nBasis
                if(isBath(j) .and. abs(getTMatEl(i,j)>eps)) then
                   ! here, we count how many 1-el integrals to the bath there are for i
                   connectionCount(i) = connectionCount(i) + 1
                endif
             enddo
          endif
       enddo
       if(SUM(connectionCount)==0) exit
       maxConnection = 0
       do j=1,nBasis
          if(connectionCount(j)>maxConnection) then
             posMax = j
             maxConnection = connectionCount(j)
          endif
       enddo
       if(posMax == 0) call stop_all(this_routine,"Error in checking bath for connections")
       isBath(posMax) = .false.
    enddo
  end subroutine verifyBath

!------------------------------------------------------------------------------------------!

  subroutine gen_excit_impurity_model(nI, ilut, nJ, ilutnJ, exFlag, IC, ExcitMat, &
                                         tParity, pGen, HElGen, part_type)
    implicit none    
    ! generate a random excitation for an impurity-type system, that is,
    ! connections between the bathsites are not taken into account
    integer, intent(in) :: nI(nel)
    integer, intent(in) :: exFlag
    integer(n_int), intent(in) :: iLut(0:niftot)
    integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(out) :: pgen
    !type(excit_gen_store_type), intent(inout), target :: store

    ! Not used
    integer(n_int), intent(out) :: ilutnJ(0:niftot)
    HElement_t(dp), intent(out) :: HElGen

    integer, intent(in), optional :: part_type

    real(dp) :: r
    character(*), parameter :: this_routine = 'gen_rand_excit'

    ! first, determine if a single or a double is created
    r = genrand_real2_dSFMT()
    if(r .lt. pSingles) then
       IC = 1
       ! for a single, check if it shall be within the impurity
       ! or between bath and impurity
       call generate_imp_single_excitation(nI,ilut,nJ,ilutnJ,ExcitMat,tParity,pGen)
    else
       IC = 2
       call generate_imp_double_excitation(nI,ilut,nJ,ilutnJ,ExcitMat,tParity,pGen)
    endif
    
  end subroutine gen_excit_impurity_model

!------------------------------------------------------------------------------------------!

  subroutine generate_imp_single_excitation(nI,ilut,nJ,ilutnJ,ExcitMat,tParity,pGen)
    implicit none    
    integer, intent(in) :: nI(nel)
    integer, intent(out) :: nJ(nel),ExcitMat(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(out) :: pGen
    integer(n_int), intent(in) :: ilut(0:niftot)
    integer(n_int), intent(out) :: ilutnJ(0:niftot)

    integer :: randImp, randDest    
    real(dp) :: r
    
    ! randomly pick an impurity orbital as all excitations contain at least
    ! one impurity orbital
    r = genrand_real2_dSFMT()
    randImp = impurityOrbitals(INT(r*nImp)+1)
    ! check if it is occupied
    r = genrand_real2_dSFMT()
    pGen = pSingles
    if(r .lt. pBath) then
       call hamiltonian_weighted_pick_single_bath(randImp,randDest,pGen,ilut,nI)
    else
       call hamiltonian_weighted_pick_single_imp(randImp,randDest,pGen,ilut,nI)
    endif
    pGen = pGen / (nImp)
    call assign_output_ilut(ExcitMat,ilut,ilutnJ,nJ,randImp,randDest)
  end subroutine generate_imp_single_excitation

!------------------------------------------------------------------------------------------!

  subroutine generate_imp_double_excitation(nI,ilut,nJ,ilutnJ,ExcitMat,tParity,pGen)
    implicit none    
    integer, intent(in) :: nI(nel)
    integer, intent(out) :: nJ(nel),ExcitMat(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(inout) :: pGen
    integer(n_int), intent(in) :: ilut(0:niftot)
    integer(n_int), intent(out) :: ilutnJ(0:niftot)

    integer :: i,j,k,l
    procedure(sumFunc), pointer :: p_getImpImpMatElDoubs
    
    pGen = 1.0_dp - pSingles
    ! pick the first electron from the impurity at random
    i = pick_random_occ_impurity(nI)
    ! then pick the second electron randomly from the remaining ones
    do 
       j = pick_random_occ_impurity(nI)
       if(i .ne. j) exit
    enddo
    ! we now have two electrons from the impurity
    ! so get two holes weighted with the matrix elements
    p_getImpImpMatElDoubs => getImpImpMatEl2Ind
    call cumulative_sum_wrapper(ilut,nI,i,0,nImp,.true.,p_getImpImpMatElDoubs,&
         pGen,k)
    p_getImpImpMatElDoubs => getImpImpMatElDoubs
    call cumulative_sum_wrapper(ilut,nI,i,0,nImp,.true.,p_getImpImpMatElDoubs,&
         pGen,l,j,k)
    call assign_output_ilut(ExcitMat,ilut,ilutnJ,nJ,i,k,j,l)
  end subroutine generate_imp_double_excitation

!------------------------------------------------------------------------------------------!

  subroutine hamiltonian_weighted_pick_single_imp(origin,destination,pGen,ilut,nI)
    use util_mod_numerical, only: binary_search_first_ge
    use sltcnd_mod, only: sltcnd_excit
    implicit none    
    ! pick an orbital from the impurity, for a single-excitation from the impurity to
    ! the impurity
    integer, intent(in) :: origin, nI(nel)
    integer, intent(out) :: destination
    real(dp), intent(inout) :: pGen
    integer(n_int), intent(in) :: ilut(0:niftot)

    procedure(sumFunc), pointer :: p_getImpImpMatEl

    p_getImpImpMatEl => getImpImpMatEl

    call cumulative_sum_wrapper(ilut,nI,origin,0,nImp,IsOcc(ilut,origin),p_getImpImpMatEl,&
         pGen,destination)
    pGen = pGen * (1.0_dp-pBath)
  end subroutine hamiltonian_weighted_pick_single_imp

!------------------------------------------------------------------------------------------!

  subroutine hamiltonian_weighted_pick_single_bath(origin,destination,pGen,ilut,nI)
    use util_mod_numerical, only: binary_search_first_ge
    implicit none    
    ! pick an orbital from the bath with probability weighted according to the 
    ! hamiltonian. This is kept general, so it can generate excitations both from
    ! the impurity to the bath and vice versa, depending if origin is occupied
    integer, intent(in) :: origin
    integer, intent(out) :: destination
    real(dp), intent(inout) :: pGen
    integer(n_int), intent(in) :: ilut(0:niftot)
    integer, intent(in) :: nI(nel)
   
    procedure(sumFunc), pointer :: p_getBathImpMatEl

    p_getBathImpMatEl => getBathImpMatEl

    ! Only count holes/parts depending on whether we are considering a part/hole in the
    ! impurity
    call cumulative_sum_wrapper(ilut,nI,origin,nImp,nBath,IsOcc(ilut,origin),p_getBathImpMatEl,&
         pGen,destination)
    pGen = pGen * pBath
  end subroutine hamiltonian_weighted_pick_single_bath

!------------------------------------------------------------------------------------------!

  subroutine cumulative_sum_wrapper(ilut,nI,origin,offset,nTarget,occupancy,sfunc,&
       pGen,destination,aux_a,aux_b)
    integer(n_int), intent(in) :: ilut(0:niftot)
    integer, intent(in) :: offset, nTarget, origin, nI(nel)
    logical, intent(in) :: occupancy
    integer, intent(in), optional :: aux_a, aux_b
    real(dp), intent(out) :: pGen
    integer, intent(out) :: destination
    procedure(sumFunc), pointer :: sfunc
    real(dp) :: cumSum(nTarget/2)
    integer :: k,l

    if(present(aux_a) .and. present(aux_b)) then
       k = aux_a
       l = aux_b
    else
       k = 0
       l = 0
    endif

    call create_cumulative_list(cumSum,ilut,nI,origin,offset,nTarget,occupancy,sfunc,&
         k,l)
    if(abs(cumSum(nTarget/2)) .le. eps) then 
       destination = origin
    else
       call pick_from_cumulative_sum(destination,pGen,cumSum,nTarget,offset,mod(origin,2))
    endif
  end subroutine cumulative_sum_wrapper

!------------------------------------------------------------------------------------------!

 subroutine create_cumulative_list(cumSum,ilut,nI,origin,offset,nTarget,occupancy,&
      sumFunc,aux_a,aux_b)
    implicit none    
    integer(n_int), intent(in) :: ilut(0:niftot)
    integer, intent(in) :: offset, nTarget, origin, nI(nel)
    real(dp), intent(out) :: cumSum(nTarget/2)
    integer, intent(in) :: aux_a, aux_b
    logical, intent(in) :: occupancy

    integer :: i, ms
    
    interface sf
       function sumFunc(i,j,k,l,nI)
         use constants, only: dp
         use SystemData, only: nel
         real(dp) :: sumFunc
         integer, intent(in) :: i,j,k,l,nI(nel)
       end function sumFunc
    end interface sf

    ms = mod(origin,2)
    
    if(occupancy .neqv. IsOcc(ilut,offset+2-ms)) then
       cumSum(1) = abs(sumFunc(origin,offset+2-ms,aux_a,aux_b,nI))
    else
       cumSum(1) = 0.0_dp
    endif
    do i = 2, nTarget/2
       if(occupancy .neqv. IsOcc(ilut,offset+2*i-ms)) then
          cumSum(i) = cumSum(i-1) + abs(sumFunc(origin,offset+2*i-ms,aux_a,aux_b,nI))
       else
          cumSum(i) = cumSum(i-1)
       endif
    end do
  end subroutine create_cumulative_list

!------------------------------------------------------------------------------------------!

  subroutine pick_from_cumulative_sum(destination,pGen,cumSum,nTargets,offset,ms)
    implicit none    
    integer, intent(out) :: destination
    integer, intent(in) :: nTargets, ms, offset
    real(dp), intent(inout) :: pGen
    real(dp), intent(in) :: cumSum(:)

    integer :: nOrbs
    real(dp) :: r

    nOrbs = nTargets/2
    r = genrand_real2_dSFMT() * cumSum(nOrbs)
    destination = binary_search_first_ge(cumSum,r)
    if(destination == 1) then
       pGen = pGen * cumSum(1)/cumSum(nOrbs)
    else
       pGen = pGen * (cumSum(destination) - cumSum(destination-1))/cumSum(nOrbs)
    endif
    destination = offset + 2*destination - ms
  end subroutine pick_from_cumulative_sum

!------------------------------------------------------------------------------------------!

  function pick_random_occ_impurity(nI) result(i)
    implicit none    
    integer, intent(in) :: nI(nel)
    integer :: i
    integer :: nOccImp
    real(dp) :: r

    r = genrand_real2_dSFMT()
    ! here, we need the number of occupied impurity orbitals
    nOccImp = binary_search_first_ge(nI,nImp)
    ! then, pick one of those at random
    i = nI(int(nOccImp*r)+1)
  end function pick_random_occ_impurity

!------------------------------------------------------------------------------------------!

  subroutine assign_output_ilut(ex,ilut,ilutnJ,nJ,i,j,k,l)
    use bit_reps, only: decode_bit_det
    implicit none
    integer, intent(in) :: i,j
    integer, intent(in), optional :: k,l
    integer(n_int), intent(in) :: ilut(0:niftot)
    integer, intent(out) :: nJ(nel), ex(2,2)
    integer(n_int), intent(out) :: ilutnJ(0:niftot)

    ex = 0
    ex(1,1) = i
    ex(1,2) = j
    ilutnJ = ilut
    clr_orb(ilutnJ,i)
    set_orb(ilutnJ,j)
    if(present(k) .and. present(l)) then
       clr_orb(ilutnJ,k)
       set_orb(ilutnJ,l)
       ex(2,1) = k
       ex(2,2) = l
    endif
    call decode_bit_det(nJ,ilutnJ)
  end subroutine assign_output_ilut

!------------------------------------------------------------------------------------------!

  function getBathImpMatEl(i,j,k,l,nI) result(mel)
    use constants, only: dp
    use SystemData, only: nel
    HElement_t(dp) :: mel
    integer, intent(in) :: i,j,k,l,nI(nel)
    
    mel = getTMatEl(i,j)
  end function getBathImpMatEl

!------------------------------------------------------------------------------------------!

  function getImpImpMatEl(i,j,k,l,nI) result(mel)
    use constants, only: dp
    use SystemData, only: nel
    HElement_t(dp) :: mel
    integer, intent(in) :: i,j,k,l,nI(nel)
    integer :: nOccImp,m

    mel = getTMatEl(i,j)
    nOccImp = binary_search_first_ge(nI,nImp)
    do m = 1,nOccImp
       mel = mel + get_umat_el(i,nI(k),j,nI(k)) - get_umat_el(i,nI(k),nI(k),i)
    enddo
  end function getImpImpMatEl

!------------------------------------------------------------------------------------------!

  function getImpImpMatElDoubs(i,j,k,l,nI) result(mel)
    use constants, only: dp
    use SystemData, only: nel
    HElement_t(dp) :: mel
    integer, intent(in) :: i,j,k,l,nI(nel)
    
    mel = get_umat_el(i,k,l,j) - get_umat_el(i,k,j,l)

  end function getImpImpMatElDoubs

!------------------------------------------------------------------------------------------!

  function getImpImpMatEl2Ind(i,j,k,l,nI) result(mel)
    use constants, only: dp
    use SystemData, only: nel
    HElement_t(dp) :: mel
    integer, intent(in) :: i,j,k,l,nI(nel)
    
    mel = get_umat_el(i,j,i,j) - get_umat_el(i,j,j,i)

  end function getImpImpMatEl2Ind

end module impurityModels
