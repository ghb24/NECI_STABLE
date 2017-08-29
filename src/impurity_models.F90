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
       call generate_imp_single_excitation(nI,ilut,nJ,ExcitMat,tParity,pGen)
    else
       IC = 2
       call generate_imp_double_excitation(nI,nJ,ExcitMat,tParity,pGen,ilut)
    endif
    
  end subroutine gen_excit_impurity_model

!------------------------------------------------------------------------------------------!

  subroutine generate_imp_double_excitation(nI,nJ,ExcitMat,tParity,pGen,ilut)
    implicit none    
    integer, intent(in) :: nI(nel)
    integer, intent(out) :: nJ(nel),ExcitMat(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(out) :: pGen
    integer(n_int), intent(in) :: ilut(0:niftot)

    real(dp) :: r, cumSum(nImp)
    integer :: i,j,k,l

    ! pick the first electron from the impurity at random
    i = pick_random_occ_impurity(nI)
    ! then pick the second electron randomly from the remaining ones
    do 
       j = pick_random_occ_impurity(nI)
       if(i .ne. j) exit
    enddo
    ! we now have two electrons from the impurity
    ! so get two holes weighted with the matrix elements
    if(IsNotOcc(ilut,1)) then
       cumSum(1) = sqrt(abs_l1(get_umat_el(i,1,i,1)))
    else
       cumSum(1) = 0
    endif
    do k = 2, nImp
       if(IsNotOcc(ilut,k)) then
          cumSum(k) = cumSum(k-1) + sqrt(abs_l1(get_umat_el(max(i,k),min(i,k),max(i,k),min(i,k))))
       else
          cumSum(k) = cumSum(k-1)
       endif
    enddo
  end subroutine generate_imp_double_excitation

!------------------------------------------------------------------------------------------!

  subroutine generate_imp_single_excitation(nI,ilut,nJ,ExcitMat,tParity,pGen)
    implicit none    
    integer, intent(in) :: nI(nel)
    integer, intent(out) :: nJ(nel),ExcitMat(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(out) :: pGen
    integer(n_int), intent(in) :: ilut(0:niftot)

    integer :: iBath, iEl
    integer :: randImp, randDest
    logical :: occ
    real(dp) :: r
    
    ! randomly pick an impurity orbital as all excitations contain at least
    ! one impurity orbital
    r = genrand_real2_dSFMT()
    randImp = impurityOrbitals(INT(r*nImp)+1)
    ! check if it is occupied
    r = genrand_real2_dSFMT()
    if(r .lt. pBath) then
       call hamiltonian_weighted_pick_single_bath(randImp,randDest,pGen,ilut)
    else
       call hamiltonian_weighted_pick_single_imp(randImp,randDest,pGen,ilut,nI)
    endif
    pGen = pGen / (nImp)

  end subroutine generate_imp_single_excitation

!------------------------------------------------------------------------------------------!

  subroutine hamiltonian_weighted_pick_single_imp(origin,destination,pGen,ilut,nI)
    use util_mod_numerical, only: binary_search_first_ge
    use sltcnd_mod, only: sltcnd_excit
    implicit none    
    ! pick an orbital from the impurity, for a single-excitation from the impurity to
    ! the impurity
    integer, intent(in) :: origin, nI(nel)
    integer, intent(out) :: destination
    real(dp), intent(out) :: pGen
    integer(n_int), intent(in) :: ilut(0:niftot)

    integer :: i
    real(dp) :: cumSum(nImp), r
    integer :: ex(2,2)

    ex = 0

    if(IsOcc(ilut,origin)) then
       ex(1,1) = origin
       if(IsNotOcc(ilut,1)) then
          ex(2,1) = 1
          cumSum(1) = abs(sltcnd_excit(nI,1,ex))
       else
          cumSum(1) = 0.0_dp
       endif
       do i = 2, nImp
          if(IsNotOcc(ilut,i)) then
             ex(2,1) = i
             cumSum(i) = cumSum(i-1) + abs(sltcnd_excit(nI,1,ex))
          else
             cumSum(i) = cumSum(i-1)   
          endif
       enddo
    else
       ex(2,1) = origin
       if(IsOcc(ilut,1)) then
          ex(2,1) = 1
          cumSum(1) = abs(sltcnd_excit(nI,1,ex))
       else
          cumSum(1) = 0.0_dp
       endif
       do i = 2, nImp
          if(IsOcc(ilut,i)) then
             ex(1,1) = i
             cumSum(i) = cumSum(i-1) + abs(sltcnd_excit(nI,1,ex))
          else
             cumSum(i) = cumSum(i-1)   
          endif
       enddo
    endif
    call pick_from_cumulative_sum(destination,pGen,cumSum,nImp)
    pGen = pGen * (1-pBath)
  end subroutine hamiltonian_weighted_pick_single_imp

!------------------------------------------------------------------------------------------!

  subroutine hamiltonian_weighted_pick_single_bath(origin,destination,pGen,ilut)
    use util_mod_numerical, only: binary_search_first_ge
    implicit none    
    ! pick an orbital from the bath with probability weighted according to the 
    ! hamiltonian. This is kept general, so it can generate excitations both from
    ! the impurity to the bath and vice versa, depending if origin is occupied
    integer, intent(in) :: origin
    integer, intent(out) :: destination
    real(dp), intent(out) :: pGen
    integer(n_int), intent(in) :: ilut(0:niftot)

    integer :: i
    real(dp) :: cumSum(nBath), r

    ! Only count holes/parts depending on whether we are considering a part/hole in the
    ! impurity
    if(IsOcc(ilut,origin)) then
       if(IsNotOcc(ilut,nImp+1)) then
          cumSum(1) = abs(getTMatEl(origin,nImp+1))    
       else
          cumSum(1) = 0.0_dp
       endif
       do i = nImp+2, nBasis
          if(IsNotOcc(ilut,i)) then
             cumSum(i) = cumSum(i-1) + abs(getTMatEl(origin,i))
          else
             cumSum(i) = cumSum(i-1)
          endif
       end do
    else
       if(IsOcc(ilut,1)) then
          cumSum(1) = abs(getTMatEl(origin,nImp+1))    
       else
          cumSum(1) = 0.0_dp
       endif
       do i = 2, nBath
          if(IsOcc(ilut,i+nImp)) then
             cumSum(i) = cumSum(i-1) + abs(getTMatEl(origin,i+nImp))
          else
             cumSum(i) = cumSum(i-1)
          endif
       end do
    endif
    call pick_from_cumulative_sum(destination,pGen,cumSum,nBath)
    pGen = pGen * pBath
  end subroutine hamiltonian_weighted_pick_single_bath

!------------------------------------------------------------------------------------------!

 subroutine create_cumulative_list(cumSum,ilut,offset,nTarget,occupancy,sumFunc)
    implicit none    
    integer(n_int), intent(in) :: ilut
    integer, intent(in) :: offset, nTarget
    real(dp), intent(out) :: cumSum
    logical, intent(in) :: occupancy
    
    interface sf
       function sumFunc(i,j)
         use constants, only: dp
         real(dp) :: sumFunc
         integer, intent(in) :: i,j
       end function sumFunc
    end interface sf
  end subroutine create_cumulative_list

!------------------------------------------------------------------------------------------!

  function pick_random_occ_impurity(nI) result(i)
    implicit none    
    integer, intent(in) :: nI(nel)
    integer :: i
    integer :: nOccImp
    real(dp) :: r

    r = genrand_real2_dSFMT()
    ! here, we need the number of occupied impurity orbitals
    nOccImp = binary_search_first_ge(nI,nImp+1)
    ! then, pick one of those at random
    i = nI(int(nOccImp*r)+1)
  end function pick_random_occ_impurity

!------------------------------------------------------------------------------------------!

  subroutine pick_from_cumulative_sum(destination,pGen,cumSum,nTargets)
    implicit none    
    integer, intent(out) :: destination
    integer, intent(in) :: nTargets
    real(dp), intent(out) :: pGen
    real(dp), intent(in) :: cumSum(:)

    real(dp) :: r

    r = genrand_real2_dSFMT() * cumSum(nTargets)
    destination = binary_search_first_ge(cumSum,r)
    if(destination == 1) then
       pGen = cumSum(1)/cumSum(nTargets)
    else
       pGen = (cumSum(destination) - cumSum(destination-1))/cumSum(nTargets)
    endif
  end subroutine pick_from_cumulative_sum

end module impurityModels
