#include "macros.h"

module impurity_models
use SystemData, only: nBasis, nel, G1
use bit_rep_data, only: NIfTot, NIfD
use procedure_pointers, only: get_umat_el
use OneEInts, only: getTMatEl
use FciMCData, only: pSingles
use dSFMT_interface, only: genrand_real2_dSFMT
use constants, only: dp, n_int, eps, bits_n_int, HEl_zero, maxExcit
use util_mod, only: abs_l1
use util_mod_numerical, only: binary_search_first_ge
use get_excit, only: make_single, make_double
use UMatCache, only: GtID
use procedure_pointers, only: get_umat_el
! This module contains utility for treating impurity models, i.e. systems 
! which are seperable into a (large) noninteracting bath and a (small) interacting
! impurity

! In isImp, isImp(I)==.true. means that the orbital I is an impurity orbital

! DISCLAIMER: 
! 1) We assume spin conservation for now. This changes for spin-orbit coupling
! 2) The impurity orbitals have to be the first ones in the orbital list
! TODO: ADD A SUBROUTINE TO INITIALIZATION THAT ORDERS THE ORBITALS, SO
! THE ORDERING DOES NOT HAVE TO BE PROVIDED

implicit none

private

public :: setupImpurityExcitgen, clearImpurityExcitgen, gen_excit_impurity_model

integer, allocatable :: connections(:,:), ImpuritySites(:), nConnects(:)
integer :: nImp

interface
   function sumFunc(i,j,k,l,ilut)
     use constants, only: dp, n_int
     use bit_rep_data, only: NIfTot
     HElement_t(dp) :: sumFunc
     integer, intent(in) :: i,j,k,l
     integer(n_int), intent(in) :: ilut(0:NIfTot)
   end function sumFunc
end interface
contains 

!------------------------------------------------------------------------------------------!

subroutine setupImpurityExcitgen()
  implicit none
  
  ! get the number of ImpuritySites
  call constructBath()
  ! first, get the number of connections for each orbital
  allocate(nConnects(nBasis))
  ! get the number of connected orbs per orb
  call constructConnections()
  ! and then the connections 
  allocate(connections(nBasis,maxval(nConnects)))
  ! get the connections
  call constructConnections()

  ! set pSingles
  call assignPSingles()

end subroutine setupImpurityExcitgen

!------------------------------------------------------------------------------------------!

subroutine constructConnections()
  implicit none
  integer :: i,j,k,l
  logical :: tconnected
  integer(n_int) :: dummy(0:NIfTot)

  dummy = 0.0_dp
  nConnects = 0
  do i = 1, nBasis
     ! for each orbital, check which ones have a nonzero one-election matrix element
     ! this defines the bath structure
     do j = 1, nBasis
        ! use a flag to store if i and j are connected
        tconnected = .false.
        if(abs(gettmatel(i,j)) > eps .and. i .ne. j) then
           tconnected = .true.
        else
           ! if the 1-electron integral is 0, check the 2-electron one
           if(i <= nImp .and. j <= nImp .and. i .ne. j) then
           ! this can only be the case for impurity orbitals
              do k = 1, nImp
                 ! the orbitals in question cannot count here, as the 
                 ! occupation of k cannot change
                 if(k .ne. i .and. k .ne. j) then
                    ! check the 2-electron integral
                    if(abs(get_umat_el(gtID(i),gtID(k),gtID(k),gtID(j))) > eps) then
                       tconnected = .true.
                       exit
                    endif
                 endif
              enddo
           endif
        endif
           if(tconnected) then
              ! increase the number of connections of i by one
              nConnects(i) = nConnects(i) + 1
              ! if allocated, add the connection
              if(allocated(connections)) &
                   connections(i,nConnects(i)) = j
           endif
     enddo
  enddo

end subroutine constructConnections

!------------------------------------------------------------------------------------------!

  subroutine constructBath()
    implicit none    
    ! This subroutine identifies the bath-sites for a given UMAT/TMAT
    logical :: isImp(nBasis)
    integer :: i,j,k,l
    character(*), parameter :: t_r = "constructBath"

    nImp = 0
    isImp = .false.
    do i=1,nBasis
        ! for each orbital, check if it is in the bath, that is, check if there are 
        ! UMAT entries for that orbital
        do j=i+1,nBasis
            do k=1,nBasis
                do l=k+1,nBasis
                    if(abs(get_umat_el(gtId(i),gtID(j),gtID(k),gtID(l)))>eps) then
                        ! as we know the first orbitals have to be the ImpuritySites,
                        ! the only information needed is their number
                        ! auxiliary flag for checking the correct ordering
                        isImp(i) = .true.
                        exit
                    endif
                    ! No double counting
                    if(isImp(i)) exit
                enddo
                if(isImp(i)) exit
            enddo
            if(isImp(i)) exit
        enddo
    enddo
    nImp = count(isImp)

    ! dummy array of impurity orbitals
    allocate(ImpuritySites(nImp))
    do i = 1, nImp
       ImpuritySites(i) = i
    enddo
    
    ! do a quick check if the orbital ordering is correct
    do i = nImp+1, nBasis
       if(isImp(i)) call stop_all(t_r, &
            "Incorrect orbital ordering. The first orbitals have to be the ImpuritySites")
    enddo
  end subroutine constructBath

!------------------------------------------------------------------------------------------!

  subroutine assignPSingles()
    use constants, only: iout
    ! set the initial values for pSingles/pDoubles
    implicit none

    pSingles = 0.9_dp
    write (iout,'(A,F14.6)') " pDoubles set to: ", 1.0_dp - pSingles
    write (iout,'(A,F14.6)') " pSingles set to: ", pSingles

  end subroutine assignPSingles

!------------------------------------------------------------------------------------------!

  ! optional, only use in the star-geometry
 
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
                if(isBath(j) .and. abs(getTMatEl(i,j))>eps) then
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
                                         tParity, pGen, HElGen, store, part_type)
    use FciMCData, only: excit_gen_store_type
    implicit none    
    ! generate a random excitation for an impurity-type system, that is,
    ! there is significant part of the system which is non-interacting (called bath)
    ! for the one-electron integrals, we have the connections cached
    integer, intent(in) :: nI(nel)
    integer, intent(in) :: exFlag
    integer(n_int), intent(in) :: iLut(0:NIfTot)
    integer, intent(out) :: nJ(nel), IC, ExcitMat(2,maxExcit)
    logical, intent(out) :: tParity
    type(excit_gen_store_type), intent(inout), target :: store
    real(dp), intent(out) :: pgen
    !type(excit_gen_store_type), intent(inout), target :: store

    integer(n_int), intent(out) :: ilutnJ(0:NIfTot)
    HElement_t(dp), intent(out) :: HElGen

    integer, intent(in), optional :: part_type

    real(dp) :: r
    integer :: nImpEls
    character(*), parameter :: this_routine = 'gen_rand_excit'

    unused_var(store)
    unused_var(part_type)
    unused_var(exFlag)

    HElGen = HEl_zero
    pSingles = 1.0_dp

    ! first, determine if a single or a double is created
    ! if there are enoguh holes and electrons in the impurity, pick randomly
    nImpEls = binary_search_first_ge(nI,nImp+1)-1
    if(nImpEls > 1 .and. (nImp - nImpEls) > 1) then
       r = genrand_real2_dSFMT()
       if(r .lt. pSingles) then
          IC = 1
          pGen = pSingles
          call generate_imp_single_excitation(nI,ilut,nJ,ilutnJ,ExcitMat,tParity,pGen)
       else
          IC = 2
          pGen = 1.0_dp - pSingles
          call generate_imp_double_excitation(nI,ilut,nJ,ilutnJ,ExcitMat,tParity,pGen)
       endif
    else
       ! else, we target a single
       IC = 1
       pGen = 1.0_dp
       call generate_imp_single_excitation(nI,ilut,nJ,ilutnJ,ExcitMat,tParity,pGen)
    endif
    
  end subroutine gen_excit_impurity_model

!------------------------------------------------------------------------------------------!

  subroutine generate_imp_single_excitation(nI,ilut,nJ,ilutnJ,ex,tParity,pGen)
    implicit none    
    integer, intent(in) :: nI(nel)
    integer, intent(out) :: nJ(nel),ex(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(inout) :: pGen
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    integer(n_int), intent(out) :: ilutnJ(0:NIfTot)
    integer, allocatable :: destPool(:)
    integer :: randSource, randDest, iElec
    real(dp) :: r

    ! first, randomly pick an orb
    randSource = pick_source_el_single_excit(nI,ilut,iElec,pGen,destPool)
    ! Randomly select from the coupled unoccupied orbitals
    r = genrand_real2_dSFMT()
    randDest = destPool(int(size(destPool)*r)+1)
   ! assign the output
    ! if the excitation is valid, generate the output
    call make_single(nI,nJ,iElec,randDest,ex,tParity)
    call assign_output_ilut(ilut,ilutnJ,randSource,randDest)
  end subroutine generate_imp_single_excitation

!------------------------------------------------------------------------------------------!

  subroutine generate_imp_double_excitation(nI,ilut,nJ,ilutnJ,ex,tParity,pGen)
      implicit none    
      integer, intent(in) :: nI(nel)
      integer, intent(out) :: nJ(nel),ex(2,2)
      logical, intent(out) :: tParity
      real(dp), intent(inout) :: pGen
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      integer(n_int), intent(out) :: ilutnJ(0:NIfTot)

      integer :: i,j,k,l,nOccImp
      procedure(sumFunc), pointer :: p_getImpImpMatElDoubs

      ! pick the first electron from the impurity at random
      i = pick_random_occ_impurity(nI,nOccImp)
      ! If an excitation is not possible, abort
      if(nOccImp == 1 .or. (nImp - nOccImp) < 2) then
          call invalidate()
          return
      endif
      ! then pick the second electron randomly from the remaining ones
      j = pick_random_occ_impurity(nI,nOccImp)
      if(i .ne. j) then
          call invalidate()
          return
      endif
      ! pgen is 1.0/nOccImp per orbital, but they can be picked in either order
      pgen = pgen*2.0/(nOccImp**2)

      ! we now have two electrons from the impurity
      ! so get two holes weighted with the matrix elements
      k = pick_random_unocc_impurity(nI,G1(i)%Ms,pgen)
      j = pick_random_unocc_impurity(nI,G1(j)%Ms,pgen)
      ! Since the spin is fixed, the holes cannot be picked in any order
      ! and write the resulting det to output
      call make_double(nI,nJ,i,j,k,l,ex,tparity)
      call assign_output_ilut(ilut,ilutnJ,i,k,j,l)

  contains

      subroutine invalidate()
          nJ = 0
          ilutnJ = 0
          pgen = 1.0
      end subroutine invalidate
  end subroutine generate_imp_double_excitation

  !------------------------------------------------------------------------------------------!
  
  function pick_source_el_single_excit(nI,ilut,iElec,pGen,pool) result(source)
      implicit none
      integer, intent(in) :: nI(nel)
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      integer, intent(out) :: iElec
      real(dp), intent(inout) :: pGen
      integer, allocatable, intent(out) :: pool(:)
      integer :: source

      logical :: tAvail(nel)
      integer :: i,j,nAvail
      integer :: nIPick(nel), poolsize(nel), pools(nel,nbasis)
      real(dp) :: r

      ! check for each electron, if there is some space to excite to
      tAvail = .false.
      nAvail = 0
      pools = 0
      poolsize = 0
      do i = 1, nel
          do j = 1, nConnects(nI(i))
              if(IsNotOcc(ilut,connections(nI(i),j))) then
                  tAvail(i) = .true.
                  poolsize(i) = poolsize(i) + 1
                  pools(i,poolsize(i)) = connections(nI(i),j)
              end if
          end do
      end do
      nAvail = count(poolsize/=0)

      nIPick = pack(nI,tAvail,vector=nI)
      r = genrand_real2_dSFMT()
      source = nIPick(int(r*nAvail)+1)
      ! each pickable electron has the same probability to be chosen
      pGen = pGen / nAvail
      ! get the electron index
      iElec = binary_search_first_ge(nI,source)
      ! return the coupled, unoccupied orbitals
      allocate(pool(1:poolsize(iElec)))
      pool = pools(iElec,1:poolsize(iElec))
    
  end function pick_source_el_single_excit
!------------------------------------------------------------------------------------------!

  function pick_random_occ_impurity(nI,nOccImp) result(i)
      implicit none    
      integer, intent(in) :: nI(nel)
      integer :: i
      integer, intent(out) :: nOccImp
      real(dp) :: r

      r = genrand_real2_dSFMT()
      ! get the number of occupied imp orbitals in nI
      nOccImp = binary_search_first_ge(nI,nImp+1)-1
      ! then, pick one of the occupied impurity orbitals at random
      i = nI(int(nOccImp*r)+1)
  end function pick_random_occ_impurity

!------------------------------------------------------------------------------------------!

  function pick_random_unocc_impurity(nI,ms,pgen) result(i)
      implicit none
      integer, intent(in) :: nI(nel)
      integer, intent(in) :: ms
      real(dp), intent(inout) :: pgen
      integer :: i

      real(dp) :: r
      integer :: pool(nImp)
      integer :: j, nEmpty
      logical :: tEmpty(nImp)

      ! First, get the unoccupied impurity sites
      do j = 1, nImp
          tEmpty(j) = (all(nI/=ImpuritySites(j)) .and. G1(ImpuritySites(j))%Ms==ms)
      end do
      pool = pack(ImpuritySites,tEmpty,vector=ImpuritySites)
      nEmpty = count(tEmpty)

      r = genrand_real2_dSFMT()
      i = pool(int(r*nEmpty)+1)

      ! Adjust pgen
      pgen = pgen * 1.0 / nEmpty
      
  end function pick_random_unocc_impurity

!------------------------------------------------------------------------------------------!

  subroutine assign_output_ilut(ilut,ilutnJ,i,j,k,l)
    implicit none
    integer, intent(in) :: i,j
    integer, intent(in), optional :: k,l
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    integer(n_int), intent(out) :: ilutnJ(0:NIfTot)

    ilutnJ = ilut
    clr_orb(ilutnJ,i)
    set_orb(ilutnJ,j)
    if(present(k) .and. present(l)) then
       clr_orb(ilutnJ,k)
       set_orb(ilutnJ,l)
    endif
  end subroutine assign_output_ilut

!------------------------------------------------------------------------------------------!

  ! at the end of the run, clear the impurity excitation generator
  subroutine clearImpurityExcitGen()
    implicit none
    
    ! deallocate the auxiliary resources used by the impurity excitation generator
    deallocate(ImpuritySites)
    deallocate(connections)
    deallocate(nConnects)
  end subroutine clearImpurityExcitGen

end module impurity_models
