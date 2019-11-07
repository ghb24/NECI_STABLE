#include "macros.h"

module impurityModels
use SystemData, only: nBasis, nel
use bit_rep_data, only: NIfTot, NIfD
use procedure_pointers, only: get_umat_el
use OneEInts, only: getTMatEl
use FciMCData, only: pSingles
use dSFMT_interface, only: genrand_real2_dSFMT
use constants, only: dp, n_int, eps, bits_n_int, HEl_zero
use util_mod, only: abs_l1
use util_mod_numerical, only: binary_search_first_ge
use get_excit, only: make_single, make_double
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
                    if(abs(getImpImpMatElDoubs(i,k,k,j,dummy)) > eps) then
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
       do j=1,nBasis
          do k=1,nBasis
             do l=1,nBasis
                if(abs(get_umat_el(i,j,k,l))>eps) then
                   ! as we know the first orbitals have to be the ImpuritySites,
                   ! the only information needed is their number
                   nImp = nImp + 1
                   ! auxiliary flag for checking the correct ordering
                   isImp(i) = .true.
                endif
             enddo
          enddo
       enddo
    enddo

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

    pSingles = 1.0_dp
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
    integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
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

    unused_variable(exFlag)
    unused_variable(part_type)
    unused_variable(store)

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

    integer :: randSource, randDest, iElec
    real(dp) :: r

    !r = genrand_real2_dSFMT()
    !iElec = INT(r*nel) + 1
    !randSource = nI(iElec)
    ! first, randomly pick an orb
    randSource = pick_source_el_single_excit(nI,ilut,pGen)
    ! get the electron index
    iElec = binary_search_first_ge(nI,randSource)
    call hamiltonian_weighted_pick_single(randSource,randDest,pGen,ilut)
    ! assign the output
    ! note that an invalid excitation (i.e. no possible excitations with the
    ! chosen electron) is passed as randDest=0
    if(randDest == 0) then
       nJ(1) = 0
       ilutnJ = ilut
       return
    endif
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
    ! then pick the second electron randomly from the remaining ones
    do
       j = pick_random_occ_impurity(nI,nOccImp)
       if(i .ne. j .or. nOccImp .eq. 1) exit
    enddo
    ! we now have two electrons from the impurity
    ! so get two holes weighted with the matrix elements
    p_getImpImpMatElDoubs => getImpImpMatEl2Ind
    call cumulative_sum_wrapper(ilut,i,nImp,ImpuritySites,p_getImpImpMatElDoubs,&
         pGen,k)
    ! only pick a second hole if k is a valid hole
    if(k .ne. 0) then
       ! pick a second hole
       p_getImpImpMatElDoubs => getImpImpMatElDoubs
       call cumulative_sum_wrapper(ilut,j,nImp,ImpuritySites,p_getImpImpMatElDoubs,&
            pGen,l,i,k)
    endif

    ! and write the resulting det to output
    call make_double(nI,nJ,i,j,k,l,ex,tparity)
    call assign_output_ilut(ilut,ilutnJ,i,k,j,l)
  end subroutine generate_imp_double_excitation

!------------------------------------------------------------------------------------------!

  subroutine hamiltonian_weighted_pick_single(origin,destination,pGen,ilut)
    use util_mod_numerical, only: binary_search_first_ge
    use sltcnd_mod, only: sltcnd_excit
    implicit none
    ! pick an orbital from the impurity, for a single-excitation from the impurity to
    ! the impurity
    integer, intent(in) :: origin
    integer, intent(out) :: destination
    real(dp), intent(inout) :: pGen
    integer(n_int), intent(in) :: ilut(0:NIfTot)

    procedure(sumFunc), pointer :: p_getImpImpMatEl

    ! we use the cumulative sum framework, with the tmat-matrix elements
    p_getImpImpMatEl => getImpImpMatEl

    call cumulative_sum_wrapper(ilut,origin,nConnects(origin),connections(origin,:),p_getImpImpMatEl,&
         pGen,destination)
  end subroutine hamiltonian_weighted_pick_single

!------------------------------------------------------------------------------------------!

  subroutine cumulative_sum_wrapper(ilut,origin,nTarget,targetOrbs,sfunc,&
       pGen,destination,aux_a,aux_b)
    ! This wrapper function takes a given range of orbitals and create a cumulative
    ! list of the matrix elements with origin, using some function sfunc
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    integer, intent(in) :: nTarget, origin
    integer, intent(in) :: targetOrbs(nTarget)
    integer, intent(in), optional :: aux_a, aux_b
    real(dp), intent(out) :: pGen
    integer, intent(out) :: destination
    procedure(sumFunc), pointer :: sfunc
    real(dp) :: cumSum(nTarget)
    integer :: k,l

    if(present(aux_a) .and. present(aux_b)) then
       k = aux_a
       l = aux_b
    else
       k = 0
       l = 0
    endif


    ! first, generate the cumulative sum
    cumSum = create_cumulative_list(ilut,origin,nTarget,targetOrbs,sfunc,&
         k,l)
    ! for an empty sum, the excitation is invalid
    if(abs(cumSum(nTarget)) .le. eps) then
       destination = 0
       pgen = 0.0_dp
       ! pGen does not change, no probabilistics involved
    else
        pgen = 1.0_dp
        call pick_from_cumulative_sum(destination,pGen,cumSum,nTarget,targetOrbs)
    endif
  end subroutine cumulative_sum_wrapper

!------------------------------------------------------------------------------------------!

 function create_cumulative_list(ilut,origin,nTarget,targetOrbs, &
      sumFunc,aux_a,aux_b) result(cumSum)
   ! This creates a list of cumulatively added function values of an input function
    implicit none
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    ! origin is the starting point and fixed argument of the input function
    ! offset is the first orbital to start with
    integer, intent(in) :: nTarget, origin
    integer, intent(in) :: targetOrbs(nTarget)
    real(dp) :: cumSum(nTarget)
    ! aux_a and aux_b are additional arguments for the input function, which might
    ! not matter
    integer, intent(in) :: aux_a, aux_b

    integer :: i, ms

! The function of which to create a cumulative list
    interface sf
       function sumFunc(i,j,k,l,ilut)
         use constants, only: dp, n_int
         use bit_rep_data, only: NIfTot
         HElement_t(dp) :: sumFunc
         integer, intent(in) :: i,j,k,l
         integer(n_int), intent(in) :: ilut(0:NIfTot)
       end function sumFunc
    end interface sf

    ! Only consider orbitals with the correct spin
    ! we assume s_z conservation here, else, we have to drop this
    ms = mod(origin,2)

    if(.not. IsOcc(ilut,targetOrbs(1))) then
       cumSum(1) = abs(sumFunc(origin,targetOrbs(1),aux_a,aux_b,ilut))
    else
       cumSum(1) = 0.0_dp
    endif
    do i = 2, nTarget
       if(.not. IsOcc(ilut,targetOrbs(i))) then
          cumSum(i) = cumSum(i-1) + abs(sumFunc(origin,targetOrbs(i),aux_a,aux_b,ilut))
       else
          cumSum(i) = cumSum(i-1)
       endif
    end do
  end function create_cumulative_list

!------------------------------------------------------------------------------------------!

  subroutine pick_from_cumulative_sum(destination,pGen,cumSum,nTargets,targetOrbs)
    implicit none
    integer, intent(out) :: destination
    integer, intent(in) :: nTargets
    integer, intent(in) :: targetOrbs(nTargets)
    real(dp), intent(inout) :: pGen
    real(dp), intent(in) :: cumSum(:)

    real(dp) :: r

    ! randomly pick one entry of cumSum
    r = genrand_real2_dSFMT() * cumSum(nTargets)
    ! with probability given by the difference to the previous entry (entry 0 is said to be 0)
    destination = binary_search_first_ge(cumSum,r)
    ! get pGen from the entries
    if(destination == 1) then
       pGen = pGen * cumSum(1)/cumSum(nTargets)
    else
       pGen = pGen * (cumSum(destination) - cumSum(destination-1))/cumSum(nTargets)
    endif
    ! map the entry to an orbital
    destination = targetOrbs(destination)
  end subroutine pick_from_cumulative_sum

!------------------------------------------------------------------------------------------!

  function pick_source_el_single_excit(nI,ilut,pGen) result(source)
    implicit none
    integer, intent(in) :: nI(nel)
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    real(dp), intent(inout) :: pGen
    integer :: source

    logical :: tAvail(nel)
    integer :: i,j,nAvail
    integer :: nIPick(nel)
    real(dp) :: r

    ! check for each electron, if there is some space to excite to
    tAvail = .false.
    nAvail = 0
    do i = 1, nel
       do j = 1, nConnects(i)
          if(IsNotOcc(ilut,connections(i,j))) then
             tAvail(i) = .true.
             nAvail = nAvail + 1
             exit
          end if
       end do
    end do


    nIPick = pack(nI,tAvail)
    r = genrand_real2_dSFMT()
    source = nIPick(int(r*nAvail)+1)
    ! each pickable electron has the same probability to be chosen
    pGen = pGen / nAvail

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

  function getBathImpMatEl(i,j,k,l,ilut) result(mel)
    use constants, only: dp, n_int
    use bit_rep_data, only: NIfTot
    HElement_t(dp) :: mel
    integer, intent(in) :: i,j,k,l
    integer(n_int), intent(in) :: ilut(0:NIfTot)

    unused_variable(k)
    unused_variable(l)
    unused_variable(ilut)

    mel = getTMatEl(i,j)
  end function getBathImpMatEl

!------------------------------------------------------------------------------------------!

  function getImpImpMatEl(i,j,k,l,ilut) result(mel)
    use constants, only: dp, n_int
    use bit_rep_data, only: NIfTot
    HElement_t(dp) :: mel
    integer, intent(in) :: i,j,k,l
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    integer :: m

    unused_variable(k)
    unused_variable(l)

    mel = getTMatEl(i,j)
    if(i <= nImp .and. j <= nImp) then
       do m = 1, nImp
          if(IsOcc(ilut,ImpuritySites(m))) &
          mel = mel + getImpImpMatElDoubs(i,ImpuritySites(m),j,ImpuritySites(m),ilut)
       enddo
    end if
  end function getImpImpMatEl

!------------------------------------------------------------------------------------------!

  function getImpImpMatElDoubs(i,j,k,l,ilut) result(mel)
    use constants, only: dp, n_int
    use bit_rep_data, only: NIfTot
    HElement_t(dp) :: mel
    integer, intent(in) :: i,j,k,l
    integer(n_int), intent(in) :: ilut(0:NIfTot)

    unused_variable(ilut)

    mel = get_umat_el(i,k,l,j) - get_umat_el(i,k,j,l)

  end function getImpImpMatElDoubs

!------------------------------------------------------------------------------------------!

  function getImpImpMatEl2Ind(i,j,k,l,ilut) result(mel)
    use constants, only: dp, n_int
    use bit_rep_data, only: NIfTot
    HElement_t(dp) :: mel
    integer, intent(in) :: i,j,k,l
    integer(n_int), intent(in) :: ilut(0:NIfTot)

    unused_variable(k)
    unused_variable(l)
    unused_variable(ilut)

    mel = get_umat_el(i,j,i,j) - get_umat_el(i,j,j,i)

  end function getImpImpMatEl2Ind

!------------------------------------------------------------------------------------------!

  ! at the end of the run, clear the impurity excitation generator
  subroutine clearImpurityExcitGen()
    implicit none

    ! deallocate the auxiliary resources used by the impurity excitation generator
    deallocate(ImpuritySites)
    deallocate(connections)
    deallocate(nConnects)
  end subroutine clearImpurityExcitGen

end module impurityModels
