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
use get_excit, only: make_single, make_double
! This module contains utility for treating impurity models, i.e. systems 
! which are seperable into a (large) noninteracting bath and a (small) interacting
! impurity

! In isBath, isBath(I)==.true. means that the orbital I is a bath orbital

implicit none 
integer, allocatable :: connections(:,:), impuritySites(:), bathSites(:)
integer :: nBath, nImp
logical, allocatable :: isImp(:)

interface
   function sumFunc(i,j,k,l,ilut)
     use constants, only: dp
     use bit_rep_data, only: NIfTot
     HElement_t(dp) :: sumFunc
     integer, intent(in) :: i,j,k,l,ilut(NIfTot)
   end function sumFunc
end interface
contains 

!------------------------------------------------------------------------------------------!

subroutine setupImpurityExcitgen()
  include lattice_mod, only: lat, aim
  implicit none
  integer :: i,j
  
  ! set up a cluster aim lattice
  lat => aim('cluster',1,1)
  ! the lattice has information on the connectivity of the sites
  ! first, allocate the direct connection cache
  allocate(connections(nBasis,lat%get_nconnect_max()))
  do i = 1, nBasis
     ! in each entry, we store the neighbors 
     connections(i) = lat%get_neighbors(i)
  enddo
  ! store the impurity and bath orbitals contiguously
  nImp = lat%get_n_imps()
  allocate(impuritySites(nImp))
  impuritySites = lat%get_impurities()
  
  nBath = lat%get_n_bath()
  allocate(bathSites(nBath))
  bathSites = lat%get_bath()

  ! also, for each orbital, we need to be able to quickly check
  ! if it is in the impurity
  allocate(isImp(nBasis))
  do i = 1, nBasis
     isImp(i) = lat%is_impurity_site(i)
  enddo

end subroutine setup_impurity_excitgen

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
    allocate(bathSites(nBasis))
    allocate(impuritySites(nBasis))
    do i = 1,nBasis
       if(isBath(i)) then
          iBath = iBath + 1
          bathSites(iBath) = i
       else
          iImp = iImp + 1
          impuritySites(iImp) = i
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
                                         tParity, pGen, HElGen, part_type)
    implicit none    
    ! generate a random excitation for an impurity-type system, that is,
    ! there is significant part of the system which is non-interacting (called bath)
    ! for the one-electron integrals, we have the connections cached
    integer, intent(in) :: nI(nel)
    integer, intent(in) :: exFlag
    integer(n_int), intent(in) :: iLut(0:niftot)
    integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(out) :: pgen
    !type(excit_gen_store_type), intent(inout), target :: store

    integer(n_int), intent(out) :: ilutnJ(0:niftot)
    HElement_t(dp), intent(out) :: HElGen

    integer, intent(in), optional :: part_type

    real(dp) :: r
    integer :: nImpEls
    character(*), parameter :: this_routine = 'gen_rand_excit'

    ! first, determine if a single or a double is created
    ! if there are enoguh holes and electrons in the impurity, pick randomly
    nImpEls = binary_search_first_ge(nI,nImp+1)
    if(nImpEls > 1 .and. (nImp - nImpEls) > 1) then
       r = genrand_real2_dSFMT()
       if(r .lt. pSingles) then
          IC = 1
          ! for a single, check if it shall be within the impurity
          ! or between bath and impurity
          pGen = pSingles
          call generate_imp_single_excitation(nI,ilut,nJ,ilutnJ,ExcitMat,tParity,pGen)
       else
          IC = 2
          pGen = 1.0_dp - pSingles
          call generate_imp_double_excitation(nI,ilut,nJ,ilutnJ,ExcitMat,tParity,pGen)
       endif
    else
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
    integer(n_int), intent(in) :: ilut(0:niftot)
    integer(n_int), intent(out) :: ilutnJ(0:niftot)

    integer :: randSource, randDest, iElec
    real(dp) :: r

    ! first, randomly pick an electron
    r = genrand_real2_dSFMT()
    iElec = INT(r*nel)+1
    randSource = nI(iElec)
    ! check if it is occupied
    call hamiltonian_weighted_pick_single(randSource,randDest,pGen,ilut,nI)
    ! each electron has the same probability to be picked
    pGen = pGen / (nel)
    ! assign the output
    ! note that an invalid excitation (i.e. no possible excitations with the
    ! chosen electron) is passed as randDest=0
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
    integer(n_int), intent(in) :: ilut(0:niftot)
    integer(n_int), intent(out) :: ilutnJ(0:niftot)

    integer :: i,j,k,l
    procedure(sumFunc), pointer :: p_getImpImpMatElDoubs
    
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

    ! and write the resulting det to output
    call make_double(nI,nJ,i,j,k,l,ex,tparity)
    call assign_output_ilut(ilut,ilutnJ,i,k,j,l)
  end subroutine generate_imp_double_excitation

!------------------------------------------------------------------------------------------!

  subroutine hamiltonian_weighted_pick_single(origin,destination,pGen,ilut,nI)
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

    ! we use the cumulative sum framework, with the tmat-matrix elements
    p_getImpImpMatEl => getImpImpMatEl

    call cumulative_sum_wrapper(ilut,nI,origin,0,nImp,IsOcc(ilut,origin),p_getImpImpMatEl,&
         pGen,destination)
    pGen = pGen * (1.0_dp-pBath)
  end subroutine hamiltonian_weighted_pick_single_imp

!------------------------------------------------------------------------------------------!

  subroutine cumulative_sum_wrapper(ilut,nI,origin,offset,nTarget,occupancy,sfunc,&
       pGen,destination,aux_a,aux_b)
    ! This wrapper function takes a given range of orbitals and create a cumulative 
    ! list of the matrix elements with origin, using some function sfunc
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
   ! This creates a list of cumulatively added function values of an input function
    implicit none    
    integer(n_int), intent(in) :: ilut(0:niftot)
    ! origin is the starting point and fixed argument of the input function
    integer, intent(in) :: offset, nTarget, origin
    ! The targetOrbs are the
    integer, intnet(in) :: targetOrbs
    real(dp), intent(out) :: cumSum(nTarget/2)
    ! aux_a and aux_b are additional arguments for the input function, which might
    ! not matter
    integer, intent(in) :: aux_a, aux_b
    logical, intent(in) :: occupancy

    integer :: i, ms
    
! The function of which to create a cumulative list
    interface sf
       function sumFunc(i,j,k,l,ilut)
         use constants, only: dp
         use bit_rep_data, only: NIfTot
         HElement_t(dp) :: sumFunc
         integer, intent(in) :: i,j,k,l
         real(dp) :: ilut(0:NIfTot)
       end function sumFunc
    end interface sf

    ms = mod(origin,2)
    
    if(occupancy .neqv. IsOcc(ilut,offset+2-ms)) then
       cumSum(1) = abs(sumFunc(origin,offset+2-ms,aux_a,aux_b,ilut))
    else
       cumSum(1) = 0.0_dp
    endif
    do i = 2, nTarget/2
       if(occupancy .neqv. IsOcc(ilut,offset+2*i-ms)) then
          cumSum(i) = cumSum(i-1) + abs(sumFunc(origin,offset+2*i-ms,aux_a,aux_b,ilut))
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

  subroutine assign_output_ilut(ilut,ilutnJ,i,j,k,l)
    implicit none
    integer, intent(in) :: i,j
    integer, intent(in), optional :: k,l
    integer(n_int), intent(in) :: ilut(0:niftot)
    integer(n_int), intent(out) :: ilutnJ(0:niftot)

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
    use constants, only: dp
    use bit_rep_data, only: NIfTot
    HElement_t(dp) :: mel
    integer, intent(in) :: i,j,k,l,ilut(NIfTot)
    
    mel = getTMatEl(i,j)
  end function getBathImpMatEl

!------------------------------------------------------------------------------------------!

  function getImpImpMatEl(i,j,k,l,ilut) result(mel)
    use constants, only: dp
    use bit_rep_data, only: NIfTot
    HElement_t(dp) :: mel
    integer, intent(in) :: i,j,k,l,ilut(NIfTot)
    integer :: nOccImp,m

    mel = getTMatEl(i,j)
    if(isImp(i) .and. isImp(j)) then       
       do m = 1, nBath
          if(isOcc(ilut,bathSites(m))) &
          mel = mel + getImpImpMatElDoubs(i,bathSites(m),j,bathSites(m),nI)
       enddo
    end if
  end function getImpImpMatEl

!------------------------------------------------------------------------------------------!

  function getImpImpMatElDoubs(i,j,k,l,ilut) result(mel)
    use constants, only: dp
    use bit_rep_data, only: NIfTot
    HElement_t(dp) :: mel
    integer, intent(in) :: i,j,k,l,ilut(NIfTot)
    
    mel = get_umat_el(i,k,l,j) - get_umat_el(i,k,j,l)

  end function getImpImpMatElDoubs

!------------------------------------------------------------------------------------------!

  function getImpImpMatEl2Ind(i,j,k,l,ilut) result(mel)
    use constants, only: dp
    use bit_rep_data, only: NIfTot
    HElement_t(dp) :: mel
    integer, intent(in) :: i,j,k,l,ilut(NIfTot)
    
    mel = get_umat_el(i,j,i,j) - get_umat_el(i,j,j,i)

  end function getImpImpMatEl2Ind

!------------------------------------------------------------------------------------------!

  ! at the end of the run, clear the impurity excitation generator
  subroutine clearImpurityExcitGen()
    implicit none
    
    ! deallocate the auxiliary resources used by the impurity excitation generator
    deallocate(bathSites)
    deallocate(impuritySites)
    deallocate(connections)
  end subroutine clearImpurityExcitGen

end module impurityModels
