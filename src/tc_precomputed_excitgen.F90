#include "macros.h"
module tc_precomputed_excitgen
  use constants
  use tc_three_body_data
  use aliasSampling
  use bit_reps, only: niftot
  use SystemData, only: nel
  use FciMCData, only: excit_gen_store_type, pDoubles, pSingles
  use dSFMT_interface , only : genrand_real2_dSFMT
  implicit none

  ! these are the samplers used for generating single excitations
  ! we have one hole sampler for each electron
  type(aliasSampler_t), allocatable :: single_hole_samplers(:)
  ! we have a single electron sampler, the electron is always chosen first
  type(aliasSampler_t) :: single_elec_sampler

contains

  !------------------------------------------------------------------------------------------!
  ! Main routines generating the excitation
  !------------------------------------------------------------------------------------------!

  ! this is the interface routine that calls the corresponding generator to create
  ! a single/double or potentially triple excitation
  subroutine gen_ran_excit_pcpp(nI, ilut, nJ, ilutnJ, exFlag, ic, ExcitMat, tParity, pGen, &
       HElGen, store, part_type)
    implicit none
    ! The interface is common to all excitation generators, see proc_ptrs.F90
    integer, intent(in) :: nI(nel), exFlag
    integer(n_int), intent(in) :: iLut(0:niftot)
    integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(out) :: pGen
    type(excit_gen_store_type), intent(inout), target :: store

    ! Not used
    integer(n_int), intent(out) :: ilutnJ(0:niftot)
    HElement_t(dp), intent(out) :: HElGen

    integer, intent(in), optional :: part_type

    real(dp) :: r

    ! decide whether to generate a single or double excitation
    r = genrand_real2_dSFMT()
    if(r < pDoubles) then
       call generate_double_pcpp(nI, ilut, nJ, excitMat, tParity, pGen)
       IC = 2
       pGen = pGen * pDoubles
    else
       call generate_single_pcpp(nI, ilut, nJ, excitMat, tParity, pGen)
       IC = 1
       pGen = pGen * pSingles
    endif
    
  end subroutine gen_ran_excit_pcpp

  !------------------------------------------------------------------------------------------!

  subroutine generate_double_pcpp(nI, ilut, nJ, excitMat, tParity, pGen)
    implicit none
    ! given the initial determinant (both as nI and ilut), create a random single excitation
    ! given by nJ/ilutnJ/excitMat with probability pGen. tParity indicates the fermi sign
    ! picked up by applying the excitation operator
    integer, intent(in) :: nI(nel)
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    integer, intent(out) :: nJ(nel)
    integer, intent(out) :: ExcitMat(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(out) :: pGen

    
  end subroutine generate_double_pcpp

  !------------------------------------------------------------------------------------------!

  subroutine generate_single_pcpp(nI, ilut, nJ, excitMat, tParity, pGen)
    implicit none
    ! given the initial determinant (both as nI and ilut), create a random double excitation
    ! given by nJ/ilutnJ/excitMat with probability pGen. tParity indicates the fermi sign
    ! picked up by applying the excitation operator
    integer, intent(in) :: nI(nel)
    integer(n_int), intent(in) :: ilut(0:NIfTot)
    integer, intent(out) :: nJ(nel)
    integer, intent(out) :: ExcitMat(2,2)
    logical, intent(out) :: tParity
    real(dp), intent(out) :: pGen

    integer :: src, tgt
    real(dp) :: pHole

    ! get a random electron
    call single_elec_sampler%sample(src,pGen)
    ! get a random associated orbital
    call single_hole_samplers(src)%sample(tgt,pHole)

    ! map the electron/hole to the current determinant
    src = map_elec_from_ref(nI, src)
    tgt = map_orb_from_ref(nI, tgt)

    call make_single(nI, nJ, src, tgt, excitMat, tParity)
    ! add the probability to find this hole from this electron
    pGen = pGen * pHole
    
  end subroutine generate_single_pcpp

  !------------------------------------------------------------------------------------------!
  ! Functions that map orbital and electron indices between reference and current determinant
  !------------------------------------------------------------------------------------------!  

  function map_elec_from_ref(nI, iElec) result(src)
    implicit none
    integer, intent(in) :: nI(nel)
    integer, intent(in) :: iElec
    integer :: src
  end function map_elec_from_ref

  function map_orb_from_ref(nI, iOrb) result(tgt)
    implicit none
    integer, intent(in) :: nI(nel)
    integer, intent(in) :: iOrb
    integer :: tgt
  end function map_orb_from_ref
  
  !------------------------------------------------------------------------------------------!
  ! Initialization routines for the pcpp excitation generator
  !------------------------------------------------------------------------------------------!

  subroutine init_pcpp_excitgen()
    implicit none

    call init_pcpp_doubles_excitgen()
    call init_pcpp_singles_excitgen()
   
  end subroutine init_pcpp_excitgen
  
  !------------------------------------------------------------------------------------------!

  subroutine init_pcpp_doubles_excitgen()
    implicit none
  end subroutine init_pcpp_doubles_excitgen

  !------------------------------------------------------------------------------------------!

  subroutine init_pcpp_singles_excitgen()
    implicit none
  end subroutine init_pcpp_singles_excitgen
  
end module tc_precomputed_excitgen
