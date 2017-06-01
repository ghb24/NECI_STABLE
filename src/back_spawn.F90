#include "macros.h"

module back_spawn

    use CalcData, only: t_back_spawn, tTruncInitiator
    use SystemData, only: nel, nbasis, G1
    use constants, only: n_int, dp
    use bit_rep_data, only: nifd
    use fcimcdata, only: projedet
    use dSFMT_interface, only: genrand_real2_dSFMT

    implicit none

    ! i need a list to indicate the virtual orbitals in the reference 
    ! determinant: the idea of the first implementation is for non-initiators
    ! to only pick electrons from these orbitals, so that the chance to 
    ! de-excite relative to the reference determinant is higher and thus to 
    ! increase the chance to hit already occupied determinants

    ! i could use a mask in the ilut format.. 
    integer(n_int), allocatable :: mask_virt_ilut(:)

    ! or i could use a list of orbitals in the nI format
    integer, allocatable :: mask_virt_ni(:)

    ! and i guess it could also be wise to do it in a spatial resolved way.
    integer, allocatable :: mask_virt_spat(:)

contains

    ! what do i need..
    subroutine init_back_spawn() 
        ! init routine
        character(*), parameter :: this_routine = "init_back_spawn"

        ! first it only makes sense if we actually use the initiator method
        if (.not. tTruncInitiator) then 
            call stop_all(this_routine, &
                "back spawning makes only sense in the initiator method!")
        end if

        ! first use the most simple implementation of an nI style 
        ! virtual orbital indication:
        if (allocated(mask_virt_ni)) deallocate(mask_virt_ni)

        allocate(mask_virt_ni(nBasis - nel))

        ! and assure that this routine is called after the first HFDET is 
        ! already assigned
        ASSERT(allocated(projedet))
        if (.not.allocated(projedet)) then 
            call stop_all(this_routine, &
                "init_back_spawn() called to early; run reference not yet setup!")
        end if

        call setup_virtual_mask()

    end subroutine init_back_spawn

    subroutine setup_virtual_mask()
        ! routine to setup the list of virtual orbitals in the current 
        ! reference determinant, these are then used to choose the electrons
        ! for non-initiator determinants
        ! for now this is only done for single-runs! not dneci, mneci for now!
        character(*), parameter :: this_routine = "setup_virtual_mask"
        integer :: i, j

        ASSERT(allocated(projedet))

        ! i guess the easiest way to do that is to loop over all the 
        ! spin-orbitals and only write an entry if this orbital is not 
        ! occupied in the reference
        j = 1
        do i = 1, nbasis
            ! if (i) is in the reference cycle
            if (any(i == projedet)) cycle
            ! otherwise fill up the virtual mask
            mask_virt_ni(j) = i
            j = j + 1
        end do

    end subroutine setup_virtual_mask

    subroutine pick_virtual_electrons_double(nI, elecs, src, sym_prod, ispn, &
                                                sum_ml, pgen)
        ! this is the important routine! 
        ! for non-initiator determinants this pick electrons only from the 
        ! virtual orbitals of the reference determinant to increase the 
        ! chance to de-excite and to spawn to an already occupied 
        ! determinant from an non-initiator!
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: elecs(2), src(2), sym_prod, ispn, sum_ml
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "pick_virtual_electrons_double"

        integer :: i, n_valid, j, ind, n_valid_pairs, ind_1, ind_2
        integer, allocatable :: virt_elecs(:)

        ! i guess for now i only want to choose uniformly from all the 
        ! available electron in the virtual orbitals of the reference

        ! what do we need here? 
        ! count all the electrons in the virtual of the reference, then 
        ! pick two random orbitals out of those! 
        ! check the routine in symrandexcit3.f90 this does the job i guess..

!         print *, "test picking 2 virtual electrons:"
!         print *, "nI: ", nI
!         print *, "mask_virt_ni: ", mask_virt_ni
        n_valid = 0

        do i = 1, nel
            if (any(nI(i) == mask_virt_ni)) then
                ! the electron is in the virtual of the 
                n_valid = n_valid + 1
            end if
        end do

!         print *, "n_valid: ", n_valid
        if (n_valid < 2) then
            ! something went wrong
            ! in this case i have to abort as no valid double excitation 
            ! could have been found
            elecs = 0
            src = 0
            pgen = 0.0_dp
            return
!             call stop_all(this_routine, & 
!                 "something went wront, did not find 2 valid virtual electrons!")
        end if

        allocate(virt_elecs(n_valid)) 

        j = 1
        do i = 1, nel
            if (any(nI(i) == mask_virt_ni)) then
                virt_elecs(j) = i
                j = j + 1
            end if
        end do

!         print *, "virt_elecs: ", virt_elecs
        ! determine how many valid pairs there are now
        n_valid_pairs = (n_valid * (n_valid - 1)) / 2

        ! and the pgen is now: 
        pgen = 1.0_dp / real(n_valid_pairs, dp)

        ! and is it now enough to do is just like in the symrandexcit3 routine:
        ind = 1 + int(n_valid_pairs * genrand_real2_dSFMT())
        ind_1 = ceiling((1 + sqrt(1 + 8*real(ind,dp))) / 2)
        ind_2 = ind - ((ind_1 - 1) * (ind_1 - 2)) / 2

        ! and retro pick the electron number from the created list? 
        elecs(1) = virt_elecs(ind_1)
        elecs(2) = virt_elecs(ind_2)

!         print *, "ind_1, ind_2: ", ind_1, ind_2

        ! hm.. test this tomorrow
        
        ! now i have to pick two random ones from the list! 
        ! all the symmetry related stuff at the end:
        src = nI(elecs)

!         print *, "src: ", src
        
        if (is_beta(src(1)) .eqv. is_beta(src(2))) then
            if (is_beta(src(1))) then
                iSpn = 1
            else
                iSpn = 3
            end if
        else
            iSpn = 2
        end if

        ! The Ml value is obtained from the orbitals
        sum_ml = sum(G1(src)%Ml)

        ! And the spatial symmetries
!         sym_prod = RandExcitSymLabelProd(SpinOrbSymLabel(src(1)), &
!                                          SpinOrbSymLabel(src(2)))

    end subroutine pick_virtual_electrons_double

    subroutine pick_virtual_electron_single(nI, elec, pgen_elec)
        ! same as above for a single excitation
        ! remember: elec is really just the number in the ilut!
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: elec
        real(dp), intent(out) :: pgen_elec
        character(*), parameter :: this_routine = "pick_virtual_electron_single"

        integer :: i, n_valid, j, ind
        integer, allocatable :: virt_elecs(:)

        ! what do we need here? 
        ! count all the electrons in the virtual of the reference, then 
        ! create a list of them and pick one uniformly
        n_valid = 0
        do i = 1, nel
            if (any(nI(i) == mask_virt_ni)) then
                ! the electron is in the virtual of the 
                n_valid = n_valid + 1
            end if
        end do

        if (n_valid == 0) then
            ! something went wrong
            call stop_all(this_routine, & 
                "something went wront, did not find valid virtual single electron!")
        end if

        allocate(virt_elecs(n_valid)) 

        j = 1
        do i = 1, nel
            if (any(nI(i) == mask_virt_ni)) then
                virt_elecs(j) = i
                j = j + 1
            end if
        end do

        ! and now pick a random number: 
        ind = 1 + floor(genrand_real2_dSFMT() * n_valid) 

        elec = virt_elecs(ind)

#ifdef __DEBUG
        print *, "test picking single virtual elec:"
        print *, "nI: ", nI 
        print *, "mask_virt_ni: ", mask_virt_ni
        print *, "virt_elecs", virt_elecs
        print *, "elec: ", elec
#endif

        pgen_elec = 1.0_dp / real(n_valid, dp)

    end subroutine pick_virtual_electron_single
    
end module back_spawn

