#include "macros.h"
module spin_project
	use SystemData, only: LMS, STOT, nel, nbasis
    use SymExcitDataMod, only: scratchsize
    use bit_reps, only: NIfD, NIfTot
	use csf, only: csf_get_yamas, get_num_csfs, csf_coeff, csf_get_random_det
    use constants, only: dp, bits_n_int, lenof_sign, n_int, end_n_int
    use FciMCData, only: TotWalkers, CurrentDets

	implicit none

contains

    pure function csf_spin_project_elem (dorder_i, dorder_j, nopen) &
                                        result (ret)

        ! Obtain the projection: \sum_Y <J|Y><Y|I>
        ! Used for spin projection.
        !
        ! n.b. I, J are determinants, and the sum performed is over all CSFs
        !      with the same spatial structure, with the overall specified
        !      value of S, Ms.
        !
        ! In: dorder_i, dorder_j - The list of alpha/beta for each spin
        !                          orbital in det i,j (only unpaired orbs)
        !     nopen              - The number of unpaired electrons

        integer, intent(in) :: nopen, dorder_i(nopen), dorder_j(nopen)
        integer :: yamas (get_num_csfs(nopen, STOT), nopen)
        integer :: i, ncsf
        real(dp) :: ret

        ! Generate the list of CSFs
        ncsf = ubound(yamas, 1)
        call csf_get_yamas (nopen, STOT, yamas, ncsf)

        ret = 0
        do i = 1, ncsf
            ret = ret + (csf_coeff (yamas(i, :), dorder_i, nopen) * &
                         csf_coeff (yamas(i, :), dorder_j, nopen))
        enddo
    end function

    subroutine generate_excit_spin_proj (nI, iLutI, nJ, iLutJ, exFlag, IC, &
                                         ex, tParity, pGen, tFilled, &
                                         scratch1, scratch2, scratch3)

        ! This returns an excitation of the source determiant (iLutI).
        !
        !   --> It does not need nI.
        !   --> New det returned in nJ
        !   --> List of unpaired electrons (in reverse order) is returned
        !       in last nopen elements of nJ.
        !   --> The element of nJ before the unpaired electrons start is
        !       marked with -1 --> don't have to count nopen again.
        integer, intent(in) :: nI(nel)
        integer(kind=n_int), intent(in) :: iLutI(0:niftot)
        integer, intent(in) :: exFlag
        integer, intent(inout) :: scratch1(scratchsize)
        integer, intent(inout) :: scratch2(scratchsize)
        integer, intent(inout) :: scratch3(scratchsize)
        integer, intent(out) :: nJ(nel) 
        integer(kind=n_int), intent(out) :: iLutJ(0:niftot)
        integer, intent(out) :: ic, ex(2,2)
        real(dp), intent(out) :: pGen
        logical, intent(inout) :: tFilled
        logical, intent(out) :: tParity

        integer :: nopen, nchoose, i, j
        integer :: nTmp(nel)
        logical :: found
        character(*), parameter :: this_routine = 'generate_excit_spin_proj'

        ! TODO: this test should end up somewhere else...
        if (LMS /= STOT) &
            call stop_all (this_routine, "STOT must equal LMS")

        ! Loop over the bit representation to find the unpaired electrons
        ! Can we store the results of this bit?
        nopen = 0
        do i = 1, nbasis-1, 2
            ! Is this an unpaired electron?
            if (IsOcc(ilutI, i) .neqv. IsOcc(ilutI, i+1)) then
                if (IsOcc(ilutI, i)) then
                    nJ(nel-nopen) = i
                else
                    nJ(nel-nopen) = i + 1
                endif
                nopen = nopen + 1
            endif
        enddo

        ! If we know that there are no possible excitations to be made
        if (nopen == STOT) then
            nJ(1) = 0
            return
        endif

        nTmp(nel-nopen+1:nel) = nJ(nel-nopen+1:nel)
        do while (.not. found)
            nchoose = csf_get_random_det (nJ, nopen, LMS)

            ! Removable for speed?
            if (nchoose == -1) &
                call stop_all (this_routine, "All possible cases here should &
                                             &have been excluded above")

            if (.not.all(nJ(nel-nopen+1:nel) == nTmp(nel-nopen+1:nel))) exit
        enddo

        ! Change the spin structure of nI (only the unpaired elecs)
        ilutJ = ilutI
        do i = nel-nopen+1, nel
            set_orb(ilutJ, nJ(i))
            clr_orb(ilutJ, ab_pair(nJ(i)))
        enddo

        ! Mark the end of the unpaired electrons section.
        if (nopen < nel) nJ(nel-nopen) = -1

        ! Generation probability, -1 as we exclude the starting det above.
        pGen = 1_dp / real(nchoose - 1, dp)
    end subroutine


















    !subroutine stochastic_spin_project ()
    !    
    !    integer :: det(nel), flags
    !    integer, dimension(lenof_sign) :: sgn
    !    integer :: j, p
    !    character(*), parameter :: this_routine = 'stochastic_spin_project'

    !    if (lenof_sign /= 1) call stop_all (this_routine, &
    !               "Spin projection does not work with complex walkers yet.")

    !    ! Loop over all the determinants in the system.
    !    ! n.b. TotWalkers /= number of walkers in the system.
    !    do j = 1, TotWalkers

    !        ! Decode determinant from (stored) bit-representation
    !        ! call extract_sign (CurrentDets(:,j), sgn)
    !        call extract_bit_rep (CurrentDets(:,j), det, sgn, flags)

    !        ! Loop over all the particles
    !        do p = 1, abs(sgn(1))

    !            ! Generate excitation

    !            ! Test overlap and spawn

    !        enddo
    !        

    !    enddo

    !    ! Annihilate

    !end subroutine
end module
