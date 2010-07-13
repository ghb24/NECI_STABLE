#include "macros.h"
module spin_project
    use SystemData, only: LMS, STOT, nel, nbasiS
    use CalcData, only: tau
    use SymExcitDataMod, only: scratchsize
    use bit_reps, only: NIfD, NIfTot
    use csf, only: csf_get_yamas, get_num_csfs, csf_coeff, random_spin_permute
    use constants, only: dp, bits_n_int, lenof_sign, n_int, end_n_int
    use FciMCData, only: TotWalkers, CurrentDets
    use DeterminantData, only: write_det
    use dSFMT_interface, only: genrand_real2_dSFMT

    implicit none

    logical :: tSpinProject
    integer :: spin_proj_interval
    real(dp) :: spin_proj_gamma, spin_proj_shift

contains

    function csf_spin_project_elem (dorder_i, dorder_j, nopen) &
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

    function csf_spin_project_elem_self (dorder, nopen) result (ret)

        ! Obtain the projection: \sum_Y <I|Y><Y|I>
        ! Used for spin projection.
        !
        ! This is essentially the same as the above routine, but for the
        ! self element, so the sum is simplied.

        integer, intent(in) :: nopen, dorder(nopen)
        integer :: yamas (get_num_csfs(nopen, STOT), nopen)
        integer :: i, ncsf
        real(dp) :: ret, tmp

        ! Generate the list of CSFs
        ncsf = ubound(yamas, 1)
        call csf_get_yamas (nopen, STOT, yamas, ncsf)

        ret = 0
        do i = 1, ncsf
            tmp = csf_coeff (yamas(i, :), dorder, nopen)
            ret = ret + (tmp * tmp)
        enddo
    end function
    
    function get_spawn_helement_spin_proj (nI, nJ, ilutI, ilutJ, ic, ex, &
                                         tParity, prob) result (hel)

        ! Calculate ( \delta_\gamma \sum_Y <J|Y><Y|I> ) / \delta_\tau

        integer, intent(in) :: nI(nel), nJ(nel)
        integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer, intent(in) :: ic, ex(2,2)
        logical, intent(in) :: tParity
        real(dp), intent(in) :: prob
        HElement_t :: hel
        
        integer :: dorder_i(nel), dorder_j(nel), nopen, nopen2, i
        character(*), parameter :: this_routine = 'get_spawn_helement_spin_proj'
        ! nJ contains a list of the unpaired electrons in ilutJ
        do i = 1, nel
            if (nJ(i) == -1) exit
            if (is_alpha(nJ(i))) then
                dorder_j(i) = 1
            else
                dorder_j(i) = 0
            endif
        enddo
        nopen = i - 1

        ! We need the dorder for nI, which we don't have as easily...
        ! This is duplacated from generate_excit_spin_proj. Do we need to
        ! preserve nI? If not, we can just pass it in.
        ! TODO: make this be passed  from generate...not recalculated.
        !       or even, can we pointerise the extraction?
        nopen2 = 0
        i = 1
        do while (i <= nel)
            if (nopen2 >= nopen) exit
            if (is_beta(nI(i)) .and. i < nel) then
                if (is_in_pair(nI(i), nI(i+1))) then
                    i = i + 2
                    cycle
                endif
            endif

            nopen2 = nopen2 + 1
            if (is_alpha(nI(i))) then
                dorder_i(nopen2) = 1
            else
                dorder_i(nopen2) = 0
            endif

            i = i + 1
        enddo
        ASSERT(nopen == nopen2)

        hel = csf_spin_project_elem (dorder_i, dorder_j, nopen)
        hel = hel * spin_proj_gamma / tau

    end function get_spawn_helement_spin_proj

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
        character(*), parameter :: this_routine = 'generate_excit_spin_proj'

        ! TODO: this test should end up somewhere else...
        if (LMS /= STOT) &
            call stop_all (this_routine, "STOT must equal LMS")

        ! Loop over the bit representation to find the unpaired electrons
        ! Can we store the results of this bit?
        ! TODO: we can store the stuff at the starte of nJ now...
        nopen = 0
        i = 1
        do while (i <= nel)
            ! Is this an unpaired electron?
            if (is_beta(nI(i)) .and. i < nel) then
                if (is_in_pair(nI(i), nI(i+1))) then
                    i = i + 2
                    cycle
                endif
            endif

            nopen = nopen + 1
            nJ(nopen) = nI(i)
            i = i + 1
        enddo

        ! If we know that there are no possible excitations to be made
        if (nopen == STOT) then
            nJ(1) = 0
            return
        endif

        nTmp(1:nopen) = nJ(1:nopen)
        do while (.true.)
            nchoose = random_spin_permute (nJ(1:nopen), LMS)

            ! Removable for speed?
            if (nchoose == -1) &
                call stop_all (this_routine, "All possible cases here should &
                                             &have been excluded above")

            ! If we have found our target, exit the loop
            if (.not.all(nJ(1:nopen) == nTmp(1:nopen))) exit
        enddo

        ! Change the spin structure of nI (only the unpaired elecs)
        ilutJ = ilutI
        do i = 1, nopen
            set_orb(ilutJ, nJ(i))
            clr_orb(ilutJ, ab_pair(nJ(i)))
        enddo

        ! Mark the end of the unpaired electrons section.
        if (nopen < nel) nJ(nopen + 1) = -1

        ! Generation probability, -1 as we exclude the starting det above.
        pGen = 1_dp / real(nchoose - 1, dp)
    end subroutine generate_excit_spin_proj

    function attempt_die_spin_proj (nI, Kii, wSign) result (ndie)

        integer, intent(in) :: nI(nel)
        integer, dimension(lenof_sign), intent(in) :: wSign
        real(dp), intent(in) :: Kii
        integer, dimension(lenof_sign) :: ndie

        real(dp) :: elem, r, rat
        integer :: dorder(nel), i, nopen

        ! Extract the dorder for determinant nI
        nopen = 0
        i = 1
        do while (i <= nel)
            if (is_beta(nI(i)) .and. i < nel) then
                if (is_in_pair(nI(i), nI(i+1))) then
                    i = i + 2
                    cycle
                endif
            endif

            nopen = nopen + 1
            if (is_alpha(nI(i))) then
                dorder(nopen) = 1
            else
                dorder(nopen) = 0
            endif

            i = i + 1
        enddo

        ! Subtract the crurrent value of the shift and multiply by
        ! delta_gamma. If there are multiple particles, scale the
        ! probability.
        elem = csf_spin_project_elem_self (dorder, nopen)
        elem = elem - 1 + spin_proj_shift
        elem = - elem * spin_proj_gamma

        do i = 1, lenof_sign
            rat = elem * abs(wSign(i))
            
            ndie(i) = int(rat)
            rat = rat - real(ndie(i), dp)

            ! Choose to die or not stochastically
            r = genrand_real2_dSFMT()
            if (abs(rat) > r) ndie(i) = ndie(i) + nint(sign(1.0_dp, rat))
        enddo

    end function attempt_die_spin_proj


end module
