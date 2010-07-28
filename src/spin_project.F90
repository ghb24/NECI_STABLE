#include "macros.h"
module spin_project
    use SystemData, only: LMS, STOT, nel, nbasis
    use CalcData, only: tau
    use SymExcitDataMod, only: scratchsize
    use bit_reps, only: NIfD, NIfTot, extract_sign
    use csf, only: csf_get_yamas, get_num_csfs, csf_coeff, random_spin_permute
    use constants, only: dp, bits_n_int, lenof_sign, n_int, end_n_int, int32
    use FciMCData, only: TotWalkers, CurrentDets, fcimc_iter_data, yama_global
    use DeterminantData, only: write_det
    use dSFMT_interface, only: genrand_real2_dSFMT
    use util_mod, only: choose, binary_search

    implicit none

    logical :: tSpinProject, spin_proj_stochastic_yama
    integer :: spin_proj_interval, spin_proj_cutoff
    real(dp) :: spin_proj_gamma, spin_proj_shift

    ! Store the data from iterations
    type(fcimc_iter_data), target :: iter_data_spin_proj

contains

    subroutine test_spin_proj (nI, ilutI)

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:nIfTot)

        integer :: dorder_i(nel), dorder_j(nel), open_el(nel)
        integer :: nopen, i, nup, orb, count_dets, ndet, pos
        integer(n_int) :: iluttmp(0:niftot)
        integer, dimension(lenof_sign) :: sgnI, sgnJ
        real*8 :: tot_cpt, elem

        ! Extract the dorder for det nI
        nopen = 0
        i = 1
        nup = 0
        do while (i <= nel)
            if (is_beta(nI(i)) .and. i < nel) then
                if (is_in_pair(nI(i), nI(i+1))) then
                    i = i + 2
                    cycle
                endif
            endif

            nopen = nopen + 1
            open_el(nopen) = nI(i)
            if (is_alpha(nI(i))) then
                dorder_i(nopen) = 0
                nup = nup + 1
            else
                dorder_i(nopen) = 1
            endif

            i = i + 1
        enddo

        if (nopen /= 0) then
            ! How many dets are there to choose from
            ndet = choose (nopen, nup)

            ! initialise iluttmp
            iluttmp = ilutI
            tot_cpt = 0

            ! Get the sign for det I
            call extract_sign (iLutI, sgnI)

            ! Generate the list of all possible determinants, one by one.
            count_dets = 0
            dorder_j(1) = -1
            call get_lexicographic (dorder_j, nopen, nup)
            do while (dorder_j(1) /= -1)

                count_dets = count_dets + 1

                ! Obtain the bit representation of determinant J
                do i = 1, nopen
                    if (dorder_j(i) == 0) then
                        orb = get_alpha(open_el(i))
                    else
                        orb = get_beta(open_el(i))
                    endif
                    set_orb(iluttmp, orb)
                    clr_orb(iluttmp, ab_pair(orb))
                enddo

                pos = binary_search (CurrentDets, iLutTmp, nIfTot+1, &
                                     int(TotWalkers,int32), nIfD+1)

                call extract_sign (CurrentDets(:,pos), sgnJ)

                elem =  csf_spin_project_elem (dorder_i, dorder_j, nopen)
                tot_cpt = tot_cpt + (sgnJ(1) * elem)

                !print*, '   -- cpt', sgnj(1), pos, elem, tot_cpt

                call get_lexicographic (dorder_j, nopen, nup)
            enddo

            if (count_dets /= ndet) then
                !print*, ndet, count_dets
                call stop_all ("Det count failure", "Here")
            endif

            !print*, 'test cpt', tot_cpt, '(', sgnI(1), ')'
        endif

        !print*, 'test cpt', 'this is closed shell...'

    end subroutine

    subroutine csf_spin_project_one_yama (nI, yama)

        ! Apply the operator Os to the determinant nI (ilutI) using only the
        ! csf Yama, and print out the resultant coefficients.
        !
        ! In: nI    - The determinant to apply Os to
        !     yama  - The Yamanouchi symbol to use

        integer, intent(in) :: nI(nel), yama(:)

        integer :: nopen, nup, count_dets, ndet, i
        integer :: dorder_i(nel), dorder_j(nel)
        real(dp) :: elem, elem_i, elem_j, tot_wgt, elem_sum, tot_sum
        real(dp) :: tot_wgt_2, tot_sum_2
        character(20) :: fmt_str, fmt_num
        character(*), parameter :: this_routine = 'csf_spin_project_one_yama'

        ! Get the dorder for nI
        nopen = 0
        nup = 0
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
                nup = nup + 1
                dorder_i(nopen) = 0
            else
                dorder_i(nopen) = 1
            endif

            i = i + 1
        enddo

        ! <Y|I>
        elem_i = csf_coeff (yama, dorder_i, nopen)

        ! How many dets are there to choose from
        ndet = choose (nopen, nup)

        ! Construct the format string
        write(fmt_num, '(i6)') nopen
        fmt_str = '(a,'//trim(adjustl(fmt_num))//'i1," - ")'

        print*, 'ELEM I', elem_i

        ! Obtain a list of possible determinants
        count_dets = 0
        dorder_j(1) = - 1
        tot_wgt = 0
        tot_sum = 0
        tot_wgt_2 = 0
        tot_sum_2 = 0
        call get_lexicographic (dorder_j, nopen, nup)
        print*, '     |J>            <J|Y><Y|I>             <J|Y>        &
                &     sum_Y <J|Y><Y|I>'
        do while (dorder_j(1) /= -1)

            count_dets = count_dets + 1

            ! <J|Y>
            elem_j = csf_coeff (yama, dorder_j, nopen)
            elem = elem_i * elem_j
            tot_wgt = tot_wgt + abs(elem)
            tot_wgt_2 = tot_wgt_2 + (elem*elem)

            elem_sum = csf_spin_project_elem (dorder_i, dorder_j, nopen)
            tot_sum = tot_sum + abs(elem_sum)
            tot_sum_2 = tot_sum_2 + elem_sum*elem_sum

            write (6, fmt_str, advance='no') 'det: ', dorder_j(1:nopen)
            write (6,*) elem, '(', elem_j, ')', elem_sum

            call get_lexicographic (dorder_j, nopen, nup)
        enddo

        print*, 'total amplitude: ', tot_wgt, tot_sum
        print*, 'total amplitude squared: ', tot_wgt_2, tot_sum_2

        ASSERT(count_dets /= ndet)

    end subroutine

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
        integer :: i, ncsf, ind
        real(dp) :: ret, r

        ! Generate the list of CSFs
        ncsf = ubound(yamas, 1)
        call csf_get_yamas (nopen, STOT, yamas, ncsf)

        if (spin_proj_stochastic_yama) then
            r = genrand_real2_dSFMT()
            ind = int(ncsf * r) + 1

            ret = ncsf * (csf_coeff (yamas(ind, :), dorder_i, nopen) * &
                          csf_coeff (yamas(ind, :), dorder_j, nopen))
        else
            ret = 0
            do i = 1, ncsf
                ret = ret + (csf_coeff (yamas(i, :), dorder_i, nopen) * &
                             csf_coeff (yamas(i, :), dorder_j, nopen))
            enddo
        endif

    end function

    function csf_spin_project_elem_self (dorder, nopen) result (ret)

        ! Obtain the projection: \sum_Y <I|Y><Y|I>
        ! Used for spin projection.
        !
        ! This is essentially the same as the above routine, but for the
        ! self element, so the sum is simplied.

        integer, intent(in) :: nopen, dorder(nopen)
        integer :: yamas (get_num_csfs(nopen, STOT), nopen)
        integer :: i, ncsf, ind
        real(dp) :: ret, tmp, r

        ! Generate the list of CSFs
        ncsf = ubound(yamas, 1)
        call csf_get_yamas (nopen, STOT, yamas, ncsf)

        if (spin_proj_stochastic_yama) then
            r = genrand_real2_dSFMT()
            ind = int(ncsf * r) + 1
            tmp = csf_coeff (yamas(ind, :), dorder, nopen)

            ret = ncsf * (tmp * tmp)
        else
            ret = 0
            do i = 1, ncsf
                tmp = csf_coeff (yamas(i, :), dorder, nopen)
                ret = ret + (tmp * tmp)
            enddo
        endif
        
    end function
    
    function get_spawn_helement_spin_proj (nI, nJ, ilutI, ilutJ, ic, ex, &
                                         tParity, prob) result (hel)

        ! Calculate ( - \delta_\gamma \sum_Y <J|Y><Y|I> ) / \delta_\tau
        !
        ! n.b. the negative sign. We need to spawn walkers of the same sign
        !      as the initial det. if element is positive. attempt_create
        !      does the opposite, so we invert the sign of the element.

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
                dorder_j(i) = 0
            else
                dorder_j(i) = 1
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
                dorder_i(nopen2) = 0
            else
                dorder_i(nopen2) = 1
            endif

            i = i + 1
        enddo
        ASSERT(nopen == nopen2)

        ! We invert the sign of the returned element, so that the
        ! the attempt_create routine creates walkers of the correct sign.
        hel = csf_spin_project_elem (dorder_i, dorder_j, nopen)
        hel = - hel * spin_proj_gamma / tau

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
        integer, dimension(lenof_sign) :: sgn_tmp
        character(*), parameter :: this_routine = 'generate_excit_spin_proj'

        call extract_sign (iLutI, sgn_tmp)
        if (sum(abs(sgn_tmp(1:lenof_sign))) < spin_proj_cutoff) then
            nJ(1) = 0
            return
        endif

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
        ! Invert sign so that a positive overlap element spawns walkers with
        ! the same sign
        pGen = 1_dp / real(nchoose - 1, dp)
        !print*, 'Generation prob', pGen, nchoose-1
    end subroutine generate_excit_spin_proj

    function attempt_die_spin_proj (nI, Kii, wSign) result (ndie)

        integer, intent(in) :: nI(nel)
        integer, dimension(lenof_sign), intent(in) :: wSign
        real(dp), intent(in) :: Kii
        integer, dimension(lenof_sign) :: ndie

        real(dp) :: elem, r, rat
        integer :: dorder(nel), i, nopen

        if (sum(abs(wSign(1:lenof_sign))) < spin_proj_cutoff) then
            ndie = 0
            return
        endif

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
                dorder(nopen) = 0
            else
                dorder(nopen) = 1
            endif

            i = i + 1
        enddo

        if (nopen == STOT) then
            ndie = 0
            return
        endif

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
            !print*, 'RAT die', rat

            ! Choose to die or not stochastically
            r = genrand_real2_dSFMT()
            if (abs(rat) > r) ndie(i) = ndie(i) + nint(sign(1.0_dp, rat))
        enddo

    end function attempt_die_spin_proj

    subroutine get_lexicographic (dorder, nopen, nup)

        ! Unlike the csf version, this uses 1 == alpha, 0 = beta.

        integer, intent(in) :: nopen, nup
        integer, intent(inout) :: dorder(nopen)
        integer :: comb(nup)
        integer :: i, j
        logical :: bInc

        ! Initialise
        if (dorder(1) == -1) then
            dorder(1:nup) = 0
            dorder(nup+1:nopen) = 1
        else
            ! Get the list of positions of the beta electrons
            j = 0
            do i = 1, nopen
                if (dorder(i) == 0) then
                    j = j + 1
                    comb(j) = i
                    
                    ! Have we reached the last possibility?
                    if (j == 1 .and. i == nopen - nup + 1) then
                        dorder(1) = -1
                        return
                    endif

                    if (j == nup) exit
                endif
            enddo

            do i = 1, nup
                bInc = .false.
                if (i == nup) then
                    bInc = .true.
                else if (i < nup) then
                    if (comb(i+1) /= comb(i) + 1) bInc = .true.
                endif

                if (bInc) then
                    comb(i) = comb(i) + 1
                    exit
                else
                    comb(i) = i
                endif
            enddo

            dorder = 1
            dorder(comb) = 0
        endif
    end subroutine


end module
