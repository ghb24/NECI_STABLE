#include "macros.h"
module spin_project
    use SystemData, only: LMS, STOT, nel, nbasis, tHPHF
    use CalcData, only: tau, tTruncInitiator, tAlLRealCoeff, &
                        tRealCoeffByExcitLevel, RealCoeffExcitThresh
    use SymExcitDataMod, only: scratchsize
    use bit_rep_data, only: extract_sign
    use bit_reps, only: NIfD, NIfTot, flag_initiator, test_flag, set_flag, &
                        get_initiator_flag_by_run
    use csf, only: csf_get_yamas, get_num_csfs, csf_coeff, random_spin_permute
    use constants, only: dp, bits_n_int, lenof_sign, n_int, end_n_int, int32, sizeof_int, &
                         inum_runs, maxExcit
    use FciMCData, only: TotWalkers, CurrentDets, fcimc_iter_data, &
                         yama_global, excit_gen_store_type, &
                         fcimc_excit_gen_store
    use DeterminantData, only: write_det, get_lexicographic
    use dSFMT_interface, only: genrand_real2_dSFMT
    use util_mod, only: choose, binary_search
    use DetBitOps, only: IsAllowedHPHF, count_open_orbs

    implicit none

    ! Logical(4) datatypes for compilation with builds of openmpi that don't
    ! have support for logical(8). Gah.
    logical :: spin_proj_spawn_initiators, spin_proj_no_death

    ! tSpinProject moved to CalcData to avoid circular dependencies
    logical :: spin_proj_stochastic_yama
    logical :: disable_spin_proj_varyshift
    integer :: spin_proj_interval, spin_proj_cutoff, spin_proj_iter_count
    integer :: spin_proj_nopen_max
    real(dp) :: spin_proj_gamma
    real(dp), target :: spin_proj_shift

    ! Store the data from iterations
    type(fcimc_iter_data), target :: iter_data_spin_proj

    ! For spin projection, pre-calculate and store the Yamanouchi
    ! symbols --> don't need to re-calculate them all the time.
    type yama_storage_type
        integer, allocatable :: yamas(:, :)
        integer :: nyama
    end type
    type(yama_storage_type), allocatable, target :: y_storage(:)

contains

    subroutine init_yama_store()

        ! Calculate all of  the allowed Yamanouchi symbols with the given
        ! values of S, Ms for all allowed unpaired electrons.

        integer :: nopen, ncsf
#ifdef DEBUG_
        character(*), parameter :: this_routine = 'init_yama_store'
#endif

        ! Just in case...
        ASSERT(STOT == LMS)

        ! Allocate the storage super-object.
        ! TODO: This could be a shared memory object with jiggling of bounds.
        allocate(y_storage(LMS:spin_proj_nopen_max))

        ! Loop over all (allowed) numbers of unpaired electrons
        do nopen = LMS, spin_proj_nopen_max, 2

            ! Obtain all of the csfs
            ncsf = get_num_csfs(nopen, STOT)
            y_storage(nopen)%nyama = ncsf
            allocate(y_storage(nopen)%yamas(ncsf, nopen))

            if (ncsf > 0 .and. nopen > 0) then
                call csf_get_yamas(nopen, STOT, y_storage(nopen)%yamas, ncsf)
            end if
        end do

    end subroutine

    subroutine clean_yama_store()

        integer :: i

        if (allocated(y_storage)) then

            do i = lbound(y_storage, 1), ubound(y_storage, 1)
                if (allocated(y_storage(i)%yamas)) &
                    deallocate(y_storage(i)%yamas)
            end do

            deallocate(y_storage)

        end if

    end subroutine

    subroutine test_spin_proj(nI, ilutI)

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:nIfTot)

        integer :: dorder_i(nel), dorder_j(nel), open_el(nel)
        integer :: nopen, i, nup, orb, orb2, count_dets, ndet, pos
        integer(n_int) :: iluttmp(0:niftot)
        real(dp) :: sgnI(lenof_sign), sgnJ(lenof_sign)
        real(dp) :: tot_cpt, elem

        ! Extract the dorder for det nI
        nopen = 0
        i = 1
        nup = 0
        do while (i <= nel)
            if (is_beta(nI(i)) .and. i < nel) then
                if (is_in_pair(nI(i), nI(i + 1))) then
                    i = i + 2
                    cycle
                end if
            end if

            nopen = nopen + 1
            open_el(nopen) = nI(i)
            if (is_alpha(nI(i))) then
                dorder_i(nopen) = 0
                nup = nup + 1
            else
                dorder_i(nopen) = 1
            end if

            i = i + 1
        end do

        if (nopen /= 0) then
            ! How many dets are there to choose from
            ndet = int(choose(nopen, nup), sizeof_int)

            ! initialise iluttmp
            iluttmp = ilutI
            tot_cpt = 0

            ! Get the sign for det I
            call extract_sign(iLutI, sgnI)

            ! Generate the list of all possible determinants, one by one.
            count_dets = 0
            dorder_j(1) = -1
            call get_lexicographic(dorder_j, nopen, nup)
            do while (dorder_j(1) /= -1)

                count_dets = count_dets + 1

                ! Obtain the bit representation of determinant J
                do i = 1, nopen
                    if (dorder_j(i) == 0) then
                        orb = get_alpha(open_el(i))
                    else
                        orb = get_beta(open_el(i))
                    end if
                    orb2 = ab_pair(orb)
                    set_orb(iluttmp, orb)
                    clr_orb(iluttmp, orb2)
                end do

                pos = binary_search(CurrentDets(:, 1:TotWalkers), iLutTmp, &
                                    nIfD + 1)

                call extract_sign(CurrentDets(:, pos), sgnJ)

                elem = csf_spin_project_elem(dorder_i, dorder_j, nopen)
                tot_cpt = tot_cpt + (sgnJ(1) * elem)

                !print*, '   -- cpt', sgnj(1), pos, elem, tot_cpt

                call get_lexicographic(dorder_j, nopen, nup)
            end do

            if (count_dets /= ndet) then
                !print*, ndet, count_dets
                call stop_all("Det count failure", "Here")
            end if

            !print*, 'test cpt', tot_cpt, '(', sgnI(1), ')'
        end if

        !print*, 'test cpt', 'this is closed shell...'

    end subroutine

    subroutine csf_spin_project_one_yama(nI, yama)

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
#ifdef DEBUG_
        character(*), parameter :: this_routine = 'csf_spin_project_one_yama'
#endif

        ! Get the dorder for nI
        nopen = 0
        nup = 0
        i = 1
        do while (i <= nel)
            if (is_beta(nI(i)) .and. i < nel) then
                if (is_in_pair(nI(i), nI(i + 1))) then
                    i = i + 2
                    cycle
                end if
            end if

            nopen = nopen + 1
            if (is_alpha(nI(i))) then
                nup = nup + 1
                dorder_i(nopen) = 0
            else
                dorder_i(nopen) = 1
            end if

            i = i + 1
        end do

        ! <Y|I>
        elem_i = csf_coeff(yama, dorder_i, nopen)

        ! How many dets are there to choose from
        ndet = int(choose(nopen, nup), sizeof_int)

        ! Construct the format string
        write(fmt_num, '(i6)') nopen
        fmt_str = '(a,'//trim(adjustl(fmt_num))//'i1," - ")'

        !print*, 'ELEM I', elem_i

        ! Obtain a list of possible determinants
        count_dets = 0
        dorder_j(1) = -1
        tot_wgt = 0
        tot_sum = 0
        tot_wgt_2 = 0
        tot_sum_2 = 0
        call get_lexicographic(dorder_j, nopen, nup)
        !print*, '     |J>            <J|Y><Y|I>             <J|Y>        &
        !        &     sum_Y <J|Y><Y|I>'
        do while (dorder_j(1) /= -1)

            count_dets = count_dets + 1

            ! <J|Y>
            elem_j = csf_coeff(yama, dorder_j, nopen)
            elem = elem_i * elem_j
            tot_wgt = tot_wgt + abs(elem)
            tot_wgt_2 = tot_wgt_2 + (elem * elem)

            elem_sum = csf_spin_project_elem(dorder_i, dorder_j, nopen)
            tot_sum = tot_sum + abs(elem_sum)
            tot_sum_2 = tot_sum_2 + elem_sum * elem_sum

            write(6, fmt_str, advance='no') 'det: ', dorder_j(1:nopen)
            write(6, *) elem, '(', elem_j, ')', elem_sum

            call get_lexicographic(dorder_j, nopen, nup)
        end do

!        print*, 'total amplitude: ', tot_wgt, tot_sum
!        print*, 'total amplitude squared: ', tot_wgt_2, tot_sum_2

        ASSERT(count_dets /= ndet)

    end subroutine

    function csf_spin_project_elem(dorder_i, dorder_j, nopen) &
        result(ret)

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
        integer, pointer :: yamas(:, :)
        integer :: i, ncsf, ind
        real(dp) :: ret, r

        ! Generate the list of CSFs
        !ncsf = ubound(yamas, 1)
        !call csf_get_yamas (nopen, STOT, yamas, ncsf)
        ncsf = y_storage(nopen)%nyama
        yamas => y_storage(nopen)%yamas

        if (spin_proj_stochastic_yama) then
            r = genrand_real2_dSFMT()
            ind = int(ncsf * r) + 1

            ret = ncsf * (csf_coeff(yamas(ind, :), dorder_i, nopen) * &
                          csf_coeff(yamas(ind, :), dorder_j, nopen))
        else
            ret = 0
            do i = 1, ncsf
                ret = ret + (csf_coeff(yamas(i, :), dorder_i, nopen) * &
                             csf_coeff(yamas(i, :), dorder_j, nopen))
            end do
        end if

    end function

    function csf_spin_project_elem_self(dorder, nopen) result(ret)

        ! Obtain the projection: \sum_Y <I|Y><Y|I>
        ! Used for spin projection.
        !
        ! This is essentially the same as the above routine, but for the
        ! self element, so the sum is simplied.

        integer, intent(in) :: nopen, dorder(nopen)
        integer, pointer :: yamas(:, :)
        integer :: i, ncsf
        real(dp) :: ret, tmp

        ! Generate the list of CSFs
        !ncsf = ubound(yamas, 1)
        !call csf_get_yamas (nopen, STOT, yamas, ncsf)
        ncsf = y_storage(nopen)%nyama
        yamas => y_storage(nopen)%yamas

        ! We cannot use stochastic yama here - otherwise we may end up with
        ! a 0 on the bottom of the term in get_spawn_helement
        !if (spin_proj_stochastic_yama) then
        !    r = genrand_real2_dSFMT()
        !    ind = int(ncsf * r) + 1
        !    tmp = csf_coeff (yamas(ind, :), dorder, nopen)

        !    ret = ncsf * (tmp * tmp)
        !else
        ret = 0
        do i = 1, ncsf
            tmp = csf_coeff(yamas(i, :), dorder, nopen)
            ret = ret + (tmp * tmp)
        end do
         end if

    end function

    function get_spawn_helement_spin_proj(nI, nJ, ilutI, ilutJ, ic, ex, &
                                          tParity, HElGen) result(hel)

        ! Calculate ( - \delta_\gamma \sum_Y <J|Y><Y|I> ) / \delta_\tau
        !
        ! n.b. the negative sign. We need to spawn walkers of the same sign
        !      as the initial det. if element is positive. attempt_create
        !      does the opposite, so we invert the sign of the element.

        integer, intent(in) :: nI(nel), nJ(nel)
        integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer, intent(in) :: ic, ex(2, ic)
        logical, intent(in) :: tParity
        HElement_t(dp), intent(in) :: HElGen
        HElement_t(dp) :: hel

        integer :: iUnused
        integer(n_int) :: iUnused2
        logical :: lUnused
        HElement_t(dp) :: hUnused
#ifdef DEBUG_
        character(*), parameter :: this_routine = 'get_spawn_helement_spin_proj'
#endif

        ! We invert the sign of the returned element, so that the
        ! the attempt_create routine creates walkers of the correct sign.
        hel = csf_spin_project_elem(fcimc_excit_gen_store%dorder_i, &
                                    fcimc_excit_gen_store%dorder_j, &
                                    fcimc_excit_gen_store%nopen)
        hel = -hel * spin_proj_gamma / tau

        ! If we are not permitting death, modify this
        if (spin_proj_no_death) &
            hel = hel / &
                  csf_spin_project_elem_self(fcimc_excit_gen_store%dorder_i, &
                                             fcimc_excit_gen_store%nopen)

        if (tHPHF) then
            ASSERT(count_open_orbs(ilutI) /= 0)
            hel = 2 * hel
        end if

        ! Avoid warnings
        lUnused = tParity; iUnused = IC; iUnused = ex(1, 1)
        iUnused2 = iLutI(0); iUnused2 = iLutJ(0); iUnused = nI(1)
        iUnused = nJ(1); hUnused = helgen

    end function get_spawn_helement_spin_proj

    subroutine generate_excit_spin_proj(nI, iLutI, nJ, iLutJ, exFlag, IC, &
                                        ex, tParity, pGen, HElGen, store, part_type)

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
        integer, intent(out) :: nJ(nel)
        integer(kind=n_int), intent(out) :: iLutJ(0:niftot)
        integer, intent(out) :: ic, ex(2, maxExcit)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        integer :: nopen, nchoose, i
        integer :: nTmp(nel), iUnused
        integer :: open_orbs(nel), open_pos(nel), orb2
        real(dp) :: sgn_tmp(lenof_sign)
        character(*), parameter :: this_routine = 'generate_excit_spin_proj'

        unused_var(part_type)

        ! Only consider determinants with a significant (specified) weight.
        call extract_sign(iLutI, sgn_tmp)
        if (sum(abs(sgn_tmp(1:lenof_sign))) < spin_proj_cutoff) then
            nJ(1) = 0
            return
        end if

        ! TODO: this test should end up somewhere else...
        if (LMS /= STOT) &
            call stop_all(this_routine, "STOT must equal LMS")

        ! Loop over the bit representation to find the unpaired electrons
        ! Can we store the results of this bit?
        ! TODO: we can store the stuff at the starte of nJ now...
        nopen = 0
        i = 1
        do while (i <= nel)
            ! Is this an unpaired electron?
            if (is_beta(nI(i)) .and. i < nel) then
                if (is_in_pair(nI(i), nI(i + 1))) then
                    i = i + 2
                    cycle
                end if
            end if

            ! We have another unpaired electron
            nopen = nopen + 1
            open_orbs(nopen) = nI(i)
            open_pos(nopen) = i

            ! Is it alpha/beta?
            if (is_alpha(nI(i))) then
                fcimc_excit_gen_store%dorder_i(nopen) = 0
            else
                fcimc_excit_gen_store%dorder_i(nopen) = 1
            end if

            ! Loop
            i = i + 1
        end do
        fcimc_excit_gen_store%nopen = nopen

        ! If we know that there are no possible excitations to be made
        if (nopen == STOT .or. nopen > spin_proj_nopen_max .or. &
            (tHPHF .and. nopen == 2)) then
            nJ(1) = 0
            return
        end if

        nTmp(1:nopen) = open_orbs(1:nopen)
        do while (.true.)
            nchoose = random_spin_permute(nTmp(1:nopen), LMS)

            ! Removable for speed?
            if (nchoose == -1) &
                call stop_all(this_routine, "All possible cases here should &
                                             &have been excluded above")

            ! In HPHF, a determinant is allowed if its last e- is alpha
            ! TODO: Is this definition general with > 64 orbitals?
            if (tHPHF) then
                if (is_beta(nTmp(nopen))) cycle
            end if

            ! If we have found our target, exit the loop
            if (.not. all(open_orbs(1:nopen) == nTmp(1:nopen))) exit
        end do

        ! Change the spin structure of nI (only the unpaired elecs)
        ilutJ = ilutI
        nJ = nI
        do i = 1, nopen
            orb2 = ab_pair(nTmp(i))
            set_orb(ilutJ, nTmp(i))
            clr_orb(ilutJ, orb2)
            nJ(open_pos(i)) = nTmp(i)

            ! Construct the dorder --> used in spawn_helement
            if (is_alpha(nTmp(i))) then
                fcimc_excit_gen_store%dorder_j(i) = 0
            else
                fcimc_excit_gen_store%dorder_j(i) = 1
            end if
        end do

        ! If we are in initiator mode, then we may want to make all of the
        ! children into initiators as well
        if (tTruncInitiator) then
            do i = 1, inum_runs
                ! We always want our particles to survive.
                call set_flag(ilutJ, get_initiator_flag_by_run(i))
            end do
        end if

        ! Generation probability, -1 as we exclude the starting det above.
        ! Invert sign so that a positive overlap element spawns walkers with
        ! the same sign
        if (tHPHF) then
            pGen = 2_dp / real(nchoose - 2, dp)
        else
            pGen = 1_dp / real(nchoose - 1, dp)
        end if

        ! Protect against compiler warnings
        ex(1, 1) = ex(1, 1); IC = IC; iUnused = exFlag; tParity = tParity
        HelGen = HelGen; iUnused = store%nopen

    end subroutine generate_excit_spin_proj

    function attempt_die_spin_proj(nI, Kii, RealwSign, WalkExcitLevel, DetPosition) result(ndie)

        integer, intent(in) :: nI(nel)
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        real(dp), intent(in) :: Kii
        real(dp), dimension(lenof_sign) :: ndie
        integer, intent(in) :: WalkExcitLevel
        integer, intent(in), optional :: DetPosition

        real(dp) :: elem, r, rat, rUnused
        integer :: i, iUnused

        unused_var(DetPosition)

        ! If we are not allowing death, or we are below the cutoff for
        ! consideration, then the particle cannot die
        if (spin_proj_no_death .or. &
            sum(abs(realwSign(1:lenof_sign))) < spin_proj_cutoff) then
            ndie = 0.0_dp
            return
        end if

        if (fcimc_excit_gen_store%nopen == STOT .or. &
            fcimc_excit_gen_store%nopen > spin_proj_nopen_max) then
            ndie = 0.0_dp
            return
        end if

        ! Subtract the crurrent value of the shift and multiply by
        ! delta_gamma. If there are multiple particles, scale the
        ! probability.
        elem = csf_spin_project_elem_self(fcimc_excit_gen_store%dorder_i, &
                                          fcimc_excit_gen_store%nopen)
        if (tHPHF) elem = 2 * elem
        elem = elem - 1 + spin_proj_shift
        elem = -elem * spin_proj_gamma

        if (tAllRealCoeff .or. (tRealCoeffByExcitLevel .and. (WalkExcitLevel <= RealCoeffExcitThresh))) then
            ndie(1) = elem * abs(realwSign(1))
        else
            do i = 1, lenof_sign
                rat = elem * abs(realwSign(i))

                ndie(i) = real(int(rat), dp)
                rat = rat - real(ndie(i), dp)
                !print*, 'RAT die', rat

                ! Choose to die or not stochastically
                r = genrand_real2_dSFMT()
                if (abs(rat) > r) ndie(i) = ndie(i) + real(nint(sign(1.0_dp, rat)), dp)
            end do
        end if

        ! Protect against compiler warnings
        rUnused = Kii; iUnused = nI(1)

    end function attempt_die_spin_proj

end module
