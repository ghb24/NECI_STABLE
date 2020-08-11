#include "macros.h"
#:include "macros.fpph"

module hist

    use DeterminantData, only: get_lexicographic, calculated_ms
    use MemoryManager
    use SystemData, only: nbasis, nel, LMS, LMS, tHPHF, &
                          tOddS_HPHF, G1, tGUGA
    use DetBitOps, only: count_open_orbs, EncodeBitDet, spatial_bit_det, &
                         DetBitEq, count_open_orbs, TestClosedShellDet, &
                         CalcOpenOrbs, IsAllowedHPHF, FindBitExcitLevel
    use load_balance_calcnodes, only: DetermineDetNode
    use CalcData, only: tFCIMC, tTruncInitiator
    use DetCalcData, only: FCIDetIndex, det
    use FciMCData, only: tFlippedSign, TotWalkers, CurrentDets, iter, &
                         all_norm_psi_squared, ilutRef
    use HPHFRandExcitMod, only: FindExcitBitDetSym
    use hphf_integrals, only: hphf_sign
    use constants, only: n_int, bits_n_int, size_n_int, lenof_sign
    use bit_rep_data, only: NIfTot, NIfD, extract_sign
    use bit_reps, only: encode_sign, extract_bit_rep, &
                        decode_bit_det, flag_initiator, test_flag, &
                        get_initiator_flag, &
                        any_run_is_initiator
    use parallel_neci
    use searching, only: BinSearchParts2
    use hist_data
    use timing_neci
    use Determinants, only: write_det
    use util_mod

    implicit none

contains

    subroutine ilut_nifd_pointer_assign(ptr, ilut)

        ! This is a useful helper function to get around some irritating
        ! behaviour (namely that pointer array slices are indexed to
        ! begin at 0, not 1 --> wrong for iluts).

        integer(n_int), intent(in), target :: ilut(0:NIfTot)
        integer(n_int), intent(out), pointer :: ptr(:)

        ptr => ilut

    end subroutine

    subroutine add_hist_spawn(ilut, sign, ExcitLevel, dProbFin)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: ExcitLevel
        real(dp), intent(in) :: dProbFin
        real(dp), intent(in) :: sign(lenof_sign)

        integer :: PartInd, open_orbs
        integer(n_int) :: ilut_sym(0:NIfTot)
        real(dp) :: delta(lenof_sign)
        logical :: tSuccess
        character(*), parameter :: t_r = 'add_hist_spawn'

        if (ExcitLevel == nel) then
            call BinSearchParts2(ilut, HistMinInd(ExcitLevel), det, PartInd, &
                                 tSuccess)
            ! CCMC doesn't sum particle contributions in order, so we must
            ! search the whole space again!
            if (tFCIMC) HistMinInd(ExcitLevel) = PartInd
        else if (ExcitLevel == 0) then
            PartInd = 1
            tSuccess = .true.
        else
            call BinSearchParts2(ilut, HistMinInd(ExcitLevel), &
                                 FCIDetIndex(ExcitLevel + 1) - 1, PartInd, &
                                 tSuccess)
        end if

        if (tSuccess) then
            delta = sign / dProbFin

            if (tHPHF) then
                call FindExcitBitDetSym(ilut, ilut_sym)
                if (.not. DetBitEq(ilut, ilut_sym)) delta = delta / sqrt(2.0_dp)
            end if

            if (tFlippedSign) delta = -delta

            Histogram(:, PartInd) = Histogram(:, PartInd) + delta
            if (tHistSpawn) &
                InstHist(:, PartInd) = InstHist(:, PartInd) + delta

            ! In HPHF we also include the spin-coupled determinant, which will
            ! have the same amplitude as the original determinant, unless it
            ! is antisymmetric.
            if (tHPHF .and. .not. DetBitEQ(ilut, ilut_sym)) then
                if (ExcitLevel == nel) then
                    call BinSearchParts2(ilut_sym, FCIDetIndex(ExcitLevel), &
                                         Det, PartInd, tSuccess)
                else if (ExcitLevel == 0) then
                    PartInd = 1
                    tSuccess = .true.
                else
                    call BinSearchParts2(ilut_sym, FCIDetIndex(ExcitLevel), &
                                         FCIDetIndex(ExcitLevel + 1) - 1, &
                                         PartInd, tSuccess)
                end if
                if (tSuccess) then
                    delta = (sign / sqrt(2.0_dp)) / dProbFin

                    call CalcOpenOrbs(ilut_sym, open_orbs)
                    if ((mod(open_orbs, 2) == 1) .neqv. tOddS_HPHF) &
                        delta = -delta

                    if (tFlippedSign .neqv. tOddS_HPHF) delta = -delta

                    Histogram(:, PartInd) = Histogram(:, PartInd) + delta
                    if (tHistSpawn) &
                        InstHist(:, PartInd) = InstHist(:, PartInd) + delta
                end if
            end if
        else
            call writebitdet(6, ilut, .true.)
            write(6, *) '***', ilut
            write(6, *) '***', ExcitLevel, HistMinInd(ExcitLevel), Det
            call stop_all(t_r, "Cannot find corresponding FCI determinant &
                                &when histogramming")
        end if

    end subroutine

    subroutine find_hist_coeff_explicit(ilut, ExcitLevel, PartInd, tSuccess)

        implicit none
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: ExcitLevel
        integer, intent(out) :: PartInd
        logical, intent(out) :: tSuccess

        tSuccess = .false.
        if (ExcitLevel == nel) then
            call BinSearchParts2(ilut, FCIDetIndex(ExcitLevel), det, PartInd, &
                                 tSuccess)
        else if (ExcitLevel == 0) then
            PartInd = 1
            tSuccess = .true.
        else
            call BinSearchParts2(ilut, FCIDetIndex(ExcitLevel), &
                                 FCIDetIndex(ExcitLevel + 1) - 1, PartInd, &
                                 tSuccess)
        end if

    end subroutine

    subroutine add_hist_energies(ilut, Sign, HDiag)

        ! This will histogram the energies of the particles, rather than the
        ! determinants themselves.

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp), dimension(lenof_sign), intent(in) :: Sign
        real(dp), intent(in) :: HDiag
        integer :: bin
        character(*), parameter :: t_r = "add_hist_energies"

        integer(n_int) :: iUnused

        bin = int(HDiag / BinRange) + 1
        if (bin > iNoBins) &
            call stop_all(t_r, "Histogramming energies higher than the &
                         &arrays can cope with. Increase iNoBins or BinRange")
        HistogramEnergy(bin) = HistogramEnergy(bin) + sum(abs(Sign))

        ! Avoid compiler warnings
        iUnused = ilut(0)

    end subroutine

    function calc_s_squared(only_init) result(ssq)

        ! Calculate the instantaneous value of S^2 for the walkers stored
        ! in CurrentDets
        !
        ! --> This could be generalised to an arbitrary list of iluts. We
        !     would also then need to calculate the value of psi_squared

        real(dp) :: ssq(inum_runs), tmp(inum_runs)
        integer :: i, run
        logical, intent(in) :: only_init
        type(timer), save :: s2_timer

        s2_timer%timer_name = 'S2 local'
        call set_timer(s2_timer)

        ssq = 0
        do i = 1, int(TotWalkers, sizeof_int)
            if ((test_flag(CurrentDets(:, i), get_initiator_flag(1)) .or. &
                 test_flag(CurrentDets(:, i), get_initiator_flag(lenof_sign))) &
                .and. .not. TestClosedShellDet(CurrentDets(:, i))) then
                ssq = ssq + ssquared_contrib(CurrentDets(:, i), only_init)
            end if
        end do

        ! Sum over all processors and normalise
        call MPISum(ssq, tmp)
        ssq = tmp

        do run = 1, inum_runs
            if (abs(all_norm_psi_squared(run)) < 1.0e-10_dp) then
                ssq(run) = 0.0_dp
            else

                ssq(run) = ssq(run) / all_norm_psi_squared(run)

                ! TODO: n.b. This is a hack. LMS appears to contain -2*Ms of the
                !            system I am somewhat astounded I haven't noticed this
                !            before...
                ssq(run) = ssq(run) &
                           + real(calculated_ms * (calculated_ms + 2), dp) / 4
            end if
        end do

        call halt_timer(s2_timer)

    end function

    function calc_s_squared_multi() result(ssq)

        integer :: max_linked, max_per_proc, max_spawned
        real(dp), dimension(inum_runs) :: ssq
        type(timer), save :: s2_timer

        s2_timer%timer_name = 'S^2'
        call set_timer(s2_timer)

        max_linked = int(choose(nel, (nel + LMS) / 2))
        max_per_proc = 2 * (max_linked / nProcessors) + 1
        max_spawned = max_per_proc * nProcessors

        ssq = calc_s_squared_multi_worker(max_per_proc, max_spawned)

        call halt_timer(s2_timer)

    end function

    function calc_s_squared_multi_worker(max_per_proc, max_spawned) &
        result(ssq)

        integer :: i, j, k, orb2, orbtmp, pos, ierr
        integer :: nI(nel), nJ(nel), proc
        integer(n_int), pointer :: detcurr(:)
        integer(n_int) :: splus(0:NIfTot), sminus(0:NIfTot)
        logical :: running, any_running
        real(dp), dimension(inum_runs) :: ssq, tmp
        integer :: max_per_proc, max_spawned, run
        real(dp) :: sgn1(lenof_sign), sgn2(lenof_sign)

        ! Could we pre-initialise all of these data structures
        !integer :: max_linked = int(choose(nel, (nel + LMS)/2))
        integer(n_int) :: det_list(0:NIfTot, max_spawned)
        integer(n_int) :: recv_dets(0:NIfTot, max_spawned)
        integer :: proc_pos(nProcessors), proc_pos_init(nProcessors)
        integer :: send_count(nProcessors), recv_count(nProcessors)
        integer(MPIArg) :: send_data(nProcessors), recv_data(nProcessors)
        integer(MPIArg) :: send_off(nProcessors), recv_off(nProcessors)

        running = (TotWalkers > 0)
        any_running = .true.
        j = 1
        ssq = 0
        forall (i=1:nProcessors) proc_pos_init(i) = (i - 1) * max_per_proc + 1

        do while (any_running)

            ! Clear transmission lists
            proc_pos = proc_pos_init

            if (running) then

                ! Generate items, add to list (and use the sgn of initial
                ! walker, so we send it to the target processor)
                call ilut_nifd_pointer_assign(detcurr, CurrentDets(0:NIfTot, j))
                call decode_bit_det(nI, detcurr)

                do i = 1, nel
                    orbtmp = get_alpha(nI(i))
                    if (is_beta(nI(i)) &
                        .and. IsNotOcc(detcurr, orbtmp)) then
                        splus = detcurr
                        clr_orb(splus, nI(i))
                        set_orb(splus, orbtmp)

                        do k = 1, nel
                            orb2 = nI(k)
                            if (k == i) orb2 = get_alpha(orb2)

                            orbtmp = get_beta(orb2)
                            if (is_alpha(orb2) &
                                .and. IsNotOcc(splus, orbtmp)) then
                                sminus = splus
                                clr_orb(sminus, orb2)
                                set_orb(sminus, orbtmp)

                                ! Store this det (n.b. contains original sgn)
                                call decode_bit_det(nJ, sminus)
                                proc = DetermineDetNode(nel, nJ, 0) + 1
                                det_list(:, proc_pos(proc)) = sminus
                                proc_pos(proc) = proc_pos(proc) + 1
                            end if
                        end do
                    end if
                end do

                ! Walk through the list. Stop when we get to the end.
                j = j + 1
                if (j > TotWalkers) running = .false.

            end if

            ! How many elements are there in each list?
            send_count = proc_pos - proc_pos_init
            if (any(send_count > max_per_proc)) &
                send_count = max_per_proc + 1

            call MPIAlltoAll(send_count, 1, recv_count, 1, ierr)

            send_off = int((proc_pos_init - 1) * (NIfTot + 1), MPIArg)
            recv_off(1) = 0
            do i = 2, nProcessors
                recv_off(i) = recv_off(i - 1) + int(recv_count(i - 1), MPIArg)
            end do
            recv_off = recv_off * int(NIfTot + 1, MPIArg)
            send_data = int(send_count * (NIfTot + 1), MPIArg)
            recv_data = int(recv_count * (NIfTot + 1), MPIArg)

            call MPIAlltoAllv(det_list, send_data, send_off, &
                              recv_dets, recv_data, recv_off, ierr)

            ! Find the det in list, and sum in its term.
            do i = 1, sum(recv_count)

                ! The sign of the source particle
                call extract_sign(recv_dets(:, i), sgn1)

                ! And the generated, connected particle
                pos = binary_search(CurrentDets(:, 1:TotWalkers), &
                                    recv_dets(:, i), NIfD + 1)
                if (pos > 0) then
                    call extract_sign(CurrentDets(:, pos), sgn2)
                    ssq = ssq + (sgn1(1) * sgn2(1))
                end if

            end do

            ! Is there anything left to do on any process?
            call MPIAllLORLogical(running, any_running)

        end do

        call MPISum(ssq, tmp)
        ssq = tmp

        do run = 1, inum_runs
            if (abs(all_norm_psi_squared(run)) < 1.0e-10_dp) then
                ssq(run) = 0.0_dp
            else
                ssq(run) = ssq(run) / all_norm_psi_squared(run)

                ! TODO: n.b. This is a hack. LMS appears to contain -2Ms of the
                !            system. I am somewhat astounded I haven't noticed
                !            this before...
                ssq(run) = ssq(run) &
                           + real(calculated_ms * (calculated_ms + 2), dp) / 4
            end if
        end do

    end function

    function calc_s_squared_star(only_init) result(ssq)

        real(dp), dimension(inum_runs) :: ssq
        integer, parameter :: max_per_proc = 1000
        integer(n_int) :: recv_dets(0:NIfTot, max_per_proc)
        integer :: proc_dets, start_pos, nsend, i, p
        integer :: bcast_tmp(2), run
        real(dp) :: sgn_tmp(lenof_sign)
        type(timer), save :: s2_timer, s2_timer_init
        real(dp), dimension(inum_runs) :: ssq_sum, psi_squared
        real(dp), dimension(inum_runs):: All_ssq_sum, All_psi_squared
        logical, intent(in) :: only_init

        s2_timer%timer_name = 'S^2 star'
        s2_timer_init%timer_name = 'S^2 init'
        if (only_init) then
            call set_timer(s2_timer_init)
        else
            call set_timer(s2_timer)
        end if

        ssq_sum = 0.0_dp
        psi_squared = 0.0_dp
        do p = 0, nProcessors - 1

            ! How many dets are on processor p
            proc_dets = int(TotWalkers, sizeof_int)
            call MPIBcast(proc_dets, iProcIndex == p)

            ! Send the dets around bit by bit
            start_pos = 1
            do while (start_pos <= proc_dets)

                if (tTruncInitiator .and. only_init) then
                    ! Loop over walkers and only add initiators to bcast list
                    nsend = 0
                    if (p == iProcIndex) then
                        do i = start_pos, int(TotWalkers, sizeof_int)
                            ! Break up the list into correctly sized chunks
                            if (nsend == max_per_proc) exit

                            if (any_run_is_initiator(CurrentDets(:, i))) then
                                nsend = nsend + 1
                                recv_dets(:, nsend) = CurrentDets(:, i)

                                ! If we are using initiators only, keep track
                                ! of the overall magnitude of the init-only
                                ! wavefunction.
                                call extract_sign(CurrentDets(:, i), sgn_tmp)
                                psi_squared = psi_squared + sum(sgn_tmp**2)
                            end if
                        end do
                        start_pos = i
                    end if
                    bcast_tmp = (/nsend, start_pos/)
                    call MPIBcast(bcast_tmp, p == iProcIndex)
                    nsend = bcast_tmp(1)
                    start_pos = bcast_tmp(2)

                else
                    ! How many dets to send in this batch?
                    if (start_pos + max_per_proc - 1 > proc_dets) then
                        nsend = proc_dets - start_pos + 1
                    else
                        nsend = max_per_proc
                    end if

                    ! Update list of received dets
                    if (p == iProcIndex) recv_dets(:, 1:nsend) = &
                        CurrentDets(:, start_pos:start_pos + nsend - 1)

                    ! Increment position
                    start_pos = start_pos + nsend

                end if

                ! Broadcast to all processors
                call MPIBcast(recv_dets(:, 1:nsend), iProcIndex == p)

                ! All processors loop over these dets, and calculate their
                ! contribution to S^2
                do i = 1, nsend
                    if (.not. TestClosedShellDet(recv_dets(:, i))) then
                        ! If nopen == 2, and tHPHF, then this can be
                        ! optimised further
                        ssq_sum = ssq_sum + ssquared_contrib(recv_dets(:, i), &
                                                             only_init)
                    end if
                end do

            end do

        end do ! Loop over processors

        ! Sum all of the s squared terms
        if (tTruncInitiator .and. only_init) then
            call MPISum(psi_squared, All_psi_squared)
            psi_squared = All_psi_squared
        else
            psi_squared = all_norm_psi_squared
        end if
        call MPISum(ssq_sum, All_ssq_sum)
        ssq_sum = All_ssq_sum

        do run = 1, inum_runs
            if (abs(psi_squared(run)) < 1.0e-10_dp) then
                ssq(run) = 0.0_dp
            else
                ssq(run) = real(ssq_sum(run), dp) / psi_squared(run)

                ! TODO: n.b. This is a hack. LMS appears to contain -2Ms of the
                !            system. I am somewhat astounded I haven't noticed
                !            this before...
                ssq(run) = ssq(run) &
                           + real(calculated_ms * (calculated_ms + 2), dp) / 4.0
            end if
        end do

        if (only_init) then
            call halt_timer(s2_timer_init)
        else
            call halt_timer(s2_timer)
        end if

    end function

    function ssquared_contrib(ilut, only_init_, n_opt, ilut_list_opt) result(ssq)

        ! Calculate the contribution to s-squared from the determinant
        ! provided (from walkers on this processor).
        !
        ! This applies the operator S-S+, returning the result:
        !
        !      <Psi(iProcIndex) | S-S+ | D_i>

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        logical, intent(in), optional :: only_init_
        integer, intent(in), optional :: n_opt
        integer(n_int), intent(in), optional :: ilut_list_opt(0:, :)
        real(dp) :: ssq
#ifdef DEBUG_
        character(*), parameter :: this_routine = "ssquared_contrib"
#endif
        integer(n_int) :: splus(0:NIfD), sminus(0:NIfD)
        integer(n_int) :: ilut_srch(0:NIfD), ilut_sym(0:NIfD)
        real(dp) :: sgn(lenof_sign), sgn2(lenof_sign), sgn_hphf
        integer :: flg, nI(nel), j, k, orb2, pos, orb_tmp
        logical :: only_init, inc
        integer :: n_states
        integer(n_int), allocatable :: ilut_list(:, :)

        if (present(only_init_)) then
            only_init = only_init_
        else
            only_init = .false.
        end if

        ! make this routine more flexible and usable not only for CurrentDets
        if (present(n_opt)) then
            ASSERT(present(ilut_list_opt))
            ASSERT(size(ilut_list_opt, 2) == n_opt)

            n_states = n_opt
            allocate(ilut_list(0:niftot, n_states), source=ilut_list_opt(0:niftot, 1:n_opt))

        else
            n_states = int(TotWalkers)
            allocate(ilut_list(0:niftot, TotWalkers), source=CurrentDets(0:niftot, 1:TotWalkers))

        end if

        ! Extract details of determinant
        call extract_bit_rep(ilut, nI, sgn, flg)

        ssq = 0.0_dp
        do j = 1, nel
            if (is_beta(nI(j)) &
                .and. IsNotOcc(ilut, get_alpha(nI(j)))) then
                splus = ilut(0:NIfD)
                orb_tmp = get_alpha(nI(j))
                clr_orb(splus, nI(j))
                set_orb(splus, orb_tmp)

                do k = 1, nel
                    orb2 = nI(k)
                    if (k == j) orb2 = get_alpha(orb2)

                    orb_tmp = get_beta(orb2)
                    if (is_alpha(orb2) &
                        .and. IsNotOcc(splus, orb_tmp)) then
                        sminus = splus
                        clr_orb(sminus, orb2)
                        set_orb(sminus, orb_tmp)

                        ! Adjust for the sign of the paired det in HPHF.
                        sgn_hphf = 1.0
                        if (tHPHF) then
                            if (IsAllowedHPHF(sminus, ilut_sym)) then
                                ilut_srch = sminus
                            else
                                ilut_srch = ilut_sym
                                sgn_hphf = hphf_sign(ilut_srch)
                            end if
                        else
                            ilut_srch = sminus
                        end if

                        ! --> sminus is an allowed result of applying S-S+
                        pos = binary_search(ilut_list(:, 1:n_states), &
                                            ilut_srch, NIfD + 1)
                        if (pos > 0) then

                            ! If we are looking for the spin of the initiator
                            ! only wavefunction, we need to ensure that we
                            ! are projecting onto an initiator...
                            inc = .true.
                            if (tTruncInitiator .and. only_init) then
                                inc = (any_run_is_initiator(ilut_list(:, pos)))
                            end if

                            call extract_sign(ilut_list(:, pos), sgn2)
                            ssq = ssq + sgn(1) * sgn2(1) * sgn_hphf

                        end if
                    end if
                end do
            end if
        end do

    end function

    subroutine init_hist_excit_tofrom()

        integer :: ierr
        character(*), parameter :: this_routine = 'init_hist_excit_tofrom'

        ! Initialise storage for these excitaitons
        allocate(hist_excit_tofrom(0:nel, 0:nel), stat=ierr)
        @:log_alloc(hist_excit_tofrom, tag_hist_excit, ierr)

        ! Zero everything
        hist_excit_tofrom = 0

        if (iProcIndex == root) then

            ! Open an output file
            excit_tofrom_unit = get_free_unit()
            open (excit_tofrom_unit, file='spawns_tofrom_excit', &
                  status='replace')

            ! Write a header in the output file
            write (excit_tofrom_unit, '("# Number of particles spawned between &
                                       &excitation levels from the Hartree--&
                                       &Fock on this update cycle")')
            write (excit_tofrom_unit, '("# iter,  0-->0,  0-->1,  0-->2, ..., &
                                      &1-->0, ...")')

        end if

    end subroutine

    subroutine clean_hist_excit_tofrom()

        character(*), parameter :: this_routine = 'clean_hist_excit_tofrom'

        ! Clean up the stored data.
        deallocate(hist_excit_tofrom)
        log_dealloc(tag_hist_excit)

        ! Close the output file
        if (iProcIndex == root) &
            close(excit_tofrom_unit)

    end subroutine

    subroutine add_hist_excit_tofrom(iluti, ilutj, child)

        integer(n_int), intent(in) :: ilutI(0:NIfTot), ilutJ(0:NIfTot)
        real(dp), intent(in) :: child(lenof_sign)
        real(dp) :: abschild
        integer :: exlevelI, exlevelJ

        character(*), parameter :: this_routine = "add_hist_excit_tofrom"
        ! We want a total count of the particle weight formed.
        abschild = sum(abs(child))

        ! Get the excitation levels of the source and target
        ASSERT(.not. tGUGA)
        exlevelI = FindBitExcitLevel(ilutRef(:, 1), ilutI, t_hphf_ic=.true.)
        exlevelJ = FindBitExcitLevel(ilutRef(:, 1), ilutJ, t_hphf_ic=.true.)

        ! And store it!
        hist_excit_tofrom(exlevelI, exlevelJ) = &
            hist_excit_tofrom(exlevelI, exlevelJ) + abschild

    end subroutine

    subroutine write_zero_hist_excit_tofrom()

        real(dp) :: all_hist(0:nel, 0:nel)
        integer :: i, j

        ! Gather the accumulator data to this process, and then zero the data
        ! for the next update cycle.
        call MPISum(hist_excit_tofrom, all_hist)
        hist_excit_tofrom = 0

        if (iProcIndex == root) then

            ! Output the current iteration
            write(excit_tofrom_unit, '(i12)', advance='no') iter

            ! And output the accumulated data
            do i = 0, nel
                do j = 0, nel
                    write(excit_tofrom_unit, '(f16.5)', advance='no') &
                        all_hist(i, j)
                end do
            end do
            write(excit_tofrom_unit, *)
        end if

    end subroutine

end module
