#include "macros.h"

module fcimc_pointed_fns

    use SystemData, only: nel, tGUGA, tGen_nosym_guga, &
                          tGen_sym_guga_mol, t_consider_diff_bias, nSpatOrbs, thub, &
                          tUEG, nBasis, tgen_guga_crude, &
                          t_k_space_hubbard, t_new_real_space_hubbard, &
                          t_trans_corr_2body, t_trans_corr_hop, &
                          t_precond_hub, uhub

    use LoggingData, only: tHistExcitToFrom, FciMCDebug

    use CalcData, only: RealSpawnCutoff, tRealSpawnCutoff, tAllRealCoeff, &
                        RealCoeffExcitThresh, tau, DiagSft, &
                        tRealCoeffByExcitLevel, InitiatorWalkNo, &
                        t_fill_frequency_hists, t_truncate_spawns, n_truncate_spawns, &
                        t_matele_cutoff, matele_cutoff, tEN2Truncated, &
                        t_hist_tau_search_option, &
                        tTruncInitiator, tSkipRef, t_truncate_unocc, t_consider_par_bias, &
                        tAdaptiveShift, LAS_Sigma, LAS_F1, LAS_F2, &
                        AAS_Thresh, AAS_Expo, AAS_Cut, &
                        tPrecond, AAS_Const, EAS_Scale, ShiftOffset, tAS_Offset
    use DetCalcData, only: FciDetIndex, det
    use procedure_pointers, only: get_spawn_helement, shiftFactorFunction
    use fcimc_helper, only: CheckAllowedTruncSpawn

    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet, count_open_orbs

    use load_balance, only: scaleFunction

    use bit_rep_data, only: NIfTot, test_flag

    use tau_search, only: log_death_magnitude, fill_frequency_histogram_nosym_diff, &
                          fill_frequency_histogram_nosym_nodiff, log_spawn_magnitude

    use bit_reps, only: get_initiator_flag, get_initiator_flag_by_run

    use rdm_general, only: calc_rdmbiasfac
    use hist, only: add_hist_excit_tofrom
    use searching, only: BinSearchParts2
    use UMatCache, only: UMatInd, gtID
    use util_mod
    use FciMCData
    use constants

    use bit_reps, only: nifguga

    use guga_matrixElements, only: calcDiagMatEleGUGA_nI

#ifdef DEBUG_
    use guga_bitRepOps, only: convert_ilut_toGUGA, write_det_guga
    use guga_excitations, only: print_excitInfo, global_excitInfo
#endif

    use real_time_data, only: runge_kutta_step, t_real_time_fciqmc

    use tau_search_hist, only: fill_frequency_histogram_4ind, &
                               fill_frequency_histogram_sd, &
                               fill_frequency_histogram

    use excit_gen_5, only: pgen_select_a_orb

    use cepa_shifts, only: t_cepa_shift, cepa_shift

    use hphf_integrals, only: hphf_diag_helement

    use Determinants, only: get_helement
    use global_det_data, only: get_tot_spawns, get_acc_spawns

    implicit none

contains

    function attempt_create_trunc_spawn(DetCurr, iLutCurr, RealwSign, nJ, &
                                        iLutnJ, prob, HElGen, ic, ex, tparity, walkExcitLevel, part_type, &
                                        AvSignCurr, AvExPerWalker, RDMBiasFacCurr, precond_fac) result(child)

        integer, intent(in) :: DetCurr(nel), nJ(nel), part_type
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2, ic), walkExcitLevel
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        real(dp), dimension(lenof_sign) :: child
        real(dp), dimension(lenof_sign), intent(in) :: AvSignCurr
        real(dp), intent(in) :: AvExPerWalker
        real(dp), intent(out) :: RDMBiasFacCurr
        real(dp), intent(in) :: precond_fac
        logical :: tAllowForEN2Calc
        HElement_t(dp), intent(inout) :: HElGen

        ! If tEN2Truncated is true, then we want to allow the otherwise
        ! truncated spawning to allow an EN2 correction to be calculated,
        ! and then we will cancel it later in Annihilation, at the point
        ! the correction is added in. However, this is currently only
        ! implemented in the full, non-initiator scheme. So if initiator
        ! is on, then ignore this.
        tAllowForEN2Calc = tEN2Truncated .and. (.not. tTruncInitiator)

        if (.not. tAllowForEN2Calc) then
            if (CheckAllowedTruncSpawn(walkExcitLevel, nJ, iLutnJ, IC)) then
                child = attempt_create_normal(DetCurr, iLutCurr, RealwSign, &
                                              nJ, iLutnJ, prob, HElGen, ic, ex, tParity, walkExcitLevel, &
                                              part_type, AvSignCurr, AvExPerWalker, RDMBiasFacCurr, precond_fac)
            else
                child = 0
            end if
        else
            child = attempt_create_normal(DetCurr, iLutCurr, RealwSign, nJ, &
                                          iLutnJ, prob, HElGen, ic, ex, tParity, walkExcitLevel, part_type, &
                                          AvSignCurr, AvExPerWalker, RDMBiasFacCurr, precond_fac)
        end if
    end function

!Decide whether to spawn a particle at nJ from DetCurr. (bit strings iLutnJ and iLutCurr respectively).
!  ic and ex specify the excitation of nJ from DetCurr, along with the sign change tParity.
!  part_type:           Is the parent real (1) or imaginary (2)
!  wSign:               wSign gives the sign of the particle we are trying to spawn from
!                          if part_type is 1, then it will only use wsign(1)
!                                          2,                       wsign(2)
!                       Only the sign, not magnitude is used.
!  prob:                prob is the generation probability of the excitation in order to unbias.
!                       The probability of spawning is divided by prob to do this.
!  HElGen:              If the matrix element has already been calculated, it is sent in here.
!  get_spawn_helement:  A function pointer for looking up or calculating the relevant matrix element.
!  walkExcitLevel:      Is Unused
!
!  child:      A lenof_sign array containing the particles spawned.
    function att_create_trunc_spawn_enc(DetCurr, &
                                        iLutCurr, RealwSign, nJ, iLutnJ, prob, HElGen, &
                                        ic, ex, tparity, walkExcitLevel, part_type, &
                                        AvSignCurr, AvExPerWalker, RDMBiasFacCurr, precond_fac) result(child)

        integer, intent(in) :: DetCurr(nel), nJ(nel), part_type
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2, ic), walkExcitLevel
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        real(dp), dimension(lenof_sign) :: child
        real(dp), dimension(lenof_sign), intent(in) :: AvSignCurr
        real(dp), intent(in) :: AvExPerWalker
        real(dp), intent(out) :: RDMBiasFacCurr
        real(dp), intent(in) :: precond_fac
        logical :: tAllowForEN2Calc
        ! make sure that HElgen is assigned on return
        HElement_t(dp), intent(inout) :: HElGen

        call EncodeBitDet(nJ, iLutnJ)

        ! If tEN2Truncated is true, then we want to allow the otherwise
        ! truncated spawning to allow an EN2 correction to be calculated,
        ! and then we will cancel it later in Annihilation, at the point
        ! the correction is added in. However, this is currently only
        ! implemented in the full, non-initiator scheme. So if initiator
        ! is on, then ignore this.
        tAllowForEN2Calc = tEN2Truncated .and. (.not. tTruncInitiator)

        if (.not. tAllowForEN2Calc) then
            if (CheckAllowedTruncSpawn(walkExcitLevel, nJ, iLutnJ, IC)) then
                child = attempt_create_normal(DetCurr, iLutCurr, RealwSign, &
                                              nJ, iLutnJ, prob, HElGen, ic, ex, tParity, walkExcitLevel, &
                                              part_type, AvSignCurr, AvExPerWalker, RDMBiasFacCurr, precond_fac)
            else
                child = 0
            end if
        else
            child = attempt_create_normal(DetCurr, iLutCurr, RealwSign, nJ, &
                                          iLutnJ, prob, HElGen, ic, ex, tParity, walkExcitLevel, part_type, &
                                          AvSignCurr, AvExPerWalker, RDMBiasFacCurr, precond_fac)
        end if
    end function

    function attempt_create_normal(DetCurr, iLutCurr, RealwSign, nJ, iLutnJ, &
                                   prob, HElGen, ic, ex, tParity, walkExcitLevel, part_type, AvSignCurr, &
                                   AvExPerWalker, RDMBiasFacCurr, precond_fac) result(child)

        use orb_idx_mod, only: SpinOrbIdx_t
        use gasci, only: operator(.contains.), GAS_specification
        use SystemData, only: tGAS

        integer, intent(in) :: DetCurr(nel), nJ(nel)
        integer, intent(in) :: part_type    ! odd = Real parent particle, even = Imag parent particle
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ(0:niftot)
        integer, intent(in) :: ic, ex(2, ic), walkExcitLevel
        real(dp), dimension(lenof_sign), intent(in) :: RealwSign
        logical, intent(in) :: tParity
        real(dp), intent(inout) :: prob
        real(dp), dimension(lenof_sign) :: child
        real(dp), dimension(lenof_sign), intent(in) :: AvSignCurr
        real(dp), intent(in) :: AvExPerWalker
        real(dp), intent(out) :: RDMBiasFacCurr
        real(dp), intent(in) :: precond_fac
        HElement_t(dp), intent(inout) :: HElGen
        character(*), parameter :: this_routine = 'attempt_create_normal'

        real(dp) :: rat, r, walkerweight, pSpawn, nSpawn, MatEl, p_spawn_rdmfac
        integer :: extracreate, tgt_cpt, component, i, iUnused
        integer :: TargetExcitLevel
        logical :: tRealSpawning
        HElement_t(dp) :: rh_used, Eii_curr
#ifdef DEBUG_
        integer :: nOpen
        integer(n_int) :: ilutTmpI(0:nifguga), ilutTmpJ(0:nifguga)
#endif
        logical :: t_par
        real(dp) :: temp_prob, pgen_a, dummy_arr(nBasis), cum_sum
        integer :: ispn, n_double

        integer :: temp_ex(2, ic)

        unused_var(AvSignCurr)
        unused_var(WalkExcitLevel)

        ! Just in case
        child = 0.0_dp

        prob = prob * AvExPerWalker
        ! In the case of using HPHF, and when tGenMatHEl is on, the matrix
        ! element is calculated at the time of the excitation generation,
        ! and returned in HElGen. In this case, get_spawn_helement simply
        ! returns HElGen, rather than recomputing the matrix element.

        temp_ex(1, :) = ex(2, :)
        temp_ex(2, :) = ex(1, :)

        ! We actually want to calculate Hji - take the complex conjugate,
        ! rather than swap around DetCurr and nJ.
        rh_used = get_spawn_helement(nJ, DetCurr, ilutnJ, iLutCurr, ic, temp_ex, &
                                     tParity, HElGen)

        ! assign the matrix element
        HElGen = abs(rh_used)

        ! essentially here i have all the information for my frequency
        ! analysis of the H_ij/pgen ratio so call the routine here
        ! but i have to remember to keep it parallel! so dont forget to
        ! sum up all the contributions from different cores!
        ! and divide prob by AvMCExcits again to get correct pgen!
        if (t_fill_frequency_hists) then
            ! use specific ones for different types of excitation gens
            if (tHUB .or. tUEG .or. &
                (t_new_real_space_hubbard .and. .not. t_trans_corr_hop) .or. &
                (t_k_space_hubbard .and. .not. t_trans_corr_2body)) then
                call fill_frequency_histogram(abs(rh_used / precond_fac), prob)

            else
                if (t_consider_par_bias) then
                    ! determine if excitation was parallel or anti-parallel
                    ! ex(1,1) and ex(1,2) are the electrons
                    if (ic == 2) then
                        t_par = (is_beta(ex(1, 1)) .eqv. is_beta(ex(1, 2)))
                    else
                        t_par = .false.
                    end if

                    call fill_frequency_histogram_4ind(abs(rh_used / precond_fac), prob, &
                                                       ic, t_par, ex)

                else if (tGen_nosym_guga) then
                    ! have to also check if diff bias is considered
                    if (t_consider_diff_bias) then
                        call fill_frequency_histogram_nosym_diff(abs(rh_used / precond_fac), prob, &
                                                                 ic, ex(1, 1), ex(1, 2))
                    else
                        call fill_frequency_histogram_nosym_nodiff(abs(rh_used / precond_fac), &
                                                                   prob, ic, ex(1, 1))
                    end if

                else
                    ! for any other excitation generator just use one histogram
                    ! for all the excitations..
                    call fill_frequency_histogram_sd(abs(rh_used / precond_fac), prob, ic)

                end if
            end if
        end if

        ! If each walker does not have exactly one spawning attempt
        ! (if AvExPerWalker /= 1.0_dp) then the probability of an excitation
        ! having been chosen, prob, must be altered accordingly.
        ! the simple preconditioner for hubbard

        ! i am not so sure about hongjuns change here...
        ! is this for GUGA only?? Ask him

        if (tGUGA .and. t_precond_hub) then
            Eii_curr = calcDiagMatEleGUGA_nI(nJ)
        end if

        tRealSpawning = .false.
        if (tAllRealCoeff) then
            tRealSpawning = .true.
        else if (tRealCoeffByExcitLevel) then
            if (tGUGA) call stop_all(this_routine, &
                                     "excit level does not work with GUGA here...")

            TargetExcitLevel = FindBitExcitLevel(iLutRef(:, 1), iLutnJ)

            if (TargetExcitLevel <= RealCoeffExcitThresh) &
                tRealSpawning = .true.
        end if

        ! Spawn to real and imaginary particles. Note that spawning from
        ! imaginary parent particles has slightly different rules:
        !       - Attempt to spawn REAL walkers with prob +AIMAG(Hij)/P
        !       - Attempt to spawn IMAG walkers with prob -REAL(Hij)/P

#if !defined(CMPLX_) && (defined(PROG_NUMRUNS_) || defined(DOUBLERUN_))
        child = 0.0_dp
        tgt_cpt = part_type
        walkerweight = sign(1.0_dp, RealwSign(part_type))
        matEl = real(rh_used, dp)
        if (t_matele_cutoff) then
            if (abs(matEl) < matele_cutoff) matel = 0.0_dp
        end if
#else
        do tgt_cpt = 1, (lenof_sign / inum_runs)

            ! Real, single run:    inum_runs=1, lenof_sign=1 --> 1 loop
            ! Real, double run:    inum_runs=2, lenof_sign=1 --> 1 loop
            ! Complex, single run: inum_runs=1, lenof_sign=2 --> 2 loops
            ! Complex, double run: inum_runs=2, lenof_sign=4 --> 2 loops
            ! Complex, multiple run: inum_runs=m, lenof_sign=2*m --> 2 loops

            ! For spawning from imaginary particles, we cross-match the
            ! real/imaginary matrix-elements/target-particles.

#if defined(CMPLX_) && (defined(PROG_NUMRUNS_) || defined(DOUBLERUN_))
            component = part_type + tgt_cpt - 1
            if (.not. btest(part_type, 0)) then
                ! even part_type => imag replica =>  map 4->3,4 ; 6->5,6 etc.
                component = part_type - tgt_cpt + 1
            end if
#else
            component = tgt_cpt
            if ((part_type == 2) .and. (inum_runs == 1)) component = 3 - tgt_cpt
#endif

            ! Get the correct part of the matrix element
            walkerweight = sign(1.0_dp, RealwSign(part_type))
            if (btest(component, 0)) then
                ! real component
                MatEl = real(rh_used, dp)
                if (t_matele_cutoff) then
                    if (abs(MatEl) < matele_cutoff) MatEl = 0.0_dp
                end if
            else
#ifdef CMPLX_
                MatEl = real(aimag(rh_used), dp)
                if (t_matele_cutoff) then
                    if (abs(MatEl) < matele_cutoff) MatEl = 0.0_dp
                end if
                ! n.b. In this case, spawning is of opposite sign.
                if (.not. btest(part_type, 0)) then
                    ! imaginary parent -> imaginary child
                    walkerweight = -walkerweight
                end if
#endif
            end if
#endif

!             ASSERT(prob /= 0.0_dp)
            nSpawn = -tau * MatEl * walkerweight / prob

#ifdef DEBUG_
            ! so.. does a |nSpawn| > 1 in my new tau-search implitly mean
            ! that the tau is to small for this kind of exciation?
            ! i guess so yeah.. so do it really brute force to start with
            if (t_truncate_spawns .and. abs(nSpawn) > n_truncate_spawns) then
                ! in debug mode i should output some additional information
                ! to analyze the type of excitations and how many open-orbitals
                ! etc. are used
                if (abs(nSpawn) > 10.0_dp) then
                    if (tGUGA) then
                        write(iout, *) "=================================================="
                        call convert_ilut_toGUGA(iLutCurr, ilutTmpI)
                        call convert_ilut_toGUGA(ilutnj, ilutTmpJ)
                        write(iout, *) "nSpawn > n_truncate_spawns!", nSpawn
                        write(iout, *) "limit the number of spawned walkers to: ", n_truncate_spawns
                        write(iout, *) "for spawn from determinant: "
                        call write_det_guga(6, ilutTmpI)
                        write(iout, *) "to: "
                        call write_det_guga(iout, ilutTmpJ)
                        nOpen = count_open_orbs(iLutCurr)
                        write(iout, *) "# of openshell orbitals: ", nOpen, count_open_orbs(ilutnj)
                        write(iout, *) "open/spatial: ", nOpen / real(nSpatOrbs, dp)
                        write(iout, *) "(t |H_ij|/pgen) / #open ratio: ", abs(nSpawn) / real(nOpen, dp)
                        write(iout, *) " H_ij, pgen: ", MatEl, prob
                        write(iout, *) "=================================================="
                        ! excitation type would be cool too.. but how do i get it
                        ! to here?? do i still have global_excitInfo??
                        call print_excitInfo(global_excitInfo)
                        call neci_flush(iout)
                    end if
                end if
            end if
#endif

            if (tGUGA .and. t_precond_hub .and. abs(Eii_curr) > 10.0_dp) then
                nSpawn = nSpawn / (1.0_dp + dble(Eii_curr))
            end if

            ! [Werner Dobrautz 4.4.2017:]
            ! apply the spawn truncation, when using histogramming tau-search
            if ((t_truncate_spawns .and. .not. t_truncate_unocc) .and. abs(nspawn) > &
                n_truncate_spawns .and. .not. tEScaleWalkers) then
                ! does not work with scaled walkers, as the scaling factor is not
                ! computed here for performance reasons (it was a huge performance bottleneck)
                ! TODO: add some additional output if this event happens
#ifdef DEBUG_
                write(iout, *) "Truncating spawn magnitude from: ", abs(nspawn), " to ", n_truncate_spawns
#endif
                truncatedWeight = truncatedWeight + abs(nSpawn) - n_truncate_spawns
                nSpawn = sign(n_truncate_spawns, nspawn)

            end if

            ! n.b. if we ever end up with |walkerweight| /= 1, then this
            !      will need to ffed further through.
            if (tSearchTau .and. (.not. tFillingStochRDMonFly)) then
                ! in the back-spawning i have to adapt the probabilites
                ! back, to be sure the time-step covers the changed
                ! non-initiators spawns!

                call log_spawn_magnitude(ic, ex, matel / precond_fac, prob)
            end if

            ! Keep track of the biggest spawn this cycle
            max_cyc_spawn = max(abs(nSpawn), max_cyc_spawn)

            if (tRealSpawning) then
                ! Continuous spawning. Add in acceptance probabilities.
                if (tRealSpawnCutoff .and. &
                    abs(nSpawn) < RealSpawnCutoff) then
                    p_spawn_rdmfac = abs(nSpawn) / RealSpawnCutoff
                    nSpawn = RealSpawnCutoff &
                             * stochastic_round(nSpawn / RealSpawnCutoff)
                else
                    p_spawn_rdmfac = 1.0_dp !The acceptance probability of some kind of child was equal to 1
                end if
            else
                if (abs(nSpawn) >= 1.0) then
                    p_spawn_rdmfac = 1.0_dp !We were certain to create a child here.
                    ! This is the special case whereby if P_spawn(j | i) > 1,
                    ! then we will definitely spawn from i->j.
                    ! I.e. the pair Di,Dj will definitely be in the SpawnedParts list.
                    ! We don't care about multiple spawns - if it's in the list, an RDM contribution will result
                    ! regardless of the number spawned - so if P_spawn(j | i) > 1, we treat it as = 1.
                else
                    p_spawn_rdmfac = abs(nSpawn)
                end if

                ! How many children should we spawn?

                ! And round this to an integer in the usual way
                ! HACK: To use the same number of random numbers for the tests.
                nSpawn = real(stochastic_round(nSpawn), dp)

            end if
            ! And create the parcticles
#ifdef CMPLX_
            child((part_type_to_run(part_type) - 1) * 2 + tgt_cpt) = nSpawn
#else
            child(tgt_cpt) = nSpawn
#endif

#if defined(CMPLX_) || !defined(PROG_NUMRUNS_) && !defined(DOUBLERUN_)
        end do
#endif

        if (tFillingStochRDMonFly) then
            if (.not. near_zero(child(part_type))) then
                !Only add in contributions for spawning events within population 1
                !(Otherwise it becomes tricky in annihilation as spawnedparents doesn't tell you which population
                !the event came from at present)
                call calc_rdmbiasfac(p_spawn_rdmfac, prob, realwSign(part_type), RDMBiasFacCurr)
            else
                RDMBiasFacCurr = 0.0_dp
            end if
        else
            ! Not filling the RDM stochastically, bias is zero.
            RDMBiasFacCurr = 0.0_dp
        end if

    end function

    !
    ! This is a null routine for encoding spawned sites
    ! --> DOES NOTHING!!!
    subroutine null_encode_child(ilutI, ilutJ, ic, ex)
        implicit none
        integer(kind=n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(in) :: ic
        integer, intent(in) :: ex(2, ic)
        integer(kind=n_int), intent(inout) :: ilutj(0:niftot)

        unused_var(ilutI); unused_var(ilutJ); unused_var(ic); unused_var(ex)

    end subroutine

    subroutine new_child_stats_hist_hamil(iter_data, iLutI, nJ, iLutJ, ic, &
                                          walkExLevel, child, parent_flags, &
                                          part_type)
        ! Based on old AddHistHamilEl. Histograms the hamiltonian matrix, and
        ! then calls the normal statistics routine.

        integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer, intent(in) :: ic, walkExLevel, parent_flags, nJ(nel)
        integer, intent(in) :: part_type
        real(dp), dimension(lenof_sign), intent(in) :: child
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = 'new_child_stats_hist_hamil'
        integer :: partInd, partIndChild, childExLevel
        logical :: tSuccess

        if (walkExLevel == nel) then
            call BinSearchParts2(iLutI, FCIDetIndex(walkExLevel), Det, &
                                 PartInd, tSuccess)
        else
            call BinSearchParts2(iLutI, FCIDetIndex(walkExLevel), &
                                 FciDetIndex(walkExLevel + 1) - 1, partInd, &
                                 tSuccess)
        end if

        if (.not. tSuccess) &
            call stop_all(this_routine, 'Cannot find determinant nI in list')

        childExLevel = FindBitExcitLevel(iLutHF, iLutJ, nel)
        if (tGUGA) call stop_all(this_routine, &
                                 "excit level does not work with GUGA here...")

        if (childExLevel == nel) then
            call BinSearchParts2(iLutJ, FCIDetIndex(childExLevel), Det, &
                                 partIndChild, tSuccess)
        else if (childExLevel == 0) then
            partIndChild = 1
            tSuccess = .true.
        else
            call BinSearchParts2(iLutJ, FCIDetIndex(childExLevel), &
                                 FciDetIndex(childExLevel + 1) - 1, &
                                 partIndChild, tSuccess)
        end if

        histHamil(partIndChild, partInd) = &
            histHamil(partIndChild, partInd) + (1.0_dp * child(1))
        histHamil(partInd, partIndChild) = &
            histHamil(partInd, partIndChild) + (1.0_dp * child(1))
        avHistHamil(partIndChild, partInd) = &
            avHistHamil(partIndChild, partInd) + (1.0_dp * child(1))
        avHistHamil(partInd, partIndChild) = &
            avHistHamil(partInd, partIndChild) + (1.0_dp * child(1))

        ! Call the normal stats routine
        call new_child_stats_normal(iter_data, iLutI, nJ, iLutJ, ic, &
                                    walkExLevel, child, parent_flags, &
                                    part_type)

    end subroutine

    subroutine new_child_stats_normal(iter_data, iLutI, nJ, iLutJ, ic, &
                                      walkExLevel, child, parent_flags, &
                                      part_type)

        integer(kind=n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer, intent(in) :: ic, walkExLevel, parent_flags, nJ(nel)
        integer, intent(in) :: part_type
        real(dp), dimension(lenof_sign), intent(in) :: child
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer(n_int) :: iUnused
        integer :: run
        integer :: i

        unused_var(nJ); unused_var(walkExLevel)

        ! Write out some debugging information if asked
        IFDEBUG(FCIMCDebug, 3) then
        write(iout, "(A)", advance='no') "Creating "
        do i = 1, lenof_sign
            write(iout, "(f10.5)", advance='no') child(i)
        end do
        write(iout, "(A)", advance='no') " particles: "
        write(iout, "(A,2I4,A)", advance='no') &
            "Parent flag: ", parent_flags, part_type
        call writebitdet(iout, ilutJ, .true.)
        call neci_flush(iout)
        end if

        ! Count the number of children born
        ! in the real-time fciqmc, it is probably good to keep track of the
        ! stats of the 2 distinct RK loops .. use a global variable for the step
        ! RT_M_Merge: Merge conflict with master due to introduction of kmneci resolved.
        ! For rmneci, the __REALTIME case will have to be adapted to inum_runs>1

        ! rmneci_setup: Added multirun support for real-time case

        if (.not. t_real_time_fciqmc .or. runge_kutta_step == 2) then
#if defined( CMPLX_)
            do run = 1, inum_runs
                NoBorn(run) = NoBorn(run) + sum(abs(child(min_part_type(run):max_part_type(run))))
                if (ic == 1) SpawnFromSing(run) = SpawnFromSing(run) + sum(abs(child(min_part_type(run):max_part_type(run))))

                ! Count particle blooms, and their sources
                if (sum(abs(child(min_part_type(run):max_part_type(run)))) > InitiatorWalkNo) then
                    bloom_count(ic) = bloom_count(ic) + 1
                    bloom_sizes(ic) = max(real(sum(abs(child(min_part_type(run):max_part_type(run)))), dp), bloom_sizes(ic))
                end if
            end do
#else
            NoBorn = NoBorn + abs(child)
            if (ic == 1) SpawnFromSing = SpawnFromSing + abs(child)

            ! Count particle blooms, and their sources
            if (abs(child(part_type)) > InitiatorWalkNo) then
                bloom_count(ic) = bloom_count(ic) + 1
                bloom_sizes(ic) = max(real((abs(child(part_type))), dp), bloom_sizes(ic))
            end if
#endif
        end if
        if (.not. tPrecond) iter_data%nborn = iter_data%nborn + abs(child)

        ! Histogram the excitation levels as required
        if (tHistExcitToFrom) &
            call add_hist_excit_tofrom(ilutI, ilutJ, child)

        ! Avoid compiler warnings
        iUnused = iLutI(0); iUnused = iLutJ(0)

        end subroutine

        function attempt_die_normal(DetCurr, Kii, realwSign, WalkExcitLevel, DetPosition) result(ndie)

            ! Should we kill the particle at determinant DetCurr.
            ! The function allows multiple births (if +ve shift), or deaths from
            ! the same particle. The returned number is the number of deaths if
            ! positive, and the
            !
            ! In:  DetCurr - The determinant to consider
            !      Kii     - The diagonal matrix element of DetCurr (-Ecore)
            !      wSign   - The sign of the determinant being considered. If
            !                |wSign| > 1, attempt to die multiple particles at
            !                once (multiply probability of death by |wSign|)
            ! Ret: ndie    - The number of deaths (if +ve), or births (If -ve).

            integer, intent(in) :: DetCurr(nel)
            real(dp), dimension(lenof_sign), intent(in) :: RealwSign
            real(dp), intent(in) :: Kii
            real(dp), dimension(lenof_sign) :: ndie
            integer, intent(in) :: WalkExcitLevel
            integer, intent(in), optional :: DetPosition
            character(*), parameter :: t_r = 'attempt_die_normal'

            integer(kind=n_int) :: iLutnI(0:niftot)
            integer :: N_double

            real(dp) :: probsign, r
            real(dp), dimension(inum_runs) :: fac
            integer :: i, run, iUnused
#ifdef CMPLX_
            real(dp) :: rat(2)
#else
            real(dp) :: rat(1)
#endif
            real(dp) :: shift
            real(dp) :: relShiftOffset ! ShiftOffset relative to Hii

            do i = 1, inum_runs
                if (t_cepa_shift) then

                    fac(i) = tau * (Kii - (DiagSft(i) - cepa_shift(i, WalkExcitLevel)))
                    call log_death_magnitude(Kii - (DiagSft(i) - cepa_shift(i, WalkExcitLevel)))

                else if (tSkipRef(i) .and. all(DetCurr == projEdet(:, i))) then
                    !If we are fixing the population of reference det, skip death/birth
                    fac(i) = 0.0
                else
                    ! rescale the shift. It is better to wait till the shift starts varying
                    if (tAdaptiveShift .and. .not. tSinglePartPhase(i)) then
                        if (tAS_Offset) then
                            !Since DiagSft is relative to Hii, ShiftOffset should also be adjusted
                            relShiftOffset = ShiftOffset(i) - Hii
                        else
                            !Each replica's shift should be scaled using its own reference
                            relShiftOffset = proje_ref_energy_offsets(i) ! i.e. E(Ref(i)) - Hii
                        end if
                      shift = relShiftOffset + (DiagSft(i) - relShiftOffset) * shiftFactorFunction(DetPosition, i, mag_of_run(RealwSign, i))
                    else
                        shift = DiagSft(i)
                    end if

                    fac(i) = tau * (Kii - shift)

                    if (t_precond_hub .and. Kii > 10.0_dp) then
                        call log_death_magnitude(((Kii - shift) / (1.0_dp + Kii)))
                    else
                        call log_death_magnitude(Kii - shift)
                    end if
                end if

                if (t_precond_hub .and. Kii > 10.0_dp) then
                    fac(i) = fac(i) / (1.0_dp + Kii)
                end if
            end do

            if (any(fac > 1.0_dp)) then
                if (any(fac > 2.0_dp)) then
                    if (tSearchTauOption .or. t_hist_tau_search_option) then
                        ! If we are early in the calculation, and are using tau
                        ! searching, then this is not a big deal. Just let the
                        ! searching deal with it
                        write(iout, '("** WARNING ** Death probability > 2: Algorithm unstable.")')
                        write(iout, '("** WARNING ** Truncating spawn to ensure stability")')
                        do i = 1, inum_runs
                            fac(i) = min(2.0_dp, fac(i))
                        end do
                    else
                        print *, "Diag sft", DiagSft
                        print *, "Death probability", fac
                        call stop_all(t_r, "Death probability > 2: Algorithm unstable. Reduce timestep.")
                    end if
                else
                    write (iout, '("** WARNING ** Death probability > 1: Creating Antiparticles. "&
                        & //"Timestep errors possible: ")', advance='no')
                    do i = 1, inum_runs
                        write(iout, '(1X,f13.7)', advance='no') fac(i)
                    end do
                    write(iout, '()')
                end if
            end if

            if ((tRealCoeffByExcitLevel .and. (WalkExcitLevel <= RealCoeffExcitThresh)) &
                .or. tAllRealCoeff) then
                do run = 1, inum_runs
                    ndie(min_part_type(run)) = fac(run) * abs(realwSign(min_part_type(run)))
#ifdef CMPLX_
                    ndie(max_part_type(run)) = fac(run) * abs(realwSign(max_part_type(run)))
#endif
                end do
            else
                do run = 1, inum_runs

                    ! Subtract the current value of the shift, and multiply by tau.
                    ! If there are multiple particles, scale the probability.

                    rat(:) = fac(run) * abs(realwSign(min_part_type(run):max_part_type(run)))

                    ndie(min_part_type(run):max_part_type(run)) = real(int(rat), dp)
                    rat(:) = rat(:) - ndie(min_part_type(run):max_part_type(run))

                    ! Choose to die or not stochastically
                    r = genrand_real2_dSFMT()
                    if (abs(rat(1)) > r) ndie(min_part_type(run)) = &
                        ndie(min_part_type(run)) + real(nint(sign(1.0_dp, rat(1))), dp)
#ifdef CMPLX_
                    r = genrand_real2_dSFMT()
                    if (abs(rat(2)) > r) ndie(max_part_type(run)) = &
                        ndie(max_part_type(run)) + real(nint(sign(1.0_dp, rat(2))), dp)
#endif
                end do
            end if

            ! Avoid compiler warnings
            iUnused = DetCurr(1)

        end function attempt_die_normal

        function attempt_die_precond(DetCurr, Kii, realwSign, WalkExcitLevel, DetPos) result(ndie)
            ! Should we kill the particle at determinant DetCurr.
            ! The function allows multiple births (if +ve shift), or deaths from
            ! the same particle. The returned number is the number of deaths if
            ! positive, and the
            !
            ! In:  DetCurr - The determinant to consider
            !      Kii     - The diagonal matrix element of DetCurr (-Ecore)
            !      wSign   - The sign of the determinant being considered. If
            !                |wSign| > 1, attempt to die multiple particles at
            !                once (multiply probability of death by |wSign|)
            !      DetPos  - Position of the spawning determinant
            ! Ret: ndie    - The number of deaths (if +ve), or births (If -ve).

            integer, intent(in) :: DetCurr(nel)
            real(dp), dimension(lenof_sign), intent(in) :: RealwSign
            real(dp), intent(in) :: Kii
            real(dp), dimension(lenof_sign) :: ndie
            integer, intent(in) :: WalkExcitLevel

            integer, intent(in), optional :: DetPos
            character(*), parameter :: t_r = 'attempt_die_normal'

            real(dp) :: probsign, r
            real(dp), dimension(inum_runs) :: fac
            integer :: i, run, iUnused
#ifdef CMPLX_
            real(dp) :: rat(2)
#else
            real(dp) :: rat(1)
#endif
            unused_var(Kii); unused_var(DetPos)

            do i = 1, inum_runs
                fac(i) = tau

                ! And for tau searching purposes
                call log_death_magnitude(1.0_dp)
            end do

            if ((tRealCoeffByExcitLevel .and. (WalkExcitLevel <= RealCoeffExcitThresh)) &
                .or. tAllRealCoeff) then
                do run = 1, inum_runs
                    ndie(min_part_type(run)) = fac(run) * abs(realwSign(min_part_type(run)))
#ifdef CMPLX_
                    ndie(max_part_type(run)) = fac(run) * abs(realwSign(max_part_type(run)))
#endif
                end do
            else
                do run = 1, inum_runs

                    ! Subtract the current value of the shift, and multiply by tau.
                    ! If there are multiple particles, scale the probability.

                    rat(:) = fac(run) * abs(realwSign(min_part_type(run):max_part_type(run)))

                    ndie(min_part_type(run):max_part_type(run)) = real(int(rat), dp)
                    rat(:) = rat(:) - ndie(min_part_type(run):max_part_type(run))

                    ! Choose to die or not stochastically
                    r = genrand_real2_dSFMT()
                    if (abs(rat(1)) > r) ndie(min_part_type(run)) = &
                        ndie(min_part_type(run)) + real(nint(sign(1.0_dp, rat(1))), dp)
#ifdef CMPLX_
                    r = genrand_real2_dSFMT()
                    if (abs(rat(2)) > r) ndie(max_part_type(run)) = &
                        ndie(max_part_type(run)) + real(nint(sign(1.0_dp, rat(2))), dp)
#endif
                end do
            end if

            ! Avoid compiler warnings
            iUnused = DetCurr(1)
            iUnused = DetPos

        end function attempt_die_precond

!------------------------------------------------------------------------------------------!

        pure function powerScaleFunction(hdiag) result(Si)
            implicit none

            real(dp), intent(in) :: hdiag
            real(dp) :: Si

            Si = 1.0 / ((sFAlpha * (hdiag) + 1)**sFBeta)
        end function powerScaleFunction

!------------------------------------------------------------------------------------------!

        pure function expScaleFunction(hdiag) result(Si)
            implicit none

            real(dp), intent(in) :: hdiag
            real(dp) :: Si

            Si = 1.0 / (sfBeta * exp(sFAlpha * hdiag))
        end function expScaleFunction

!------------------------------------------------------------------------------------------!

        pure function expCOScaleFunction(hdiag) result(Si)
            implicit none

            real(dp), intent(in) :: hdiag
            real(dp) :: Si

            Si = (1 - sFbeta) / (exp(sFAlpha * hdiag)) + sFbeta
        end function expCOScaleFunction

!------------------------------------------------------------------------------------------!

        pure function negScaleFunction(hdiag) result(Si)
            implicit none

            real(dp), intent(in) :: hdiag
            real(dp) :: Si

            unused_var(hdiag)

            Si = -1

        end function negScaleFunction

!------------------------------------------------------------------------------------------!

        pure function expShiftFactorFunction(pos, run, pop) result(f)
            implicit none
            ! Exponential scale function for the shift
            ! Input: pos - position of given determinant in CurrentDets
            ! Input: run - run for which the factor is needed
            ! Input: pop - population of given determinant
            ! Output: f - scaling factor for the shift
            integer, intent(in) :: pos
            integer, intent(in) :: run
            real(dp), intent(in) :: pop
            real(dp) :: f
#ifdef DEBUG_
            ! Disable compiler warnings
            real(dp) :: dummy
            dummy = pos
            dummy = run
#endif

            f = 1.0 - exp(-pop / EAS_Scale)

        end function expShiftFactorFunction

!------------------------------------------------------------------------------------------!

        pure function constShiftFactorFunction(pos, run, pop) result(f)
            implicit none
            ! Dummy scale function for the shift: S' = S
            ! Input: pos - position of given determinant in CurrentDets
            ! Input: run - run for which the factor is needed
            ! Input: pop - population of given determinant
            ! Output: f - scaling factor for the shift
            integer, intent(in) :: pos
            integer, intent(in) :: run
            real(dp), intent(in) :: pop
            real(dp) :: f
#ifdef WARNING_WORKAROUND_
            ! Disable compiler warnings
            real(dp) :: dummy
            dummy = pos
            dummy = run
            dummy = pop
#endif
            f = 1.0
        end function constShiftFactorFunction

!------------------------------------------------------------------------------------------!

        pure function linearShiftFactorFunction(pos, run, pop) result(f)
            implicit none
            ! Piecewise-linear scale function for the shift
            ! Input: pos - position of given determinant in CurrentDets
            ! Input: run - run for which the factor is needed
            ! Input: pop - population of given determinant
            ! Output: f - scaling factor for the shift
            integer, intent(in) :: pos
            integer, intent(in) :: run
            real(dp), intent(in) :: pop
            real(dp) :: f, slope
#ifdef WARNING_WORKAROUND_
            ! Disable compiler warnings
            real(dp) :: dummy
            dummy = pos
            dummy = run
#endif

            if (pop > InitiatorWalkNo) then
                f = 1.0
            else if (pop < LAS_Sigma) then
                f = 0.0
            else
                if (InitiatorWalkNo.isclose.LAS_Sigma) then
                    !In this case the slope is ill-defined.
                    !Since initiators are strictly large than InitiatorWalkNo, set shift to zero
                    f = 0.0
                else
                    !Apply linear modifcation that equals F1 at Sigma and F2 at InitatorWalkNo
                    slope = (LAS_F2 - LAS_F1) / (InitiatorWalkNo - LAS_Sigma)
                    f = (LAS_F1 + (pop - LAS_Sigma) * slope)
                end if
            end if
        end function linearShiftFactorFunction
!------------------------------------------------------------------------------------------!

        pure function autoShiftFactorFunction(pos, run, pop) result(f)
            implicit none
            ! Scale function for the shift based on the ratio of reject spawns
            ! Input: pos - position of given determinant in CurrentDets
            ! Input: run - run for which the factor is needed
            ! Input: pop - population of given determinant
            ! Output: f - scaling factor for the shift
            integer, intent(in) :: pos
            integer, intent(in) :: run
            real(dp), intent(in) :: pop
            real(dp) :: f, tot, acc, tmp

            unused_var(pop)

            tot = get_tot_spawns(pos, run)
            acc = get_acc_spawns(pos, run)

            if (test_flag(CurrentDets(:, pos), get_initiator_flag_by_run(run))) then
                tmp = 1.0
            else if (tot > AAS_Thresh) then
                tmp = acc / tot
            else
                tmp = 0.0
            end if
            !The factor is actually never zero, because at least the parent is occupied
            !As a heuristic, we use the connectivity of HF
            if (tmp < AAS_Cut) then
                tmp = AAS_Cut
            end if
            tmp = (tmp + AAS_Const) / (1 + AAS_Const)
            f = tmp**AAS_Expo

        end function autoShiftFactorFunction

    end module
