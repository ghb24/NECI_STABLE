#include "macros.h"
module fcimc_output

    use SystemData, only: nel, tHPHF, tFixLz, tMolpro, tMolproMimic, MolproID, &
                          tGen_4ind_weighted, tGen_4ind_2, tGUGA, tGen_sym_guga_mol, &
                          tgen_guga_crude, t_new_real_space_hubbard, t_no_ref_shift

    use LoggingData, only: tLogComplexPops, tMCOutput, tCalcInstantS2, &
                           tCalcInstantS2Init, instant_s2_multiplier_init, &
                           instant_s2_multiplier, tPrintFCIMCPsi, &
                           iWriteHistEvery, tDiagAllSpaceEver, OffDiagMax, &
                           OffDiagBinRange, tCalcVariationalEnergy, &
                           iHighPopWrite, tLogEXLEVELStats, StepsPrint, &
                           maxInitExLvlWrite, AllInitsPerExLvl, t_force_replica_output

    use hist_data, only: Histogram, AllHistogram, InstHist, AllInstHist, &
                         BeforeNormHist, iNoBins, BinRange, HistogramEnergy, &
                         AllHistogramEnergy

    use CalcData, only: tTruncInitiator, tTrialWavefunction, tReadPops, &
                        DiagSft, tSpatialOnlyHash, tOrthogonaliseReplicas, &
                        StepsSft, tPrintReplicaOverlaps, tStartTrialLater, &
                        tEN2, &
                        tSemiStochastic, t_truncate_spawns, tLogGreensfunction

    use DetBitOps, only: FindBitExcitLevel, count_open_orbs, EncodeBitDet, &
                         TestClosedShellDet

    use IntegralsData, only: frozen_orb_list, frozen_orb_reverse_map, &
                             nel_pre_freezing

    use DetCalcData, only: det, fcidets, ReIndex, NDet, NRow, HAMIL, LAB

    use bit_rep_data, only: test_flag, extract_sign

    use bit_reps, only: decode_bit_det, get_initiator_flag

    use semi_stoch_procs, only: global_most_populated_states, GLOBAL_RUN, core_space_weight

    use bit_rep_data, only: niftot, nifd, flag_initiator

    use hist, only: calc_s_squared_star, calc_s_squared

    use fcimc_helper, only: LanczosFindGroundE

    use hphf_integrals, only: hphf_diag_helement
    use Determinants, only: get_helement, writeDetBit
    use DeterminantData, only: write_det
    use adi_data, only: AllCoherentDoubles, AllIncoherentDets, nRefs, &
         ilutRefAdi, tAdiActive, nConnection, AllConnection

    use rdm_data, only: en_pert_main

    use Parallel_neci

    use FciMCData

    use constants

    use sort_mod

    use util_mod

    use real_time_data, only:  gf_count, &
                              normsize, snapShotOrbs, &
                              current_overlap, t_real_time_fciqmc, elapsedRealTime, &
                              elapsedImagTime, overlap_real, overlap_imag, dyn_norm_psi,&
                              dyn_norm_red, real_time_info, allPopSnapshot, numSnapshotOrbs

    use tc_three_body_data, only: tLMatCalc, lMatCalcStatsIters, &
                                  lMatCalcHit,   lMatCalcTot,    lMatCalcHUsed,   lMatCalcHSize
    use fortran_strings, only: str

    use guga_matrixElements, only: calcDiagMatEleGUGA_nI

    use matmul_mod, only: my_hpsi

    implicit none

contains

    SUBROUTINE WriteFciMCStatsHeader()
        integer :: j, k, run, offset
        integer(int64) :: i
        character(256) label
        character(32) tchar_r, tchar_i, tchar_j, tchar_k
        character(17) trunc_caption
        character(38) validExCaption

        call getProjEOffset()

        IF(iProcIndex == root) THEN
!Print out initial starting configurations
            write(stdout,*) ""
            IF(tTruncInitiator) THEN
               write(initiatorstats_unit,"(A2,A17,16A23)", advance = 'no') &
                    "# ","1.Step","2.TotWalk","3.Annihil","4.Died", &
                    & "5.Born","6.TotUniqDets",&
                    &               "7.InitDets","8.NonInitDets","9.InitWalks","10.NonInitWalks","11.AbortedWalks", &
                    "12. Removed Dets",  "13. Initiator Proj.E", "14. CoreNonInits"
               offset = 14
               if(tTrialWavefunction .or. tStartTrialLater) then
                  write(initiatorstats_unit,"(A)", advance = 'no') &
                  "15. TrialNumerators (inits)   16. TrialDenom (inits)"
                  offset = 16
               end if
                do k = 1, maxInitExLvlWrite
                   write(tchar_k,*) k+offset
                   write(tchar_r,*) k
                   tchar_r = trim(adjustl(tchar_k))//'. Inits on ex. lvl '//trim(adjustl(tchar_r))
                   write(initiatorstats_unit,'(1x,a)', advance = 'no') &
                        trim(adjustl(tchar_r))
                end do
                write(initiatorstats_unit,'()', advance = 'yes')
            end if
            IF(tLogComplexPops) THEN
                write(complexstats_unit,"(A)") '#   1.Step  2.Shift     3.RealShift     4.ImShift   5.TotParts      " &
                & //"6.RealTotParts      7.ImTotParts'
            end if
            if (tLogEXLEVELStats) then
                write(EXLEVELStats_unit, '(a)', advance='no') '# 1.Step'
                k = 1
                do run = 1, inum_runs
                    tchar_r = ''
                    if (inum_runs>1) then
                        write(tchar_r,*)run
                        tchar_r = '(run='//trim(adjustl(tchar_r))//')'
                    end if
                    do i = 0, 2
                        write(tchar_i,*)i
                        do j = 0, NEl
                            k = k + 1
                            write(tchar_j,*)j
                            write(tchar_k,*)k
                            write(EXLEVELStats_unit, '(1x,a)', &
                                   advance='no') trim(adjustl(tchar_k))// &
                                   &'.W'//trim(adjustl(tchar_j))//'^'// &
                                   trim(adjustl(tchar_i))// &
                                   trim(adjustl(tchar_r))
                        end do ! j
                    end do ! i
                end do ! run
                write(EXLEVELStats_unit, '()', advance='yes')
            end if ! tLogEXLEVELStats

#ifdef CMPLX_
            if(tMCOutput) then
                write(stdout, '(a)') "       Step     Shift      WalkerCng(Re)  &
                       &WalkerCng(Im)    TotWalkers(Re)   TotWalkers(Im)    &
                       &Proj.E(Re)   ProjE(Im)     Proj.E.ThisCyc(Re)  &
                       &Proj.E.ThisCyc(Im)   NoatHF(Re)   NoatHF(Im)   &
                       &NoatDoubs      AccRat     UniqueDets   NumDetsSpawned   IterTime"
            end if
            write(fcimcstats_unit, "(a,i4,a,l1,a,l1,a,l1)") &
                   "# FCIMCStats VERSION 2 - COMPLEX : NEl=", nel, &
                   " HPHF=", tHPHF, ' Lz=', tFixLz, &
                   ' Initiator=', tTruncInitiator
            write(fcimcstats_unit, "(a)", advance = 'no') &
                   "#     1.Step   2.Shift    3.WalkerCng(Re)  &
                   &4.WalkerCng(Im)   5.TotWalkers(Re)  6.TotWalkers(Im)  &
                   &7.Proj.E(Re)   8.Proj.E(Im)   9.Proj.E.ThisCyc(Re)  &
                   &10.Proj.E.ThisCyc(Im)  11.Tot-Proj.E.ThisCyc(Re)  12.NoatHF(Re)  &
                   &13.NoatHF(Im)  14.NoatDoubs  15.AccRat  16.UniqueDets  17.IterTime &
                   &18.FracSpawnFromSing  19.WalkersDiffProc  20.TotImagTime &
                   &  21.HFInstShift  22.TotInstShift  &
                   &23.HFContribtoE(Both)  &
                   &24.NumContribtoE(Re)  &
                   &25.NumContribtoE(Im)  26.HF weight   27.|Psi|    &
                   &28.Inst S^2  29.PartsDiffProc   30.MaxCycSpawn"
! Dongxia comment 28-30 off because they are not printed out
! 28.SpawnedParts  29.MergedParts  30.Zero elems   31.PartsDiffProc   32.MaxCycSpawn"
            if (tTrialWavefunction .or. tStartTrialLater) then
                   write(fcimcstats_unit, "(A)", advance = 'no') &
                   "  31.TrialNumerator(Re)  32.TrialNumerator(Im)  33.TrialDenom(Re)  &
                   &  34.TrialDenom(Im)  35.TrialOverlap  36.TrialProjE(Re)  37.TrialProjE(Im)"
            end if

            write(fcimcstats_unit, "()", advance = 'yes')

#elif defined(DOUBLERUN_)
            write(fcimcstats_unit2, "(a,i4,a,l1,a,l1,a,l1)") &
                  "# FCIMCStats VERSION 2 - REAL : NEl=", nel, &
                  " HPHF=", tHPHF, ' Lz=', tFixLz, &
                  ' Initiator=', tTruncInitiator
            write(fcimcstats_unit2, "(A)", advance = 'no') &
                  "#     1.Step   2.Shift    3.WalkerCng  4.GrowRate     &
                  &5.TotWalkers  6.Annihil  7.NoDied  8.NoBorn  &
                  &9.Proj.E       10.Av.Shift 11.Proj.E.ThisCyc  12.NoatHF &
                  &13.NoatDoubs  14.AccRat  15.UniqueDets  16.IterTime &
                  &17.FracSpawnFromSing  18.WalkersDiffProc  19.TotImagTime  &
                  &20.ProjE.ThisIter  21.HFInstShift  22.TotInstShift  &
                  &23.Tot-Proj.E.ThisCyc   24.HFContribtoE  25.NumContribtoE &
                  &26.HF weight    27.|Psi|     28.Inst S^2 29.Inst S^2 30.AbsProjE &
                  &31.PartsDiffProc    32.|Semistoch|/|Psi|     33.MaxCycSpawn"
           if (tTrialWavefunction .or. tStartTrialLater) then
                  write(fcimcstats_unit2, "(A)", advance = 'no') &
                  "  34.TrialNumerator  35.TrialDenom  36.TrialOverlap"
              trunc_caption = "  37. TruncWeight"
           else
              trunc_caption = "  34. TruncWeight"
           end if
           if(t_truncate_spawns) write(fcimcstats_unit2, "(A)", advance = 'no') &
                trunc_caption

           write(fcimcstats_unit2, "()", advance = 'yes')
#endif
#ifndef CMPLX_
            if(tMCOutput) then
                write(stdout, "(A)", advance = 'no') "        Step    Shift           &
                      &WalkerCng       GrowRate        TotWalkers      Annihil         &
                      &NoDied          NoBorn          Proj.E          Av.Shift        &
                      &Proj.E.Cyc"
                if (tTrialWavefunction .or. tStartTrialLater) write(stdout, "(A)", advance = 'no') &
                      "    Trial.E.Cyc "
                write(stdout, "(A)", advance = 'yes') "      NoatHF          NoatDoubs       &
                &AccRat        UniqueDets    NumDetsSpawned   IterTime"
            end if
            write(fcimcstats_unit, "(a,i4,a,l1,a,l1,a,l1)") &
                  "# FCIMCStats VERSION 2 - REAL : NEl=", nel, &
                  " HPHF=", tHPHF, ' Lz=', tFixLz, &
                  ' Initiator=', tTruncInitiator
            write(fcimcstats_unit, "(A)", advance = 'no') &
                  "#     1.Step   2.Shift    3.WalkerCng  4.GrowRate     &
                  &5.TotWalkers  6.Annihil  7.NoDied  8.NoBorn  &
                  &9.Proj.E       10.Av.Shift 11.Proj.E.ThisCyc  12.NoatHF &
                  &13.NoatDoubs  14.AccRat  15.UniqueDets  16.IterTime &
                  &17.FracSpawnFromSing  18.WalkersDiffProc  19.TotImagTime  &
                  &20.ProjE.ThisIter  21.HFInstShift  22.TotInstShift  &
                  &23.Tot-Proj.E.ThisCyc   24.HFContribtoE  25.NumContribtoE &
                  &26.HF weight    27.|Psi|     28.Inst S^2 &
                  &29.Inst S^2   30.AbsProjE   31.PartsDiffProc &
                  &32.|Semistoch|/|Psi|  33.MaxCycSpawn "
           if (tTrialWavefunction .or. tStartTrialLater) then
              write(fcimcstats_unit, "(A)", advance = 'no') &
                   "  34.TrialNumerator  35.TrialDenom  36.TrialOverlap"
              validExCaption = "  37.InvalidExcits  38. ValidExcits  "
              trunc_caption = "  39. TruncWeight"
           else
              trunc_caption = "  36. TruncWeight"
              validExCaption = "  34.InvalidExcits  35. ValidExcits  "
           end if
           write(fcimcstats_unit, "(A)", advance = 'no') validExCaption
           if(t_truncate_spawns) write(fcimcstats_unit, "(A)", advance = 'no') &
                trunc_caption

           write(fcimcstats_unit, "()", advance = 'yes')

#endif

        end if

    END SUBROUTINE WriteFciMCStatsHeader

    subroutine WriteFCIMCStats()

        INTEGER :: i, j, run
        real(dp),dimension(inum_runs) :: FracFromSing
        real(dp) :: E_ref_tmp(inum_runs)

        ! get the offset for the projected energy (i.e. reference energy)
        call getProjEOffset()

        ! What is the current value of S2
        if (tCalcInstantS2) then
            if (mod(iter / StepsPrint, instant_s2_multiplier) == 0) then
                if (tSpatialOnlyhash) then
                    curr_S2 = calc_s_squared (.false.)
                else
                    curr_S2 = calc_s_squared_star (.false.)
                end if
            end if
        else
            curr_S2 = -1
        end if

        ! What is the current value of S2 considering only initiators
        if (tCalcInstantS2Init) then
            if (mod(iter / StepsPrint, instant_s2_multiplier_init) == 0) then
                if (tSpatialOnlyhash) then
                    curr_S2_init = calc_s_squared (.true.)
                else
                    curr_S2_init = calc_s_squared_star (.true.)
                end if
            end if
        else
            curr_S2_init = -1
        end if

        !To prevent /0 problems
        do run=1,inum_runs
            if(.not. near_zero(AllNoBorn(run))) then
                FracFromSing(run)=real(AllSpawnFromSing(run),dp) / real(AllNoBorn(run),dp)
            else
                FracFromSing(run)=0.0_dp
            end if

            if(t_no_ref_shift)then
                if (tHPHF) then
                    E_ref_tmp(run) = hphf_diag_helement (ProjEDet(:,run), iLutRef(:,run))
                else if(tguga)then
                    E_ref_tmp(run) = calcDiagMatEleGUGA_nI(ProjEDet(:,run))
                else
                    E_ref_tmp(run) = get_helement (ProjEDet(:,run), ProjEDet(:,run), 0)
                end if
            else
                E_ref_tmp(run) = 0.0_dp
            end if


        end do

        if (iProcIndex == root) then

#ifdef CMPLX_
            write(fcimcstats_unit,"(I12,5G16.7,8G18.9e3,&
                                  &G13.5,I12,G13.5,G17.5,I13,G13.5,8G18.9e3,I13,&

                                  &g16.7)",advance='no') &
                Iter + PreviousCycles, &                !1.
                DiagSft + E_ref_tmp, &                  !2.
                AllTotParts(1) - AllTotPartsLastOutput(1), &   !3.
                AllTotParts(2) - AllTotPartsLastOutput(2), &   !4.
                AllTotParts(1), &                       !5.
                AllTotParts(2), &                       !6.
                real(ProjectionE, dp), &                !7.     real \sum[ nj H0j / n0 ]
                aimag(projectionE), &                   !8.     Im   \sum[ nj H0j / n0 ]
                real(proje_iter, dp), &                 !9.
                aimag(proje_iter), &                    !10.
                real(proje_iter,dp) + OutputHii, &            !11.
                AllNoatHF(1), &                         !12.
                AllNoatHF(2), &                         !13.
                AllNoatDoubs, &                         !14.
                AccRat, &                               !15.
                AllTotWalkers, &                        !16.
                IterTime, &                             !17.
                FracFromSing(1), &                      !18.
                WalkersDiffProc, &                           !19.
                TotImagTime, &                               !20.
                HFShift, &                                   !21.
                InstShift, &                                 !22.
                real((AllHFOut*conjg(AllHFOut)),dp), &     !23 |n0|^2  denominator for both calcs
                real((AllENumOut*conjg(AllHFOut)),dp), &   !24. Re[\sum njH0j]xRe[n0]+Im[\sum njH0j]xIm[n0]   No div by StepsPrint
                aimag(AllENumOut*conjg(AllHFOut)), &       !25.Im[\sum njH0j]xRe[n0]-Re[\sum njH0j]xIm[n0]   since no physicality
                sqrt(sum(AllNoatHF**2)) / norm_psi, & !26
                norm_psi, &                           !27
                curr_S2, &                            !28
                PartsDiffProc, &                      !29
                all_max_cyc_spawn                     !30
                if (tTrialWavefunction .or. tStartTrialLater) then
                    write(fcimcstats_unit, "(7(1X,es18.11))", advance = 'no') &
                    (tot_trial_numerator(1) / StepsPrint), &              ! 31. 32
                    (tot_trial_denom(1) / StepsPrint), &                  ! 33. 34
                    abs((tot_trial_denom(1) / (norm_psi(1)*StepsPrint))), &  ! 35.
                    tot_trial_numerator(1)/tot_trial_denom(1)           ! 36. 37.
                end if
                write(fcimcstats_unit, "()", advance = 'yes')

            if(tMCOutput) then
                write(stdout, "(I12,13G16.7,2I12,G13.5)") &
                    Iter + PreviousCycles, &
                    DiagSft + E_ref_tmp, &
                    AllTotParts(1) - AllTotPartsLastOutput(1), &
                    AllTotParts(2) - AllTotPartsLastOutput(2), &
                    AllTotParts(1), &
                    AllTotParts(2), &
                    real(ProjectionE, dp), &
                    aimag(ProjectionE), &
                    real(proje_iter, dp), &
                    aimag(proje_iter), &
                    AllNoatHF(1), &
                    AllNoatHF(2), &
                    AllNoatDoubs, &
                    AccRat, &
                    AllTotWalkers, &
                    nspawned_tot, &
                    IterTime
            end if
            if (tTruncInitiator) then
               write(initiatorstats_unit,"(I12,4G16.7,3I20,4G16.7,F16.9,2G16.7,2I20,1G16.7,1I20)", &
                    advance = 'no')&
                   Iter + PreviousCycles, sum(AllTotParts), &
                   AllAnnihilated(1), AllNoDied(1), AllNoBorn(1), AllTotWalkers,&
                   AllNoInitDets(1), AllNoNonInitDets(1), AllNoInitWalk(1), &
                   AllNoNonInitWalk(1),AllNoAborted(1), AllNoRemoved(1), &
                   inits_proje_iter(1) + Hii, all_n_core_non_init
               if(tTrialWavefunction .or. tStartTrialLater) &
                    write(initiatorstats_unit, "(2G16.7)", advance = 'no') &
                    tot_init_trial_numerator(1)/StepsPrint, tot_init_trial_denom(1)/StepsPrint
               do j = 1, maxInitExLvlWrite
                  write(initiatorstats_unit,'(1I20)', advance ='no') AllInitsPerExLvl(j)/StepsPrint
               end do
               write(initiatorstats_unit,'()', advance = 'yes')
            end if
            if (tLogComplexPops) then
                write(complexstats_unit,"(I12,6G16.7)") &
                    Iter + PreviousCycles, DiagSft, DiagSftRe, DiagSftIm, &
                    sum(AllTotParts), AllTotParts(1), AllTotParts(lenof_sign)
            end if
#elif defined(DOUBLERUN_)
            write(fcimcstats_unit2,"(i12,7g16.7,5g18.9e3,g13.5,i12,g13.5,g17.5,&
                                   &i13,g13.5,4g18.9e3,1X,2(es18.11,1X),5g18.9e3,&
                                   &i13,2g16.7)",advance = 'no') &
                Iter + PreviousCycles, &                   ! 1.
                DiagSft(2) + E_ref_tmp(2), &                              ! 2.
                AllTotParts(2) - AllTotPartsLastOutput(2), &      ! 3.
                AllGrowRate(2), &                          ! 4.
                AllTotParts(2), &                          ! 5.
                AllAnnihilated(2), &                       ! 6.
                AllNoDied(2), &                            ! 7.
                AllNoBorn(2), &                            ! 8.
                ProjectionE(2), &                          ! 9.
                AvDiagSft(2), &                            ! 10.
                proje_iter(2), &                           ! 11.
                AllNoatHF(2), &                            ! 12.
                AllNoatDoubs(2), &                         ! 13.
                AccRat(2), &                               ! 14.
                AllTotWalkers, &                           ! 15.
                IterTime, &                                ! 16.
                FracFromSing(2), &                         ! 17.
                WalkersDiffProc, &                         ! 18.
                TotImagTime, &                             ! 19.
                0.0_dp, &                                  ! 20.
                HFShift(2), &                              ! 21.
                InstShift(2), &                            ! 22.
                proje_iter(2) + OutputHii, &                     ! 23.
                (AllHFOut(2) / StepsPrint), &                ! 24.
                (AllENumOut(2) / StepsPrint), &              ! 25.
                AllNoatHF(2) / norm_psi(2), &              ! 26.
                norm_psi(2), &                             ! 27.
                curr_S2(2), curr_S2_init(2), &             ! 28, 29.
                AbsProjE(2), &                             ! 30.
                PartsDiffProc, &                           ! 31.
                norm_semistoch(2)/norm_psi(2), &           ! 32.
                all_max_cyc_spawn                          ! 33.
                if (tTrialWavefunction .or. tStartTrialLater) then
                    write(fcimcstats_unit2, "(3(1X,es17.10))", advance = 'no') &
                    (tot_trial_numerator(2) / StepsPrint), &
                    (tot_trial_denom(2) / StepsPrint), &
                    abs(tot_trial_denom(2) / (norm_psi(2)*StepsPrint))
                end if
                if(t_truncate_spawns) then
                   write(fcimcstats_unit2, "(1X,es18.11)", advance = 'no') AllTruncatedWeight
                end if

                write(fcimcstats_unit2, "()", advance = 'yes')
#endif
#ifndef CMPLX_

            write(fcimcstats_unit,"(i12,7g16.7,5g18.9e3,g13.5,i12,g13.5,g17.5,&
                                   &i13,g13.5,4g18.9e3,1X,2(es18.11,1X),5g18.9e3,&
                                   &i13,4g16.7)",advance = 'no') &
                Iter + PreviousCycles, &                   ! 1.
                DiagSft(1) + E_ref_tmp(1), &                              ! 2.
                AllTotParts(1) - AllTotPartsLastOutput(1), &      ! 3.
                AllGrowRate(1), &                          ! 4.
                AllTotParts(1), &                          ! 5.
                AllAnnihilated(1)/StepsPrint, &            ! 6.
                AllNoDied(1)/StepsPrint, &                 ! 7.
                AllNoBorn(1)/StepsPrint, &                 ! 8.
                ProjectionE(1), &                          ! 9.
                AvDiagSft(1), &                            ! 10.
                proje_iter(1), &                           ! 11.
                AllNoatHF(1), &                            ! 12.
                AllNoatDoubs(1), &                         ! 13.
                AccRat(1), &                               ! 14.
                AllTotWalkers, &                           ! 15.
                IterTime, &                                ! 16.
                FracFromSing(1), &                         ! 17.
                WalkersDiffProc, &                         ! 18.
                TotImagTime, &                             ! 19.
                0.0_dp, &                                  ! 20.
                HFShift(1), &                              ! 21.
                InstShift(1), &                            ! 22.
                proje_iter(1) + OutputHii, &                     ! 23.
                (AllHFOut(1) / StepsPrint), &                ! 24.
                (AllENumOut(1) / StepsPrint), &              ! 25.
                AllNoatHF(1) / norm_psi(1), &              ! 26.
                norm_psi(1), &                             ! 27.
                curr_S2(1), curr_S2_init(1), &             ! 28, 29.
                AbsProjE(1), &                             ! 30.
                PartsDiffProc, &                           ! 31.
                norm_semistoch(1)/norm_psi(1), &           ! 32.
                all_max_cyc_spawn                          ! 33
                if (tTrialWavefunction .or. tStartTrialLater) then
                    write(fcimcstats_unit, "(3(1X,es18.11))", advance = 'no') &
                    (tot_trial_numerator(1) / StepsPrint), &              ! 34.
                    (tot_trial_denom(1) / StepsPrint), &                  ! 35.
                    abs((tot_trial_denom(1) / (norm_psi(1)*StepsPrint)))  ! 36.
                 end if
                 write(fcimcstats_unit, "(2g16.7)", advance = 'no') &
                      allNInvalidExcits, & ! 34/37.
                      allNValidExcits      ! 35/38.
                if(t_truncate_spawns) then
                   write(fcimcstats_unit, "(1X,es18.11)", advance = 'no') AllTruncatedWeight
                end if
                write(fcimcstats_unit, "()", advance = 'yes')

            if(tMCOutput) then
                write(stdout, "(I12,10G16.7)", advance = 'no') &
                    Iter + PreviousCycles, &
                    DiagSft(1)+E_ref_tmp(1), &
                    AllTotParts(1) - AllTotPartsLastOutput(1), &
                    AllGrowRate(1), &
                    AllTotParts(1), &
                    AllAnnihilated(1), &
                    AllNoDied(1), &
                    AllNoBorn(1), &
                    ProjectionE(1), &
                    AvDiagSft(1), &
                    proje_iter(1)
                if (tTrialWavefunction) then
                     write(stdout, "(G20.11)", advance = 'no') &
                         (tot_trial_numerator(1)/tot_trial_denom(1))
                else if (tStartTrialLater) then
                     write(stdout, "(G20.11)", advance = 'no') 0.0_dp
                end if
                write(stdout, "(3G16.7,2I12,G13.5)", advance = 'yes') &
                    AllNoatHF(1), &
                    AllNoatDoubs(1), &
                    AccRat(1), &
                    AllTotWalkers, &
                    nspawned_tot, &
                    IterTime
            end if
            if (tTruncInitiator) then
               write(initiatorstats_unit,"(I12,4G16.7,3I20,4G16.7,F16.9,2G16.7,2I20,1G16.7,1I20)", &
                    advance = 'no')&
                   Iter + PreviousCycles, AllTotParts(1), &
                   AllAnnihilated(1), AllNoDied(1), AllNoBorn(1), AllTotWalkers,&
                   AllNoInitDets(1), AllNoNonInitDets(1), AllNoInitWalk(1), &
                   AllNoNonInitWalk(1),AllNoAborted(1), AllNoRemoved(1), &
                   inits_proje_iter(1) + Hii, all_n_core_non_init
               if(tTrialWavefunction .or. tStartTrialLater) &
                    write(initiatorstats_unit, "(2G16.7)", advance = 'no') &
                    tot_init_trial_numerator(1)/StepsPrint, tot_init_trial_denom(1)/StepsPrint
               do j = 1, maxInitExLvlWrite
                  write(initiatorstats_unit,'(1I20)', advance ='no') AllInitsPerExLvl(j)
               end do
               write(initiatorstats_unit,'()', advance = 'yes')
            end if
#endif
            if (tLogEXLEVELStats) then
                write(EXLEVELStats_unit, '(i12)', advance='no') &
                      &Iter + PreviousCycles
                do run = 1, inum_runs
                    do i = 0, 2
                        do j = 0, NEl
                            write(EXLEVELStats_unit, '(1x,G18.9e3)', &
                                   advance='no') AllEXLEVEL_WNorm(i,j,run)
                        end do ! j
                    end do ! i
                end do ! run
                write(EXLEVELStats_unit, '()', advance='yes')
            end if ! tLogEXLEVELStats

            if (tMCOutput .and. tLMatCalc .and. mod(Iter, lMatCalcStatsIters) == 0) then
                write(stdout, *) "============ LMatCalc Caching Stats ==============="
                write(stdout, *) "LMatCalc Cache Fill Ratio: ", &
                    real(lMatCalcHUsed,dp)/real(lMatCalcHSize,dp)
                write(stdout, *) "LMatCalc Cache Hit Rate  : ", lMatCalcHit/real(lMatCalcTot)
                lMatCalcHit = 0
                lMatCalcTot = 0
                write(stdout, *) "==================================================="
            end if

            if(tMCOutput) then
                call neci_flush(stdout)
            end if
            call neci_flush(fcimcstats_unit)
            if (inum_runs.eq.2) call neci_flush(fcimcstats_unit2)
            if (tLogEXLEVELStats) call neci_flush(EXLEVELStats_unit)

        end if

    end subroutine WriteFCIMCStats


    subroutine open_create_stats(stem,funit)

        integer, intent(in) :: funit
        character(*), intent(in) :: stem
        character(*), parameter :: t_r = 'open_create_fciqmc_stats'

        character(30) :: filename
        character(43) :: filename2
        character(12) :: num

        ! If we are using Molpro, then append the molpro ID to uniquely
        ! identify the output
        if (tMolpro .and. .not. tMolproMimic) then
            filename = stem // adjustl(MolproID)
        else
            filename = stem
        end if


        if (tReadPops) then

            ! If we are reading from a POPSFILE, then we want to continue an
            ! existing fciqmc_stats file if it exists.
            open(funit, file=filename, status='unknown', position='append')
        else

           call open_new_file(funit, filename)

        end if

    end subroutine open_create_stats

    subroutine write_fcimcstats2(iter_data, initial)

        ! Write output to our FCIMCStats file.
        ! This should be done in a nice general way, so that we make merges
        ! somewhat easier.
        !
        ! --> The column numbers are no longer as well fixed (oh well...)

        type(fcimc_iter_data), intent(in) :: iter_data
        logical, intent(in), optional :: initial

        ! Use a state type to keep things compact and tidy below.
        type(write_state_t), save :: state
        type(write_state_t), save :: state_i
        logical, save :: inited = .false.
        character(5) :: tmpc, tmpc2, tmgf
        integer :: p, q, iGf, run
        logical :: init
        real(dp) :: l1_norm

! Is in the interface to refactor the procedure lateron.
        unused_var(iter_data)

        call getProjEOffset()
        ! Provide default 'initial' option
        if (present(initial)) then
            state%init = initial
            if (tTruncInitiator) state_i%init = initial
        else
            state%init = .false.
            if (tTruncInitiator) state_i%init = .false.
        end if

        ! If the output file hasn't been opened yet, then create it.
        if (iProcIndex == Root .and. .not. inited) then
           call open_state_file('fciqmc_stats',state)
           ! For the initiator stats file here:
           if (tTruncInitiator) call open_state_file('initiator_stats',state_i)

           inited = .true.
        end if

        ! What is the current value of S2
        if (tCalcInstantS2) then
            if (mod(iter / StepsSft, instant_s2_multiplier) == 0) then
                if (tSpatialOnlyhash) then
                    curr_S2 = calc_s_squared (.false.)
                else
                    curr_S2 = calc_s_squared_star (.false.)
                end if
            end if
        else
            curr_S2 = -1
        end if

        ! What is the current value of S2 considering only initiators
        if (tCalcInstantS2Init) then
            if (mod(iter / StepsSft, instant_s2_multiplier_init) == 0) then
                if (tSpatialOnlyhash) then
                    curr_S2_init = calc_s_squared (.true.)
                else
                    curr_S2_init = calc_s_squared_star (.true.)
                end if
            end if
        else
            curr_S2_init = -1
        end if

        ! ------------------------------------------------
        ! This is where any calculation that needs multiple nodes should go
        ! ------------------------------------------------
        ! ------------------------------------------------

        if (iProcIndex == root) then

            ! Only do the actual outputting on the head node.
            call write_padding_init(state)
            call write_padding_init(state_i)

            ! And output the actual data!
            state%cols = 0
            state%cols_mc = 0
            state%mc_out = tMCOutput

            if(t_real_time_fciqmc) then
                l1_norm = 0.0
                do run = 1, inum_runs
                    l1_norm = l1_norm + mag_of_run(AllTotParts, run)
                end do
            end if
            call stats_out(state,.true., iter + PreviousCycles, 'Iter.')
            if (.not. tOrthogonaliseReplicas) then
               ! note that due to the averaging, the printed value is not necessarily
               ! an integer
                call stats_out(state,.true., sum(abs(AllTotParts))/inum_runs, &
                     'Tot. parts real')
                if(t_real_time_fciqmc) then
                    call stats_out(state,.true., real_time_info%time_angle,'Time rot. angle')
                    call stats_out(state,.false., l1_norm/inum_runs ,'L1 Norm')
                else
                    call stats_out(state,.true., sum(abs(AllNoatHF))/inum_runs, 'Tot. ref')
                end if
            end if

            if(.not. t_real_time_fciqmc) then
#ifdef CMPLX_
                call stats_out(state,.true., real(proje_iter_tot), 'Re Proj. E')
                call stats_out(state,.true., aimag(proje_iter_tot), 'Im Proj. E')
#else
                call stats_out(state,.true., proje_iter_tot, 'Proj. E (cyc)')
#endif
            end if
            call stats_out(state,.true., sum(DiagSft)/inum_runs, 'Shift. (cyc)')

            if(t_real_time_fciqmc) &
                call stats_out(state, .true., real(sum(dyn_norm_psi))/normsize, '|psi|^2')

            call stats_out(state,.false., sum(AllNoBorn), 'No. born')
            call stats_out(state,.false., sum(AllNoInitDets), 'No. Inits')
            if(t_real_time_fciqmc) then
                call stats_out(state,.false., TotImagTime, 'Elapsed complex time')
                call stats_out(state,.false., real_time_info%damping, 'eta')
                call stats_out(state,.false., IterTime, 'Iter. time')
            else
                call stats_out(state,.false., sum(AllAnnihilated), 'No. annihil')
            end if

            call stats_out(state,.false., sum(AllSumWalkersCyc), 'SumWalkersCyc')
            call stats_out(state,.false., sum(AllNoAborted), 'No aborted')
#ifdef CMPLX_
            call stats_out(state,.true., real(proje_iter_tot) + OutputHii, &
                           'Tot. Proj. E')
            call stats_out(state,.false.,allDoubleSpawns,'Double spawns')
#else
            call stats_out(state,.true., proje_iter_tot + OutputHii, &
                           'Tot. Proj. E')
#endif
            call stats_out(state,.true., AllTotWalkers, 'Dets occ.')
            call stats_out(state,.true., nspawned_tot, 'Dets spawned')
            call stats_out(state,.false., Hii, 'reference energy')

            if(t_real_time_fciqmc) then
                call stats_out(state,.false., real(sum(dyn_norm_red(:,1))/normsize),'GF normalization')
            else
                call stats_out(state,.true., IterTime, 'Iter. time')
            end if
            if(t_real_time_fciqmc) then
                call stats_out(state, .true., elapsedRealTime, 'Re. time')
                call stats_out(state, .true., elapsedImagTime, 'Im. time')
            else
                call stats_out(state,.false., TotImagTime, 'Im. time')
            end if

            ! Put the conditional columns at the end, so that the column
            ! numbers of the data are as stable as reasonably possible (for
            ! people who want to use gnuplot/not analyse column headers too
            ! frequently).
            ! This also makes column contiguity on resumes as likely as
            ! possible.
            if(t_real_time_fciqmc .or. tLogGreensfunction) then
                ! also output the overlaps and norm..
                do iGf = 1, gf_count
                   write(tmgf, '(i5)') iGf
                   call stats_out(state,.true., overlap_real(iGf), 'Re. <y_i(0)|y(t)> (i=' // &
                        trim(adjustl(tmgf)) // ')' )
                   call stats_out(state,.true., overlap_imag(iGf), 'Im. <y_i(0)|y(t)> (i=' // &
                        trim(adjustl(tmgf)) // ')' )
                end do
                do iGf = 1, gf_count
                   write(tmgf, '(i5)') iGf
                   do p = 1, normsize
                      write(tmpc, '(i5)') p
                      call stats_out(state,.false.,real(current_overlap(p,iGf)), 'Re. <y(0)|y(t)>(rep ' // &
                           trim(adjustl(tmpc)) // ' i=' // trim(adjustl(tmgf)) //  ')')
                      call stats_out(state,.false.,aimag(current_overlap(p,iGf)), 'Im. <y(0)|y(t)>(rep ' // &
                           trim(adjustl(tmpc)) // ' i=' // trim(adjustl(tmgf)) //')')
                   end do
                end do
                if(t_real_time_fciqmc) then
                    do p = 1, numSnapshotOrbs
                        ! if any orbitals are monitored, output their population
                        write(tmpc, '(i5)') snapShotOrbs(p)
                        call stats_out(state,.false.,allPopSnapshot(p),'Population of ' &
                            // trim(adjustl(tmpc)))
                    end do
                end if
            end if

            ! if we truncate walkers, print out the total truncated weight here
            if(t_truncate_spawns) call stats_out(state, .false., AllTruncatedWeight, &
                 'trunc. Weight')
            ! If we are running multiple (replica) simulations, then we
            ! want to record the details of each of these
#ifdef PROG_NUMRUNS_
            if(.not. t_real_time_fciqmc) then
                do p = 1, inum_runs
                    write(tmpc, '(i5)') p
                    call stats_out (state, .false., AllTotParts(p), &
                                    'Parts (' // trim(adjustl(tmpc)) // ')')
                    call stats_out (state, .false., AllNoatHF(p), &
                                    'Ref (' // trim(adjustl(tmpc)) // ')')
                    call stats_out(state, .false., proje_ref_energy_offsets(p), &
                                    'ref. energy offset('//trim(adjustl(tmpc))// ')')
                    call stats_out (state, .false., DiagSft(p) + Hii, &
                                    'Shift (' // trim(adjustl(tmpc)) // ')')
#ifdef CMPLX_
                    call stats_out (state, .false., real(proje_iter(p) + OutputHii), &
                                    'Tot ProjE real (' // trim(adjustl(tmpc)) // ')')
                    call stats_out (state, .false., aimag(proje_iter(p) + OutputHii), &
                                    'Tot ProjE imag (' // trim(adjustl(tmpc)) // ')')

                    call stats_out (state, .false., real(AllHFOut(p) / StepsPrint), &
                                    'ProjE Denom real (' // trim(adjustl(tmpc)) // ")")
                    call stats_out (state, .false., aimag(AllHFOut(p) / StepsPrint), &
                                    'ProjE Denom imag (' // trim(adjustl(tmpc)) // ")")

                    call stats_out (state, .false., &
                                    real((AllENumOut(p) + OutputHii*AllHFOut(p))) / StepsPrint,&
                                    'ProjE Num real (' // trim(adjustl(tmpc)) // ")")
                    call stats_out (state, .false., &
                                    aimag((AllENumOut(p) + OutputHii*AllHFOut(p))) / StepsPrint,&
                                    'ProjE Num imag (' // trim(adjustl(tmpc)) // ")")
                    if (tTrialWavefunction .or. tStartTrialLater) then
                        call stats_out (state, .false., &
                                        real(tot_trial_numerator(p) / StepsPrint), &
                                        'TrialE Num real (' // trim(adjustl(tmpc)) // ")")
                        call stats_out (state, .false., &
                                        aimag(tot_trial_numerator(p) / StepsPrint), &
                                        'TrialE Num imag (' // trim(adjustl(tmpc)) // ")")

                        call stats_out (state, .false., &
                                        real(tot_trial_denom(p) / StepsPrint), &
                                        'TrialE Denom real (' // trim(adjustl(tmpc)) // ")")
                        call stats_out (state, .false., &
                                        aimag(tot_trial_denom(p) / StepsPrint), &
                                        'TrialE Denom imag (' // trim(adjustl(tmpc)) // ")")
                    end if
#else
                    call stats_out (state, .false., proje_iter(p) + OutputHii, &
                                    'Tot ProjE (' // trim(adjustl(tmpc)) // ")")
                    call stats_out (state, .false., AllHFOut(p) / StepsPrint, &
                                    'ProjE Denom (' // trim(adjustl(tmpc)) // ")")
                    call stats_out (state, .false., &
                                    (AllENumOut(p) + OutputHii*AllHFOut(p)) / StepsPrint,&
                                    'ProjE Num (' // trim(adjustl(tmpc)) // ")")
                    if (tTrialWavefunction .or. tStartTrialLater) then
                        call stats_out (state, .false., &
                                        tot_trial_numerator(p) / StepsPrint, &
                                        'TrialE Num (' // trim(adjustl(tmpc)) // ")")
                        call stats_out (state, .false., &
                                        tot_trial_denom(p) / StepsPrint, &
                                        'TrialE Denom (' // trim(adjustl(tmpc)) // ")")
                    end if
#endif


                    call stats_out (state, .false., &
                                    AllNoBorn(p), &
                                    'Born (' // trim(adjustl(tmpc)) // ')')
                    call stats_out (state, .false., &
                                    AllNoDied(p), &
                                    'Died (' // trim(adjustl(tmpc)) // ')')
                    call stats_out (state, .false., &
                                    AllAnnihilated(p), &
                                    'Annihil (' // trim(adjustl(tmpc)) // ')')
                    call stats_out (state, .false., &
                                    AllNoAtDoubs(p), &
                                    'Doubs (' // trim(adjustl(tmpc)) // ')')
                end do

                call stats_out(state,.false.,all_max_cyc_spawn, &
                     'MaxCycSpawn')

                ! Print overlaps between replicas at the end.
                do p = 1, inum_runs
                    write(tmpc, '(i5)') p
                        if (tPrintReplicaOverlaps) then
                        do q = p+1, inum_runs
                            write(tmpc2, '(i5)') q
#ifdef CMPLX_
                            call stats_out(state, .false.,  replica_overlaps_real(p, q),&
                                           '<psi_' // trim(adjustl(tmpc)) // '|' &
                                           // 'psi_' // trim(adjustl(tmpc2)) &
                                           // '> (real)')
                            call stats_out(state, .false.,  replica_overlaps_imag(p, q),&
                                           '<psi_' // trim(adjustl(tmpc)) // '|' &
                                           // 'psi_' // trim(adjustl(tmpc2)) &
                                           // '> (imag)')

#else
                            call stats_out(state, .false.,  replica_overlaps_real(p, q),&
                                 '<psi_' // trim(adjustl(tmpc)) // '|' &
                                 // 'psi_' // trim(adjustl(tmpc2)) &
                                 // '>')
#endif

                        end do
                    end if
                end do
            end if
#endif
            if (tEN2) call stats_out(state,.true., en_pert_main%ndets_all, 'EN2 Dets.')

            if (tTruncInitiator) then
                call stats_out(state_i, .false., Iter + PreviousCycles, 'Iter.')
                call stats_out(state_i, .false., AllTotWalkers, 'TotDets.')
                do p = 1, inum_runs
                    write(tmpc, '(i5)') p
                    call stats_out(state_i, .false., AllTotParts(p), 'TotWalk. (' // trim(adjustl(tmpc)) // ")")
                    call stats_out(state_i, .false., AllAnnihilated(p), 'Annihil. (' // trim(adjustl(tmpc)) // ")")
                    call stats_out(state_i, .false., AllNoBorn(p), 'Born (' // trim(adjustl(tmpc)) // ")")
                    call stats_out(state_i, .false., AllNoDied(p), 'Died (' // trim(adjustl(tmpc)) // ")")
                    call stats_out(state_i, .false., AllNoRemoved(p), 'Removed Dets (' // trim(adjustl(tmpc)) // ")")
                    call stats_out(state_i, .false., AllNoAborted(p), 'AbortedWalks (' // trim(adjustl(tmpc)) // ")")
                    call stats_out(state_i, .false., AllNoInitDets(p), 'InitDets (' // trim(adjustl(tmpc)) // ")")
                    call stats_out(state_i, .false., AllNoNonInitDets(p), 'NonInitDets (' // trim(adjustl(tmpc)) // ")")
                    call stats_out(state_i, .false., AllNoInitWalk(p), 'InitWalks (' // trim(adjustl(tmpc)) // ")")
                    call stats_out(state_i, .false., AllNoNonInitWalk(p), 'NonInitWalks (' // trim(adjustl(tmpc)) // ")")
                end do
            end if
            ! And we are done
            write(state%funit, *)
            if (tTruncInitiator) write(state_i%funit, *)
            if (tMCOutput) write(stdout, *)

            call neci_flush(state%funit)
            if (tTruncInitiator) call neci_flush(state_i%funit)
            call neci_flush(stdout)

        end if

    end subroutine write_fcimcstats2

    subroutine open_state_file(filename,state)
      implicit none
      character(*), intent(in) :: filename
      type(write_state_t), intent(inout) :: state
      ! mini-subroutine for opening a file and assigning it to a state
      state%funit = get_free_unit()
      call open_create_stats(filename,state%funit)

    end subroutine open_state_file

    subroutine write_padding_init(state)
      implicit none
      type(write_state_t), intent(inout) :: state

      ! Don't treat the header line as data. Add padding to align the
      ! other columns. We also add a # to the first line of data, so
      ! that there aren't repeats if starting from POPSFILES
      if (state%init .or. state%prepend) then
         write(state%funit, '("#")', advance='no')
         if (tMCOutput) write(stdout, '("#")', advance='no')
         state%prepend = state%init
      else if (.not. state%prepend) then
         write(state%funit, '(" ")', advance='no')
         if (tMCOutput) write(stdout, '(" ")', advance='no')
      end if
    end subroutine write_padding_init

!Similar to WriteHistogram, but will only print out in order of maximum component, and only the averaged wavefunction
    SUBROUTINE PrintFCIMCPsi()
        use DetCalcData , only : FCIDets
        INTEGER :: i,nI(NEl),ExcitLevel,j, iunit
        real(dp) :: norm2
        real(dp), dimension(lenof_sign) :: norm1, norm

        CALL MPISumAll(Histogram,AllHistogram)
        norm1=0.0_dp
        do i=1,Det
            do j=1,lenof_sign
                norm1(j)=norm1(j)+AllHistogram(j,i)**2
            end do
        end do
#ifdef CMPLX_
        norm2=SQRT(sum(norm1))
#else
        norm1=SQRT(norm1)
#endif
        write(stdout,*) "Total FCIMC Wavefuction normalisation:",norm1
        do i=1,Det
            do j=1,lenof_sign
#ifdef CMPLX_
                AllHistogram(j,i)=AllHistogram(j,i)/norm2
#else
                AllHistogram(j,i)=AllHistogram(j,i)/norm1(j)
#endif
            end do
        end do

        iunit = 0
        IF(tPrintFCIMCPsi) THEN
!Order and print wavefunction


            IF(iProcIndex.eq.0) THEN

                ! We now want to order AllHistogram, taking the corresponding
                ! element(s) of FCIDets with it...
                call sort (AllHistogram, FCIDets)

                open(iunit,FILE='FCIMCPsi',STATUS='UNKNOWN')

                norm=0.0_dp
                do i=1,Det
                    do j=1,lenof_sign
                        norm(j)=norm(j)+AllHistogram(j,i)**2
                    end do
!write out FCIMC Component weight (normalised), current normalisation, excitation level
                    ExcitLevel = FindBitExcitLevel(iLutHF, FCIDets(:,i), nel)
                    CALL decode_bit_det(nI,FCIDets(0:NIfTot,i))
#ifdef CMPLX_
                    write(iunit,"(I13,G25.16,I6,G20.10)",advance='no') i,AllHistogram(1,i),ExcitLevel,sum(norm)
#else
                    write(iunit,"(I13,G25.16,I6,G20.10)",advance='no') i,AllHistogram(1,i),ExcitLevel,norm(1)
#endif
                    do j=1,NEl-1
                        write(iunit,"(I5)",advance='no') nI(j)
                    end do
                    write(iunit,"(I5)") nI(NEl)
                end do

                close(iunit)

            end if
        end if

    END SUBROUTINE PrintFCIMCPsi


!This routine will write out the average wavevector from the spawning run up until now.
    SUBROUTINE WriteHistogram()
        INTEGER :: i,IterRead, io1, io2, io3, Tot_No_Unique_Dets,ierr,val,j
        integer :: posn, FinalPop
        real(dp) :: ShiftRead,AllERead,NumParts,DDOT
        real(dp),dimension(lenof_sign) :: norm,norm1,norm2,norm3
        real(dp) :: norm_c,norm1_c,norm2_c,norm3_c
        real(dp) :: AvVarEnergy, VarEnergy, GroundE_Ever, ProjGroundE
        real(dp) , allocatable :: CKN(:),HOrderedHist(:),HOrderedInstHist(:)
        integer , allocatable :: ExpandedWalkerDets(:,:)
        CHARACTER(len=22) :: abstr,abstr2
        LOGICAL :: exists
        character(len=*), parameter :: t_r='WriteHistogram'

!This will open a file called SpawnHist-"Iter" on unit number 17.
        abstr = 'SpawnHist-'//str(Iter)
        IF(iProcIndex.eq.0) THEN
            write(stdout,*) "Writing out the average wavevector up to iteration number: ", Iter
            CALL neci_flush(stdout)
        end if

        IF(iProcIndex.eq.0) THEN
            AllHistogram(:,:)=0.0_dp
            AllInstHist(:,:)=0.0_dp
            AllAvAnnihil(:,:)=0.0_dp
            AllInstAnnihil(:,:)=0.0_dp
        end if

        CALL MPIReduce(Histogram,MPI_SUM,AllHistogram)
        CALL MPIReduce(InstHist,MPI_SUM,AllInstHist)
        CALL MPIReduce(InstAnnihil,MPI_SUM,AllInstAnnihil)
        CALL MPIReduce(AvAnnihil,MPI_SUM,AllAvAnnihil)

        BeforeNormHist(:) = AllHistogram(1,:)

        IF(iProcIndex.eq.0) THEN

            norm=0.0_dp
            norm1=0.0_dp
            norm2=0.0_dp
            norm3=0.0_dp
            do i=1,Det
                do j=1,lenof_sign
                    norm(j)=norm(j)+AllHistogram(j,i)**2
                    norm1(j)=norm1(j)+AllInstHist(j,i)**2
                    norm2(j)=norm2(j)+AllInstAnnihil(j,i)**2
                    norm3(j)=norm3(j)+AllAvAnnihil(j,i)**2
                end do
            end do
#ifdef CMPLX_
            norm_c=SQRT(sum(norm))
            norm1_c=SQRT(sum(norm1))
            norm2_c=SQRT(sum(norm2))
            norm3_c=SQRT(sum(norm3))
#else
            norm=SQRT(norm)
            norm1=SQRT(norm1)
            norm2=SQRT(norm2)
            norm3=SQRT(norm3)
#endif
            do i=1,Det
                do j=1,lenof_sign
#ifdef CMPLX_
                    AllHistogram(j,i)=AllHistogram(j,i)/norm_c
                    AllInstHist(j,i)=AllInstHist(j,i)/norm1_c
                    IF(.not. near_zero(norm2_c)) THEN
                        AllInstAnnihil(j,i)=AllInstAnnihil(j,i)/norm2_c
                    end if
                    IF(.not. near_zero(norm3_c)) THEN
                        AllAvAnnihil(j,i)=AllAvAnnihil(j,i)/norm3_c
                    end if
#else
                    AllHistogram(j,i)=AllHistogram(j,i)/norm(j)
                    AllInstHist(j,i)=AllInstHist(j,i)/norm1(j)
                    IF(.not. near_zero(norm2(j))) THEN
                        AllInstAnnihil(j,i)=AllInstAnnihil(1,i)/norm2(j)
                    end if
                    IF(.not. near_zero(norm3(j))) THEN
                    AllAvAnnihil(j,i)=AllAvAnnihil(j,i)/norm3(j)
                    end if
#endif
                end do
            end do

            io1 = get_free_unit()
            open(io1,FILE=abstr,STATUS='UNKNOWN')

            abstr = 'Energies-'//str(Iter - iWriteHistEvery)
            abstr2 = 'Energies-'//str(Iter)

            io2 = get_free_unit()
            open(io2,FILE=abstr2,STATUS='UNKNOWN')

            INQUIRE(FILE=abstr,EXIST=exists)
            IF(exists) THEN
                io3 = get_free_unit()
                open(io3,FILE=abstr,STATUS='OLD',POSITION='REWIND',ACTION='READ')
                do while(.true.)
                    read(io3,"(I13,3G25.16)",END=99) IterRead,ShiftRead,AllERead,NumParts
                    write(io2,"(I13,3G25.16)") IterRead,ShiftRead,AllERead,NumParts
                end do
99              CONTINUE
#ifdef CMPLX_
                IF(near_zero(AllHFOut(1))) then
                    write(io2,"(I13,3G25.16)") Iter,DiagSft,AllERead,SUM(AllTotPartsLastOutput)
                ELSE
                    write(io2,"(I13,3G25.16)") Iter,DiagSft,AllENumOut/AllHFOut,SUM(AllTotPartsLastOutput)
                end if
#else
                IF(near_zero(AllHFOut(1))) THEN
                    write(io2,"(I13,3G25.16)") Iter,DiagSft(1),AllERead,AllTotPartsLastOutput(1)
                ELSE
                    write(io2,"(I13,3G25.16)") Iter,DiagSft(1),AllENumOut(1)/AllHFOut(1),AllTotPartsLastOutput(1)
                end if
#endif
                close(io2)
                close(io3)

            ELSE
                open(io2,FILE=abstr2,STATUS='UNKNOWN')
#ifdef CMPLX_
                write(io2,"(I13,3G25.16)") Iter,DiagSft,AllENumOut/AllHFOut,SUM(AllTotPartsLastOutput)
#else
                write(io2,"(I13,3G25.16)") Iter,DiagSft(1),AllENumOut(1)/AllHFOut(1),AllTotPartsLastOutput(1)
#endif
                close(io2)
            end if


            norm=0.0_dp
            norm1=0.0_dp
            Tot_No_Unique_Dets = 0
            do i=1,Det

                posn=binary_search_ilut(CurrentDets(:,1:TotWalkers), FCIDETS(:,i), NifD+1)

                if (posn.lt.0) then
                    FinalPop = 0
                else
                    FinalPop = int(CurrentDets(NifD+1,posn))
                end if

                do j=1,lenof_sign
                    norm(j)=norm(j)+AllHistogram(j,i)**2
                    norm1(j)=norm1(j)+AllAvAnnihil(j,i)**2
                end do
                IF(lenof_sign.eq.1) THEN
                    write(io1,"(I13,6G25.16,I13,G25.16)") i, &
                          AllHistogram(1,i), norm, AllInstHist(1,i), &
                          AllInstAnnihil(1,i), AllAvAnnihil(1,i), norm1, &
                          FinalPop, BeforeNormHist(i)
                ELSE
#ifdef CMPLX_
                    write(io1,"(I13,6G25.16)") i, AllHistogram(1,i), sum(norm), &
                          AllInstHist(1,i), AllInstAnnihil(1,i), &
                          AllAvAnnihil(1,i), sum(norm1)
#else
                    write(io1,"(I13,6G25.16)") i, AllHistogram(1,i), norm(1), &
                          AllInstHist(1,i), AllInstAnnihil(1,i), &
                          AllAvAnnihil(1,i), norm1(1)
#endif
                end if
                IF(.not. near_zero(AllHistogram(1,i))) Tot_No_Unique_Dets = Tot_No_Unique_Dets + 1
            end do
            if(tCalcVariationalEnergy) then
                !Calculate the variational FCIMC energy

                !Check that Det eq NDet
                if(Det.ne.NDet) call stop_all(t_r,'Error')
                allocate(CKN(1:Det))
                CKN = 0.0_dp
                !We need to reorder the vectors so that they are in the same order as the hamiltonian
                allocate(HOrderedHist(1:Det))
                HOrderedHist(:) = AllHistogram(1,:)
                allocate(HOrderedInstHist(1:Det))
                HOrderedInstHist(:) = AllInstHist(1,:)
                call sort(ReIndex(1:Det),HOrderedHist(1:Det),HOrderedInstHist(1:Det))

#ifdef CMPLX_
                call stop_all(t_r, "does not work for complex")
#else
                call my_hpsi(Det,1,NROW,LAB,HAMIL,HOrderedHist,CKN,.true.)
#endif
                AvVarEnergy = DDOT(Det,HOrderedHist,1,CKN,1)

                CKN = 0.0_dp
#ifdef CMPLX_
                call stop_all(t_r, "does not work for complex")
#else
                call my_hpsi(Det,1,NROW,LAB,HAMIL,HOrderedInstHist,CKN,.true.)
#endif
                VarEnergy = DDOT(Det,HOrderedInstHist,1,CKN,1)

                deallocate(CKN)
                deallocate(HOrderedHist,HOrderedInstHist)
                if(tDiagAllSpaceEver) then

                    allocate(ExpandedWalkerDets(NEl,Tot_No_Unique_Dets),stat=ierr)
                    val=1
                    do i=1,Det
                        if(.not. near_zero(AllHistogram(1,i))) then
                            call decode_bit_det(ExpandedWalkerDets(:,val),FCIDets(0:NIfTot,i))
                            val=val+1
                        end if
                    end do
                    if((val-1).ne.Tot_No_Unique_Dets) call stop_all(t_r,'Wrong counting')

                    call LanczosFindGroundE(ExpandedWalkerDets,Tot_No_Unique_Dets,GroundE_Ever,ProjGroundE,.true.)
                    if(abs(GroundE_Ever-ProjGroundE).gt.1.0e-7_dp) call stop_all(t_r,'Why do these not agree?!')
                    write(Tot_Unique_Dets_Unit,"(2I14,3G25.15)") Iter,Tot_No_Unique_Dets,AvVarEnergy,VarEnergy,GroundE_Ever
                else
                    write(Tot_Unique_Dets_Unit,"(2I14,2G25.15)") Iter,Tot_No_Unique_Dets,AvVarEnergy,VarEnergy
                end if
            else
                if(tDiagAllSpaceEver) then

                    allocate(ExpandedWalkerDets(NEl,Tot_No_Unique_Dets),stat=ierr)
                    val=1
                    do i=1,Det
                        if(.not. near_zero(AllHistogram(1,i))) then
                            call decode_bit_det(ExpandedWalkerDets(:,val),FCIDets(0:NIfTot,i))
                            val=val+1
                        end if
                    end do
                    if((val-1).ne.Tot_No_Unique_Dets) call stop_all(t_r,'Wrong counting')

                    call LanczosFindGroundE(ExpandedWalkerDets,Tot_No_Unique_Dets,GroundE_Ever,ProjGroundE,.true.)
                    if(abs(GroundE_Ever-ProjGroundE).gt.1.0e-7_dp) call stop_all(t_r,'Why do these not agree?!')
                    write(Tot_Unique_Dets_Unit,"(2I14,3G25.15)") Iter,Tot_No_Unique_Dets,GroundE_Ever
                else
                    write(Tot_Unique_Dets_Unit,"(2I14)") Iter, Tot_No_Unique_Dets
                end if
            end if

            close(io1)
        end if
        InstHist(:,:)=0.0_dp
        InstAnnihil(:,:)=0.0_dp

    END SUBROUTINE WriteHistogram

    ! TODO: Move to hist.F90
    SUBROUTINE WriteHistogramEnergies()
        INTEGER :: i, io(8)
        real(dp) :: Norm,EnergyBin

        IF(iProcIndex.eq.Root) THEN
            AllHistogramEnergy(:)=0.0_dp
            AllAttemptHist(:)=0.0_dp
            AllSpawnHist(:)=0.0_dp
            AllDoublesHist(:)=0.0_dp
            AllDoublesAttemptHist(:)=0.0_dp
            AllSinglesHist(:)=0.0_dp
            AllSinglesAttemptHist(:)=0.0_dp
            AllSinglesHistOccOcc(:)=0.0_dp
            AllSinglesHistOccVirt(:)=0.0_dp
            AllSinglesHistVirtOcc(:)=0.0_dp
            AllSinglesHistVirtVirt(:)=0.0_dp
        end if
        CALL MPIReduce(HistogramEnergy,MPI_SUM,AllHistogramEnergy)
        CALL MPIReduce(AttemptHist,MPI_SUM,AllAttemptHist)
        CALL MPIReduce(SpawnHist,MPI_SUM,AllSpawnHist)
        CALL MPIReduce(SinglesHist,MPI_SUM,AllSinglesHist)
        CALL MPIReduce(SinglesAttemptHist,MPI_SUM,AllSinglesAttemptHist)
        CALL MPIReduce(DoublesHist,MPI_SUM,AllDoublesHist)
        CALL MPIReduce(DoublesAttemptHist,MPI_SUM,AllDoublesAttemptHist)
        CALL MPIReduce(SinglesHistOccOcc,MPI_SUM,AllSinglesHistOccOcc)
        CALL MPIReduce(SinglesHistOccVirt,MPI_SUM,AllSinglesHistOccVirt)
        CALL MPIReduce(SinglesHistVirtOcc,MPI_SUM,AllSinglesHistVirtOcc)
        CALL MPIReduce(SinglesHistVirtVirt,MPI_SUM,AllSinglesHistVirtVirt)

        IF(iProcIndex.eq.Root) THEN
            AllHistogramEnergy=AllHistogramEnergy/sum(AllHistogramEnergy)
            AllAttemptHist=AllAttemptHist/sum(AllAttemptHist)
            AllSpawnHist=AllSpawnHist/sum(AllSpawnHist)
            AllSinglesAttemptHist=AllSinglesAttemptHist/sum(AllSinglesAttemptHist)
            AllDoublesHist=AllDoublesHist/sum(AllDoublesHist)
            Norm=sum(AllDoublesAttemptHist)
            AllDoublesAttemptHist=AllDoublesAttemptHist/Norm
            do i=1,iOffDiagNoBins
                AllDoublesAttemptHist(i)=AllDoublesAttemptHist(i)/Norm
            end do
            Norm=0.0_dp
            do i=1,iOffDiagNoBins
                Norm=Norm+AllSinglesHist(i)
            end do
!            write(stdout,*) "AllSinglesHistNorm = ",Norm
            do i=1,iOffDiagNoBins
                AllSinglesHist(i)=AllSinglesHist(i)/Norm
            end do

!            Norm=0.0_dp
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistOccOcc(i)
!            end do
            do i=1,iOffDiagNoBins
                AllSinglesHistOccOcc(i)=AllSinglesHistOccOcc(i)/Norm
            end do
!            Norm=0.0_dp
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistOccVirt(i)
!            end do
            do i=1,iOffDiagNoBins
                AllSinglesHistOccVirt(i)=AllSinglesHistOccVirt(i)/Norm
            end do
!            Norm=0.0_dp
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistVirtOcc(i)
!            end do
            do i=1,iOffDiagNoBins
                AllSinglesHistVirtOcc(i)=AllSinglesHistVirtOcc(i)/Norm
            end do
!            Norm=0.0_dp
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistVirtVirt(i)
!            end do
            do i=1,iOffDiagNoBins
                AllSinglesHistVirtVirt(i)=AllSinglesHistVirtVirt(i)/Norm
            end do


            io(1) = get_free_unit()
            open(io(1),FILE='EVERYENERGYHIST',STATUS='UNKNOWN')
            io(2) = get_free_unit()
            open(io(2),FILE='ATTEMPTENERGYHIST',STATUS='UNKNOWN')
            io(3) = get_free_unit()
            open(io(3),FILE='SPAWNENERGYHIST',STATUS='UNKNOWN')

            EnergyBin=BinRange/2.0_dp
            do i=1,iNoBins
                IF(AllHistogramEnergy(i).gt.0.0_dp) write(io(1),*) EnergyBin, AllHistogramEnergy(i)
                IF(AllAttemptHist(i).gt.0.0_dp) write(io(2),*) EnergyBin, AllAttemptHist(i)
                IF(AllSpawnHist(i).gt.0.0_dp) write(io(3),*) EnergyBin, AllSpawnHist(i)
                EnergyBin=EnergyBin+BinRange
            end do
            close(io(1))
            close(io(2))
            close(io(3))
            open(io(1),FILE='SINGLESHIST',STATUS='UNKNOWN')
            open(io(2),FILE='ATTEMPTSINGLESHIST',STATUS='UNKNOWN')
            open(io(3),FILE='DOUBLESHIST',STATUS='UNKNOWN')
            io(4) = get_free_unit()
            open(io(4),FILE='ATTEMPTDOUBLESHIST',STATUS='UNKNOWN')
            io(5) = get_free_unit()
            open(io(5),FILE='SINGLESHISTOCCOCC',STATUS='UNKNOWN')
            io(6) = get_free_unit()
            open(io(6),FILE='SINGLESHISTOCCVIRT',STATUS='UNKNOWN')
            io(7) = get_free_unit()
            open(io(7),FILE='SINGLESHISTVIRTOCC',STATUS='UNKNOWN')
            io(8) = get_free_unit()
            open(io(8),FILE='SINGLESHISTVIRTVIRT',STATUS='UNKNOWN')

            EnergyBin=-OffDiagMax+OffDiagBinRange/2.0_dp
            do i=1,iOffDiagNoBins
                IF(AllSinglesHist(i).gt.0.0_dp) write(io(1),*) EnergyBin, AllSinglesHist(i)
                IF(AllSinglesAttemptHist(i).gt.0.0_dp) write(io(2),*) EnergyBin, AllSinglesAttemptHist(i)
                IF(AllDoublesHist(i).gt.0.0_dp) write(io(3),*) EnergyBin, AllDoublesHist(i)
                IF(AllDoublesAttemptHist(i).gt.0.0_dp) write(io(4),*) EnergyBin, AllDoublesAttemptHist(i)
                IF(AllSinglesHistOccOcc(i).gt.0.0_dp) write(io(5),*) EnergyBin, AllSinglesHistOccOcc(i)
                IF(AllSinglesHistOccVirt(i).gt.0.0_dp) write(io(6),*) EnergyBin, AllSinglesHistOccVirt(i)
                IF(AllSinglesHistVirtOcc(i).gt.0.0_dp) write(io(7),*) EnergyBin, AllSinglesHistVirtOcc(i)
                IF(AllSinglesHistVirtVirt(i).gt.0.0_dp) write(io(8),*) EnergyBin, AllSinglesHistVirtVirt(i)
                EnergyBin=EnergyBin+OffDiagBinRange
!                write(stdout,*) i
            end do

            close(io(1))
            close(io(2))
            close(io(3))
            close(io(4))
            close(io(5))
            close(io(6))
            close(io(7))
            close(io(8))
        end if

    END SUBROUTINE WriteHistogramEnergies

    !Routine to print the highest populated determinants at the end of a run
    SUBROUTINE PrintHighPops()
        use adi_references, only: update_ref_signs, print_reference_notification, nRefs
        use adi_data, only: tSetupSIs
        use guga_data, only: ExcitationInformation_t
        use guga_bitrepops, only: identify_excitation, write_guga_list
        real(dp), dimension(lenof_sign) :: SignCurr
        integer :: ierr,i,j,counter,ExcitLev,nopen
        integer :: full_orb, run
        real(dp) :: HighSign, norm
        integer(n_int) , allocatable :: GlobalLargestWalkers(:, :)
        integer, allocatable :: GlobalProc(:), tmp_ni(:)
        real(dp), allocatable :: GlobalHdiag(:)
        character(100) :: bufEnd, bufStart
        integer :: lenEnd, lenStart
        character(len=*), parameter :: t_r='PrintHighPops'

        character(1024) :: header
        character(11), allocatable :: walker_string(:)
        character(13), allocatable :: amplitude_string(:)
        character(9), allocatable :: init_string(:)

        integer :: lenof_out, this_run, offset
        logical :: t_replica_resolved_output

        real(dp), parameter :: eps_high = 1.0e-7_dp

        type(ExcitationInformation_t) :: excitInfo

        allocate(GlobalLargestWalkers(0:NIfTot,iHighPopWrite), source=0_n_int)
        allocate(GlobalHdiag(iHighPopWrite), source=0.0_dp)
        allocate(GlobalProc(iHighPopWrite), source=0)

        ! Decide if each replica shall have its own output
        t_replica_resolved_output = tOrthogonaliseReplicas .or. t_force_replica_output
        if(t_replica_resolved_output) then
            lenof_out = rep_size
        else
            lenof_out = lenof_sign
        end if

        ! Walkers(replica) Amplitude(replica) Init?(replica)
        allocate(walker_string(lenof_out))
        allocate(init_string(lenof_out))

        do run = 1, inum_runs
            ! If t_replica_resolved_output is set:
            ! Execute this once per run with run instead of GLOBAL_RUN -> prints the highest
            ! determinants for each replica
            if(t_replica_resolved_output) then
                write(stdout,*) "============================================================="
                write(stdout,*) "Reference and leading determinants for replica",run
                write(stdout,*) "============================================================="
            end if
            if(t_replica_resolved_output) then
                this_run = run
                offset = rep_size * (run - 1)
            else
                this_run = GLOBAL_RUN
                offset = 0
            end if
            call global_most_populated_states(iHighPopWrite, this_run, GlobalLargestWalkers, &
                norm, rank_of_largest=GlobalProc, hdiag_largest=GlobalHdiag)

            ! This has to be done by all procs
            if(tAdiActive) call update_ref_signs()

            if(iProcIndex.eq.Root) then
                !Now print out the info contained in GlobalLargestWalkers and GlobalProc

                counter=0
                do i=1,iHighPopWrite
                    !How many non-zero determinants do we actually have?
                    call extract_sign(GlobalLargestWalkers(:,i),SignCurr)
                    HighSign = core_space_weight(SignCurr,this_run)
                    if (HighSign > eps_high) counter = counter + 1
                end do


                write(stdout,*) ""
                if (tReplicaReferencesDiffer) then
                    write(stdout,'(A)') "Current references: "
                    call write_det(stdout, ProjEDet(:,run), .true.)
                    call writeDetBit(stdout, ilutRef(:, run), .true.)
                else
                    write(stdout,'(A)') "Current reference: "
                    call write_det (stdout, ProjEDet(:,1), .true.)
                    if(tSetupSIs) call print_reference_notification(&
                        1,nRefs,"Used Superinitiator",.true.)
                    write(stdout,*) "Number of superinitiators", nRefs
                end if

                write(stdout,*)
                write(stdout,'("Input DEFINEDET line (includes frozen orbs):")')
                write(stdout,'("definedet ")', advance='no')
                if (allocated(frozen_orb_list)) then
                    allocate(tmp_ni(nel_pre_freezing))
                    tmp_ni(1:nel) = frozen_orb_reverse_map(ProjEDet(:,run))
                    if (nel /= nel_pre_freezing) &
                        tmp_ni(nel+1:nel_pre_freezing) = frozen_orb_list
                    call sort(tmp_ni)
                    call writeDefDet(tmp_ni, nel_pre_freezing)
                    !                    do i = 1, nel_pre_freezing
                    !                        write(stdout, '(i3," ")', advance='no') tmp_ni(i)
                    !                    end do
                    deallocate(tmp_ni)
                else
                    call writeDefDet(ProjEDet(:,run), nel)
                    !                    do i = 1, nel
                    !                        write(stdout, '(i3," ")', advance='no') ProjEDet(i, run)
                    !                    end do
                end if
                do i = 1, nel
                    full_orb = ProjEDet(i, run)
                    if (allocated(frozen_orb_list)) &
                        full_orb = full_orb  + count(frozen_orb_list <= ProjEDet(i, run))
                end do
                write(stdout,*)

                write(stdout,*) ""
                write(stdout,"(A,I10,A)") "Most occupied ",counter," determinants as excitations from reference: "
                write(stdout,*)
                if(lenof_sign.eq.1) then
                    if(tHPHF) then
                        write(stdout,"(A)") " Excitation   ExcitLevel   Seniority    Walkers    Amplitude    Init?   <D|H|D>  Proc  Spin-Coup?"
                    else
                        write(stdout,"(A)") " Excitation   ExcitLevel   Seniority    Walkers    Amplitude    Init?   <D|H|D>  Proc"
                    end if
                else
#ifdef CMPLX_
                    if(tHPHF) then
                        write(stdout,"(A)") " Excitation   ExcitLevel Seniority  Walkers(Re)   Walkers(Im)  Weight   &
                            &Init?(Re)   Init?(Im)   <D|H|D>  Proc  Spin-Coup?"
                    else
                        write(stdout,"(A)") " Excitation   ExcitLevel Seniority   Walkers(Re)   Walkers(Im)  Weight   &
                            &Init?(Re)   Init?(Im)   <D|H|D>  Proc"
                    end if
#else
                    ! output the weight of every replica, and do not only assume
                    ! it is a complex run

                    do i = 1, lenof_out
                        write(walker_string(i), '(a,i0,a)') "Walkers(", i, ")"
                        !                     write(amplitude_string(i), '(a,i0,a)') "Amplitude(", i, ")"
                        write(init_string(i), '(a,i0,a)') "Init?(", i, ")"
                    end do

                    block
                        character(:), allocatable :: fmt_str
                        fmt_str = '(3a11,' // str(lenof_out) // 'a11, a13,' // str(lenof_out) // 'a9,1x,16a,1x,a)'
                        write(header, fmt_str) "Excitation ", "ExcitLevel ", "Seniority ", &
                            walker_string, "Amplitude ", init_string, "<D|H|D>", "Proc "
                    end block

                    if (tHPHF) then
                        header = trim(header) // " Spin-Coup?"
                    end if

                    write(stdout, '(a)') trim(header)

#endif
                end if
                do i=1,iHighPopWrite
                    call extract_sign(GlobalLargestWalkers(:,i),SignCurr)
                    HighSign = core_space_weight(SignCurr,this_run)
                    if(HighSign < eps_high) cycle
                    call WriteDetBit(stdout,GlobalLargestWalkers(:,i),.false.)
                    if (tGUGA) then
                        excitInfo = identify_excitation(iLutRef(:,1), GlobalLargestWalkers(:,i))
                        excitLev = excitInfo%excitLvl
                    else
                        Excitlev=FindBitExcitLevel(iLutRef(:,1),GlobalLargestWalkers(:,i),nEl,.true.)
                    end if
                    write(stdout,"(I5)",advance='no') Excitlev
                    nopen=count_open_orbs(GlobalLargestWalkers(:,i))
                    write(stdout,"(I5)",advance='no') nopen
                    do j=1,lenof_out
                        write(stdout,"(G16.7)",advance='no') SignCurr(j+offset)
                    end do
                    if(tHPHF.and.(.not.TestClosedShellDet(GlobalLargestWalkers(:,i)))) then
                        !Weight is proportional to (nw/sqrt(2))**2
                        write(stdout,"(F9.5)",advance='no') ((HighSign/sqrt(2.0_dp))/norm )
                    else
                        write(stdout,"(F9.5)",advance='no') (HighSign/norm)
                    end if
                    do j=1,lenof_out
                        if(.not.tTruncInitiator) then
                            write(stdout,"(A3)",advance='no') 'Y'
                        else
                            if(test_flag(GlobalLargestWalkers(:,i),get_initiator_flag(j+offset))) then
                                write(stdout,"(A3)",advance='no') 'Y'
                            else
                                write(stdout,"(A3)",advance='no') 'N'
                            end if
                        end if
                    end do
                    write(stdout,"(1x,es16.8,1x)",advance='no') GlobalHdiag(i)
                    if(tHPHF.and.(.not.TestClosedShellDet(GlobalLargestWalkers(:,i)))) then
                        write(stdout,"(I7)",advance='no') GlobalProc(i)
                        write(stdout,"(A3)") "*"
                    else
                        write(stdout,"(I7)") GlobalProc(i)
                    end if
                end do
                ! Keep the reference weight in a separate output variable
                ! which can be accessed from the library wrappers
                call extract_sign(GlobalLargestWalkers(:,1),SignCurr)
                fciqmc_run_ref_weight = SignCurr(1)
                if(tHPHF) then
                    write(stdout,"(A)") " * = Spin-coupled function implicitly has time-reversed determinant with same weight."
                end if

                write(stdout,*) ""
            end if
            ! Only continue if printing the replica-resolved output
            if(.not. t_replica_resolved_output) exit
        end do

        deallocate(GlobalLargestWalkers,GlobalProc,GlobalHdiag)
        deallocate(walker_string, init_string)

        contains

            subroutine writeDefDet(defdet, numEls)
                implicit none
                integer, intent(in) :: numEls
                integer, intent(in) :: defdet(:)
                logical :: nextInRange, previousInRange

                do i = 1, numEls
                    ! if the previous orbital is in the same contiguous range
                    if (i == 1) then
                        ! for the first one, there is no previous one
                        previousInRange = .false.
                    else
                        previousInRange = defdet(i) == defdet(i-1) + 1
                    end if

                    ! if the following orbital is in the same contiguous range
                    if(i.eq.numEls) then
                        ! there is no following orbital
                        nextInRange = .false.
                    else
                        nextInRange = defdet(i).eq.defdet(i+1)-1
                    end if
                    ! there are three cases that need output:

                    ! the last orbital of a contigous range of orbs
                    if(previousInRange .and. .not.nextInRange) then
                        write(bufEnd,'(i3)') defdet(i)
                        lenEnd = len_trim(bufEnd)
                        bufStart(lenStart+2:lenStart+lenEnd+1) = adjustl(trim(bufEnd))
                        write(stdout,'(A7)',advance='no') trim(adjustl(bufStart))
                        ! the first orbital of a contiguous range of orbs
                    else if(.not.previousInRange .and. nextInRange) then
                        write(bufStart,'(i3)') defdet(i)
                        lenStart = len_trim(bufStart)
                        bufStart(lenStart+1:lenStart+1) = "-"
                        ! and an orbital not in any range
                    else if(.not.previousInRange .and. .not.nextInRange) then
                        write(stdout,'(i3," ")', advance='no') defdet(i)
                    end if
                end do

            end subroutine writeDefDet

    END SUBROUTINE PrintHighPops

!------------------------------------------------------------------------------------------!

    !> Print out an already genereated 2d-histogram to disk
    !! The histogram is written with the two axes as first rows, then the data as a 2d-matrix
    !> @param[in] filename  name of the file to write to
    !> @param[in] label1  label of the first axis
    !> @param[in] label2  label of the second axis
    !> @param[in] hists  array of integers containing the histogram data for each pair of bins
    !> @param[in] bins1  bins of the first axis
    !> @param[in] bins2  bins of the second axis
    subroutine print_2d_hist(filename, label1, label2, hist, bins1, bins2)
      implicit none
      character(len=*), intent(in) :: filename, label1, label2
      real(dp), intent(in) :: bins1(:), bins2(:)
      integer, intent(in) :: hist(:,:)
      integer :: hist_unit
      integer :: i, j
      character, parameter :: tab = char(9)

         if(iProcIndex == root) then

            ! output the histogram
            hist_unit = get_free_unit()
            open(hist_unit, file = filename, status = 'unknown')
            write(hist_unit,"(A, A)") "# Boundaries of the bins of the first (vertical) dimension - ", label1

            do j = 1, size(bins1)
               write(hist_unit, '(G17.5)', advance = 'no') bins1(j)
            end do
            write(hist_unit, '()', advance = 'yes')

            write(hist_unit, "(A, A)") "# Boundaries of the bins of the second (horizontal) dimension - ", label2

            do j = 1, size(bins2)
               write(hist_unit, '(G17.5)', advance = 'no') bins2(j)
            end do
            write(hist_unit, '()', advance = 'yes')

            write(hist_unit, "(A)") "# Histogram - Note: Values laying exactly on a bounday are binned to the left"
            do i = 1, size(bins1)-1
                do j = 1, size(bins2)-1
                  write(hist_unit, '(G17.5)', advance = 'no') hist(i,j)
               end do
               write(hist_unit, '()', advance = 'yes')
            end do

            close(hist_unit)
         end if

    end subroutine print_2d_hist

    !> Wrapper function to create a 2d-histogram of the shift scale factors over energy
    !> @param[in] EnergyBinsNum  resolution of the energy axis (number of bins)
    !> @param[in] FValBinsNum  resolution of the factor axis (number of bins)
    subroutine print_fval_energy_hist(EnergyBinsNum, FValBinsNum)
      use CalcData, only: tAutoAdaptiveShift
      implicit none
      integer, intent(in) :: EnergyBinsNum, FValBinsNum ! number of points in the histogram's axes
      integer, allocatable :: hist(:,:), allHist(:,:)
      real(dp), allocatable :: EnergyBins(:), FValBins(:)

      if(tAutoAdaptiveShift) then
         ! allocate the buffers
         allocate(EnergyBins(EnergyBinsNum+1))
         allocate(FValBins(FValBinsNum+1))
         allocate(hist(EnergyBinsNum,FValBinsNum))
         allocate(allHist(EnergyBinsNum, FValBinsNum))
         ! generate the histogram

         call generate_fval_energy_hist(hist, EnergyBins, FvalBins, EnergyBinsNum, &
              FvalBinsNum ,allHist)

         call print_2d_hist("FValsEnergyHist", "Energy", "FVal", allHist, EnergyBins, FValBins)

         ! deallocate the buffers
         if(allocated(allHist)) deallocate(allHist)
         if(allocated(hist)) deallocate(hist)
         if(allocated(EnergyBins)) deallocate(EnergyBins)
         if(allocated(FValBins)) deallocate(FValBins)
      end if
    end subroutine print_fval_energy_hist

    !> Wrapper function to create a 2d-histogram of the shift scale factors over population
    !> @param[in] PopBinsNum  resolution of the population axis (number of bins)
    !> @param[in] FValBinsNum  resolution of the factor axis (number of bins)
    subroutine print_fval_pop_hist(PopBinsNum, FValBinsNum)
      use CalcData, only: tAutoAdaptiveShift
      implicit none
      integer, intent(in) :: PopBinsNum, FValBinsNum ! number of points in the histogram's axes
      integer, allocatable :: hist(:,:), allHist(:,:)
      real(dp), allocatable :: PopBins(:), FValBins(:)

      if(tAutoAdaptiveShift) then
         ! allocate the buffers
         allocate(PopBins(PopBinsNum+1))
         allocate(FValBins(FValBinsNum+1))
         allocate(hist(PopBinsNum,FValBinsNum))
         allocate(allHist(PopBinsNum, FValBinsNum))
         ! generate the histogram

         call generate_fval_pop_hist(hist, PopBins, FvalBins, PopBinsNum, &
              FvalBinsNum ,allHist)

         call print_2d_hist("FValsPopHist", "Population", "FVal", allHist, PopBins, FValBins)

         ! deallocate the buffers
         if(allocated(allHist)) deallocate(allHist)
         if(allocated(hist)) deallocate(hist)
         if(allocated(PopBins)) deallocate(PopBins)
         if(allocated(FValBins)) deallocate(FValBins)
      end if
    end subroutine print_fval_pop_hist
!------------------------------------------------------------------------------------------!

    !> For a given value, get the position in a histogram of given window sizes
    !> @param[in] val  value to get the position
    !> @param[in] minVal  smallest value appearing in the histogram
    !> @param[in] nPoints  number of bins in the histogram
    !> @param[in] windowSize  size of each bin
    !> @return ind  index of val in the histogram
    function getHistIndex(val, minVal, nPoints, windowSize) result(ind)
        implicit none
        real(dp), intent(in) :: val, minVal, windowSize
        integer, intent(in) :: nPoints
        integer :: ind

        if(val <= minVal) then
            ind = 1
        else if(val > minVal + nPoints*windowSize) then
            ind = nPoints
        else
            ind = ceiling((val - minVal) / windowSize)
        end if
    end function getHistIndex

    !> Create the data written out in the histogram of shift factor over energy.
    !! The generated data can be passed to print_2d_hist. This is a synchronizing routine.
    !> @param[out] hist  on return, histogram data of this proc only
    !> @param[out] histEnergy  on return, energy axis of the histogram
    !> @param[out] histAccRate  on return, shift factor axis of the histogram
    !> @param[in] enPoints  number of bins on the energy axis
    !> @param[in] accRatePoints  number of bins on the shift factor axis
    !> @param[out] allHist  on return, histogram data over all procs
    subroutine generate_fval_energy_hist(hist, histEnergy, histAccRate, &
        enPoints, accRatePoints, allHist)
        use global_det_data, only: det_diagH, get_acc_spawns, get_tot_spawns
        ! count the acceptance ratio per energy and create
        ! a histogram hist(:,:) with axes histEnergy(:), histAccRate(:)
        ! entries of hist(:,:) : numer of occurances
        ! entries of histEnergy(:) : energies of histogram entries (with tolerance, first dimension)
        ! entries of histAccRate(:) : acc. rates of histogram entries (second dimension)
        implicit none
        ! number of energy/acc rate windows in the histogram
        integer, intent(in) :: enPoints, accRatePoints
        integer, intent(out) :: hist(enPoints,accRatePoints)
        integer, intent(out) :: allHist(:,:)
        real(dp), intent(out) :: histEnergy(enPoints), histAccRate(accRatePoints)
        integer :: i, run
        integer :: enInd, arInd
        real(dp) :: minEn, maxEn, enWindow, arWindow, locMinEn, locMaxEn, totSpawn

        ! get the energy window size (the acc. rate window size is just 1.0/(accRatePoints+1))
        locMinEn = det_diagH(1)
        locMaxEn = det_diagH(1)
        do i = 2, int(TotWalkers)
            if(det_diagH(i) > locMaxEn) locMaxEn = det_diagH(i)
            if(det_diagH(i) < locMinEn) locMinEn = det_diagH(i)
        end do
        ! communicate the energy window
        call MPIAllReduce(locMinEn, MPI_MIN, minEn)
        call MPIAllReduce(locMaxEn, MPI_MAX, maxEn)
        enWindow = (maxEn - minEn) / real(enPoints,dp)

        if(enWindow > eps) then
            ! set up the histogram axes
            arWindow = 1.0_dp / real(accRatePoints,dp)
            do i = 1, accRatePoints+1
                histAccRate(i) = (i-1)*arWindow
            end do

            do i = 1, enPoints+1
                histEnergy(i) = minEn + (i-1)*enWindow
            end do

            ! then, fill the histogram itself
            hist = 0
            do i = 1, int(TotWalkers)
                do run = 1, inum_runs
                    totSpawn = get_tot_spawns(i,run)
                    if(abs(totSpawn) > eps) then
                        enInd = getHistIndex(det_diagH(i), minEn, enPoints, enWindow)
                        arInd = getHistIndex(get_acc_spawns(i,run)/totSpawn, &
                            0.0_dp, accRatePoints, arWindow)
                        hist(enInd, arInd) = hist(enInd, arInd) + 1
                    end if
                end do
            end do

            ! communicate the histogram
            call MPISum(hist, allHist)
        else
            write(stdout,*) "WARNING: Empty energy histogram of acceptance rates"
        end if
    end subroutine generate_fval_energy_hist

    !> Create the data written out in the histogram of shift factor over population.
    !! The generated data can be passed to print_2d_hist. This is a synchronizing routine.
    !> @param[out] hist  on return, histogram data of this proc only
    !> @param[out] histPop  on return, population axis of the histogram
    !> @param[out] histAccRate  on return, shift factor axis of the histogram
    !> @param[in] popPoints  number of bins on the population axis
    !> @param[in] accRatePoints  number of bins on the shift factor axis
    !> @param[out] allHist  on return, histogram data over all procs
    subroutine generate_fval_pop_hist(hist, histPop, histAccRate, &
        popPoints, accRatePoints, allHist)
        use global_det_data, only: get_acc_spawns, get_tot_spawns
        ! count the acceptance ratio per energy and create
        ! a histogram hist(:,:) with axes histEnergy(:), histAccRate(:)
        ! entries of hist(:,:) : numer of occurances
        ! entries of histPop(:) : populations of histogram entries (with tolerance, first dimension)
        ! entries of histAccRate(:) : acc. rates of histogram entries (second dimension)
        implicit none
        ! number of energy/acc rate windows in the histogram
        integer, intent(in) :: popPoints, accRatePoints
        integer, intent(out) :: hist(popPoints,accRatePoints)
        integer, intent(out) :: allHist(:,:)
        real(dp), intent(out) :: histPop(popPoints), histAccRate(accRatePoints)
        integer :: i, run
        integer :: popInd, arInd
        real(dp) :: maxPop, locMaxPop, totSpawn, pop
        real(dp), dimension(lenof_sign) :: sgn

        ! The bins of the population has width of 1, except the last bin which
        ! constains all populations larger than popPoints-1
        hist = 0
        locMaxPop = 0.0
        do i = 1, int(TotWalkers)
            call extract_sign(CurrentDets(:,i), sgn)
            do run = 1, inum_runs
                pop = mag_of_run(sgn,run)
                if(pop > locMaxPop) locMaxPop = pop
                totSpawn = get_tot_spawns(i,run)
                if(abs(totSpawn) > eps) then
                    popInd = getHistIndex(pop, 0.0_dp, popPoints, 1.0_dp)
                    arInd = getHistIndex(get_acc_spawns(i,run)/totSpawn, &
                        0.0_dp, accRatePoints, 1.0_dp/real(accRatePoints,dp))
                    hist(popInd, arInd) = hist(popInd, arInd) + 1
                end if
            end do
        end do

        call MPIAllReduce(locMaxPop, MPI_MAX, maxPop)
        do i = 1, popPoints+1
            histPop(i) = real(i-1)
        end do
        if(maxPop>histPop(popPoints+1))  histPop(popPoints+1)= maxPop

        do i = 1, accRatePoints+1
            histAccRate(i) = (i-1)/real(accRatePoints,dp)
        end do

        ! communicate the histogram
        call MPISum(hist, allHist)

    end subroutine generate_fval_pop_hist

!------------------------------------------------------------------------------------------!

    subroutine end_iteration_print_warn (totWalkersNew)

        ! Worker function for PerformFciMCycPar. Prints warnings about
        ! particle blooms and memory usage.
        integer, intent(in) :: totWalkersNew
        integer :: i
        real(dp) :: rat

        ! Too many particles?
        rat = real(TotWalkersNew,dp) / real(MaxWalkersPart,dp)
        if (rat > 0.95_dp) then
#ifdef DEBUG_
            if(tMolpro) then
                write(stderr, '(a)') '*WARNING* - Number of particles/determinants &
                                 &has increased to over 95% of allotted memory. &
                                 &Errors imminent. Increase MEMORYFACWALKERS, or reduce rate of growth.'
            else
                write(stderr, '(a)') '*WARNING* - Number of particles/determinants &
                                 &has increased to over 95% of allotted memory. &
                                 &Errors imminent. Increase MEMORYFACPART, or reduce rate of growth.'
            end if
#else
            if(tMolpro) then
                write(stderr,*) '*WARNING* - Number of particles/determinants &
                                 &has increased to over 95% of allotted memory on task ', iProcIndex, '. &
                                 &Errors imminent. Increase MEMORYFACWALKERS, or reduce rate of growth.'
            else
                write(stderr,*) '*WARNING* - Number of particles/determinants &
                                 &has increased to over 95% of allotted memory on task ', iProcIndex, '. &
                                 &Errors imminent. Increase MEMORYFACPART, or reduce rate of growth.'
            end if
#endif
            call neci_flush(stderr)
        end if

        ! Are ony of the sublists near the end of their alloted space?
        if (nNodes > 1) then
            do i = 0, nNodes-1
                rat = real(ValidSpawnedList(i) - InitialSpawnedSlots(i),dp) /&
                             real(InitialSpawnedSlots(1), dp)
                if (rat > 0.95_dp) then
#ifdef DEBUG_
                    if(tMolpro) then
                        write(stderr, '(a)') '*WARNING* - Highest processor spawned &
                                         &particles has reached over 95% of allotted memory.&
                                         &Errors imminent. Increase MEMORYFACSPAWNED, or reduce spawning rate.'
                    else
                        write(stderr, '(a)') '*WARNING* - Highest processor spawned &
                                         &particles has reached over 95% of allotted memory.&
                                         &Errors imminent. Increase MEMORYFACSPAWN, or reduce spawning rate.'
                    end if
#else
                    if(tMolpro) then
                        write(stderr,*) '*WARNING* - Highest processor spawned &
                                         &particles has reached over 95% of allotted memory on task ',iProcIndex,' .&
                                         &Errors imminent. Increase MEMORYFACSPAWNED, or reduce spawning rate.'
                    else
                        write(stderr,*) '*WARNING* - Highest processor spawned &
                                         &particles has reached over 95% of allotted memory on task ',iProcIndex,' .&
                                         &Errors imminent. Increase MEMORYFACSPAWN, or reduce spawning rate.'
                    end if
#endif
                    call neci_flush(stderr)
                end if
            end do
        else
            rat = real(ValidSpawnedList(0), dp) / real(MaxSpawned, dp)
            if (rat > 0.95_dp) then
#ifdef DEBUG_
                if(tMolpro) then
                    write(stderr, '(a)') '*WARNING* - Highest processor spawned &
                                     &particles has reached over 95% of allotted memory.&
                                     &Errors imminent. Increase MEMORYFACSPAWNED, or reduce spawning rate.'
                else
                    write(stderr, '(a)') '*WARNING* - Highest processor spawned &
                                     &particles has reached over 95% of allotted memory.&
                                     &Errors imminent. Increase MEMORYFACSPAWN, or reduce spawning rate.'
                end if
#else
                if(tMolpro) then
                    write(stderr,*) '*WARNING* - Highest processor spawned &
                                     &particles has reached over 95% of allotted memory on task ',iProcIndex,' .&
                                     &Errors imminent. Increase MEMORYFACSPAWNED, or reduce spawning rate.'
                else
                    write(stderr,*) '*WARNING* - Highest processor spawned &
                                     &particles has reached over 95% of allotted memory on task ',iProcIndex,' .&
                                     &Errors imminent. Increase MEMORYFACSPAWN, or reduce spawning rate.'
                end if
#endif
                call neci_flush(stderr)
            end if
         end if

    end subroutine end_iteration_print_warn


    subroutine getProjEOffset()
        ! get the offset of the projected energy versus the total energy,
        ! which is the reference energy

        implicit none
        ! if the reference energy is used as an offset to the hamiltonian (default behaviour)
        ! just get it
        if(.not.tZeroRef) then
            OutputHii = Hii
        ! else, calculate the reference energy
        else if (tHPHF) then
            OutputHii = hphf_diag_helement (ProjEDet(:,1), iLutRef(:,1))
        else if (tGUGA) then
            OutputHii = calcDiagMatEleGUGA_nI(ProjEDet(:,1))
        else
            OutputHii = get_helement (ProjEDet(:,1), ProjEDet(:,1), 0)
        end if


    end subroutine getProjEOffset

end module
