#include "macros.h"
module fcimc_output

    use SystemData, only: nel, tHPHF, tFixLz, tMolpro, tMolproMimic, MolproID, &
                          tGen_4ind_weighted, tGen_4ind_2, tGUGA, tGen_sym_guga_mol, &
                          tGen_nosym_guga, t_consider_diff_bias, tgen_guga_crude, &
                          t_new_real_space_hubbard, t_no_ref_shift

    use LoggingData, only: tLogComplexPops, tMCOutput, tCalcInstantS2, &
                           tCalcInstantS2Init, instant_s2_multiplier_init, &
                           instant_s2_multiplier, tPrintFCIMCPsi, &
                           iWriteHistEvery, tDiagAllSpaceEver, OffDiagMax, &
                           OffDiagBinRange, tCalcVariationalEnergy, &
                           iHighPopWrite, tLogEXLEVELStats, StepsPrint, &
                           maxInitExLvlWrite, AllInitsPerExLvl

    use hist_data, only: Histogram, AllHistogram, InstHist, AllInstHist, &
                         BeforeNormHist, iNoBins, BinRange, HistogramEnergy, &
                         AllHistogramEnergy

    use CalcData, only: tTruncInitiator, tTrialWavefunction, tReadPops, &
                        DiagSft, tSpatialOnlyHash, tOrthogonaliseReplicas, &
                        StepsSft, tPrintReplicaOverlaps, tStartTrialLater, &
                        frq_step_size, tEN2, &
                        frequency_bins_singles, frequency_bins_para, &
                        frequency_bins_anti, frequency_bins_doubles, &
                        frequency_bins_type2, frequency_bins_type3, &
                        frequency_bins_type4, frequency_bins_type2_diff, &
                        frequency_bins_type3_diff, n_frequency_bins, &
                        tSemiStochastic, t_truncate_spawns, tLogGreensfunction

    use DetBitOps, only: FindBitExcitLevel, count_open_orbs, EncodeBitDet, &
                         TestClosedShellDet

    use IntegralsData, only: frozen_orb_list, frozen_orb_reverse_map, &
                             nel_pre_freezing

    use DetCalcData, only: det, fcidets, ReIndex, NDet, NRow, HAMIL, LAB

    use bit_reps, only: decode_bit_det, test_flag, extract_sign, get_initiator_flag

    use semi_stoch_procs, only: return_most_populated_states, GLOBAL_RUN

    use bit_rep_data, only: niftot, nifd, flag_initiator

    use hist, only: calc_s_squared_star, calc_s_squared

    use fcimc_helper, only: LanczosFindGroundE

    use hphf_integrals, only: hphf_diag_helement
    use Determinants, only: write_det, get_helement
    use adi_data, only: AllCoherentDoubles, AllIncoherentDets, nRefs, &
         ilutRefAdi, tAdiActive, nConnection, AllConnection

    use rdm_data, only: en_pert_main

    use Parallel_neci

    use FciMCData

    use constants

    use sort_mod

    use util_mod

    use tau_search, only: comm_frequency_histogram, comm_frequency_histogram_spec

    use real_time_data, only:  gf_count, &
                              normsize, snapShotOrbs, &
                              current_overlap, t_real_time_fciqmc, elapsedRealTime, &
                              elapsedImagTime, overlap_real, overlap_imag, dyn_norm_psi,&
                              dyn_norm_red, real_time_info, allPopSnapshot, numSnapshotOrbs

    use tc_three_body_data, only: tLMatCalc, lMatCalcStatsIters, &
                                  lMatCalcHit,   lMatCalcTot,    lMatCalcHUsed,   lMatCalcHSize
    use fortran_strings, only: str

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

        IF(iProcIndex.eq.root) THEN
!Print out initial starting configurations
            WRITE(iout,*) ""
            IF(tTruncInitiator) THEN
               WRITE(initiatorstats_unit,"(A2,A17,16A23)", advance = 'no') &
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
            ENDIF
            IF(tLogComplexPops) THEN
                WRITE(complexstats_unit,"(A)") '#   1.Step  2.Shift     3.RealShift     4.ImShift   5.TotParts      " &
                & //"6.RealTotParts      7.ImTotParts'
            ENDIF
            if (tLogEXLEVELStats) then
                write (EXLEVELStats_unit, '(a)', advance='no') '# 1.Step'
                k = 1
                do run = 1, inum_runs
                    tchar_r = ''
                    if (inum_runs>1) then
                        write(tchar_r,*)run
                        tchar_r = '(run='//trim(adjustl(tchar_r))//')'
                    endif
                    do i = 0, 2
                        write(tchar_i,*)i
                        do j = 0, NEl
                            k = k + 1
                            write(tchar_j,*)j
                            write(tchar_k,*)k
                            write (EXLEVELStats_unit, '(1x,a)', &
                                   advance='no') trim(adjustl(tchar_k))// &
                                   &'.W'//trim(adjustl(tchar_j))//'^'// &
                                   trim(adjustl(tchar_i))// &
                                   trim(adjustl(tchar_r))
                        enddo ! j
                    enddo ! i
                enddo ! run
                write (EXLEVELStats_unit, '()', advance='yes')
            endif ! tLogEXLEVELStats

#ifdef CMPLX_
            if(tMCOutput) then
                write(iout, '(a)') "       Step     Shift      WalkerCng(Re)  &
                       &WalkerCng(Im)    TotWalkers(Re)   TotWalkers(Im)    &
                       &Proj.E(Re)   ProjE(Im)     Proj.E.ThisCyc(Re)  &
                       &Proj.E.ThisCyc(Im)   NoatHF(Re)   NoatHF(Im)   &
                       &NoatDoubs      AccRat     UniqueDets   NumDetsSpawned   IterTime"
            endif
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
                write(iout, "(A)", advance = 'no') "        Step    Shift           &
                      &WalkerCng       GrowRate        TotWalkers      Annihil         &
                      &NoDied          NoBorn          Proj.E          Av.Shift        &
                      &Proj.E.Cyc"
                if (tTrialWavefunction .or. tStartTrialLater) write(iout, "(A)", advance = 'no') &
                      "    Trial.E.Cyc "
                write(iout, "(A)", advance = 'yes') "      NoatHF          NoatDoubs       &
                &AccRat        UniqueDets    NumDetsSpawned   IterTime"
            endif
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

        ENDIF

    END SUBROUTINE WriteFciMCStatsHeader

    subroutine WriteFCIMCStats()
        use guga_matrixElements, only: calcDiagMatEleGUGA_nI

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
        endif

        !To prevent /0 problems
        do run=1,inum_runs
            if(.not. near_zero(AllNoBorn(run))) then
                FracFromSing(run)=real(AllSpawnFromSing(run),dp) / real(AllNoBorn(run),dp)
            else
                FracFromSing(run)=0.0_dp
            endif

            if(t_no_ref_shift)then
                if(.not. tguga)then
                    E_ref_tmp(run) = get_helement (ProjEDet(:,run), ProjEDet(:,run), 0)
                else
                    E_ref_tmp(run) = calcDiagMatEleGUGA_nI(ProjEDet(:,run))
                end if
            else
                E_ref_tmp(run) = 0.0_dp
            end if


        enddo

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
                write (iout, "(I12,13G16.7,2I12,G13.5)") &
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
            endif
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
            endif
            if (tLogComplexPops) then
                write (complexstats_unit,"(I12,6G16.7)") &
                    Iter + PreviousCycles, DiagSft, DiagSftRe, DiagSftIm, &
                    sum(AllTotParts), AllTotParts(1), AllTotParts(lenof_sign)
            endif
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
                endif

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
                endif
                write(fcimcstats_unit, "()", advance = 'yes')

            if(tMCOutput) then
                write (iout, "(I12,10G16.7)", advance = 'no') &
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
                     write(iout, "(G20.11)", advance = 'no') &
                         (tot_trial_numerator(1)/tot_trial_denom(1))
                else if (tStartTrialLater) then
                     write(iout, "(G20.11)", advance = 'no') 0.0_dp
                end if
                write (iout, "(3G16.7,2I12,G13.5)", advance = 'yes') &
                    AllNoatHF(1), &
                    AllNoatDoubs(1), &
                    AccRat(1), &
                    AllTotWalkers, &
                    nspawned_tot, &
                    IterTime
            endif
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
            endif
#endif
            if (tLogEXLEVELStats) then
                write (EXLEVELStats_unit, '(i12)', advance='no') &
                      &Iter + PreviousCycles
                do run = 1, inum_runs
                    do i = 0, 2
                        do j = 0, NEl
                            write (EXLEVELStats_unit, '(1x,G18.9e3)', &
                                   advance='no') AllEXLEVEL_WNorm(i,j,run)
                        enddo ! j
                    enddo ! i
                enddo ! run
                write (EXLEVELStats_unit, '()', advance='yes')
            endif ! tLogEXLEVELStats

            if (tMCOutput .and. tLMatCalc .and. mod(Iter, lMatCalcStatsIters) == 0) then
                write(iout, *) "============ LMatCalc Caching Stats ==============="
                write(iout, *) "LMatCalc Cache Fill Ratio: ", &
                    real(lMatCalcHUsed,dp)/real(lMatCalcHSize,dp)
                write(iout, *) "LMatCalc Cache Hit Rate  : ", lMatCalcHit/real(lMatCalcTot)
                lMatCalcHit = 0
                lMatCalcTot = 0
                write(iout, *) "==================================================="
            end if

            if(tMCOutput) then
                call neci_flush(iout)
            endif
            call neci_flush(fcimcstats_unit)
            if (inum_runs.eq.2) call neci_flush(fcimcstats_unit2)
            if (tLogEXLEVELStats) call neci_flush(EXLEVELStats_unit)

        endif

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
        endif

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
            endif
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
                endif
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
            endif

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
            endif
            if(t_real_time_fciqmc) then
                call stats_out(state, .true., elapsedRealTime, 'Re. time')
                call stats_out(state, .true., elapsedImagTime, 'Im. time')
            else
                call stats_out(state,.false., TotImagTime, 'Im. time')
            endif

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
                enddo
                do iGf = 1, gf_count
                   write(tmgf, '(i5)') iGf
                   do p = 1, normsize
                      write(tmpc, '(i5)') p
                      call stats_out(state,.false.,real(current_overlap(p,iGf)), 'Re. <y(0)|y(t)>(rep ' // &
                           trim(adjustl(tmpc)) // ' i=' // trim(adjustl(tmgf)) //  ')')
                      call stats_out(state,.false.,aimag(current_overlap(p,iGf)), 'Im. <y(0)|y(t)>(rep ' // &
                           trim(adjustl(tmpc)) // ' i=' // trim(adjustl(tmgf)) //')')
                   enddo
                enddo
                if(t_real_time_fciqmc) then
                    do p = 1, numSnapshotOrbs
                        ! if any orbitals are monitored, output their population
                        write(tmpc, '(i5)') snapShotOrbs(p)
                        call stats_out(state,.false.,allPopSnapshot(p),'Population of ' &
                            // trim(adjustl(tmpc)))
                    end do
                endif
            endif

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
            if (tMCOutput) write(iout, *)

            call neci_flush(state%funit)
            if (tTruncInitiator) call neci_flush(state_i%funit)
            call neci_flush(iout)

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
         if (tMCOutput) write(iout, '("#")', advance='no')
         state%prepend = state%init
      else if (.not. state%prepend) then
         write(state%funit, '(" ")', advance='no')
         if (tMCOutput) write(iout, '(" ")', advance='no')
      end if
    end subroutine write_padding_init

    subroutine writeMsWalkerCountsAndCloseUnit()
        integer :: ms, tempni(1:nel)
        integer(int64) :: i
        real(dp) :: totWalkPopByMsReal(nel+1), totWalkPopByMsImag(nel+1), &
                    tempSign(lenof_sign)

        do i=1,TotWalkers
            call extract_sign(WalkVecDets(:,i),TempSign)
            call decode_bit_det(TempnI, WalkVecDets(:,i))
            ms = sum(get_spin_pn(Tempni(1:nel)))
            walkPopByMsReal(1+nel/2+ms/2) = walkPopByMsReal(1+nel/2+ms/2)+abs(TempSign(1))
#ifdef CMPLX_
            walkPopByMsImag(1+nel/2+ms/2) = walkPopByMsImag(1+nel/2+ms/2)+abs(TempSign(2))
#endif
            write(mswalkercounts_unit,*) ms, TempSign
        enddo

        totWalkPopByMsReal = walkPopByMsReal
        totWalkPopByMsImag = walkPopByMsImag

        ! sum the populations from all processors
        call MPISumAll(walkPopByMsReal, totWalkPopByMsReal)
        call MPISumAll(walkPopByMsImag, totWalkPopByMsImag)
        ms = -1*nel
        do i =1,nel+1
            write(mswalkercounts_unit,*) ms, totWalkPopByMsReal(i),&
                totWalkPopByMsImag(i), (totWalkPopByMsReal(i)**2+totWalkPopByMsImag(i)**2)**(0.5)
            ms = ms+2
        enddo
        close(mswalkercounts_unit)

    endsubroutine writeMsWalkerCountsAndCloseUnit





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
            enddo
        enddo
#ifdef CMPLX_
        norm2=SQRT(sum(norm1))
#else
        norm1=SQRT(norm1)
#endif
        WRITE(iout,*) "Total FCIMC Wavefuction normalisation:",norm1
        do i=1,Det
            do j=1,lenof_sign
#ifdef CMPLX_
                AllHistogram(j,i)=AllHistogram(j,i)/norm2
#else
                AllHistogram(j,i)=AllHistogram(j,i)/norm1(j)
#endif
            enddo
        enddo

        iunit = 0
        IF(tPrintFCIMCPsi) THEN
!Order and print wavefunction


            IF(iProcIndex.eq.0) THEN

                ! We now want to order AllHistogram, taking the corresponding
                ! element(s) of FCIDets with it...
                call sort (AllHistogram, FCIDets)

                OPEN(iunit,FILE='FCIMCPsi',STATUS='UNKNOWN')

                norm=0.0_dp
                do i=1,Det
                    do j=1,lenof_sign
                        norm(j)=norm(j)+AllHistogram(j,i)**2
                    enddo
!write out FCIMC Component weight (normalised), current normalisation, excitation level
                    ExcitLevel = FindBitExcitLevel(iLutHF, FCIDets(:,i), nel)
                    CALL decode_bit_det(nI,FCIDets(0:NIfTot,i))
#ifdef CMPLX_
                    WRITE(iunit,"(I13,G25.16,I6,G20.10)",advance='no') i,AllHistogram(1,i),ExcitLevel,sum(norm)
#else
                    WRITE(iunit,"(I13,G25.16,I6,G20.10)",advance='no') i,AllHistogram(1,i),ExcitLevel,norm(1)
#endif
                    do j=1,NEl-1
                        WRITE(iunit,"(I5)",advance='no') nI(j)
                    enddo
                    WRITE(iunit,"(I5)") nI(NEl)
                enddo

                CLOSE(iunit)

            ENDIF
        ENDIF

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
            WRITE(iout,*) "Writing out the average wavevector up to iteration number: ", Iter
            CALL neci_flush(iout)
        ENDIF

        IF(iProcIndex.eq.0) THEN
            AllHistogram(:,:)=0.0_dp
            AllInstHist(:,:)=0.0_dp
            AllAvAnnihil(:,:)=0.0_dp
            AllInstAnnihil(:,:)=0.0_dp
        ENDIF

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
                enddo
            enddo
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
                    ENDIF
                    IF(.not. near_zero(norm3_c)) THEN
                        AllAvAnnihil(j,i)=AllAvAnnihil(j,i)/norm3_c
                    ENDIF
#else
                    AllHistogram(j,i)=AllHistogram(j,i)/norm(j)
                    AllInstHist(j,i)=AllInstHist(j,i)/norm1(j)
                    IF(.not. near_zero(norm2(j))) THEN
                        AllInstAnnihil(j,i)=AllInstAnnihil(1,i)/norm2(j)
                    ENDIF
                    IF(.not. near_zero(norm3(j))) THEN
                    AllAvAnnihil(j,i)=AllAvAnnihil(j,i)/norm3(j)
                    ENDIF
#endif
                enddo
            enddo

            io1 = get_free_unit()
            OPEN(io1,FILE=abstr,STATUS='UNKNOWN')

            abstr = 'Energies-'//str(Iter - iWriteHistEvery)
            abstr2 = 'Energies-'//str(Iter)

            io2 = get_free_unit()
            OPEN(io2,FILE=abstr2,STATUS='UNKNOWN')

            INQUIRE(FILE=abstr,EXIST=exists)
            IF(exists) THEN
                io3 = get_free_unit()
                OPEN(io3,FILE=abstr,STATUS='OLD',POSITION='REWIND',ACTION='READ')
                do while(.true.)
                    READ(io3,"(I13,3G25.16)",END=99) IterRead,ShiftRead,AllERead,NumParts
                    WRITE(io2,"(I13,3G25.16)") IterRead,ShiftRead,AllERead,NumParts
                enddo
99              CONTINUE
#ifdef CMPLX_
                IF(near_zero(AllHFOut(1))) then
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft,AllERead,SUM(AllTotPartsLastOutput)
                ELSE
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft,AllENumOut/AllHFOut,SUM(AllTotPartsLastOutput)
                ENDIF
#else
                IF(near_zero(AllHFOut(1))) THEN
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft(1),AllERead,AllTotPartsLastOutput(1)
                ELSE
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft(1),AllENumOut(1)/AllHFOut(1),AllTotPartsLastOutput(1)
                ENDIF
#endif
                CLOSE(io2)
                CLOSE(io3)

            ELSE
                OPEN(io2,FILE=abstr2,STATUS='UNKNOWN')
#ifdef CMPLX_
                WRITE(io2,"(I13,3G25.16)") Iter,DiagSft,AllENumOut/AllHFOut,SUM(AllTotPartsLastOutput)
#else
                WRITE(io2,"(I13,3G25.16)") Iter,DiagSft(1),AllENumOut(1)/AllHFOut(1),AllTotPartsLastOutput(1)
#endif
                CLOSE(io2)
            ENDIF


            norm=0.0_dp
            norm1=0.0_dp
            Tot_No_Unique_Dets = 0
            do i=1,Det

                posn=binary_search(CurrentDets(:,1:TotWalkers), FCIDETS(:,i), NifD+1)

                if (posn.lt.0) then
                    FinalPop = 0
                else
                    FinalPop = int(CurrentDets(NifD+1,posn))
                endif

                do j=1,lenof_sign
                    norm(j)=norm(j)+AllHistogram(j,i)**2
                    norm1(j)=norm1(j)+AllAvAnnihil(j,i)**2
                enddo
                IF(lenof_sign.eq.1) THEN
                    WRITE(io1,"(I13,6G25.16,I13,G25.16)") i, &
                          AllHistogram(1,i), norm, AllInstHist(1,i), &
                          AllInstAnnihil(1,i), AllAvAnnihil(1,i), norm1, &
                          FinalPop, BeforeNormHist(i)
                ELSE
#ifdef CMPLX_
                    WRITE(io1,"(I13,6G25.16)") i, AllHistogram(1,i), sum(norm), &
                          AllInstHist(1,i), AllInstAnnihil(1,i), &
                          AllAvAnnihil(1,i), sum(norm1)
#else
                    WRITE(io1,"(I13,6G25.16)") i, AllHistogram(1,i), norm(1), &
                          AllInstHist(1,i), AllInstAnnihil(1,i), &
                          AllAvAnnihil(1,i), norm1(1)
#endif
                ENDIF
                IF(.not. near_zero(AllHistogram(1,i))) Tot_No_Unique_Dets = Tot_No_Unique_Dets + 1
            enddo
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

                call my_hpsi(Det,1,NROW,LAB,HAMIL,HOrderedHist,CKN,.true.)
                AvVarEnergy = DDOT(Det,HOrderedHist,1,CKN,1)

                CKN = 0.0_dp
                call my_hpsi(Det,1,NROW,LAB,HAMIL,HOrderedInstHist,CKN,.true.)
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
                        endif
                    enddo
                    if((val-1).ne.Tot_No_Unique_Dets) call stop_all(t_r,'Wrong counting')

                    call LanczosFindGroundE(ExpandedWalkerDets,Tot_No_Unique_Dets,GroundE_Ever,ProjGroundE,.true.)
                    if(abs(GroundE_Ever-ProjGroundE).gt.1.0e-7_dp) call stop_all(t_r,'Why do these not agree?!')
                    write(Tot_Unique_Dets_Unit,"(2I14,3G25.15)") Iter,Tot_No_Unique_Dets,AvVarEnergy,VarEnergy,GroundE_Ever
                else
                    write(Tot_Unique_Dets_Unit,"(2I14,2G25.15)") Iter,Tot_No_Unique_Dets,AvVarEnergy,VarEnergy
                endif
            else
                if(tDiagAllSpaceEver) then

                    allocate(ExpandedWalkerDets(NEl,Tot_No_Unique_Dets),stat=ierr)
                    val=1
                    do i=1,Det
                        if(.not. near_zero(AllHistogram(1,i))) then
                            call decode_bit_det(ExpandedWalkerDets(:,val),FCIDets(0:NIfTot,i))
                            val=val+1
                        endif
                    enddo
                    if((val-1).ne.Tot_No_Unique_Dets) call stop_all(t_r,'Wrong counting')

                    call LanczosFindGroundE(ExpandedWalkerDets,Tot_No_Unique_Dets,GroundE_Ever,ProjGroundE,.true.)
                    if(abs(GroundE_Ever-ProjGroundE).gt.1.0e-7_dp) call stop_all(t_r,'Why do these not agree?!')
                    write(Tot_Unique_Dets_Unit,"(2I14,3G25.15)") Iter,Tot_No_Unique_Dets,GroundE_Ever
                else
                    write(Tot_Unique_Dets_Unit,"(2I14)") Iter, Tot_No_Unique_Dets
                endif
            endif

            CLOSE(io1)
        ENDIF
        InstHist(:,:)=0.0_dp
        InstAnnihil(:,:)=0.0_dp

    END SUBROUTINE WriteHistogram

!This routine will write out the average hamiltonian from the spawning run up until now.
    SUBROUTINE WriteHamilHistogram()
        INTEGER :: i,j
        integer :: iunit
        CHARACTER(len=22) :: abstr

!This will open a file called HamilHist-"Iter" on unit number 17.
        abstr = 'HamilHist-'//str(Iter)
        IF(iProcIndex.eq.0) THEN
            WRITE(iout,*) "Writing out the average hamiltonian up to iteration number: ", Iter
            CALL neci_flush(iout)
        ENDIF

        IF(iProcIndex.eq.0) THEN
            AllHistHamil(:,:)=0.0_dp
            AllAvHistHamil(:,:)=0.0_dp
        ENDIF

        CALL MPIReduce(HistHamil,MPI_SUM,AllHistHamil)
        CALL MPIReduce(AvHistHamil,MPI_SUM,AllAvHistHamil)

        IF(iProcIndex.eq.0) THEN
!How do we normalise this!
            iunit = get_free_unit()
            OPEN(iunit,FILE=abstr,STATUS='UNKNOWN')
            do i=1,Det
                do j=1,Det
                    WRITE(iunit,*) j,i,AllAvHistHamil(j,i),AllHistHamil(j,i)
                enddo
                WRITE(iunit,*) ""
            enddo
            CLOSE(iunit)
        ENDIF
        HistHamil(:,:)=0.0_dp

    END SUBROUTINE WriteHamilHistogram

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
        ENDIF
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
            enddo
            Norm=0.0_dp
            do i=1,iOffDiagNoBins
                Norm=Norm+AllSinglesHist(i)
            enddo
!            WRITE(iout,*) "AllSinglesHistNorm = ",Norm
            do i=1,iOffDiagNoBins
                AllSinglesHist(i)=AllSinglesHist(i)/Norm
            enddo

!            Norm=0.0_dp
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistOccOcc(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistOccOcc(i)=AllSinglesHistOccOcc(i)/Norm
            enddo
!            Norm=0.0_dp
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistOccVirt(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistOccVirt(i)=AllSinglesHistOccVirt(i)/Norm
            enddo
!            Norm=0.0_dp
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistVirtOcc(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistVirtOcc(i)=AllSinglesHistVirtOcc(i)/Norm
            enddo
!            Norm=0.0_dp
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistVirtVirt(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistVirtVirt(i)=AllSinglesHistVirtVirt(i)/Norm
            enddo


            io(1) = get_free_unit()
            OPEN(io(1),FILE='EVERYENERGYHIST',STATUS='UNKNOWN')
            io(2) = get_free_unit()
            OPEN(io(2),FILE='ATTEMPTENERGYHIST',STATUS='UNKNOWN')
            io(3) = get_free_unit()
            OPEN(io(3),FILE='SPAWNENERGYHIST',STATUS='UNKNOWN')

            EnergyBin=BinRange/2.0_dp
            do i=1,iNoBins
                IF(AllHistogramEnergy(i).gt.0.0_dp) WRITE(io(1),*) EnergyBin, AllHistogramEnergy(i)
                IF(AllAttemptHist(i).gt.0.0_dp) WRITE(io(2),*) EnergyBin, AllAttemptHist(i)
                IF(AllSpawnHist(i).gt.0.0_dp) WRITE(io(3),*) EnergyBin, AllSpawnHist(i)
                EnergyBin=EnergyBin+BinRange
            enddo
            CLOSE(io(1))
            CLOSE(io(2))
            CLOSE(io(3))
            OPEN(io(1),FILE='SINGLESHIST',STATUS='UNKNOWN')
            OPEN(io(2),FILE='ATTEMPTSINGLESHIST',STATUS='UNKNOWN')
            OPEN(io(3),FILE='DOUBLESHIST',STATUS='UNKNOWN')
            io(4) = get_free_unit()
            OPEN(io(4),FILE='ATTEMPTDOUBLESHIST',STATUS='UNKNOWN')
            io(5) = get_free_unit()
            OPEN(io(5),FILE='SINGLESHISTOCCOCC',STATUS='UNKNOWN')
            io(6) = get_free_unit()
            OPEN(io(6),FILE='SINGLESHISTOCCVIRT',STATUS='UNKNOWN')
            io(7) = get_free_unit()
            OPEN(io(7),FILE='SINGLESHISTVIRTOCC',STATUS='UNKNOWN')
            io(8) = get_free_unit()
            OPEN(io(8),FILE='SINGLESHISTVIRTVIRT',STATUS='UNKNOWN')

            EnergyBin=-OffDiagMax+OffDiagBinRange/2.0_dp
            do i=1,iOffDiagNoBins
                IF(AllSinglesHist(i).gt.0.0_dp) WRITE(io(1),*) EnergyBin, AllSinglesHist(i)
                IF(AllSinglesAttemptHist(i).gt.0.0_dp) WRITE(io(2),*) EnergyBin, AllSinglesAttemptHist(i)
                IF(AllDoublesHist(i).gt.0.0_dp) WRITE(io(3),*) EnergyBin, AllDoublesHist(i)
                IF(AllDoublesAttemptHist(i).gt.0.0_dp) WRITE(io(4),*) EnergyBin, AllDoublesAttemptHist(i)
                IF(AllSinglesHistOccOcc(i).gt.0.0_dp) WRITE(io(5),*) EnergyBin, AllSinglesHistOccOcc(i)
                IF(AllSinglesHistOccVirt(i).gt.0.0_dp) WRITE(io(6),*) EnergyBin, AllSinglesHistOccVirt(i)
                IF(AllSinglesHistVirtOcc(i).gt.0.0_dp) WRITE(io(7),*) EnergyBin, AllSinglesHistVirtOcc(i)
                IF(AllSinglesHistVirtVirt(i).gt.0.0_dp) WRITE(io(8),*) EnergyBin, AllSinglesHistVirtVirt(i)
                EnergyBin=EnergyBin+OffDiagBinRange
!                WRITE(6,*) i
            enddo

            CLOSE(io(1))
            CLOSE(io(2))
            CLOSE(io(3))
            CLOSE(io(4))
            CLOSE(io(5))
            CLOSE(io(6))
            CLOSE(io(7))
            CLOSE(io(8))
        ENDIF

    END SUBROUTINE WriteHistogramEnergies

    !Routine to print the highest populated determinants at the end of a run
    SUBROUTINE PrintHighPops()
      use adi_references, only: update_ref_signs, print_reference_notification, nRefs
      use adi_data, only: tSetupSIs
        real(dp), dimension(lenof_sign) :: SignCurr, LowSign
        integer :: ierr,i,j,counter,ExcitLev,SmallestPos,HighPos,nopen
        integer :: full_orb, run
        real(dp) :: HighSign,reduce_in(1:2),reduce_out(1:2),Norm,AllNorm
        integer(n_int) , allocatable :: LargestWalkers(:,:)
        integer(n_int) , allocatable :: GlobalLargestWalkers(:,:)
        integer(n_int) :: HighestDet(0:NIfTot)
        integer, allocatable :: GlobalProc(:), tmp_ni(:)
        character(100) :: bufEnd, bufStart
        integer :: lenEnd, lenStart
        character(len=*), parameter :: t_r='PrintHighPops'

        character(1024) :: header
        character(25) :: format_string
        character(11), allocatable :: walker_string(:)
        character(13), allocatable :: amplitude_string(:)
        character(9), allocatable :: init_string(:)

        !Allocate memory to hold highest iHighPopWrite determinants
        allocate(LargestWalkers(0:NIfTot,iHighPopWrite),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"error allocating here")

        ! Return the most populated states in CurrentDets on *this* processor only.
        call return_most_populated_states(iHighPopWrite, GLOBAL_RUN, &
            LargestWalkers, norm = norm)

        call MpiSum(norm,allnorm)
        if(iProcIndex.eq.Root) norm=sqrt(allnorm)

        !Now have sorted list of the iHighPopWrite largest weighted determinants on the process
        if(iProcIndex.eq.Root) then
            allocate(GlobalLargestWalkers(0:NIfTot,iHighPopWrite),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"error allocating here")
            GlobalLargestWalkers(:,:)=0
            allocate(GlobalProc(iHighPopWrite),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"error allocating here")
            GlobalProc(:)=0
        endif

        do i=1,iHighPopWrite

            ! Find highest sign on each processor. Since all lists are
            ! sorted, this is just the first nonzero value.
           HighSign = 0.0_dp
           HighPos = 1
            do j=iHighPopWrite,1,-1
                call extract_sign (LargestWalkers(:,j), SignCurr)
                if (any(LargestWalkers(:,j) /= 0)) then

#ifdef CMPLX_
                    HighSign = sqrt(sum(abs(SignCurr(1::2)))**2 + sum(abs(SignCurr(2::2)))**2)
#else
                    HighSign = sum(real(abs(SignCurr),dp))
#endif

                    ! We have the largest sign
                    HighPos = j
                    exit
                end if
            end do
            reduce_in=(/ HighSign,real(iProcIndex,dp) /)
            call MPIAllReduceDatatype(reduce_in,1,MPI_MAXLOC,MPI_2DOUBLE_PRECISION,reduce_out)
            !Now, reduce_out(2) has the process of the largest weighted determinant - broadcast it!
            if(iProcIndex.eq.nint(reduce_out(2))) then
                HighestDet(0:NIfTot) = LargestWalkers(:,HighPos)
            else
                HighestDet(0:NIfTot) = 0
            endif
            call MPIBCast(HighestDet(:),NIfTot+1,nint(reduce_out(2)))
            if(iProcIndex.eq.Root) then
                GlobalLargestWalkers(0:NIfTot,i) = HighestDet(:)
                GlobalProc(i) = nint(reduce_out(2))
            endif

            !Now delete this highest determinant from the list on the corresponding processor
            if(iProcIndex.eq.nint(reduce_out(2))) then
                LargestWalkers(:,HighPos) = 0
                !No need to resort any more
!                call sort(LargestWalkers(:,1:iHighPopWrite), sign_lt, sign_gt)
            endif
        enddo

        ! This has to be done by all procs
        if(tAdiActive) call update_ref_signs()
        if(iProcIndex.eq.Root) then
            !Now print out the info contained in GlobalLargestWalkers and GlobalProc

            counter=0
            do i=1,iHighPopWrite
                !How many non-zero determinants do we actually have?
                call extract_sign(GlobalLargestWalkers(:,i),SignCurr)
#ifdef CMPLX_
                HighSign = sqrt(sum(abs(SignCurr(1::2)))**2 + sum(abs(SignCurr(2::2)))**2)
#else
                HighSign=sum(real(abs(SignCurr),dp))
#endif
                if (HighSign > 1.0e-7_dp) counter = counter + 1
            enddo


            write(iout,*) ""
            if (tReplicaReferencesDiffer) then
                write(iout,'(A)') "Current references: "
                do run = 1, inum_runs
                    call write_det(iout, ProjEDet(:,run), .true.)
                    call writeDetBit(iout, ilutRef(:, run), .true.)
                end do
            else
                write(iout,'(A)') "Current reference: "
                call write_det (iout, ProjEDet(:,1), .true.)
                if(tSetupSIs) call print_reference_notification(&
                     1,nRefs,"Used Superinitiator",.true.)
                write(iout,*) "Number of superinitiators", nRefs
            end if

            write(iout,*)
            write(iout,'("Input DEFINEDET line (includes frozen orbs):")')
            do run = 1, inum_runs
                write(6,'("definedet ")', advance='no')
                if (allocated(frozen_orb_list)) then
                    allocate(tmp_ni(nel_pre_freezing))
                    tmp_ni(1:nel) = frozen_orb_reverse_map(ProjEDet(:,run))
                    if (nel /= nel_pre_freezing) &
                        tmp_ni(nel+1:nel_pre_freezing) = frozen_orb_list
                    call sort(tmp_ni)
                    call writeDefDet(tmp_ni, nel_pre_freezing)
!                    do i = 1, nel_pre_freezing
!                        write(6, '(i3," ")', advance='no') tmp_ni(i)
!                    end do
                    deallocate(tmp_ni)
                else
                   call writeDefDet(ProjEDet(:,run), nel)
!                    do i = 1, nel
!                        write(6, '(i3," ")', advance='no') ProjEDet(i, run)
!                    end do
                end if
                do i = 1, nel
                    full_orb = ProjEDet(i, run)
                    if (allocated(frozen_orb_list)) &
                        full_orb = full_orb  + count(frozen_orb_list <= ProjEDet(i, run))
                end do
                write(iout,*)

                if (.not. tReplicaReferencesDiffer) exit
            end do

            write(iout,*) ""
            write(iout,"(A,I10,A)") "Most occupied ",counter," determinants as excitations from reference: "
            write(iout,*)
            if(lenof_sign.eq.1) then
                if(tHPHF) then
                    write(iout,"(A)") " Excitation   ExcitLevel   Seniority    Walkers    Amplitude    Init?   Proc  Spin-Coup?"
                else
                    write(iout,"(A)") " Excitation   ExcitLevel   Seniority    Walkers    Amplitude    Init?   Proc"
                endif
            else
#ifdef CMPLX_
                if(tHPHF) then
                    write(iout,"(A)") " Excitation   ExcitLevel Seniority  Walkers(Re)   Walkers(Im)  Weight   &
                                        &Init?(Re)   Init?(Im)   Proc  Spin-Coup?"
                else
                    write(iout,"(A)") " Excitation   ExcitLevel Seniority   Walkers(Re)   Walkers(Im)  Weight   &
                                        &Init?(Re)   Init?(Im)   Proc"
                endif
#else
                ! output the weight of every replica, and do not only assume
                ! it is a complex run
                write(format_string, '(a,i0,a,a,i0,a)') &
                    '(3a11,', lenof_sign, 'a11,', 'a13,', lenof_sign,'a9,a)'
                ! Walkers(replica) Amplitude(replica) Init?(replica)
                allocate(walker_string(lenof_sign))
!                 allocate(amplitude_string(lenof_sign))
                allocate(init_string(lenof_sign))

                do i = 1, lenof_sign
                    write(walker_string(i), '(a,i0,a)') "Walkers(", i, ")"
!                     write(amplitude_string(i), '(a,i0,a)') "Amplitude(", i, ")"
                    write(init_string(i), '(a,i0,a)') "Init?(", i, ")"
                end do

                write(header, format_string) "Excitation ", "ExcitLevel ", "Seniority ", &
                    walker_string, "Amplitude ", init_string, "Proc "

                if (tHPHF) then
                    header = trim(header) // " Spin-Coup?"
                end if

                write(iout, '(a)') trim(header)

#endif
            endif
            do i=1,counter
!                call WriteBitEx(iout,iLutRef,GlobalLargestWalkers(:,i),.false.)
                call WriteDetBit(iout,GlobalLargestWalkers(:,i),.false.)
                Excitlev=FindBitExcitLevel(iLutRef(:,1),GlobalLargestWalkers(:,i),nEl,.true.)
                write(iout,"(I5)",advance='no') Excitlev
                nopen=count_open_orbs(GlobalLargestWalkers(:,i))
                write(iout,"(I5)",advance='no') nopen
                call extract_sign(GlobalLargestWalkers(:,i),SignCurr)
                do j=1,lenof_sign
                    write(iout,"(G16.7)",advance='no') SignCurr(j)
                enddo
#ifdef CMPLX_
                HighSign = sqrt(sum(abs(SignCurr(1::2)))**2 + sum(abs(SignCurr(2::2)))**2)
#else
                HighSign=sum(real(abs(SignCurr),dp))
#endif
                if(tHPHF.and.(.not.TestClosedShellDet(GlobalLargestWalkers(:,i)))) then
                    !Weight is proportional to (nw/sqrt(2))**2
                    write(iout,"(F9.5)",advance='no') ((HighSign/sqrt(2.0_dp))/norm )
                else
                    write(iout,"(F9.5)",advance='no') (HighSign/norm)
                endif
                do j=1,lenof_sign
                    if(.not.tTruncInitiator) then
                        write(iout,"(A3)",advance='no') 'Y'
                    else
                        if(test_flag(GlobalLargestWalkers(:,i),get_initiator_flag(j))) then
                            write(iout,"(A3)",advance='no') 'Y'
                        else
                            write(iout,"(A3)",advance='no') 'N'
                        endif
                    endif
                enddo
                if(tHPHF.and.(.not.TestClosedShellDet(GlobalLargestWalkers(:,i)))) then
                    write(iout,"(I7)",advance='no') GlobalProc(i)
                    write(iout,"(A3)") "*"
                else
                    write(iout,"(I7)") GlobalProc(i)
                endif
            enddo

            if(tHPHF) then
                write(iout,"(A)") " * = Spin-coupled function implicitly has time-reversed determinant with same weight."
            endif

            deallocate(GlobalLargestWalkers,GlobalProc)
            write(iout,*) ""
        endif

        deallocate(LargestWalkers)

        contains

          subroutine writeDefDet(defdet, numEls)
            implicit none
            integer, intent(in) :: numEls
            integer, intent(in) :: defdet(:)
            logical :: nextInRange, previousInRange

            do i = 1, numEls
               ! if the previous orbital is in the same contiguous range
               if(i.eq.1) then
                  ! for the first one, there is no previous one
                  previousInRange = .false.
               else
                  previousInRange = defdet(i).eq.defdet(i-1)+1
               endif

               ! if the following orbital is in the same contiguous range
               if(i.eq.numEls) then
                  ! there is no following orbital
                  nextInRange = .false.
               else
                  nextInRange = defdet(i).eq.defdet(i+1)-1
               endif
               ! there are three cases that need output:

               ! the last orbital of a contigous range of orbs
               if(previousInRange .and. .not.nextInRange) then
                  write(bufEnd,'(i3)') defdet(i)
                  lenEnd = len_trim(bufEnd)
                  bufStart(lenStart+2:lenStart+lenEnd+1) = adjustl(trim(bufEnd))
                  write(iout,'(A7)',advance='no') trim(adjustl(bufStart))
               ! the first orbital of a contiguous range of orbs
               else if(.not.previousInRange .and. nextInRange) then
                  write(bufStart,'(i3)') defdet(i)
                  lenStart = len_trim(bufStart)
                  bufStart(lenStart+1:lenStart+1) = "-"
               ! and an orbital not in any range
               else if(.not.previousInRange .and. .not.nextInRange) then
                  write(iout,'(i3," ")', advance='no') defdet(i)
               endif
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
         endif

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
      endif
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
      endif
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
        endif
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
            write(iout,*) "WARNING: Empty energy histogram of acceptance rates"
        endif
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
                write (iout, '(a)') '*WARNING* - Number of particles/determinants &
                                 &has increased to over 95% of allotted memory. &
                                 &Errors imminent. Increase MEMORYFACWALKERS, or reduce rate of growth.'
            else
                write (iout, '(a)') '*WARNING* - Number of particles/determinants &
                                 &has increased to over 95% of allotted memory. &
                                 &Errors imminent. Increase MEMORYFACPART, or reduce rate of growth.'
            endif
#else
            if(tMolpro) then
                write (iout,*) '*WARNING* - Number of particles/determinants &
                                 &has increased to over 95% of allotted memory on task ', iProcIndex, '. &
                                 &Errors imminent. Increase MEMORYFACWALKERS, or reduce rate of growth.'
            else
                write (iout,*) '*WARNING* - Number of particles/determinants &
                                 &has increased to over 95% of allotted memory on task ', iProcIndex, '. &
                                 &Errors imminent. Increase MEMORYFACPART, or reduce rate of growth.'
            endif
#endif
            call neci_flush(iout)
        end if

        ! Are ony of the sublists near the end of their alloted space?
        if (nNodes > 1) then
            do i = 0, nNodes-1
                rat = real(ValidSpawnedList(i) - InitialSpawnedSlots(i),dp) /&
                             real(InitialSpawnedSlots(1), dp)
                if (rat > 0.95_dp) then
#ifdef DEBUG_
                    if(tMolpro) then
                        write (iout, '(a)') '*WARNING* - Highest processor spawned &
                                         &particles has reached over 95% of allotted memory.&
                                         &Errors imminent. Increase MEMORYFACSPAWNED, or reduce spawning rate.'
                    else
                        write (iout, '(a)') '*WARNING* - Highest processor spawned &
                                         &particles has reached over 95% of allotted memory.&
                                         &Errors imminent. Increase MEMORYFACSPAWN, or reduce spawning rate.'
                    endif
#else
                    if(tMolpro) then
                        write (iout,*) '*WARNING* - Highest processor spawned &
                                         &particles has reached over 95% of allotted memory on task ',iProcIndex,' .&
                                         &Errors imminent. Increase MEMORYFACSPAWNED, or reduce spawning rate.'
                    else
                        write (iout,*) '*WARNING* - Highest processor spawned &
                                         &particles has reached over 95% of allotted memory on task ',iProcIndex,' .&
                                         &Errors imminent. Increase MEMORYFACSPAWN, or reduce spawning rate.'
                    endif
#endif
                    call neci_flush(iout)
                endif
            enddo
        else
            rat = real(ValidSpawnedList(0), dp) / real(MaxSpawned, dp)
            if (rat > 0.95_dp) then
#ifdef DEBUG_
                if(tMolpro) then
                    write (iout, '(a)') '*WARNING* - Highest processor spawned &
                                     &particles has reached over 95% of allotted memory.&
                                     &Errors imminent. Increase MEMORYFACSPAWNED, or reduce spawning rate.'
                else
                    write (iout, '(a)') '*WARNING* - Highest processor spawned &
                                     &particles has reached over 95% of allotted memory.&
                                     &Errors imminent. Increase MEMORYFACSPAWN, or reduce spawning rate.'
                endif
#else
                if(tMolpro) then
                    write (iout,*) '*WARNING* - Highest processor spawned &
                                     &particles has reached over 95% of allotted memory on task ',iProcIndex,' .&
                                     &Errors imminent. Increase MEMORYFACSPAWNED, or reduce spawning rate.'
                else
                    write (iout,*) '*WARNING* - Highest processor spawned &
                                     &particles has reached over 95% of allotted memory on task ',iProcIndex,' .&
                                     &Errors imminent. Increase MEMORYFACSPAWN, or reduce spawning rate.'
                endif
#endif
                call neci_flush(iout)
            endif
         endif

    end subroutine end_iteration_print_warn

    subroutine print_frequency_histogram_spec
        ! this routine is the adapted version to print out the single and
        ! double(para/anti) histograms and a combined one of all of them
        integer :: iunit, i, max_size, old_size
        character(255) :: filename
        integer :: all_frequency_bins_s(n_frequency_bins)
        integer :: all_frequency_bins_d(n_frequency_bins)
        integer :: all_frequency_bins_p(n_frequency_bins)
        integer :: all_frequency_bins_a(n_frequency_bins)
        integer :: all_frequency_bins(n_frequency_bins)
        integer :: all_frequency_bins_2(n_frequency_bins)
        integer :: all_frequency_bins_2_d(n_frequency_bins)
        integer :: all_frequency_bins_3(n_frequency_bins)
        integer :: all_frequency_bins_3_d(n_frequency_bins)
        integer :: all_frequency_bins_4(n_frequency_bins)
        real(dp) :: step_size, norm
!         real(dp), allocatable :: all_frequency_bounds(:)
        character(*), parameter :: this_routine = "print_frequency_histogram_spec"

        integer(int64) :: sum_all

        ! this is only called in the 4ind weighted or GUGA case so singles
        ! are always there so do them first
        ! why does it crash here? and not on my laptop... compiler issue i
        ! guess.. too much memory requested... lol
        all_frequency_bins_s = 0
        call MPIReduce(frequency_bins_singles, MPI_SUM, all_frequency_bins_s)

        if (iProcIndex == 0) then
            max_size = size(all_frequency_bins_s)

            iunit = get_free_unit()
            call get_unique_filename("frequency_histogram_singles",.true.,.true.,&
                1, filename)
            open(iunit, file = filename, status = "unknown")

            do i = 1, max_size
                if (all_frequency_bins_s(i) == 0) cycle
                write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                write(iunit, "(i12)") all_frequency_bins_s(i)
            end do

            close(iunit)

!             deallocate(all_frequency_bounds)

        end if

!         deallocate(all_frequency_bins)

        ! then dependent if it is guga or 4ind print out remaining
        if (tGen_4ind_weighted .or. tGen_4ind_2) then
            ! do para first
            call comm_frequency_histogram_spec(size(frequency_bins_para), &
                frequency_bins_para, all_frequency_bins_p)

            if (iProcIndex == 0) then
                max_size = size(all_frequency_bins_p)

                iunit = get_free_unit()
                call get_unique_filename("frequency_histogram_para", .true., &
                    .true., 1, filename)
                open(iunit, file = filename, status = "unknown")

                do i = 1, max_size
                    if (all_frequency_bins_p(i) == 0) cycle
                    write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                    write(iunit, "(i12)") all_frequency_bins_p(i)
                end do
                close(iunit)

!                 deallocate(all_frequency_bounds)
            end if

            ! then anti
!             deallocate(all_frequency_bins)

            call comm_frequency_histogram_spec(size(frequency_bins_anti), &
                frequency_bins_anti, all_frequency_bins_a)

            if (iProcIndex == 0) then
                max_size = size(all_frequency_bins_a)

                iunit = get_free_unit()
                call get_unique_filename("frequency_histogram_anti", .true., &
                    .true., 1, filename)
                open(iunit, file = filename, status = "unknown")

                do i = 1, max_size
                    if (all_frequency_bins_a(i) == 0) cycle
                    write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                    write(iunit, "(i12)") all_frequency_bins_a(i)
                end do
                close(iunit)

!                 deallocate(all_frequency_bounds)
            end if
!             deallocate(all_frequency_bins)

            ! i also want to add up all the bins for a final output
            if (iProcIndex == 0) then
                max_size = max(size(all_frequency_bins_s), size(all_frequency_bins_p), &
                    size(all_frequency_bins_a))

!                 allocate(all_frequency_bins(max_size))
                all_frequency_bins = 0
                all_frequency_bins(1:size(all_frequency_bins_s)) = all_frequency_bins_s

                all_frequency_bins(1:size(all_frequency_bins_p)) = &
                    all_frequency_bins(1:size(all_frequency_bins_p)) + &
                    all_frequency_bins_p

                all_frequency_bins(1:size(all_frequency_bins_a)) = &
                    all_frequency_bins(1:size(all_frequency_bins_a)) + &
                    all_frequency_bins_a

                ! and also need the max bounds

                if (.not. any(all_frequency_bins < 0)) then
                    iunit = get_free_unit()
                    call get_unique_filename("frequency_histogram", .true., .true., &
                        1, filename)
                    open(iunit, file = filename, status = "unknown")

                    do i = 1, max_size
                        if (all_frequency_bins(i) == 0) cycle
                        write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                        write(iunit, "(i12)") all_frequency_bins(i)
                    end do
                    close(iunit)
                else
                    write(iout,*) "Integer overflow in all_frequency_bins!"
                    write(iout,*) "Do no print it!"
                end if

                ! also print out a normed frequency histogram to better
                ! compare runs with different length
                sum_all = sum(all_frequency_bins)
                if (.not. sum_all < 0) then
                    ! we have a int overflow..
                    ! how to deal with that?? hm..
                    norm = real(sum_all,dp)

                    iunit = get_free_unit()
                    call get_unique_filename("frequency_histogram_normed", .true., &
                        .true., 1, filename)
                    open(iunit, file = filename, status = "unknown")

                    ! and change x and y axis finally
                    do i = 1, max_size
                        if (all_frequency_bins(i) == 0) cycle
                        write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                        write(iunit, "(f16.7)") real(all_frequency_bins(i),dp) / norm
                    end do
                    close(iunit)

                else
                    write(iout,*) "Integer overflow in normed frequency histogram!"
                    write(iout,*) "Do no print it!"

                end if
!                 deallocate(all_frequency_bins)
            end if

        else if (tGen_sym_guga_mol .or. (tgen_guga_crude .and. .not. t_new_real_space_hubbard)) then
            ! do only doubles for now in the guga case

            call comm_frequency_histogram_spec(size(frequency_bins_doubles), &
                frequency_bins_doubles, all_frequency_bins_d)

            if (iProcIndex == 0) then
                max_size = size(all_frequency_bins_d)

                iunit = get_free_unit()
                call get_unique_filename("frequency_histogram_doubles", .true., &
                    .true., 1, filename)
                open(iunit, file = filename, status = "unknown")

                do i = 1, max_size
                    if (all_frequency_bins_d(i) == 0) cycle
                    write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                    write(iunit, "(i12)") all_frequency_bins_d(i)
                end do
                close(iunit)

!                 deallocate(all_frequency_bounds)

                ! and also do the add up from singles and doubles
                max_size = max(size(all_frequency_bins_s), size(all_frequency_bins_d))

!                 allocate(all_frequency_bins(max_size))
                all_frequency_bins = 0

                all_frequency_bins(1:size(all_frequency_bins_s)) = all_frequency_bins_s

                all_frequency_bins(1:size(all_frequency_bins_d)) = &
                    all_frequency_bins(1:size(all_frequency_bins_d)) + &
                    all_frequency_bins_d

                if (.not. any(all_frequency_bins < 0)) then

                    iunit = get_free_unit()
                    call get_unique_filename("frequency_histogram", .true., .true., &
                        1, filename)
                    open(iunit, file = filename, status = "unknown")

                    do i = 1, max_size
                        if (all_frequency_bins(i) == 0) cycle
                        write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                        write(iunit, "(i12)") all_frequency_bins(i)
                    end do
                    close(iunit)
                else
                    write(iout,*) "Integer overflow in all_frequency_bins!"
                    write(iout,*) "Do no print it!"
                end if

                ! also print out a normed frequency histogram to better
                ! compare runs with different length
                sum_all = sum(all_frequency_bins)
                if (.not. sum_all < 0) then

                    norm = real(sum_all,dp)

                    iunit = get_free_unit()
                    call get_unique_filename("frequency_histogram_normed", .true., &
                        .true., 1, filename)
                    open(iunit, file = filename, status = "unknown")

                    ! and change x and y axis finally
                    do i = 1, max_size
                        if (all_frequency_bins(i) == 0) cycle
                        write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                        write(iunit, "(f16.7)") real(all_frequency_bins(i),dp) / norm
                    end do
                    close(iunit)

                else
                    write(iout,*) "Integer overflow in normed frequency histogram!"
                    write(iout,*) "Do no print it!"
                end if

!                 deallocate(all_frequency_bins)

            end if

        else if (tGen_nosym_guga) then
            ! here a lot of stuff has to be printed and also dependent on
            ! t_consider_diff_bias..
            ! singles are already printed
            call comm_frequency_histogram_spec(size(frequency_bins_type2), &
                frequency_bins_type2, all_frequency_bins_2)

            if (iProcIndex == 0) then
                max_size = size(all_frequency_bins_2)

                iunit = get_free_unit()
                call get_unique_filename("frequency_histogram_type2", .true., &
                    .true., 1, filename)
                open(iunit, file = filename, status = "unknown")

                do i = 1, max_size
                    if (all_frequency_bins_2(i) == 0) cycle
                    write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                    write(iunit, "(i12)") all_frequency_bins_2(i)
                end do
                close(iunit)

!                 deallocate(all_frequency_bounds)
            end if

            call comm_frequency_histogram_spec(size(frequency_bins_type3), &
                frequency_bins_type3, all_frequency_bins_3)

            if (iProcIndex == 0) then
                max_size = size(all_frequency_bins_3)

                iunit = get_free_unit()
                call get_unique_filename("frequency_histogram_type3", .true., &
                    .true., 1, filename)
                open(iunit, file = filename, status = "unknown")

                do i = 1, max_size
                    if (all_frequency_bins_3(i) == 0) cycle
                    write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                    write(iunit, "(i12)") all_frequency_bins_3(i)
                end do
                close(iunit)
            end if

            call comm_frequency_histogram_spec(size(frequency_bins_type4), &
                frequency_bins_type4, all_frequency_bins_4)

            if (iProcIndex == 0) then
                max_size = size(all_frequency_bins_4)

                iunit = get_free_unit()
                call get_unique_filename("frequency_histogram_type4", .true., &
                    .true., 1, filename)
                open(iunit, file = filename, status = "unknown")

                do i = 1, max_size
                    if (all_frequency_bins_4(i) == 0) cycle
                    write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                    write(iunit, "(i12)") all_frequency_bins_4(i)
                end do
                close(iunit)

                ! determine the final max size for the summed up histogram
                max_size = max(size(all_frequency_bins_s), size(all_frequency_bins_2), &
                    size(all_frequency_bins_3), size(all_frequency_bins_4))

            end if

            if (t_consider_diff_bias) then
                call comm_frequency_histogram_spec(size(frequency_bins_type2_diff), &
                    frequency_bins_type2_diff, all_frequency_bins_2_d)

                if (iProcIndex == 0) then
                    max_size = size(all_frequency_bins_2_d)

                    iunit = get_free_unit()
                    call get_unique_filename("frequency_histogram_type2_diff", &
                        .true., .true., 1, filename)
                    open(iunit, file = filename, status = "unknown")

                    do i = 1, max_size
                        if (all_frequency_bins_2_d(i) == 0) cycle
                        write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                        write(iunit, "(i12)") all_frequency_bins_2_d(i)
                    end do
                    close(iunit)

                end if

                call comm_frequency_histogram_spec(size(frequency_bins_type3_diff), &
                    frequency_bins_type3_diff, all_frequency_bins_3_d)

                if (iProcIndex == 0) then
                    max_size = size(all_frequency_bins_3_d)

                    iunit = get_free_unit()
                    call get_unique_filename("frequency_histogram_type3_diff", &
                        .true., .true., 1, filename)
                    open(iunit, file = filename, status = "unknown")

                    do i = 1, max_size
                        if (all_frequency_bins_3_d(i) == 0) cycle
                        write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                        write(iunit, "(i12)") all_frequency_bins_3_d(i)
                    end do
                    close(iunit)
                end if

                max_size = max(size(all_frequency_bins_s), size(all_frequency_bins_2), &
                    size(all_frequency_bins_2_d), size(all_frequency_bins_3), &
                    size(all_frequency_bins_3_d), size(all_frequency_bins_4))

            end if

            if (iProcIndex == 0) then
!                 allocate(all_frequency_bins(max_size))
                all_frequency_bins = 0

                all_frequency_bins(1:size(all_frequency_bins_s)) = all_frequency_bins_s

                all_frequency_bins(1:size(all_frequency_bins_2)) = &
                    all_frequency_bins(1:size(all_frequency_bins_2)) + &
                    all_frequency_bins_2

                all_frequency_bins(1:size(all_frequency_bins_3)) = &
                    all_frequency_bins(1:size(all_frequency_bins_3)) + &
                    all_frequency_bins_3

                all_frequency_bins(1:size(all_frequency_bins_4)) = &
                    all_frequency_bins(1:size(all_frequency_bins_4)) + &
                    all_frequency_bins_4

                if (t_consider_diff_bias) then

                    all_frequency_bins(1:size(all_frequency_bins_2_d)) = &
                        all_frequency_bins(1:size(all_frequency_bins_2_d)) + &
                        all_frequency_bins_2_d

                    all_frequency_bins(1:size(all_frequency_bins_3_d)) = &
                        all_frequency_bins(1:size(all_frequency_bins_3_d)) + &
                        all_frequency_bins_3_d

                end if

                if (.not. any(all_frequency_bins < 0)) then
                    iunit = get_free_unit()
                    call get_unique_filename("frequency_histogram", .true., .true., &
                        1, filename)
                    open(iunit, file = filename, status = "unknown")

                    do i = 1, max_size
                        if (all_frequency_bins(i) == 0) cycle
                        write(iunit, "(f16.7)") frq_step_size * i
                        write(iunit, "(i12)", advance = "no") all_frequency_bins(i)
                    end do
                    close(iunit)
                else
                    write(iout,*) "Integer overflow in all_frequency_bins"
                    write(iout,*) "Do not print it!"
                end if

                ! also print out a normed frequency histogram to better
                ! compare runs with different length
                sum_all = sum(all_frequency_bins)

                if (.not. sum_all < 0) then
                    norm = real(sum_all,dp)

                    iunit = get_free_unit()
                    call get_unique_filename("frequency_histogram_normed", .true., &
                        .true., 1, filename)
                    open(iunit, file = filename, status = "unknown")

                    ! and change x and y axis finally
                    do i = 1, max_size
                        if (all_frequency_bins(i) == 0) cycle
                        write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                        write(iunit, "(f16.7)") real(all_frequency_bins(i),dp) / norm
                    end do
                    close(iunit)
                else
                    write(iout,*) "Integer overflow in normed frequency histogram!"
                    write(iout,*) "Do not print it!"
                end if

!                 deallocate(all_frequency_bins)
            end if
        end if
!
    end subroutine print_frequency_histogram_spec
!
    subroutine print_frequency_histogram
        ! routine to write a file with the H_ij/pgen ratio frequencies
        integer :: iunit, i, max_size, old_size
        character(255) :: filename
        integer :: all_frequency_bins(n_frequency_bins)
        real(dp) :: step_size
        real(dp), allocatable :: all_frequency_bounds(:)
        real(dp) :: norm
        integer :: sum_all
        ! i can test how to deal with the MPI stuff here to get the same
        ! results as in the single runs
        ! first i need the maximum length of all processors
        call comm_frequency_histogram(all_frequency_bins)

        if (iProcIndex == 0) then
            max_size = size(all_frequency_bins)

            iunit = get_free_unit()
            call get_unique_filename("frequency_histogram",.true.,.true.,1,filename)
            open(iunit, file = filename, status = "unknown")


            do i = 1, max_size
                write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                write(iunit, "(i12)") all_frequency_bins(i)
            end do

            close(iunit)

            ! and print out normed frequency histogram if possible
            sum_all = sum(all_frequency_bins)
            if (.not. sum_all < 0) then
                ! we have a int overflow..
                ! how to deal with that?? hm..
                norm = real(sum_all,dp)

                iunit = get_free_unit()
                call get_unique_filename("frequency_histogram_normed", .true., &
                    .true., 1, filename)
                open(iunit, file = filename, status = "unknown")

                ! and change x and y axis finally
                do i = 1, max_size
                    write(iunit, "(f16.7)", advance = "no") frq_step_size * i
                    write(iunit, "(f16.7)") real(all_frequency_bins(i),dp) / norm
                end do
                close(iunit)

            else
                write(iout,*) "Integer overflow in normed frequency histogram!"
                write(iout,*) "Do no print it!"

            end if

!             deallocate(all_frequency_bounds)
        end if

    end subroutine print_frequency_histogram

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
         OutputHii = hphf_diag_helement (ProjEDet(:,1), &
              iLutRef(:,1))
      else
         OutputHii = get_helement (ProjEDet(:,1), &
              ProjEDet(:,1), 0)
      endif


    end subroutine getProjEOffset

end module

