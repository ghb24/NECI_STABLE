#include "macros.h"
module fcimc_output

    use SystemData, only: nel, tHPHF, tFixLz, tMolpro, tMolproMimic, MolproID
    use LoggingData, only: tLogComplexPops, tMCOutput, tCalcInstantS2, &
                           tCalcInstantS2Init, instant_s2_multiplier_init, &
                           instant_s2_multiplier, tPrintFCIMCPsi, &
                           iWriteHistEvery, tDiagAllSpaceEver, OffDiagMax, &
                           OffDiagBinRange, tCalcVariationalEnergy, &
                           iHighPopWrite, tLogEXLEVELStats, tWriteConflictLvls, &
                           maxInitExLvlWrite, AllInitsPerExLvl
    use hist_data, only: Histogram, AllHistogram, InstHist, AllInstHist, &
                         BeforeNormHist, iNoBins, BinRange, HistogramEnergy, &
                         AllHistogramEnergy
    use CalcData, only: tTruncInitiator, tTrialWavefunction, tReadPops, &
                        DiagSft, tSpatialOnlyHash, tOrthogonaliseReplicas, &
                        StepsSft, tPrintReplicaOverlaps, tStartTrialLater, tEN2, &
                        tGlobalInitFlag, t_truncate_spawns, tTimedDeaths
    use DetBitOps, only: FindBitExcitLevel, count_open_orbs, EncodeBitDet, &
                         TestClosedShellDet
    use IntegralsData, only: frozen_orb_list, frozen_orb_reverse_map, &
                             nel_pre_freezing
    use DetCalcData, only: det, fcidets, ReIndex, NDet, NRow, HAMIL, LAB
    use bit_reps, only: decode_bit_det, test_flag, extract_sign, get_initiator_flag
    use semi_stoch_procs, only: return_most_populated_states
    use bit_rep_data, only: niftot, nifd, flag_initiator
    use hist, only: calc_s_squared_star, calc_s_squared
    use fcimc_helper, only: LanczosFindGroundE
    use Determinants, only: write_det
    use adi_data, only: AllCoherentDoubles, AllIncoherentDets, nRefs, &
         ilutRefAdi, tAdiActive, nConnection, AllConnection
    use rdm_data, only: en_pert_main
    use Parallel_neci
    use FciMCData
    use constants
    use sort_mod
    use util_mod

    implicit none

contains

    SUBROUTINE WriteFciMCStatsHeader()
        integer i, j, k, run, offset
        character(256) label
        character(32) tchar_r, tchar_i, tchar_j, tchar_k
        character(17) trunc_caption

        IF(iProcIndex.eq.root) THEN
!Print out initial starting configurations
            WRITE(iout,*) ""
            IF(tTruncInitiator) THEN
               WRITE(initiatorstats_unit,"(A2,A17,15A23)", advance = 'no') &
                    "# ","1.Step","2.TotWalk","3.Annihil","4.Died", &
                    & "5.Born","6.TotUniqDets",&
                    &               "7.InitDets","8.NonInitDets","9.InitWalks","10.NonInitWalks","11.AbortedWalks", &
                    "12. Removed Dets",  "13. Initiator Proj.E"
               offset = 13            
               if(tTrialWavefunction .or. tStartTrialLater) then
                  write(initiatorstats_unit,"(A)", advance = 'no') &
                  "14. TrialNumerators (inits)   15. TrialDenom (inits)"
                  offset = 15
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

#ifdef __CMPLX
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

#elif defined(__DOUBLERUN)
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
                  &31.|Semistoch|/|Psi|   32.PartsDiffProc"
           if (tTrialWavefunction .or. tStartTrialLater) then 
                  write(fcimcstats_unit2, "(A)", advance = 'no') &
                  "  33.TrialNumerator  34.TrialDenom  35.TrialOverlap"
           end if

           write(fcimcstats_unit2, "()", advance = 'yes')
#endif
#ifndef __CMPLX
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
                  &32.|Semistoch|/|Psi|  33.MaxCycSpawn"
           if (tTrialWavefunction .or. tStartTrialLater) then 
              write(fcimcstats_unit, "(A)", advance = 'no') &
                   "  34.TrialNumerator  35.TrialDenom  36.TrialOverlap"
              trunc_caption = "  37. TruncWeight"
           else
              trunc_caption = "  34. TruncWeight"
           end if
           if(t_truncate_spawns) write(fcimcstats_unit, "(A)", advance = 'no') &
                trunc_caption

           write(fcimcstats_unit, "()", advance = 'yes')

#endif
            
        ENDIF

    END SUBROUTINE WriteFciMCStatsHeader

    subroutine WriteFCIMCStats()
        INTEGER :: i, j, run
        real(dp),dimension(inum_runs) :: FracFromSing
        real(dp) :: projE(inum_runs)

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

        !To prevent /0 problems
        do run=1,inum_runs
            if(AllNoBorn(run).ne.0) then
                FracFromSing(run)=real(AllSpawnFromSing(run),dp) / real(AllNoBorn(run),dp)
            else
                FracFromSing(run)=0.0_dp
            endif
        enddo

        if (iProcIndex == root) then

#ifdef __CMPLX
            write(fcimcstats_unit,"(I12,5G16.7,8G18.9e3,&
                                  &G13.5,I12,G13.5,G17.5,I13,G13.5,8G18.9e3,I13,&

                                  &g16.7)",advance='no') &
                Iter + PreviousCycles, &                !1.
                DiagSft, &                              !2.
                AllTotParts(1) - AllTotPartsOld(1), &   !3.
                AllTotParts(2) - AllTotPartsOld(2), &   !4.
                AllTotParts(1), &                       !5.
                AllTotParts(2), &                       !6.
                real(ProjectionE, dp), &                !7.     real \sum[ nj H0j / n0 ]
                aimag(projectionE), &                   !8.     Im   \sum[ nj H0j / n0 ]
                real(proje_iter, dp), &                 !9.     
                aimag(proje_iter), &                    !10.
                real(proje_iter,dp) + Hii, &            !11.
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
                real((AllHFCyc*conjg(AllHFCyc)),dp), &     !23 |n0|^2  denominator for both calcs
                real((AllENumCyc*conjg(AllHFCyc)),dp), &   !24. Re[\sum njH0j]xRe[n0]+Im[\sum njH0j]xIm[n0]   No div by StepsSft
                aimag(AllENumCyc*conjg(AllHFCyc)), &       !25.Im[\sum njH0j]xRe[n0]-Re[\sum njH0j]xIm[n0]   since no physicality
                sqrt(sum(AllNoatHF**2)) / norm_psi, & !26
                norm_psi, &                           !27
                curr_S2, &                            !28
                PartsDiffProc, &                      !29
                all_max_cyc_spawn                     !30.
                if (tTrialWavefunction .or. tStartTrialLater) then
                    write(fcimcstats_unit, "(7(1X,es18.11))", advance = 'no') &
                    (tot_trial_numerator(1) / StepsSft), &              ! 31. 32
                    (tot_trial_denom(1) / StepsSft), &                  ! 33. 34
                    abs((tot_trial_denom(1) / (norm_psi(1)*StepsSft))), &  ! 35.
                    tot_trial_numerator(1)/tot_trial_denom(1)           ! 36. 37.
                end if
                write(fcimcstats_unit, "()", advance = 'yes')

            if(tMCOutput) then
                write (iout, "(I12,13G16.7,2I12,G13.5)") &
                    Iter + PreviousCycles, &
                    DiagSft, &
                    AllTotParts(1) - AllTotPartsOld(1), &
                    AllTotParts(2) - AllTotPartsOld(2), &
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
               write(initiatorstats_unit,"(I12,4G16.7,3I20,7G16.7,2I20,1G16.7)", &
                    advance = 'no')&
                   Iter + PreviousCycles, sum(AllTotParts), &
                   AllAnnihilated(1), AllNoDied(1), AllNoBorn(1), AllTotWalkers,&
                   AllNoInitDets(1), AllNoNonInitDets(1), AllNoInitWalk(1), &
                   AllNoNonInitWalk(1),AllNoAborted(1), AllNoRemoved(1), &
                   inits_proje_iter(1) + Hii
               if(tTrialWavefunction .or. tStartTrialLater) &
                    write(initiatorstats_unit, "(2G16.7)", advance = 'no') &
                    tot_init_trial_numerator(1)/StepsSft, tot_init_trial_denom(1)/StepsSft
               do j = 1, maxInitExLvlWrite
                  write(initiatorstats_unit,'(1I20)', advance ='no') AllInitsPerExLvl(j)
               end do
               write(initiatorstats_unit,'()', advance = 'yes')
            endif
            if (tLogComplexPops) then
                write (complexstats_unit,"(I12,6G16.7)") &
                    Iter + PreviousCycles, DiagSft, DiagSftRe, DiagSftIm, &
                    sum(AllTotParts), AllTotParts(1), AllTotParts(lenof_sign)
            endif
#elif defined(__DOUBLERUN)
            write(fcimcstats_unit2,"(i12,7g16.7,5g18.9e3,g13.5,i12,g13.5,g17.5,&
                                   &i13,g13.5,4g18.9e3,1X,2(es18.11,1X),5g18.9e3,&
                                   &i13,2g16.7)",advance = 'no') &
                Iter + PreviousCycles, &                   ! 1.
                DiagSft(2), &                              ! 2.
                AllTotParts(2) - AllTotPartsOld(2), &      ! 3.
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
                proje_iter(2) + Hii, &                     ! 23.
                (AllHFCyc(2) / StepsSft), &                ! 24.
                (AllENumCyc(2) / StepsSft), &              ! 25.
                AllNoatHF(2) / norm_psi(2), &              ! 26.
                norm_psi(2), &                             ! 27.
                curr_S2(2), curr_S2_init(2), &             ! 28, 29.
                AbsProjE(2), &                             ! 30.
                PartsDiffProc, &                           ! 31.
                norm_semistoch(2)/norm_psi(2), &           ! 32.
                all_max_cyc_spawn                          ! 33.
                if (tTrialWavefunction .or. tStartTrialLater) then
                    write(fcimcstats_unit2, "(3(1X,es17.10))", advance = 'no') &
                    (tot_trial_numerator(2) / StepsSft), &
                    (tot_trial_denom(2) / StepsSft), &
                    abs(tot_trial_denom(2) / (norm_psi(2)*StepsSft))
                end if
                
                write(fcimcstats_unit2, "()", advance = 'yes')
#endif
#ifndef __CMPLX

            write(fcimcstats_unit,"(i12,7g16.7,5g18.9e3,g13.5,i12,g13.5,g17.5,&
                                   &i13,g13.5,4g18.9e3,1X,2(es18.11,1X),5g18.9e3,&
                                   &i13,2g16.7)",advance = 'no') &
                Iter + PreviousCycles, &                   ! 1.
                DiagSft(1), &                              ! 2.
                AllTotParts(1) - AllTotPartsOld(1), &      ! 3.
                AllGrowRate(1), &                          ! 4.
                AllTotParts(1), &                          ! 5.
                AllAnnihilated(1), &                       ! 6.
                AllNoDied(1), &                            ! 7.
                AllNoBorn(1), &                            ! 8.
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
                proje_iter(1) + Hii, &                     ! 23.
                (AllHFCyc(1) / StepsSft), &                ! 24.
                (AllENumCyc(1) / StepsSft), &              ! 25.
                AllNoatHF(1) / norm_psi(1), &              ! 26.
                norm_psi(1), &                             ! 27.
                curr_S2(1), curr_S2_init(1), &             ! 28, 29.
                AbsProjE(1), &                             ! 30.
                PartsDiffProc, &                           ! 31.
                norm_semistoch(1)/norm_psi(1), &           ! 32.
                all_max_cyc_spawn                          ! 33.
                if (tTrialWavefunction .or. tStartTrialLater) then
                    write(fcimcstats_unit, "(3(1X,es18.11))", advance = 'no') &
                    (tot_trial_numerator(1) / StepsSft), &              ! 34.
                    (tot_trial_denom(1) / StepsSft), &                  ! 35.
                    abs((tot_trial_denom(1) / (norm_psi(1)*StepsSft)))  ! 36.
                end if
                if(t_truncate_spawns) then
                   write(fcimcstats_unit, "(1X,es18.11)", advance = 'no') AllTruncatedWeight
                endif
                write(fcimcstats_unit, "()", advance = 'yes')

            if(tMCOutput) then
                write (iout, "(I12,10G16.7)", advance = 'no') &
                    Iter + PreviousCycles, &
                    DiagSft(1), &
                    AllTotParts(1) - AllTotPartsOld(1), &
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
               write(initiatorstats_unit,"(I12,4G16.7,3I20,7G16.7,2I20,1G16.7)", &
                    advance = 'no')&
                   Iter + PreviousCycles, AllTotParts(1), &
                   AllAnnihilated(1), AllNoDied(1), AllNoBorn(1), AllTotWalkers,&
                   AllNoInitDets(1), AllNoNonInitDets(1), AllNoInitWalk(1), &
                   AllNoNonInitWalk(1),AllNoAborted(1), AllNoRemoved(1), &
                   inits_proje_iter(1) + Hii
               if(tTrialWavefunction .or. tStartTrialLater) &
                    write(initiatorstats_unit, "(2G16.7)", advance = 'no') &
                    tot_init_trial_numerator(1)/StepsSft, tot_init_trial_denom(1)/StepsSft
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
        integer :: i, ierr
        logical :: exists

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

            ! If we are doing a normal calculation, move existing fciqmc_stats
            ! files so that they are not overwritten, and then create a new one
            inquire(file=filename, exist=exists)
            if (exists) then

                ! Loop until we find an available spot to move the existing
                ! file to.
                i = 1
                do while(exists)
                    write(num, '(i12)') i
                    filename2 = trim(adjustl(filename)) // "." // &
                                trim(adjustl(num))
                    inquire(file=filename2, exist=exists)
                    if (i > 10000) &
                        call stop_all(t_r, 'Error finding free fciqmc_stats.*')
                    i = i + 1
                end do

                ! Move the file
                call rename(filename, filename2)

            end if

            ! And finally open the file
            open(funit, file=filename, status='unknown', iostat=ierr)

        end if

    end subroutine open_create_stats

    subroutine write_unoccstats(initial)
      logical, intent(in), optional :: initial
      
      type(write_state_t), save :: state_ud
      integer :: p
      logical, save :: inited = .false.
      character(5) :: tmpc

      if(present(initial)) then
         state_ud%init = initial
      else
         state_ud%init = .false.
      endif

      ! only root prints the info on the unocc dets    
      if(iProcIndex == root) then
         if(.not. inited) then
            call open_state_file('unoccupied_stats',state_ud)
            inited = .true.
         endif

         call write_padding_init(state_ud)

         call stats_out(state_ud, .false., Iter + PreviousCycles, 'Iter')
         call stats_out(state_ud, .false., AllNUnoccDets, 'Unocc Dets')
         do p = 1, maxHoleExLvlWrite
            ! write the number of conflicts of this excitation lvl
            write(tmpc,('(i5)')) p
            call stats_out(state_ud, .false., AllHolesByExLvl(p), 'nUnocc (ex = '//&
                 trim(adjustl(tmpc)) // ")")
         end do

         write(state_ud%funit,*)
         ! flush output
         call neci_flush(state_ud%funit)
      end if
      
    end subroutine write_unoccstats

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
        type(write_state_t), save :: state_cl
        logical, save :: inited = .false.
        character(5) :: tmpc, tmpc2
        integer :: p, q
        logical :: init
       
        ! Provide default 'initial' option
        if (present(initial)) then
            state%init = initial
            if (tTruncInitiator) state_i%init = initial
            if(tWriteConflictLvls) state_cl%init = initial
        else
            state%init = .false.
            if (tTruncInitiator) state_i%init = .false.
            if(tWriteConflictLvls) state_cl%init = .false.
        end if

        ! If the output file hasn't been opened yet, then create it.
        if (iProcIndex == Root .and. .not. inited) then
           
           call open_state_file('fciqmc_stats',state)
           ! For the initiator stats file here:
           if (tTruncInitiator) call open_state_file('initiator_stats',state_i)

           if(tWriteConflictLvls) call open_state_file('conflict_stats',state_cl)

           inited = .true.
        end if

        ! ------------------------------------------------
        ! This is where any calculation that needs multiple nodes should go
        ! ------------------------------------------------
        ! ------------------------------------------------
    
        if (iProcIndex == root) then

            ! Only do the actual outputting on the head node.

            call write_padding_init(state)
            call write_padding_init(state_i)
            call write_padding_init(state_cl)

            ! And output the actual data!
            state%cols = 0
            state%cols_mc = 0
            state%mc_out = tMCOutput
            call stats_out(state,.true., iter + PreviousCycles, 'Iter.')
            if (.not. tOrthogonaliseReplicas) then
                call stats_out(state,.true., sum(abs(AllTotParts)), 'Tot. parts')
                call stats_out(state,.true., sum(abs(AllNoatHF)), 'Tot. ref')
#ifdef __CMPLX
                call stats_out(state,.true., real(proje_iter_tot), 'Re Proj. E')
                call stats_out(state,.true., aimag(proje_iter_tot), 'Im Proj. E')
#else
                call stats_out(state,.true., proje_iter_tot, 'Proj. E (cyc)')
#endif
                call stats_out(state,.true., sum(DiagSft / inum_runs), 'Shift. (cyc)')
                call stats_out(state,.false., sum(AllNoBorn), 'No. born')
                call stats_out(state,.false., sum(AllNoDied), 'No. died')
                call stats_out(state,.false., sum(AllAnnihilated), 'No. annihil')
!!            call stats_out(state,.false., AllGrowRate(1), 'Growth fac.')
!!            call stats_out(state,.false., AccRat(1), 'Acc. rate')
#ifdef __CMPLX
                call stats_out(state,.true., real(proje_iter_tot) + Hii, &
                               'Tot. Proj. E')
                call stats_out(state,.true., aimag(proje_iter_tot) + Hii, &
                               'Tot. Proj. E')
#else
                call stats_out(state,.true., proje_iter_tot + Hii, &
                               'Tot. Proj. E')
#endif
            end if
            call stats_out(state,.true., AllTotWalkers, 'Dets occ.')
            call stats_out(state,.true., nspawned_tot, 'Dets spawned')

            call stats_out(state,.true., IterTime, 'Iter. time')
            call stats_out(state,.false., TotImagTime, 'Im. time')

            ! Put the conditional columns at the end, so that the column
            ! numbers of the data are as stable as reasonably possible (for
            ! people who want to use gnuplot/not analyse column headers too
            ! frequently).
            ! This also makes column contiguity on resumes as likely as
            ! possible.
            
            ! if we truncate walkers, print out the total truncated weight here
            if(t_truncate_spawns) call stats_out(state, .false., AllTruncatedWeight, &
                 'trunc. Weight')

            ! If we are running multiple (replica) simulations, then we
            ! want to record the details of each of these
#ifdef __PROG_NUMRUNS
            do p = 1, inum_runs
                write(tmpc, '(i5)') p
                call stats_out (state, .false., AllTotParts(p), &
                                'Parts (' // trim(adjustl(tmpc)) // ")")
                call stats_out (state, .false., AllNoatHF(p), &
                                'Ref (' // trim(adjustl(tmpc)) // ")")
                call stats_out (state, .false., DiagSft(p) + Hii, &
                                'Shift (' // trim(adjustl(tmpc)) // ")")
#ifdef __CMPLX
                call stats_out (state, .false., real(proje_iter(p) + Hii), &
                                'Tot ProjE real (' // trim(adjustl(tmpc)) // ")")
                call stats_out (state, .false., aimag(proje_iter(p) + Hii), &
                                'Tot ProjE imag (' // trim(adjustl(tmpc)) // ")")

                call stats_out (state, .false., real(AllHFCyc(p) / StepsSft), &
                                'ProjE Denom real (' // trim(adjustl(tmpc)) // ")")
                call stats_out (state, .false., aimag(AllHFCyc(p) / StepsSft), &
                                'ProjE Denom imag (' // trim(adjustl(tmpc)) // ")")

                call stats_out (state, .false., &
                                real((AllENumCyc(p) + Hii*AllHFCyc(p))) / StepsSft,&
                                'ProjE Num real (' // trim(adjustl(tmpc)) // ")")
                call stats_out (state, .false., &
                                aimag((AllENumCyc(p) + Hii*AllHFCyc(p))) / StepsSft,&
                                'ProjE Num imag (' // trim(adjustl(tmpc)) // ")")
                if (tTrialWavefunction .or. tStartTrialLater) then
                    call stats_out (state, .false., &
                                    real(tot_trial_numerator(p) / StepsSft), &
                                    'TrialE Num real (' // trim(adjustl(tmpc)) // ")")
                    call stats_out (state, .false., &
                                    aimag(tot_trial_numerator(p) / StepsSft), &
                                    'TrialE Num imag (' // trim(adjustl(tmpc)) // ")")

                    call stats_out (state, .false., &
                                    real(tot_trial_denom(p) / StepsSft), &
                                    'TrialE Denom real (' // trim(adjustl(tmpc)) // ")")
                    call stats_out (state, .false., &
                                    aimag(tot_trial_denom(p) / StepsSft), &
                                    'TrialE Denom imag (' // trim(adjustl(tmpc)) // ")")
                end if
#else
                call stats_out (state, .false., proje_iter(p) + Hii, &
                                'Tot ProjE (' // trim(adjustl(tmpc)) // ")")
                call stats_out (state, .false., AllHFCyc(p) / StepsSft, &
                                'ProjE Denom (' // trim(adjustl(tmpc)) // ")")
                call stats_out (state, .false., &
                                (AllENumCyc(p) + Hii*AllHFCyc(p)) / StepsSft,&
                                'ProjE Num (' // trim(adjustl(tmpc)) // ")")
                if (tTrialWavefunction .or. tStartTrialLater) then
                    call stats_out (state, .false., &
                                    tot_trial_numerator(p) / StepsSft, &
                                    'TrialE Num (' // trim(adjustl(tmpc)) // ")")
                    call stats_out (state, .false., &
                                    tot_trial_denom(p) / StepsSft, &
                                    'TrialE Denom (' // trim(adjustl(tmpc)) // ")")
                end if
#endif


                call stats_out (state, .false., &
                                AllNoBorn(p), &
                                'Born (' // trim(adjustl(tmpc)) // ")")
                call stats_out (state, .false., &
                                AllNoDied(p), &
                                'Died (' // trim(adjustl(tmpc)) // ")")
                call stats_out (state, .false., &
                                AllAnnihilated(p), &
                                'Annihil (' // trim(adjustl(tmpc)) // ")")
                call stats_out (state, .false., &
                                AllNoAtDoubs(p), &
                                'Doubs (' // trim(adjustl(tmpc)) // ")")
            end do

            ! Print overlaps between replicas at the end.
            do p = 1, inum_runs
                write(tmpc, '(i5)') p
                if (tPrintReplicaOverlaps) then
                    do q = p+1, inum_runs
                        write(tmpc2, '(i5)') q
#ifdef __CMPLX
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

#endif

            if (tEN2) call stats_out(state,.true., en_pert_main%ndets_all, 'EN2 Dets.')

            if (tTruncInitiator) then
                call stats_out(state_i, .false., Iter + PreviousCycles, 'Iter.')
                call stats_out(state_i, .false., AllTotWalkers, 'TotDets.')
                call stats_out(state_i, .false., AllNoInitsConflicts/inum_runs,&
                     'Inc. Inits (normal)')
                call stats_out(state_i, .false., AllNoSIInitsConflicts/inum_runs,&
                     'Inc. Inits (SI)')
                call stats_out(state_i, .false., AllAvSigns, 'Replica-averaged Sign')
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

            ! gather sign conflict statistics
            if(tWriteConflictLvls) then
               call stats_out(state_cl, .false., Iter + PreviousCycles, 'Iter')
               call stats_out(state_cl, .false., AllNoConflicts, 'confl. Dets')
               do p = 1, maxConflictExLvl
                  ! write the number of conflicts of this excitation lvl
                  write(tmpc,('(i5)')) p
                  call stats_out(state_cl, .false., AllConflictExLvl(p), 'confl. (ex = '//&
                       trim(adjustl(tmpc)) // ")")
               end do
            endif
     
            ! And we are done
            write(state%funit, *)
            if (tTruncInitiator) write(state_i%funit, *)
            if(tWriteConflictLvls) write(state_cl%funit,*)
            if (tMCOutput) write(iout, *)
            call neci_flush(state%funit)
            if (tTruncInitiator) call neci_flush(state_i%funit)
            if(tWriteConflictLvls) call neci_flush(state_cl%funit)
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
        integer :: i, ms, tempni(1:nel)
        real(dp) :: totWalkPopByMsReal(nel+1), totWalkPopByMsImag(nel+1), &
                    tempSign(lenof_sign)
        
        do i=1,TotWalkers
            call extract_sign(WalkVecDets(:,i),TempSign)
            call decode_bit_det(TempnI, WalkVecDets(:,i))
            ms = sum(get_spin_pn(Tempni(1:nel)))
            walkPopByMsReal(1+nel/2+ms/2) = walkPopByMsReal(1+nel/2+ms/2)+abs(TempSign(1))
#ifdef __CMPLX
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
#ifdef __CMPLX
        norm2=SQRT(sum(norm1))
#else
        norm1=SQRT(norm1)
#endif
        WRITE(iout,*) "Total FCIMC Wavefuction normalisation:",norm1
        do i=1,Det
            do j=1,lenof_sign
#ifdef __CMPLX
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
#ifdef __CMPLX
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
        abstr=''
        write(abstr,'(I12)') Iter
        abstr='SpawnHist-'//adjustl(abstr)
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
#ifdef __CMPLX
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
#ifdef __CMPLX
                    AllHistogram(j,i)=AllHistogram(j,i)/norm_c
                    AllInstHist(j,i)=AllInstHist(j,i)/norm1_c
                    IF(norm2_c.ne.0.0_dp) THEN
                        AllInstAnnihil(j,i)=AllInstAnnihil(j,i)/norm2_c
                    ENDIF
                    IF(norm3_c.ne.0.0_dp) THEN
                        AllAvAnnihil(j,i)=AllAvAnnihil(j,i)/norm3_c
                    ENDIF
#else
                    AllHistogram(j,i)=AllHistogram(j,i)/norm(j)
                    AllInstHist(j,i)=AllInstHist(j,i)/norm1(j)
                    IF(norm2(j).ne.0.0_dp) THEN
                        AllInstAnnihil(j,i)=AllInstAnnihil(1,i)/norm2(j)
                    ENDIF
                    IF(norm3(j).ne.0.0_dp) THEN
                    AllAvAnnihil(j,i)=AllAvAnnihil(j,i)/norm3(j)
                    ENDIF
#endif
                enddo
            enddo
            
            io1 = get_free_unit()
            OPEN(io1,FILE=abstr,STATUS='UNKNOWN')

            abstr=''
            write(abstr,'(I12)') Iter-iWriteHistEvery
            abstr='Energies-'//adjustl(abstr)

            abstr2=''
            write(abstr2,'(I12)') Iter
            abstr2='Energies-'//adjustl(abstr2)
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
#ifdef __CMPLX
                IF(AllHFCyc(1).eq.0.0_dp) THEN
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft,AllERead,SUM(AllTotPartsOld)
                ELSE
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft,AllENumCyc/AllHFCyc,SUM(AllTotPartsOld)
                ENDIF
#else
                IF(AllHFCyc(1).eq.0.0_dp) THEN
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft(1),AllERead,AllTotPartsOld(1)
                ELSE
                    WRITE(io2,"(I13,3G25.16)") Iter,DiagSft(1),AllENumCyc(1)/AllHFCyc(1),AllTotPartsOld(1)
                ENDIF
#endif
                CLOSE(io2)
                CLOSE(io3)

            ELSE
                OPEN(io2,FILE=abstr2,STATUS='UNKNOWN')
#ifdef __CMPLX
                WRITE(io2,"(I13,3G25.16)") Iter,DiagSft,AllENumCyc/AllHFCyc,SUM(AllTotPartsOld)
#else
                WRITE(io2,"(I13,3G25.16)") Iter,DiagSft(1),AllENumCyc(1)/AllHFCyc(1),AllTotPartsOld(1)
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
#ifdef __CMPLX
                    WRITE(io1,"(I13,6G25.16)") i, AllHistogram(1,i), sum(norm), &
                          AllInstHist(1,i), AllInstAnnihil(1,i), &
                          AllAvAnnihil(1,i), sum(norm1)
#else
                    WRITE(io1,"(I13,6G25.16)") i, AllHistogram(1,i), norm(1), &
                          AllInstHist(1,i), AllInstAnnihil(1,i), &
                          AllAvAnnihil(1,i), norm1(1)
#endif
                ENDIF
                IF(AllHistogram(1,i).ne.0.0_dp) Tot_No_Unique_Dets = Tot_No_Unique_Dets + 1
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

!                write(6,*) ReIndex(1:10)    !Must be in increasing order
!                write(6,*) HOrderedHist(:)

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
                        if(AllHistogram(1,i).ne.0.0_dp) then
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
                        if(AllHistogram(1,i).ne.0.0_dp) then
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
        abstr=''
        write(abstr,'(I12)') Iter
        abstr='HamilHist-'//adjustl(abstr)
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
        character(len=*), parameter :: t_r='PrintHighPops'

        !Allocate memory to hold highest iHighPopWrite determinants
        allocate(LargestWalkers(0:NIfTot,iHighPopWrite),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"error allocating here")

        ! Return the most populated states in CurrentDets on *this* processor only.
        call return_most_populated_states(iHighPopWrite, LargestWalkers, norm)

        call MpiSum(norm,allnorm)
        if(iProcIndex.eq.Root) norm=sqrt(allnorm)

!        write(iout,*) "Highest weighted dets on this process:"
!        do i=1,iHighPopWrite
!            write(iout,*) LargestWalkers(:,i)
!        enddo

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
                    
#ifdef __CMPLX
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
#ifdef __CMPLX
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
                    do i = 1, nel_pre_freezing
                        write(6, '(i3," ")', advance='no') tmp_ni(i)
                    end do
                    deallocate(tmp_ni)
                else
                    do i = 1, nel
                        write(6, '(i3," ")', advance='no') ProjEDet(i, run)
                    end do
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
                if(tHPHF) then
                    write(iout,"(A)") " Excitation   ExcitLevel Seniority  Walkers(Re)   Walkers(Im)  Weight   &
                                        &Init?(Re)   Init?(Im)   Proc  Spin-Coup?"
                else
                    write(iout,"(A)") " Excitation   ExcitLevel Seniority   Walkers(Re)   Walkers(Im)  Weight   &
                                        &Init?(Re)   Init?(Im)   Proc"
                endif
            endif
            do i=1,counter
!                call WriteBitEx(iout,iLutRef,GlobalLargestWalkers(:,i),.false.)
                call WriteDetBit(iout,GlobalLargestWalkers(:,i),.false.)
                Excitlev=FindBitExcitLevel(iLutRef(:,1),GlobalLargestWalkers(:,i),nEl)
                write(iout,"(I5)",advance='no') Excitlev
                nopen=count_open_orbs(GlobalLargestWalkers(:,i))
                write(iout,"(I5)",advance='no') nopen
                call extract_sign(GlobalLargestWalkers(:,i),SignCurr)
                do j=1,lenof_sign
                    write(iout,"(G16.7)",advance='no') SignCurr(j)
                enddo
#ifdef __CMPLX
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

    END SUBROUTINE PrintHighPops
            
    subroutine end_iteration_print_warn (totWalkersNew)
        
        ! Worker function for PerformFciMCycPar. Prints warnings about 
        ! particle blooms and memory usage.
        integer, intent(in) :: totWalkersNew
        integer :: i
        real(dp) :: rat

        ! Too many particles?
        rat = real(TotWalkersNew,dp) / real(MaxWalkersPart,dp)
        if (rat > 0.95_dp) then
#ifdef __DEBUG
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
                write (*,*) '*WARNING* - Number of particles/determinants &
                                 &has increased to over 95% of allotted memory on task ', iProcIndex, '. &
                                 &Errors imminent. Increase MEMORYFACWALKERS, or reduce rate of growth.'
            else
                write (*,*) '*WARNING* - Number of particles/determinants &
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
#ifdef __DEBUG
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
                        write (*,*) '*WARNING* - Highest processor spawned &
                                         &particles has reached over 95% of allotted memory on task ',iProcIndex,' .&
                                         &Errors imminent. Increase MEMORYFACSPAWNED, or reduce spawning rate.'
                    else
                        write (*,*) '*WARNING* - Highest processor spawned &
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
#ifdef __DEBUG
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
                    write (*,*) '*WARNING* - Highest processor spawned &
                                     &particles has reached over 95% of allotted memory on task ',iProcIndex,' .&
                                     &Errors imminent. Increase MEMORYFACSPAWNED, or reduce spawning rate.'
                else
                    write (*,*) '*WARNING* - Highest processor spawned &
                                     &particles has reached over 95% of allotted memory on task ',iProcIndex,' .&
                                     &Errors imminent. Increase MEMORYFACSPAWN, or reduce spawning rate.'
                endif
#endif
                call neci_flush(iout)
            endif
        endif

    end subroutine end_iteration_print_warn 

end module
