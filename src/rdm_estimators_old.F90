module rdm_estimators_old

    ! Routines related to the calculation of properties from RDMs, such as the
    ! RDM energy, spin, and dipole moments.

    use bit_rep_data, only: NIfTot
    use constants

    implicit none

contains

    subroutine rdm_output_wrapper_old(rdm, rdm_label, est)

        use FciMCData, only: tFinalRDMEnergy, Iter, IterRDMStart, PreviousCycles
        use LoggingData, only: tRDMInstEnergy, tWriteMultRDMs, IterWriteRDMs
        use LoggingData, only: tWrite_RDMs_to_read, tWrite_normalised_RDMs
        use LoggingData, only: tWriteSpinFreeRDM
        use Parallel_neci, only: iProcIndex, MPISumAll
        use rdm_data, only: rdm_t, rdm_estimates_old_t, tOpenShell
        use rdm_temp, only: Finalise_2e_RDM, calc_2e_norms, Write_out_2RDM
        use rdm_temp, only: Write_spinfree_RDM

        type(rdm_t), intent(inout) :: rdm
        integer, intent(in) :: rdm_label
        type(rdm_estimates_old_t), intent(inout) :: est

        ! Normalise, make Hermitian, etc.
        call Finalise_2e_RDM(rdm)

        if (iProcIndex == 0) then
            ! Calculate the normalisations.
            call calc_2e_norms(rdm, est%Norm_2RDM_Inst, est%Norm_2RDM, est%Trace_2RDM)

            ! There's no need to explicitly make the RDM hermitian here, as the
            ! integrals are already hermitian -- when we calculate the energy,
            ! it comes out in the wash.
            
            ! Print out the relevant 2-RDMs.
            if (tFinalRDMEnergy .or. (tWriteMultRDMs .and. (mod((Iter+PreviousCycles-IterRDMStart)+1,IterWriteRDMs) .eq. 0))) then

                ! Only ever want to print the 2-RDMs (for reading in) at the end.
                if (tFinalRDMEnergy .and. tWrite_RDMs_to_read) call Write_out_2RDM(rdm, rdm_label, est%Norm_2RDM, .false., .false.)

                ! This writes out the normalised, hermitian 2-RDMs.
                ! IMPORTANT NOTE: We assume that we want tMake_Herm=.true. here.
                if (tWrite_normalised_RDMs) call Write_out_2RDM(rdm, rdm_label, est%Norm_2RDM, .true., .true.)

                if (tWriteSpinFreeRDM) call Write_spinfree_RDM(rdm, rdm_label, est%Norm_2RDM)

             end if

            call Calc_Energy_from_RDM(rdm, est%Norm_2RDM, est%Norm_2RDM_Inst, est%Trace_2RDM_normalised, est%RDMEnergy, &
                                      est%RDMEnergy1, est%RDMEnergy2, est%RDMEnergy_Inst)

            ! Calculate the instantaneous estimate of <S^2> using the 2RDM.
            est%spin_est = calc_2rdm_spin_estimate(rdm, est%Norm_2RDM_Inst)
        end if

        ! Zero all of the instantaneous RDM arrays and data.
        if (tRDMInstEnergy) then
            rdm%aaaa(:,:) = 0.0_dp
            rdm%abab(:,:) = 0.0_dp
            rdm%abba(:,:) = 0.0_dp

            if (tOpenShell)then
                rdm%bbbb(:,:) = 0.0_dp
                rdm%baba(:,:) = 0.0_dp
                rdm%baab(:,:) = 0.0_dp
            end if
        end if

    end subroutine rdm_output_wrapper_old

    subroutine write_rdm_estimates_old(est, est_old)

        use FciMCData, only: tFinalRDMEnergy, Iter, PreviousCycles
        use LoggingData, only: tRDMInstEnergy
        use rdm_data, only: rdm_estimates_t, rdm_estimates_unit
        use rdm_data, only: rdm_estimates_old_t
        use util_mod, only: int_fmt

        type(rdm_estimates_t), intent(in) :: est
        type(rdm_estimates_old_t), intent(in) :: est_old(:)

        integer :: i

        if (tRDMInstEnergy) then
            write(rdm_estimates_unit, '(1x,i13)', advance='no') Iter+PreviousCycles
            do i = 1, est%nrdms
                write(rdm_estimates_unit, '(6(3x,es20.13))', advance='no') &
                    est_old(i)%RDMEnergy_Inst, est_old(i)%spin_est, 1.0_dp/est_old(i)%Norm_2RDM_Inst, &
                    est%energy_tot_num(i), est%spin_num(i), est%norm(i)
            end do
            write(rdm_estimates_unit,'()')

        else

            write(rdm_estimates_unit, '(1x,i13)') Iter+PreviousCycles
            do i = 1, est%nrdms
                write(rdm_estimates_unit, '(3(3x,es20.13))', advance='no') &
                    est_old(i)%RDMEnergy, est_old(i)%spin_est, 1.0_dp/est_old(i)%Norm_2RDM_Inst
            end do
            write(rdm_estimates_unit, '()')
        end if

        call neci_flush(rdm_estimates_unit)

        if (tFinalRDMEnergy) then
            do i = 1, est%nrdms
                write(6,'(1x,"FINAL ESTIMATES FOR RDM",1X,'//int_fmt(i)//',":",)') i
                write(6,'(1x,"Trace of 2-el-RDM before normalisation:",1x,es17.10)') est%trace(i)
                write(6,'(1x,"Trace of 2-el-RDM after normalisation:",1x,es17.10)') est%trace(i)/est%norm(i)
                write(6,'(1x,"Energy contribution from the 1-RDM:",1x,es17.10)') est%energy_1_num(i)/est%norm(i)
                write(6,'(1x,"Energy contribution from the 2-RDM:",1x,es17.10)') est%energy_2_num(i)/est%norm(i)
                write(6,'(1x,"*TOTAL ENERGY* CALCULATED USING THE *REDUCED DENSITY MATRICES*:",1x,es20.13,/)') &
                    est%energy_tot_num(i)/est%norm(i)
            end do
            close(rdm_estimates_unit)
        end if

    end subroutine write_rdm_estimates_old

    subroutine Calc_Energy_from_RDM(rdm, Norm_2RDM, Norm_2RDM_Inst, Trace_2RDM_normalised, RDMEnergy, &
                                    RDMEnergy1, RDMEnergy2, RDMEnergy_Inst)

        ! This routine takes the 1 electron and 2 electron reduced density
        ! matrices and calculated the energy they give.
        ! The equation for the energy is as follows:
        !
        !   E = Tr(h1 1RDM) + 1/2 Tr(h2 2RDM)
        !
        ! where h1 are the 2 index integrals, and h2 the 4 index integrals.
        ! The traces, are given by:
        !
        !   Tr(h1 1RDM) = Sum_i,j [ h1(i,j) 1RDM(j,i) ]
        !   Tr(h2 2RDM) = Sum_i,j;k,l [ h2(i,j;k,l) 2RDM(k,l;i,j) ]
        use FciMCData, only: tFinalRDMEnergy
        use global_utilities, only: set_timer, halt_timer
        use IntegralsData, only: umat
        use LoggingData, only: tRDMInstEnergy
        use rdm_data, only: RDMEnergy_Time, tOpenShell
        use rdm_data, only: rdm_t
        use RotateOrbsData, only: SpatOrbs
        use SystemData, only: tStoreSpinOrbs, ecore
        use UMatCache, only: UMatInd

        type(rdm_t), intent(inout) :: rdm
        real(dp), intent(in) :: Norm_2RDM, Norm_2RDM_Inst
        real(dp), intent(out) :: Trace_2RDM_normalised, RDMEnergy, RDMEnergy1, RDMEnergy2, RDMEnergy_Inst

        integer :: i, j, a, b, Ind1_aa, Ind1_ab, Ind2_aa, Ind2_ab
        integer :: iSpin, jSpin
        real(dp) :: Coul, Exch
        real(dp) :: Coul_aaaa, Coul_bbbb, Coul_abab, Coul_baba
        real(dp) :: Exch_aaaa, Exch_bbbb, Exch_abba, Exch_baab

        call set_timer(RDMEnergy_Time, 30)

        Trace_2RDM_normalised = 0.0_dp
        RDMEnergy_Inst = 0.0_dp
        RDMEnergy1 = 0.0_dp
        RDMEnergy2 = 0.0_dp
        RDMEnergy = 0.0_dp
    
        do i = 1, SpatOrbs
            iSpin = 2 * i

            do j = i, SpatOrbs
                jSpin = 2 * j

                Ind1_aa = ( ( (j-2) * (j-1) ) / 2 ) + i
                Ind1_ab = ( ( (j-1) * j ) / 2 ) + i

                do a = 1, SpatOrbs

                    ! Adding in contributions effectively from the 1-RDM (although these are calculated 
                    ! from the 2-RDM).
                    call calc_1RDM_and_1RDM_energy(rdm, i, j, a, iSpin,jSpin, Norm_2RDM, &
                                                   RDMEnergy_Inst, RDMEnergy1)
                    
                    do b = a, SpatOrbs

                        Ind2_aa = ( ( (b-2) * (b-1) ) / 2 ) + a
                        Ind2_ab = ( ( (b-1) * b ) / 2 ) + a

                        ! UMAT in chemical notation.
                        ! In spin or spatial orbitals.
                        Coul = real(UMAT(UMatInd(i, j, a, b)), dp)
                        Exch = real(UMAT(UMatInd(i, j, b, a)), dp)
                        
                        if ((i .ne. j) .and. (a .ne. b)) then
                            ! Cannot get i=j or a=b contributions in aaaa.

                            if (tStoreSpinOrbs)then
                                Coul_aaaa = real(UMAT(UMatInd(2*i, 2*j, 2*a, 2*b)),dp)
                                Coul_bbbb = real(UMAT(UMatInd(2*i-1, 2*j-1, 2*a-1, 2*b-1)),dp)
                                Exch_aaaa = real(UMAT(UMatInd(2*i, 2*j, 2*b, 2*a)),dp)
                                Exch_bbbb = real(UMAT(UMatInd(2*i-1, 2*j-1, 2*b-1, 2*a-1)),dp)     

                                if (tRDMInstEnergy) then 
                                    RDMEnergy_Inst = RDMEnergy_Inst + (rdm%aaaa(Ind1_aa,Ind2_aa) &
                                                      *( Coul_aaaa - Exch_aaaa ) )
                                    RDMEnergy_Inst = RDMEnergy_Inst + (rdm%bbbb(Ind1_aa,Ind2_aa) &
                                                      *  ( Coul_bbbb - Exch_bbbb ) )
                                end if

                                RDMEnergy2 = RDMEnergy2 + ( rdm%aaaa_full(Ind1_aa,Ind2_aa) &
                                                      * Norm_2RDM * ( Coul_aaaa - Exch_aaaa ) )
                                RDMEnergy2 = RDMEnergy2 + ( rdm%bbbb_full(Ind1_aa,Ind2_aa) &
                                                      * Norm_2RDM * ( Coul_bbbb - Exch_bbbb ) )

                            else 

                                if (tRDMInstEnergy) then
                                    RDMEnergy_Inst = RDMEnergy_Inst + (rdm%aaaa(Ind1_aa,Ind2_aa) &
                                                                * ( Coul - Exch ) )
                                end if

                                RDMEnergy2 = RDMEnergy2 + ( rdm%aaaa_full(Ind1_aa,Ind2_aa) &
                                                        * Norm_2RDM * ( Coul - Exch ) ) 

                                if (tOpenShell) then
                                    if (tRDMInstEnergy) then 
                                        RDMEnergy_Inst = RDMEnergy_Inst + (rdm%bbbb(Ind1_aa,Ind2_aa) &
                                                      * ( Coul - Exch ) )
                                    end if

                                    RDMEnergy2 = RDMEnergy2 + ( rdm%bbbb_full(Ind1_aa,Ind2_aa) &
                                                      * Norm_2RDM * ( Coul - Exch ) )
                                end if 

                            end if  


                            if (Ind1_aa .eq. Ind2_aa) then
                                Trace_2RDM_normalised = Trace_2RDM_normalised + &
                                                    rdm%aaaa_full(Ind1_aa,Ind2_aa) * Norm_2RDM
                                if (tOpenShell) Trace_2RDM_normalised = Trace_2RDM_normalised + &
                                                    rdm%bbbb_full(Ind1_aa,Ind2_aa) * Norm_2RDM
                            end if

                            ! For abab cases, coul element will be non-zero, exchange zero.

                            if (tStoreSpinOrbs) then
                                Coul_abab = real(UMAT(UMatInd(2*i, 2*j-1, 2*a, 2*b-1)), dp)
                                Coul_baba = real(UMAT(UMatInd(2*i-1, 2*j, 2*a-1, 2*b)), dp)

                                call neci_flush(6)

                                if (tRDMInstEnergy) then
                                    RDMEnergy_Inst = RDMEnergy_Inst + ( rdm%abab(Ind1_ab,Ind2_ab) &
                                                                    *  Coul_abab )
                                end if
                                RDMEnergy2 = RDMEnergy2 + ( rdm%abab_full(Ind1_ab,Ind2_ab) &
                                                        * Norm_2RDM * Coul_abab ) 

                                if (tRDMInstEnergy) then
                                    RDMEnergy_Inst = RDMEnergy_Inst + ( rdm%baba(Ind1_ab,Ind2_ab) &
                                                                    *  Coul_baba )
                                end if
                                RDMEnergy2 = RDMEnergy2 + ( rdm%baba_full(Ind1_ab,Ind2_ab) &
                                                        * Norm_2RDM * Coul_baba ) 

                            else 
                                if (tRDMInstEnergy) then
                                    RDMEnergy_Inst = RDMEnergy_Inst + ( rdm%abab(Ind1_ab,Ind2_ab) &
                                                                    *  Coul )
                                    if (tOpenShell) RDMEnergy_Inst = RDMEnergy_Inst + ( rdm%baba(Ind1_ab,Ind2_ab) &
                                                                    *  Coul )
                                end if

                                RDMEnergy2 = RDMEnergy2 + ( rdm%abab_full(Ind1_ab,Ind2_ab) &
                                                        * Norm_2RDM * Coul ) 
                                if (tOpenShell) RDMEnergy2 = RDMEnergy2 + ( rdm%baba_full(Ind1_ab,Ind2_ab) &
                                                        * Norm_2RDM * Coul) 

                            end if 

                            if (Ind1_ab .eq. Ind2_ab) then
                                Trace_2RDM_normalised = Trace_2RDM_normalised + &
                                                    rdm%abab_full(Ind1_ab,Ind2_ab) * Norm_2RDM
                                if (tOpenShell) Trace_2RDM_normalised = Trace_2RDM_normalised + &
                                                    rdm%baba_full(Ind1_ab,Ind2_ab) * Norm_2RDM
                            end if

                            ! For abba cases, coul element will be zero, exchange non-zero.

                            if (tStoreSpinOrbs) then
                                Exch_abba = real(UMAT(UMatInd(2*i, 2*j-1, 2*b, 2*a-1)), dp)
                                Exch_baab = real(UMAT(UMatInd(2*i-1, 2*j, 2*b-1, 2*a)), dp)

                                if (tRDMInstEnergy) then 
                                    RDMEnergy_Inst = RDMEnergy_Inst - ( rdm%abba(Ind1_aa,Ind2_aa) &
                                                                    * Exch_abba )
                                    RDMEnergy_Inst = RDMEnergy_Inst - ( rdm%baab(Ind1_aa,Ind2_aa) &
                                                                    * Exch_baab )
                                end if

                                RDMEnergy2 = RDMEnergy2 - ( rdm%abba_full(Ind1_aa,Ind2_aa) &
                                                        * Norm_2RDM * Exch_abba ) 
                                RDMEnergy2 = RDMEnergy2 - ( rdm%baab_full(Ind1_aa,Ind2_aa) &
                                                        * Norm_2RDM * Exch_baab ) 

                            else 

                                if (tRDMInstEnergy) then 
                                    RDMEnergy_Inst = RDMEnergy_Inst - ( rdm%abba(Ind1_aa,Ind2_aa) &
                                                                    * Exch )
                                    if (tOpenShell) RDMEnergy_Inst = RDMEnergy_Inst - ( rdm%baab(Ind1_aa,Ind2_aa) &
                                                                    *  Exch )
                                end if

                                RDMEnergy2 = RDMEnergy2 - ( rdm%abba_full(Ind1_aa,Ind2_aa) &
                                                        * Norm_2RDM * Exch ) 
                                if (tOpenShell) RDMEnergy2 = RDMEnergy2 - ( rdm%baab_full(Ind1_aa,Ind2_aa) &
                                                        * Norm_2RDM * Exch ) 

                           end if 

                       else if ( (i .eq. j) .or. (a .eq. b) )then
                            ! i = j or a = b
                            ! abab has both abab and abba elements in them effectively.
                            ! half will have non-zero coul, and half non-zero exchange.
                            ! For abab/baba Exch = 0, and for abba/baab Coul=0
                            ! abba/baab saved in abab/baba. Changes the sign. 
                            if (tStoreSpinOrbs) then
                                Coul_abab = real(UMAT(UMatInd(2*i, 2*j-1, 2*a, 2*b-1)), dp)
                                Coul_baba = real(UMAT(UMatInd(2*i-1, 2*j, 2*a-1, 2*b)), dp)

                                if ( (i .eq. j) .and. (a .eq. b) ) then
                                    ! This term is saved in abab only
                                    Exch_abba = real(UMAT(UMatInd(2*i, 2*j-1, 2*b, 2*a-1)), dp)

                                    if (tRDMInstEnergy) then
                                        RDMEnergy_Inst = RDMEnergy_Inst + 0.5_dp * rdm%abab(Ind1_ab,Ind2_ab) &
                                                                          * (Coul_abab+Exch_abba)
                               
                                    end if 

                                    RDMEnergy2 = RDMEnergy2 + 0.5_dp * rdm%abab_full(Ind1_ab,Ind2_ab) &
                                                              * Norm_2RDM * (Coul_abab+Exch_abba)

                                else if (i .eq. j) then
                                    ! i = j : Swap first indeces to get abba/baab terms
                                    ! abba saved in baba, baab saved in abab (sign changes)
                                    Exch_abba = real(UMAT(UMatInd(2*j, 2*i-1, 2*b, 2*a-1)), dp)
                                    Exch_baab = real(UMAT(UMatInd(2*j-1, 2*i, 2*b-1, 2*a)), dp)

                                    if (tRDMInstEnergy) then
                                        RDMEnergy_Inst = RDMEnergy_Inst + 0.5_dp * rdm%abab(Ind1_ab,Ind2_ab) &
                                                                          *  (Coul_abab+Exch_baab)
                                        RDMEnergy_Inst = RDMEnergy_Inst + 0.5_dp * rdm%baba(Ind1_ab,Ind2_ab) &
                                                                          *  (Coul_baba+Exch_abba)
                                    end if 

                                    RDMEnergy2 = RDMEnergy2 +  0.5_dp * rdm%abab_full(Ind1_ab,Ind2_ab) &
                                                               * Norm_2RDM * (Coul_abab+Exch_baab)


                                    RDMEnergy2 = RDMEnergy2 +  0.5_dp * rdm%baba_full(Ind1_ab,Ind2_ab) &
                                                               * Norm_2RDM * (Coul_baba+Exch_abba)


                                else if (a .eq. b) then
                                    ! a = b : Swap last indeces to get abba/baab terms
                                    ! abba saved in abab, baab saved in baba (sign changes)
                                    Exch_abba = real(UMAT(UMatInd(2*i, 2*j-1, 2*a, 2*b-1)), dp)
                                    Exch_baab = real(UMAT(UMatInd(2*i-1, 2*j, 2*a-1, 2*b)), dp)

                                    if (tRDMInstEnergy) then
                                        RDMEnergy_Inst = RDMEnergy_Inst + 0.5_dp * rdm%abab(Ind1_ab,Ind2_ab) &
                                                                    *(Coul_abab+Exch_abba)
                                        RDMEnergy_Inst = RDMEnergy_Inst + 0.5_dp * rdm%baba(Ind1_ab,Ind2_ab) &
                                                                    * (Coul_baba+Exch_baab)
                                    end if 

                                    RDMEnergy2 = RDMEnergy2 + 0.5_dp * rdm%abab_full(Ind1_ab,Ind2_ab) &
                                                            * Norm_2RDM * (Coul_abab+Exch_abba)


                                    RDMEnergy2 = RDMEnergy2 + 0.5_dp * rdm%baba_full(Ind1_ab,Ind2_ab) &
                                                            * Norm_2RDM * (Coul_baba+Exch_baab)
                                end if

                            else

                                if (tRDMInstEnergy) then 
                                    RDMEnergy_Inst = RDMEnergy_Inst + 0.5_dp * rdm%abab(Ind1_ab,Ind2_ab) &
                                                                    * (Coul+Exch)
                                    if (tOpenShell) RDMEnergy_Inst = RDMEnergy_Inst  &
                                                                    + 0.5_dp *rdm%baba(Ind1_ab,Ind2_ab) &
                                                                    *  (Coul+Exch)
                                end if 

                                RDMEnergy2 = RDMEnergy2 + 0.5_dp * rdm%abab_full(Ind1_ab,Ind2_ab) &
                                                        * Norm_2RDM * (Coul+Exch)


                                if (tOpenShell) RDMEnergy2 = RDMEnergy2 &
                                                        + 0.5_dp * rdm%baba_full(Ind1_ab,Ind2_ab) &
                                                        * Norm_2RDM * (Coul+Exch)

                            end if 

                            if (Ind1_ab .eq. Ind2_ab) then
                                Trace_2RDM_normalised = Trace_2RDM_normalised + &
                                                    rdm%abab_full(Ind1_ab,Ind2_ab) * Norm_2RDM
                                if (tOpenShell) then
                                    Trace_2RDM_normalised = Trace_2RDM_normalised + &
                                                    rdm%baba_full(Ind1_ab,Ind2_ab) * Norm_2RDM
                                end if  
                            end if

                        end if

                   end do

                end do
            end do
        end do

        ! Finally, add in the core energy:

        ! The total energy from the 'instantaneous' RDM.
        if (tRDMInstEnergy) RDMEnergy_Inst = RDMEnergy_Inst + Ecore/Norm_2RDM_Inst
        ! The total energy from the 'full' RDM.
        RDMEnergy = RDMEnergy1 + RDMEnergy2 + Ecore

        call halt_timer(RDMEnergy_Time)

    end subroutine Calc_Energy_from_RDM
    
    function calc_2rdm_spin_estimate(rdm, Norm_2RDM_Inst) result(spin_est)

        ! Return the (unnormalised) estimate of <S^2> from the instantaneous
        ! 2RDM estimates.

        use rdm_data, only: rdm_t, tOpenShell
        use RotateOrbsData, only: SpatOrbs
        use SystemData, only: nel
 
        type(rdm_t), intent(in) :: rdm
        real(dp), intent(in) :: Norm_2RDM_Inst

        integer :: i, j
        integer :: Ind1_aa, Ind1_ab
        real(dp) :: spin_est

        spin_est = 0.0_dp

        do i = 1, SpatOrbs

            do j = i+1, SpatOrbs

                Ind1_aa = ( ( (j-2) * (j-1) ) / 2 ) + i
                Ind1_ab = ( ( (j-1) * j ) / 2 ) + i

                spin_est = spin_est + 2*rdm%aaaa(Ind1_aa, Ind1_aa) &
                                    - 2*rdm%abab(Ind1_ab, Ind1_ab) &
                                    + 4*rdm%abba(Ind1_aa, Ind1_aa)

                if (tOpenShell) then
                    spin_est = spin_est + 2*rdm%bbbb(Ind1_aa, Ind1_aa) &
                                        - 2*rdm%baba(Ind1_ab, Ind1_ab) &
                                        + 4*rdm%baab(Ind1_aa, Ind1_aa)
                end if

            end do

            ! i = j term.
            Ind1_ab = ( ( (i-1) * i ) / 2 ) + i

            spin_est = spin_est - 6*rdm%abab(Ind1_ab, Ind1_ab)
            if (tOpenShell) spin_est = spin_est - 6*rdm%baba(Ind1_ab, Ind1_ab)

        end do 

        spin_est = spin_est + 3*real(nel,dp)/Norm_2RDM_Inst

        spin_est = spin_est/4.0_dp

    end function calc_2rdm_spin_estimate

    subroutine Calc_Lagrangian_from_RDM(rdm, Norm_1RDM, Norm_2RDM)

        ! This routine takes the 1 electron and 2 electron reduced density
        ! matrices and calculated the Lagrangian term, X, required for the
        ! calculation of forces.

        ! The equation for X is as follows:
        !
        !   X_pq = Sum_r[h_pr 1RDM_qr] + 0.5*Sum_rst[(pr|st)[2RDM_qrst + 2RDM_rqst]]
        !
        !   where 2RDM is defined in chemical notation sense:
        !
        !   2RDM_ijkl = <Psi| a_i+ a_k+ a_l a_j|Psi>

        use IntegralsData, only: UMAT
        use OneEInts, only: TMAT2D
        use Parallel_neci, only: iProcIndex
        use rdm_data, only: rdm_t, tOpenShell
        use rdm_temp, only: Find_Spatial_2RDM_Chem
        use RotateOrbsMod, only: SymLabelListInv_rot, SpatOrbs
        use UMatCache, only: UMatInd

        type(rdm_t), intent(inout) :: rdm
        real(dp), intent(in) :: Norm_2RDM
        real(dp), intent(in) :: Norm_1RDM

        integer :: p, q, r, s, t, ierr
        integer :: pSpin, rSpin
        real(dp) :: Coul
        real(dp) :: qrst, rqst
        real(dp) :: Max_Error_Hermiticity, Sum_Error_Hermiticity, Temp

        allocate(rdm%Lagrangian(SpatOrbs,SpatOrbs), stat=ierr)
        rdm%Lagrangian(:,:) = 0.0_dp
        
        ! We will begin by calculating the Lagrangian in chemical notation - we will explicitely calculate
        ! both halves (X_pq and X_qp) in order to check if it is symmetric or not.
        !  - a symmetric Lagrangian is important in allowing us to use the derivative overlap matrix rather than the 
        !    coupled-perturbed coefficients when calculating the Forces later on 
        !    (see Sherrill, Analytic Gradients of CI Energies eq38)

        write(6,'(/,"Calculating the Lagrangian, X, from the final density matrices.")')

        ! Calculating the Lagrangian X in terms of spatial orbitals.
        if (iProcIndex .eq. 0) then

            do p = 1, SpatOrbs  !Run over spatial orbitals
                pSpin = 2*p    ! Picks out beta component
                do q = 1, SpatOrbs  !Want to calculate X(p,q) separately from X(q,p) for now to see if we're symmetric
                    do r = 1, SpatOrbs
                        rSpin = 2*r
                        ! Adding in contributions effectively from the 1-RDM.
                        ! We made sure earlier that the 1RDM is contructed, so we can call directly from this.
                        if (tOpenShell) then
                            ! Include both aa and bb contributions 
                            rdm%Lagrangian(p,q) = rdm%Lagrangian(p,q) + &
                                                      (rdm%matrix(SymLabelListInv_rot(2*q),SymLabelListInv_rot(2*r)))*&
                                                      real(TMAT2D(pSpin,rSpin),8)*Norm_1RDM
                            rdm%Lagrangian(p,q) = rdm%Lagrangian(p,q) + &
                                                      (rdm%matrix(SymLabelListInv_rot(2*q-1),SymLabelListInv_rot(2*r-1)))*&
                                                      real(TMAT2D(pSpin-1,rSpin-1),8)*Norm_1RDM
                        else
                            ! We will be here most often (?)
                            rdm%Lagrangian(p,q) = rdm%Lagrangian(p,q) + rdm%matrix(SymLabelListInv_rot(q),SymLabelListInv_rot(r))* &
                                                                                  real(TMAT2D(pSpin,rSpin),8)*Norm_1RDM
                        end if

                        do s = 1, SpatOrbs
                
                            ! Here we're looking to start adding on the contributions from the 2-RDM and the 2-el integrals
                            ! For X(p,q), these have the form 0.5*Sum_rst[(pr|st)[2RDM_qrst + 2RDM_rqst]]
                            
                            ! NOTE: In some notations, the factor of a half goes *inside* the 2RDM - 
                            ! consistent with Yamaguchi etc (see eq 11 in Sherrill, Analytic Gradients
                            ! of CI energies).  We will keep it outside for now, consistent with the storage
                            ! within neci

                            do t = 1, SpatOrbs
                                
                                !Integral (pr|st) = <ps|rt>
                                !Give indices in PHYSICAL NOTATION
                                !NB, FCIDUMP is labelled in chemical notation
                                Coul = real(UMAT(UMatInd(p, s, r, t)),dp)

                                qrst = Find_Spatial_2RDM_Chem(rdm, q, r, s, t, Norm_2RDM)
                                rqst = Find_Spatial_2RDM_Chem(rdm, r, q, s, t, Norm_2RDM)
                                
                                rdm%Lagrangian(p,q) = rdm%Lagrangian(p,q) + 0.5_dp*Coul*(qrst+rqst)
                            end do
                        end do
                    end do
                end do
            end do
       
            !! Now symmetrise (make hermitian, such that X_pq = X_qp) the Lagrangian X

            Max_Error_Hermiticity = 0.0_dp
            Sum_Error_Hermiticity = 0.0_dp
            do p = 1, SpatOrbs
                do q = p, SpatOrbs
                    if (abs(rdm%Lagrangian(p,q) - rdm%Lagrangian(q,p)) .gt. Max_Error_Hermiticity) &
                        Max_Error_Hermiticity = abs(rdm%Lagrangian(p,q) - rdm%Lagrangian(q,p))

                    Sum_Error_Hermiticity = Sum_Error_Hermiticity+abs(rdm%Lagrangian(p,q) - rdm%Lagrangian(q,p))

                    Temp = (rdm%Lagrangian(p,q) + rdm%Lagrangian(q,p))/2.0_dp
                    rdm%Lagrangian(p,q) = Temp
                    rdm%Lagrangian(q,p) = Temp
                end do
            end do

            ! Output the hermiticity errors.
            write(6,'(1X,"MAX ABS ERROR IN Lagrangian HERMITICITY",F30.20)') Max_Error_Hermiticity
            write(6,'(1X,"SUM ABS ERROR IN Lagrangian HERMITICITY",F30.20)') Sum_Error_Hermiticity

        end if

    end subroutine Calc_Lagrangian_from_RDM

    subroutine calc_1RDM_and_1RDM_energy(rdm, i, j, a, iSpin, jSpin, Norm_2RDM, &
                                         RDMEnergy_Inst, RDMEnergy1)

        ! This routine calculates the 1-RDM part of the RDM energy, and
        ! also constructs the 1-RDM if required, both from the 2-RDM.

        ! gamma(i,j) = [1/(NEl - 1)] * SUM_a Gamma(i,a,j,a) 
        ! want to calculate:    gamma(i,j) * h_ij
        ! h_ij => TMAT2D(iSpin,jSpin)
        ! iSpin = 2*i, jSpin = 2*j -> alpha orbs

        use FciMCData, only: tFinalRDMEnergy
        use OneEInts, only: TMAT2D
        use LoggingData, only: tDiagRDM, tDumpForcesInfo, tDipoles, tRDMInstEnergy, tPrint1RDM
        use rdm_data, only: rdm_t, tOpenShell
        use RotateOrbsMod, only: SymLabelListInv_rot
        use SystemData, only: nel

        type(rdm_t), intent(inout) :: rdm
        integer, intent(in) :: i, j, a, iSpin, jSpin
        real(dp), intent(in) :: Norm_2RDM
        real(dp), intent(inout) :: RDMEnergy_Inst, RDMEnergy1
        real(dp) :: Parity_Factor, fac_doublecount

        integer :: Ind1_1e_ab, Ind2_1e_ab
        integer :: Ind1_1e_aa, Ind2_1e_aa
        integer :: iSpin_abab, iSpin_baba
        integer :: jSpin_abab, jSpin_baba
        integer :: iSpin_abba, iSpin_baab
        integer :: jSpin_abba, jSpin_baab
        logical :: t_abab_only, t_opposite_contri

        ! for i a -> j a excitation, when lined up as min max -> min max, 
        ! if a's are aligned, only a b a b arrays contain single excitations, 
        ! if a's not aligned, a b b a.
        ! all a a a a will contain single excitations.

        ! abab & baba terms
        if (((i .le. a) .and. (j .le. a)) .or. ((i .ge. a) .and. (j .ge. a))) then

            Ind1_1e_ab = ( ( (max(i,a)-1) * max(i,a) ) / 2 ) + min(i,a)
            Ind2_1e_ab = ( ( (max(j,a)-1) * max(j,a) ) / 2 ) + min(j,a)
            if (Ind1_1e_ab.ne.Ind2_1e_ab) then
            ! For Gamma elements corresponding to 1-RDMs ( Gamma(i,a,j,a) ), 
            ! we're only considering i =< j 
            ! therefore we need to sum in the opposite contribution too if i ne j.
                t_opposite_contri = .true.
            else  ! no opposite contribution from i = j term
                t_opposite_contri = .false.
            end if

            fac_doublecount = 1.0_dp
            if ( (i .lt. a) .or. (j .lt. a) )then
                ! i a j a  ->  i & j alpha for abab and beta for baba
                iSpin_abab = iSpin
                jSpin_abab = jSpin
                iSpin_baba = iSpin-1
                jSpin_baba = jSpin-1
                t_abab_only = .false. 
            else if ((i .gt. a) .or. (j .gt. a) )then
                ! a i a j->  i & j alpha for baba and beta for abab
                iSpin_abab = iSpin-1
                jSpin_abab = jSpin-1
                iSpin_baba = iSpin
                jSpin_baba = jSpin
                t_abab_only = .false. 
            else if ((i .eq. a) .and. (j .eq. a) )then
                ! a a a a -> i = a = j  abab and baba saved in abab array only! 
                ! -> count twice for close shell systems (fac_doublecount)
                iSpin_abab = iSpin
                jSpin_abab = jSpin
                iSpin_baba = iSpin-1
                jSpin_baba = jSpin-1
                t_abab_only = .true. 
                if (.not. tOpenShell) fac_doublecount = 2.0_dp
            end if

            if (tRDMInstEnergy) then
                RDMEnergy_Inst = RDMEnergy_Inst + &
                           fac_doublecount*( (rdm%abab(Ind1_1e_ab,Ind2_1e_ab) ) &
                                               * real(TMAT2D(iSpin_abab,jSpin_abab),dp) &
                                               * (1.0_dp / real(NEl - 1,dp)) )

                if (t_opposite_contri) RDMEnergy_Inst = RDMEnergy_Inst + &
                           ( (rdm%abab(Ind2_1e_ab, Ind1_1e_ab) ) &
                                               * real(TMAT2D(jSpin_abab,iSpin_abab),dp) &
                                               * (1.0_dp / real(NEl - 1,dp)) )
                if (tOpenShell) then ! add baba terms
                    if (.not. t_abab_only) then
                        RDMEnergy_Inst = RDMEnergy_Inst + &
                                       ( (rdm%baba(Ind1_1e_ab,Ind2_1e_ab) ) &
                                       * real(TMAT2D(iSpin_baba,jSpin_baba),dp) &
                                       * (1.0_dp / real(NEl - 1,dp)) )

                        if (t_opposite_contri)  RDMEnergy_Inst = RDMEnergy_Inst + &
                                       ( (rdm%baba(Ind2_1e_ab,Ind1_1e_ab) ) &
                                       * real(TMAT2D(jSpin_baba,iSpin_baba),dp) &
                                       * (1.0_dp / real(NEl - 1,dp)) )
                    else   ! i=j=a -> rdm%baba saved in rdm%abab & t_opposite_contri = false
                        RDMEnergy_Inst = RDMEnergy_Inst + &
                                       ( (rdm%abab(Ind1_1e_ab,Ind2_1e_ab) ) &
                                       * real(TMAT2D(iSpin_baba,jSpin_baba),dp) &
                                       * (1.0_dp / real(NEl - 1,dp)) )
                    end if
                end if
            end if 

            RDMEnergy1 = RDMEnergy1 + fac_doublecount*( (rdm%abab_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM) &
                                               * real(TMAT2D(iSpin_abab,jSpin_abab),dp) &
                                               * (1.0_dp / real(NEl - 1,dp)) )

            if (t_opposite_contri) RDMEnergy1 = RDMEnergy1 + ( (rdm%abab_full(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM) &
                                               * real(TMAT2D(jSpin_abab,iSpin_abab),dp) &
                                               * (1.0_dp / real(NEl - 1,dp)) )

            if (tOpenShell) then  ! add baba terms
                if (.not. t_abab_only) then
                    RDMEnergy1 = RDMEnergy1 + ( (rdm%baba_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM) &
                                 * real(TMAT2D(iSpin_baba,jSpin_baba),dp) &
                                 * (1.0_dp / real(NEl - 1,dp)) )

                     if (t_opposite_contri) RDMEnergy1 = RDMEnergy1 + ( (rdm%baba_full(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM) &
                                 * real(TMAT2D(jSpin_baba,iSpin_baba),dp) &
                                 * (1.0_dp / real(NEl - 1,dp)) )
                else ! rdm%baba saved in rdm%abab & t_opposite_contri = false
                    RDMEnergy1 = RDMEnergy1 + ( (rdm%abab_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM) &
                                 * real(TMAT2D(iSpin_baba,jSpin_baba),dp) &
                                 * (1.0_dp / real(NEl - 1,dp)) )
                end if

            end if


            if ((tDiagRDM .or. tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) .and. tFinalRDMEnergy) then
               
                if (.not. tOpenShell) then
                 ! Spatial orbitals
                    rdm%matrix(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = &
                             rdm%matrix(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) &
                                         + fac_doublecount*( rdm%abab_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM &
                                                 * (1.0_dp / real(NEl - 1,dp)) )

                    if (t_opposite_contri) rdm%matrix(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) = &
                             rdm%matrix(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) &
                                         + ( rdm%abab_full(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM &
                                                 * (1.0_dp / real(NEl - 1,dp)) )

                else                      
                   rdm%matrix(SymLabelListInv_rot(iSpin_abab),SymLabelListInv_rot(jSpin_abab)) = &
                           rdm%matrix(SymLabelListInv_rot(iSpin_abab),SymLabelListInv_rot(jSpin_abab)) &
                                       + ( rdm%abab_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM &
                                       * (1.0_dp / real(NEl - 1,dp)) )

                   if (t_opposite_contri)  rdm%matrix(SymLabelListInv_rot(jSpin_abab),SymLabelListInv_rot(iSpin_abab)) = &
                           rdm%matrix(SymLabelListInv_rot(jSpin_abab),SymLabelListInv_rot(iSpin_abab)) &
                                       + ( rdm%abab_full(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM &
                                       * (1.0_dp / real(NEl - 1,dp)) ) 

                   if (.not. t_abab_only) then
                       rdm%matrix(SymLabelListInv_rot(iSpin_baba),SymLabelListInv_rot(jSpin_baba)) = &
                            rdm%matrix(SymLabelListInv_rot(iSpin_baba),SymLabelListInv_rot(jSpin_baba)) &
                                      + ( rdm%baba_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM &
                                      * (1.0_dp / real(NEl - 1,dp)) )

                       if (t_opposite_contri) rdm%matrix(SymLabelListInv_rot(jSpin_baba),SymLabelListInv_rot(iSpin_baba)) = &
                            rdm%matrix(SymLabelListInv_rot(jSpin_baba),SymLabelListInv_rot(iSpin_baba)) &
                                      + ( rdm%baba_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM &
                                      * (1.0_dp / real(NEl - 1,dp)) )
                   else ! i = j = a -> baba saved in abab & t_opposite_contri = false
                       rdm%matrix(SymLabelListInv_rot(iSpin_baba),SymLabelListInv_rot(jSpin_baba)) = &
                          rdm%matrix(SymLabelListInv_rot(iSpin_baba),SymLabelListInv_rot(jSpin_baba)) &
                                      + ( rdm%abab_full(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM &
                                      * (1.0_dp / real(NEl - 1,dp)) )
                   end if
                end if

           end if
 
       end if ! abab & baba terms

       ! abba & baab & aaaa & bbbb terms
       if ((i.ne.a).and.(j.ne.a)) then   
           Ind1_1e_aa = ( ( (max(i,a)-2) * (max(i,a)-1) ) / 2 ) + min(i,a)
           Ind2_1e_aa = ( ( (max(j,a)-2) * (max(j,a)-1) ) / 2 ) + min(j,a)

           if (Ind1_1e_aa.ne.Ind2_1e_aa) then
           ! For Gamma elements corresponding to 1-RDMs (eg Gamma(i,a,a,j) ), 
           ! we're only considering i =< j 
           ! therefore we need to sum in the opposite contribution too.
               t_opposite_contri = .true.
           else  ! no opposite contribution from i = j term
               t_opposite_contri = .false.
           end if

           !abba & baab terms
           if ((i.ne.j).and.((i.lt.a).and.(j.gt.a)).or.((i.gt.a).and.(j.lt.a))) then 
               if ((i.lt.a).and.(j.gt.a))then
                   ! i a a j -> i & j alpha for abba and beta for baab
                   iSpin_abba = iSpin
                   jSpin_abba = jSpin
                   iSpin_baab = iSpin-1
                   jSpin_baab = jSpin-1
               else if ((i.gt.a).and.(j.lt.a))then
                   ! a i j a  -> i & j beta for abba and alpha for baab
                   iSpin_abba = iSpin-1
                   jSpin_abba = jSpin-1
                   iSpin_baab = iSpin
                   jSpin_baab = jSpin
               end if

               if (tRDMInstEnergy) then
                   RDMEnergy_Inst = RDMEnergy_Inst - &
                                                  ( (rdm%abba(Ind1_1e_aa,Ind2_1e_aa) ) &
                                                  * real(TMAT2D(iSpin_abba,jSpin_abba),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

                   if (t_opposite_contri .and. (.not. tOpenShell) ) then
                       RDMEnergy_Inst = RDMEnergy_Inst - &
                                                 ( (rdm%abba(Ind2_1e_aa,Ind1_1e_aa) ) &
                                                  * real(TMAT2D(jSpin_abba,iSpin_abba),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )
                   end if

                   if (tOpenShell) then
                       ! abba becomes baab when i and j are swapped for the opposite contribution
                       if (t_opposite_contri ) RDMEnergy_Inst = RDMEnergy_Inst - &
                                                 ( (rdm%baab(Ind2_1e_aa,Ind1_1e_aa) ) &
                                                  * real(TMAT2D(jSpin_abba,iSpin_abba),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

                       RDMEnergy_Inst = RDMEnergy_Inst - &
                                                  ( (rdm%baab(Ind1_1e_aa,Ind2_1e_aa) ) &
                                                  * real(TMAT2D(iSpin_baab,jSpin_baab),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

                       if (t_opposite_contri) RDMEnergy_Inst = RDMEnergy_Inst - &
                                                 ( (rdm%abba(Ind2_1e_aa,Ind1_1e_aa) ) &
                                                  * real(TMAT2D(jSpin_baab,iSpin_baab),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )
                   end if

               end if 

               RDMEnergy1 = RDMEnergy1 - ( (rdm%abba_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM) &
                                                  * real(TMAT2D(iSpin_abba,jSpin_abba),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

               if (t_opposite_contri .and. (.not. tOpenShell) ) then
                   RDMEnergy1 = RDMEnergy1 - ( (rdm%abba_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                  * real(TMAT2D(jSpin_abba,iSpin_abba),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )
               end if

               if (tOpenShell) then  ! add baab terms.
                   ! abba becomes baab when i and j are swapped for the opposite contribution.
                   if (t_opposite_contri) RDMEnergy1 = RDMEnergy1 - ( (rdm%baab_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                  * real(TMAT2D(jSpin_abba,iSpin_abba),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

                   RDMEnergy1 = RDMEnergy1 - ( (rdm%baab_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM) &
                                                  * real(TMAT2D(iSpin_baab,jSpin_baab),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

                   if (t_opposite_contri) RDMEnergy1 = RDMEnergy1 - ( (rdm%abba_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                  * real(TMAT2D(jSpin_baab,iSpin_baab),dp) &
                                                  * (1.0_dp / real(NEl - 1,dp)) )
               end if
                                                  
               if ((tDiagRDM .or. tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) .and. tFinalRDMEnergy) then
                   if (.not. tOpenShell) then
                   ! Spatial orbitals
                       rdm%matrix(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = &
                           rdm%matrix(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) &
                                      - ( rdm%abba_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                              * (1.0_dp / real(NEl - 1,dp)) )

                       if (t_opposite_contri) rdm%matrix(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) = &
                          rdm%matrix(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) &
                                      - ( rdm%abba_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                              * (1.0_dp / real(NEl - 1,dp)) )
                   else
                       rdm%matrix(SymLabelListInv_rot(iSpin_abba),SymLabelListInv_rot(jSpin_abba)) = &
                              rdm%matrix(SymLabelListInv_rot(iSpin_abba),SymLabelListInv_rot(jSpin_abba)) &
                                          - ( rdm%abba_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                  * (1.0_dp / real(NEl - 1,dp)) ) 

                       if (t_opposite_contri) rdm%matrix(SymLabelListInv_rot(jSpin_abba),SymLabelListInv_rot(iSpin_abba)) = &
                              rdm%matrix(SymLabelListInv_rot(jSpin_abba),SymLabelListInv_rot(iSpin_abba)) &
                                          - ( rdm%baab_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                                  * (1.0_dp / real(NEl - 1,dp)) ) 

                       rdm%matrix(SymLabelListInv_rot(iSpin_baab),SymLabelListInv_rot(jSpin_baab)) = &
                              rdm%matrix(SymLabelListInv_rot(iSpin_baab),SymLabelListInv_rot(jSpin_baab)) &
                                          - ( rdm%baab_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                  * (1.0_dp / real(NEl - 1,dp)) )

                       if (t_opposite_contri) rdm%matrix(SymLabelListInv_rot(jSpin_baab),SymLabelListInv_rot(iSpin_baab)) = &
                              rdm%matrix(SymLabelListInv_rot(jSpin_baab),SymLabelListInv_rot(iSpin_baab)) &
                                          - ( rdm%abba_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                                  * (1.0_dp / real(NEl - 1,dp)) )
                   end if
               end if                   

           end if ! end of abba and baab terms

           ! aaaa & bbbb terms
           if (((i .lt. a) .and. (j .lt. a)) .or. ((i .gt. a) .and. (j .gt. a))) then
               Parity_Factor = 1.0_dp
           else
               Parity_Factor = -1.0_dp
           end if
          
           if (tRDMInstEnergy) then
               RDMEnergy_Inst = RDMEnergy_Inst + &
                          ( (rdm%aaaa(Ind1_1e_aa,Ind2_1e_aa) ) &
                                                * real(TMAT2D(iSpin,jSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )
               if (t_opposite_contri) RDMEnergy_Inst = RDMEnergy_Inst + &
                          ( (rdm%aaaa(Ind2_1e_aa,Ind1_1e_aa) ) &
                                                * real(TMAT2D(jSpin,iSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )
               if (tOpenShell) then
                   RDMEnergy_Inst = RDMEnergy_Inst +&
                          ( (rdm%bbbb(Ind1_1e_aa,Ind2_1e_aa)) &
                                                * real(TMAT2D(iSpin-1,jSpin-1),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor)

                   if (t_opposite_contri) RDMEnergy_Inst = RDMEnergy_Inst +&
                          ( (rdm%bbbb(Ind2_1e_aa,Ind1_1e_aa)) &
                                                * real(TMAT2D(jSpin-1,iSpin-1),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor)
               end if

           end if 

           RDMEnergy1 = RDMEnergy1 + &
                           ( (rdm%aaaa_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM) &
                                                * real(TMAT2D(iSpin,jSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )

           if (t_opposite_contri) RDMEnergy1 = RDMEnergy1 + &
                           ( (rdm%aaaa_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                * real(TMAT2D(jSpin,iSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )

           if (tOpenShell) then
               RDMEnergy1 = RDMEnergy1 + &
                           ( (rdm%bbbb_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM) &
                                                * real(TMAT2D(iSpin-1,jSpin-1),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )
               if (t_opposite_contri) RDMEnergy1 = RDMEnergy1 + &
                           ( (rdm%bbbb_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                * real(TMAT2D(jSpin-1,iSpin-1),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )
           end if

           if ((tDiagRDM .or. tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) .and. tFinalRDMEnergy) then
               if ( .not. tOpenShell) then
                   rdm%matrix(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = &
                            rdm%matrix(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) &
                                        + ( rdm%aaaa_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 

                   if (t_opposite_contri) rdm%matrix(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) = &
                            rdm%matrix(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) &
                                        + ( rdm%aaaa_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 

               else

                   rdm%matrix(SymLabelListInv_rot(iSpin),SymLabelListInv_rot(jSpin)) = &
                            rdm%matrix(SymLabelListInv_rot(iSpin),SymLabelListInv_rot(jSpin)) &
                                        + ( rdm%aaaa_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 

                   if (t_opposite_contri) rdm%matrix(SymLabelListInv_rot(jSpin),SymLabelListInv_rot(iSpin)) = &
                            rdm%matrix(SymLabelListInv_rot(jSpin),SymLabelListInv_rot(iSpin)) &
                                        + ( rdm%aaaa_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 

                   rdm%matrix(SymLabelListInv_rot(iSpin-1),SymLabelListInv_rot(jSpin-1)) = & 
                            rdm%matrix(SymLabelListInv_rot(iSpin-1),SymLabelListInv_rot(jSpin-1)) &
                                        + ( rdm%bbbb_full(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 

                   if (t_opposite_contri) rdm%matrix(SymLabelListInv_rot(jSpin-1),SymLabelListInv_rot(iSpin-1)) = & 
                            rdm%matrix(SymLabelListInv_rot(jSpin-1),SymLabelListInv_rot(iSpin-1)) &
                                        + ( rdm%bbbb_full(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 
               end if
           end if

       end if 

    end subroutine calc_1RDM_and_1RDM_energy

    subroutine convert_mats_Molpforces(rdm, Norm_1RDM, Norm_2RDM)

        use GenRandSymExcitNUMod, only: ClassCountInd, RandExcitSymLabelProd
        use Parallel_neci, only: iProcIndex
        use rdm_data, only: rdm_t, tOpenShell
        use rdm_temp, only: Find_Spatial_2RDM_Chem
        use RotateOrbsData, only: SpatOrbs, SymLabelListInv_rot
        use sym_mod
        use SymData, only: Sym_Psi, nSymLabels
        use SymExcitDataMod, only: SpinOrbSymLabel,SymLabelCounts2
        use SystemData, only: nEl, LMS

        type(rdm_t), intent(in) :: rdm

        integer :: iblkq, iseccr, istat1, isyref, ms2
        integer :: posn1, posn2
        integer :: i, j, k, l
        integer :: myname, ifil, intrel, iout
        integer :: Sym_i, Sym_j, Sym_ij
        integer :: Sym_k, Sym_l, Sym_kl
        integer, dimension(8) :: iact, ldact !iact(:) # of active orbs per sym, 
                                             !ldact(:) - # Pairs of orbs that multiply to give given sym
        integer, dimension(8) :: icore, iclos 
        integer, dimension(nSymLabels):: blockstart1, blockstart2
        integer, dimension(nSymLabels) :: elements_assigned1, elements_assigned2
        integer :: FC_Lag_Len  ! Length of the Frozen Core Lagrangian
        integer :: Len_1RDM, Len_2RDM, FCLag_Len
        real(dp) :: ijkl, jikl
        real(dp), intent(in) :: Norm_1RDM, Norm_2RDM
        real(dp), allocatable :: SymmetryPacked2RDM(:), SymmetryPacked1RDM(:) !SymPacked2RDM in CHEMICAL Notation
        real(dp), allocatable :: SymmetryPackedLagrangian(:)
        real(dp), allocatable :: FC_Lagrangian(:)
        integer :: iWfRecord, iWfSym
#ifdef MOLPRO
        character(len=*), parameter :: t_r='convert_mats_Molpforces'
#endif

#ifdef __int64
        intrel = 1
#else
        intrel = 2
#endif

        if (iProcIndex .eq. 0) then
            ! Calculating header information needed for the molpro dump routine.
            ! Considering all electron for now.
            icore(:) = 0  ! Number of frozen orbitals per symmetry
            iclos(:) = 0  ! icore(i) + Number of 'closed' orbitals per symmetry (iclos is the same as icore for us).
            call molpro_get_reference_info(iWfRecord, iWfSym)
            iblkq = iWfRecord ! record number for orbitals
            iseccr = 0          ! record number for core orbitals
            istat1 = 1        !ground state
            isyref = Sym_Psi + 1  !spatial symmetry of the wavefunction
#ifdef MOLPRO
            if (isyref .ne. iWfSym) call stop_all(t_r,"NECI and common/cref do not agree on irrep of wave function")
#endif
            ms2 = LMS  !2 * M_s
            myname = 5001 !Arbitrary file names
            ifil = 1
            ! ^- cgk: might want to use nexfre(igrsav) for these two.
            iout = molpro_get_iout()
            ! ^- thise one should come from common/tapes. 
            ldact(:) = 0
            iact(:) = 0
            Len_1RDM = 0
            Len_2RDM = 0
            blockstart1(:) = 0
            blockstart2(:) = 0
            elements_assigned1(:) = 0
            elements_assigned2(:) = 0
            
            ! Find out the number of orbital pairs that multiply to a given sym (ldact).
            do i = 1, SpatOrbs  
                ! Run over spatial orbitals.
                do j = 1,  i    ! i .ge. j
                    Sym_i=SpinOrbSymLabel(2*i)  ! Consider only alpha orbitals.
                    Sym_j=SpinOrbSymLabel(2*j)
                    Sym_ij=RandExcitSymLabelProd(Sym_i, Sym_j)
                    ldact(Sym_ij+1)=ldact(Sym_ij+1)+1
                end do
            end do

            ! Calculate lengths of arrays, and where each sym block starts.
            do i = 0, nSymLabels-1

                ! CMO: Check if Sym_i goes 0-->7 or 1-->8.
                ! Find position of each symmetry block in sym-packed forms of RDMS 1 & 2.
                blockstart1(i+1)=Len_1RDM+1 ! N.B. Len_1RDM still being updated in this loop.
                blockstart2(i+1)=Len_2RDM+1 ! N.B. Len_2RDM still being updated in this loop.
                
                ! Count the number of active orbitals of the given symmetry.
                iact(i+1)=SymLabelCounts2(2,ClassCountInd(1,i,0)) !Count the number of active orbitals of the given symmetry.
                Len_1RDM=Len_1RDM+(iact(i+1)*(iact(i+1)+1)/2) ! add on # entries in sym-packed 1RDM for sym i.
                
                Len_2RDM = Len_2RDM + (ldact(i+1))**2 ! Assumes no frozen orbitals.
            end do

            FCLag_Len = SpatOrbs**2  ! Arbitrarily set this for now - we will not be printing it whilst nfrozen=0
            
            ! Allocate arrays accordingly.
            allocate(SymmetryPacked1RDM(Len_1RDM))
            allocate(SymmetryPackedLagrangian(Len_1RDM))
            allocate(SymmetryPacked2RDM(Len_2RDM))
            allocate(FC_Lagrangian(FCLag_Len))

            FC_Lagrangian(:) = 0  ! Frozen-core Lagrandian -- whilst we do all electron calcs.
            
            ! Constructing the Symmetry Packed arrays.
            ! We convert our 1RDM, Lagrangian and  2RDM into the required Molpro
            ! symmetry-packed format 2RDM is stored as a full square matrix,
            ! separated into symmetry blocks. We store D_ijkl (chemical notation,
            ! spatial orbs) where (i .ge. j) and (k .ge. l). For each symmetry X,
            ! there will be a block where (i,j) and (k,l) both have symmetry X
            ! making (ij,kl) totally symmetric, and D_ijkl (potentially) non-zero.
            ! 1RDM and Lagrangian are stored as upper triangles, separated by
            ! symmetry block

            SymmetryPacked2RDM(:) = 0.0_dp
            SymmetryPacked1RDM(:) = 0.0_dp
            SymmetryPackedLagrangian(:) = 0.0_dp

            do i = 1, SpatOrbs  !run over spatial orbitals, ALL ELECTRON ONLY
                do j = 1,  i    ! i .ge. j
                    Sym_i = SpinOrbSymLabel(2*i)  !Consider only alpha orbitals
                    Sym_j = SpinOrbSymLabel(2*j)
                    Sym_ij = RandExcitSymLabelProd(Sym_i, Sym_j)
                    if (Sym_ij .eq. 0) then
                        posn1 = blockstart1(Sym_i+1) + elements_assigned1(Sym_i+1)
                        ! Add pre-symmetrised contribution to the symmetry-packed 1-RDM.

                        if (tOpenShell) then
                            ! Include both aa and bb contributions.
                            SymmetryPacked1RDM(posn1)=&
                                  (rdm%matrix(SymLabelListInv_rot(2*i),SymLabelListInv_rot(2*j))&
                                 + rdm%matrix(SymLabelListInv_rot(2*i-1),SymLabelListInv_rot(2*j-1)))*Norm_1RDM
                        else
                            SymmetryPacked1RDM(posn1) = rdm%matrix(SymLabelListInv_rot(i),SymLabelListInv_rot(j))*Norm_1RDM
                        end if

                        ! Add in the symmetrised Lagrangian contribution to the
                        ! sym-packed Lagrangian
                        SymmetryPackedLagrangian(posn1) = rdm%Lagrangian(i,j)
                        elements_assigned1(Sym_i+1) = elements_assigned1(Sym_i+1) + 1
                    end if

                    do k = 1, SpatOrbs
                        do l = 1, k
                            Sym_k = SpinOrbSymLabel(2*k)  ! Consider only alpha orbitals
                            Sym_l = SpinOrbSymLabel(2*l)
                            Sym_kl = RandExcitSymLabelProd(Sym_k, Sym_l)

                            if (Sym_kl .eq. Sym_ij) then
                                posn2=blockstart2(Sym_ij+1) + elements_assigned2(Sym_ij+1)
                            
                                ! Extracts spat orb chemical notation term from spin
                                ! separated physical notation RDMs 
                                ijkl = Find_Spatial_2RDM_Chem(rdm, i, j, k, l, Norm_2RDM)
                                jikl = Find_Spatial_2RDM_Chem(rdm, j, i, k, l, Norm_2RDM)

                                SymmetryPacked2RDM(posn2) = 0.5*(ijkl+jikl)
                                elements_assigned2(Sym_ij+1) = elements_assigned2(Sym_ij+1) + 1
                            end if
                        end do
                    end do
                end do
            end do
            
            call molpro_dump_mcscf_dens_for_grad(myname,ifil, &
                icore, iclos, iact, nEL, isyref, ms2, iblkq,&
                iseccr, istat1, SymmetryPacked1RDM, Len_1RDM, &
                SymmetryPacked2RDM, Len_2RDM, SymmetryPackedLagrangian, &
                Len_1RDM, FC_Lagrangian, FC_Lag_Len, &
                iout, intrel)
            
        end if
        
    end subroutine convert_mats_Molpforces

    subroutine molpro_set_igrdsav(igrsav_)

        integer, intent(in) :: igrsav_

#ifdef MOLPRO
        include "common/cwsave"
        ! tell gradient programs where gradient information is stored
        igrsav = igrsav_
#else
        character(len=*), parameter :: t_r = 'molpro_set_igrdsav'
        integer :: iunused
        call stop_all(t_r,'Should be being called from within MOLPRO')
        iunused = igrsav_
#endif

    end subroutine molpro_set_igrdsav
    
    subroutine molpro_get_reference_info(iWfRecord, iWfSym)

        integer :: iWfSym,iWfRecord
#ifdef MOLPRO
        include "common/code"
        include "common/cref"
        ! wf(1) should contain the record number, as integer, of the last used orbital set.
        iWfRecord = wf(1)
        iWfSym = isyref
#else
        integer :: iunused
        iWfRecord = 21402

        ! ELiminate warnings
        iunused = iWfSym
#endif

    end subroutine molpro_get_reference_info
    
    function molpro_get_iout() result(iiout)

        integer :: iiout
#ifdef MOLPRO
        include "common/tapes"
        ! output i/o unit for main output (usually 6, but might be different
        ! if in logfile mode, e.g., during geometry optimization)
        iiout = iout
#else
        iiout = 6
#endif

    end function molpro_get_iout

    subroutine molpro_dump_mcscf_dens_for_grad(name,ifil, &
        icore,iclos,iact,nelec,isyref,ms2,iblkq,iseccr,istat1, &
        den1,lden1, &
        den2,lden2, &
        eps,leps, &
        epsc,lepsc,iout,intrel)

        !  Use the following subroutine to dump the information needed for molpro
        !  forces calculations, such that it is in the exact format that molpro
        !  would have dumped from a MCSCF calculation.  When interfaced with Molpro,
        !  I understand that the molpro equivalent of this routine should be called. 
    
        integer,                            intent(inout) :: name
        integer,                            intent(inout) :: ifil
        integer, dimension(8),              intent(in)    :: icore
        integer, dimension(8),              intent(in)    :: iclos
        integer, dimension(8),              intent(in)    :: iact
        integer,                            intent(in)    :: nelec
        integer,                            intent(in)    :: isyref
        integer,                            intent(in)    :: ms2
        integer,                            intent(in)    :: iblkq
        integer,                            intent(in)    :: iseccr
        integer,                            intent(in)    :: istat1
        integer,                            intent(in)    :: lden1
        real(dp), dimension(lden1), intent(in)    :: den1
        integer,                            intent(in)    :: lden2
        real(dp), dimension(lden2), intent(in)    :: den2
        integer,                            intent(in)    :: leps
        real(dp), dimension(leps),  intent(in)    :: eps
        integer,                            intent(in)    :: lepsc
        real(dp), dimension(lepsc), intent(in)    :: epsc
        integer,                            intent(in)    :: iout
        integer,                            intent(in)    :: intrel
        integer                                           :: igrsav

        !> Write mcscf density in format needed for gradient program.
        !> \param[in,out] name record number to be written
        !> \param[in,out] ifil file number to be written
        !> \param[in] icore numbers of frozen core orbitals in each symmetry
        !> \param[in] iclos plus numbers of closed-shell orbitals in each symmetry
        !> \param[in] iact numbers of active orbitals in each symmetry
        !> \param[in] nelec number of active electrons -- DONE
        !> \param[in] isyref spatial symmetry of wavefunction -- DONE
        !> \param[in] ms2 spin quauntum number times 2 -- DONE
        !> \param[in] iblkq record number * 10 + file number for orbitals (typically 21402) -- DONE
        !> \param[in] iseccr record number * 10 + file number for frozen orbitals (typically 21002) -- DONE
        !> \param[in] istat1 state number (1 for ground state) -- DONE
        !> \param[in] den1 1-particle density matrix
        !> \param[in] lden1 size of den1
        !> \param[in] den2 2-particle density matrix
        !> \param[in] lden2 size of den2
        !> \param[in] eps Lagrangian
        !> \param[in] leps size of leps
        !> \param[in] epsc Lagrangian for frozen core
        !> \param[in] lepsc size of epsc
        !> \param[in] iout unit for output, eg. 6
        !> \param[in] intrel byte size ratio of integer to double precision, normally 1 or 2
        !> \param[out] igrsav where the record was written (needed in common/cwsave)

        integer, dimension(30) :: header
        integer :: i, ncore
        integer :: lhead

        real(dp) :: runused

#ifdef MOLPRO
        character(len=6), parameter :: label='MCGRAD'
#else
        character(*), parameter :: t_r = 'molpro_dump_mcscf_dens_for_grad'
        call stop_all(t_r,'Should not be here if not running through molpro')
#endif

        ncore = 0

        do i = 1, 8
            ncore = ncore + icore(i)
            header(1-1+i) = icore(i)
            header(1+7+i) = iclos(i)
            header(1+15+i) = iact(i)
        end do

        header(1+24) = nelec
        header(1+25) = isyref
        header(1+26) = ms2
        header(1+27) = iabs(iblkq)
        header(1+28) = iabs(iseccr)
        header(1+29) = iabs(istat1)
        lhead = 30/intrel
        igrsav = 10*(name-1) + ifil
      
        name = igrsav/10
        ifil = igrsav-10*name

#ifdef MOLPRO
        call reserv(lhead+lden1+lden2+leps+lepsc,ifil,name,-1)
        call writem(header,lhead,ifil,name,0,label)
        call writem(den1,lden1,ifil,name,lhead,label)
        call writem(den2,lden2,ifil,name,lhead+lden1,label)
        call writem(eps,leps,ifil,name,lhead+lden1+lden2,label)
        if (ncore.ne.0) call writem(epsc,lepsc,ifil,name,&
                      lhead+lden1+lden2+leps,label)
#endif

        write(iout,20) istat1,isyref,name,ifil
!       call setf(2,recmc,1)
      ! ^- cgk: not sure what this does.
20      format(/' Gradient information for state',i2,'.',i1,&
                  ' saved on record  ',i8,'.',i1)
        call molpro_set_igrdsav(igrsav)

        !igrsav=nexfre(igrsav)
        !name=igrsav/10
        !ifil=igrsav-10*name
        !call reserv(lhead+lden1+lden2+leps+lepsc,ifil,name,-1)
!        open (unit = ifil, file = "fciqmc_forces_info", form='UNFORMATTED', access='sequential')
!        write(ifil) header, den1, den2, eps
        !call writem(den2,lden2,ifil,name,lhead+lden1,label)
        !call writem(eps,leps,ifil,name,lhead+lden1+lden2,label)
        !if (ncore.ne.0) call writem(epsc,lepsc,ifil,name,
        !>                          lhead+lden1+lden2+leps,label)
!        write(iout,*) "istat1, isyref, name, ifil", istat1,isyref,name,ifil
!        write(iout,*) "header, den1, den2, eps", header, den1, den2, eps
!20      format(/' Gradient information for state',i2,'.',i1,
        !>        ' saved on record  ',i8,'.',i1)
        !return

        ! Eliminate compiler warnings
        runused = den1(1); runused = den2(1); runused = eps(1); runused = epsc(1)


    end subroutine molpro_dump_mcscf_dens_for_grad

    ! To calculate the dipole moments, the casscf routine has to be called from molpro. i.e.

    ! {rhf;save,2103.2}
    ! {casscf,maxit=0;occ,10,4,4,1;closed,0,0,0,0;iprint,density,civector}
    ! gexpec,dm
    ! {fciqmc,iterations=10000,timestep=0.05,targetwalkers=10000,2RDMonFly,dipoles;core;orbital,2103.2}

    ! Notes: 
    !  o The orbitals used by the fciqmc module have to be the RHF orbitals, not the casscf natural ones.
    !  o The occ directive has to encompass the *whole* space
    !  o If core orbitals want to be frozen in the subsequent fciqmc calclation, then include these orbitals as 'closed'
    !        in the casscf call, and remove the 'core' command from fciqmc
    !  o If the system is too large to do a FCI casscf (which will hopefully normally be the case), then you can restrict
    !        the space by using the 'restict' directive directly after the 'wf' directive in the casscf, i.e.

    ! memory,128,m
    ! geometry={C;O,C,r};r=2.1316 bohr
    ! basis,VDZ;
    ! {rhf;save,2103.2}
    ! {casscf,maxit=0;occ,14,6,6,2;frozen,0,0,0,0;closed,2,0,0,0;
    ! wf,14,1,0;
    ! restrict,2,2,3.1,4.1;  !This defines a set of orbitals which must be doubly occupied in all configurations to cut down the
    ! restrict,2,2,1.2,1.3;  !size of the space, but ensure that the integrals are still calculated over the whole active space
    ! restrict,0,0,10.1,11.1,12.1,13.1,14.1;  !These are orbitals which must remain unoccupied to further 
                                              !reduce the size of the space
    ! restrict,0,0,3.2,4.2,5.2,6.2;
    ! restrict,0,0,3.3,4.3,5.3,6.3;
    ! iprint,density,civector}
    ! gexpec,dm
    ! {fciqmc,ITERATIONS=20,MAXATREF=50000,targetWALKERS=10000000;
    !   orbital,2103.2 }         
         !Note no 'core' directive included, since the core orbitals are 'closed' in the casscf and so will be ignored.

    subroutine CalcDipoles(rdm, Norm_1RDM)

        use rdm_data, only: rdm_t
#ifdef MOLPRO
        use GenRandSymExcitNUMod, only: RandExcitSymLabelProd, ClassCountInd
        use outputResult
        use RotateOrbsData, only: SpatOrbs
        use SymData, only: Sym_Psi,nSymLabels
        use SymExcitDataMod, only: SpinOrbSymLabel,SymLabelCounts2

        integer, dimension(nSymLabels) :: elements_assigned1, blockstart1
        real(dp) :: dipmom(3),znuc(3),zcor(3)
        real(dp), allocatable :: SymmetryPacked1RDM(:),zints(:,:)
        integer :: i, j, ipr, Sym_i, Sym_j, Sym_ij,
        integer :: posn1, isize, isyref, mxv, iout, ierr
        integer :: nt_frz(8), ntd_frz(8)
#endif
        ! The only thing needed is the 1RDM (normalized)
        type(rdm_t), intent(in) :: rdm
        real(dp), intent(in) :: Norm_1RDM

        character(len=*), parameter :: t_r='CalcDipoles'
        real(dp) :: runused

#ifdef MOLPRO

        if (iProcIndex .eq. 0) then
            iout = molpro_get_iout()
            ! We need to work out
            ! a) how molpro symmetry-packs UHF integrals (ROHF would be fine though)
            ! b) Ensure that the 1RDM is correctly calculated for UHF (It is always allocated as spatorbs)
            ! c) Modify this routine for contracting over spin-orbitals
            if (tOpenShell) call stop_all(t_r,'Not working for ROHF/UHF')

            isyref = Sym_Psi + 1 ! Spatial symmetry of the wavefunction.

            ! Size of symmetry packed arrays (spatial).
            isize = 0
            blockstart1(:) = 0
            do i = 0, nSymLabels-1
                ! Find position of each symmetry block in sym-packed forms of RDMS 1 & 2.
                blockstart1(i+1) = isize + 1 ! N.B. Len_1RDM still being updated in this loop.

                isize = isize + (SymLabelCounts2(2,ClassCountInd(1,i,0))* &
                    (SymLabelCounts2(2,ClassCountInd(1,i,0))+1))/2 ! Counting alpha orbitals.
            end do

            nt_frz(:) = 0
            ntd_frz(:) = 0
            do i = 0,nSymLabels-1
                nt_frz(i+1) = SymLabelCounts2(2,ClassCountInd(1,i,0))
            end do

            do i = 2, 8
                ntd_frz(i) = ntd_frz(i-1) + (nt_frz(i-1)*(nt_frz(i-1)+1))/2
            end do

            elements_assigned1(:) = 0
            allocate(SymmetryPacked1RDM(isize))
            SymmetryPacked1RDM(:) = 0.0_dp
            do i = 1, SpatOrbs ! Run over spatial orbitals, ALL ELECTRON ONLY.
                do j = 1, i ! i .ge. j
                    Sym_i = SpinOrbSymLabel(2*i)  ! Consider only alpha orbitals.
                    Sym_j = SpinOrbSymLabel(2*j)
                    Sym_ij = RandExcitSymLabelProd(Sym_i, Sym_j)
                    if (Sym_ij .eq. 0) then
                        if ((Sym_i+1) .gt. nSymLabels) call stop_all(t_r,'Error')
                        posn1 = blockstart1(Sym_i+1) + elements_assigned1(Sym_i+1)
                        if ((posn1 .gt. isize) .or. (posn1 .lt. 1)) then
                            call stop_all(t_r,'Error filling rdm')
                        end if

                        SymmetryPacked1RDM(posn1) = rdm%matrix(SymLabelListInv_rot(i),SymLabelListInv_rot(j))*Norm_1RDM
                        if (i .ne. j) then
                            ! Double the off-diagonal elements of the 1RDM, so
                            ! that when we contract over the symmetry packed
                            ! representation of the 1RDM, it is as if we are
                            ! also including the other half of the matrix
                            SymmetryPacked1RDM(posn1) = 2.0_dp*SymmetryPacked1RDM(posn1)
                        end if
                        elements_assigned1(Sym_i+1) = elements_assigned1(Sym_i+1) + 1
                    end if

                end do
            end do
                
            write(6,*) "Size of symmetry packed array: ", isize
            write(6,*) "Symmetry packed 1RDM: ", SymmetryPacked1RDM(:)
            dipmom(:) = 0.0_dp

            call clearvar('DMX')
            call clearvar('DMY')
            call clearvar('DMZ')

            allocate(zints(isize,3),stat=ierr)
            if (ierr .ne. 0) then
                write(6,*) "Alloc failed: ",ierr
                call stop_all(t_r,'Alloc failed')
            end if
            ! This now goes through an F77 wrapper file so that we can access the
            ! common blocks and check that the size of the symmetry packed arrays
            ! is correct.
            call GetDipMomInts(zints, isize, znuc, zcor, nt_frz, ntd_frz)

            do ipr = 1, 3
!                call pget(zints,ipr,znuc,zcor)

                ! Now, contract.
                do i = 1, isize
                    dipmom(ipr) = dipmom(ipr) - zints(i,ipr)*SymmetryPacked1RDM(i)
                end do
                dipmom(ipr) = dipmom(ipr) + znuc(ipr) - zcor(ipr)
            end do

            write(iout,'(/,"DIPOLE MOMENT:",1X,3f15.8,/)') dipmom(1:3)
            call output_result('FCIQMC','Dipole moment', dipmom(1:3), 1, isyref, numberformat='3f15.8', debye=.TRUE.)
            mxv = 1
            call setvar('DMX', dipmom(1), 'AU', 1, 1, mxv, -1)
            call setvar('DMY', dipmom(2), 'AU', 1, 1, mxv, -1)
            call setvar('DMZ', dipmom(3), 'AU', 1, 1, mxv, -1)
            deallocate(zints, SymmetryPacked1RDM)
        end if

#else
        call warning_neci(t_r, 'Cannot compute dipole moments if not running within molpro. Exiting...')
#endif

        ! Eliminate compiler warnings
        runused = norm_1rdm

    end subroutine CalcDipoles

end module rdm_estimators_old
