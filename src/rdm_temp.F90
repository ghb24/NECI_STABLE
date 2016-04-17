module rdm_temp

    ! This is a temporary module to prevent circular dependencies. These should
    ! be sorted later when altering the call structure of the routines.

    use constants

    implicit none

contains

    subroutine Finalise_2e_RDM(rdm)

        ! This routine sums, normalises, hermitian-ises, and prints the 2-RDMs.
        ! This may be called multiple times if we want to print multiple 2-RDMs.

        use FciMCData, only: Iter, PreviousCycles, IterRDMStart, tFinalRDMEnergy, VaryShiftIter
        use LoggingData, only: IterRDMonFly
        use LoggingData, only: tRDMInstEnergy, RDMEnergyIter
        use Parallel_neci, only: iProcIndex, MPISumAll
        use rdm_data_old, only: AllNodes_RDM_small, AllNodes_RDM_large, rdm_t
        use rdm_data, only: tCalc_RDMEnergy, tOpenShell
        use RotateOrbsData, only: SpatOrbs

        type(rdm_t), intent(inout) :: rdm

        integer :: i, ierr

        ! If Iter = 0, this means we have just read in the TwoRDM_POPS_a*** matrices into a***_RDM_full, and 
        ! just want to calculate the old energy.
        ! Don't need to do all this stuff here, because a***_RDM will be empty.

        if (((Iter+PreviousCycles) .ne. 0) .and. ((.not. tFinalRDMEnergy) .or. &
            ((.not. tCalc_RDMEnergy) .or. ((Iter - VaryShiftIter(1)) .le. IterRDMonFly) &
                      & .or. ((Iter-VaryShiftIter(inum_runs)) .le. IterRDMonFly) &
                      & .or. (mod((Iter+PreviousCycles-IterRDMStart)+1, RDMEnergyIter) .ne. 0)))) then

            ! Two different sizes arrays are needed, due to *_abab and *_baba
            ! RDMs being a bit bigger.
            allocate(AllNodes_RDM_small(((SpatOrbs*(SpatOrbs-1))/2), ((SpatOrbs*(SpatOrbs-1))/2)), stat=ierr)
            allocate(AllNodes_RDM_large(((SpatOrbs*(SpatOrbs+1))/2), ((SpatOrbs*(SpatOrbs+1))/2)), stat=ierr)
            
            ! The RDM arrays may be either inst or full, depending on whether
            ! we are calculating instantaneous energies or not.
            call MPISumAll(rdm%aaaa(:,:), AllNodes_RDM_small(:,:))
            rdm%aaaa(:,:) = AllNodes_RDM_small(:,:)

            call MPISumAll(rdm%abba(:,:), AllNodes_RDM_small(:,:))
            rdm%abba(:,:) = AllNodes_RDM_small(:,:)

            call MPISumAll(rdm%abab(:,:), AllNodes_RDM_large(:,:))
            rdm%abab(:,:) = AllNodes_RDM_large(:,:)

            if (tOpenShell) then
                call MPISumAll(rdm%bbbb(:,:), AllNodes_RDM_small(:,:))
                rdm%bbbb(:,:) = AllNodes_RDM_small(:,:)

                call MPISumAll(rdm%baab(:,:), AllNodes_RDM_small(:,:))
                rdm%baab(:,:) = AllNodes_RDM_small(:,:)

                call MPISumAll(rdm%baba(:,:), AllNodes_RDM_large(:,:))
                rdm%baba(:,:) = AllNodes_RDM_large(:,:)
            end if

            deallocate(AllNodes_RDM_small)
            deallocate(AllNodes_RDM_large)

            ! The TwoElRDM on the root is now the sum of all 'instantaneous' RDMs
            ! (summed over the energy update cycle). Whereas TwoElRDM_full is
            ! accumulated over the entire run.
            if (tRDMInstEnergy .and. (iProcIndex .eq. 0)) then

                ! We have to add the RDMs column by column. Adding the whole
                ! arrays in one go causes ifort to crash for large arrays,
                ! presumably because of large arrays being created on the
                ! stack, which ifort never likes.

                do i = 1, size(rdm%aaaa, 2)
                    rdm%aaaa_full(:,i) = rdm%aaaa_full(:,i) + rdm%aaaa(:,i)
                end do

                do i = 1, size(rdm%abba, 2)
                    rdm%abba_full(:,i) = rdm%abba_full(:,i) + rdm%abba(:,i)
                end do

                do i = 1, size(rdm%abab, 2)
                    rdm%abab_full(:,i) = rdm%abab_full(:,i) + rdm%abab(:,i)
                end do

                if (tOpenShell) then

                    do i = 1, size(rdm%bbbb, 2)
                        rdm%bbbb_full(:,i) = rdm%bbbb_full(:,i) + rdm%bbbb(:,i)
                    end do

                    do i = 1, size(rdm%baab, 2)
                        rdm%baab_full(:,i) = rdm%baab_full(:,i) + rdm%baab(:,i)
                    end do

                    do i = 1, size(rdm%baba, 2)
                        rdm%baba_full(:,i) = rdm%baba_full(:,i) + rdm%baba(:,i)
                    end do

                end if
            end if
        end if

    end subroutine Finalise_2e_RDM

    subroutine calc_2e_norms(rdm, Norm_2RDM_Inst, Norm_2RDM, Trace_2RDM)

        ! We want to 'normalise' the reduced density matrices. These are not
        ! even close to being normalised at the moment, because of the way
        ! they are calculated on the fly. They should be calculated from a
        ! normalised wavefunction.

        ! We also know that the trace of the two electron reduced density
        ! matrix must be equal to the number of electron pairs in the
        ! system = 1/2 N ( N - 1), so we can do the same for the 2RDM.

        use LoggingData, only: tRDMInstEnergy
        use rdm_data, only: tOpenShell
        use rdm_data_old, only: rdm_t
        use RotateOrbsData, only: SpatOrbs
        use SystemData, only: nel

        type(rdm_t), intent(inout) :: rdm
        real(dp), intent(out) :: Norm_2RDM_Inst, Norm_2RDM
        real(dp), intent(out) :: Trace_2RDM

        integer :: i
        real(dp) :: Trace_2RDM_Inst

        ! Find the current, unnormalised trace of each matrix.
        ! TODO: This can be merged into the spin averaging when everything is working.

        Trace_2RDM_Inst = 0.0_dp
        Trace_2RDM = 0.0_dp

        do i = 1, ((SpatOrbs*(SpatOrbs+1))/2)
            if (i .le. ((SpatOrbs*(SpatOrbs-1))/2)) then
                if (tRDMInstEnergy) then
                    Trace_2RDM_Inst = Trace_2RDM_Inst + rdm%aaaa(i,i)
                    if (tOpenShell) Trace_2RDM_Inst = Trace_2RDM_Inst + rdm%bbbb(i,i)
                end if

                Trace_2RDM = Trace_2RDM + rdm%aaaa_full(i,i)
                 if (tOpenShell) Trace_2RDM = Trace_2RDM + rdm%bbbb_full(i,i)
            end if

            if (tRDMInstEnergy) then
                Trace_2RDM_Inst = Trace_2RDM_Inst + rdm%abab(i,i)
                if (tOpenShell) Trace_2RDM_Inst = Trace_2RDM_Inst + rdm%baba(i,i)
            end if

            Trace_2RDM = Trace_2RDM + rdm%abab_full(i,i)
            if (tOpenShell) Trace_2RDM = Trace_2RDM + rdm%baba_full(i,i)
        end do

        Norm_2RDM_Inst = 0.0_dp
        Norm_2RDM = 0.0_dp

        if (tRDMInstEnergy) Norm_2RDM_Inst = ( (0.50_dp * (real(NEl,dp) * (real(NEl,dp) - 1.0_dp))) / Trace_2RDM_Inst )
        Norm_2RDM = ( (0.50_dp * (real(NEl,dp) * (real(NEl,dp) - 1.0_dp))) / Trace_2RDM )

    end subroutine calc_2e_norms

    subroutine Write_out_2RDM(rdm, rdm_label, Norm_2RDM, tNormalise, tMake_Herm)

        ! Writes out the 2-RDMs. If tNormalise is true, we print the normalised
        ! (hermitian) matrix. Otherwise we print the unnormalised 2-RDMs, and
        ! we print (in binary) both 2-RDM(Ind1,Ind2) and 2-RDM(Ind2,Ind1)
        ! because this matrix wont be hermitian.

        ! While, for instance, the TwoRDM_aaaa so far has actually been a sum
        ! of the aaaa elements and the bbbb elements.  We only want to print
        ! the aaaa elements.

        use FciMCData, only: tFinalRDMEnergy, Iter, PreviousCycles
        use LoggingData, only: tWriteMultRDMs
        use rdm_data_old, only: rdm_t
        use rdm_data, only: tOpenShell
        use RotateOrbsData, only: SpatOrbs
        use util_mod, only: get_unique_filename, get_free_unit, int_fmt

        type(rdm_t), intent(inout) :: rdm
        integer, intent(in) :: rdm_label
        real(dp), intent(in) :: Norm_2RDM
        logical, intent(in) :: tNormalise, tMake_Herm

        real(dp) :: Divide_Factor
        integer :: i, j, a, b, Ind1_aa, Ind1_ab, Ind2_aa, Ind2_ab
        integer :: No_Herm_Elements
        integer :: aaaa_RDM_unit, abab_RDM_unit, abba_RDM_unit
        integer :: bbbb_RDM_unit, baba_RDM_unit, baab_RDM_unit
        character(20) :: stem
        character(255) :: TwoRDM_aaaa_name, TwoRDM_abab_name, TwoRDM_abba_name
        character(255) :: TwoRDM_bbbb_name, TwoRDM_baba_name, TwoRDM_baab_name
        real(dp) :: Max_Error_Hermiticity, Sum_Error_Hermiticity, Sum_Herm_Percent 
        real(dp) :: Temp

        ! Eliminate compiler warnings (no functional impact)
        bbbb_RDM_unit = 0
        baba_RDM_unit = 0
        baab_RDM_unit = 0

        if (tNormalise) then
            write(6,'(1X,"Writing out the *normalised* 2 electron density matrix to a file...")')
            call neci_flush(6)
            ! This takes the TwoRDM_aaaa, and if tWriteMultPops is true (and given that we've put 
            ! .true. in the 3rd position, it'll find the next unused TwoRDM_aaaa.X file name.
#ifdef _MOLCAS_
            aaaa_RDM_unit = get_free_unit()
            call molcas_open(aaaa_RDM_unit,'TWORDM1')

            abab_RDM_unit = get_free_unit()
            call molcas_open(abab_RDM_unit,'TWORDM2')

            abba_RDM_unit = get_free_unit()
            call molcas_open(abba_RDM_unit,'TWORDM3')
#else
            write(stem, '("TwoRDM_aaaa_old.",'//int_fmt(rdm_label,0)//')') rdm_label
            call get_unique_filename(trim(stem), tWriteMultRDMs, .true., 1, TwoRDM_aaaa_name)
            aaaa_RDM_unit = get_free_unit()
            open(aaaa_RDM_unit, file=TwoRDM_aaaa_name, status='unknown')

            write(stem, '("TwoRDM_abab_old.",'//int_fmt(rdm_label,0)//')') rdm_label
            call get_unique_filename(trim(stem), tWriteMultRDMs, .true., 1, TwoRDM_abab_name)
            abab_RDM_unit = get_free_unit()
            open(abab_RDM_unit, file=TwoRDM_abab_name, status='unknown')

            write(stem, '("TwoRDM_abba_old.",'//int_fmt(rdm_label,0)//')') rdm_label
            call get_unique_filename(trim(stem), tWriteMultRDMs, .true., 1, TwoRDM_abba_name)
            abba_RDM_unit = get_free_unit()
            open(abba_RDM_unit, file=TwoRDM_abba_name, status='unknown')
#endif
            if (tOpenShell) then
                write(stem, '("TwoRDM_bbbb_old.",'//int_fmt(rdm_label,0)//')') rdm_label
                call get_unique_filename(trim(stem), tWriteMultRDMs, .true., 1, TwoRDM_bbbb_name)
                bbbb_RDM_unit = get_free_unit()
                open(bbbb_RDM_unit,file=TwoRDM_bbbb_name, status='unknown')

                write(stem, '("TwoRDM_baba_old.",'//int_fmt(rdm_label,0)//')') rdm_label
                call get_unique_filename(trim(stem), tWriteMultRDMs, .true., 1, TwoRDM_baba_name)
                baba_RDM_unit = get_free_unit()
                open(baba_RDM_unit, file=TwoRDM_baba_name, status='unknown')

                write(stem, '("TwoRDM_baab_old.",'//int_fmt(rdm_label,0)//')') rdm_label
                call get_unique_filename(trim(stem), tWriteMultRDMs, .true., 1, TwoRDM_baab_name)
                baab_RDM_unit = get_free_unit()
                open(baab_RDM_unit, file=TwoRDM_baab_name, status='unknown')
            end if

        else
            write(6,'("Writing out the *unnormalised* 2 electron density matrix to file for reading in")')
            call neci_flush(6)

            aaaa_RDM_unit = get_free_unit()
            write(stem, '("TwoRDM_POPS_aaaa_old.",'//int_fmt(rdm_label,0)//')') rdm_label
            open(aaaa_RDM_unit, file=trim(stem), status='unknown', form='unformatted')

            abab_RDM_unit = get_free_unit()
            write(stem, '("TwoRDM_POPS_abab_old.",'//int_fmt(rdm_label,0)//')') rdm_label
            open(abab_RDM_unit, file=trim(stem), status='unknown', form='unformatted')

            abba_RDM_unit = get_free_unit()
            write(stem, '("TwoRDM_POPS_abba_old.",'//int_fmt(rdm_label,0)//')') rdm_label
            open(abba_RDM_unit, file=trim(stem), status='unknown', form='unformatted')

            if (tOpenShell) then
                bbbb_RDM_unit = get_free_unit()
                write(stem, '("TwoRDM_POPS_bbbb_old.",'//int_fmt(rdm_label,0)//')') rdm_label
                open(bbbb_RDM_unit, file=trim(stem), status='unknown', form='unformatted')

                baba_RDM_unit = get_free_unit()
                write(stem, '("TwoRDM_POPS_baba_old.",'//int_fmt(rdm_label,0)//')') rdm_label
                open(baba_RDM_unit, file=trim(stem), status='unknown', form='unformatted')

                baab_RDM_unit = get_free_unit()
                write(stem, '("TwoRDM_POPS_baab_old.",'//int_fmt(rdm_label,0)//')') rdm_label
                open(baab_RDM_unit, file=trim(stem), status='unknown', form='unformatted')
            end if
        end if
       
        Max_Error_Hermiticity = 0.0_dp
        Sum_Error_Hermiticity = 0.0_dp
        Sum_Herm_Percent = 0.0_dp
        No_Herm_Elements = 0

        do i = 1, SpatOrbs

            do j = i, SpatOrbs

                Ind1_aa = ( ( (j-2) * (j-1) ) / 2 ) + i
                Ind1_ab = ( ( (j-1) * j ) / 2 ) + i

                do a = 1, SpatOrbs

                    do b = a, SpatOrbs

                        Ind2_aa = ( ( (b-2) * (b-1) ) / 2 ) + a
                        Ind2_ab = ( ( (b-1) * b ) / 2 ) + a

                        ! usually each element will have two contributions (from aaaa and bbbb).
                        ! we then need to divide each by 2.
                        ! but in cases where i and j, and a and b, are in the same spatial 
                        ! orbital, there will be only one contribution.
                        if (((i .eq. j) .and. (a .eq. b)) .or. tOpenShell) then
                            Divide_Factor = 1.0_dp
                        else
                            Divide_Factor = 2.0_dp
                        end if
                        
                        if ((i .ne. j) .and. (a .ne. b)) then

                            if ( (abs(rdm%aaaa_full(Ind1_aa,Ind2_aa)) > 1.0e-12_dp) .or. &
                                 (abs(rdm%aaaa_full(Ind2_aa,Ind1_aa)) > 1.0e-12_dp) ) then
                                ! If we're normalising (and have made the matrix hermitian) we only 
                                ! need to write out Ind1 < Ind2.
                                ! Otherwise we print out Ind1, Ind2 and Ind2, Ind1 so we can 
                                ! find the hermiticity error in the final matrix (after all runs).
                                if (tNormalise .and. (Ind1_aa .le. Ind2_aa)) then
                                    
                                    if ((abs((rdm%aaaa_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                - (rdm%aaaa_full(Ind2_aa,Ind1_aa)*Norm_2RDM))) .gt. Max_Error_Hermiticity) &
                                        Max_Error_Hermiticity = abs((rdm%aaaa_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                    - (rdm%aaaa_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                    Sum_Error_Hermiticity = Sum_Error_Hermiticity + &
                                                            abs((rdm%aaaa_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                    - (rdm%aaaa_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                    Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                        (abs((rdm%aaaa_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                - (rdm%aaaa_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / &
                                                        (abs((rdm%aaaa_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                + (rdm%aaaa_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / 2.0_dp) )
                                    No_Herm_Elements = No_Herm_Elements + 1                                                      

                                    if (tMake_Herm) then                                                            
                                        Temp = (rdm%aaaa_full(Ind1_aa,Ind2_aa) + rdm%aaaa_full(Ind2_aa,Ind1_aa)) / 2.0_dp

                                        rdm%aaaa_full(Ind1_aa,Ind2_aa) = Temp
                                        rdm%aaaa_full(Ind2_aa,Ind1_aa) = Temp
                                    end if

                                    if (tFinalRDMEnergy) then
                                        ! For the final calculation, the 2-RDMs will have been made hermitian.
                                        write(aaaa_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                                ( rdm%aaaa_full(Ind1_aa,Ind2_aa) * Norm_2RDM ) / Divide_Factor
                                    else
                                        ! If we're printing the 2-RDMs early (using writeRDMSEVERY), the actual 
                                        ! matrix will not be hermitian, but we want to print a hermitian version.
                                        ! Average the values here.
                                        write(aaaa_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                            ( ((rdm%aaaa_full(Ind1_aa,Ind2_aa) + rdm%aaaa_full(Ind2_aa,Ind1_aa))/2.0_dp) &
                                                                * Norm_2RDM ) / Divide_Factor
                                    end if
                                else if (.not. tNormalise) then
                                    ! For the popsfiles, print everything to binary.
                                    ! No divide factor, we just read them in as is.
                                    write(aaaa_RDM_unit) i, j, a, b, rdm%aaaa_full(Ind1_aa,Ind2_aa)
                                end if ! tNormalise
                            end if ! rdm%aaaa_full

                            if (tOpenShell) then

                                if ( (abs(rdm%bbbb_full(Ind1_aa,Ind2_aa)) > 1.0e-12_dp) .or. &
                                     (abs(rdm%bbbb_full(Ind2_aa,Ind1_aa)) > 1.0e-12_dp) ) then

                                    ! If we're normalising (and have made the matrix hermitian) we only 
                                    ! need to write out Ind1 < Ind2.
                                    ! Otherwise we print out Ind1, Ind2 and Ind2, Ind1 so we can 
                                    ! find the hermiticity error in the final matrix (after all runs).
                                    if (tNormalise .and. (Ind1_aa .le. Ind2_aa)) then

                                        if ((abs((rdm%bbbb_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                    - (rdm%bbbb_full(Ind2_aa,Ind1_aa)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                            Max_Error_Hermiticity = abs((rdm%bbbb_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                        - (rdm%bbbb_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                                abs((rdm%bbbb_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                        - (rdm%bbbb_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                            (abs((rdm%bbbb_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                    - (rdm%bbbb_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / &
                                                            (abs((rdm%bbbb_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                    + (rdm%bbbb_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / 2.0_dp) )
                                        No_Herm_Elements = No_Herm_Elements + 1

                                        if (tMake_Herm) then
                                            Temp = (rdm%bbbb_full(Ind1_aa,Ind2_aa) + rdm%bbbb_full(Ind2_aa,Ind1_aa)) / 2.0_dp

                                            rdm%bbbb_full(Ind1_aa,Ind2_aa) = Temp
                                            rdm%bbbb_full(Ind2_aa,Ind1_aa) = Temp
                                end if

                                        if (tFinalRDMEnergy) then
                                            ! For the final calculation, the 2-RDMs will have been made hermitian.
                                            write(bbbb_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                                    ( rdm%bbbb_full(Ind1_aa,Ind2_aa) * Norm_2RDM ) / Divide_Factor
                                        else
                                            ! If we're printing the 2-RDMs early (using writeRDMSEVERY), the actual 
                                            ! matrix will not be hermitian, but we want to print a hermitian version.
                                            ! Average the values here.
                                            write(bbbb_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                                ( ((rdm%bbbb_full(Ind1_aa,Ind2_aa) + rdm%bbbb_full(Ind2_aa,Ind1_aa))/2.0_dp) &
                                                                    * Norm_2RDM ) / Divide_Factor
                            end if
                                    else if (.not. tNormalise) then
                                        ! For the popsfiles, print everything to binary.
                                        ! no divide factor, we just read them in as is.
                                        write(bbbb_RDM_unit) i, j, a, b, rdm%bbbb_full(Ind1_aa,Ind2_aa)
                                    end if  ! tNormalise
                                end if  ! rdm%bbbb_full 
                            end if ! OpenShell

                            if (tOpenShell) then

                                if ( (abs(rdm%abba_full(Ind1_aa,Ind2_aa)) > 1.0e-12_dp) .or. &
                                     (abs(rdm%baab_full(Ind2_aa,Ind1_aa)) > 1.0e-12_dp) ) then

                                    if (tNormalise .and. (Ind1_aa .le. Ind2_aa)) then

                                        if ((abs((rdm%abba_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                    - (rdm%baab_full(Ind2_aa,Ind1_aa)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                            Max_Error_Hermiticity = abs((rdm%abba_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                    - (rdm%baab_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Error_Hermiticity = Sum_Error_Hermiticity + &
                                                                abs((rdm%abba_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                        - (rdm%baab_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Herm_Percent = Sum_Herm_Percent + &
                                                            (abs((rdm%abba_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                            - (rdm%baab_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / &
                                                            (abs((rdm%abba_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                            + (rdm%baab_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / 2.0_dp) )
                                        No_Herm_Elements = No_Herm_Elements + 1 


                                        if (tMake_Herm) then                                                            
                                            Temp = (rdm%abba_full(Ind1_aa,Ind2_aa) + rdm%baab_full(Ind2_aa,Ind1_aa)) / 2.0_dp
                                            rdm%abba_full(Ind1_aa,Ind2_aa) = Temp
                                            rdm%baab_full(Ind2_aa,Ind1_aa) = Temp
                                        end if

                                        if (tFinalRDMEnergy) then
                                            write(abba_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                                ( rdm%abba_full(Ind1_aa,Ind2_aa) * Norm_2RDM ) / Divide_Factor
                                        else
                                            write(abba_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                                ( ((rdm%abba_full(Ind1_aa,Ind2_aa) + rdm%abba_full(Ind2_aa,Ind1_aa))/2.0_dp) &
                                                                        * Norm_2RDM ) / Divide_Factor
                                        end if
                                    else if (.not.tNormalise) then
                                        write(abba_RDM_unit) i, j, a, b, rdm%abba_full(Ind1_aa,Ind2_aa) 
                                    end if
                                end if ! rdm%abba_full rdm%baab_full

                                if ( (abs(rdm%baab_full(Ind1_aa,Ind2_aa)) > 1.0e-12_dp) .or. &
                                     (abs(rdm%abba_full(Ind2_aa,Ind1_aa)) > 1.0e-12_dp) ) then

                                if (tNormalise .and. (Ind1_aa .le. Ind2_aa)) then

                                        if ((abs((rdm%baab_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                    - (rdm%abba_full(Ind2_aa,Ind1_aa)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                            Max_Error_Hermiticity = abs((rdm%baab_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                    - (rdm%abba_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                                abs((rdm%baab_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                        - (rdm%abba_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                            (abs((rdm%baab_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                            - (rdm%abba_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / &
                                                            (abs((rdm%baab_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                            + (rdm%abba_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / 2.0_dp) )
                                        No_Herm_Elements = No_Herm_Elements + 1

                                        if (tMake_Herm) then
                                            Temp = (rdm%baab_full(Ind1_aa,Ind2_aa) + rdm%abba_full(Ind2_aa,Ind1_aa)) / 2.0_dp
                                            rdm%baab_full(Ind1_aa,Ind2_aa) = Temp
                                            rdm%abba_full(Ind2_aa,Ind1_aa) = Temp  
                                        end if ! tMake_Herm = .true.

                                        if (tFinalRDMEnergy) then
                                            write(baab_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                                ( rdm%baab_full(Ind1_aa,Ind2_aa) * Norm_2RDM ) / Divide_Factor
                                        else
                                            write(baab_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                                ( ((rdm%baab_full(Ind1_aa,Ind2_aa) + rdm%baab_full(Ind2_aa,Ind1_aa))/2.0_dp) &
                                                                        * Norm_2RDM ) / Divide_Factor
                                        end if  ! tFinalRDMEnergy = .true./.false.

                                    else if (.not.tNormalise) then
                                        write(baab_RDM_unit) i, j, a, b, rdm%baab_full(Ind1_aa,Ind2_aa) 
                                    end if ! tNormalise
                                end if ! rdm%baab_full rdm%abba_full

                            else ! not tOpenShell

                                if ( (abs(rdm%abba_full(Ind1_aa,Ind2_aa)) > 1.0e-12_dp) .or. &
                                     (abs(rdm%abba_full(Ind2_aa,Ind1_aa)) > 1.0e-12_dp) ) then

                                    if (tNormalise .and. (Ind1_aa .le. Ind2_aa)) then
                                        if ((abs((rdm%abba_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                - (rdm%abba_full(Ind2_aa,Ind1_aa)*Norm_2RDM))) .gt. Max_Error_Hermiticity) &
                                        Max_Error_Hermiticity = abs((rdm%abba_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                - (rdm%abba_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                            abs((rdm%abba_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                    - (rdm%abba_full(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                        Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                        (abs((rdm%abba_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                        - (rdm%abba_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / &
                                                        (abs((rdm%abba_full(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                        + (rdm%abba_full(Ind2_aa,Ind1_aa)*Norm_2RDM)) / 2.0_dp) )
                                        No_Herm_Elements = No_Herm_Elements + 1

                                        if (tMake_Herm) then                                                            
                                            Temp = (rdm%abba_full(Ind1_aa,Ind2_aa) + rdm%abba_full(Ind2_aa,Ind1_aa)) / 2.0_dp
                                            rdm%abba_full(Ind1_aa,Ind2_aa) = Temp
                                            rdm%abba_full(Ind2_aa,Ind1_aa) = Temp
                                        end if

                                        if (tFinalRDMEnergy) then
                                            write(abba_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                                ( rdm%abba_full(Ind1_aa,Ind2_aa) * Norm_2RDM ) / Divide_Factor
                                        else
                                            write(abba_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                                ( ((rdm%abba_full(Ind1_aa,Ind2_aa) + rdm%abba_full(Ind2_aa,Ind1_aa))/2.0_dp) &
                                                                    * Norm_2RDM ) / Divide_Factor
                                        end if
                                    else if (.not. tNormalise) then
                                        write(abba_RDM_unit) i, j, a, b, rdm%abba_full(Ind1_aa,Ind2_aa) 
                                    end if
                                end if
                            end if 

                        end if  ! (i.ne.j) .and. (a.ne.b) 

                        if ( (abs(rdm%abab_full(Ind1_ab,Ind2_ab)) > 1.0e-12_dp) .or. &
                             (abs(rdm%abab_full(Ind2_ab,Ind1_ab)) > 1.0e-12_dp) ) then

                            if ((abs((rdm%abab_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                        - (rdm%abab_full(Ind2_ab,Ind1_ab)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                Max_Error_Hermiticity = abs((rdm%abab_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                        - (rdm%abab_full(Ind2_ab,Ind1_ab)*Norm_2RDM))

                            Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                    abs((rdm%abab_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                        - (rdm%abab_full(Ind2_ab,Ind1_ab)*Norm_2RDM))

                            Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                (abs((rdm%abab_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                    - (rdm%abab_full(Ind2_ab,Ind1_ab)*Norm_2RDM)) / &
                                                (abs((rdm%abab_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                   + (rdm%abab_full(Ind2_ab,Ind1_ab)*Norm_2RDM)) / 2.0_dp) )
                            No_Herm_Elements = No_Herm_Elements + 1                                                        

                            if (tMake_Herm) then                                                            
                                Temp = (rdm%abab_full(Ind1_ab,Ind2_ab) + rdm%abab_full(Ind2_ab,Ind1_ab)) / 2.0_dp
                                rdm%abab_full(Ind1_ab,Ind2_ab) = Temp
                                rdm%abab_full(Ind2_ab,Ind1_ab) = Temp
                            end if

                            if (tNormalise .and. (Ind1_ab .le. Ind2_ab)) then
                                if (tFinalRDMEnergy) then
                                    write(abab_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                        ( rdm%abab_full(Ind1_ab,Ind2_ab) * Norm_2RDM ) / Divide_Factor
                                else
                                    write(abab_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                        ( ((rdm%abab_full(Ind1_ab,Ind2_ab) + rdm%abab_full(Ind2_ab,Ind1_ab))/2.0_dp) &
                                                                * Norm_2RDM ) / Divide_Factor
                                end if
                            else if (.not. tNormalise) then
                                write(abab_RDM_unit) i, j, a, b, rdm%abab_full(Ind1_ab,Ind2_ab) 
                            end if
                        end if

                        if (tOpenShell .and. ( (a .ne. b) .or. (i .ne. j) ))then

                            if ( (abs(rdm%baba_full(Ind1_ab,Ind2_ab)) > 1.0e-12_dp).or. &
                                 (abs(rdm%baba_full(Ind2_ab,Ind1_ab)) > 1.0e-12_dp) ) then

                                if ((abs((rdm%baba_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                            - (rdm%baba_full(Ind2_ab,Ind1_ab)*Norm_2RDM))) .gt. Max_Error_Hermiticity) &
                                    Max_Error_Hermiticity = abs((rdm%baba_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                            - (rdm%baba_full(Ind2_ab,Ind1_ab)*Norm_2RDM))

                                Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                        abs((rdm%baba_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                            - (rdm%baba_full(Ind2_ab,Ind1_ab)*Norm_2RDM))

                                Sum_Herm_Percent = Sum_Herm_Percent + &
                                                    (abs((rdm%baba_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                        - (rdm%baba_full(Ind2_ab,Ind1_ab)*Norm_2RDM)) / &
                                                    (abs((rdm%baba_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                        + (rdm%baba_full(Ind2_ab,Ind1_ab)*Norm_2RDM)) / 2.0_dp) )
                                No_Herm_Elements = No_Herm_Elements + 1

                                if (tMake_Herm) then                                                            
                                    Temp = (rdm%baba_full(Ind1_ab,Ind2_ab) + rdm%baba_full(Ind2_ab,Ind1_ab)) / 2.0_dp
                                    rdm%baba_full(Ind1_ab,Ind2_ab) = Temp
                                    rdm%baba_full(Ind2_ab,Ind1_ab) = Temp
                                end if

                                if (tNormalise .and. (Ind1_ab .le. Ind2_ab)) then
                                    if (tFinalRDMEnergy) then
                                        write(baba_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                            ( rdm%baba_full(Ind1_ab,Ind2_ab) * Norm_2RDM ) / Divide_Factor
                                    else
                                        write(baba_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                            ( ((rdm%baba_full(Ind1_ab,Ind2_ab) + rdm%baba_full(Ind2_ab,Ind1_ab))/2.0_dp) &
                                                                    * Norm_2RDM ) / Divide_Factor
                                    end if
                                else if (.not. tNormalise) then
                                    write(baba_RDM_unit) i, j, a, b, rdm%baba_full(Ind1_ab,Ind2_ab)
                                end if
                            end if

                         else if (tOpenShell) then !a=b & i=j -> baba term saved in abab

                            if ( (abs(rdm%abab_full(Ind1_ab,Ind2_ab)) > 1.0e-12_dp) .or. &
                                (abs(rdm%abab_full(Ind2_ab,Ind1_ab)) > 1.0e-12_dp) ) then

                                if ((abs((rdm%abab_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                            - (rdm%abab_full(Ind2_ab,Ind1_ab)*Norm_2RDM))) .gt. Max_Error_Hermiticity) &
                                    Max_Error_Hermiticity = abs((rdm%abab_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                            - (rdm%abab_full(Ind2_ab,Ind1_ab)*Norm_2RDM))

                                Sum_Error_Hermiticity = Sum_Error_Hermiticity + &
                                                        abs((rdm%abab_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                            - (rdm%abab_full(Ind2_ab,Ind1_ab)*Norm_2RDM))

                                Sum_Herm_Percent = Sum_Herm_Percent + &
                                                    (abs((rdm%abab_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                        - (rdm%abab_full(Ind2_ab,Ind1_ab)*Norm_2RDM)) / &
                                                    (abs((rdm%abab_full(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                        + (rdm%abab_full(Ind2_ab,Ind1_ab)*Norm_2RDM)) / 2.0_dp) )
                                No_Herm_Elements = No_Herm_Elements + 1                                                        

                                if (tMake_Herm) then                                                             
                                    Temp = (rdm%abab_full(Ind1_ab,Ind2_ab) + rdm%abab_full(Ind2_ab,Ind1_ab)) / 2.0_dp
                                    rdm%abab_full(Ind1_ab,Ind2_ab) = Temp
                                    rdm%abab_full(Ind2_ab,Ind1_ab) = Temp
                                end if

                                if (tNormalise .and. (Ind1_ab .le. Ind2_ab)) then
                                    if (tFinalRDMEnergy) then
                                        write(baba_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                            ( rdm%abab_full(Ind1_ab,Ind2_ab) * Norm_2RDM ) / Divide_Factor
                                    else
                                        write(baba_RDM_unit,"(4I6,G25.17)") i, j, a, b, &
                                            ( ((rdm%abab_full(Ind1_ab,Ind2_ab) + rdm%abab_full(Ind2_ab,Ind1_ab))/2.0_dp) &
                                                                    * Norm_2RDM ) / Divide_Factor
                                    end if
                                end if
                            end if

                         end if ! tOpenShell

                    end do  
                end do  

            end do  
        end do 


        call neci_flush(aaaa_RDM_unit)
        call neci_flush(abab_RDM_unit)
        call neci_flush(abba_RDM_unit)
        close(aaaa_RDM_unit)
        close(abab_RDM_unit)
        close(abba_RDM_unit)
        if (tOpenShell) then
            close(bbbb_RDM_unit)
            close(baba_RDM_unit)
            close(baab_RDM_unit)
        end if

        if (tNormalise) then
            write(6,'(1X,"Stochastic error measures for RDM",1X,'//int_fmt(rdm_label)//',":")') rdm_label
            write(6,'(1X,I15,F30.20,5X,A41)') Iter+PreviousCycles, Max_Error_Hermiticity, &
                                            '(Iteration, MAX ABS ERROR IN HERMITICITY)'
            write(6,'(1X,I15,F30.20,5X,A41)') Iter+PreviousCycles, Sum_Error_Hermiticity, &
                                            '(Iteration, SUM ABS ERROR IN HERMITICITY)'
            write(6,'(1X,I15,F30.20,5X,A53,/)') Iter+PreviousCycles, Sum_Herm_Percent/real(No_Herm_Elements,dp), &
                                            '(Iteration, AVERAGE ABS PERCENTAGE HERMITICITY ERROR)'
        end if

    end subroutine Write_out_2RDM

    subroutine Write_spinfree_RDM(rdm, rdm_label, Norm_2RDM)

        ! Write out the spinfree 2RDM for MPQC.

        ! The final file will write out the spin-free 2RDM in the following way:
        !
        ! i j k l = \sum_{s,t} < a^+_{i,s} a^+_{j,t} a_{l,t} a_{k,s} >
        !
        ! where s and t are spin indices which are summed over.
        !

        use rdm_data_old, only: rdm_t
        use RotateOrbsData, only: SpatOrbs
        use util_mod, only: get_free_unit, int_fmt

        type(rdm_t), intent(in) :: rdm
        integer, intent(in) :: rdm_label
        real(dp), intent(in) :: Norm_2RDM

        integer :: i, j, a, b
        integer :: spinfree_RDM_unit
        character(255) :: RDM_filename

        write(6, '(1X,"Writing out the spinfree RDM for state",'//int_fmt(rdm_label,1)//')') rdm_label

        spinfree_RDM_unit = get_free_unit()
        write(RDM_filename, '("spinfree_TwoRDM_old.",'//int_fmt(rdm_label,0)//')') rdm_label
        open(spinfree_RDM_unit, file=trim(RDM_filename), status="replace")
        
        do j = 1, SpatOrbs
            do b = 1, SpatOrbs
                do a = 1, SpatOrbs
                    do i = 1, SpatOrbs
                        if (abs(Find_Spatial_2RDM_Chem(rdm, i, j, a, b, Norm_2RDM)) > 1.0e-12) then

                            write(spinfree_RDM_unit,"(4I15,F30.20)") i, a, j, b, &
                                   Find_Spatial_2RDM_Chem(rdm, i, j, a, b, Norm_2RDM)
                        end if
                    end do
                end do
            end do
        end do

        write(spinfree_RDM_unit, "(4I15,F30.20)") -1, -1, -1, -1, -1.0_dp
        close(spinfree_RDM_unit)

    end subroutine Write_spinfree_RDM

    function Find_Spatial_2RDM_Chem(rdm, p, q, r, s, Norm_2RDM) result(pqrs)

        ! This routine calculates the spatial orbital, chemical notation 2RDM component
        !                  D_pqrs = <Psi | p+ r+ s q | Psi>
        ! This is achieved by looking up the D_pr,qs component in the various spin-separated
        ! versions of the 2RDM that are currently stored with spatial orb numbering, in 
        ! physical notation (e.g. rdm%aaaa_full etc).

        ! To convert from spin to spatial orbitals, we need to apply the following:
        ! D_pr,qs = D_pr,qs(aaaa) + D_pr,qs(bbbb) + D_pr,qs(abab) + D_pr,qs(baba) (Eq. ***)

        ! We note now the following quirks of the rdm%aaaa_full-type arrays for the manner in
        ! which they store these components
        !     1. In most cases the current RDMs store the *sum* of the spin-inverted terms
        !          - ie, rdm%aaaa_full(pr,qs) contains the sum of the aaaa and bbbb contributions
        !     2. When p=r and q=s, there is only one contribution generated in NECI
        !          - ie, rdm%abab_full(pp,qq) contains only one of the two identical abab and baba contributions
        !          - Terms of this kind but be explicitly multiplied by two to satisfy Eq. *** above
        !          - This is stored in the "Mult_Factor"
        !     3. The existing 2RDMs only store terms with r>=p and s>=q
        !          - If we wish to look up a term with a different order to this, we must swap the
        !            order of the indices, considering the swapped spin and introducing appropriate signs
        !          - ie if p>r and s>q, D_pr,qs(abab) is found by looking up -D_rp,qs(abba)

        use rdm_data, only: tOpenShell
        use rdm_data_old, only: rdm_t
         
        type(rdm_t), intent(in) :: rdm
        integer, intent(in) :: p, q, r, s
        real(dp), intent(in) :: Norm_2RDM

        real(dp) :: pqrs
        real(dp) :: Mult_Factor
        integer :: Ind1_aa, Ind1_ab, Ind2_aa, Ind2_ab 

        pqrs = 0.0_dp

        if ((p .eq. r) .and. (q .eq. s)) then
            Mult_Factor = 2.0_dp
        else
            Mult_Factor = 1.0_dp
        end if
        
        if ((r .ge. p) .and. (s .ge. q)) then ! D_pr,qs correctly ordered.
            Ind1_aa = ( ( (r-2) * (r-1) ) / 2 ) + p
            Ind1_ab = ( ( (r-1) * r ) / 2 ) + p
            Ind2_aa = ( ( (s-2) * (s-1) ) / 2 ) + q
            Ind2_ab = ( ( (s-1) * s ) / 2 ) + q

            pqrs = pqrs + Mult_Factor*rdm%abab_full(Ind1_ab,Ind2_ab)*Norm_2RDM
            if (tOpenShell) pqrs = pqrs + Mult_Factor*rdm%baba_full(Ind1_ab,Ind2_ab)*Norm_2RDM

            if ((p .ne. r) .and. (q .ne. s)) then 
                pqrs = pqrs + Mult_Factor*rdm%aaaa_full(Ind1_aa,Ind2_aa)*Norm_2RDM
                if (tOpenShell) pqrs = pqrs + Mult_Factor*rdm%bbbb_full(Ind1_aa,Ind2_aa)*Norm_2RDM
            end if

        else if ((p .gt. r) .and. (s .gt. q)) then ! Need to reorder D_pr,qs to -D_rp,qs.
            Ind1_aa = ( ( (p-2) * (p-1) ) / 2 ) + r
            Ind2_aa = ( ( (s-2) * (s-1) ) / 2 ) + q

            pqrs = pqrs - Mult_Factor*rdm%abba_full(Ind1_aa,Ind2_aa)*Norm_2RDM
            if (tOpenShell) pqrs = pqrs - Mult_Factor*rdm%baab_full(Ind1_aa,Ind2_aa)*Norm_2RDM

            pqrs = pqrs - Mult_Factor*rdm%aaaa_full(Ind1_aa,Ind2_aa)*Norm_2RDM
            if (tOpenShell) pqrs = pqrs - Mult_Factor*rdm%bbbb_full(Ind1_aa,Ind2_aa)*Norm_2RDM

        else if ((r .gt. p) .and. (q .gt. s)) then ! Need to reorder D_pr,qs to -D_pr,sq.
            Ind1_aa = ( ( (r-2) * (r-1) ) / 2 ) + p
            Ind2_aa = ( ( (q-2) * (q-1) ) / 2 ) + s

            pqrs = pqrs - Mult_Factor*rdm%abba_full(Ind1_aa,Ind2_aa)*Norm_2RDM
            if (tOpenShell) pqrs = pqrs - Mult_Factor*rdm%baab_full(Ind1_aa,Ind2_aa)*Norm_2RDM

            pqrs = pqrs - Mult_Factor*rdm%aaaa_full(Ind1_aa,Ind2_aa)*Norm_2RDM
            if (tOpenShell) pqrs = pqrs - Mult_Factor*rdm%bbbb_full(Ind1_aa,Ind2_aa)*Norm_2RDM

        else ! Must need to reorder D_pr,qs to D_rp,sq.
            Ind1_aa = ( ( (p-2) * (p-1) ) / 2 ) + r
            Ind1_ab = ( ( (p-1) * p ) / 2 ) + r
            Ind2_aa = ( ( (q-2) * (q-1) ) / 2 ) + s
            Ind2_ab = ( ( (q-1) * q ) / 2 ) + s

            pqrs = pqrs + Mult_Factor*rdm%abab_full(Ind1_ab,Ind2_ab)*Norm_2RDM
            if (tOpenShell) pqrs = pqrs + Mult_Factor*rdm%baba_full(Ind1_ab,Ind2_ab)*Norm_2RDM

            if ((p .ne. r) .and. (q .ne. s)) then
                pqrs = pqrs + Mult_Factor*rdm%aaaa_full(Ind1_aa,Ind2_aa)*Norm_2RDM
                if (tOpenShell) pqrs = pqrs + Mult_Factor*rdm%bbbb_full(Ind1_aa,Ind2_aa)*Norm_2RDM
            end if

        end if
        
    end function Find_Spatial_2RDM_Chem

end module rdm_temp
