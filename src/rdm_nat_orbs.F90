module rdm_nat_orbs

    ! This contains routines related to natural orbitals calculations from the
    ! stochastically sampled reduced density matrices.

    use constants
    use RotateOrbsData, only: SymLabelCounts2_rot,SymLabelList2_rot
    use RotateOrbsData, only: SymLabelListInv_rot, NoOrbs, SpatOrbs
    use RotateOrbsData, only: SymLabelCounts2_rotTag, SymLabelList2_rotTag
    use RotateOrbsData, only: SymLabelListInv_rotTag

    implicit none

contains

    subroutine find_nat_orb_occ_numbers(rdm, irdm)

        use LoggingData, only: tPrintRODump
        use MemoryManager, only: LogMemAlloc
        use Parallel_neci, only: iProcIndex
        use rdm_data, only: one_rdm_t
        use RotateOrbsMod, only: FourIndInts, FourIndIntsTag
        use SystemData, only: tROHF, nbasis, G1, ARR, BRR

        ! Diagonalises the 1-RDM (rdm%matrix) so that, after this routine,
        ! rdm%matrix holds the eigenfunctions of the 1-RDM (the matrix
        ! transforming the MO's into the NOs). This also gets the NO
        ! occupation numbers (evaluse) and correlation entropy.

        type(one_rdm_t), intent(inout) :: rdm
        integer, intent(in) :: irdm

        integer :: ierr
        real(dp) :: SumDiag
        character(len=*), parameter :: t_r = 'find_nat_orb_occ_numbers'

        if (iProcIndex .eq. 0) then
            
            ! Diagonalises the 1-RDM. rdm%matrix goes in as the 1-RDM, comes out
            ! as the eigenvector of the 1-RDM (the matrix transforming the MO's
            ! into the NOs).
            call DiagRDM(rdm, SumDiag)

            ! Writes out the NO occupation numbers and evectors to files.
            call write_evales_and_transform_mat(rdm, irdm, SumDiag)

            if (tPrintRODump .and. tROHF) then               
                write(6,*) 'ROFCIDUMP not implemented for ROHF. Skip generation of ROFCIDUMP file.'
            else if (tPrintRODump) then
                allocate(FourIndInts(NoOrbs, NoOrbs, NoOrbs, NoOrbs), stat=ierr)
                call LogMemAlloc('FourIndInts',(NoOrbs**4), 8, t_r, &
                                                        FourIndIntsTag, ierr)
                if (ierr .ne. 0) call Stop_All(t_r, 'Problem allocating FourIndInts array,')

                ! Then, transform2ElInts.
                write(6,*) ''
                write(6,*) 'Transforming the four index integrals'
                call Transform2ElIntsMemSave_RDM(rdm%matrix, rdm%sym_list_no)

                write(6,*) 'Re-calculating the fock matrix'
                call CalcFOCKMatrix_RDM(rdm)

                write(6,*) 'Refilling the UMAT and TMAT2D'

                ! The ROFCIDUMP is also printed out in here.
                call RefillUMATandTMAT2D_RDM(rdm%matrix, rdm%sym_list_no)

                call neci_flush(6)

                call writebasis(6, G1, nbasis, ARR, BRR)

            end if
        end if

    end subroutine find_nat_orb_occ_numbers

    subroutine write_evales_and_transform_mat(rdm, irdm, SumDiag)

        use LoggingData, only: tNoNOTransform
        use rdm_data, only: tOpenShell, one_rdm_t
        use SystemData, only: nbasis, nel, BRR
        use UMatCache, only: gtID
        use util_mod, only: get_free_unit, int_fmt

        type(one_rdm_t), intent(in) :: rdm
        integer, intent(in) :: irdm
        real(dp), intent(in) :: SumDiag

        integer :: i, j, Evalues_unit, NatOrbs_unit, jSpat, jInd, NO_Number
        integer :: i_no,i_normal
        real(dp) :: Corr_Entropy, Norm_Evalues, SumN_NO_Occ
        logical :: tNegEvalue, tWrittenEvalue
        character(20) :: evals_filename, transform_filename

        if (tOpenShell) then
            Norm_Evalues = SumDiag/real(NEl,dp)
        else
            Norm_Evalues = 2.0_dp*(SumDiag/real(NEl,dp))
        end if

        ! Write out normalised Evalues to file and calculate the correlation
        ! entropy.
        Corr_Entropy = 0.0_dp

        Evalues_unit = get_free_unit()
        write(evals_filename, '("NO_OCC_NUMEBRS.",'//int_fmt(irdm,0)//')') irdm
        open(Evalues_unit, file=trim(evals_filename), status='unknown')

        write(Evalues_unit,'(A)') '# NOs (natural orbitals) ordered by occupation number' 
        write(Evalues_unit,'(A)') '# MOs (HF orbitals) ordered by energy' 
        write(Evalues_unit,'(A1,A5,A30,A20,A30)') '#','NO','NO OCCUPATION NUMBER','MO','MO OCCUPATION NUMBER'

        tNegEvalue = .false.
        SumN_NO_Occ = 0.0_dp
        NO_Number = 1

        do i = 1, NoOrbs
            if (tOpenShell) then
                write(Evalues_unit,'(I6,G35.17,I15,G35.17)') i, rdm%evalues(i)/Norm_Evalues, &
                                                                BRR(i), rdm%Rho_ii(i)
                if (rdm%evalues(i) .gt. 0.0_dp) then
                    Corr_Entropy = Corr_Entropy - ( abs(rdm%evalues(i)/ Norm_Evalues) &
                                                    * LOG(abs(rdm%evalues(i)/ Norm_Evalues)) )
                else
                    tNegEvalue = .true.
                end if
                if (i .le. NEl) SumN_NO_Occ = SumN_NO_Occ + (rdm%evalues(i)/Norm_Evalues)
            else
                write(Evalues_unit,'(I6,G35.17,I15,G35.17)') (2*i)-1,rdm%evalues(i)/Norm_Evalues, &
                                                            BRR((2*i)-1), rdm%Rho_ii(i)/2.0_dp
                if (rdm%evalues(i).gt.0.0_dp) then
                    Corr_Entropy = Corr_Entropy - (2.0_dp * ( abs(rdm%evalues(i)/Norm_Evalues) &
                                                    * LOG(abs(rdm%evalues(i)/Norm_Evalues)) ) )
                else
                    tNegEvalue = .true.
                end if
                write(Evalues_unit,'(I6,G35.17,I15,G35.17)') 2*i,rdm%evalues(i)/Norm_Evalues, &
                                                            BRR(2*i), rdm%Rho_ii(i)/2.0_dp
                if (i .le. (NEl/2)) SumN_NO_Occ = SumN_NO_Occ + (2.0_dp * (rdm%evalues(i)/Norm_Evalues))
            end if
        end do
        close(Evalues_unit)

        write(6,'(1X,A45,F30.20)') 'SUM OF THE N LARGEST NO OCCUPATION NUMBERS: ',SumN_NO_Occ
   
        write(6,'(1X,A20,F30.20)') 'CORRELATION ENTROPY', Corr_Entropy
        write(6,'(1X,A33,F30.20)') 'CORRELATION ENTROPY PER ELECTRON', Corr_Entropy / real(NEl,dp) 
        if (tNegEvalue) write(6,'(1X,"WARNING: Negative NO occupation numbers found.")')

        ! Write out the evectors to file.
        ! This is the matrix that transforms the molecular orbitals into the
        ! natural orbitals. rdm%evalues(i) corresponds to Evector NatOrbsMat(1:nBasis,i)
        ! We just want the rdm%evalues in the same order as above, but the
        ! 1:nBasis part (corresponding to the molecular orbitals), needs to
        ! refer to the actual orbital labels. Want these orbitals to preferably
        ! be in order, run through the orbital, need the position to find the
        ! corresponding NatOrbs element, use rdm%sym_list_inv_no.

        associate(arr_ind => rdm%sym_list_inv_no)

            if (.not. tNoNOTransform) then
                NatOrbs_unit = get_free_unit()
                write(transform_filename, '("NO_TRANSFORM.",'//int_fmt(irdm,0)//')') irdm
                open(NatOrbs_unit, file=trim(transform_filename), status='unknown')
                write(NatOrbs_unit,'(2A6,2A30)') '#   MO', 'NO', 'Transform Coeff', 'NO OCC NUMBER'
                ! write out in terms of spin orbitals, all alpha then all beta.
                NO_Number = 1
                do i_normal = 1, NoOrbs
                    i_no = arr_ind(i_normal)
                    tWrittenEvalue = .false.
                    do j = 1, nBasis
                        ! Here i corresponds to the natural orbital, and j to the
                        ! molecular orbital. i is actually the spin orbital in
                        ! this case.
                        if (tOpenShell) then
                            jInd = j
                        else
                            if (mod(j,2).ne.0) then
                                jInd = gtID(j)
                            else
                                cycle
                            end if
                        end if

                        if (tWrittenEvalue) then
                            if (rdm%matrix(arr_ind(jInd),i_no) .ne. 0.0_dp) &
                                write(NatOrbs_unit,'(2I6,G35.17)') j, NO_Number, rdm%matrix(arr_ind(jInd),i_no)
                        else
                            if (rdm%matrix(arr_ind(jInd),i_no) .ne. 0.0_dp) then
                                write(NatOrbs_unit,'(2I6,2G35.17)') j, NO_Number, rdm%matrix(arr_ind(jInd),i_no), &
                                                                    rdm%evalues(i_no)/Norm_Evalues
                                tWrittenEvalue = .true.
                            end if
                        end if
                    end do

                    NO_Number = NO_Number + 1
                    if (.not.tOpenShell) then
                        tWrittenEvalue = .false.
                        do j = 2, nBasis, 2
                            ! Here i corresponds to the natural orbital, and j to
                            ! the molecular orbital. i is actually the spin orbital
                            ! in this case.
                            jSpat = gtID(j)
                            if (tWrittenEvalue) then
                                if (rdm%matrix(arr_ind(jSpat),i_no) .ne. 0.0_dp) &
                                    write(NatOrbs_unit,'(2I6,G35.17)') j, NO_Number, &
                                                                    rdm%matrix(arr_ind(jSpat),i_no)
                            else
                                if (rdm%matrix(arr_ind(jSpat),i_no) .ne. 0.0_dp) then
                                    write(NatOrbs_unit,'(2I6,2G35.17)') j, NO_Number, rdm%matrix(arr_ind(jSpat),i_no), &
                                                                        rdm%evalues(i_no)/Norm_Evalues
                                    tWrittenEvalue = .true.
                                end if
                            end if

                        end do
                        NO_Number = NO_Number + 1
                    end if
                end do
                close(NatOrbs_unit)
            end if

        end associate

    end subroutine write_evales_and_transform_mat

    subroutine DiagRDM(rdm, SumTrace)

        ! The diagonalisation routine reorders the orbitals in such a way that
        ! the corresponding orbital labels are lost. In order to keep the spin
        ! and spatial symmetries, each symmetry must be fed into the
        ! diagonalisation routine separately. The best way to do this is to
        ! order the orbitals so that all the alpha orbitals follow all the beta
        ! orbitals, with the occupied orbitals first, in terms of symmetry, and
        ! the virtual second, also ordered by symmetry. This gives us
        ! flexibility w.r.t rotating only the occupied or only virtual and 
        ! looking at high spin states.

        use rdm_data, only: tOpenShell, one_rdm_t
        use MemoryManager, only: LogMemAlloc, LogMemDealloc
        use SystemData, only: G1, tUseMP2VarDenMat, tFixLz, iMaxLz

        type(one_rdm_t), intent(inout) :: rdm
        real(dp), intent(out) :: SumTrace

        real(dp) :: SumDiagTrace
        real(dp), allocatable :: WORK2(:), EvaluesSym(:), NOMSym(:,:)
        integer :: ierr, i, j, spin, Sym, LWORK2, WORK2Tag, SymStartInd, NoSymBlock
        integer :: EvaluesSymTag, NOMSymTag, k, MaxSym
        logical :: tDiffSym, tDiffLzSym
        character(len=*), parameter :: t_r = 'DiagRDM'

        ! Test that we're not breaking symmetry.
        ! And calculate the trace at the same time.
        SumTrace = 0.0_dp

        do i = 1, NoOrbs
            do j = 1, NoOrbs
                tDiffSym = .false.
                tDiffLzSym = .false.
                if (tOpenShell) then
                    if ((int(G1(SymLabelList2_rot(i))%sym%S).ne.&
                        int(G1(SymLabelList2_rot(j))%sym%S))) tDiffSym = .true.
                    if ((int(G1(SymLabelList2_rot(i))%Ml).ne.&
                        int(G1(SymLabelList2_rot(j))%Ml))) tDiffLzSym = .true.
                else
                    if ((int(G1(2*SymLabelList2_rot(i))%sym%S).ne.&
                        int(G1(2*SymLabelList2_rot(j))%sym%S))) tDiffSym = .true.
                    if ((int(G1(2*SymLabelList2_rot(i))%Ml).ne.&
                        int(G1(2*SymLabelList2_rot(j))%Ml))) tDiffLzSym = .true.
                end if
                if (tDiffSym) then
                    if (abs(rdm%matrix(i,j)).ge.1.0E-15_dp) then
                        write(6,'(6A8,A20)') 'i','j','Label i','Label j','Sym i',&
                                                                'Sym j','Matrix value'
                        if (tOpenShell) then                                                              
                            write(6,'(6I3,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j),&
                                    int(G1(SymLabelList2_rot(i))%sym%S),&
                                    int(G1(SymLabelList2_rot(j))%sym%S),rdm%matrix(i,j)
                        else
                            write(6,'(6I3,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j),&
                                    int(G1(2*SymLabelList2_rot(i))%sym%S),&
                                    int(G1(2*SymLabelList2_rot(j))%sym%S),rdm%matrix(i,j)
                        end if
                        if (tUseMP2VarDenMat) then
                            write(6,*) '**WARNING** - There is a non-zero 1-RDM &
                            &value between orbitals of different symmetry.'
                            write(6,*) 'These elements will be ignored, and the symmetry &
                            &maintained in the final transformation matrix.'
                        else
                            write(6,*) 'k,SymLabelList2_rot(k),SymLabelListInv_rot(k)'
                            do k = 1,NoOrbs
                                write(6,*) k,SymLabelList2_rot(k),SymLabelListInv_rot(k)
                            end do
                            call neci_flush(6)
                            call Stop_All(t_r,'Non-zero rdm%matrix value between &
                            &different symmetries.')
                        end if
                    end if
                    rdm%matrix(i,j) = 0.0_dp
                end if
                if (tDiffLzSym) then
                    if (abs(rdm%matrix(i,j)).ge.1.0E-15_dp) then
                        write(6,'(6A8,A40)') 'i','j','Label i','Label j','Lz i',&
                                                                'Lz j','Matrix value'
                        if (tOpenShell) then                                                              
                            write(6,'(6I8,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j),&
                                    int(G1(SymLabelList2_rot(i))%Ml),&
                                    int(G1(SymLabelList2_rot(j))%Ml),rdm%matrix(i,j)
                        else
                            write(6,'(6I8,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j),&
                                    int(G1(2*SymLabelList2_rot(i))%Ml),&
                                    int(G1(2*SymLabelList2_rot(j))%Ml),rdm%matrix(i,j)
                        end if
                        write(6,'(A)') ' **WARNING** - There is a non-zero 1-RDM element &
                        &between orbitals of different Lz symmetry.'
                    end if
                end if
            end do
            SumTrace = SumTrace + rdm%matrix(i,i)
        end do

        write(6,*) ''
        write(6,*) 'Calculating eigenvectors and eigenvalues of the 1-RDM'
        call neci_flush(6)

        ! If we want to maintain the symmetry, we cannot have all the orbitals
        ! jumbled up when the  diagonaliser reorders the eigenvectors. Must
        ! instead feed each symmetry block in separately. This means that
        ! although the transformed orbitals are jumbled within the symmetry blocks, 
        ! the symmetry labels are all that are relevant and these are unaffected.
        Sym = 0
        LWORK2 = -1
        if (tOpenShell) then
            if (tFixLz) then
                MaxSym = (16 * ( ( 2 * iMaxLz ) + 1 ) ) - 1
            else
                MaxSym = 15
            end if
        else
            if (tFixLz) then
                MaxSym = (8 * ( ( 2 * iMaxLz ) + 1 ) ) - 1
            else
                MaxSym = 7
            end if
        end if
        do while (Sym .le. MaxSym)

            NoSymBlock = SymLabelCounts2_rot(2,Sym+1)

            SymStartInd = SymLabelCounts2_rot(1,Sym+1) - 1
            ! This is one less than the index that the symmetry starts, so that when we 
            ! run through i=1,..., we can start at SymStartInd+i.

            if (NoSymBlock .gt. 1) then

                allocate(NOMSym(NoSymBlock,NoSymBlock),stat=ierr)
                call LogMemAlloc('NOMSym',NoSymBlock**2,8,t_r,NOMSymTag,ierr)
                if (ierr .ne. 0) call Stop_All(t_r,"Problem allocating NOMSym.")
                allocate(EvaluesSym(NoSymBlock),stat=ierr)
                call LogMemAlloc('EvaluesSym',NoSymBlock,8,t_r,EvaluesSymTag,ierr)
                if (ierr .ne. 0) call Stop_All(t_r,"Problem allocating EvaluesSym.")

                LWORK2 = 3*NoSymBlock + 1
                allocate(WORK2(LWORK2),stat=ierr)
                call LogMemAlloc('WORK2',LWORK2,8,t_r,WORK2Tag,ierr)
                if (ierr .ne. 0) call Stop_All(t_r,"Problem allocating WORK2.")

                do j = 1,NoSymBlock
                    do i = 1,NoSymBlock
                        NOMSym(i,j)=rdm%matrix(SymStartInd+i,SymStartInd+j)
                    end do
                end do

                call dsyev('V', 'L', NoSymBlock, NOMSym, NoSymBlock, EvaluesSym, WORK2, LWORK2, ierr)
                ! NOMSym goes in as the original NOMSym, comes out as the 
                ! eigenvectors (Coefficients).
                ! EvaluesSym comes out as the eigenvalues in ascending order.

                do i = 1, NoSymBlock
                    rdm%evalues(SymStartInd+i) = EvaluesSym(NoSymBlock-i+1)
                end do

                ! CAREFUL if eigenvalues are put in ascending order, this may not be 
                ! correct, with the labelling system.
                ! may be better to just take coefficients and transform TMAT2DRot 
                ! in transform2elints.
                ! a check that comes out as diagonal is a check of this routine anyway.

                do j = 1, NoSymBlock
                    do i = 1, NoSymBlock
                        rdm%matrix(SymStartInd+i,SymStartInd+j) = NOMSym(i,NoSymBlock-j+1)
                    end do
                end do
                ! Directly fill the coefficient matrix with the eigenvectors from 
                ! the diagonalization.

                deallocate(WORK2)
                call LogMemDealloc(t_r, WORK2Tag)

                deallocate(NOMSym)
                call LogMemDealloc(t_r, NOMSymTag)

                deallocate(EvaluesSym)
                call LogMemDealloc(t_r, EvaluesSymTag)

            else if (NoSymBlock .eq. 1) then
                ! The eigenvalue is the lone value, while the eigenvector is 1.

                rdm%evalues(SymStartInd+1) = rdm%matrix(SymStartInd+1, SymStartInd+1)
                rdm%matrix(SymStartInd+1, SymStartInd+1) = 1.0_dp
            end if

            Sym = Sym + 1
        end do

        write(6,*) 'Matrix diagonalised'
        call neci_flush(6)

        SumDiagTrace = 0.0_dp
        do i = 1, NoOrbs
            SumDiagTrace = SumDiagTrace + rdm%evalues(i)
        end do

        if ((abs(SumDiagTrace-SumTrace)).gt.1.0_dp) then
            write(6,*) 'Sum of diagonal 1-RDM elements : ',SumTrace
            write(6,*) 'Sum of eigenvalues : ',SumDiagTrace
            write(6,*) 'WARNING : &
            &The trace of the 1RDM matrix before diagonalisation is '
            write(6,*) 'not equal to that after.'
        end if

        ! The MO's still correspond to SymLabelList2_rot.
        ! Although the NO's are slightly jumbled, they are only jumbled within
        ! their symmetry blocks. They still correspond to the symmetries of
        ! SymLabelList2_rot, which is the important part.

        ! But in order to look at the output, it is easier to consider them in
        ! terms of highest occupied to lowest occupied - i.e. in terms of the
        ! NO eigenvalues (occupation numbers).
        call order_one_rdm(rdm)

    end subroutine DiagRDM

    subroutine order_one_rdm(rdm)

        ! Here, if symmetry is kept, we are going to have to reorder the
        ! eigenvectors according to the size of the eigenvalues, while taking
        ! the orbital labels (and therefore symmetries) with them. This will
        ! be put back into MP2VDM from MP2VDMTemp.

        ! Want to reorder the eigenvalues from largest to smallest, taking the
        ! eigenvectors with them and the symmetry as well. If using spin
        ! orbitals, do this for the alpha spin and then the beta.

        ! The newly sorted SymLabelList2_rot and SymLabelListInv_rot will be
        ! stored in rdm%sym_list_no and rdm%sym_list_inv_no (as they are
        ! specific to this RDM).

        use MemoryManager, only: LogMemAlloc
        use rdm_data, only: tOpenShell, one_rdm_t
        use sort_mod, only: sort
        use SystemData, only: nbasis

        type(one_rdm_t), intent(inout) :: rdm

        integer :: spin, i, j, ierr, StartSort, EndSort
        character(len=*), parameter :: t_r = 'order_one_rdm'
        integer, allocatable :: SymLabelList_temp(:)
        real(dp), allocatable :: one_rdm_Temp(:,:), EvaluesTemp(:)
        integer :: Orb, New_Pos

        ! Temporary arrays.
        allocate(one_rdm_Temp(NoOrbs,NoOrbs), stat=ierr)
        allocate(SymLabelList_temp(NoOrbs), stat=ierr)
        allocate(EvaluesTemp(NoOrbs), stat=ierr)

        ! This is just a temporary array for some sorting.
        SymLabelList_temp = SymLabelList2_rot

        ! This object will contain the updated version of SymLabelList2_rot,
        ! specifically for this RDM, and ordered with the 1-RDM eigenvalues.
        rdm%sym_list_no = SymLabelList2_rot

        StartSort = 1
        EndSort = SpatOrbs

        ! Unfortunately this sort routine orders the orbitals in ascending
        ! order... which is not quite what we want.  Just remember this when
        ! printing out the Evalues.
        call sort (rdm%EValues(startSort:endSort), &
                   rdm%matrix(1:NoOrbs, startSort:endSort), &
                   rdm%sym_list_no(startSort:endSort))

        if (tOpenShell) then                  
            StartSort = SpatOrbs + 1
            EndSort = nBasis

            call sort(rdm%EValues(startSort:endSort), &
                      rdm%matrix(1:NoOrbs, startSort:endSort), &
                      rdm%sym_list_no(startSort:endSort))

        end if                       

        ! We now have the NO's ordered according to the size of their Evalues
        ! (occupation numbers).  This will have jumbled up their symmetries.
        ! Want to reorder the MO's to match this ordering (so that we only
        ! have one SymLabelList array).

        ! Need a new SymLabelListInv_rot too. This will be stored in
        ! rdm%sym_list_inv_no.
        rdm%sym_list_inv_no = 0
        do i = 1, NoOrbs
            rdm%sym_list_inv_no(rdm%sym_list_no(NoOrbs-i+1)) = i
        end do

        one_rdm_Temp = rdm%matrix
        rdm%matrix = 0.0_dp

        do i = 1, NoOrbs
            do j = 1, NoOrbs

                ! In position j, the MO orbital Orb is currently there.
                Orb = SymLabelList_temp(j)

                ! Want to move it to the position the NO's are in.
                New_Pos = rdm%sym_list_inv_no(Orb)

                ! But we also want to reverse the order of everything... 
                rdm%matrix(New_Pos, NoOrbs-i+1) = one_rdm_Temp(j,i)
            end do
        end do

        SymLabelList_temp = rdm%sym_list_no
        EvaluesTemp = rdm%evalues
        do i = 1, NoOrbs
            rdm%sym_list_no(i) = SymLabelList_temp(NoOrbs-i+1)
            rdm%evalues(i) = EvaluesTemp(NoOrbs-i+1)
        end do

        deallocate(one_rdm_Temp)
        deallocate(SymLabelList_temp)
        deallocate(EvaluesTemp)

    end subroutine order_one_rdm

    subroutine Transform2ElIntsMemSave_RDM(one_rdm, sym_list)

        ! This is an M^5 transform, which transforms all the two-electron
        ! integrals into the new basis described by the Coeff matrix.
        ! This is very memory inefficient and currently does not use any
        ! spatial symmetry information.

        use IntegralsData, only: umat
        use MemoryManager, only: LogMemAlloc, LogMemDealloc
        use RotateOrbsMod, only: FourIndInts
        use UMatCache, only: UMatInd

        real(dp), intent(in) :: one_rdm(:,:)
        integer, intent(in) :: sym_list(:)

        integer :: i,j,k,l,a,b,g,d,ierr,Temp4indintsTag,a2,d2,b2,g2
        real(dp), allocatable :: Temp4indints(:,:)
        character(len=*), parameter :: t_r= 'Transform2ElIntsMemSave_RDM'
#ifdef __CMPLX
        call stop_all('Transform2ElIntsMemSave_RDM', &
                    'Rotating orbitals not implemented for complex orbitals.')
#endif
        
        ! Zero arrays from previous transform.

        allocate(Temp4indints(NoOrbs,NoOrbs),stat=ierr)
        call LogMemAlloc('Temp4indints',NoOrbs**2,8,&
                            'Transform2ElIntsMemSave_RDM',Temp4indintsTag,ierr)
        if (ierr .ne. 0) call Stop_All('Transform2ElIntsMemSave_RDM',&
                                    'Problem allocating memory to Temp4indints.')
 
        FourIndInts(:,:,:,:) = 0.0_dp

        ! Calculating the two-transformed, four index integrals.

        ! The untransformed <alpha beta | gamma delta> integrals are found from 
        ! UMAT(UMatInd(i,j,k,l)
        ! All our arrays are in spin orbitals - if tStoreSpinOrbs is false,
        ! UMAT will be in spatial orbitals - need to account for this.

        ! Running through 1,NoOrbs - the actual orbitals corresponding to that
        ! index are given by sym_list.

        do b = 1, NoOrbs
            b2 = sym_list(b)
            do d = 1, NoOrbs
                d2 = sym_list(d)
                do a = 1, NoOrbs
                    a2 = sym_list(a)
                    do g = 1, NoOrbs
                        g2 = sym_list(g)

                        ! UMatInd in physical notation, but FourIndInts in
                        ! chemical (just to make it more clear in these
                        ! transformations). This means that here, a and g are
                        ! interchangable, and so are b and d.
                        FourIndInts(a,g,b,d) = real(UMAT(UMatInd(a2,b2,g2,d2)),dp)
                    end do
                end do

                Temp4indints(:,:) = 0.0_dp
                call dgemm('T', 'N', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, one_rdm, NoOrbs, &
                            FourIndInts(1:NoOrbs,1:NoOrbs,b,d), NoOrbs, 0.0_dp,&
                            Temp4indints(1:NoOrbs,1:NoOrbs), NoOrbs)

                ! Temp4indints(i,g) comes out of here, so to transform g to k, 
                ! we need the transpose of this.

                call dgemm('T', 'T', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, one_rdm, NoOrbs, &
                            Temp4indints(1:NoOrbs,1:NoOrbs), NoOrbs, 0.0_dp,&
                            FourIndInts(1:NoOrbs,1:NoOrbs,b,d), NoOrbs)
                ! Get Temp4indits02(i,k).
            end do
        end do
        
        ! Calculating the 3 transformed, 4 index integrals.
        ! 01=a untransformed,02=b,03=g,04=d.
        do i = 1, NoOrbs
            do k = 1, NoOrbs

                Temp4indints(:,:) = 0.0_dp
                call dgemm('T', 'N', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, one_rdm, NoOrbs, &
                            FourIndInts(i,k,1:NoOrbs,1:NoOrbs), NoOrbs, 0.0_dp,&
                            Temp4indints(1:NoOrbs,1:NoOrbs), NoOrbs)

                call dgemm('T', 'T', NoOrbs, NoOrbs, NoOrbs, 1.0_dp, one_rdm, &
                            NoOrbs,Temp4indints(1:NoOrbs,1:NoOrbs), NoOrbs, 0.0_dp,&
                            FourIndInts(i,k,1:NoOrbs,1:NoOrbs), NoOrbs)
            end do
        end do

        deallocate(Temp4indints)
        call LogMemDeAlloc('Transform2ElIntsMemSave_RDM',Temp4indintsTag)
 
    end subroutine Transform2ElIntsMemSave_RDM

    subroutine CalcFOCKMatrix_RDM(rdm)

        ! Calculate the fock matrix in the natural orbital basis.

        use MemoryManager, only: LogMemAlloc, LogMemDealloc
        use rdm_data, only: one_rdm_t
        use SystemData, only: nbasis, ARR, BRR, tStoreSpinOrbs

        type(one_rdm_t), intent(in) :: rdm

        integer :: i, j, k, l, a, b, ierr, ArrDiagNewTag
        real(dp) :: FOCKDiagSumHF, FOCKDiagSumNew
        character(len=*), parameter :: t_r= 'CalcFOCKMatrix_RDM'
        real(dp), allocatable :: ArrDiagNew(:)

        ! This subroutine calculates and writes out the fock matrix for the
        ! transformed orbitals.
        ! ARR is originally the fock matrix in the HF basis.
        ! ARR(:,1) - ordered by energy, ARR(:,2) - ordered by spin-orbital index.

        ! When transforming the orbitals into approximate natural orbitals, we
        ! want to save memory, so don't bother calculating the whole matrix,
        ! just the diagonal elements that we actually need.

        allocate(ArrDiagNew(nBasis),stat=ierr)
        if (ierr .ne. 0) call Stop_All(t_r,'Problem allocating ArrDiagNew array,')
        call LogMemAlloc('ArrDiagNew',nBasis,8,t_r,ArrDiagNewTag,ierr)
        ArrDiagNew(:)=0.0_dp                     

        ! First calculate the sum of the diagonal elements, ARR.
        ! Check if this is already being done.
        FOCKDiagSumHF = 0.0_dp
        do a = 1, nBasis        
            FOCKDiagSumHF = FOCKDiagSumHF+Arr(a,2)
        end do

        ! Then calculate the fock matrix in the transformed basis, and the sum
        ! of the new diagonal elements.

        FOCKDiagSumNew = 0.0_dp
        do j = 1, NoOrbs
            l = rdm%sym_list_no(j)
            if (tStoreSpinOrbs) then
                ArrDiagNew(l) = 0.0_dp
            else
                ArrDiagNew(2*l) = 0.0_dp
                ArrDiagNew((2*l)-1) = 0.0_dp
            end if
            do a = 1, NoOrbs
                b = rdm%sym_list_no(a)
                if (tStoreSpinOrbs) then
                    ArrDiagNew(l)=ArrDiagNew(l)+(rdm%matrix(a,j)*ARR(b,2)*rdm%matrix(a,j))
                else
                    ArrDiagNew(2*l)=ArrDiagNew(2*l)+(rdm%matrix(a,j)*ARR(2*b,2)*rdm%matrix(a,j))
                    ArrDiagNew((2*l)-1)=ArrDiagNew((2*l)-1)+(rdm%matrix(a,j)*ARR((2*b)-1,2)*rdm%matrix(a,j))
                end if
            end do
            if (tStoreSpinOrbs) then
                FOCKDiagSumNew = FOCKDiagSumNew + (ArrDiagNew(l))
            else
                FOCKDiagSumNew = FOCKDiagSumNew + (ArrDiagNew(2*l))
                FOCKDiagSumNew = FOCKDiagSumNew + (ArrDiagNew((2*l)-1))
            end if
        end do
        ! If we are truncation the virtual space, only the unfrozen entries will 
        ! be transformed.

        ! Refill ARR(:,1) (ordered in terms of energies), and ARR(:,2)
        ! (ordered in terms of orbital number).
        ! ARR(:,2) needs to be ordered in terms of symmetry and then energy
        ! (like SymLabelList), so currently this ordering will not be correct
        ! when reading in qchem INTDUMPs as the orbital number ordering is by energy.

        do j = 1,nBasis
            ARR(j,2) = ArrDiagNew(j)
            ARR(j,1) = ArrDiagNew(BRR(j))
        end do

        deallocate(ArrDiagNew)
        call LogMemDealloc(t_r, ArrDiagNewTag)

    end subroutine CalcFOCKMatrix_RDM

    subroutine RefillUMATandTMAT2D_RDM(one_rdm, sym_list)

        ! This routine refills these to more easily write out the ROFCIDUMP,
        ! and originally to be able to continue a calculation (although I doubt
        ! this works at the moment).

        ! UMat is in spin or spatial orbitals, TMAT2D only spin.

        use IntegralsData, only: umat
        use MemoryManager, only: LogMemAlloc, LogMemDealloc
        use OneEInts, only: TMAT2D
        use rdm_data, only: tRotatedNOs
        use RotateOrbsMod, only: FourIndInts
        use SystemData, only: nbasis, tStoreSpinOrbs
        use UMatCache, only: UMatInd

        real(dp), intent(in) :: one_rdm(:,:)
        integer, intent(in) :: sym_list(:)

        integer :: l, k, j, i, a, b, g, d, c
        integer :: nBasis2, TMAT2DPartTag, ierr
        real(dp) :: NewTMAT
        real(dp), allocatable :: TMAT2DPart(:,:)
        character(len=*), parameter :: t_r = 'RefillUMATandTMAT2D_RDM'

#ifdef __CMPLX
        call stop_all('RefillUMATandTMAT2D_RDM', &
                    'Rotating orbitals not implemented for complex orbitals.')
#endif

        ! TMAT2D is always in spin orbitals.
        allocate(TMAT2DPart(nBasis,nBasis),stat=ierr)
        if (ierr .ne. 0) call Stop_All(t_r,'Problem allocating TMAT2DPart array,')
        call LogMemAlloc('TMAT2DPart',nBasis*nBasis,8,&
                                    'RefillUMAT_RDM',TMAT2DPartTag,ierr)
        TMAT2DPart = 0.0_dp

        ! Make the UMAT elements the four index integrals.
        ! These are calculated by transforming the HF orbitals using the
        ! coefficients that have been found.
        do l = 1, NoOrbs
            d = sym_list(l)
            do k = 1, NoOrbs
                b = sym_list(k)
                do j = 1, NoOrbs
                    g = sym_list(j)
                    do i = 1, NoOrbs
                        a = sym_list(i)
                        ! The FourIndInts are in chemical notation, the UMatInd
                        ! in physical.
                        UMAT(UMatInd(a,b,g,d)) = FourIndInts(i,j,k,l)
                    end do
                end do
            end do
        end do

        ! Also calculate the 2 index integrals, and make these the elements
        ! of the TMAT2D matrix. TMAT2D is in spin orbitals.

        do a = 1, nBasis
            do k = 1, NoOrbs
                i = sym_list(k)
                NewTMAT = 0.0_dp
                do b = 1,NoOrbs
                    d = sym_list(b)
                    if (tStoreSpinOrbs) then
                        NewTMAT = NewTMAT + (one_rdm(b,k)*real(TMAT2D(d,a),dp))
                    else
                        NewTMAT = NewTMAT + (one_rdm(b,k)*real(TMAT2D(2*d,a),dp))
                    end if
                end do
                if (tStoreSpinOrbs) then
                    TMAT2DPart(i,a) = NewTMAT
                else
                    if (mod(a,2).eq.0) then
                        TMAT2DPart(2*i,a) = NewTMAT
                    else
                        TMAT2DPart((2*i)-1,a) = NewTMAT
                    end if
                end if
            end do
        end do

        do k = 1, nBasis
            do l = 1, NoOrbs
                j = sym_list(l)
                NewTMAT = 0.0_dp
                do a = 1, NoOrbs
                    c = sym_list(a)
                    if (tStoreSpinOrbs) then
                        NewTMAT = NewTMAT+(one_rdm(a,l)*TMAT2DPart(k,c))
                    else
                        NewTMAT = NewTMAT+(one_rdm(a,l)*TMAT2DPart(k,2*c))
                    end if
                end do
                if (tStoreSpinOrbs) then
                    TMAT2D(k,j) = NewTMAT
                else
                    if (mod(k,2) .eq. 0) then
                        TMAT2D(k,2*j) = NewTMAT
                    else
                        TMAT2D(k,(2*j)-1) = NewTMAT
                    end if
                end if
            end do
        end do

        deallocate(TMAT2DPart)
        call LogMemDeAlloc('RefillUMAT_RDM',TMAT2DPartTag)

        if (.not. tRotatedNOs) then
            call PrintROFCIDUMP_RDM("ROFCIDUMP")
        end if

    end subroutine RefillUMATandTMAT2D_RDM

    subroutine PrintROFCIDUMP_RDM(filename)

        ! This prints out a new FCIDUMP file in the same format as the old one.

        use IntegralsData, only: umat
        use LoggingData, only: tBrokenSymNOs
        use OneEInts, only: TMAT2D
        use rdm_data, only: tRotatedNOs
        use SystemData, only: nel, G1, tFixLz, tStoreSpinOrbs, lms, ARR, ecore
        use UMatCache, only: UMatInd
        use util_mod, only: get_free_unit

        integer :: i, j, k, l, iunit, orb
        character(len=9) :: filename

!        PrintROFCIDUMP_Time%timer_name='PrintROFCIDUMP'
!        call set_timer(PrintROFCIDUMP_Time,30)

        iunit = get_free_unit()
        open(iunit,file=filename,status='unknown') !'ROFCIDUMP',status='unknown')
        
        write(iunit,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=', NoOrbs, ',NELEC=', NEl, ',MS2=', LMS, ','
        write(iunit,'(A9)',advance='no') 'ORBSYM='
        do i=1,NoOrbs
            if (tStoreSpinOrbs) then
                write(iunit,'(I1,A1)',advance='no') int(G1(i)%sym%S)+1,','
            else
                if (tRotatedNOs .and. tBrokenSymNOs) then
                    write(iunit,'(I1,A1)',advance='no') 1,','
                else
                    write(iunit,'(I1,A1)',advance='no') int(G1(i*2)%sym%S)+1,','
                end if
            end if
        end do

        write(iunit,*) ''

        if (tStoreSpinOrbs) then
            write(iunit,'(A7,I1,A11)') 'ISYM=',1,' UHF=.TRUE.'
        else
            write(iunit,'(A7,I1,A12)') 'ISYM=',1,' UHF=.FALSE.'
        end if

        if (tFixLz) then
            write(iunit,'(A7)',advance='no') 'SYML='
            do i = 1, NoOrbs
                if (i .eq. NoOrbs) then
                    write(iunit,'(I3,A1)') -20,','
                else
                    write(iunit,'(I3,A1)',advance='no') -20,','
                end if
            end do
            write(iunit,'(A8)',advance='no') 'SYMLZ='
            do i = 1, NoOrbs
                orb = i
                if (.not. tStoreSpinOrbs) orb = 2 * orb
                write(iunit, '(i3,",")', advance='no') int(g1(orb)%ml)
            end do
            write(iunit,*)
        end if

        write(iunit,'(A5)') '&end'
       
        do i = 1, NoOrbs
            do j = 1, NoOrbs
                do l = 1, j
                    ! Potential to put symmetry in here, have currently taken it out, 
                    ! because when we're only printing non-zero values, it is kind 
                    ! of unnecessary - although it may be used to speed things up.
                    do k = 1, i
                        ! UMatInd is in physical notation <ij|kl>, but the indices
                        ! printed in the FCIDUMP are in chemical notation (ik|jl).
                        if ((abs(real(UMat(UMatInd(i,j,k,l)),dp))).ne.0.0_dp) &
                                write(iunit,'(F21.12,4I5)') &
                                real(UMat(UMatInd(i,j,k,l)),dp), i, k, j, l 
 
                    end do
                end do
           end do
        end do

        ! TMAT2D stored as spin orbitals.
        do i=1,NoOrbs
            ! Symmetry?
            do k=1,NoOrbs
                if (tStoreSpinOrbs) then
                    if ((real(TMAT2D(i,k),dp)).ne.0.0_dp) write(iunit,'(F21.12,4I5)') &
                                                        real(TMAT2D(i,k),dp),i,k,0,0
                else
                    if ((real(TMAT2D(2*i,2*k),dp)).ne.0.0_dp) write(iunit,'(F21.12,4I5)') &
                                                        real(TMAT2D(2*i,2*k),dp),i,k,0,0
                end if
            end do
        end do

        ! ARR has the energies of the orbitals (eigenvalues).
        ! ARR(:,2) has ordering we want.
        ! ARR is stored as spin orbitals.

        do k=1,NoOrbs
            if (tStoreSpinOrbs) then
                write(iunit,'(F21.12,4I5)') Arr(k,2), k, 0, 0, 0
            else
                write(iunit,'(F21.12,4I5)') Arr(2*k,2), k, 0, 0, 0
            end if
        end do

        write(iunit,'(F21.12,4I5)') ECore, 0, 0, 0, 0
        
        call neci_flush(iunit)

        close(iunit)

!        call halt_timer(PrintROFCIDUMP_Time)

    end subroutine PrintROFCIDUMP_RDM

    subroutine BrokenSymNO(evalues, occ_numb_diff)

        ! This rouine finds natural orbitals (NOs) whose occupation
        ! numbers differ by a small relative threshold (occ_numb_diff) and
        ! rotates them by calling the Rotate2Orbs routine in order to
        ! break symmetry and maximally localise the NOs
        ! We'd like to keep both the original NOs (Natural Orbitals) 
        ! and the broken maximally localised NOs for checking.
        ! This is not very (time-) efficient at the moment.

        use IntegralsData, only: umat
        use LoggingData, only: tBreakSymNOs, local_cutoff, tagRotNOs, RotNOs
        use LoggingData, only: rottwo, rotthree, rotfour
        use MemoryManager, only: LogMemDealloc
        use Parallel_neci, only: iProcIndex
        use rdm_data, only: tRotatedNOs
        use SystemData, only: nel, tStoreSpinOrbs
        use UMatCache, only: UMatInd

        real(dp), intent(in) :: evalues(:)
        real(dp), intent(in) :: occ_numb_diff

        real(dp) :: diffnorm, SumDiag, sum_old, sum_new, selfint_old
        real(dp), allocatable :: one_rdm(:,:), trans_2orbs_coeffs(:,:)
        real(dp), allocatable :: selfint(:)
        integer, allocatable :: rotate_list(:,:), rotorbs(:,:)
        integer, allocatable :: sym_list(:)
        integer :: l1, l2, l3, l4, l5, l6, m, n
        integer :: iumat,jumat
        logical :: partnerfound, localdelocal

        allocate(one_rdm(NoOrbs,NoOrbs))
        allocate(sym_list(NoOrbs))
        allocate(trans_2orbs_coeffs(2,2))
        allocate(rotorbs(6,2))
        allocate(rotate_list((NoOrbs*(NoOrbs-1)),4))
        allocate(selfint(NoOrbs))

        if (iProcIndex .eq. 0) then

            ! No symmetry ordering is applied.
            do l1 = 1, NoOrbs
                sym_list(l1) = l1
            end do

            ! Normalisation.
            SumDiag = sum(evalues)

            if (tStoreSpinOrbs) then
                diffnorm = SumDiag/dble(NEl)
            else
                diffnorm = 2.0_dp*(SumDiag/dble(NEl))
            end if
            
            trotatedNOs = .true.
     
            if (tStoreSpinorbs) then
                call Stop_all("BrokenSymNO", "Broken symmetry NOs currently not implemented for UHF")
            end if

            write(6,*) '------------------------------------------------------------------------------'
            write(6,*) 'Localising NOs whose occupation numbers differ by less than threshold'
            write(6,*) '------------------------------------------------------------------------------'

            if (tBreakSymNOs) then
                write(6,*) 'Rotating specified NOs'
            else
                write(6,*) 'Threshold for orbitals to rotate:', occ_numb_diff
            end if

            ! Self-interactions.
            selfint(:) = 0.0_dp
            do l1 = 1, NoOrbs
                selfint(l1) = Umat(UmatInd(l1,l1,l1,l1))
            end do

            write(6,*) 'Self-interactions for NOs:'
            do l1 = 1, NoOrbs
                write(6,'(I3,3X,G25.12)') l1, selfint(l1)
            end do
            write(6,*) 'Sum of NO selfinteractions:',sum(selfint)
            selfint_old = sum(selfint)

            ! If the NOs to be rotated are specified in the input file.
            if (tBreakSymNOs) then
                rotate_list(:,:) = 0
                rotorbs(:,:) = 0
                m = 0
                do l1 = 2,(2*rottwo),2
                    m = m + 1
                    rotate_list(m,1) = RotNOs(l1-1)
                    rotate_list(m,2) = RotNOs(l1)
                end do
                do l1 = ((2*rottwo)+3),((2*rottwo)+(3*rotthree)),3
                    m = m + 1
                    rotate_list(m,1) = RotNOs(l1-2)
                    rotate_list(m,2) = RotNOs(l1-1)
                    rotate_list(m,3) = RotNOs(l1)
                end do
                do l1 = ((2*rottwo)+(3*rotthree)+4),((2*rottwo)+(3*rotthree)+(4*rotfour)),4
                    m = m + 1
                    rotate_list(m,1) = RotNOs(l1-3)
                    rotate_list(m,2) = RotNOs(l1-2)
                    rotate_list(m,3) = RotNOs(l1-1)
                    rotate_list(m,4) = RotNOs(l1)
                end do
            else
                ! If the threshold is used to generate a list of NOs to be
                ! rotated.

                ! Generate the list of orbitals which are rotated amongst each
                ! other.
                rotate_list(:,:) = 0
                rotorbs(:,:) = 0

                ! Need to account for spatial and spin orbital representations
                ! since orbitals of different spin cannot be mixed.
                ! List contains the NOs which are rotated.
                ! It can deal with a maximum of four NOs which are mixed.
                m = 0
                n = 1
                do l1 = 1, NoOrbs
                    if ((m .ne. 0) .and. (l1 .le. rotate_list(m,n))) cycle
                        partnerfound = .false.
                        n = 1
                        do l2 = (l1+1), NoOrbs

                            if ((abs((evalues(l1)/diffnorm)-(evalues(l2)/diffnorm))/abs((evalues(l2)/diffnorm)))&
                                & .lt. occ_numb_diff) then
                            if (.not. partnerfound) then
                                m = m + 1
                                n = n + 1
                                rotate_list(m,1) = l1
                                rotate_list(m,2) = l2
                                partnerfound = .true.
                            else if (partnerfound) then
                                n = n + 1
                                ! this is for up to 2-fold degenearcy
!                                if (n.gt.2) then
!                                    n = 2
!                                    write(6,*) '***Warning***'
!                                    write(6,*) 'Threshold generated more than 2-fold degeneracy'
!                                    write(6,*) 'NOs around:',l2
!                                    cycle  ! don't want to rotate more than 2 orbitals
!                                end if
                                ! this is for up to four-fold degeneracy
                                if (n.gt.4) then
                                    n = 4
                                    write(6,*) '***Warning***'
                                    write(6,*) 'Threshold generated more than 4-fold degeneracy'
                                    write(6,*) 'NOs around:',l2
                                    cycle  ! don't want to rotate more than 4 orbitals
                                end if
                                rotate_list(m,n) = l2
                            end if
                        end if
                    end do
                end do
            end if

            write(6,*) 'The following pairs of orbitals will be rotated:'
            do l1 = 1, m
                write(6,'(I3,3X,4(I3))') l1, rotate_list(l1,:)
            end do

            one_rdm(:,:) = 0.0_dp
            do l1 = 1, NoOrbs
                one_rdm(l1,l1) = 1.0_dp
            end do

            ! Rotate two-fold degenerate pairs first.
            do l1 = 1,m
                ! If only two orbitals have the same occupation numbers.
                if (rotate_list(l1,3).eq.0) then
                    write(6,'(A20,4(I3))') 'Rotating NOs:', rotate_list(l1,:)
                    iumat = rotate_list(l1,1)
                    jumat = rotate_list(l1,2)
                    if (jumat.le.local_cutoff) then
                        localdelocal = .false.
                    else if (jumat.gt.local_cutoff) then
                        localdelocal = .true.
                    end if
                    call Rotate2Orbs(iumat,jumat,trans_2orbs_coeffs,selfint(iumat),&
                        &selfint(jumat),localdelocal)

                    ! The new NOs are 
                    ! phi_{i'} = cos a p_{i} + sin a p_{j}
                    ! phi_{j'} = -sin a p_{i} + cos a p_{j}
                    one_rdm(iumat,iumat) = trans_2orbs_coeffs(1,1)
                    one_rdm(jumat,iumat) = trans_2orbs_coeffs(2,1)
                    one_rdm(iumat,jumat) = trans_2orbs_coeffs(1,2)
                    one_rdm(jumat,jumat) = trans_2orbs_coeffs(2,2)

                    write(6,*) 'Sum of rotated NO self-interactions:',sum(selfint)

                    selfint_old = sum(selfint)
                end if
           end do

           ! Transform integral corresponding to rotated NO.
           ! These are required for not printing out RPFCIDUMP or 
           ! BSFCIDUMP every time.
           trotatedNOs = .true.
           call Transform2ElIntsMemSave_RDM(one_rdm, sym_list)
           call RefillUMATandTMAT2D_RDM(one_rdm, sym_list)

           ! If three orbitals are degenerate.
           do l1 = 1, m
                 if ((rotate_list(l1,3) .ne. 0) .and. (rotate_list(l1,4) .eq. 0)) then

                        sum_new = sum(selfint)
                        rotorbs(1,1) = 1
                        rotorbs(1,2) = 2
                        rotorbs(2,1) = 1
                        rotorbs(2,2) = 3
                        rotorbs(3,1) = 2
                        rotorbs(3,2) = 3

                        ! These have to be done self-consistently since all
                        ! three orbitals can intermix.
                        do
                            sum_old = sum_new
                            write(6,'(A20,4(I3))') 'Rotating NOs:', rotate_list(l1,:)
                            do l3 = 1,3
                                iumat = rotate_list(l1,rotorbs(l3,1))
                                jumat = rotate_list(l1,rotorbs(l3,2))
                                one_rdm = 0.0_dp
                                do l4=1,NoOrbs
                                    if ((l4 .ne. iumat).and.(l4 .ne. jumat)) then
                                        one_rdm = 1.0_dp
                                    end if
                                end do
                                if (jumat .le. local_cutoff) then
                                    localdelocal = .false.
                                else if (jumat .gt. local_cutoff) then
                                    localdelocal = .true.
                                end if
                                call Rotate2Orbs(iumat,jumat,&
                                    &trans_2orbs_coeffs,selfint(iumat),selfint(jumat)&
                                    &,localdelocal)

                                ! The new NOs are 
                                ! phi_{i'} = cos a p_{i} + sin a p_{j}
                                ! phi_{j'} = -sin a p_{i} + cos a p_{j}
                                one_rdm(iumat,iumat) = trans_2orbs_coeffs(1,1)
                                one_rdm(jumat,iumat) = trans_2orbs_coeffs(2,1)
                                one_rdm(iumat,jumat) = trans_2orbs_coeffs(1,2)
                                one_rdm(jumat,jumat) = trans_2orbs_coeffs(2,2)

                                ! Transform integral corresponding to rotated NO.
                                ! These are required for not printing out
                                ! RPFCIDUMP or BSFCIDUMP every time.
                                trotatedNOs = .true.
                                call Transform2ElIntsMemSave_RDM(one_rdm, sym_list)
                                call RefillUMATandTMAT2D_RDM(one_rdm, sym_list)
                            end do

                            ! Check for convergence.
                            sum_new = sum(selfint)
                            write(6,'(A50,2G20.12)') 'Current and previous selfinteraction:',&
                                &sum_new,sum_old
                            if (abs(sum_new-sum_old).lt.1e-12_dp) then
                                exit
                            end if
                        end do

                    write(6,*) 'Sum of rotated NO self-interactions:', sum(selfint)

                    selfint_old = sum(selfint)

                else if ((rotate_list(l1,3).ne.0).and.(rotate_list(l1,4).ne.0)) then

                        sum_new = sum(selfint)
                        rotorbs(1,1) = 1
                        rotorbs(1,2) = 2
                        rotorbs(2,1) = 3
                        rotorbs(2,2) = 4
                        rotorbs(3,1) = 1
                        rotorbs(3,2) = 3
                        rotorbs(4,1) = 2
                        rotorbs(4,2) = 4
                        rotorbs(5,1) = 1
                        rotorbs(5,2) = 4
                        rotorbs(6,1) = 2
                        rotorbs(6,2) = 3

                        ! These have to be done self-consistently since all
                        ! three orbitals can intermix.
                        do
                            sum_old = sum_new
                            write(6,'(A20,4(I3))') 'Rotating NOs:',rotate_list(l1,:)
                            do l3 = 1, 3
                                one_rdm(:,:) = 0.0_dp
                                do l4 = 1, NoOrbs
                                    if ((l4 .ne. rotate_list(l1,1)) .and. (l4 .ne. rotate_list(l1,2))&
                                        & .and. (l4.ne.rotate_list(l1,3)) .and. (l4 .ne. rotate_list(l1,4))) then
                                            one_rdm(l4,l4) = 1.0_dp
                                    end if
                                end do
                                ! Rotate these two independently.
                                do l5 = 0, 1
                                    iumat = rotate_list(l1,rotorbs(((2*l3)-l5),1))
                                    jumat = rotate_list(l1,rotorbs(((2*l3)-l5),2))
                                    if (jumat.le.local_cutoff) then
                                        localdelocal = .false.
                                    else if (jumat.gt.local_cutoff) then
                                        localdelocal = .true.
                                    end if
                                    call Rotate2Orbs(iumat, jumat, trans_2orbs_coeffs, &
                                                     selfint(iumat), selfint(jumat), localdelocal)

                                    ! The new NOs are 
                                    ! phi_{i'} = cos a p_{i} + sin a p_{j}
                                    ! phi_{j'} = -sin a p_{i} + cos a p_{j}
                                    one_rdm(iumat,iumat) = trans_2orbs_coeffs(1,1)
                                    one_rdm(jumat,iumat) = trans_2orbs_coeffs(2,1)
                                    one_rdm(iumat,jumat) = trans_2orbs_coeffs(1,2)
                                    one_rdm(jumat,jumat) = trans_2orbs_coeffs(2,2)
                                end do

                                ! Transform integral corresponding to rotated NO.
                                ! These are required for not printing out
                                ! RPFCIDUMP or BSFCIDUMP every time.
                                trotatedNOs = .true.
                                call Transform2ElIntsMemSave_RDM(one_rdm, sym_list)
                                call RefillUMATandTMAT2D_RDM(one_rdm, sym_list)
                            end do
                            ! Check for convergence.
                            sum_new = sum(selfint)

                            write(6,"(A50,2G20.12)") 'Current and previous selfinteractions:',&
                                &sum_new,sum_old

                            if (abs(sum_new-sum_old) .lt. 1e-12_dp) then
                                exit
                            end if
                        end do

                    write(6,*) 'Sum of rotated NO self-interactions:',sum(selfint)

                    selfint_old = sum(selfint)
                end if
           end do

            write(6,*) 'Final self-interactions for rotated NOs:'
            do l1 = 1,NoOrbs
                write(6,'(I3,3X,G25.12)') l1, selfint(l1)
            end do
            write(6,*) 'Sum of rotated NO self-interactions:', sum(selfint)

            write(6,*) '------------------------------------------------------'
            write(6,*) 'Writing out BSFCIDUMP...'
            write(6,*) '------------------------------------------------------'
            
            call PrintROFCIDUMP_RDM("BSFCIDUMP")
 
        end if

        deallocate(one_rdm)
        deallocate(sym_list)
        deallocate(trans_2orbs_coeffs)
        deallocate(rotate_list)
        deallocate(rotorbs)
        deallocate(selfint)
        
        if (tBreakSymNOs) then
            deallocate(RotNOs)
            call LogMemDealloc('BrokenSymNO',tagRotNOs)
        end if

    end subroutine BrokenSymNO

    subroutine Rotate2Orbs(i, j, trans_2orbs_coeffs, selfintorb1, selfintorb2, localdelocal)

        ! This routine takes two orbitals i,j, and rotates them in order to
        ! maximally localise these. It employs an Edminston-Ruedenberg type
        ! localisation which maximises the self-interaction
        ! \sum_{i=1}^{2} \sum_{r,s,u,v} (c_{ir})*(c_{is})*c_{iu}c_{iv} <p_{i}p_{i}|u|p_{i}p_{i}>
        ! where p_{i} are the original NOs.
        ! The coefficients c are given by the following matrix:
        ! c =  cos a  sin a
        !     -sin a  cos a
        ! Then angle a is found by differentiating and setting it equal to 0
        ! which gives the following analytical expression of the form
        ! tan a = -x/y
        ! where x and y are sums of the original NO four index integrals.

        use IntegralsData, only: umat
        use UMatCache, only: UMatInd

        real(dp), allocatable, intent(inout) :: trans_2orbs_coeffs(:,:)
        real(dp), intent(inout) :: selfintorb1,selfintorb2

        real(dp) :: alpha2(17)
        real(dp) :: selfinteractions(17)
        real(dp) :: coeffcos,coeffsin,maxint
        integer :: maxangle(1)
        integer :: indicesij(2)
        integer, intent(in) :: i,j
        integer :: l1,l2,l3,l4,l5
        logical, intent(in) :: localdelocal

        indicesij(1) = i
        indicesij(2) = j
        trans_2orbs_coeffs(:,:) = 0.0_dp

        ! Umat(UMatInd(i,j,k,l)) contains the four-index integrals
        ! <ij|kl> (physical notation) in the NO basis

        coeffcos = Umat(UmatInd(i,i,i,j)) + Umat(UmatInd(i,i,j,i)) + Umat(UmatInd(i,j,i,i)) &
            & - Umat(UmatInd(i,j,j,j)) + Umat(UmatInd(j,i,i,i)) - Umat(UmatInd(j,i,j,j)) &
            & - Umat(UmatInd(j,j,i,j)) - Umat(UmatInd(j,j,j,i))

        coeffsin = -Umat(UmatInd(i,i,i,i)) + Umat(UmatInd(i,i,j,j)) + Umat(UmatInd(i,j,i,j)) &
            & + Umat(UmatInd(i,j,j,i)) + Umat(UmatInd(j,i,i,j)) + Umat(UmatInd(j,i,j,i)) &
            & + Umat(UmatInd(j,j,i,i)) - Umat(UmatInd(j,j,j,j))

        ! atan return a value in [-pi/2,pi/2]
        ! because of the 4*alpha in the equation there are 8 distinct solutions
        ! i.e. in the range 0,2*pi
        ! i.e. possible solutions are separated by (2*pi/8)=pi/4
        ! for safety 16 solutions are evaluated.
        alpha2(9) = atan((-coeffcos/coeffsin))
        alpha2(9) = alpha2(9)/4.0_dp
        do l1 = 8, 1, -1
            alpha2(l1) = alpha2(l1+1) - (pi/4.0_dp)
        end do
        do l1 = 10, 17
            alpha2(l1) = alpha2(l1-1) + (pi/4.0_dp)
        end do
        
        ! second derivatives to find maximum (necessary since the minimum, i.e. fully delocalised
        ! orbitals satisfy the same conditions
        !secondderiv(1) = (4.0_dp*coeffsin*cos(4.0_dp*alpha(1))) - (4.0_dp*coeffcos*sin(4.0_dp*alpha(1)))
        !secondderiv(2) = (4.0_dp*coeffsin*cos(4.0_dp*alpha(2))) - (4.0_dp*coeffcos*sin(4.0_dp*alpha(2)))

        ! Compute selfinteractions to check which one is largest.
        ! This is a better measure than the second derivatives.
        selfinteractions(:) = 0.0_dp 

        do l1 = 1, 17
            trans_2orbs_coeffs(1,1) = cos(alpha2(l1))
            trans_2orbs_coeffs(2,1) = sin(alpha2(l1))
            trans_2orbs_coeffs(1,2) = -sin(alpha2(l1))
            trans_2orbs_coeffs(2,2) = cos(alpha2(l1))

            do l2 = 1, 2
                do l3 = 1, 2
                    do l4 = 1, 2
                        do l5 = 1, 2
                            selfinteractions(l1) = selfinteractions(l1) + trans_2orbs_coeffs(l2,1)&
                                &*trans_2orbs_coeffs(l3,1)*&
                                &trans_2orbs_coeffs(l4,1)*trans_2orbs_coeffs(l5,1)*&
                                &Umat(UmatInd(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5)))
                        end do
                    end do
                end do
            end do
            do l2 = 1, 2
                do l3 = 1, 2
                    do l4 = 1, 2
                        do l5 = 1, 2
                            selfinteractions(l1) = selfinteractions(l1) + trans_2orbs_coeffs(l2,2)&
                                &*trans_2orbs_coeffs(l3,2)*&
                                &trans_2orbs_coeffs(l4,2)*trans_2orbs_coeffs(l5,2)*&
                                &Umat(UmatInd(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5)))
                        end do
                    end do
                end do
            end do
        end do

        ! Choose the angle which maximises the self interactions.
        if (.not.localdelocal) then
            ! Maximally delocalised.
            maxangle = minloc(selfinteractions)
            maxint = minval(selfinteractions)
        else if (localdelocal) then
            ! Maximally localised.
            maxangle = maxloc(selfinteractions)
            maxint = maxval(selfinteractions)
        end if


        ! Return transformatin coefficients.
        trans_2orbs_coeffs(1,1) =  cos(alpha2(maxangle(1)))
        trans_2orbs_coeffs(2,1) =  sin(alpha2(maxangle(1)))
        trans_2orbs_coeffs(1,2) = -sin(alpha2(maxangle(1)))
        trans_2orbs_coeffs(2,2) =  cos(alpha2(maxangle(1)))
 
        ! New self-interactions for transformed orbitals.
        selfintorb1 = 0.0_dp
        selfintorb2 = 0.0_dp

        do l2 = 1, 2
            do l3 = 1, 2
                do l4 = 1, 2
                    do l5 = 1, 2
                        selfintorb1 = selfintorb1 + trans_2orbs_coeffs(l2,1)&
                            &*trans_2orbs_coeffs(l3,1)*&
                            &trans_2orbs_coeffs(l4,1)*trans_2orbs_coeffs(l5,1)*&
                            &Umat(UmatInd(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5)))

                        selfintorb2 = selfintorb2 + trans_2orbs_coeffs(l2,2)&
                            &*trans_2orbs_coeffs(l3,2)*&
                            &trans_2orbs_coeffs(l4,2)*trans_2orbs_coeffs(l5,2)*&
                            &Umat(UmatInd(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5)))
                    end do
                end do
            end do
        end do

    end subroutine Rotate2Orbs

end module rdm_nat_orbs
