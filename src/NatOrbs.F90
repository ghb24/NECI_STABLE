! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
MODULE NatOrbsMod
! This file is primarily concerned with finding the one electron reduced density matrix, from the wavefunction
! constructed by a previous spawning calculation.
! Diagonalisation of this density matrix gives a set of eigenvectors which rotate the HF orbitals into the 
! CI natural orbitals within the given excitation level.
! Once these eigenvectors have been obtained, the relevant routines from RotateOrbs are called to transform the 
! integrals and produce a ROFCIDUMP file in the natural orbital basis.
        
        USE Global_utilities
        USE Parallel_neci
        USE IntegralsData , only : UMAT
        USE UMatCache , only : UMatInd
        USE SystemData , only : NEl,nBasis,G1,ARR,BRR,lNoSymmetry,LMS,tStoreSpinOrbs,nOccAlpha,&
                                nOccBeta,tSeparateOccVirt, tRotateOccOnly,tRotateVirtOnly,&
                                tFindCINatOrbs,tUseMP2VarDenMat,nBasisMax,ALAT,iSpinSkip
        use bit_reps, only: NIfY, NIfTot
        USE RotateOrbsData , only : SymLabelList2_rot,SymLabelCounts2_rot,SymLabelCounts2_rotTag,&
                                    SymLabelListInv_rot,NoOrbs,SpatOrbs,FillOneRDM_time,&
                                    FillMP2VDM_Time,DiagNatOrbMat_Time,OrderCoeff_Time,FillCoeff_Time,&
                                    NoFrozenVirt,SymLabelList3_rot
        use sort_mod
        use bit_reps, only: decode_bit_det
        use MemoryManager, only: TagIntType
        use util_mod, only: get_free_unit
        IMPLICIT NONE
        INTEGER(TagIntType) :: NoSpinCyc,SymOrbs_rotTempTag
        real(dp) , ALLOCATABLE :: NatOrbMat(:,:),Evalues(:)
        INTEGER , ALLOCATABLE :: SymOrbs_rotTemp(:)
        INTEGER ierr
        INTEGER(TagIntType) :: NatOrbMatTag,EvaluesTag

    contains
    
    SUBROUTINE FindNatOrbs()
        IMPLICIT NONE
        INTEGER :: ierr
        CHARACTER(len=*) , PARAMETER :: this_routine='FindNatOrbs'

! Fed into this routine will be the wavefunction, Psi, and its amplitudes within the given excitation level.    

! First need to set up the orbital labels and symmetries etc.
! This is done slightly differently for spin and spatial and whether or not we are truncating the virtual 
! space when writing out the final ROFCIDUMP file.

! Allocate the matrix used to find the natural orbitals.
        
        ALLOCATE(NatOrbMat(NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('NatOrbMat',NoOrbs**2,8,this_routine,NatOrbMatTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for NatOrbMat failed.")
        NatOrbMat(:,:)=0.0_dp

        ALLOCATE(Evalues(NoOrbs),stat=ierr)
        CALL LogMemAlloc('Evalues',NoOrbs,8,this_routine,EvaluesTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for Evalues failed.")
        Evalues(:)=0.0_dp

! First need to fill the relevant matrix for calculating the type of natural orbitals we want.
        IF(tFindCINatOrbs) THEN

            ! For the CISD, CISDT etc natural orbitals, the relevant matrix is the one electron reduced 
            ! density matrix from the previous spawning calculation (trucated at a certain excitation).
            CALL FillOneRDM()

        ELSEIF(tUseMP2VarDenMat) THEN

            ! For the MP2 natural orbitals, the natural orbital matrix is the MP2 variational density matrix.
            CALL FillMP2VDM()

        ENDIF

! Then need to diagonalise this, maintaining the various symmetries (spin and spatial).    
        CALL DiagNatOrbMat()

! We then need to put the resulting eigenvectors back into the ordering we want, and copy these over to CoeffT1.
        CALL OrderCoeffT1()

        
    END SUBROUTINE FindNatOrbs



    SUBROUTINE SetupNatOrbLabels()
        use MemoryManager, only: TagIntType
        IMPLICIT NONE
        INTEGER :: x,i,j,ierr,NoOcc,StartFill01,StartFill02,Symi,SymCurr,Prev,EndFill01,EndFill02
        CHARACTER(len=*) , PARAMETER :: this_routine='SetupNatOrbLabels'
        INTEGER , ALLOCATABLE :: LabVirtOrbs(:),LabOccOrbs(:),SymVirtOrbs(:),SymOccOrbs(:)
        INTEGER(TagIntType) :: LabVirtOrbsTag,LabOccOrbsTag,SymVirtOrbsTag,SymOccOrbsTag
        integer :: lo, hi
        

! The earlier test should pick this up, if it crashes here, will want to put in an earlier test so that we don't 
! get all the way to this stage.
        IF((LMS.ne.0).and.(.not.tStoreSpinOrbs)) CALL Stop_All("FindNatOrbs","Open shell system, and UMAT is &
                                                                &not being stored as spin orbitals.")

! We now need two slightly different sets of orbital labels for the case of spin orbitals and spatial orbitals. 
! When using spin orbitals we want all the beta spin followed by all the alpha spin. 
! Then we want two values for the number of occupied orbitals to allow for high spin cases.
! With spatial, it is equivalent to just keeping the beta spin.
        IF(tStoreSpinOrbs) THEN
            NoSpinCyc=2
        ELSE
            NoSpinCyc=1
        ENDIF

        do x=1,NoSpinCyc
            IF(.not.tSeparateOccVirt) THEN
                NoOcc=0
            ELSE
                IF(x.eq.1) THEN
                    IF(tStoreSpinOrbs) THEN
                        NoOcc=nOccBeta
                    ELSE
                        NoOcc=NEl/2
                    ENDIF
                ENDIF
                IF(x.eq.2) NoOcc=nOccAlpha
            ENDIF

            IF(tSeparateOccVirt) THEN
                ALLOCATE(LabOccOrbs(NoOcc),stat=ierr)
                CALL LogMemAlloc('LabOccOrbs',(NoOcc),4,this_routine,LabOccOrbsTag,ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for LabOccOrbs failed.")
                LabOccOrbs(:)=0
                ALLOCATE(SymOccOrbs(NoOcc),stat=ierr)
                CALL LogMemAlloc('SymOccOrbs',(NoOcc),4,this_routine,SymOccOrbsTag,ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for SymOccOrbs failed.")
                SymOccOrbs(:)=0
            ENDIF

            ALLOCATE(LabVirtOrbs(SpatOrbs-NoOcc),stat=ierr)
            CALL LogMemAlloc('LabVirtOrbs',(SpatOrbs-NoOcc),4,this_routine,LabVirtOrbsTag,ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for LabVirtOrbs failed.")
            LabVirtOrbs(:)=0
            ALLOCATE(SymVirtOrbs(SpatOrbs-NoOcc),stat=ierr)
            CALL LogMemAlloc('SymVirtOrbs',(SpatOrbs-NoOcc),4,this_routine,SymVirtOrbsTag,ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for SymVirtOrbs failed.")
            SymVirtOrbs(:)=0

! First fill SymLabelList2_rot.

! Brr has the orbital numbers in order of energy... i.e Brr(2) = the orbital index with the second lowest energy.

!            do i=1,(SpatOrbs*2)
!                WRITE(6,*) BRR(i)
!            enddo
!            CALL neci_flush(6)
!            CALL Stop_All('','')

! this picks out the NoOcc lowest energy orbitals from BRR as these will be the occupied.
! these are then ordered according to symmetry, and the same done to the virtual.
            do i=1,NoOcc
                IF(x.eq.1) THEN
                    IF(tStoreSpinOrbs) THEN
                        LabOccOrbs(i)=BRR(2*i)-1
                        SymOccOrbs(i)=INT(G1(LabOccOrbs(i))%sym%S,4)
                    ELSE
                        LabOccOrbs(i)=BRR(2*i)/2
                        SymOccOrbs(i)=INT(G1(LabOccOrbs(i)*2)%sym%S,4)
                    ENDIF
                ELSEIF(x.eq.2) THEN
                    LabOccOrbs(i)=BRR(2*i)
                    SymOccOrbs(i)=INT(G1(LabOccOrbs(i))%sym%S,4)
                ENDIF
            enddo
            IF(tSeparateOccVirt) call sort (SymOccOrbs, LabOccOrbs)
            ! Sorts LabOrbs according to the order of SymOccOrbs (i.e. in terms of symmetry). 


            ! Same for the virtual
            do i=1,SpatOrbs-NoOcc
                IF(x.eq.1) THEN
                    IF(tStoreSpinOrbs) THEN
                        IF(tSeparateOccVirt) THEN
                            LabVirtOrbs(i)=BRR((2*i)+(NoOcc*2))-1
                            SymVirtOrbs(i)=INT(G1(LabVirtOrbs(i))%sym%S,4)
                        ELSE
                            LabVirtOrbs(i)=BRR((2*i))-1
                            SymVirtOrbs(i)=INT(G1(LabVirtOrbs(i))%sym%S,4)
                        ENDIF
                    ELSE
                        IF(tSeparateOccVirt) THEN
                            LabVirtOrbs(i)=BRR((2*i)+(NoOcc*2))/2
                            SymVirtOrbs(i)=INT(G1(LabVirtOrbs(i)*2)%sym%S,4)
                        ELSE
                            LabVirtOrbs(i)=BRR((2*i))/2
                            SymVirtOrbs(i)=INT(G1(LabVirtOrbs(i)*2)%sym%S,4)
                        ENDIF
                    ENDIF
                ELSEIF(x.eq.2) THEN
                    IF(tSeparateOccVirt) THEN
                        LabVirtOrbs(i)=BRR((2*i)+(NoOcc*2))
                        SymVirtOrbs(i)=INT(G1(LabVirtOrbs(i))%sym%S,4)
                    ELSE
                        LabVirtOrbs(i)=BRR((2*i))
                        SymVirtOrbs(i)=INT(G1(LabVirtOrbs(i))%sym%S,4)
                    ENDIF
                ENDIF
            enddo

!            WRITE(6,*) 'before sort'
!            do i=1,NoOrbs
!                WRITE(6,*) LabVirtOrbs(i),SymVirtOrbs(i)
!            enddo

            call sort (SymVirtOrbs, LabVirtOrbs)
            ! Sorts LabOrbs according to the order of SymOccOrbs (i.e. in terms of symmetry). 

!            WRITE(6,*) 'after sort'
!            do i=1,NoOrbs
!                WRITE(6,*) LabVirtOrbs(i),SymVirtOrbs(i)
!            enddo

!            CALL neci_flush(6)
!            stop

 
! SymLabelList2_rot is then filled with the symmetry ordered occupied then virtual arrays for each spin.        
            IF(x.eq.1) THEN
                StartFill01=1
                StartFill02=NoOcc+1
                EndFill01=NoOcc
                EndFill02=SpatOrbs
            ELSEIF(x.eq.2) THEN
                StartFill01=SpatOrbs+1
                StartFill02=SpatOrbs+NoOcc+1
                EndFill01=SpatOrbs+NoOcc
                EndFill02=NoOrbs
            ENDIF

            j=0
            do i=StartFill01,EndFill01
                j=j+1
                SymLabelList2_rot(i)=LabOccOrbs(j)
            enddo
            j=0
            do i=StartFill02,EndFill02
                j=j+1
                SymLabelList2_rot(i)=LabVirtOrbs(j)
            enddo


!************
! Second fill SymLabelCounts2_rot.
! - the first 8 places of SymLabelCounts2_rot(1,:) and SymLabelCounts2_rot(2,:) refer to the occupied orbitals 
! - and the second 8 to the virtuals.

            IF(lNoSymmetry) THEN
                ! if we are ignoring symmetry, all orbitals essentially have symmetry 0.
                IF(x.eq.1) THEN
                    SymLabelCounts2_rot(1,1)=1
                    SymLabelCounts2_rot(1,9)=NoOcc+1
                    SymLabelCounts2_rot(2,1)=NoOcc
                    SymLabelCounts2_rot(2,9)=SpatOrbs-NoOcc
                ELSEIF(x.eq.2) THEN
                    SymLabelCounts2_rot(1,17)=1
                    SymLabelCounts2_rot(1,25)=NoOcc+1
                    SymLabelCounts2_rot(2,17)=NoOcc
                    SymLabelCounts2_rot(2,25)=SpatOrbs-NoOcc
                ENDIF
                
            ELSE 
                ! otherwise we run through the occupied orbitals, counting the number with each symmetry
                ! and noting where in SymLabelList2_rot each symmetry block starts.
                IF(x.eq.1) THEN
                    StartFill01=1
                    StartFill02=9
                    Prev=0
                ELSEIF(x.eq.2) THEN
                    StartFill01=17
                    StartFill02=25
                    Prev=SpatOrbs
                ENDIF
                SymCurr=0
                SymLabelCounts2_rot(1,StartFill01)=1+Prev
                do i=1,NoOcc
                    IF(tStoreSpinOrbs) THEN
                        Symi=INT(G1(SymLabelList2_rot(i+Prev))%sym%S,4)
                    ELSE
                        Symi=INT(G1((SymLabelList2_rot(i+Prev)*2))%sym%S,4)
                    ENDIF
                    SymLabelCounts2_rot(2,(Symi+StartFill01))=SymLabelCounts2_rot(2,(Symi+StartFill01))+1
                    IF(Symi.ne.SymCurr) THEN
                        SymLabelCounts2_rot(1,(Symi+StartFill01))=i+Prev
                        SymCurr=Symi
                    ENDIF
                enddo
                ! the same is then done for the virtuals.
                SymCurr=0
                SymLabelCounts2_rot(1,StartFill02)=NoOcc+1+Prev
                do i=NoOcc+1,SpatOrbs
                    IF(tStoreSpinOrbs) THEN
                        Symi=INT(G1(SymLabelList2_rot(i+Prev))%sym%S,4)
                    ELSE
                        Symi=INT(G1((SymLabelList2_rot(i+Prev)*2))%sym%S,4)
                    ENDIF
                    SymLabelCounts2_rot(2,(Symi+StartFill02))=SymLabelCounts2_rot(2,(Symi+StartFill02))+1
                    IF(Symi.ne.SymCurr) THEN
                        SymLabelCounts2_rot(1,(Symi+StartFill02))=i+Prev
                        SymCurr=Symi
                    ENDIF
                enddo
            ENDIF
     
            ! Go through each symmetry group, making sure the orbital pairs are ordered lowest to highest.
            IF(x.eq.1) THEN
                do i=1,16
                    IF(SymLabelCounts2_rot(2,i).ne.0) THEN
                        lo = SymLabelCounts2_rot(1, i)
                        hi = lo + SymLabelCounts2_rot(2, i) - 1
                        call sort (SymLabelList2_rot (lo:hi))
                    ENDIF
                enddo
            ELSEIF(x.eq.2) THEN
                do i=17,32
                    IF(SymLabelCounts2_rot(2,i).ne.0) THEN
                        lo = SymLabelCounts2_rot(1, i)
                        hi = lo + SymLabelCounts2_rot(2, i) - 1
                        call sort (SymLabelList2_rot (lo:hi))
                    ENDIF
                enddo
            ENDIF

! Deallocate the arrays just used in this routine.

            IF(tSeparateOccVirt) THEN
                DEALLOCATE(SymOccOrbs)
                CALL LogMemDealloc(this_routine,SymOccOrbsTag)

                DEALLOCATE(LabOccOrbs)
                CALL LogMemDealloc(this_routine,LabOccOrbsTag)
            ENDIF

            DEALLOCATE(SymVirtOrbs)
            CALL LogMemDealloc(this_routine,SymVirtOrbsTag)

            DEALLOCATE(LabVirtOrbs)
            CALL LogMemDealloc(this_routine,LabVirtOrbsTag)

        enddo

        do i=1,NoOrbs
            SymLabelListInv_rot(SymLabelList2_rot(i))=i
        enddo

        do i=1,NoOrbs
            SymLabelList3_rot(i) = SymLabelList2_rot(i)
        enddo

        IF(.not.tSeparateOccVirt) THEN
            ! basically we treat all the orbitals as virtuals and set NoOcc to zero in each routine. 
            tRotateVirtOnly=.true.
        ENDIF
        

!        WRITE(6,*) 'Sym Label Counts'
!        do i=1,32
!            WRITE(6,*) i,SymLabelCounts2_rot(1,i),SymLabelCounts2_rot(2,i)
!        enddo
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and their symmetries according to G1'
!        do i=1,NoOrbs
!            IF(tStoreSpinOrbs) THEN
!                WRITE(6,*) i,SymLabelList2_rot(i),INT(G1(SymLabelList2_rot(i))%sym%S,4)
!            ELSE
!                WRITE(6,*) i,SymLabelList2_rot(i),INT(G1(SymLabelList2_rot(i)*2)%sym%S,4)
!            ENDIF
!        enddo
!        WRITE(6,*) 'i','ARR(SymLabelList2_rot(i),1)','ARR(SymLabelList2_rot(i),2)','Sym'
!        do i=1,NoOrbs
!            IF(tStoreSpinOrbs) THEN
!                WRITE(6,*) i,ARR(SymLabelList2_rot(i),1),ARR(SymLabelList2_rot(i),2),&
!                                    INT(G1(SymLabelList2_rot(i))%sym%S,4)
!            ENDIF
!        enddo
!
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and its inverse'
!        do i=1,NoOrbs
!            WRITE(6,*) SymLabelList2_rot(i),SymLabelListInv_rot(i)
!        enddo
!        CALL neci_flush(6)
!        CALL Stop_All('SetupNatOrbLabels','Checking orbital labelling.')


    END SUBROUTINE SetupNatOrbLabels



    SUBROUTINE FillOneRDM()
        USE DetCalcData , only : Det,FCIDets,FCIDetIndex,ICILevel
        use DetBitOps, only: FindBitExcitLevel
        use bit_reps, only: decode_bit_det
! Det is the number of determinants in FCIDets.
! FCIDets contains the list of all determinants in the system in bit string representation, FCIDets(0:NIfTot,1:Det) 
! ICILevel is the max excitation level of the calculation - as in EXCITE ICILevel.
! FCIDetIndex(1:NEl) contains the index of FCIDets where each excitation level starts.
! As in FCIDetIndex(1) = 2 always I think - Excitation level 1 starts at the second determinant (after HF).
! Pretty sure FCIDetIndex always goes from 1:NEl even from truncated excite calculations.
        use hist_data, only: AllHistogram
! The elements of AllHistogram correspond to the rows of FCIDets - i.e to each determinant in the system.
! AllHistogram contains the final (normalised) amplitude of the determinant - with sign.
        IMPLICIT NONE
        INTEGER :: excit,i,j,Starti,Endi,Startj,Endj,ExcitLevel,Ex(2,1),Orbi,Orbj,nJ(NEl),Orbk,k,nI(NEl),MaxExcit
        INTEGER :: Spins
        LOGICAL :: tSign
        real(dp) :: SignDet

! Density matrix    = D_pq = < Psi | a_q+ a_p | Psi > 
!                   = sum_ij [ c_i* c_j < D_i | a_q+ a_p | D_j > ]
! Where a_p is the annihilation, and a_q+ the creation operators.
! In other words, < D_i | a_q+ a_p | D_j > will only be non-zero if D_i and D_j are connected by 
! an annihilation in p and a creation in q.

! Want to run through all determinants D_i in the final wavefunction, Psi.  For each, find all determinants, D_j
! that are connected to D_i by a single excitation - i.e. those which differ by just one orbital.
! Only need to consider those of the same excitation level or one above or one below.
! Find the differing orbitals - these will be p and q.
! Sum in the occupations of D_i and D_j (c_i x c_j) to the matrix element D_pq.
! Take, for instance, p always =< q.

! Will get the orbitals in the original labelling system - convert it to this system.

! Get a list of the wavefunction with amplitudes in order of excitation.

! Depending on the type of reduced density matrix want to:
! Run through the determinants with excitation level one less, the same and one more.

        
        FillOneRDM_Time%timer_name='FillOneRDM'
        CALL set_timer(FillOneRDM_Time,30)


!        WRITE(6,*) 'ICIIndex',ICILevel
!        WRITE(6,*) 'Det',Det
!        WRITE(6,*) 'FCIDetIndex'
!        do i=1,NEl
!            WRITE(6,*) FCIDetIndex(i)
!        enddo

!        WRITE(6,*) 'FCIDets'
!        do i=1,Det
!            WRITE(6,*) FCIDets(0:NIfTot,i),AllHistogram(i)
!        enddo
!        CALL neci_flush(6)
!        stop
        WRITE(6,*) '*** The weight of the HF determinant is : ', AllHistogram(1,1)

        WRITE(6,*) 'Beginning to fill the one-electron reduced density matrix.'

        IF(ICILevel.eq.0) THEN
            MaxExcit=NEl
        ELSE
            MaxExcit=ICILevel
        ENDIF

        do excit=0,MaxExcit         
        ! Run through all determinants D_i, in the final wavefunction, Psi. 
        ! If this is done by excitation block, we then don't have to check 
        ! the excitation level of the determinant each time.

            IF(tRotateVirtOnly.and.tSeparateOccVirt.and.(excit.eq.0)) CYCLE      
            ! The HF only involves 'occupied' orbitals - these are not required if only rotating virt.

! This next bit is a bit messy because there is no row in FCIDetIndex for the HF - there is probably a 
! tidier way to achieve the same thing, but it does the trick for now.
            IF(excit.eq.0) THEN         ! i is the HF det.
                Starti=1
                Endi=1
                Startj=1
                Endj=MIN((FCIDetIndex(2)-1),Det)   ! If i is the HF det, just run over singly excited j. 
            ELSEIF(excit.eq.MaxExcit) THEN
                Starti=FCIDetIndex(excit)
                Endi=Det
                Startj=FCIDetIndex(excit-1)
                Endj=Det
            ELSE
                Starti=FCIDetIndex(excit)
                Endi=FCIDetIndex(excit+1)-1
                IF(excit.eq.1) THEN
                    Startj=1
                    IF(NEl.lt.3) THEN
                        Endj=Det
                    ELSE
                        Endj=FCIDetIndex(3)-1
                    ENDIF
                ELSEIF(excit.eq.(MaxExcit-1)) THEN
                    Startj=FCIDetIndex(excit-1)
                    Endj=Det
                ELSE
                    Startj=FCIDetIndex(excit-1)
                    Endj=FCIDetIndex(excit+2)-1
                ENDIF
            ENDIF
!            WRITE(6,*) 'Starti',Starti
!            WRITE(6,*) 'Endi',Endi
!            WRITE(6,*) 'Startj',Startj
!            WRITE(6,*) 'Endj',Endj

            do i=Starti,Endi
            ! Then run through the determinants in that excitation level.

                do j=Startj,i
!                do j=Startj,Endj
                ! Run through all determinants D_j, with the potential to be connected to i 
                ! by a single excitation, i.e from one excitation lower to one excitation higher.
                    IF((i.gt.Det).or.(j.gt.Det)) THEN
                        CALL Stop_All('FillOneRDM',&
                            'Running through i or j larger than the number of determinants.')
                    ENDIF

                    ExcitLevel = FindBitExcitLevel(FCIDets(:,i), &
                                                   FCIDets(:,j),2)
                    ! Need to find the excitation level between D_i and D_j. 
                    ! If this is 1 - go on to add their contributions to the OneRDM.

                    IF(ExcitLevel.eq.1) THEN
                        Ex(:,:)=0
                        Ex(1,1)=ExcitLevel

                        CALL GetBitExcitation(FCIDets(:,i),FCIDets(:,j),Ex,tSign)
                        ! Gives the orbitals involved in the excitation Ex(1,1) 
                        ! in i -> Ex(2,1) in j (in spin orbitals).

                        IF(tStoreSpinOrbs) THEN
                            ! OneRDM will be in spin orbitals - simply add the orbitals involved.
                            Orbi=SymLabelListInv_rot(Ex(1,1))
                            Orbj=SymLabelListInv_rot(Ex(2,1))
                            Spins=1
                        ELSE
                            Orbi=SymLabelListInv_rot(CEILING(REAL(Ex(1,1))/2.0_dp))
                            Orbj=SymLabelListInv_rot(CEILING(REAL(Ex(2,1))/2.0_dp))
                            Spins=2
                        ENDIF
                        IF(tSign) THEN
                            SignDet=(-1.0_dp)
                        ELSE
                            SignDet=1.0_dp
                        ENDIF

                        NatOrbMat(Orbi,Orbj)=NatOrbMat(Orbi,Orbj) &
                                                + (SignDet*AllHistogram(1,i)*AllHistogram(1,j))
                        NatOrbMat(Orbj,Orbi)=NatOrbMat(Orbj,Orbi) &
                                                + (SignDet*AllHistogram(1,i)*AllHistogram(1,j))

                        ! AllHistogram are the normalised amplitudes of the determinants.
!                        IF(((AllHistogram(i)*AllHistogram(j).ne.0.0_dp).and.&
!                            (INT(G1(SymLabelList2_rot(Orbi)*2)%sym%S,4).ne.&
!                            INT(G1(SymLabelList2_rot(Orbj)*2)%sym%S,4)))&
!                        &.or.(Ex(1,1).gt.(SpatOrbs*2)).or.(Ex(2,1).gt.(SpatOrbs*2))) THEN

                        IF((AllHistogram(1,i)*AllHistogram(1,j).ne.0.0_dp).and. &
                         (INT(G1(SymLabelList2_rot(Orbi)*Spins)%sym%S,4).ne.&
                         INT(G1(SymLabelList2_rot(Orbj)*Spins)%sym%S,4))) THEN
                            WRITE(6,*) 'ERROR in symmetries'
                            WRITE(6,*) 'Ex,',Ex(1,1),Ex(2,1)
                            WRITE(6,*) CEILING(REAL(Ex(1,1)/2.0_dp)),CEILING(REAL(Ex(2,1)/2.0_dp))
                            WRITE(6,*) 'Orbi,',Orbi,'Orbj,',Orbj
                            WRITE(6,*) 'Sym(Orbi)',INT(G1(SymLabelList2_rot(Orbi)*Spins)%sym%S,4),'Sym(Orbj)', &
                                INT(G1(SymLabelList2_rot(Orbj)*Spins)%sym%S,4)
                            call decode_bit_det (nI, FCIDets(0:NIfTot,i))
                            WRITE(6,*) 'i',nI
                            call decode_bit_det (nJ, FCIDets(0:NIfTot,j))
                            WRITE(6,*) 'j',nJ
                            WRITE(6,*) 'AllHistogram(1,i)',AllHistogram(1,i)
                            WRITE(6,*) 'AllHistogram(1,j)',AllHistogram(1,j)
                            CALL neci_flush(6)
                            CALL Stop_All('FillOneRDM','Non-zero element between different symmetries.')
                        ENDIF

                    ELSEIF(ExcitLevel.eq.0) THEN
                        CALL Decode_Bit_Det(nJ,FCIDets(0:NIfTot,j))
                        do k=1,NEl
!                            WRITE(6,*) 'k',k
                            IF(tStoreSpinOrbs) THEN
                                Orbk=SymLabelListInv_rot(nJ(k))
                            ELSE
                                Orbk=SymLabelListInv_rot(CEILING(REAL(nJ(k))/2.0_dp))
                            ENDIF
                            NatOrbMat(Orbk,Orbk)=NatOrbMat(Orbk,Orbk)+(AllHistogram(1,j)**2)
!                            NatOrbMat(Orbk,Orbk)=NatOrbMat(Orbk,Orbk)+(0.5 * (AllHistogram(j)**2))
                            ! 0.5 x because this will be added twice since we are not currently 
                            ! restricting i<k or anything.
                        enddo
                    ENDIF
                        
                enddo

            enddo
        enddo

        WRITE(6,*) 'Filled OneRDM'
!        do i=1,NoOrbs
!            do j=1,NoOrbs
!                WRITE(6,'(F10.5)',advance='no') NatOrbMat(j,i)
!            enddo
!            WRITE(6,*) ''
!        enddo
!        CALL neci_flush(6)
!        stop

        CALL halt_timer(FillOneRDM_Time)



    END SUBROUTINE FillOneRDM



    SUBROUTINE FillMP2VDM()
! In this routine, the natural orbital matrix is calculated from the MP2 variational density matrix.

! MP2VDM = D2_ab = sum_ijc [ t_ij^ac ( 2 t_ij^bc - t_ji^bc ) ]
! Where :  t_ij^ac = - < ab | ij > / ( E_a - E_i + E_b - Ej )
! Ref : J. Chem. Phys. 131, 034113 (2009) - note: in Eqn 1, the cb indices are the wrong way round (should be bc).
        USE Integrals_neci , only : GetUMatEl
        USE SystemData , only : tUEG
        use constants, only: dp
        INTEGER :: a,b,c,i,j,a2,b2,c2,i2,j2,x,y,z,w
        INTEGER :: Startab,Endab,NoOcc,NoOccC,Startc,Endc,Starti,Endi,Startj,Endj
        real(dp) :: MP2VDMSum
        CHARACTER(len=*), PARAMETER :: this_routine='FillMP2VDM'
        HElement_t :: HEl01,HEl02

#ifdef __CMPLX
         call stop_all('FillMP2VDM', 'Natural Orbitals not implemented for complex orbitals.')
#endif
! Calculating the MP2VDM (D2_ab) matrix whose eigenvectors become the transformation matrix.        
! This goes in the natural orbital matrix of this module.
! The eigenvalues are the occupation numbers of the new orbitals.  These should decrease exponentially so that 
! when we remove the orbitals with small occupation numbers we should have little affect on the energy.


! For the MP2VDM, we always only rotate the virtual orbitals - denomonator term of the above expression would 
! be 0 if a and b were occupied.
! The orbital labels are ordered occupied then virtual if spatial orbitals are being used,
! otherwise they go occupied beta, virtual beta, occupied alpha, virtual alpha.
! This is so the alpha and beta spins can be diagonalised separately and we can keep track of which is which 
! when the evectors are reordered and maintain spin symmetry.


        WRITE(6,*) 'Filling MP2VDM nat orb matrix'
        CALL neci_flush(6)
        
        FillMP2VDM_Time%timer_name='FillMP2VDM'
        CALL set_timer(FillMP2VDM_Time,30)

!        WRITE(6,*) 'nOccBeta',nOccBeta
!        WRITE(6,*) 'nOccAlpha',nOccAlpha
!        WRITE(6,*) 'tStoreSpinOrbs',tStoreSpinOrbs

        do x=1,NoSpinCyc
            IF(x.eq.1) THEN
                IF(tStoreSpinOrbs) THEN
                    NoOcc=nOccBeta
                ELSE
                    NoOcc=NEl/2
                ENDIF
                Startab=NoOcc+1
                Endab=SpatOrbs
            ELSEIF(x.eq.2) THEN
                NoOcc=nOccAlpha
                Startab=SpatOrbs+NoOcc+1
                Endab=NoOrbs
            ENDIF

            ! a and b must be of the same spin to mix, so only need to run over both beta then both alpha.
            do a2=Startab,Endab
                a=SymLabelList2_rot(a2)
!                do b2=Startab,Endab
                do b2=Startab,a2

                    b=SymLabelList2_rot(b2)

                    MP2VDMSum=0.0_dp
!                    WRITE(6,*) 'a',a,'b',b,'a2',a2,'b2',b2

                    ! when a and b beta, run over both alpha and beta virtual for c, then both alpha 
                    ! and beta virtual for both i and j etc. 

                    do y=1,NoSpinCyc
                        IF(y.eq.1) THEN
                            IF(tStoreSpinOrbs) THEN
                                NoOccC=nOccBeta
!                                NoOccC=MAX(nOccBeta,nOccAlpha)
                            ELSE
                                NoOccC=NEl/2
                            ENDIF
                            Startc=NoOccC+1
                            Endc=SpatOrbs
!                            Startc=1
!                            Endc=SpatOrbs
                        ELSEIF(y.eq.2) THEN
                            NoOccC=nOccAlpha
!                            NoOccC=MAX(nOccBeta,nOccAlpha)
                            Startc=SpatOrbs+NoOccC+1
                            Endc=NoOrbs
!                            Startc=SpatOrbs+1
!                            Endc=NoOrbs
                        ENDIF


                        do c2=Startc,Endc
                            c=SymLabelList2_rot(c2)

                            do z=1,NoSpinCyc
                                IF(z.eq.1) THEN
                                    Starti=1
                                    IF(tStoreSpinOrbs) THEN
                                        Endi=nOccBeta
!                                        Endi=MAX(nOccBeta,nOccAlpha)
                                    ELSE
                                        Endi=NEl/2
                                    ENDIF
                                ELSEIF(z.eq.2) THEN
                                    Starti=1+SpatOrbs
                                    Endi=SpatOrbs+nOccAlpha
!                                    Endi=SpatOrbs+MAX(nOccBeta,nOccAlpha)
                                ENDIF

                                do i2=Starti,Endi
                                    i=SymLabelList2_rot(i2)

                                    do w=1,NoSpinCyc
                                        IF(w.eq.1) THEN
                                            Startj=1
                                            IF(tStoreSpinOrbs) THEN
                                                Endj=nOccBeta
!                                                Endj=MAX(nOccBeta,nOccAlpha)
                                            ELSE
                                                Endj=NEl/2
                                            ENDIF
                                        ELSEIF(w.eq.2) THEN
                                            Startj=1+SpatOrbs
                                            Endj=SpatOrbs+nOccAlpha
!                                            Endj=SpatOrbs+MAX(nOccBeta,nOccAlpha)
                                        ENDIF


                                        do j2=Startj,Endj
                                            j=SymLabelList2_rot(j2)

                                            IF(tUEG) THEN
                                                HEl01=GETUMATEL(a,c,i,j)
                                                HEl02=GETUMATEL(b,c,i,j)
                                                MP2VDMSum=MP2VDMSum+&
                                                    &(( (REAL(HEl01,dp)) * (2.0_dp*(REAL(HEl02,dp))) )/&
                                                    &( (ARR(2*i,2)+ARR(2*j,2)-ARR(2*a,2)-ARR(2*c,2)) &
                                                    &* (ARR(2*i,2)+ARR(2*j,2)-ARR(2*b,2)-ARR(2*c,2)) ) )

                                                HEl02=GETUMATEL(c,b,i,j)
                                                MP2VDMSum=MP2VDMSum-&
                                                    &(( (REAL(HEl01,dp)) * (REAL(HEl02,dp)) )/&
                                                    &( (ARR(2*i,2)+ARR(2*j,2)-ARR(2*a,2)-ARR(2*c,2)) * &
                                                    &(ARR(2*i,2)+ARR(2*j,2)-ARR(2*c,2)-ARR(2*b,2)) ) )

                                            ELSEIF(tStoreSpinOrbs) THEN
                                                IF((ARR(i,2)+ARR(j,2)-ARR(a,2)-ARR(c,2)).eq.0.0_dp) THEN
                                                    IF((REAL(UMAT(UMatInd(a,c,i,j,0,0)),dp)).ne.0.0_dp) THEN
                                                        WRITE(6,*) i,j,a,c,REAL(UMAT(UMatInd(a,c,i,j,0,0)),dp)
                                                        CALL Stop_All(this_routine,"Dividing a non-zero by zero.")
                                                    ENDIF
                                                ENDIF
                                                MP2VDMSum=MP2VDMSum+&
                                                   (((REAL(UMAT(UMatInd(a,c,i,j,0,0)),dp)) & 
                                                   * (2.0_dp*(REAL(UMAT(UMatInd(b,c,i,j,0,0)),dp))))/&
                                                   ( (ARR(i,2)+ARR(j,2)-ARR(a,2)-ARR(c,2)) &
                                                   * (ARR(i,2)+ARR(j,2)-ARR(b,2)-ARR(c,2)) ) )
                                                MP2VDMSum=MP2VDMSum-&
                                                    (( (REAL(UMAT(UMatInd(a,c,i,j,0,0)),dp)) &
                                                    * (REAL(UMAT(UMatInd(c,b,i,j,0,0)),dp)) )/ &
                                                    ( (ARR(i,2)+ARR(j,2)-ARR(a,2)-ARR(c,2)) &
                                                    * (ARR(i,2)+ARR(j,2)-ARR(c,2)-ARR(b,2)) ) )
                                            ELSE
                                                MP2VDMSum=MP2VDMSum+&
                                                    (( (REAL(UMAT(UMatInd(a,c,i,j,0,0)),dp)) &
                                                    * (2.0_dp*(REAL(UMAT(UMatInd(b,c,i,j,0,0)),dp))) )/&
                                                    ((ARR(2*i,2)+ARR(2*j,2)-ARR(2*a,2)-ARR(2*c,2)) &
                                                    * (ARR(2*i,2)+ARR(2*j,2)-ARR(2*b,2)-ARR(2*c,2))))
                                               MP2VDMSum=MP2VDMSum-&
                                                    (( (REAL(UMAT(UMatInd(a,c,i,j,0,0)),dp)) &
                                                    * (REAL(UMAT(UMatInd(c,b,i,j,0,0)),dp)) )/&
                                                    ( (ARR(2*i,2)+ARR(2*j,2)-ARR(2*a,2)-ARR(2*c,2)) &
                                                    * (ARR(2*i,2)+ARR(2*j,2)-ARR(2*c,2)-ARR(2*b,2))))
                                            ENDIF

                                        enddo
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
!                    WRITE(6,*) 'MP2VDMSum',MP2VDMSum
                    NatOrbMat(a2,b2)=MP2VDMSum
                    NatOrbMat(b2,a2)=MP2VDMSum
                enddo
            enddo
        enddo

        WRITE(6,*) 'Finished filling MP2VDM'

!        do i=1,NoOrbs
!            do j=1,NoOrbs
!                WRITE(6,*) NatOrbMat(i,j)
!            enddo
!        enddo
!        CALL neci_flush(6)
!        CALL Stop_All('','')

        CALL halt_timer(FillMP2VDM_Time)


    END SUBROUTINE FillMP2VDM




    SUBROUTINE DiagNatOrbMat()
! The diagonalisation routine reorders the orbitals in such a way that the corresponding orbital 
! labels are lost.
! In order to keep the spin and spatial symmetries, each symmetry must be fed into the diagonalisation 
! routine separately.
! The best way to do this is to order the orbitals so that all the alpha orbitals follow all the 
! beta orbitals, with the occupied orbitals first, in terms of symmetry, and the virtual second, 
! also ordered by symmetry.
! This gives us flexibility w.r.t rotating only the occupied or only virtual and looking at high spin states.
        use MemoryManager, only: TagIntType
        IMPLICIT NONE
        real(dp) :: SumTrace,SumDiagTrace
        real(dp) , ALLOCATABLE :: WORK2(:),EvaluesSym(:),NOMSym(:,:)
        INTEGER :: ierr,i,j,x,z,Sym,LWORK2,SymStartInd,NoSymBlock,PrevSym,StartOccVirt,EndOccVirt,Prev,NoOcc
        INTEGER(TagIntType) :: EvaluesSymTag,NOMSymTag,WORK2Tag
        CHARACTER(len=*), PARAMETER :: this_routine='DiagNatOrbMat'

 
        DiagNatOrbMat_Time%timer_name='DiagNatOrb'
        CALL set_timer(DiagNatOrbMat_Time,30)

        do x=1,NoSpinCyc
            IF(tSeparateOccVirt) THEN
                IF(x.eq.1) THEN
                    IF(tStoreSpinOrbs) THEN
                        NoOcc=nOccBeta
                    ELSE
                        NoOcc=NEl/2
                    ENDIF
                    Prev=0
                ELSEIF(x.eq.2) THEN
                    NoOcc=nOccAlpha
                    Prev=SpatOrbs
                ENDIF
            ELSE
                NoOcc=0
            ENDIF
            IF(tRotateVirtOnly) THEN
                do i=1,NoOcc
                    do j=1,SpatOrbs
                        NatOrbMat(i+Prev,j+Prev)=0.0_dp
                        NatOrbMat(j+Prev,i+Prev)=0.0_dp
                        IF(i.eq.j) NatOrbMat(i+Prev,j+Prev)=1.0_dp
                    enddo
                    Evalues(i+Prev)=1.0_dp
                enddo
            ELSEIF(tRotateOccOnly) THEN
                do i=NoOcc+1,SpatOrbs
                    do j=1,SpatOrbs
                        NatOrbMat(i+Prev,j+Prev)=0.0_dp
                        NatOrbMat(j+Prev,i+Prev)=0.0_dp
                        IF(i.eq.j) NatOrbMat(i+Prev,j+Prev)=1.0_dp
                    enddo
                    Evalues(i+Prev)=1.0_dp
                enddo
            ELSEIF(tSeparateOccVirt) THEN
                do i=1,NoOcc
                    do j=NoOcc+1,SpatOrbs
                        NatOrbMat(i+Prev,j+Prev)=0.0_dp
                        NatOrbMat(j+Prev,i+Prev)=0.0_dp
                    enddo
                enddo
            ENDIF
        enddo

! Test that we're not breaking symmetry.
        do i=1,NoOrbs
            do j=1,NoOrbs
                IF(tStoreSpinOrbs) THEN
!                    WRITE(6,*) INT(G1(SymLabelList2_rot(i))%sym%S,4),INT(G1(SymLabelList2_rot(j))%sym%S,4),NatOrbMat(i,j)
                    IF((INT(G1(SymLabelList2_rot(i))%sym%S,4).ne.INT(G1(SymLabelList2_rot(j))%sym%S,4))) THEN
                        IF(ABS(NatOrbMat(i,j)).ge.1.0E-15) THEN
                            WRITE(6,'(6A8,A20)') 'i','j','Label i','Label j','Sym i','Sym j','Matrix value'
                            WRITE(6,'(6I3,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j), &
                                INT(G1(SymLabelList2_rot(i))%sym%S,4), &
                                INT(G1(SymLabelList2_rot(j))%sym%S,4),NatOrbMat(i,j)
                            IF(tUseMP2VarDenMat) THEN
                                WRITE(6,*) '**WARNING** - There is a non-zero NatOrbMat value between " &
                                 & //"orbitals of different symmetry.'
                                WRITE(6,*) 'These elements will be ignored, and the symmetry maintained " &
                                 & //"in the final transformation matrix.'
                            ELSE
                                CALL Stop_All(this_routine,'Non-zero NatOrbMat value between different symmetries.')
                            ENDIF
                        ENDIF
                        NatOrbMat(i,j)=0.0_dp
                    ENDIF
                ELSE
!                    WRITE(6,*) INT(G1(SymLabelList2_rot(i)*2)%sym%S,4),INT(G1(SymLabelList2_rot(j)*2)%sym%S,4),NatOrbMat(i,j)
                    IF((INT(G1(SymLabelList2_rot(i)*2)%sym%S,4).ne.INT(G1(SymLabelList2_rot(j)*2)%sym%S,4))) THEN
                        IF(ABS(NatOrbMat(i,j)).ge.1.0E-15) THEN
                            WRITE(6,'(6A8,A20)') 'i','j','Label i','Label j','Sym i','Sym j','Matrix value'
                            WRITE(6,'(6I3,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j), &
                             INT(G1(SymLabelList2_rot(i)*2)%sym%S,4),INT(G1(SymLabelList2_rot(j)*2)%sym%S,4),NatOrbMat(i,j)
                            IF(tUseMP2VarDenMat) THEN
                                WRITE(6,*) '**WARNING** - There is a non-zero NatOrbMat value between orbitals of "&
                                & //"different symmetry.'
                                WRITE(6,*) 'These elements will be ignored, and the symmetry maintained in the " &
                                & //"final transformation matrix.'
                            ELSE
                                CALL Stop_All(this_routine,'Non-zero NatOrbMat value between different symmetries.')
                            ENDIF
                        ENDIF
                        NatOrbMat(i,j)=0.0_dp
                    ENDIF
                ENDIF
            enddo
        enddo

        SumTrace=0.0_dp
        do i=1,NoOrbs
            SumTrace=SumTrace+NatOrbMat(i,i)
        enddo

        WRITE(6,*) 'Calculating eigenvectors and eigenvalues of NatOrbMat'
        CALL neci_flush(6)

        ! If we are using spin orbitals, need to feed in the alpha and beta spins separately.
        ! Otherwise these jumble up and the final ordering is uncorrect. 
        ! There should be no non-zero values between these, but can put a check in for this.

        do x=1,NoSpinCyc

! If we want to maintain the symmetry, we cannot have all the orbitals jumbled up when the diagonaliser 
! reorders the eigenvectors.
! Must instead feed each symmetry block in separately.
! This means that although the transformed orbitals are jumbled within the symmetry blocks, the symmetry 
! labels are all that are relevant and these are unaffected.
            StartOccVirt=1
            EndOccVirt=2
            IF(tRotateVirtOnly) StartOccVirt=2
            IF(tRotateOccOnly) EndOccVirt=1

            do z=StartOccVirt,EndOccVirt
                IF(x.eq.1) THEN
                    IF(z.eq.1) PrevSym=1
                    IF(z.eq.2) PrevSym=9
                ELSEIF(x.eq.2) THEN
                    IF(z.eq.1) PrevSym=17
                    IF(z.eq.2) PrevSym=25
                ENDIF

                Sym=0
                LWORK2=-1
                do while (Sym.le.7)

                    NoSymBlock=SymLabelCounts2_rot(2,Sym+PrevSym)

                    SymStartInd=SymLabelCounts2_rot(1,Sym+PrevSym)-1
                    ! This is one less than the index that the symmetry starts, so that when we run through i=1,..., 
                    ! we can start at SymStartInd+i.

                    IF(NoSymBlock.gt.1) THEN

                        ALLOCATE(NOMSym(NoSymBlock,NoSymBlock),stat=ierr)
                        CALL LogMemAlloc('NOMSym',NoSymBlock**2,8,this_routine,NOMSymTag,ierr)
                        IF(ierr.ne.0) CALL Stop_All(this_routine,"Problem allocating NOMSym.")
                        ALLOCATE(EvaluesSym(NoSymBlock),stat=ierr)
                        CALL LogMemAlloc('EvaluesSym',NoSymBlock,8,this_routine,EvaluesSymTag,ierr)
                        IF(ierr.ne.0) CALL Stop_All(this_routine,"Problem allocating EvaluesSym.")

                        LWORK2=3*NoSymBlock+1
                        ALLOCATE(WORK2(LWORK2),stat=ierr)
                        CALL LogMemAlloc('WORK2',LWORK2,8,this_routine,WORK2Tag,ierr)
                        IF(ierr.ne.0) CALL Stop_All(this_routine,"Problem allocating WORK2.")

                        do j=1,NoSymBlock
                            do i=1,NoSymBlock
                                NOMSym(i,j)=NatOrbMat(SymStartInd+i,SymStartInd+j)
                            enddo
                        enddo


                        WRITE(6,*) '*****'
                        WRITE(6,*) 'Symmetry ',Sym, 'with x ',x,' and z ',z,' has ',NoSymBlock,' orbitals.'
                        WRITE(6,*) 'The NatOrbMat for this symmetry block is '
                        do i=1,NoSymBlock
                            do j=1,NoSymBlock
                                WRITE(6,'(F20.10)',advance='no') NOMSym(j,i)
                            enddo
                            WRITE(6,*) ''
                        enddo

                        CALL DSYEV('V','L',NoSymBlock,NOMSym,NoSymBlock,EvaluesSym,WORK2,LWORK2,ierr)
                        ! NOMSym goes in as the original NOMSym, comes out as the eigenvectors (Coefficients).
                        ! EvaluesSym comes out as the eigenvalues in ascending order.

                        WRITE(6,*) 'After diagonalization, the e-vectors (diagonal elements) of this matrix are ,'
                        do i=1,NoSymBlock
                            WRITE(6,'(F20.10)',advance='no') EvaluesSym(i)
                        enddo
                        WRITE(6,*) ''
                        WRITE(6,*) 'These go from orbital ,',SymStartInd+1,' to ',SymStartInd+NoSymBlock

                        do i=1,NoSymBlock
                            Evalues(SymStartInd+i)=EvaluesSym(i)
                        enddo

                        ! CAREFUL if eigenvalues are put in ascending order, this may not be correct, 
                        ! with the labelling system. may be better to just take coefficients and transform 
                        ! TMAT2DRot in transform2elints. a check that comes out as diagonal is a check of 
                        ! this routine anyway.

                        WRITE(6,*) 'The eigenvectors (coefficients) for symmetry block ',Sym
                        do i=1,NoSymBlock
                            do j=1,NoSymBlock
                                WRITE(6,'(F20.10)',advance='no') NOMSym(j,i)
                            enddo
                            WRITE(6,*) ''
                        enddo
                 
                        do j=1,NoSymBlock
                            do i=1,NoSymBlock
                                NatOrbMat(SymStartInd+i,SymStartInd+j)=NOMSym(i,j)
                            enddo
                        enddo
                        ! Directly fill the coefficient matrix with the eigenvectors from the diagonalization.

                        DEALLOCATE(WORK2)
                        CALL LogMemDealloc(this_routine,WORK2Tag)

                        DEALLOCATE(NOMSym)
                        CALL LogMemDealloc(this_routine,NOMSymTag)

                        DEALLOCATE(EvaluesSym)
                        CALL LogMemDealloc(this_routine,EvaluesSymTag)

                    ELSEIF(NoSymBlock.eq.1) THEN
                        ! The eigenvalue is the lone value, while the eigenvector is 1.

                        Evalues(SymStartInd+1)=NatOrbMat(SymStartInd+1,SymStartInd+1)
                        NatOrbMat(SymStartInd+1,SymStartInd+1)=1.0_dp
                        WRITE(6,*) '*****'
                        WRITE(6,*) 'Symmetry ',Sym,' has only one orbital.'
                        WRITE(6,*) 'Copying diagonal element ,',SymStartInd+1,'to NatOrbMat'
                    ENDIF

                    Sym=Sym+1
                enddo
            enddo
        enddo

        WRITE(6,*) 'Matrix diagonalised'
        CALL neci_flush(6)

        SumDiagTrace=0.0_dp
        do i=1,NoOrbs
            SumDiagTrace=SumDiagTrace+Evalues(i)
        enddo
        IF((ABS(SumDiagTrace-SumTrace)).gt.10.0_dp) THEN
            WRITE(6,*) 'Sum of diagonal NatOrbMat elements : ',SumTrace
            WRITE(6,*) 'Sum of eigenvalues : ',SumDiagTrace
!            CALL Stop_All(this_routine,'The trace of the 1RDM matrix before diagonalisation is 
!                                            not equal to that after.')
            WRITE(6,*) 'WARNING, The trace of the 1RDM matrix before diagonalisation is not equal to that after.'
        ENDIF

        CALL halt_timer(DiagNatOrbMat_Time)

    END SUBROUTINE DiagNatOrbMat



    SUBROUTINE OrderCoeffT1()
        USE RotateOrbsData , only : SymLabelList3_rot
        USE Logging , only : tTruncRODump
        IMPLICIT NONE
        INTEGER :: x,i,ierr,StartSort,EndSort,NoOcc
        CHARACTER(len=*), PARAMETER :: this_routine='OrderCoeffT1'
        

! Here, if symmetry is kept, we are going to have to reorder the eigenvectors according to the size of the 
! eigenvalues, while taking the orbital labels (and therefore symmetries) with them. This will be put back 
! into MP2VDM from MP2VDMTemp.

! Want to reorder the eigenvalues from largest to smallest, taking the eigenvectors with them and the 
! symmetry as well.  
! If using spin orbitals, do this for the alpha spin and then the beta.
 
        OrderCoeff_Time%timer_name='OrderCoeff'
        CALL set_timer(OrderCoeff_Time,30)


        IF(tTruncRODump) THEN
            ! If we are truncating, the orbitals stay in this order, so we want to take their symmetries with them.
            ALLOCATE(SymOrbs_rotTemp(NoOrbs),stat=ierr)
            CALL LogMemAlloc('SymOrbs_rotTemp',NoOrbs,4,this_routine,SymOrbs_rotTempTag,ierr)
            SymOrbs_rotTemp(:)=0

            IF(tStoreSpinOrbs) THEN
                do i=1,NoOrbs
                    SymOrbs_rotTemp(i)=INT(G1(SymLabelList2_rot(i))%sym%S,4)
                enddo
            ELSE 
                do i=1,NoOrbs
                    SymOrbs_rotTemp(i)=INT(G1(SymLabelList2_rot(i)*2)%sym%S,4)
                enddo
            ENDIF

            do x=1,NoSpinCyc

                IF(x.eq.1) THEN
                    IF(tSeparateOccVirt) THEN
                        IF(tStoreSpinOrbs) THEN
                            NoOcc=nOccBeta
                        ELSE
                            NoOcc=NEl/2
                        ENDIF
                    ELSE
                        NoOcc=0
                    ENDIF
                    StartSort=1
                    EndSort=SpatOrbs
                    IF(tRotateVirtOnly) StartSort=NoOcc+1
                    IF(tRotateOccOnly) EndSort=NoOcc
                ELSEIF(x.eq.2) THEN
                    IF(tSeparateOccVirt) THEN
                        NoOcc=nOccAlpha
                    ELSE
                        NoOcc=0
                    ENDIF
                    StartSort=SpatOrbs+1
                    EndSort=NoOrbs
                    IF(tRotateVirtOnly) StartSort=SpatOrbs+NoOcc+1
                    IF(tRotateOccOnly) EndSort=NoOcc+SpatOrbs
                ENDIF

                call sort (Evalues(startSort:endSort), &
                           natOrbMat(startSort:endSort, startSort:endSort), &
                           SymOrbs_rotTemp(startSort:endSort))

            enddo
               
        ELSE
            ! If we are not truncating, the orbitals get put back into their original order, so the symmetry 
            ! information is still correct, no need for the SymOrbs_rot array.
            ! Instead, just take the labels of SymLabelList3_rot with them.

            do x=1,NoSpinCyc

                IF(x.eq.1) THEN
                    IF(tSeparateOccVirt) THEN
                        IF(tStoreSpinOrbs) THEN
                            NoOcc=nOccBeta
                        ELSE
                            NoOcc=NEl/2
                        ENDIF
                    ELSE
                        NoOcc=0
                    ENDIF
                    StartSort=1
                    EndSort=SpatOrbs
                    IF(tRotateOccOnly) EndSort=NoOcc
                    IF(tRotateVirtOnly) StartSort=NoOcc+1

                ELSEIF(x.eq.2) THEN
                    IF(tSeparateOccVirt) THEN
                        NoOcc=nOccAlpha
                    ELSE
                        NoOcc=0
                    ENDIF
                    StartSort=SpatOrbs+1
                    EndSort=NoOrbs
                    IF(tRotateOccOnly) EndSort=SpatOrbs+NoOcc
                    IF(tRotateVirtOnly) StartSort=SpatOrbs+NoOcc+1
                ENDIF

                call sort (EValues(startSort:endSort), &
                           NatOrbMat(startSort:endSort, startSort:endSort), &
                           SymLabelList3_rot(startSort:endSort))
            enddo 
            
        ENDIF

        CALL halt_timer(OrderCoeff_Time)

        WRITE(6,*) 'Eigen-values: '
        do i=1,NoOrbs
            WRITE(6,*) Evalues(i)
        enddo
       

    END SUBROUTINE OrderCoeffT1



    SUBROUTINE FillCoeffT1
        USE RotateOrbsData , only : CoeffT1,SymLabelList3_rot,SymOrbs_rot,SymOrbs_rotTag,&
                                    TruncEval,NoRotOrbs,EvaluesTrunc,EvaluesTruncTag
        USE Logging , only : tTruncRODump,tTruncDumpbyVal
        IMPLICIT NONE
        INTEGER :: l,k,i,j,NoRotAlphBet, io1, io2
        CHARACTER(len=*), PARAMETER :: this_routine='FillCoeffT1'
        CHARACTER(len=5) :: Label
        CHARACTER(len=20) :: LabelFull
        real(dp) :: OccEnergies(1:NoRotOrbs)
  
        FillCoeff_Time%timer_name='FillCoeff'
        CALL set_timer(FillCoeff_Time,30)

        WRITE(6,*) 'NatOrbMat'
        do i = 1, nBasis
            do j = 1, nBasis
                WRITE(6,'(F10.6)',advance='no') NatOrbMat(j,i)
            enddo
            WRITE(6,*) ''
        enddo

        IF(tTruncRODump) THEN

            IF(tTruncDumpbyVal) THEN
                NoFrozenVirt=0
                IF(tStoreSpinOrbs) THEN
                    do i=SpatOrbs,1,-1
                        IF(Evalues(i).gt.TruncEval) EXIT
                        IF(Evalues(i+SpatOrbs).gt.TruncEval) EXIT
                        NoFrozenVirt=NoFrozenVirt+2
                    enddo
                    IF(NoFrozenVirt.ge.(NoOrbs-NEl)) CALL Stop_All(this_routine,&
                                                                'Freezing all virtual orbitals.')
                ELSE
                    do i=SpatOrbs,1,-1
                        IF(Evalues(i).gt.TruncEval) EXIT
                        NoFrozenVirt=NoFrozenVirt+1
                    enddo
                    IF(NoFrozenVirt.ge.(SpatOrbs-(NEl/2))) CALL Stop_All(this_routine,&
                                                                'Freezing all virtual orbitals.')
                ENDIF
                NoRotOrbs=NoOrbs-NoFrozenVirt
            ENDIF

            ALLOCATE(SymOrbs_rot(NoOrbs),stat=ierr)
            CALL LogMemAlloc('SymOrbs_rot',NoOrbs,4,this_routine,SymOrbs_rotTag,ierr)
            SymOrbs_rot(:)=0

            ALLOCATE(EvaluesTrunc(NoOrbs-NoFrozenVirt),stat=ierr)
            CALL LogMemAlloc('EvaluesTrunc',NoOrbs-NoFrozenVirt,4,this_routine,EvaluesTruncTag,ierr)
            EvaluesTrunc(:)=0.0_dp

            IF(tStoreSpinOrbs) THEN
                NoRotAlphBet=SpatOrbs-(NoFrozenVirt/2)
            ELSE 
                NoRotAlphBet=NoOrbs-NoFrozenVirt
            ENDIF


            IF(tStoreSpinOrbs) THEN                                            
! As we reorder these so that they are truncated, we also need to pair up symmetries.

! Get the beta symmetry - find the alpha to match etc.
! 
!                IF(nOccBeta.ge.nOccAlpha) THEN
!                    k=1
!                    do i=1,NoRotAlphBet
!                        tSymFound=.false.
!                        SymFirst=SymOrbs_rotTemp(i)
!                        CoeffT1(:,k)=NatOrbMat(:,i)
!                        EvaluesTrunc(k)=Evalues(i)
!                        SymOrbs_rot(k)=SymOrbs_rotTemp(i)
!                        do j=SpatOrbs+1,SpatOrbs+NoRotAlphBet
!                            IF(SymOrbs_rotTemp(j).eq.SymFirst) THEN
!                                SymOrbs_rot(k+1)=SymOrbs_rotTemp(j)
!                                CoeffT1(:,k+1)=NatOrbMat(:,j)
!                                EvaluesTrunc(k+1)=Evalues(j)
!                                SymOrbs_rotTemp(j)=9
!                                tSymFound=.true.
!                                EXIT
!                            ENDIF
!                        enddo
!                        IF(.not.tSymFound) THEN
!                            do j=SpatOrbs+1,SpatOrbs+NoRotAlphBet
!                                IF(SymOrbs_rotTemp(j).lt.9) THEN
!                                    SymOrbs_rot(k+1)=SymOrbs_rotTemp(j)
!                                    CoeffT1(:,k+1)=NatOrbMat(:,j)
!                                    EvaluesTrunc(k+1)=Evalues(j)
!                                    SymOrbs_rotTemp(j)=9
!                                    EXIT
!                                ENDIF
!                            enddo
!                        ENDIF
!                        k=k+2
!                        IF(k.gt.(NoRotAlphBet*2)) EXIT
!                    enddo
!                ELSE
!                    k=2
!                    do i=SpatOrbs+1,SpatOrbs+NoRotAlphBet
!                        tSymFound=.false.
!                        SymFirst=SymOrbs_rotTemp(i)
!                        CoeffT1(:,k)=NatOrbMat(:,i)
!                        EvaluesTrunc(k)=Evalues(i)
!                        SymOrbs_rot(k)=SymOrbs_rotTemp(i)
!                        do j=1,NoRotAlphBet
!                            IF(SymOrbs_rotTemp(j).eq.SymFirst) THEN
!                                SymOrbs_rot(k-1)=SymOrbs_rotTemp(j)
!                                CoeffT1(:,k-1)=NatOrbMat(:,j)
!                                EvaluesTrunc(k-1)=Evalues(j)
!                                SymOrbs_rotTemp(j)=9
!                                tSymFound=.true.
!                                EXIT
!                            ENDIF
!                        enddo

!                        IF(.not.tSymFound) THEN
!                            do j=1,NoRotAlphBet
!                                IF(SymOrbs_rotTemp(j).lt.9) THEN
!                                    SymOrbs_rot(k-1)=SymOrbs_rotTemp(j)
!                                    CoeffT1(:,k-1)=NatOrbMat(:,j)
!                                    EvaluesTrunc(k-1)=Evalues(j)
!                                    SymOrbs_rotTemp(j)=9
!                                    EXIT
!                                ENDIF
!                            enddo
!                        ENDIF
! 
!                        k=k+2
!                        IF(k.gt.(NoRotAlphBet*2)) EXIT
!                    enddo
!                ENDIF

                CALL CalcOccEnergies(OccEnergies)

!First nOccBeta, then nOccAlpha.
                do i=1,(2*nOccBeta),2
                    k=1
                    do while(OccEnergies(k).eq.0.0_dp)
                        k=k+2
                    enddo
                    do j=1,(2*nOccBeta),2
                        IF((OccEnergies(j).lt.OccEnergies(k)).and.(OccEnergies(j).ne.0.0_dp)) k=j
                    enddo
                    l=CEILING(REAL(k)/2.0_dp)
                    CoeffT1(:,i)=NatOrbMat(:,l)
                    EvaluesTrunc(i)=Evalues(l)
                    SymOrbs_rot(i)=SymOrbs_rotTemp(l)
                    OccEnergies(k)=0.0_dp
                enddo
                do i=2,(2*nOccAlpha),2
                    k=2
                    do while(OccEnergies(k).eq.0.0_dp)
                        k=k+2
                    enddo
                    do j=2,(2*nOccAlpha),2
                        IF((OccEnergies(j).lt.OccEnergies(k)).and.(OccEnergies(j).ne.0.0_dp)) k=j
                    enddo
                    l=(k/2)+SpatOrbs
                    CoeffT1(:,i)=NatOrbMat(:,l)
                    EvaluesTrunc(i)=Evalues(l)
                    SymOrbs_rot(i)=SymOrbs_rotTemp(l)
                    OccEnergies(k)=0.0_dp
                enddo
                
!Need to fill coeffT1 so that it goes alpha beta alpha beta.
                k=(2*nOccBeta)+1
                do i=nOccBeta+1,NoRotAlphBet
                    CoeffT1(:,k)=NatOrbMat(:,i)
                    EvaluesTrunc(k)=Evalues(i)
                    SymOrbs_rot(k)=SymOrbs_rotTemp(i)
                    k=k+2
                enddo
                k=(2*nOccAlpha)+2
                do i=SpatOrbs+nOccAlpha+1,SpatOrbs+NoRotAlphBet
                    CoeffT1(:,k)=NatOrbMat(:,i)
                    SymOrbs_rot(k)=SymOrbs_rotTemp(i)
                    EvaluesTrunc(k)=Evalues(i)
                    k=k+2
                enddo
 
            ELSE

!Order occupied in terms of energy again - this makes sure freezing etc doesn't get screwed up.                    
                CALL CalcOccEnergies(OccEnergies)
!OccEnergies has the orbital energies as they are ordered currently - need to put NatOrbMat into 
!CoeffT1 so that this goes from lowest energy to highest. 

                do i=1,NEl/2
                    k=1
                    do while(OccEnergies(k).eq.0.0_dp)
                        k=k+1
                    enddo
                    do j=1,NEl/2
                        IF((OccEnergies(j).lt.OccEnergies(k)).and.(OccEnergies(j).ne.0.0_dp)) k=j
                    enddo
                    CoeffT1(:,i)=NatOrbMat(:,k)
                    EvaluesTrunc(i)=Evalues(k)
                    SymOrbs_rot(i)=SymOrbs_rotTemp(k)
                    OccEnergies(k)=0.0_dp
                enddo

                do i=(NEl/2)+1,NoRotAlphBet
                    CoeffT1(:,i)=NatOrbMat(:,i)
                    EvaluesTrunc(i)=Evalues(i)
                    SymOrbs_rot(i)=SymOrbs_rotTemp(i)
                enddo
            ENDIF

!            WRITE(6,*) SymOrbs(:)
!            CALL neci_flush(6)
!            CALL Stop_All('','')

        ELSE


            do i=1,NoOrbs
                CoeffT1(:,i)=NatOrbMat(:,i)
            enddo

        ENDIF

!        do i=1,NoOrbs
!            WRITE(6,*) Evalues(i)
!        enddo
!        do i=1,NoOrbs
!            WRITE(6,*) NatOrbMat(:,i)
!        enddo
!        CALL neci_flush(6)
!        stop

        IF(tTruncRODump) THEN

            Label=''
            LabelFull=''
            WRITE(Label,'(I5)') NoFrozenVirt
            LabelFull='EVALUES-TRUNC-'//adjustl(Label)

            io1 = get_free_unit()
            OPEN(io1,FILE=LabelFull,status='unknown')
            IF(tStoreSpinOrbs) THEN
                WRITE(io1,*) NoOrbs-NoFrozenVirt
                do i=1,NoOrbs-NoFrozenVirt,2
                    WRITE(io1,'(I5,ES20.10,I5,A5,I5,ES20.10,I5)') i,EvaluesTrunc(i),SymOrbs_rot(i),&
                        '  *  ',i+1, EvaluesTrunc(i+1),SymOrbs_rot(i+1)
                enddo
            ELSE
                WRITE(io1,*) NoOrbs-NoFrozenVirt
                do i=1,NoOrbs-NoFrozenVirt
                    WRITE(io1,'(ES20.10,I5)') EvaluesTrunc(i),SymOrbs_rot(i)
                enddo
            ENDIF
            CLOSE(io1) 
        ELSE
            io2 = get_free_unit()
            OPEN(io2,FILE='EVALUES',status='unknown')
            WRITE(io2,*) NoOrbs
            IF(tStoreSpinOrbs) THEN
                k=0
                do i=1,NoOrbs,2
                    k=k+1
                    IF(tTruncRODump) THEN
                        WRITE(io2,'(2I5,ES20.10,I5,A5,I5,ES20.10,I5)') (NoOrbs-i+1),i,Evalues(k),SymOrbs_rot(i),&
                            '  *  ', i+1,Evalues(k+SpatOrbs),SymOrbs_rot(i+1)
                    ELSE
                        WRITE(io2,'(2I5,ES20.10,I5,A5,I5,ES20.10,I5)') (NoOrbs-i+1),i,Evalues(k), &
                                INT(G1(SymLabelList3_rot(k))%Sym%S,4),'  *  ',&
                             i+1,Evalues(k+SpatOrbs),INT(G1(SymLabelList3_rot(k+SpatOrbs))%Sym%S,4)
                    ENDIF
                enddo
            ELSE
                do i=1,SpatOrbs
                    WRITE(io2,'(3I5,ES20.10)') i,NoOrbs-i+1,(NoOrbs-i+1)*2,Evalues(i)
                enddo
            ENDIF
            CLOSE(io2)
        ENDIF

        CALL HistNatOrbEvalues()

!        WRITE(6,*) 'NatOrbMat matrix'
!        do i=1,NoOrbs
!            WRITE(6,*) NatOrbMat(:,i)
!        enddo

!        OPEN(io1,FILE='TRANSFORMMAT',status='unknown')
!        do i=1,NoOrbs
!            do j=1,NoOrbs-NoFrozenVirt
!                WRITE(io1,*) i,j,CoeffT1(i,j)
!            enddo
!        enddo
!        CLOSE(io1)

        CALL halt_timer(FillCoeff_Time)
 

    ENDSUBROUTINE FillCoeffT1




    SUBROUTINE HistNatOrbEvalues()
        USE Logging , only : tTruncRODump
        USE RotateOrbsData , only : CoeffT1,EvaluesTrunc
        IMPLICIT NONE
        INTEGER :: i,k,a,b,NoOcc,io1, io2
        real(dp) :: EvalueEnergies(1:NoOrbs),OrbEnergies(1:NoOrbs)
        real(dp) :: SumEvalues
        
        io1 = get_free_unit()
        NoOcc = NEl/2   !Is this correct in all cases?!

        OPEN(io1,FILE='EVALUES-PLOTRAT',status='unknown')
        IF(tStoreSpinOrbs) THEN
            k=0
            do i=1,SpatOrbs
                k=k+2
                WRITE(io1,'(F20.10,ES20.10)') REAL(k-1)/REAL(NoOrbs),Evalues(i)
                WRITE(io1,'(F20.10,ES20.10)') REAL(k)/REAL(NoOrbs),Evalues(SpatOrbs+i)
            enddo
        ELSEIF(tRotateOccOnly) THEN
            k=0
            do i=1,NoOcc
                k=k+1
                WRITE(io1,'(F20.10,ES20.10)') REAL(k)/REAL(NoOcc),Evalues(i)
            enddo
        ELSEIF(tRotateVirtOnly) THEN
            k=NoOcc
            do i=NoOcc+1,NoOrbs
                k=k+1
                WRITE(io1,'(F20.10,ES20.10)') REAL(k-NoOcc)/REAL(NoOrbs-NoOcc),Evalues(i)
            enddo
        ELSE
            k=0
            do i=1,SpatOrbs
                k=k+1
                WRITE(io1,'(F20.10,ES20.10)') REAL(k)/REAL(NoOrbs),Evalues(i)
            enddo
        ENDIF
        CLOSE(io1)

!        OPEN(io2,FILE='EVALUES-plot',status='unknown')
!        EvaluesCount(:,:)=0.0_dp

!        do x=1,NoSpinCyc

!            IF(tSeparateOccVirt) THEN
!                IF(x.eq.1) THEN
!                    IF(tStoreSpinOrbs) THEN
!                        NoOcc=nOccBeta
!                    ELSE
!                        NoOcc=NEl/2
!                    ENDIF
!                ELSEIF(x.eq.2) THEN
!                    NoOcc=nOccAlpha
!                ENDIF
!            ELSE
!                NoOcc=0
!            ENDIF

!            k=1
!            EvaluesCount(k,1)=Evalues(1)
!            EvaluesCount(k,2)=1.0_dp
!            do i=2,NoOrbs
!                IF((ABS(Evalues(i)-Evalues(i-1))).ge.(1E-10)) THEN
!                    k=k+1
!                    EvaluesCount(k,1)=Evalues(i)
!                    EvaluesCount(k,2)=1.0_dp
!                ELSE
!                    EvaluesCount(k,2)=EvaluesCount(k,2)+1.0_dp
!                ENDIF
!            enddo
!            NoEvalues=k

!            do i=1,NoEvalues
!                WRITE(io2,*) EvaluesCount(i,1),Evaluescount(i,2)
!            enddo


!        enddo

!        CLOSE(io2)

! Want to write out the eigenvectors in order of the energy of the new orbitals - so that we can see 
! the occupations of the type of orbital.
! For now, keep this separate to the transformation of ARR - even though it is equivalent.

!        WRITE(6,*) 'ARR'
!        do i=1,NoOrbs
!            WRITE(6,*) ARR(2*i,2)
!        enddo
!        WRITE(6,*) 'Evalues'
!        do i=1,NoOrbs
!            WRITE(6,*) Evalues(i)
!        enddo

        OrbEnergies(:)=0.0_dp
        EvalueEnergies(:)=0.0_dp
        SumEvalues=0.0_dp
        do i=1,NoOrbs
            IF(tStoreSpinOrbs) THEN
                SumEvalues=SumEvalues+Evalues(i)
            ELSE
                SumEvalues=SumEvalues+(2*Evalues(i))
            ENDIF
            IF(tTruncRODump) THEN
                EvalueEnergies(i)=EvaluesTrunc(i)
            ELSE
                EvalueEnergies(i)=Evalues(i)
            ENDIF
! We are only interested in the diagonal elements.            
            do a=1,NoOrbs
                b=SymLabelList2_rot(a)
                IF(tStoreSpinOrbs) THEN
                    OrbEnergies(i)=OrbEnergies(i)+(CoeffT1(a,i)*ARR(b,2)*CoeffT1(a,i))
                ELSE
                    OrbEnergies(i)=OrbEnergies(i)+(CoeffT1(a,i)*ARR(2*b,2)*CoeffT1(a,i))
                ENDIF
            enddo
        enddo
! ARR(:,1) - ordered by energy, ARR(:,2) - ordered by spin-orbital index.

!        WRITE(6,*) 'OrbEnergies'
!        do i=1,NoOrbs
!            WRITE(6,*) OrbEnergies(i)
!        enddo
!        WRITE(6,*) 'EvalueEnergies'
!        do i=1,NoOrbs
!            WRITE(6,*) EvalueEnergies(i)
!        enddo

!        CALL SortEvecbyEval(NoOrbs,OrbEnergies(1:NoOrbs),1,EvalueEnergies(1,1:NoOrbs))
        call sort (orbEnergies(1:noOrbs), EvalueEnergies(1:noOrbs))

        io2 = get_free_unit()
        OPEN(io2,FILE='EVALUES-ENERGY',status='unknown')
        do i=1,NoOrbs
            WRITE(io2,*) OrbEnergies(NoOrbs-i+1),EvalueEnergies(NoOrbs-i+1)
        enddo
        WRITE(io2,*) 'The sum of the occupation numbers (eigenvalues) = ',SumEvalues
        WRITE(io2,*) 'The number of electrons = ',NEl
        CALL neci_flush(io2)
        CLOSE(io2)
        CALL neci_flush(6)

        CALL PrintOccTable()


    END SUBROUTINE HistNatOrbEvalues


    SUBROUTINE CalcOccEnergies(OccEnergies)
        USE RotateOrbsData , only : CoeffT1,NoRotOrbs
        real(dp) :: OccEnergies(1:NoRotOrbs)
        INTEGER :: i,a,b,NoOcc,x,Prev,k

        OccEnergies(:)=0.0_dp
        IF(tStoreSpinOrbs) THEN
            do x=1,2
                IF(x.eq.1) THEN
                    NoOcc=nOccBeta
                    Prev=0
                    k=1
                ELSE
                    NoOcc=nOccAlpha
                    Prev=SpatOrbs
                    k=2
                ENDIF
                do i=1+Prev,NoOcc+Prev
! We are only interested in the diagonal elements.            
                    do a=1,NoOrbs
                        b=SymLabelList2_rot(a)
                        OccEnergies(k)=OccEnergies(k)+(NatOrbMat(a,i)*ARR(b,2)*NatOrbMat(a,i))
                    enddo
                    k=k+2
                enddo
            enddo
        ELSE
            NoOcc=NEl/2
            do i=1,NoOcc
! We are only interested in the diagonal elements.            
                do a=1,NoOrbs
                    b=SymLabelList2_rot(a)
                    OccEnergies(i)=OccEnergies(i)+(NatOrbMat(a,i)*ARR(2*b,2)*NatOrbMat(a,i))
                enddo
            enddo
        ENDIF

    END SUBROUTINE CalcOccEnergies



    SUBROUTINE PrintOccTable()
        USE Logging , only : tTruncRODump
        USE RotateOrbsData , only : CoeffT1,EvaluesTrunc
        USE SystemData , only : tUseHFOrbs
        INTEGER x,i,a,b, io2

        io2 = get_free_unit()
        OPEN(io2,FILE='OccupationTable',status='unknown')
        x=1
        do while (x.le.NoOrbs)
            WRITE(io2,'(A16,A5)',advance='no') 'HF Orb En    ','Sym'
            IF(.not.tUseHFOrbs) THEN
                do i=x,x+9
                    IF(i.gt.NoOrbs) THEN
                        WRITE(io2,*) ''
                        EXIT
                    ENDIF
                    IF(tTruncRODump) THEN
                        WRITE(io2,'(ES16.6)',advance='no') EvaluesTrunc(i)
                    ELSE
                        WRITE(io2,'(ES16.6)',advance='no') Evalues(i)
                    ENDIF
                enddo
            ENDIF
            WRITE(io2,*) ''

            do a=1,NoOrbs
                b=SymLabelListInv_rot(a)
                IF(tStoreSpinOrbs) THEN
                    WRITE(io2,'(F16.10,I5)',advance='no') ARR(a,1),INT(G1(a)%sym%S,4)
                ELSE
                    WRITE(io2,'(F16.10,I5)',advance='no') ARR(2*a,1),INT(G1(2*a)%sym%S,4)
                ENDIF
                do i=x,x+9
                    IF(i.gt.NoOrbs) THEN
                        WRITE(io2,*) ''
                        EXIT
                    ENDIF
!                    WRITE(io2,'(F16.10)',advance='no') NatOrbMat(b,i)
                    WRITE(io2,'(F16.10)',advance='no') CoeffT1(b,i)
                enddo
                WRITE(io2,*) ''
            enddo
            WRITE(io2,*) ''
            x=x+10
        enddo
        CALL neci_flush(io2)
        CLOSE(io2)
        CALL neci_flush(6)
        

    END SUBROUTINE PrintOccTable


    SUBROUTINE PrintOrbOccs(OrbOccs)
! This routine takes whatever orbital basis we're using and is called at the end of a spawn to find 
! the contribution of each orbital to the final wavefunction.    
! This is done by histogramming the determinant populations, and then running over these adding the 
! coefficients of each determinant to the orbitals occupied.
! This is essentially < Psi | a_p+ a_p | Psi > - the diagonal terms of the one electron reduced density matrix.
!        USE Logging , only : OrbOccs
        IMPLICIT NONE
        real(dp) :: Norm,OrbOccs(nBasis),AllOrbOccs(nBasis)
        INTEGER :: i,error, iunit
        LOGICAL :: tWarning

        AllOrbOccs = 0.0_dp

        call MPIReduce(OrbOccs,MPI_SUM,AllOrbOccs)

! Want to normalise the orbital contributions for convenience.        
        tWarning=.false.
        IF(iProcIndex.eq.0) THEN
            Norm=0.0_dp
            do i=1,nBasis
                Norm=Norm+AllOrbOccs(i)
                IF((AllOrbOccs(i).lt.0).or.(Norm.lt.0)) THEN
                    WRITE(6,*) 'WARNING: Integer overflow when calculating the orbital occupations.'
                    tWarning=.true.
                ENDIF
            enddo
            IF(Norm.ne.0.0_dp) THEN
                do i=1,nBasis
                    AllOrbOccs(i)=AllOrbOccs(i)/Norm
                enddo
            ENDIF

            iunit = get_free_unit()
            OPEN(iunit,FILE='ORBOCCUPATIONS',STATUS='UNKNOWN')
            WRITE(iunit,'(A15,A30)') '# Orbital no.','Normalised occupation'
            IF(tWarning) WRITE(iunit,*) 'WARNING: INTEGER OVERFLOW OCCURRED WHEN CALCULATING THESE OCCUPATIONS'
            do i=1,nBasis
                WRITE(iunit,'(I15,F30.10)') i,AllOrbOccs(i)
            enddo
            CLOSE(iunit)
        ENDIF

    END SUBROUTINE PrintOrbOccs

    SUBROUTINE PrintDoubUEGOccs(OrbOccs)
! Based on PrintOrbOccs (above), but for PrintDoubsUEG
! Histogram determinant populations for all doubles
! This hopefully prints it all out
        IMPLICIT NONE
        real(dp) :: Norm,OrbOccs(nEl,nEl,nBasis,4),AllOrbOccs(nEl,nEl,nBasis,4)
        INTEGER :: i,i2,i3,error, iunit
        LOGICAL :: tWarning

        AllOrbOccs = 0.0_dp

        call MPISum(OrbOccs,AllOrbOccs)
!#ifdef PARALLEL
!        CALL MPI_Reduce(OrbOccs,AllOrbOccs,nEl*nEl*nBasis*4,MPI_DOUBLE_PRECISION,& 
!                                    MPI_SUM,0,MPI_COMM_WORLD,error)
!#else
!        AllOrbOccs=OrbOccs
!#endif

! Want to normalise the orbital contributions for convenience.        
        tWarning=.false.
        IF(iProcIndex.eq.0) THEN

!            Norm=0.0_dp
!            do i2=1,nEl
!                do i3=1,nEl
!                    do i=1,nBasis
!                        Norm=Norm+AllOrbOccs(i2,i3,i,1)
!                        !No need for this test at the moment
!                        !IF((Norm.lt.0)) THEN
!                        !    WRITE(6,*) 'WARNING: Integer overflow when calculating &
!                                               & the orbital occupations.'
!                        !    tWarning=.true.
!                        !ENDIF
!                    enddo
!                enddo
!            enddo
!            IF(Norm.ne.0.0_dp) THEN
!                do i2=1,nEl
!                    do i3=1,nEl
!                        do i=1,nBasis
!                            AllOrbOccs(i2,i3,i,1)=AllOrbOccs(i2,i3,i,1)/Norm
!                        enddo
!                    enddo
!                enddo
!            ENDIF

            iunit = get_free_unit()
            OPEN(iunit,FILE='DOUBOCCUPATIONS',STATUS='UNKNOWN')
!            WRITE(iunit,'(A15,A30)') '# Orbital no.','Normalised occupation'
            IF(tWarning) WRITE(iunit,*) &
                        'WARNING: INTEGER OVERFLOW OCCURRED WHEN CALCULATING THESE OCCUPATIONS'
            do i2=1,nEl
                do i3=1,nEl
                    do i=1,nBasis
                        WRITE(iunit,'(I15,I15,I15,F30.10,F30.10,F30.10,F30.10)') i2,i3,i,&
                                            AllOrbOccs(i2,i3,i,1), AllOrbOccs(i2,i3,i,2),&
                                            AllOrbOccs(i2,i3,i,3),AllOrbOccs(i2,i3,i,4)
                    enddo
                enddo
            enddo
            CLOSE(iunit)
        ENDIF

    END SUBROUTINE PrintDoubUEGOccs


    SUBROUTINE DeallocateNatOrbs()
        USE Logging , only : tTruncRODump
        IMPLICIT NONE
        CHARACTER(len=*), PARAMETER :: this_routine='DeallocateNatOrbs'

! Deallocate the natural orbitals matrix.    

        IF(tTruncRODump) THEN
            DEALLOCATE(SymOrbs_rotTemp)
            CALL LogMemDeAlloc(this_routine,SymOrbs_rotTempTag)
        ENDIF
        DEALLOCATE(NatOrbMat)
        CALL LogMemDeAlloc(this_routine,NatOrbMatTag)
        DEALLOCATE(Evalues)
        CALL LogMemDeAlloc(this_routine,EvaluesTag)

    END SUBROUTINE DeallocateNatOrbs




!This file was primarily concerned with the creation of natural orbitals from a rotation 
!of the previous orbitals.
!The 1-electron Reduced density matrix was inputted, and the natural orbitals constructed. 
!From there, the 1 and 2 electron integrals were transformed and replaced into UMat.
    SUBROUTINE FindNatOrbsOld()
        IMPLICIT NONE
        INTEGER :: i,j, iunit

        iunit = get_free_unit()
        OPEN(iunit,FILE='ONEEL-RDM',STATUS='UNKNOWN')
        do i=1,nBasis
            do j=1,nBasis
                WRITE(iunit,"(F18.7)",advance='no') NatOrbMat(i,j)
            enddo
            WRITE(iunit,*) ""
        enddo
        CLOSE(iunit)

        CALL Stop_All('FindNatOrbsOld',&
            'This is the old routine for finding the natural orbitals - likely buggy.')

!First, diagonalize the 1-RDM...
        CALL Diag1RDMOld()

!Setup symmetry information needed...
!    CALL SetupSymNO()

!Find orbitals energy...
!    CALL FindOrbEnergies()

!Transform integrals...
!    CALL TransRotInts()

    END SUBROUTINE FindNatOrbsOld


    SUBROUTINE Diag1RDMOld()
       use MemoryManager, only: TagIntType
! The diagonalisation routine reorders the orbitals in such a way that the corresponding 
! orbital labels are lost.
! In order to keep the spin and spatial symmetries, each symmetry must be fed into the 
! diagonalisation routine separately.
! The best way to do this is to order the orbitals so that all the alpha orbitals follow all the beta 
! orbitals, with the occupied orbitals first, in terms of symmetry, and the virtual second, also 
! ordered by symmetry. This gives us flexibility w.r.t rotating only the occupied or only virtual and 
! looking at high spin states.
        IMPLICIT NONE
        real(dp) , ALLOCATABLE :: NOccNums(:),Work(:)
        INTEGER(TagIntType) :: nOccNumsTag=0,WorkTag=0
        INTEGER :: iErr,WorkSize,WorkCheck,i
        CHARACTER(len=*), PARAMETER :: this_routine='Diag1RDM'

        ALLOCATE(NOccNums(nBasis),stat=ierr)
        CALL LogMemAlloc('NOccNums',nBasis,8,this_routine,NOccNumsTag,iErr)

!Find desired optimal workspace size
        WorkSize=-1
        WorkCheck=3*nBasis+1
!    CALL DSYEV('V','U',nBasis,NatOrbMat,nBasis,NOccNums,WorkCheck,WorkSize,iErr)
!    IF(iErr.ne.0) THEN
!        CALL Stop_All(this_routine,"Error in finding scratch space size")
!    ENDIF
        WorkSize=WorkCheck

        WRITE(6,*) "Optimal size of scratch space for diagonalization = ",WorkCheck

        ALLOCATE(Work(WorkSize),stat=ierr)
        CALL LogMemAlloc('Work',WorkSize,8,this_routine,WorkTag,iErr)

        WRITE(6,"(A)",advance='no') "Diagonalizing 1-RDM to find approximate natural orbitals..."
        CALL neci_flush(6)

!Diagonalize... Matrix must be symmetric
        CALL DSYEV('V','U',nBasis,NatOrbMat,nBasis,NOccNums,Work,WorkSize,iErr)
        IF(iErr.eq.0) THEN
            WRITE(6,"(A)") "DONE!"
        ELSE
            WRITE(6,"(A)") "FAILED!"
            CALL Stop_All(this_routine,"Diagonalization of 1-RDM failed...")
        ENDIF

        DEALLOCATE(Work)
        CALL LogMemDealloc(this_routine,WorkTag)

        WRITE(6,*) "Occupation numbers of approximate natural spin-orbitals are: "
        do i=1,nBasis
            WRITE(6,*) i,NOccNums(i)    !Symmetry of orbital...?
        enddo
        CALL neci_flush(6)

    END SUBROUTINE Diag1RDMOld


END MODULE NatOrbsMod
