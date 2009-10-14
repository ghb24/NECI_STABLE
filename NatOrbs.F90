MODULE NatOrbsMod
! This file is primarily concerned with finding the one electron reduced density matrix, from the wavefunction
! constructed by a previous spawning calculation.
! Diagonalisation of this density matrix gives a set of eigenvectors which rotate the HF orbitals into the 
! CI natural orbitals within the given excitation level.
! Once these eigenvectors have been obtained, the relevant routines from RotateOrbs are called to transform the 
! integrals and produce a ROFCIDUMP file in the natural orbital basis.
        
        USE Global_utilities
        USE Parallel
        USE IntegralsData , only : UMAT
        USE UMatCache , only : UMatInd
        USE SystemData , only : NEl,nBasis,NIfD,G1,ARR,BRR,lNoSymmetry,LMS,tStoreSpinOrbs,nOccAlpha,nOccBeta,tSeparateOccVirt
        USE SystemData , only : tRotateOccOnly,tRotateVirtOnly,tFindCINatOrbs,tUseMP2VarDenMat
        USE RotateOrbsData , only : SymLabelList2,SymLabelCounts2,SymLabelCounts2Tag,SymLabelListInv,NoOrbs,SpatOrbs,FillOneRDM_time
        USE RotateOrbsData , only : FillMP2VDM_Time,DiagNatOrbMat_Time,OrderCoeff_Time,FillCoeff_Time,NoFrozenVirt
        IMPLICIT NONE
        INTEGER :: NoSpinCyc,SymOrbsTempTag
        REAL*8 , ALLOCATABLE :: NatOrbMat(:,:),Evalues(:)
        INTEGER , ALLOCATABLE :: SymOrbsTemp(:)
        INTEGER :: NatOrbMatTag,ierr,EvaluesTag

    contains
    
    SUBROUTINE FindNatOrbs()
        IMPLICIT NONE
        INTEGER :: ierr
        CHARACTER(len=*) , PARAMETER :: this_routine='FindNatOrbs'

! Fed into this routine will be the wavefunction, Psi, and its amplitudes within the given excitation level.    

! First need to set up the orbital labels and symmetries etc.
! This is done slightly differently for spin and spatial and whether or not we are truncating the virtual space when
! writing out the final ROFCIDUMP file.

! Allocate the matrix used to find the natural orbitals.
        
        ALLOCATE(NatOrbMat(NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('NatOrbMat',NoOrbs**2,8,this_routine,NatOrbMatTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for NatOrbMat failed.")
        NatOrbMat(:,:)=0.D0

        ALLOCATE(Evalues(NoOrbs),stat=ierr)
        CALL LogMemAlloc('Evalues',NoOrbs,8,this_routine,EvaluesTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for Evalues failed.")
        Evalues(:)=0.D0

! First need to fill the relevant matrix for calculating the type of natural orbitals we want.
        IF(tFindCINatOrbs) THEN

            ! For the CISD, CISDT etc natural orbitals, the relevant matrix is the one electron reduced density matrix
            ! from the previous spawning calculation (trucated at a certain excitation).
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
        IMPLICIT NONE
        INTEGER :: x,i,j,ierr,NoOcc,StartFill01,StartFill02,Symi,SymCurr,Prev,EndFill01,EndFill02
        CHARACTER(len=*) , PARAMETER :: this_routine='SetupNatOrbLabels'
        INTEGER , ALLOCATABLE :: LabVirtOrbs(:),LabOccOrbs(:),SymVirtOrbs(:),SymOccOrbs(:)
        INTEGER :: LabVirtOrbsTag,LabOccOrbsTag,SymVirtOrbsTag,SymOccOrbsTag
        

! The earlier test should pick this up, if it crashes here, will want to put in an earlier test so that we don't 
! get all the way to this stage.
        IF((LMS.ne.0).and.(.not.tStoreSpinOrbs)) CALL Stop_All("FindNatOrbs","Open shell system, and UMAT is not being &
                                                            &stored as spin orbitals.")

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

! First fill SymLabelList2.

! Brr has the orbital numbers in order of energy... i.e Brr(2) = the orbital index with the second lowest energy.

!            do i=1,(SpatOrbs*2)
!                WRITE(6,*) BRR(i)
!            enddo
!            CALL FLUSH(6)
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
            IF(tSeparateOccVirt) CALL NECI_SORT2I(NoOcc,SymOccOrbs,LabOccOrbs)
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

            CALL NECI_SORT2I((SpatOrbs-NoOcc),SymVirtOrbs,LabVirtOrbs)
            ! Sorts LabOrbs according to the order of SymOccOrbs (i.e. in terms of symmetry). 

!            WRITE(6,*) 'after sort'
!            do i=1,NoOrbs
!                WRITE(6,*) LabVirtOrbs(i),SymVirtOrbs(i)
!            enddo

!            CALL FLUSH(6)
!            stop

 
! SymLabelList2 is then filled with the symmetry ordered occupied then virtual arrays for each spin.        
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
                SymLabelList2(i)=LabOccOrbs(j)
            enddo
            j=0
            do i=StartFill02,EndFill02
                j=j+1
                SymLabelList2(i)=LabVirtOrbs(j)
            enddo


!************
! Second fill SymLabelCounts2.
! - the first 8 places of SymLabelCounts2(1,:) and SymLabelCounts2(2,:) refer to the occupied orbitals 
! - and the second 8 to the virtuals.

            IF(lNoSymmetry) THEN
                ! if we are ignoring symmetry, all orbitals essentially have symmetry 0.
                IF(x.eq.1) THEN
                    SymLabelCounts2(1,1)=1
                    SymLabelCounts2(1,9)=NoOcc+1
                    SymLabelCounts2(2,1)=NoOcc
                    SymLabelCounts2(2,9)=SpatOrbs-NoOcc
                ELSEIF(x.eq.2) THEN
                    SymLabelCounts2(1,17)=1
                    SymLabelCounts2(1,25)=NoOcc+1
                    SymLabelCounts2(2,17)=NoOcc
                    SymLabelCounts2(2,25)=SpatOrbs-NoOcc
                ENDIF
                
            ELSE 
                ! otherwise we run through the occupied orbitals, counting the number with each symmetry
                ! and noting where in SymLabelList2 each symmetry block starts.
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
                SymLabelCounts2(1,StartFill01)=1+Prev
                do i=1,NoOcc
                    IF(tStoreSpinOrbs) THEN
                        Symi=INT(G1(SymLabelList2(i+Prev))%sym%S,4)
                    ELSE
                        Symi=INT(G1((SymLabelList2(i+Prev)*2))%sym%S,4)
                    ENDIF
                    SymLabelCounts2(2,(Symi+StartFill01))=SymLabelCounts2(2,(Symi+StartFill01))+1
                    IF(Symi.ne.SymCurr) THEN
                        SymLabelCounts2(1,(Symi+StartFill01))=i+Prev
                        SymCurr=Symi
                    ENDIF
                enddo
                ! the same is then done for the virtuals.
                SymCurr=0
                SymLabelCounts2(1,StartFill02)=NoOcc+1+Prev
                do i=NoOcc+1,SpatOrbs
                    IF(tStoreSpinOrbs) THEN
                        Symi=INT(G1(SymLabelList2(i+Prev))%sym%S,4)
                    ELSE
                        Symi=INT(G1((SymLabelList2(i+Prev)*2))%sym%S,4)
                    ENDIF
                    SymLabelCounts2(2,(Symi+StartFill02))=SymLabelCounts2(2,(Symi+StartFill02))+1
                    IF(Symi.ne.SymCurr) THEN
                        SymLabelCounts2(1,(Symi+StartFill02))=i+Prev
                        SymCurr=Symi
                    ENDIF
                enddo
            ENDIF
     
            ! Go through each symmetry group, making sure the orbital pairs are ordered lowest to highest.
            IF(x.eq.1) THEN
                do i=1,16
                    IF(SymLabelCounts2(2,i).ne.0) THEN
                        CALL NECI_SORTI(SymLabelCounts2(2,i),SymLabelList2(SymLabelCounts2(1,i):(SymLabelCounts2(1,i)+SymLabelCounts2(2,i)-1)))
                    ENDIF
                enddo
            ELSEIF(x.eq.2) THEN
                do i=17,32
                    IF(SymLabelCounts2(2,i).ne.0) THEN
                        CALL NECI_SORTI(SymLabelCounts2(2,i),SymLabelList2(SymLabelCounts2(1,i):(SymLabelCounts2(1,i)+SymLabelCounts2(2,i)-1)))
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
            SymLabelListInv(SymLabelList2(i))=i
        enddo

        IF(.not.tSeparateOccVirt) THEN
            ! basically we treat all the orbitals as virtuals and set NoOcc to zero in each routine. 
            tRotateVirtOnly=.true.
        ENDIF
        

!        WRITE(6,*) 'Sym Label Counts'
!        do i=1,32
!            WRITE(6,*) i,SymLabelCounts2(1,i),SymLabelCounts2(2,i)
!        enddo
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and their symmetries according to G1'
!        do i=1,NoOrbs
!            IF(tStoreSpinOrbs) THEN
!                WRITE(6,*) i,SymLabelList2(i),INT(G1(SymLabelList2(i))%sym%S,4)
!            ELSE
!                WRITE(6,*) i,SymLabelList2(i),INT(G1(SymLabelList2(i)*2)%sym%S,4)
!            ENDIF
!        enddo
!        WRITE(6,*) 'i','ARR(SymLabelList2(i),1)','ARR(SymLabelList2(i),2)','Sym'
!        do i=1,NoOrbs
!            IF(tStoreSpinOrbs) THEN
!                WRITE(6,*) i,ARR(SymLabelList2(i),1),ARR(SymLabelList2(i),2),INT(G1(SymLabelList2(i))%sym%S,4)
!            ENDIF
!        enddo
!
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and its inverse'
!        do i=1,NoOrbs
!            WRITE(6,*) SymLabelList2(i),SymLabelListInv(i)
!        enddo
!        CALL FLUSH(6)
!        CALL Stop_All('SetupNatOrbLabels','Checking orbital labelling.')


    END SUBROUTINE SetupNatOrbLabels



    SUBROUTINE FillOneRDM()
        USE DetCalc , only : Det,FCIDets,FCIDetIndex,ICILevel
! Det is the number of determinants in FCIDets.
! FCIDets contains the list of all determinants in the system in bit string representation, FCIDets(0:nBasis/32,1:Det) 
! ICILevel is the max excitation level of the calculation - as in EXCITE ICILevel.
! FCIDetIndex(1:NEl) contains the index of FCIDets where each excitation level starts.
! As in FCIDetIndex(1) = 2 always I think - Excitation level 1 starts at the second determinant (after HF).
! Pretty sure FCIDetIndex always goes from 1:NEl even from truncated excite calculations.
        USE FciMCData , only : AllHistogram
! The elements of AllHistogram correspond to the rows of FCIDets - i.e to each determinant in the system.
! AllHistogram contains the final (normalised) amplitude of the determinant - with sign.
        IMPLICIT NONE
        INTEGER :: x,excit,i,j,NoOcc,Starti,Endi,Startj,Endj,ExcitLevel,Ex(2,1),Ex2(2,1),Orbi,Orbj,nJ(NEl),Orbk,k,nI(NEl),MaxExcit
        INTEGER :: FCIDetIndex2(0:(NEl+1)),Spins
        LOGICAL :: tSign
        REAL*8 :: SignDet

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
!            WRITE(6,*) FCIDets(0:NIfD,i),AllHistogram(i)
!        enddo
!        CALL FLUSH(6)
!        stop
        WRITE(6,*) '*** The weight of the HF determinant is : ', AllHistogram(1)

        WRITE(6,*) 'Beginning to fill the one-electron reduced density matrix.'

        IF(ICILevel.eq.0) THEN
            MaxExcit=NEl
        ELSE
            MaxExcit=ICILevel
        ENDIF

        do excit=0,MaxExcit         
        ! Run through all determinants D_i, in the final wavefunction, Psi. 
        ! If this is done by excitation block, we then don't have to check the excitation level of the determinant each time.
            IF(tRotateVirtOnly.and.tSeparateOccVirt.and.(excit.eq.0)) CYCLE      ! The HF only involves 'occupied' orbitals - these are not required if only rotating virt.

! This next bit is a bit messy because there is no row in FCIDetIndex for the HF - there is probably an tidier way to achieve the same thing, but it does the trick for now.
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
                ! Run through all determinants D_j, with the potential to be connected to i by a single excitation, i.e from one excitation
!               ! lower to one excitation higher.
                    IF((i.gt.Det).or.(j.gt.Det)) CALL Stop_All('FillOneRDM','Running through i or j larger than the number of determinants.')

                    CALL FindBitExcitLevel(FCIDets(0:NIfD,i),FCIDets(0:NIfD,j),NIfD,ExcitLevel,2)
                    ! Need to find the excitation level between D_i and D_j. If this is 1 - go on to add their contributions to the OneRDM.

                    IF(ExcitLevel.eq.1) THEN
                        Ex(:,:)=0
                        Ex(1,1)=ExcitLevel

                        CALL GetBitExcitation(FCIDets(0:NIfD,i),FCIDets(0:NIfD,j),NIfD,NEl,Ex,tSign)
                        ! Gives the orbitals involved in the excitation Ex(1,1) in i -> Ex(2,1) in j (in spin orbitals).

                        IF(tStoreSpinOrbs) THEN
                            ! OneRDM will be in spin orbitals - simply add the orbitals involved.
                            Orbi=SymLabelListInv(Ex(1,1))
                            Orbj=SymLabelListInv(Ex(2,1))
                            Spins=1
                        ELSE
                            Orbi=SymLabelListInv(CEILING(REAL(Ex(1,1))/2.D0))
                            Orbj=SymLabelListInv(CEILING(REAL(Ex(2,1))/2.D0))
                            Spins=2
                        ENDIF
                        IF(tSign) THEN
                            SignDet=(-1.D0)
                        ELSE
                            SignDet=1.D0
                        ENDIF

                        NatOrbMat(Orbi,Orbj)=NatOrbMat(Orbi,Orbj)+(SignDet*AllHistogram(i)*AllHistogram(j))
                        NatOrbMat(Orbj,Orbi)=NatOrbMat(Orbj,Orbi)+(SignDet*AllHistogram(i)*AllHistogram(j))

                        ! AllHistogram are the normalised amplitudes of the determinants.
!                        IF(((AllHistogram(i)*AllHistogram(j).ne.0.D0).and.(INT(G1(SymLabelList2(Orbi)*2)%sym%S,4).ne.INT(G1(SymLabelList2(Orbj)*2)%sym%S,4)))&
!                        &.or.(Ex(1,1).gt.(SpatOrbs*2)).or.(Ex(2,1).gt.(SpatOrbs*2))) THEN

                        IF((AllHistogram(i)*AllHistogram(j).ne.0.D0).and.(INT(G1(SymLabelList2(Orbi)*Spins)%sym%S,4).ne.INT(G1(SymLabelList2(Orbj)*Spins)%sym%S,4))) THEN
                            WRITE(6,*) 'ERROR in symmetries'
                            WRITE(6,*) 'Ex,',Ex(1,1),Ex(2,1)
                            WRITE(6,*) CEILING(REAL(Ex(1,1)/2.D0)),CEILING(REAL(Ex(2,1)/2.D0))
                            WRITE(6,*) 'Orbi,',Orbi,'Orbj,',Orbj
                            WRITE(6,*) 'Sym(Orbi)',INT(G1(SymLabelList2(Orbi)*Spins)%sym%S,4),'Sym(Orbj)',INT(G1(SymLabelList2(Orbj)*Spins)%sym%S,4)
                            CALL DecodeBitDet(nI,FCIDets(0:NIfD,i),NEl,NIfD)
                            WRITE(6,*) 'i',nI
                            CALL DecodeBitDet(nJ,FCIDets(0:NIfD,j),NEl,NIfD)
                            WRITE(6,*) 'j',nJ
                            WRITE(6,*) 'AllHistogram(i)',AllHistogram(i)
                            WRITE(6,*) 'AllHistogram(j)',AllHistogram(j)
                            CALL FLUSH(6)
                            CALL Stop_All('FillOneRDM','Non-zero element between different symmetries.')
                        ENDIF

                    ELSEIF(ExcitLevel.eq.0) THEN
                        CALL DecodeBitDet(nJ,FCIDets(0:NIfD,j),NEl,NIfD)
                        do k=1,NEl
!                            WRITE(6,*) 'k',k
                            IF(tStoreSpinOrbs) THEN
                                Orbk=SymLabelListInv(nJ(k))
                            ELSE
                                Orbk=SymLabelListInv(CEILING(REAL(nJ(k))/2.D0))
                            ENDIF
                            NatOrbMat(Orbk,Orbk)=NatOrbMat(Orbk,Orbk)+(AllHistogram(j)**2)
!                            NatOrbMat(Orbk,Orbk)=NatOrbMat(Orbk,Orbk)+(0.5 * (AllHistogram(j)**2))
                            ! 0.5 x because this will be added twice since we are not currently restricting i<k or anything.
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
!        CALL FLUSH(6)
!        stop

        CALL halt_timer(FillOneRDM_Time)



    END SUBROUTINE FillOneRDM



    SUBROUTINE FillMP2VDM()
! In this routine, the natural orbital matrix is calculated from the MP2 variational density matrix.

! MP2VDM = D2_ab = sum_ijc [ t_ij^ac ( 2 t_ij^bc - t_ji^bc ) ]
! Where :  t_ij^ac = - < ab | ij > / ( E_a - E_i + E_b - Ej )
! Ref : J. Chem. Phys. 131, 034113 (2009) - note: in Eqn 1, the cb indices are the wrong way round (should be bc).
        INTEGER :: a,b,c,i,j,a2,b2,c2,i2,j2,x,y,z,w
        INTEGER :: Startab,Endab,NoOcc,NoOccC,Startc,Endc,Starti,Endi,Startj,Endj
        REAL*8 :: MP2VDMSum
        CHARACTER(len=*), PARAMETER :: this_routine='FillMP2VDM'

! Calculating the MP2VDM (D2_ab) matrix whose eigenvectors become the transformation matrix.        
! This goes in the natural orbital matrix of this module.
! The eigenvalues are the occupation numbers of the new orbitals.  These should decrease exponentially so that when we remove the 
! orbitals with small occupation numbers we should have little affect on the energy.


! For the MP2VDM, we always only rotate the virtual orbitals - denomonator term of the above expression would be 0 if a and b were occupied.
! The orbital labels are ordered occupied then virtual if spatial orbitals are being used,
! otherwise they go occupied beta, virtual beta, occupied alpha, virtual alpha.
! This is so the alpha and beta spins can be diagonalised separately and we can keep track of which is which when the evectors are reordered 
! and maintain spin symmetry.


        WRITE(6,*) 'Filling MP2VDM nat orb matrix'
        CALL FLUSH(6)
        
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
                a=SymLabelList2(a2)
!                do b2=Startab,Endab
                do b2=Startab,a2

                    b=SymLabelList2(b2)

                    MP2VDMSum=0.D0

                    ! when a and b beta, run over both alpha and beta virtual for c, then both alpha and beta virtual for both i and j etc. 

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
                            c=SymLabelList2(c2)

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
                                    i=SymLabelList2(i2)

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
                                            j=SymLabelList2(j2)

                                            IF(tStoreSpinOrbs) THEN
                                                IF((ARR(i,2)+ARR(j,2)-ARR(a,2)-ARR(c,2)).eq.0.D0) THEN
                                                    IF((REAL(UMAT(UMatInd(a,c,i,j,0,0))%v,8)).ne.0.D0) THEN
                                                        WRITE(6,*) i,j,a,c,REAL(UMAT(UMatInd(a,c,i,j,0,0))%v,8)
                                                        CALL Stop_All(this_routine,"Dividing a non-zero by zero.")
                                                    ENDIF
                                                ENDIF
                                                MP2VDMSum=MP2VDMSum+&
                                                            &(( (REAL(UMAT(UMatInd(a,c,i,j,0,0))%v,8)) * (2.D0*(REAL(UMAT(UMatInd(b,c,i,j,0,0))%v,8))) )/&
                                                            &( (ARR(i,2)+ARR(j,2)-ARR(a,2)-ARR(c,2)) * (ARR(i,2)+ARR(j,2)-ARR(b,2)-ARR(c,2)) ) )
                                                MP2VDMSum=MP2VDMSum-&                                            
                                                            &(( (REAL(UMAT(UMatInd(a,c,i,j,0,0))%v,8)) * (REAL(UMAT(UMatInd(c,b,i,j,0,0))%v,8)) )/&
                                                            &( (ARR(i,2)+ARR(j,2)-ARR(a,2)-ARR(c,2)) * (ARR(i,2)+ARR(j,2)-ARR(c,2)-ARR(b,2)) ) )
                 
                                            ELSE
                                                MP2VDMSum=MP2VDMSum+&
                                                            &(( (REAL(UMAT(UMatInd(a,c,i,j,0,0))%v,8)) * (2.D0*(REAL(UMAT(UMatInd(b,c,i,j,0,0))%v,8))) )/&
                                                            &( (ARR(2*i,2)+ARR(2*j,2)-ARR(2*a,2)-ARR(2*c,2)) * (ARR(2*i,2)+ARR(2*j,2)-ARR(2*b,2)-ARR(2*c,2)) ) )
                                                MP2VDMSum=MP2VDMSum-&                                            
                                                            &(( (REAL(UMAT(UMatInd(a,c,i,j,0,0))%v,8)) * (REAL(UMAT(UMatInd(c,b,i,j,0,0))%v,8)) )/&
                                                            &( (ARR(2*i,2)+ARR(2*j,2)-ARR(2*a,2)-ARR(2*c,2)) * (ARR(2*i,2)+ARR(2*j,2)-ARR(2*c,2)-ARR(2*b,2)) ) )
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
!        CALL FLUSH(6)
!        CALL Stop_All('','')

        CALL halt_timer(FillMP2VDM_Time)


    END SUBROUTINE FillMP2VDM




    SUBROUTINE DiagNatOrbMat()
! The diagonalisation routine reorders the orbitals in such a way that the corresponding orbital labels are lost.
! In order to keep the spin and spatial symmetries, each symmetry must be fed into the diagonalisation routine separately.
! The best way to do this is to order the orbitals so that all the alpha orbitals follow all the beta orbitals, with the 
! occupied orbitals first, in terms of symmetry, and the virtual second, also ordered by symmetry.
! This gives us flexibility w.r.t rotating only the occupied or only virtual and looking at high spin states.
        IMPLICIT NONE
        REAL*8 :: SumTrace,SumDiagTrace
        REAL*8 , ALLOCATABLE :: Work(:),WORK2(:),EvaluesSym(:),NOMSym(:,:)
        INTEGER :: ierr,i,j,x,z,Sym,LWORK2,WORK2Tag,SymStartInd,NoSymBlock,PrevSym,StartOccVirt,EndOccVirt,Prev,NoOcc
        INTEGER :: EvaluesSymTag,NOMSymTag
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
                        NatOrbMat(i+Prev,j+Prev)=0.D0
                        NatOrbMat(j+Prev,i+Prev)=0.D0
                        IF(i.eq.j) NatOrbMat(i+Prev,j+Prev)=1.D0
                    enddo
                    Evalues(i+Prev)=1.D0
                enddo
            ELSEIF(tRotateOccOnly) THEN
                do i=NoOcc+1,SpatOrbs
                    do j=1,SpatOrbs
                        NatOrbMat(i+Prev,j+Prev)=0.D0
                        NatOrbMat(j+Prev,i+Prev)=0.D0
                        IF(i.eq.j) NatOrbMat(i+Prev,j+Prev)=1.D0
                    enddo
                    Evalues(i+Prev)=1.D0
                enddo
            ELSEIF(tSeparateOccVirt) THEN
                do i=1,NoOcc
                    do j=NoOcc+1,SpatOrbs
                        NatOrbMat(i+Prev,j+Prev)=0.D0
                        NatOrbMat(j+Prev,i+Prev)=0.D0
                    enddo
                enddo
            ENDIF
        enddo

! Test that we're not breaking symmetry.
        do i=1,NoOrbs
            do j=1,NoOrbs
                IF(tStoreSpinOrbs) THEN
!                    WRITE(6,*) INT(G1(SymLabelList2(i))%sym%S,4),INT(G1(SymLabelList2(j))%sym%S,4),NatOrbMat(i,j)
                    IF((INT(G1(SymLabelList2(i))%sym%S,4).ne.INT(G1(SymLabelList2(j))%sym%S,4))) THEN
                        IF(ABS(NatOrbMat(i,j)).ge.1.0E-15) THEN
                            WRITE(6,'(6A8,A20)') 'i','j','Label i','Label j','Sym i','Sym j','Matrix value'
                            WRITE(6,'(6I3,F40.20)') i,j,SymLabelList2(i),SymLabelList2(j),INT(G1(SymLabelList2(i))%sym%S,4),INT(G1(SymLabelList2(j))%sym%S,4),NatOrbMat(i,j)
                            IF(tUseMP2VarDenMat) THEN
                                WRITE(6,*) '**WARNING** - There is a non-zero NatOrbMat value between orbitals of different symmetry.'
                                WRITE(6,*) 'These elements will be ignored, and the symmetry maintained in the final transformation matrix.'
                            ELSE
                                CALL Stop_All(this_routine,'Non-zero NatOrbMat value between different symmetries.')
                            ENDIF
                        ENDIF
                        NatOrbMat(i,j)=0.D0
                    ENDIF
                ELSE
!                    WRITE(6,*) INT(G1(SymLabelList2(i)*2)%sym%S,4),INT(G1(SymLabelList2(j)*2)%sym%S,4),NatOrbMat(i,j)
                    IF((INT(G1(SymLabelList2(i)*2)%sym%S,4).ne.INT(G1(SymLabelList2(j)*2)%sym%S,4))) THEN
                        IF(ABS(NatOrbMat(i,j)).ge.1.0E-15) THEN
                            WRITE(6,'(6A8,A20)') 'i','j','Label i','Label j','Sym i','Sym j','Matrix value'
                            WRITE(6,'(6I3,F40.20)') i,j,SymLabelList2(i),SymLabelList2(j),INT(G1(SymLabelList2(i)*2)%sym%S,4),INT(G1(SymLabelList2(j)*2)%sym%S,4),NatOrbMat(i,j)
                            CALL Stop_All(this_routine,'Non-zero NatOrbMat value between different symmetries.')
                        ENDIF
                        NatOrbMat(i,j)=0.D0
                    ENDIF
                ENDIF
            enddo
        enddo

        SumTrace=0.D0
        do i=1,NoOrbs
            SumTrace=SumTrace+NatOrbMat(i,i)
        enddo

        WRITE(6,*) 'Calculating eigenvectors and eigenvalues of NatOrbMat'
        CALL FLUSH(6)

        ! If we are using spin orbitals, need to feed in the alpha and beta spins separately.
        ! Otherwise these jumble up and the final ordering is uncorrect. 
        ! There should be no non-zero values between these, but can put a check in for this.

        do x=1,NoSpinCyc

! If we want to maintain the symmetry, we cannot have all the orbitals jumbled up when the diagonaliser reorders the eigenvectors.
! Must instead feed each symmetry block in separately.
! This means that although the transformed orbitals are jumbled within the symmetry blocks, the symmetry labels are all that are relevant and these are unaffected.
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

                    NoSymBlock=SymLabelCounts2(2,Sym+PrevSym)

                    SymStartInd=SymLabelCounts2(1,Sym+PrevSym)-1
                    ! This is one less than the index that the symmetry starts, so that when we run through i=1,..., we can
                    ! start at SymStartInd+i.

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

                        ! CAREFUL if eigenvalues are put in ascending order, this may not be correct, with the labelling system.
                        ! may be better to just take coefficients and transform TMAT2DRot in transform2elints.
                        ! a check that comes out as diagonal is a check of this routine anyway.

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
                        NatOrbMat(SymStartInd+1,SymStartInd+1)=1.D0
                        WRITE(6,*) '*****'
                        WRITE(6,*) 'Symmetry ',Sym,' has only one orbital.'
                        WRITE(6,*) 'Copying diagonal element ,',SymStartInd+1,'to NatOrbMat'
                    ENDIF

                    Sym=Sym+1
                enddo
            enddo
        enddo

        WRITE(6,*) 'Matrix diagonalised'
        CALL FLUSH(6)

        SumDiagTrace=0.D0
        do i=1,NoOrbs
            SumDiagTrace=SumDiagTrace+Evalues(i)
        enddo
        IF((ABS(SumDiagTrace-SumTrace)).gt.1E-10) THEN
            WRITE(6,*) 'Sum of diagonal NatOrbMat elements : ',SumTrace
            WRITE(6,*) 'Sum of eigenvalues : ',SumDiagTrace
            CALL Stop_All(this_routine,'The trace of the 1RDM matrix before diagonalisation is not equal to that after.')
        ENDIF

        CALL halt_timer(DiagNatOrbMat_Time)

    END SUBROUTINE DiagNatOrbMat



    SUBROUTINE OrderCoeffT1()
        USE RotateOrbsData , only : SymLabelList3
        USE Logging , only : tTruncRODump
        IMPLICIT NONE
        INTEGER :: x,i,ier,ierr,StartSort,EndSort,NoOcc
        CHARACTER(len=*), PARAMETER :: this_routine='OrderCoeffT1'
        

! Here, if symmetry is kept, we are going to have to reorder the eigenvectors according to the size of the eigenvalues, while taking
! the orbital labels (and therefore symmetries) with them. This will be put back into MP2VDM from MP2VDMTemp.

! Want to reorder the eigenvalues from largest to smallest, taking the eigenvectors with them and the symmetry as well.  
! If using spin orbitals, do this for the alpha spin and then the beta.
 
        OrderCoeff_Time%timer_name='OrderCoeff'
        CALL set_timer(OrderCoeff_Time,30)


        IF(tTruncRODump) THEN
            ! If we are truncating, the orbitals stay in this order, so we want to take their symmetries with them.
            ALLOCATE(SymOrbsTemp(NoOrbs),stat=ierr)
            CALL LogMemAlloc('SymOrbsTemp',NoOrbs,4,this_routine,SymOrbsTempTag,ierr)
            SymOrbsTemp(:)=0

            IF(tStoreSpinOrbs) THEN
                do i=1,NoOrbs
                    SymOrbsTemp(i)=INT(G1(SymLabelList2(i))%sym%S,4)
                enddo
            ELSE 
                do i=1,NoOrbs
                    SymOrbsTemp(i)=INT(G1(SymLabelList2(i)*2)%sym%S,4)
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

                CALL SortEvecbyEvalPlus1(((EndSort-StartSort)+1),Evalues(StartSort:EndSort),((EndSort-StartSort)+1),NatOrbMat(StartSort:EndSort,&
                                            &StartSort:EndSort),SymOrbsTemp(StartSort:EndSort))

            enddo
               
        ELSE
            ! If we are not truncating, the orbitals get put back into their original order, so the symmetry information is still 
            ! correct, no need for the SymOrbs array.
            ! Instead, just take the labels of SymLabelList3 with them.

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

                CALL SortEvecbyEvalPlus1(((EndSort-StartSort)+1),Evalues(StartSort:EndSort),((EndSort-StartSort)+1),&
                                            &NatOrbMat(StartSort:EndSort,StartSort:EndSort),SymLabelList3(StartSort:EndSort))
            enddo 
            
        ENDIF

        CALL halt_timer(OrderCoeff_Time)

        WRITE(6,*) 'Eigen-values: '
        do i=1,NoOrbs
            WRITE(6,*) Evalues(i)
        enddo
       

    END SUBROUTINE OrderCoeffT1



    SUBROUTINE FillCoeffT1
        USE RotateOrbsData , only : CoeffT1,SymLabelList3,SymOrbs,SymOrbsTag,TruncEval,NoRotOrbs
        USE Logging , only : tTruncRODump,tTruncDumpbyVal
        IMPLICIT NONE
        INTEGER :: k,i,j,NoRotAlphBet,SymFirst
        REAL*8 :: EvaluesTrunc(NoOrbs-NoFrozenVirt)
        CHARACTER(len=*), PARAMETER :: this_routine='FillCoeffT1'
        CHARACTER(len=5) :: Label
        CHARACTER(len=20) :: LabelFull
        LOGICAL :: tSymFound
  
        FillCoeff_Time%timer_name='FillCoeff'
        CALL set_timer(FillCoeff_Time,30)

        IF(tTruncRODump) THEN

            IF(tTruncDumpbyVal) THEN
                NoFrozenVirt=0
                IF(tStoreSpinOrbs) THEN
                    do i=SpatOrbs,1,-1
                        IF(Evalues(i).gt.TruncEval) EXIT
                        IF(Evalues(i+SpatOrbs).gt.TruncEval) EXIT
                        NoFrozenVirt=NoFrozenVirt+2
                    enddo
                    IF(NoFrozenVirt.ge.(NoOrbs-NEl)) CALL Stop_All(this_routine,'Freezing all virtual orbitals.')
                ELSE
                    do i=SpatOrbs,1,-1
                        IF(Evalues(i).gt.TruncEval) EXIT
                        NoFrozenVirt=NoFrozenVirt+1
                    enddo
                    IF(NoFrozenVirt.ge.(SpatOrbs-(NEl/2))) CALL Stop_All(this_routine,'Freezing all virtual orbitals.')
                ENDIF
                NoRotOrbs=NoOrbs-NoFrozenVirt
            ENDIF

            ALLOCATE(SymOrbs(NoOrbs),stat=ierr)
            CALL LogMemAlloc('SymOrbs',NoOrbs,4,this_routine,SymOrbsTag,ierr)
            SymOrbs(:)=0

            EvaluesTrunc(:)=0.D0

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
!                        SymFirst=SymOrbsTemp(i)
!                        CoeffT1(:,k)=NatOrbMat(:,i)
!                        EvaluesTrunc(k)=Evalues(i)
!                        SymOrbs(k)=SymOrbsTemp(i)
!                        do j=SpatOrbs+1,SpatOrbs+NoRotAlphBet
!                            IF(SymOrbsTemp(j).eq.SymFirst) THEN
!                                SymOrbs(k+1)=SymOrbsTemp(j)
!                                CoeffT1(:,k+1)=NatOrbMat(:,j)
!                                EvaluesTrunc(k+1)=Evalues(j)
!                                SymOrbsTemp(j)=9
!                                tSymFound=.true.
!                                EXIT
!                            ENDIF
!                        enddo
!                        IF(.not.tSymFound) THEN
!                            do j=SpatOrbs+1,SpatOrbs+NoRotAlphBet
!                                IF(SymOrbsTemp(j).lt.9) THEN
!                                    SymOrbs(k+1)=SymOrbsTemp(j)
!                                    CoeffT1(:,k+1)=NatOrbMat(:,j)
!                                    EvaluesTrunc(k+1)=Evalues(j)
!                                    SymOrbsTemp(j)=9
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
!                        SymFirst=SymOrbsTemp(i)
!                        CoeffT1(:,k)=NatOrbMat(:,i)
!                        EvaluesTrunc(k)=Evalues(i)
!                        SymOrbs(k)=SymOrbsTemp(i)
!                        do j=1,NoRotAlphBet
!                            IF(SymOrbsTemp(j).eq.SymFirst) THEN
!                                SymOrbs(k-1)=SymOrbsTemp(j)
!                                CoeffT1(:,k-1)=NatOrbMat(:,j)
!                                EvaluesTrunc(k-1)=Evalues(j)
!                                SymOrbsTemp(j)=9
!                                tSymFound=.true.
!                                EXIT
!                            ENDIF
!                        enddo

!                        IF(.not.tSymFound) THEN
!                            do j=1,NoRotAlphBet
!                                IF(SymOrbsTemp(j).lt.9) THEN
!                                    SymOrbs(k-1)=SymOrbsTemp(j)
!                                    CoeffT1(:,k-1)=NatOrbMat(:,j)
!                                    EvaluesTrunc(k-1)=Evalues(j)
!                                    SymOrbsTemp(j)=9
!                                    EXIT
!                                ENDIF
!                            enddo
!                        ENDIF
! 
!                        k=k+2
!                        IF(k.gt.(NoRotAlphBet*2)) EXIT
!                    enddo
!                ENDIF

                k=1
                do i=1,NoRotAlphBet
                    CoeffT1(:,k)=NatOrbMat(:,i)
                    EvaluesTrunc(k)=Evalues(i)
                    SymOrbs(k)=SymOrbsTemp(i)
                    k=k+2
                enddo
                k=2
                do i=SpatOrbs+1,SpatOrbs+NoRotAlphBet
                    CoeffT1(:,k)=NatOrbMat(:,i)
                    SymOrbs(k)=SymOrbsTemp(i)
                    EvaluesTrunc(k)=Evalues(i)
                    k=k+2
                enddo

            ELSE
                do i=1,NoRotAlphBet
                    CoeffT1(:,i)=NatOrbMat(:,i)
                    EvaluesTrunc(i)=Evalues(i)
                    SymOrbs(i)=SymOrbsTemp(i)
                enddo
            ENDIF

!            WRITE(6,*) SymOrbs(:)
!            CALL FLUSH(6)
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
!        CALL FLUSH(6)
!        stop

        IF(tTruncRODump) THEN

            Label=''
            LabelFull=''
            WRITE(Label,'(I5)') NoFrozenVirt
            LabelFull='EVALUES-TRUNC-'//adjustl(Label)

            OPEN(74,FILE=LabelFull,status='unknown')
            IF(tStoreSpinOrbs) THEN
                WRITE(74,*) NoOrbs-NoFrozenVirt
                do i=1,NoOrbs-NoFrozenVirt,2
                    WRITE(74,'(I5,ES20.10,I5,A5,I5,ES20.10,I5)') i,EvaluesTrunc(i),SymOrbs(i),'  *  ',i+1,EvaluesTrunc(i+1),SymOrbs(i+1)
                enddo
            ELSE
                WRITE(74,*) NoOrbs-NoFrozenVirt
                do i=1,NoOrbs-NoFrozenVirt
                    WRITE(74,'(ES20.10,I5)') EvaluesTrunc(i),SymOrbs(i)
                enddo
            ENDIF
            CLOSE(74) 
        ELSE
            OPEN(73,FILE='EVALUES',status='unknown')
            WRITE(73,*) NoOrbs
            IF(tStoreSpinOrbs) THEN
                k=0
                do i=1,NoOrbs,2
                    k=k+1
                    IF(tTruncRODump) THEN
                        WRITE(73,'(2I5,ES20.10,I5,A5,I5,ES20.10,I5)') (NoOrbs-i+1),i,Evalues(k),SymOrbs(i),'  *  ',i+1,Evalues(k+SpatOrbs),SymOrbs(i+1)
                    ELSE
                        WRITE(73,'(2I5,ES20.10,I5,A5,I5,ES20.10,I5)') (NoOrbs-i+1),i,Evalues(k),INT(G1(SymLabelList3(k))%Sym%S,4),'  *  ',&
                                                             &i+1,Evalues(k+SpatOrbs),INT(G1(SymLabelList3(k+SpatOrbs))%Sym%S,4)
                    ENDIF
                enddo
            ELSE
                do i=1,SpatOrbs
                    WRITE(73,'(3I5,ES20.10)') i,NoOrbs-i+1,(NoOrbs-i+1)*2,Evalues(i)
                enddo
            ENDIF
            CLOSE(73)
        ENDIF

        CALL HistNatOrbEvalues()

!        WRITE(6,*) 'NatOrbMat matrix'
!        do i=1,NoOrbs
!            WRITE(6,*) NatOrbMat(:,i)
!        enddo

!        OPEN(74,FILE='TRANSFORMMAT',status='unknown')
!        do i=1,NoOrbs
!            do j=1,NoOrbs-NoFrozenVirt
!                WRITE(74,*) i,j,CoeffT1(i,j)
!            enddo
!        enddo
!        CLOSE(74)

        CALL halt_timer(FillCoeff_Time)
 

    ENDSUBROUTINE FillCoeffT1




    SUBROUTINE HistNatOrbEvalues()
        IMPLICIT NONE
        INTEGER :: i,k,x,NoEvalues,a,b,NoOcc
        REAL*8 :: EvaluesCount(NoOrbs,2),OrbEnergies(1:NoOrbs),EvalueEnergies(1:NoOrbs)
        REAL*8 :: SumEvalues


        OPEN(74,FILE='EVALUES-PLOTRAT',status='unknown')
        IF(tStoreSpinOrbs) THEN
            k=0
            do i=1,SpatOrbs
                k=k+2
                WRITE(74,'(F20.10,ES20.10)') REAL(k-1)/REAL(NoOrbs),Evalues(i)
                WRITE(74,'(F20.10,ES20.10)') REAL(k)/REAL(NoOrbs),Evalues(SpatOrbs+i)
            enddo
        ELSEIF(tRotateOccOnly) THEN
            k=0
            do i=1,NoOcc
                k=k+1
                WRITE(74,'(F20.10,ES20.10)') REAL(k)/REAL(NoOcc),Evalues(i)
            enddo
        ELSEIF(tRotateVirtOnly) THEN
            k=NoOcc
            do i=NoOcc+1,NoOrbs
                k=k+1
                WRITE(74,'(F20.10,ES20.10)') REAL(k-NoOcc)/REAL(NoOrbs-NoOcc),Evalues(i)
            enddo
        ELSE
            k=0
            do i=1,SpatOrbs
                k=k+1
                WRITE(74,'(F20.10,ES20.10)') REAL(k)/REAL(NoOrbs),Evalues(i)
            enddo
        ENDIF
        CLOSE(74)

!        OPEN(73,FILE='EVALUES-plot',status='unknown')
!        EvaluesCount(:,:)=0.D0

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
!            EvaluesCount(k,2)=1.D0
!            do i=2,NoOrbs
!                IF((ABS(Evalues(i)-Evalues(i-1))).ge.(1E-10)) THEN
!                    k=k+1
!                    EvaluesCount(k,1)=Evalues(i)
!                    EvaluesCount(k,2)=1.D0
!                ELSE
!                    EvaluesCount(k,2)=EvaluesCount(k,2)+1.D0
!                ENDIF
!            enddo
!            NoEvalues=k

!            do i=1,NoEvalues
!                WRITE(73,*) EvaluesCount(i,1),Evaluescount(i,2)
!            enddo


!        enddo

!        CLOSE(73)

! Want to write out the eigenvectors in order of the energy of the new orbitals - so that we can see the occupations 
! of the type of orbital.
! For now, keep this separate to the transformation of ARR - even though it is equivalent.

!        WRITE(6,*) 'ARR'
!        do i=1,NoOrbs
!            WRITE(6,*) ARR(2*i,2)
!        enddo
!        WRITE(6,*) 'Evalues'
!        do i=1,NoOrbs
!            WRITE(6,*) Evalues(i)
!        enddo

        OrbEnergies(:)=0.D0
        EvalueEnergies(:)=0.D0
        SumEvalues=0.D0
        do i=1,NoOrbs
            IF(tStoreSpinOrbs) THEN
                SumEvalues=SumEvalues+Evalues(i)
            ELSE
                SumEvalues=SumEvalues+(2*Evalues(i))
            ENDIF
            EvalueEnergies(i)=Evalues(i)
! We are only interested in the diagonal elements.            
            do a=1,NoOrbs
                b=SymLabelList2(a)
                IF(tStoreSpinOrbs) THEN
                    OrbEnergies(i)=OrbEnergies(i)+(NatOrbMat(a,i)*ARR(b,2)*NatOrbMat(a,i))
                ELSE
                    OrbEnergies(i)=OrbEnergies(i)+(NatOrbMat(a,i)*ARR(2*b,2)*NatOrbMat(a,i))
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
        CALL Sort2Real(NoOrbs,OrbEnergies(1:NoOrbs),EvalueEnergies(1:NoOrbs))

        OPEN(73,FILE='EVALUES-ENERGY',status='unknown')
        do i=1,NoOrbs
            WRITE(73,*) OrbEnergies(NoOrbs-i+1),EvalueEnergies(NoOrbs-i+1)
        enddo
        WRITE(73,*) 'The sum of the occupation numbers (eigenvalues) = ',SumEvalues
        WRITE(73,*) 'The number of electrons = ',NEl
        CALL FLUSH(73)
        CLOSE(73)
        CALL FLUSH(6)


        OPEN(73,FILE='OccupationTable',status='unknown')
        x=1
        do while (x.le.NoOrbs)
            WRITE(73,'(A16,A5)',advance='no') 'HF Orb En    ','Sym'
            do i=x,x+9
                IF(i.gt.NoOrbs) THEN
                    WRITE(73,*) ''
                    EXIT
                ENDIF
                WRITE(73,'(ES16.6)',advance='no') Evalues(i)
            enddo
            WRITE(73,*) ''

            do a=1,NoOrbs
                b=SymLabelListInv(a)
                IF(tStoreSpinOrbs) THEN
                    WRITE(73,'(F16.10,I5)',advance='no') ARR(a,1),INT(G1(a)%sym%S,4)
                ELSE
                    WRITE(73,'(F16.10,I5)',advance='no') ARR(2*a,1),INT(G1(2*a)%sym%S,4)
                ENDIF
                do i=x,x+9
                    IF(i.gt.NoOrbs) THEN
                        WRITE(73,*) ''
                        EXIT
                    ENDIF
                    WRITE(73,'(F16.10)',advance='no') NatOrbMat(b,i)
                enddo
                WRITE(73,*) ''
            enddo
            WRITE(73,*) ''
            x=x+10
        enddo
        CALL FLUSH(73)
        CLOSE(73)
        CALL FLUSH(6)
        

    END SUBROUTINE HistNatOrbEvalues



    SUBROUTINE DeallocateNatOrbs()
        USE Logging , only : tTruncRODump
        IMPLICIT NONE
        CHARACTER(len=*), PARAMETER :: this_routine='DeallocateNatOrbs'

! Deallocate the natural orbitals matrix.    

        IF(tTruncRODump) THEN
            DEALLOCATE(SymOrbsTemp)
            CALL LogMemDeAlloc(this_routine,SymOrbsTempTag)
        ENDIF
        DEALLOCATE(NatOrbMat)
        CALL LogMemDeAlloc(this_routine,NatOrbMatTag)
        DEALLOCATE(Evalues)
        CALL LogMemDeAlloc(this_routine,EvaluesTag)

    END SUBROUTINE DeallocateNatOrbs




!This file was primarily concerned with the creation of natural orbitals from a rotation of the previous orbitals.
!The 1-electron Reduced density matrix was inputted, and the natural orbitals constructed. From there, the
!1 and 2 electron integrals were transformed and replaced into UMat.
    SUBROUTINE FindNatOrbsOld()
        IMPLICIT NONE
        INTEGER :: i,j

        OPEN(12,FILE='ONEEL-RDM',STATUS='UNKNOWN')
        do i=1,nBasis
            do j=1,nBasis
                WRITE(12,"(F18.7)",advance='no') NatOrbMat(i,j)
            enddo
            WRITE(12,*) ""
        enddo
        CLOSE(12)

        CALL Stop_All('FindNatOrbsOld','This is the old routine for finding the natural orbitals - likely buggy.')

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
! The diagonalisation routine reorders the orbitals in such a way that the corresponding orbital labels are lost.
! In order to keep the spin and spatial symmetries, each symmetry must be fed into the diagonalisation routine separately.
! The best way to do this is to order the orbitals so that all the alpha orbitals follow all the beta orbitals, with the 
! occupied orbitals first, in terms of symmetry, and the virtual second, also ordered by symmetry.
! This gives us flexibility w.r.t rotating only the occupied or only virtual and looking at high spin states.
        IMPLICIT NONE
        REAL*8 , ALLOCATABLE :: NOccNums(:),Work(:)
        INTEGER :: nOccNumsTag=0,iErr,WorkSize,WorkCheck,WorkTag=0,i
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
        CALL FLUSH(6)

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
        CALL FLUSH(6)

    END SUBROUTINE Diag1RDMOld

END MODULE NatOrbsMod
