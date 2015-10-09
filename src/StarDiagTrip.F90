    MODULE StarDiagTripMod
      use constants, only: dp
!      USE global_utilities
!      use SystemData , only : NEl
!      USE Determinants , only : FDet


!      IMPLICIT NONE

!Try to maintain same name system as StarDiagMod, therefore
!.. ExcitInfo(0,...) corresponds to J=I
!.. ExcitInfo(J,0) = RHOJJ
!.. ExcitInfo(J,1) = RHOIJ
!.. ExcitInfo(J,2) = HIJ
!      HElement_t(dp) , ALLOCATABLE :: ExcitInfo(:,:)
!      INTEGER :: ExcitInfoTag=0
!
!!Triples stores info on the number of triples from each double
!      INTEGER , ALLOCATABLE :: Triples(:)
!      INTEGER :: TriplesTag=0
!
!!TripsInfo stores the star information for the excited stars the same way as ExcitInfo
!      HElement_t(dp) , ALLOCATABLE :: TripsInfo(:,:)
!      INTEGER :: TripsInfoTag=0
!
!!Incase we want to do a complete diagonalisation of the triples star
!      real(dp) , ALLOCATABLE :: ExcitMat(:,:)
!      real(dp) , ALLOCATABLE :: HamMat(:)
!      INTEGER :: ExcitMatTag=0,HamMatTag=0
!
!!Vals and first element of the eigenvectors from the prediagonalised excited stars
!      real(dp) , ALLOCATABLE :: Vals(:)
!      real(dp) , ALLOCATABLE :: Vecs(:)
!      INTEGER :: ValsTag=0,VecsTag=0
!      
!      HElement_t(dp) :: rhii
!      HElement_t(dp) :: Hii

      Contains

! This routine creates stars of all triples from all doubles. These stars of triples are diagonalised, and the
! eigenvector components reattached to the HF determinant. Should scale as N^3 M^3.
! Currently, only PolyMax works with the routine (lowest energy state), and there is a forced counting of matrix elements.
    SUBROUTINE StarDiagTrips(Energyxw,Weight)
!        use SystemData, only: Alat,Beta,Brr,ECore,G1,nBasis,nBasisMax,nMsh,Arr
!        use CalcData , only : i_P,NWHTay,RhoEps,dBeta,TFullDiag
!        use IntegralsData, only : fck,nMax,UMat,nTay
!        USE Determinants , only : GetHElement2
!        USE LoggingData , only : iLogging
!        USE StarDiagMod , only : GetValsnVecs
!        IMPLICIT NONE
!        CHARACTER(len=*), PARAMETER :: this_routine='StarDiagTrips'
!        real(dp) , ALLOCATABLE :: TempDiags(:),Work(:),Vals2(:)
!        INTEGER :: TempDiagsTag=0,WorkTag=0,Vals2Tag=0
!        real(dp) :: Energy
        real(dp) :: Weight,Energyxw
!        HElement_t(dp) :: rh,rhjj,rhij,Norm,OffDiagNorm,Hij,rhjk
!        INTEGER , ALLOCATABLE :: nExcit(:),nExcit2(:)
!        INTEGER :: nExcitTag,ierr,exFlagHF,exFlagDoub,iMaxExcit,ExcitCount
!        INTEGER :: Meth,nStore(6),nExcitMemLen,nJ(NEl),iExcit,nStore2(6)
!        type(timer), save :: proc_timerTrips
!        INTEGER :: iMaxExcit2,nExcitMemLen2,nK(NEl),ExcitCountDoubs,nExcitTag2
!        INTEGER :: i,j,k,nRoots,Vert,IC,iGetExcitLevel,DoubIndex,TotElem,Info
!        LOGICAL :: TCountExcits,iExcit2
    
        ! Disable warnings
        weight = weight
        energyxw = energyxw

        CALL Stop_All("StarDiagTrips","This code has now been commented out.")

!        IF(HElement_t_size.gt.1) STOP 'Only Real orbitals allowed in StarDiagTrips so far'
!    
!        proc_timerTrips%timer_name='StarDiagTrips'
!        call set_timer(proc_timerTrips)
!        nExcitTag=0
!        nExcitTag2=0
!
!        Meth=NWHTay(2,2)
!        TCountExcits=BTEST(Meth,8)
!!Allow single and double excitations from HF
!        exFlagHF=3
!!Allow only single excitations from doubles (giving singles and triples)
!        exFlagDoub=1
!
!        IF(.NOT.TCountExcits) THEN
!            WRITE(6,*) "Currently, you have to count the excitations in StarDiag Trips"
!            STOP "Currently, you have to count the excitations in StarDiag Trips"
!        ENDIF
!
!!Find rhii
!        CALL CalcRho2(FDet,FDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhii,nTay,0,ECore)
!
!!Counting excitations
!        IF(TCountExcits) THEN
!            nStore(1)=0
!            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlagHF)
!            ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
!            CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag)
!            nExcit(1)=0
!            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlagHF)
!
!!We have to provide memory to store the number of triples from each double
!            ALLOCATE(Triples(iMaxExcit),stat=ierr)
!            CALL LogMemAlloc('Triples',iMaxExcit,4,this_routine,TriplesTag)
!            Triples(1:iMaxExcit)=0
!
!!Run through all possible excitations from HF...
!            ExcitCount=0
!            ExcitCountDoubs=0
!       lp2: do while(.true.)
!                CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,exFlagHF)
!                IF(nJ(1).eq.0) exit lp2
!                CALL CalcRho2(FDet,nJ,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,iExcit,ECore)
!                IF(abs(rh).gt.RhoEps) THEN
!!Count contribution from double excitation...
!                    ExcitCount=ExcitCount+1
!                    ExcitCountDoubs=ExcitCountDoubs+1
!
!!Now, setup excitation generators from each double excitation
!                    nStore2(1)=0
!                    CALL GenSymExcitIt2(nJ,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen2,nK,iMaxExcit2,0,nStore2,exFlagDoub)
!                    ALLOCATE(nExcit2(nExcitMemLen2),stat=ierr)
!                    CALL LogMemAlloc('nExcit2',nExcitMemLen2,4,this_routine,nExcitTag2)
!                    nExcit2(1)=0
!                    CALL GenSymExcitIt2(nJ,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit2,nK,iMaxExcit,0,nStore2,exFlagDoub)
!                    
!                    i=0
!!Run through all single excitations from doubles
!               lp1: do while(.true.)
!                        CALL GenSymExcitIt2(nJ,NEl,G1,nBasis,nBasisMax,.false.,nExcit2,nK,iExcit2,0,nStore2,exFlagDoub)
!                        IF(nK(1).eq.0) exit lp1
!!                        CALL CalcRho2(nJ,nK,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,iExcit2,ECore)
!                        rh=(0.0_dp)
!!Uncomment this code if you want to only allow triple excitations from each double
!!                        IC=iGetExcitLevel(FDet,nK,NEl)
!!                        IF(IC.ne.3) THEN
!!                            rh=(0.0_dp)
!!                        ENDIF
!                        IF(abs(rh).gt.RhoEps) THEN
!                            i=i+1
!                            ExcitCount=ExcitCount+1
!                        ENDIF
!                    enddo lp1
!
!!Store the number of triples from each double
!                    Triples(ExcitCountDoubs)=i
!
!!Destroy doubles excitation generators
!                    DEALLOCATE(nExcit2)
!                    CALL LogMemDealloc(this_routine,nExcitTag2)
!
!                ENDIF
!        
!            enddo lp2
!
!            WRITE(6,*) "Double Excitations found: ",ExcitCountDoubs
!            WRITE(6,*) "Including ALL single excitations from each double as a separate star..."
!            WRITE(6,*) "Total Excitations found: ",ExcitCount
!        
!!Destroy excitation generator from HF
!            DEALLOCATE(nExcit)
!            CALL LogMemDealloc(this_routine,nExcitTag)
!    
!        ENDIF
!
!!Allocate memory for calculation
!        IF(TFullDiag) THEN
!            ALLOCATE(ExcitMat(ExcitCount+1,ExcitCount+1),stat=ierr)
!            CALL LogMemAlloc('ExcitMat',(ExcitCount+1)**2,8,this_routine,ExcitMatTag)
!            ExcitMat=0.0_dp
!
!            ExcitMat(1,1)=1.0_dp
!
!            ALLOCATE(HamMat(ExcitCount+1),stat=ierr)
!            CALL LogMemAlloc('HamMat',ExcitCount+1,8,this_routine,HamMatTag)
!            HamMat=0.0_dp
!            Hii=GetHElement2(FDet,FDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)
!            HamMat(1)=Hii
!            
!        ENDIF
!
!        ALLOCATE(ExcitInfo(0:ExcitCount,0:2),stat=ierr)
!        CALL LogMemAlloc('ExcitInfo',(ExcitCount+1)*3,8*HElement_t_size,this_routine,ExcitInfoTag)
!        ExcitInfo=(0.0_dp)
!
!!Still divide all elements by rhii
!        ExcitInfo(0,0)=(1.0_dp)
!        ExcitInfo(0,1)=(1.0_dp)
!        ExcitInfo(0,2)=GetHElement2(FDet,FDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)
!        
!!Reinitialise HF excitation generator
!        nStore(1)=0
!        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlagHF)
!        ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
!        CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag)
!        nExcit(1)=0
!        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlagHF)
!
!!Run through all possible excitations from HF. j counts the double excitations, while Vert counts the 
!vertex which is being reattached to HF
!        j=0
!        Vert=0
!        TotElem=1
!   lp3: do while(.true.)
!            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,exFlagHF)
!            IF(nJ(1).eq.0) exit lp3
!            CALL CalcRho2(FDet,nJ,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhij,nTay,iExcit,ECore)
!            IF(abs(rhij).gt.RhoEps) THEN
!                j=j+1
!
!!Find matrix elements for double excitation
!                CALL CalcRho2(nJ,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhjj,nTay,0,ECore)
!                Hij=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,iExcit,ECore)
!                IF(TFullDiag) THEN
!                    TotElem=TotElem+1
!                    ExcitMat(1,TotElem)=(rhij)/(rhii)
!                    ExcitMat(TotElem,1)=(rhij)/(rhii)
!                    ExcitMat(TotElem,TotElem)=(rhjj)/(rhii)
!                    HamMat(TotElem)=(Hij)
!                    DoubIndex=TotElem
!                ENDIF
!
!!Allocate memory for triples star
!                ALLOCATE(TripsInfo(0:Triples(j),0:1),stat=ierr)
!                CALL LogMemAlloc('TripsInfo',(Triples(j)+1)*2,8*HElement_t_size,this_routine,TripsInfoTag)
!                TripsInfo=(0.0_dp)
!
!!Need to divide everything by the rhjj element. The final eigenvalues will need to be multiplied by them at the end
!                TripsInfo(0,0)=(1.0_dp)
!                TripsInfo(0,1)=(1.0_dp)
!
!!Now, setup excitation generators from each double excitation
!                nStore2(1)=0
!                CALL GenSymExcitIt2(nJ,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen2,nK,iMaxExcit2,0,nStore2,exFlagDoub)
!                ALLOCATE(nExcit2(nExcitMemLen2),stat=ierr)
!                CALL LogMemAlloc('nExcit2',nExcitMemLen2,4,this_routine,nExcitTag2)
!                nExcit2(1)=0
!                CALL GenSymExcitIt2(nJ,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit2,nK,iMaxExcit2,0,nStore2,exFlagDoub)
!
!                IF(iMaxExcit2.lt.Triples(j)) THEN
!                    WRITE(6,*) "Problem in counting here..."
!                    STOP 'Problem in counting here...'
!                ENDIF
!                
!                i=0
!!Run through all single excitations from doubles, which i counts
!           lp4: do while(.true.)
!                    CALL GenSymExcitIt2(nJ,NEl,G1,nBasis,nBasisMax,.false.,nExcit2,nK,iExcit2,0,nStore2,exFlagDoub)
!                    IF(nK(1).eq.0) exit lp4
!                    rhjk=(0.0_dp)
!!                    CALL CalcRho2(nJ,nK,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhjk,nTay,iExcit2,ECore)
!!Uncomment this code if you want to only allow triple excitations from each double
!!                    IC=iGetExcitLevel(FDet,nK,NEl)
!!                    IF(IC.ne.3) THEN
!!                        rh=(0.0_dp)
!!                    ENDIF
!                    IF(abs(rhjk).gt.RhoEps) THEN
!                        i=i+1
!                        TripsInfo(i,1)=rhjk/rhjj
!                        
!!Calculate diagonal element
!                        CALL CalcRho2(nK,nK,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,0,ECore)
!                        TripsInfo(i,0)=rh/rhjj
!
!                        IF(TFullDiag) THEN
!                            TotElem=TotElem+1
!                            ExcitMat(DoubIndex,TotElem)=(rhjk)/(rhii)
!                            ExcitMat(TotElem,DoubIndex)=(rhjk)/(rhii)
!                            ExcitMat(TotElem,TotElem)=(rh)/(rhii)
!                        ENDIF
!
!                    ENDIF
!
!                enddo lp4
!
!                IF(i.ne.Triples(j)) THEN
!                    WRITE(6,*) "Problem in counting here..."
!                    STOP 'Problem in counting here...'
!                ENDIF
!
!!Deallocate the excitation generators for the triples
!                DEALLOCATE(nExcit2)
!                CALL LogMemDealloc(this_routine,nExcitTag2)
!                nStore2(1:6)=0
!
!!Allocate memory to store all eigenvalues, and first elements of eigenvectors of matrix
!                ALLOCATE(Vals(i+1),stat=ierr)
!                CALL LogMemAlloc('Vals',i+1,8,this_routine,ValsTag)
!                ALLOCATE(Vecs(i+1),stat=ierr)
!                CALL LogMemAlloc('Vecs',i+1,8,this_routine,VecsTag)
!
!!Next, we need to diagonalise this excited star matrix. ALL eigenvalues, and the first 
!element of all the eigenvectors needed...bummer...
!                IF(.NOT.BTEST(Meth,0)) THEN
!!This will diagonalise each excited star fully - v. slow - order N^3
!!For some annoying reason, DiagRhos needs to be of type real(dp)...
!                    
!                    ALLOCATE(TempDiags(i+1),stat=ierr)
!                    CALL LogMemAlloc('TempDiags',i+1,8,this_routine,TempDiagsTag)
!                    do k=1,i+1
!                        TempDiags(k)=TripsInfo(k-1,0)
!                    enddo
!
!                    CALL GetValsnVecs(i+1,TempDiags(1:i+1),TripsInfo(1:i,1),Vals,Vecs)
!
!                    DEALLOCATE(TempDiags)
!                    CALL LogMemDealloc(this_routine,TempDiagsTag)
!
!                ELSE
!!Diagonalise excited stars using polynomial root searching algorithm - this is order N^2
!                    STOP 'Polynomial diagonalisation is not implemented yet - use diag'
!                    
!                ENDIF
!
!!Excitation information about excited stars can now be discarded...
!                DEALLOCATE(TripsInfo)
!                CALL LogMemDealloc(this_routine,TripsInfoTag)
!
!!Since the excited star matrix elements were divided by rhjj, the eigenvalues need to be multiplied by this value.
!!In the final star matrix, all the elements will also be divided by rhii, so this needs to be done for all the elements.
!                Norm=rhjj/rhii
!                OffDiagNorm=rhij/rhii
!
!!The eigenvectors from the excited stars now need to be reattached to the original root...
!                do k=1,i+1
!                    Vert=Vert+1
!                    ExcitInfo(Vert,0)=(Vals(k))*Norm
!                    ExcitInfo(Vert,1)=(Vecs(k))*OffDiagNorm
!                    ExcitInfo(Vert,2)=(Vecs(k))*Hij
!                enddo
!
!                DEALLOCATE(Vecs)
!                CALL LogMemDealloc(this_routine,VecsTag)
!                DEALLOCATE(Vals)
!                CALL LogMemDealloc(this_routine,ValsTag)
!
!            ENDIF
!
!            DEALLOCATE(Vecs)
!            CALL LogMemDealloc(this_routine,VecsTag)
!            DEALLOCATE(Vals)
!            CALL LogMemDealloc(this_routine,ValsTag)
!
!        enddo lp3
!
!!Deallocate the excitation generators for the doubles 
!        DEALLOCATE(nExcit)
!        CALL LogMemDealloc(this_routine,nExcitTag)
!        DEALLOCATE(Triples)
!        CALL LogMemDealloc(this_routine,TriplesTag)
!        nStore(1:6)=0
!        
!        IF(TFullDiag.and.(TotElem.ne.(ExcitCount+1))) THEN
!            WRITE(6,*) "We have a counting problem"
!            STOP "We have a counting problem"
!        ENDIF
!        IF(j.ne.ExcitCountDoubs) THEN
!            WRITE(6,*) "We have a counting problem"
!            STOP "We have a counting problem"
!        ENDIF
!
!        IF(Vert.ne.ExcitCount) THEN
!            WRITE(6,*) "We have a counting problem"
!            STOP "We have a counting problem"
!        ENDIF
!
!!We now have a large star matrix, with its information in ExcitInfo, which we need to diagonalise.
!!We need to diagonalise this large final star matrix. Only largest eigenvector needed
!            WRITE(65,*) "RHOMAT is "
!            DO I=1,ExcitCount+1
!                DO J=1,ExcitCount+1
!                    IF(I.eq.1) THEN
!                        WRITE(65,"(E14.6)",advance='no') ExcitInfo(J-1,1)
!                    ELSEIF(J.eq.1) THEN
!                        WRITE(65,"(E14.6)",advance='no') ExcitInfo(I-1,1)
!                    ELSEIF(J.eq.I) THEN
!                        WRITE(65,"(E14.6)",advance='no') ExcitInfo(I-1,0)
!                    ELSE
!                        WRITE(65,"(E14.6)",advance='no') 0.0_dp
!                    ENDIF
!                ENDDO
!                WRITE(65,*) ""
!            ENDDO
!            WRITE(65,*) "**************"
!            DO I=1,ExcitCount+1
!                WRITE(65,"(E14.6)") ExcitInfo(I-1,2)
!            ENDDO
!!        IF(.NOT.BTEST(Meth,0)) THEN
!!This will diagonalise each excited star fully - v. slow - order N^3
!!            WRITE(6,*) "Beginning Complete Star Tridiagonalization"
!            CALL StarDiag(ExcitCount+1,ExcitInfo,ExcitCount+1,i_P,Weight,dBeta,Energyxw)
!
!!        ELSE
!!Use polynomial diagonalisation, order N
!!            WRITE(6,*) "Beginning Polynomial Star Diagonalization"
!!            nRoots=ExcitCount
!!            IF(BTEST(Meth,1)) THEN
!!                nRoots=1
!!                WRITE(6,*) "Searching for 1 root"
!!            ELSEIF(BTEST(Meth,2)) THEN
!!                WRITE(6,*) "Searching for enough roots to converge sum"
!!                nRoots=ExcitCount+1
!!            ELSEIF(BTEST(Meth,6)) THEN
!!                WRITE(6,*) "Searching for enough roots to converge sum - Method 2"
!!                nRoots=ExcitCount+2
!!            ELSE
!!                WRITE(6,*) "Searching for all roots"
!!            ENDIF
!
!
!!            CALL SORT3RN(ExcitCount,ExcitInfo(1:ExcitCount,0),ExcitInfo(1:ExcitCount,1),ExcitInfo(1:ExcitCount,2),HElement_t_size)
!
!
!!            CALL StarDiag2(ExcitCount+1,ExcitInfo,ExcitCount+1,i_P,Weight,dBeta,Energyxw,nRoots)
!!
!!        ENDIF
!
!        IF(TFullDiag) THEN
!!Need to diagonalise the full matrix - this could take a while!
!            ALLOCATE(Work(3*(ExcitCount+1)),stat=ierr)
!            CALL LogMemAlloc('Work',3*(ExcitCount+1),8,this_routine,WorkTag)
!            ALLOCATE(Vals2(ExcitCount+1),stat=ierr)
!            CALL LogMemAlloc('Vals2',ExcitCount+1,8,this_routine,Vals2Tag)
!            
!            WRITE(66,*) "RHOMAT is "
!            DO I=1,ExcitCount+1
!                DO J=1,ExcitCount+1
!                    WRITE(66,"(E14.6)",advance='no') ExcitMat(I,J)
!                ENDDO
!                WRITE(66,*) ""
!            ENDDO
!            WRITE(66,*) "**************"
!            DO I=1,ExcitCount+1
!                WRITE(66,"(E14.6)") HamMat(I)
!            ENDDO
!
!
!            
!            CALL DSYEV('V','U',ExcitCount+1,ExcitMat,ExcitCount+1,Vals2,Work,3*(ExcitCount+1),INFO)
!            IF(INFO.ne.0) THEN
!                WRITE(6,*) "DYSEV error in AddTrips: ",INFO
!                STOP
!            ENDIF
!            WRITE(66,*) "**************"
!            DO I=1,ExcitCount+1
!                WRITE(66,*) Vals2(i) 
!            ENDDO
!            WRITE(66,*) "**************"
!            DO I=1,ExcitCount+1
!                WRITE(66,*) ABS(ExcitMat(I,ExcitCount+1))
!            ENDDO
!                
!            
!            Energy=0.0_dp
!            do i=2,ExcitCount+1
!                Energy=Energy+HamMat(i)*ABS(ExcitMat(i,ExcitCount+1))
!            enddo
!            Energy=Energy/(ABS(ExcitMat(1,ExcitCount+1)))
!            Energy=Energy+(Hii)
!
!            WRITE(6,"(A,G25.17)") "From complete diagonalisation of the rho matrix, the energy is given as: ",Energy
!
!            DEALLOCATE(HamMat)
!            CALL LogMemDealloc(this_routine,HamMatTag)
!            DEALLOCATE(ExcitMat)
!            CALL LogMemDealloc(this_routine,ExcitMatTag)
!            DEALLOCATE(Work)
!            CALL LogMemDealloc(this_routine,WorkTag)
!            DEALLOCATE(Vals2)
!            CALL LogMemDealloc(this_routine,Vals2Tag)
!
!        ENDIF
!
!        DEALLOCATE(ExcitInfo)
!        CALL LogMemDealloc(this_routine,ExcitInfoTag)
!
!        IF(iErr.ne.0) THEN
!            WRITE(6,*) "Error in allocation/deallocation"
!            STOP "Error in allocation/deallocation"
!        ENDIF
!
!        call halt_timer(proc_timerTrips)
!        
!        RETURN

    END SUBROUTINE StarDiagTrips

    END MODULE StarDiagTripMod
