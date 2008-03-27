!This code it designed to generate a random graph, with Ndets distinct determinants in it (including HF), and improve it in an iterative fashion.
!Once a graph is found, all possible excitations from each determinant  in the graph are found, and then the rho matrix from these additional
!excitations operate on the original eigenvector of the graph (null padded). This produces a much larger vector (not eigenvector), from which
!excitations which have a strong connection can be chosen stochastically. This process is then repeated with the new graph, 
!chosen from the vector (ensuring connectivity). This will improve the connections between the determinants, 
!which will hopefully lower the energy after only a few iterations.

MODULE GraphMorph
    
    USE System , only : NEl
    USE Determinants , only : FDet
!Iters is the number of interations of morphing the graph. Nd is the number of determinants in the graph.
    USE Calc , only : Iters,NDets,GraphBias,Biasing
    USE MemoryManager , only : LogMemAlloc,LogMemDealloc
    USE HElem
    IMPLICIT NONE
    SAVE
!This will hold all the determinants in the graph
    INTEGER , ALLOCATABLE :: GraphDets(:,:)
    INTEGER :: GraphDetsTag=0

!This will hold all the Hamiltononian elements in the graph (zero if > double excit) (H_11) = HF
    TYPE(HElement) , ALLOCATABLE :: HamElems(:)
    INTEGER :: HamElemsTag=0
    
!This will hold all the determinants of the excitations of the graph
    INTEGER , ALLOCATABLE :: ExcitsDets(:,:)
    INTEGER :: ExcitsDetsTag=0

!This is the number of possible excitations from each determinant in the graph, stored
!in a cumulative way
    INTEGER , ALLOCATABLE :: NoExcits(:)
    INTEGER :: NoExcitsTag=0
!The total number of excitations from the graph
    INTEGER :: TotExcits

!Holds the rho elements between each determinant and its excitations
    TYPE(HElement) , ALLOCATABLE :: ConnectionsToExcits(:)
    INTEGER :: ConnectionsToExcitsTag=0

!The full Rho Matrix for the graph
    TYPE(HElement) , ALLOCATABLE :: GraphRhoMat(:,:)
    INTEGER :: GraphRhoMatTag=0
!The Largest Eigenvector for the graph
    TYPE(HElement) , ALLOCATABLE :: Eigenvector(:)
    INTEGER :: EigenvectorTag=0
!The largest Eigenvalue of the graph
    REAL*8 :: Eigenvalue

!This is the vector of propensity to move towards the excited determinants of the graph
    TYPE(HElement) , ALLOCATABLE :: ExcitsVector(:)
    INTEGER :: ExcitsVectorTag=0

!The rho matrix and Hamiltonian element for the HF determinant
    TYPE(HElement) :: rhii,Hii

!These are the weight and energy of the current graph respectively
    REAL*8 :: SI,DLWDB

!This is the seed for the random numbers needed in the routines
    INTEGER :: Seed

!Various stats for printing
    REAL*8 :: PStay,Orig_Graph,SucRat,MeanExcit

    contains

    SUBROUTINE MorphGraph(Weight,Energyxw)
        USE System, only: Alat,Beta,Brr,ECore,G1,nBasis,nBasisMax,Arr
        USE Calc , only : i_P,G_VMC_Seed
        USE Integrals, only : fck,nMax,nMsh,UMat,nTay
        USE Determinants , only : GetHElement2
        IMPLICIT NONE
        TYPE(HDElement) :: Weight,Energyxw
        INTEGER :: ierr,i
        CHARACTER(len=*), PARAMETER :: this_routine='MorphGraph'

        OPEN(63,file='MorphStats',Status='unknown')
        
        IF(HElementSize.ne.1) STOP 'Only real orbitals allowed in GraphMorph so far'

!Initialise random number generator
        Seed=G_VMC_Seed

!Find rho_ii value, which all rho elements will be divided by, and Hii value
        CALL CalcRho2(FDet,FDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhii,nTay,0,ECore)
        Hii=GetHElement2(FDet,FDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)

!It is first necessary to construct the original graph to work with.
        WRITE(6,"(A,I5,A)") "Constructing random initial graph to morph from with ", NDets, " vertices."
        CALL ConstructInitialGraph()
        
!Allocate space for largest Eigenvector for various graphs
        ALLOCATE(Eigenvector(NDets),stat=ierr)
        CALL LogMemAlloc('Eigenvector',NDets,8*HElementSize,this_routine,EigenvectorTag)

        WRITE(63,*) "Iteration  Energy  P_stay  Orig_Graph  SucRat   MeanExcit"
        IF(Biasing) THEN
            WRITE(6,"(A,F10.7)") "Graph growing biased towards original determinants with probability ", GraphBias
        ENDIF

!Once the graph is found, the loop over the morphing iterations can begin
        do i=1,Iters

            WRITE(6,*) ""
            WRITE(6,"(A,I4,A)") "Starting Iteration ", i, " ..."
            WRITE(6,*) ""
            CALL FLUSH(6)

!The graph is first diagonalised, and the energy of the graph found, along with the largest eigenvector and eigenvalues.
            CALL DiagGraphMorph()
            WRITE(6,"(A,2G20.12)") "Weight and Energy of current graph is: ", SI, DLWDB
            CALL FLUSH(6)

!Excitation generators are initialised for each of the determinants in the graph, and the total number of possible
!connected determinants calculated. Memory allocated for the ensuing calculation.
            CALL CountExcits()
!            WRITE(6,*) "Fraction of space which is space of excitations is: ", TotExcits/(TotExcits+NDets)
            WRITE(6,*) "Total number of determinants available from current graph is: ",TotExcits

!Run through each determinant in the graph, calculating the excitations, storing them, and the rho elements to them.
            CALL FindConnections()

!Multiply the correct rho element, with the correct element of the eigenvector, to create the value for the 
!tendancy to move to that excited determinant - store this in ExcitsVector
            CALL CreateExcitsVector()

!Add the original eigenvector*eigenvalue to the list of determinants - now have vector (contained in Eigenvector*eigenvalue and ExcitsVector) with all
!determinants in graph, and all possible determinants to excite to. This needs normalising.
            CALL NormaliseVector()

!Pick NDets new excitations stocastically from normalised list of determinants with probability |c|^2. Ensure connections,
!allocate and create rho matrix for new graph. Deallocate info for old graph.
            WRITE(6,*) "Choosing new graph stochastically from previous graph and its excitations..."
            CALL FLUSH(6)
            CALL PickNewDets()

!Write out stats
            WRITE(63,"(I10,5G20.12)") i,DLWDB,PStay,Orig_Graph,SucRat,MeanExcit
            CALL FLUSH(63)

!Once graph is fully constructed, the next iteration can begin.
        enddo
        WRITE(6,*) ""

!Diagonalise final graph, and find weight & energy. 
!There is no beta-dependance, so only largest eigenvector and value needed.
        CALL DiagGraphMorph()
        WRITE(6,"(A,2G20.12)") "Weight and Energy of final graph is: ", SI, DLWDB

!Deallocate info...
        DEALLOCATE(Eigenvector)
        CALL LogMemDealloc(this_routine,EigenvectorTag)
        DEALLOCATE(GraphDets)
        CALL LogMemDealloc(this_routine,GraphDetsTag)
        DEALLOCATE(HamElems)
        CALL LogMemDealloc(this_routine,HamElemsTag)

        Weight=HDElement(SI-1.D0)
        Energyxw=HDElement((DLWDB*SI)-Hii%v)
!        Weight=SI-1.D0
!        Energy=DLWDB-(Hii%v)

        Close(63)
        
        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in ConstructInitialGraph'

    END SUBROUTINE MorphGraph

    SUBROUTINE ConstructInitialGraph()
        USE System , only : G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,Arr
        USE Calc , only : i_P,RhoEps,G_VMC_Pi
        USE Integrals , only : fck,nMax,nMsh,UMat,nTay
        USE Determinants , only : GetHElement2
        IMPLICIT NONE
        INTEGER :: ierr,nStore(6),nExcitMemLen,iMaxExcit,nJ(NEl),nExcitTag,iExcit
        INTEGER :: iPathTag,XijTag,RhoiiTag,RhoijTag,HijsTag,iSubInit,i,j,diff
        INTEGER , ALLOCATABLE :: iPath(:,:)
        REAL*8 , ALLOCATABLE :: Xij(:,:)
#if defined(POINTER8)
        INTEGER*8 :: ExcitGen(0:NDets)
#else
        INTEGER :: ExcitGen(0:NDets)
#endif
        TYPE(HDElement) , ALLOCATABLE :: Rhoii(:)
        TYPE(HElement) , ALLOCATABLE :: Rhoij(:,:),Hijs(:)
        TYPE(HElement) :: rh
        INTEGER , ALLOCATABLE :: nExcit(:)
        REAL*8 :: PGen,OldImport
        CHARACTER(len=*), PARAMETER :: this_routine='ConstructInitialGraph'
        
        CALL TISET('ConsInitGraph',iSubInit)

!        WRITE(6,*) "FDET is ",FDet(:)

!Set the importance parameter to be equal to 1 if we want random double excitation connected star graphs as initial graphs.
        OldImport=G_VMC_Pi
        G_VMC_Pi=1.D0
        
!Set Tags to zero, so we know when they are allocated/deallocated
        nExcitTag=0
        iPathTag=0
        XijTag=0
        RhoiiTag=0
        RhoijTag=0
        HijsTag=0
!Setup excitation generator
        CALL IAZZERO(nStore,6)
        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
        ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
        CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag)
        CALL IAZZERO(nExcit,nExcitMemLen)
        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,3)
        CALL GenRandSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,nExcit,nJ,Seed,iExcit,0,UMat,nMax,PGen)

!Allocate memory for graph
        ALLOCATE(iPath(NEl,0:NDets),stat=ierr)
        CALL LogMemAlloc('iPath',NEl*(NDets+1),4,this_routine,iPathTag)
        CALL IAZZERO(iPath,NEl*(NDets+1))
        ALLOCATE(Xij(0:NDets-1,0:NDets-1),stat=ierr)
        CALL LogMemAlloc('Xij',NDets*NDets,8,this_routine,XijTag)
        CALL AZZERO(Xij,NDets*NDets)
        ALLOCATE(Rhoii(0:NDets),stat=ierr)
        CALL LogMemAlloc('Rhoii',NDets+1,8,this_routine,RhoiiTag)
        CALL AZZERO(Rhoii,NDets+1)
        ALLOCATE(Rhoij(0:NDets,0:NDets),stat=ierr)
        CALL LogMemAlloc('Rhoij',(NDets+1)*(NDets+1),8*HElementSize,this_routine,RhoijTag)
        CALL AZZERO(Rhoij,(NDets+1)*(NDets+1)*HElementSize)
        ALLOCATE(Hijs(0:NDets),stat=ierr)
        CALL LogMemAlloc('Hijs',NDets+1,8*HElementSize,this_routine,HijsTag)
        CALL AZZERO(Hijs,(NDets+1)*HElementSize)

!The first and last determinants of iPath want to be the FDet...
        do i=1,NEl
            iPath(i,0)=FDet(i)
            iPath(i,NDets)=FDet(i)
        enddo

!ExcitGen is an array of pointers which points to the excitation generators for each vertex of the graph.
!Is used, and the memory is deallocated in the graph generation algorithm, so do not need to worry about it.
        ExcitGen(:)=0

        !Generate the initial graph...(Can speed this up as do not need to know prob, and are recalculating the rho matrix)
        !Arr sent instead of NMax since getting values straight from modules (not passed through)
        CALL Fmcpr4d2GenGraph(FDet,NEl,Beta,i_P,iPath,NDets,Xij,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,Alat,UMat,nTay,RhoEps,Rhoii,Rhoij,ECore,Seed,Hijs,nExcit,0,ExcitGen)

!Store initial graph determinants in the GraphDets array
        ALLOCATE(GraphDets(NDets,NEl),stat=ierr)
        CALL LogMemAlloc('GraphDets',NDets*NEl,4,this_routine,GraphDetsTag)
        CALL IAZZERO(GraphDets,NDets*NEl)
        do i=1,NDets
            do j=1,NEl
                GraphDets(i,j)=iPath(j,i-1)
                IF(GraphDets(i,j).eq.0) STOP 'Error in creating initial graph'
            enddo
        enddo
        DEALLOCATE(iPath)
        CALL LogMemDealloc(this_routine,iPathTag)

!Xij is the probability matrix - we do not need this information...
        DEALLOCATE(Xij)
        CALL LogMemDealloc(this_routine,XijTag)

!Rho matrix is stored as paths, and so do not need last column and row (should be same as first) - same with H elements
        ALLOCATE(GraphRhoMat(NDets,NDets),stat=ierr)
        CALL LogMemAlloc('GraphRhoMat',NDets*NDets,8*HElementSize,this_routine,GraphRhoMatTag)
        CALL AZZERO(GraphRhoMat,NDets*NDets*HElementSize)

!Trust rho matrix from Fmcpr4d2GenGraph - no need to recalculate...
        do i=1,NDets
            do j=1,NDets
                GraphRhoMat(i,j)=Rhoij(i-1,j-1)
            enddo
        enddo

!To make sure, put in diagonal elements separatly
        do i=1,NDets
            GraphRhoMat(i,i)=Rhoii(i-1)
        enddo

!First element is not put in for you...
        GraphRhoMat(1,1)=rhii

!        do i=1,NDets
!            do j=i,NDets
!                IF(i.eq.j) THEN
!                    diff=0
!                ELSE
!                    diff=-1
!                ENDIF
!                CALL CalcRho2(GraphDets(i,:),GraphDets(j,:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,diff,ECore)
!                GraphRhoMat(i,j)=rh
!                GraphRhoMat(j,i)=rh
!            enddo
!        enddo

        IF(DREAL(GraphRhoMat(1,1)).ne.(rhii%v)) THEN
!This could be because the elements haven't been divided by rhii - check value of rhii & then divide by it...
            STOP 'Rho matrix elements for initial graph incorrect'
        ENDIF

!Write out rho matrix - debugging
!        do i=1,NDets
!            do j=1,NDets
!                WRITE(6,"(E14.6)",advance='no') GraphRhoMat(i,j)
!            enddo
!            WRITE(6,*) ""
!        enddo

        ALLOCATE(HamElems(NDets),stat=ierr)
        CALL LogMemAlloc('HamElems',NDets,8*HElementSize,this_routine,HamElemsTag)
        CALL AZZERO(HamElems,NDets*HElementSize)

!Trust the fmcpr4d2gengraph Hij values. Recalculated ones are identical
        do i=1,NDets
            HamElems(i)=Hijs(i-1)
        enddo
        HamElems(1)=Hii

!        do i=1,NDets
!            HamElems(i)=GetHElement2(FDet,GraphDets(i,:),NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,-1,ECore)
!        enddo
        
        IF((HamElems(1)%v).ne.(Hii%v)) THEN
            STOP 'H elements for initial graph incorrect'
        ENDIF

!Return G_VMC_Pi to original value (Just in case it wants to be used later)
        G_VMC_Pi=OldImport
        
        DEALLOCATE(Rhoii)
        CALL LogMemDealloc(this_routine,RhoiiTag)
        DEALLOCATE(Rhoij)
        CALL LogMemDealloc(this_routine,RhoijTag)
        DEALLOCATE(Hijs)
        CALL LogMemDealloc(this_routine,HijsTag)
        DEALLOCATE(nExcit)
        CALL LogMemDealloc(this_routine,nExcitTag)

        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in ConstructInitialGraph'
        CALL TIHALT('ConsInitGraph',iSubInit)

    END SUBROUTINE ConstructInitialGraph

!This routine stocastically picks NDets new determinants stochastically from the vector obtained.
    SUBROUTINE PickNewDets()
        USE System , only : G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,Arr
        USE Calc , only : i_P,RhoEps
        USE Integrals , only : fck,nMax,nMsh,UMat,nTay
        USE Determinants , only : GetHElement2
        IMPLICIT NONE
        INTEGER :: ierr,i,j,k,NoVerts,AttemptDet(NEl),GrowGraphTag,IC,dist
        INTEGER :: Success,Failure,OrigDets,ExcitDets,Tries,iGetExcitLevel
        INTEGER , ALLOCATABLE :: GrowGraph(:,:)
        REAL*8 :: r,Ran2
        LOGICAL :: Attach,OriginalPicked
        TYPE(HElement) :: rh
        CHARACTER(len=*), PARAMETER :: this_routine='PickNewDets'

        GrowGraphTag=0

!Allocate Rho Matrix for growing graph
        ALLOCATE(GraphRhoMat(NDets,NDets),stat=ierr)
        CALL LogMemAlloc('GraphRhoMat',NDets*NDets,8*HElementSize,this_routine,GraphRhoMatTag)
        CALL AZZERO(GraphRhoMat,NDets*NDets*HElementSize)
        CALL AZZERO(HamElems,NDets*HElementSize)

!Initial determinant has to be HF determinant - ensure that this is stored in GraphDets(:,:)
        ALLOCATE(GrowGraph(NDets,NEl),stat=ierr)
        CALL LogMemAlloc('GrowGraph',NDets*NEl,4,this_routine,GrowGraphTag)
        CALL IAZZERO(GrowGraph,NEl*NDets)
        do i=1,NEl
            GrowGraph(1,i)=FDet(i)
        enddo
        GraphRhoMat(1,1)=rhii
        HamElems(1)=Hii

!NoVerts indicates the number of vertices in the graph so far
        NoVerts=1

!For interesting statistics, keep a record of successful/unsuccessful attachments, and 
!mean number of excitations away from HF
        Success=0
        Failure=0
        MeanExcit=0.D0
!Also keep a ratio of determinants which were in the graph before, and new ones added.
        OrigDets=0
        ExcitDets=0
!Allow a maximum average of ten attempts for every successfully attached determinant
        Tries=NDets*100000
        k=0

!Continue trying to build graph until fully constructed
        do while ((NoVerts.lt.NDets).and.(k.le.Tries))
            
            k=k+1
!Stochastically pick determinant to add to the graph from vectors
            r=RAN2(Seed)
!Set i=1, since we do not want to choose the first determinant of the original graph
            i=1

!First look through the normalised eigenvector*eigenvalue
            do while ((r.gt.0.D0).and.(i.lt.NDets))
                i=i+1
                r=r-((Eigenvector(i)%v)*(Eigenvector(i)%v))
            enddo

!If still not found, then look in Vector of excitations
            IF(r.gt.0.D0) THEN
                i=0
                do while ((r.ge.0.D0).and.(i.lt.TotExcits))
                    i=i+1
                    r=r-((ExcitsVector(i)%v)*(ExcitsVector(i)%v))
                enddo

!Error - cannot find determinant to attach in original graph, or its excitations
                IF(r.gt.0.D0) THEN
                    WRITE(6,*) "k = ", k
                    WRITE(6,*) "r = ", r
                    STOP 'Error in stochastic sampling in PickNewDets'
                ENDIF

!Determinant picked is from excitations...
                OriginalPicked=.false.
                do j=1,NEl
                    AttemptDet(j)=ExcitsDets(i,j)
                enddo

            ELSE
!Determinant picked is from original graph...
                OriginalPicked=.true.
                do j=1,NEl
                    AttemptDet(j)=GraphDets(i,j)
                enddo

            ENDIF

!Need to test whether graph is accepted - if not, they cycle around for another attempt at finding a valid determinant
            Attach=.false.
            DO i=1,NoVerts

!Check the number of excitations away from each other vertex in graph
                IC=IGetExcitLevel(AttemptDet(:),GrowGraph(i,:),NEl)
!                dist=IGetExcitLevel(FDet(:),AttemptDet(:),NEl)
!                IF(dist.gt.2) THEN
!                    WRITE(6,*) "Higher excitation selected ", dist
!                ENDIF
                IF(IC.eq.0) THEN
!Determinant is already in the graph - exit loop - determinant not valid
                    Attach=.false.
                    EXIT
                ENDIF

                !Determine connectivity to other determinants in the graph, if no connection yet found
                IF(.not.Attach) THEN
                    CALL CalcRho2(AttemptDet(:),GrowGraph(i,:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
                    IF(rh.agt.RhoEps) Attach=.true. 
                ENDIF

            ENDDO

!Attach new determinant to the list if it passes tests
            IF(Attach) THEN

!First check on attachment to HF, since wants to be stored in HamElems
                IC=IGetExcitLevel(AttemptDet,FDet,NEl)
!                IF(IC.gt.2) WRITE(6,*) "Higher excitation attached ", IC
     
                MeanExcit=MeanExcit+IC

!Only double excitations (or single? <- include for completness) of HF contribute
                IF(IC.eq.2) THEN!.or.(IC.eq.1)) THEN
                    HamElems(NoVerts+1)=GetHElement2(FDet,AttemptDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,IC,ECore)
                    CALL CalcRho2(FDet,AttemptDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
!Add rhoij contribution. No real need to set a rho epsilon here? We are diagonalising anyway...
                    GraphRhoMat(1,NoVerts+1)=rh
                    GraphRhoMat(NoVerts+1,1)=rh

                ELSE
!Added determinant not connected to root
                    HamElems(NoVerts+1)=HElement(0.D0)
                    GraphRhoMat(1,NoVerts+1)=HElement(0.D0)
                    GraphRhoMat(NoVerts+1,1)=HElement(0.D0)
                ENDIF

!Run through rest of excitations in the graph testing contributions and adding to rho matrix
                do i=2,NoVerts
                    CALL CalcRho2(GrowGraph(i,:),AttemptDet(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,-1,ECore)
                    GraphRhoMat(i,NoVerts+1)=rh
                    GraphRhoMat(NoVerts+1,i)=rh
                enddo

!Include diagonal rho matrix element for chosen determinant
                CALL CalcRho2(AttemptDet,AttemptDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,0,ECore)
                GraphRhoMat(NoVerts+1,NoVerts+1)=rh

!Add determinant to growing graph
                do i=1,NEl
                    GrowGraph(NoVerts+1,i)=AttemptDet(i)
                enddo
                NoVerts=NoVerts+1
!                WRITE(6,"A,I5") "Vertex Added - ",NoVerts
!                CALL FLUSH(6)

                Success=Success+1
                IF(OriginalPicked) THEN
                    OrigDets=OrigDets+1
                ELSE
                    ExcitDets=ExcitDets+1
                ENDIF

            ELSE

!Determinant chosen not attached to growing graph
                Failure=Failure+1

            ENDIF

!Attempt to find another determinant to attach
        enddo

!Find average number of excitations away from HF determinant
        MeanExcit=MeanExcit/(NDets-1)

!Test if graph growing was successful. If not, possible reasons are that the space is too small for
!such a large graph to be grown easily - reduce NDets, or increase time to search for determinants.
        IF(NoVerts.ne.NDets) STOP 'Error in attaching determinants to new graph'
        IF(Success.ne.(NDets-1)) STOP 'Error in attaching determinants to new graph 2'

!Determine success ratio for attachment of determinants into new graph
        SucRat=(Success+0.D0)/(Success+Failure+0.D0)
        Orig_Graph=((OrigDets+0.D0)/(ExcitDets+OrigDets+0.D0))*100
        WRITE(6,"(A,F9.5)") "New graph created. Success ratio for adding determinants: ",SucRat
        WRITE(6,"(A,F9.4)") "Percentage of determinants found in original graph: ",Orig_Graph

!Replace original list of determinants in graph with new list
        CALL IAZZERO(GraphDets,NDets*NEl)
        do i=1,NDets
            do j=1,NEl
                IF(GrowGraph(i,j).eq.0) STOP 'Error in allocating determinants to graph'
                GraphDets(i,j)=GrowGraph(i,j)
            enddo
        enddo

!Deallocate temporary list of determinants in graph, list of all excitations, and the ExcitsVector vector
        DEALLOCATE(ExcitsVector)
        CALL LogMemDealloc(this_routine,ExcitsVectorTag)
        DEALLOCATE(ExcitsDets)
        CALL LogMemDealloc(this_routine,ExcitsDetsTag)
        DEALLOCATE(GrowGraph)
        CALL LogMemDealloc(this_routine,GrowGraphTag)

        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in PickNewDets'

    END SUBROUTINE PickNewDets

!This is used to normalise the vector of current and excited determinants from which we are going to pick the next graph.
!The vector is spread between the arrays for the original eigenvector (which needs to be multiplied by corresponding eigenvalue)
!and the ExcitsVector. The first element of the eigenvector should not be included in the normalisation, as
!it cannot be picked - the HF is always in each graph.
    SUBROUTINE NormaliseVector()
        IMPLICIT NONE
        INTEGER :: i
        REAL*8 :: Stay,Move
        TYPE(HElement) :: Norm,Norm1,Norm2

        Biasing=.true.

!The bias towards the determinants already in the graph is given by the largest eigenvector, multiplied by its eigenvalue.
!Since we no longer need the largest eigenvector, we can multiply its elements by its eigenvalue
        do i=2,NDets
            Eigenvector(i)=Eigenvector(i)*HElement(Eigenvalue)
        enddo

!An addititional bias against the determinants already in the graph can be designed.
        IF(Biasing) THEN
!GraphBias is the probability of picking a determinant which is already in the graph.
!Remember that this is not actually a strict probability, since the way that the graph is grown,
!means that some determinants will be automatically preferred over others.
            
            IF((GraphBias.gt.1.D0).or.(GraphBias.lt.0.D0)) THEN
                STOP 'Value for Graphbias must be between 1 and 0'
            ENDIF

            Norm1=HElement(0.D0)
            do i=2,NDets
                Norm1=Norm1+(Eigenvector(i)*Eigenvector(i))
            enddo
            Norm1=HElement(SQRT((Norm1%v)/GraphBias))

!Divide elements of the eigenvector by the normalisation
            Stay=0.D0
            do i=2,NDets
                Eigenvector(i)=Eigenvector(i)/Norm1
                Stay=Stay+((Eigenvector(i)%v)**2)
            enddo

!Now find normalisation for the excitations
            Norm2=HElement(0.D0)
            do i=1,TotExcits
                Norm2=Norm2+(ExcitsVector(i)*ExcitsVector(i))
            enddo
            Norm2=HElement(SQRT((Norm2%v)/(1.D0-GraphBias)))

!Divide elements of ExcitsVector by new normalisation
            Move=0.D0
            do i=1,TotExcits
                ExcitsVector(i)=ExcitsVector(i)/Norm2
                Move=Move+((ExcitsVector(i)%v)**2)
            enddo

        ELSE
!No biasing towards excitations

!We need to find the normalisation constant, given by the sum of the squares of all the elements
            Norm=HElement(0.D0)
!First, sum the squares of the original determinants in the graph
            do i=2,NDets
                Norm=Norm+(Eigenvector(i)*Eigenvector(i))
            enddo

!Then, sum the squares of the vector for the excitations
            do i=1,TotExcits
                Norm=Norm+(ExcitsVector(i)*ExcitsVector(i))
            enddo

            Norm=HElement(SQRT(Norm%v))

!Once the normalisation is found, all elements need to be divided by it.
!Stay is the total probability of staying with original graph
            Stay=0.D0
            Move=0.D0
            do i=2,NDets
                Eigenvector(i)=Eigenvector(i)/Norm
                Stay=Stay+((Eigenvector(i)%v)**2)
            enddo
            do i=1,TotExcits
                ExcitsVector(i)=ExcitsVector(i)/Norm
                Move=Move+((ExcitsVector(i)%v)**2)
            enddo

        ENDIF

        PStay=Stay
!        WRITE(6,*) "Probability of staying at original determinants: ", Stay
        WRITE(6,*) "Probability of Moving: ", Move
        WRITE(6,*) "Total Probability: ", Stay+Move
!        WRITE(6,*) "Normalisation constant for propagation vector: ", Norm

    END SUBROUTINE NormaliseVector

!This routine creates a new vector (ExcitsVector), which is formed by the matrix product of the rho matrix for the new excitations
!(ConnectionsToExcits) and the original eigenvector.
!This is the rho operator for the excitations acting on the original eigenstate to show the probabilities of being in the new
!excitations. Since each excitation can only be attached to one vertex of the graph, there is only one element per
!row of the operator, and only one multiplication is needed per excitation.
    SUBROUTINE CreateExcitsVector()
        IMPLICIT NONE
        INTEGER :: ierr,NoExcitsAttached,Element,i,j
        CHARACTER(len=*), PARAMETER :: this_routine='CreateExcitsVector'

!Allocate space to hold this new vector - could get away without allocating any more memory by simply multiplying Connections
!array by needed element of eigenvector...?
        ALLOCATE(ExcitsVector(TotExcits),stat=ierr)
        CALL LogMemAlloc('ExcitsVector',TotExcits,8*HElementSize,this_routine,ExcitsVectorTag)
        CALL AZZERO(ExcitsVector,TotExcits*HElementSize)

!Separatly deal with first vertex in graph for clarity
        do j=1,NoExcits(1)
!Multiply rho element to excitation by eigenvector component corresponding to determinant being excited
            ExcitsVector(j)=ConnectionsToExcits(j)*Eigenvector(1)
        enddo

!Cycle through all remaining vertices in graph
        do i=2,NDets

            NoExcitsAttached=NoExcits(i)-NoExcits(i-1)
!Cycle through all determinants attached to vertex i
            do j=1,NoExcitsAttached

!Corresponding element in 'ConnectionsToExcits' and 'ExcitsVector' is calculated
                Element=j+NoExcits(i-1)
                ExcitsVector(Element)=ConnectionsToExcits(Element)*Eigenvector(i)

            enddo

        enddo

!Deallocate NoExcits Array and ConnectionsToExcits array
        DEALLOCATE(ConnectionsToExcits)
        CALL LogMemDealloc(this_routine,ConnectionsToExcitsTag)
        DEALLOCATE(NoExcits)
        CALL LogMemDealloc(this_routine,NoExcitsTag)

        IF(Element.ne.TotExcits) THEN
            STOP 'Error in counting in CreateExcitsVector'
        ENDIF

        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in CreateExcitsVector'
    END SUBROUTINE CreateExcitsVector

!This subroutine will go through all determinants, find all excitations, and then calculate the rho element to
!the determinant in the graph to which it is connected.
    SUBROUTINE FindConnections()
        USE System, only: G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,Arr
        USE Calc , only : i_P
        USE Integrals, only: fck,nMax,nMsh,UMat,nTay
        IMPLICIT NONE
        TYPE(HElement) :: rh
        INTEGER , SAVE :: iSubConns
        INTEGER :: ierr,i,j,DetCurr(NEl),nJ(NEl),nStore(6),iMaxExcit,nExcitMemLen
        INTEGER :: nExcitTag,iExcit,ExcitCurr,dist,iGetExcitLevel
        INTEGER , ALLOCATABLE :: nExcit(:)
        CHARACTER(len=*), PARAMETER :: this_routine='FindConnections'

        CALL TISET('FindConns',iSubConns)

        nExcitTag=0

!Allocate memory to hold connections, and form of the determinants of the excitations
        ALLOCATE(ConnectionsToExcits(TotExcits),stat=ierr)
        CALL LogMemAlloc('ConnectionsToExcits',TotExcits,8*HElementSize,this_routine,ConnectionsToExcitsTag)
        CALL AZZERO(ConnectionsToExcits,TotExcits*HElementSize)
        ALLOCATE(ExcitsDets(TotExcits,NEl),stat=ierr)
        CALL LogMemAlloc('ExcitsDets',TotExcits*NEl,4,this_routine,ExcitsDetsTag)
        CALL IAZZERO(ExcitsDets,TotExcits*NEl)

!Cycle over all vertices in graph
        ExcitCurr=0
        do i=1,NDets

!Find determinant form for current vertex
            do j=1,NEl
                DetCurr(j)=GraphDets(i,j)
            enddo

!Setup excitation generators again for this determinant
            iMaxExcit=0
            CALL IAZZERO(nStore,6)
            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
            ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
            CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag)
            CALL IAZZERO(nExcit,nExcitMemLen)
            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,3)

            IF(i.eq.1) THEN
                IF(iMaxExcit.ne.NoExcits(i)) STOP 'Error in counting in FindConnections'
            ELSE
                IF(iMaxExcit.ne.(NoExcits(i)-NoExcits(i-1))) STOP 'Error in counting in FindConnections'
            ENDIF

!Cycle through all excitations of each determinant
        lp: do while(.true.)
                CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,iExcit,0,nStore,3)
                IF(nJ(1).eq.0) exit lp
                ExcitCurr=ExcitCurr+1
                CALL CalcRho2(DetCurr,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,iExcit,ECore)
!                Dist=IGetExcitLevel(FDet,nJ,NEl)
!                IF(Dist.gt.4) THEN
!                    WRITE(6,*) "Higher than double excitation found - attached to det:", DetCurr
!                    WRITE(6,*) "Determinant is: ", nJ
!                    WRITE(6,*) Dist," fold excitation"
!                ENDIF

!Store path of determinant in ExcitsDets, and the rho elements in ConnectionsToExcits
                ConnectionsToExcits(ExcitCurr)=rh
                do j=1,NEl
                    ExcitsDets(ExcitCurr,j)=nJ(j)
                enddo
            enddo lp

            IF(ExcitCurr.ne.NoExcits(i)) STOP 'Incorrect counting in FindConnections'

!Deallocate excitation generators
            DEALLOCATE(nExcit)
            CALL LogMemDealloc(this_routine,nExcitTag)

        enddo

        IF(ExcitCurr.ne.TotExcits) STOP 'Incorrect counting in FondConnections'

        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in FindConnections'
        CALL TIHALT('FindConns',iSubConns)

    END SUBROUTINE FindConnections

!This routine initialises and then destroys excitation generators, in order to calculate the space required for all excitations.
!This could improve on space if all excitations were explicitly calculated, though would take longer.
    SUBROUTINE CountExcits()

        USE System, only : G1,nBasis,nBasisMax
        IMPLICIT NONE
        INTEGER :: ierr,i,j,nStore(6),DetCurr(NEl),nJ(NEl),iMaxExcit,nExcitMemLen
        INTEGER :: nExcitTag
        CHARACTER(len=*), PARAMETER :: this_routine='CountExcits'
        INTEGER , ALLOCATABLE :: nExcit(:)

        TotExcits=0
        nExcitTag=0
!Allocate Memory for NoExcits array
        ALLOCATE(NoExcits(NDets),stat=ierr)
        CALL LogMemAlloc('NoExcits',NDets,4,this_routine,NoExcitsTag)
        CALL IAZZERO(NoExcits,NDets)

!Cycle over all vertices in graph
        do i=1,NDets
            
!Find determinant form for current vertex
            do j=1,NEl
                DetCurr(j)=GraphDets(i,j)
            enddo

!Create excitation generator
            iMaxExcit=0
            CALL IAZZERO(nStore,6)
            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
            ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
            CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag)
            CALL IAZZERO(nExcit,nExcitMemLen)
            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,3)
            
!Store number of excitations from each determinant cumulativly.
            IF(i.eq.1) THEN
                NoExcits(i)=iMaxExcit
            ELSE
                NoExcits(i)=NoExcits(i-1)+iMaxExcit
            ENDIF

            TotExcits=TotExcits+iMaxExcit

!Destroy excitation generator
            DEALLOCATE(nExcit)
            CALL LogMemDealloc(this_routine,nExcitTag)
            CALL IAZZERO(nStore,6)
        
        enddo

!Perform check that all excitations accounted for
        IF(NoExcits(NDets).ne.TotExcits) THEN
            STOP 'Error in counting excits in CountExcits'
        ENDIF

        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in CountExcits'
    END SUBROUTINE CountExcits

!This is a routine to find the energy of the graph by diagonalisation, and return as well its eigenvalues and largest eigenvector. Deallocate the RhoMatrix when done.
    SUBROUTINE DiagGraphMorph()
        IMPLICIT NONE
        INTEGER, SAVE :: iSubDiag
        INTEGER :: Info,ierr,i
        REAL*8 , ALLOCATABLE :: Work(:),Eigenvalues(:)
        INTEGER :: WorkTag,EigenvaluesTag
        CHARACTER(len=*), PARAMETER :: this_routine='DiagGraphMorph'

        CALL TISET('DiagGraphMorph',iSubDiag)

        WorkTag=0
        EigenvaluesTag=0

        ALLOCATE(Eigenvalues(NDets),stat=ierr)
        CALL LogMemAlloc('Eigenvalues',NDets,8,this_routine,EigenvaluesTag)
        CALL AZZERO(Eigenvalues,NDets)

!Workspace needed for diagonaliser
        ALLOCATE(Work(3*NDets),stat=ierr)
        CALL LogMemAlloc('WORK',NDets*3,8,this_routine,WorkTag)

!Diagonalise...
!Eigenvalues now stored in ascending order
        CALL DSYEV('V','U',NDets,GraphRhoMat,NDets,Eigenvalues,Work,3*NDets,Info)
        IF(Info.ne.0) THEN
            WRITE(6,*) "DYSEV error in DiagGraphMorph: ",Info
            STOP
        ENDIF

        DEALLOCATE(Work)
        CALL LogMemDealloc(this_routine,WorkTag)

!Store largest eigenvector - last column of GraphRhoMat (zero it)
        CALL AZZERO(Eigenvector,NDets*HElementSize)

        do i=1,NDets
            Eigenvector(i)=GraphRhoMat(i,NDets)
        enddo

!Also need to save the largest Eigenvalue
        Eigenvalue=Eigenvalues(NDets)

!Deallocate RhoMatrix and Eigenvalues (not needed)
        DEALLOCATE(Eigenvalues)
        CALL LogMemDealloc(this_routine,EigenvaluesTag)
        DEALLOCATE(GraphRhoMat)
        CALL LogMemDealloc(this_routine,GraphRhoMatTag)

!Find weight and energy of graph.
!There is no beta-dependance, so only largest eigenvector needed.
!Note - since we are not having a beta-dependance, 'weight' takes on a slightly different
!meaning. Here, the weight is simply the square of the first element of the largest eigenvector,
!i.e. the magnitude of the projection of the graph back onto the HF...
        SI=(Eigenvector(1)%v)*(Eigenvector(1)%v)
        DLWDB=0.D0
        do i=2,NDets
            DLWDB=DLWDB+(HamElems(i)%v)*(Eigenvector(i)%v)
        enddo
        DLWDB=(DLWDB/(Eigenvector(1)%v))+(HamElems(1)%v)
        
        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in DiagGraphMorph'
        CALL TIHALT('DiagGraphMorph',iSubDiag)

    END SUBROUTINE DiagGraphMorph

END MODULE GraphMorph
