!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!   TO DO   !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!Fix Lanczos code
!Fix MoveDets code
!New initial graph growing routine which adds dets in order until complete (not stochastically)
!Go to larger graphs, and look at possible convergence
!Look into MC sampling of rest of excitation space to aid convergence

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
    USE Calc , only : Iters,NDets,GraphBias,TBiasing,NoMoveDets,TMoveDets,TInitStar
    USE Calc , only : TNoSameExcit,TLanczos,TMaxExcit,iMaxExcitLevel,TOneExcitConn
    USE Calc , only : TSinglesExcitSpace,TMCExcitSpace,NoMCExcits,TGrowInitGraph
    USE Logging , only : TDistrib
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

!This is needed if we are using the MoveDets graph growing algorithm to store a copy of the rho matrix in
    REAL*8 , ALLOCATABLE :: CopyRhoMat(:,:)
    INTEGER :: CopyRhoMatTag=0

!If TDistrib is on, then this will show the distribution of determinants among the excitations
    INTEGER , ALLOCATABLE :: Distribs(:,:)
    INTEGER :: DistribsTag=0

!If TNoSameExcit is on, then GraphExcitLevel stores the Excitation levels of the determinants in the graph
!This is because we do not allow connections between determinants of the same excitation level.
!We also require the array for TOneExcitConn for the same reasons
    INTEGER , ALLOCATABLE :: GraphExcitLevel(:)
    INTEGER :: GraphExcitLevelTag=0

!The rho matrix and Hamiltonian element for the HF determinant
    TYPE(HElement) :: rhii,Hii

!These are the weight and energy of the current graph respectively
    REAL*8 :: SI,DLWDB

!This is the seed for the random numbers needed in the routines
    INTEGER :: Seed

!This is the iteration of the GraphMorph that we are currently on...
    INTEGER :: Iteration

!This is needed in case a disconnected graph is generated - if so, a full new graph will be generated in the first MC
    LOGICAL :: ReturntoTMoveDets

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
        REAL*8 :: LowestE,BestSI
        INTEGER :: ierr,i
        CHARACTER(len=*), PARAMETER :: this_routine='MorphGraph'

        OPEN(63,file='MorphStats',Status='unknown')
        IF(TDistrib) THEN
            OPEN(64,file='MorphDistrib',Status='unknown')
        ENDIF
        
        IF(HElementSize.ne.1) STOP 'Only real orbitals allowed in GraphMorph so far'

        IF(TMoveDets) THEN
            WRITE(6,*) "***************************************"
            WRITE(6,*) "*** WARNING - MOVEDETS MAY NOT WORK ***"
            WRITE(6,*) "***************************************"
        ENDIF
        IF(TMoveDets.and.TOneExcitConn) THEN
            WRITE(6,*) "TOneExcitConn is not yet implimented in TMoveDets"
            STOP "TOneExcitConn is not yet implimented in TMoveDets"
        ENDIF
        IF(TLanczos) THEN
            WRITE(6,*) "Sorry, but Lanczos diagonalisation not yet working for GraphMorph"
            STOP "Sorry, but Lanczos diagonalisation not yet working for GraphMorph"
        ENDIF
        IF(TMCExcitSpace) THEN
            WRITE(6,"(A,I10,A)") "Searching Excitations of each graph stochastically with ",NoMCExcits, " per determinant"
        ENDIF

        ReturntoTMoveDets=.false.
!Initialise random number generator
        Seed=G_VMC_Seed

!Find rho_ii value, which all rho elements will be divided by, and Hii value
        CALL CalcRho2(FDet,FDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhii,nTay,0,ECore)
        Hii=GetHElement2(FDet,FDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)

!It is first necessary to construct the original graph to work with.
        IF(TInitStar) THEN
!Let the initial graph be a complete connected star graph.
            WRITE(6,"(A)") "Constructing random initial star graph to morph from."
            CALL ConstructInitialStarGraph()
        ELSEIF(TGrowInitGraph) THEN
            WRITE(6,"(A,I8,A)") "Constructing initial graph in increasing excitation space with ",NDets, " vertices."
            CALL ConstructExcitsInitGraph()
        ELSE
            WRITE(6,"(A,I8,A)") "Constructing random fully connected initial graph to morph from with ", NDets, " vertices."
            CALL ConstructInitialGraph()
        ENDIF

        IF(TDistrib) THEN
            WRITE(64,"(A8)",advance='no') "1"
            do i=1,NEl
                WRITE(64,"(I8)",advance='no') Distribs(i,1)
            enddo
            WRITE(64,*) ""
        ENDIF
        
!Allocate space for largest Eigenvector for various graphs
        ALLOCATE(Eigenvector(NDets),stat=ierr)
        CALL LogMemAlloc('Eigenvector',NDets,8*HElementSize,this_routine,EigenvectorTag)

        WRITE(63,*) "Iteration  Energy  P_stay  Orig_Graph  SucRat   MeanExcit"
        IF(TBiasing) THEN
            WRITE(6,"(A,F10.7)") "Graph growing biased towards original determinants with probability ", GraphBias
        ELSEIF(TMoveDets) THEN
            WRITE(6,"(A,I3,A)") "Choosing new graph by moving ", NoMoveDets," determinants in previous graph to their excitations..."
        ENDIF

!Once the graph is found, the loop over the morphing iterations can begin
        do Iteration=1,Iters

            WRITE(6,*) ""
            WRITE(6,"(A,I4,A)") "Starting Iteration ", Iteration, " ..."
            WRITE(6,*) ""
            CALL FLUSH(6)

            IF((Iteration.eq.2).and.ReturntoTMoveDets) THEN
!If one iteration of regrowing graphs is used, then return to the moving dets algorithm                
                WRITE(6,*) "Returning to MoveDets algorithm"
                TMoveDets=.true.
                ReturntoTMoveDets=.false.
                TBiasing=.false.
            ENDIF

!The graph is first diagonalised, and the energy of the graph found, along with the largest eigenvector and eigenvalues.
            IF(TLanczos) THEN
                CALL DiagGraphLanc()
            ELSE
                CALL DiagGraphMorph()
            ENDIF

            WRITE(6,"(A,2G20.12)") "Weight and Energy of current graph is: ", SI, DLWDB
            IF(Iteration.eq.1) THEN
                LowestE=DLWDB
                BestSI=SI
            ELSE
                IF(DLWDB.lt.LowestE) THEN
                    LowestE=DLWDB
                    BestSI=SI
                ENDIF
            ENDIF
            CALL FLUSH(6)

!Write out stats
            IF(Iteration.eq.1) THEN
                WRITE(63,"(I10,G20.12,3A20,G20.12)") Iteration,DLWDB,"N/A","N/A","N/A",MeanExcit
            ELSE
                WRITE(63,"(I10,5G20.12)") Iteration,DLWDB,PStay,Orig_Graph,SucRat,MeanExcit
            ENDIF
            CALL FLUSH(63)

!Excitation generators are initialised for each of the determinants in the graph, and the total number of possible
!connected determinants calculated. Memory allocated for the ensuing calculation.
            CALL CountExcits()
!            WRITE(6,*) "Fraction of space which is space of excitations is: ", TotExcits/(TotExcits+NDets)
            WRITE(6,*) "Total number of determinants available from current graph is: ",TotExcits
!            CALL FLUSH(6)

!Run through each determinant in the graph, calculating the excitations, storing them, and the rho elements to them.
            CALL FindConnections()

!Multiply the correct rho element, with the correct element of the eigenvector, to create the value for the 
!tendancy to move to that excited determinant - store this in ExcitsVector
            CALL CreateExcitsVector()

            IF(TMoveDets) THEN
!MoveDets means that a new graph-growing algorithm is used, where NoMoveDets determinants are selected from the current graph, and replaced
!by the same number of determinants from its excitations.

!New arrays created - the inverse array from the original set of determinants, and the normalised array for the excitations
                CALL NormaliseVectorSep()

                CALL MoveDetsGraph()

            ELSE

!Add the original eigenvector*eigenvalue to the list of determinants - now have vector (contained in Eigenvector*eigenvalue and ExcitsVector) with all
!determinants in graph, and all possible determinants to excite to. This needs normalising.
                CALL NormaliseVector()

!Pick NDets new excitations stocastically from normalised list of determinants with probability |c|^2. Ensure connections,
!allocate and create rho matrix for new graph. Deallocate info for old graph.
                WRITE(6,*) "Choosing new graph stochastically from previous graph and its excitations..."
                CALL FLUSH(6)
                CALL PickNewDets()

            ENDIF
        
            IF(TDistrib) THEN
                WRITE(64,"(I8)",advance='no') Iteration+1
                do i=1,NEl
                    WRITE(64,"(I8)",advance='no') Distribs(i,Iteration+1)
                enddo
                WRITE(64,*) ""
                CALL FLUSH(64)
            ENDIF

!Once graph is fully constructed, the next iteration can begin.
        enddo
        WRITE(6,*) ""

!Diagonalise final graph, and find weight & energy. 
!There is no beta-dependance, so only largest eigenvector and value needed.
        CALL DiagGraphMorph()
        WRITE(6,"(A,2G20.12)") "Weight and Energy of final graph is: ", SI, DLWDB
        IF(DLWDB.lt.LowestE) THEN
            LowestE=DLWDB
            BestSI=SI
        ENDIF

!Deallocate info...
        DEALLOCATE(Eigenvector)
        CALL LogMemDealloc(this_routine,EigenvectorTag)
        DEALLOCATE(GraphDets)
        CALL LogMemDealloc(this_routine,GraphDetsTag)
        DEALLOCATE(HamElems)
        CALL LogMemDealloc(this_routine,HamElemsTag)
        IF(TDistrib) THEN
            DEALLOCATE(Distribs)
            CALL LogMemDealloc(this_routine,DistribsTag)
        ENDIF
        IF(TNoSameExcit.or.TOneExcitConn) THEN
            DEALLOCATE(GraphExcitLevel)
            CALL LogMemDealloc(this_routine,GraphExcitLevelTag)
        ENDIF

        Weight=HDElement(BestSI-1.D0)
        Energyxw=HDElement((LowestE*BestSI)-Hii%v)

        Close(63)
        IF(TDistrib) THEN
            Close(64)
        ENDIF
        
        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in ConstructInitialGraph'

    END SUBROUTINE MorphGraph

!This routine constructs an initial graph in a non-stocastic manner by adding determinants from increasing excitation levels
    SUBROUTINE ConstructExcitsInitGraph()
        USE System , only : G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,Arr
        USE Calc , only : i_P,RhoEps
        USE Integrals , only : fck,nMax,nMsh,UMat,nTay
        USE Determinants , only : GetHElement2
        IMPLICIT NONE
        INTEGER :: nStore(6),exFlag,nExcitMemLen,iMaxExcit,nJ(NEl)
        INTEGER , ALLOCATABLE :: nExcit(:)
        INTEGER :: nExcitTag=0
        TYPE(HElement) :: rh
        INTEGER :: iExcit,excitcount,i,j,k,IC,ICRoot
        CHARACTER(len=*), PARAMETER :: this_routine='ConstructExcitsInitGraph'
        INTEGER :: ierr,iSubInitExcit,Root,RootDet(NEl),iGetExcitLevel
        LOGICAL :: SameDet,Connection

        CALL TISET('InitExcitGraph',iSubInitExcit)

!Allow single and double excitations
        exFlag=3

        IF(TDistrib) THEN
!We want to record the number of each excitation level for each determinant in the graph
            ALLOCATE(Distribs(NEl,Iters+1),stat=ierr)
            CALL LogMemAlloc('Distribs',NEl*(Iters+1),4,this_routine,DistribsTag)
            CALL IAZZERO(Distribs,NEl*(Iters+1))
        ENDIF

!Allocate space for graph determinants and rho matrix
        ALLOCATE(GraphDets(NDets,NEl),stat=ierr)
        CALL LogMemAlloc('GraphDets',NDets*NEl,4,this_routine,GraphDetsTag)
        CALL IAZZERO(GraphDets,NDets*NEl)
        ALLOCATE(GraphRhoMat(NDets,NDets),stat=ierr)
        CALL LogMemAlloc('GraphRhoMat',NDets*NDets,8*HElementSize,this_routine,GraphRhoMatTag)
        CALL AZZERO(GraphRhoMat,NDets*NDets*HElementSize)
        ALLOCATE(HamElems(NDets),stat=ierr)
        CALL LogMemAlloc('HamElems',NDets,8*HElementSize,this_routine,HamElemsTag)
        CALL AZZERO(HamElems,NDets)
        IF(TNoSameExcit.or.TOneExcitConn) THEN
            ALLOCATE(GraphExcitLevel(NDets),stat=ierr)
            CALL LogMemAlloc('GraphExcitLevel',NDets,4,this_routine,GraphExcitLevelTag)
            CALL IAZZERO(GraphExcitLevel,NDets)
        ENDIF

!Put HF Determinant into first element of GraphDets
        do i=1,NEl
            GraphDets(1,i)=FDet(i)
        enddo
        GraphRhoMat(1,1)=rhii
        HamElems(1)=Hii
        MeanExcit=0.D0

        i=1
        Root=1
        do while(i.lt.NDets)
!Cycle through all excitations consecutivly, adding them where possible

            IF(Root.gt.i) THEN
                WRITE(6,*) "Error - trying to make an unavailable determinant root"
                STOP "Error - trying to make an unavailable determinant root"
            ENDIF

!Let root of excitation generator be next available determinant
            do j=1,NEl
                RootDet(j)=GraphDets(Root,j)
                IF(RootDet(j).eq.0) THEN
                    WRITE(6,*) "Incorrect assignment of root determinant"
                    STOP "Incorrect assignment of root determinant"
                ENDIF
            enddo

!Setup excitation generator for new root
            CALL IAZZERO(nStore,6)
            CALL GenSymExcitIt2(RootDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
            ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
            CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag)
            nExcit(1)=0
            CALL GenSymExcitIt2(RootDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)

            do while(i.lt.NDets)

                CALL GenSymExcitIt2(RootDet,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,iExcit,0,nStore,exFlag)
                IF(nJ(1).eq.0) THEN
!In the next sweep, look at the next determinant being the root
                    Root=Root+1
!Deallocate excitation generators
                    DEALLOCATE(nExcit)
                    CALL LogMemDealloc(this_routine,nExcitTag)
                    EXIT
                ENDIF

                Connection=.false.
                do j=1,i
!Look through determinants already in graph to check for double counting and connections
                    SameDet=.true.
                    do k=1,NEl
                        IF(nJ(k).ne.GraphDets(j,k)) THEN
                            SameDet=.false.
                        ENDIF
                    enddo
!If we have found the determinant somewhere else in the graph, then discount it
                    IF(SameDet) EXIT

                    IF(.not.Connection) THEN
                        CALL CalcRho2(GraphDets(j,:),nJ(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,-1,ECore)
                        IF(rh.agt.0.D0) THEN
!A connection has been found - we do not need to look for any others
                            Connection=.true.
                        ENDIF
                    ENDIF
                enddo

!                IF(.not.Connection) WRITE(6,*) "A determinant has been created and not attached"

                IF(Connection.and.(.not.SameDet)) THEN
!A valid determinant has been found - add it
                    i=i+1
                    do j=1,NEl
                        GraphDets(i,j)=nJ(j)
                    enddo

!Find diagonal rho matrix element, and connection to HF for HElement
                    CALL CalcRho2(nJ,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,0,ECore)
                    GraphRhoMat(i,i)=rh
                    IF(Root.eq.1) THEN
                        ICRoot=iExcit
                    ELSE
                        ICRoot=iGetExcitLevel(FDet,nJ,NEl)
                    ENDIF
                    MeanExcit=MeanExcit+(ICRoot+0.D0)
                    IF(ICRoot.gt.2) THEN
                        HamElems(i)=HElement(0.D0)
                        GraphRhoMat(i,1)=HElement(0.D0)
                        GraphRhoMat(1,i)=HElement(0.D0)
                    ELSE
                        HamElems(i)=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,ICRoot,ECore)
                        CALL CalcRho2(FDet(:),nJ(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,ICRoot,ECore)
                        GraphRhoMat(i,1)=rh
                        GraphRhoMat(1,i)=rh
                    ENDIF
!Store the excitation level from HF if we are removing excitations from the same level
                    IF(TNoSameExcit.or.TOneExcitConn) GraphExcitLevel(i)=ICRoot

                    IF(TDistrib) THEN
!Need to add excitation level of determinant. Also, store the excitation level of the determinants in the current graph
                        Distribs(ICRoot,1)=Distribs(ICRoot,1)+1
                    ENDIF

                    do j=2,i-1
!Need to find connection to all other determinants
            
                        IF(TNoSameExcit) THEN
!Don't allow connections between excitations of the same level

                            IF(GraphExcitLevel(j).ne.ICRoot) THEN
                                IC=iGetExcitLevel(GraphDets(j,:),nJ(:),NEl)
                                CALL CalcRho2(GraphDets(j,:),nJ(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
                                GraphRhoMat(i,j)=rh
                                GraphRhoMat(j,i)=rh
                            ELSE
                                GraphRhoMat(i,j)=HElement(0.D0)
                                GraphRhoMat(j,i)=HElement(0.D0)
                            ENDIF

                        ELSEIF(TOneExcitConn) THEN
!Only allow connections between determinants which differ in excitation level by one.

                            IF((ABS(GraphExcitLevel(j)-ICRoot)).eq.1) THEN
                                IC=iGetExcitLevel(GraphDets(j,:),nJ(:),NEl)
                                CALL CalcRho2(GraphDets(j,:),nJ(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
                                GraphRhoMat(i,j)=rh
                                GraphRhoMat(j,i)=rh
                            ELSE
                                GraphRhoMat(i,j)=HElement(0.D0)
                                GraphRhoMat(j,i)=HElement(0.D0)
                            ENDIF

                        ELSE
                            
                            IC=iGetExcitLevel(GraphDets(j,:),nJ(:),NEl)
!Fully connect the graph
                            CALL CalcRho2(GraphDets(j,:),nJ(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
                            GraphRhoMat(i,j)=rh
                            GraphRhoMat(j,i)=rh

                        ENDIF

!End of searching for connections
                    enddo
!Whether we connect the determinant created to the graph
                ENDIF
!End of loop for a given root
            enddo
!End of searching for excitation - all found
        enddo
                    
!Deallocate excitation generators
        DEALLOCATE(nExcit)
        CALL LogMemDealloc(this_routine,nExcitTag)

        IF(i.ne.NDets) THEN
            WRITE(6,*) "Graph Growing routine has not found all determinant requested"
            STOP "Graph Growing routine has not found all determinant requested"
        ENDIF

        WRITE(6,"(A,I4)") "Number of roots needed to create initial graph = ", Root

        MeanExcit=MeanExcit/(NDets-1)

        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in ConstructExcitsInitGraph'
        CALL TIHALT('InitExcitGraph',iSubInitExcit)

    END SUBROUTINE ConstructExcitsInitGraph


!This routine constructs the complete connected star graph as the initial graph to begin morphing from
    SUBROUTINE ConstructInitialStarGraph()
        USE System , only : G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,Arr
        USE Calc , only : i_P,RhoEps
        USE Integrals , only : fck,nMax,nMsh,UMat,nTay
        USE Determinants , only : GetHElement2
        IMPLICIT NONE
        INTEGER :: ierr,iSubInitStar,nStore(6),exFlag,nExcitMemLen,iMaxExcit,nJ(NEl)
        INTEGER , ALLOCATABLE :: nExcit(:)
        INTEGER :: nExcitTag=0
        TYPE(HElement) :: rh
        INTEGER :: iExcit,excitcount,i,j
        CHARACTER(len=*), PARAMETER :: this_routine='ConstructInitialStarGraph'
        
        CALL TISET('InitStarGraph',iSubInitStar)
        
        CALL IAZZERO(nStore,6)
!Having exFlag=2 means that only double excitations are generated
        exFlag=2
        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
        ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
        CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag)
        CALL IAZZERO(nExcit,nExcitMemLen)
        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
!First run through all excitations, so we do not store more excitations than necessary
        
        excitcount=0
        do while(.true.)
            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,iExcit,0,nStore,exFlag)
            IF(nJ(1).eq.0) exit
            CALL CalcRho2(FDet,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,iExcit,ECore)
            IF(rh.agt.0.D0) excitcount=excitcount+1
        enddo

        DEALLOCATE(nExcit)
        CALL IAZZERO(nStore,6)
        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
        ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
        CALL IAZZERO(nExcit,nExcitMemLen)
        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
        iMaxExcit=excitcount
!NDets is equal to the number of excitations plus 1 (for root)
        NDets=excitcount+1

        WRITE(6,"(A,I10)") "Total number of determinants in GraphMorph is ", NDets
        
        IF(TDistrib) THEN
!If TDistribs is true, then we are recording the change in distributions of the graphs over the iterations
            ALLOCATE(Distribs(NEl,Iters+1),stat=ierr)
            CALL LogMemAlloc('Distribs',NEl*(Iters+1),4,this_routine,DistribsTag)
            CALL IAZZERO(Distribs,NEl*(Iters+1))
        ENDIF
        IF(TNoSameExcit.or.TOneExcitConn) THEN
!We want to store the excitation levels of the determinants in the graph
            ALLOCATE(GraphExcitLevel(NDets),stat=ierr)
            CALL LogMemAlloc('GraphExcitLevel',NDets,4,this_routine,GraphExcitLevelTag)
!This will automatically set the first element (HF) to zero
            CALL IAZZERO(GraphExcitLevel,NDets)
        ENDIF

!Allocate memory needed for calculations
        ALLOCATE(GraphDets(NDets,NEl),stat=ierr)
        CALL LogMemAlloc('GraphDets',NDets*NEl,4,this_routine,GraphDetsTag)
        CALL IAZZERO(GraphDets,NDets*NEl)
        ALLOCATE(GraphRhoMat(NDets,NDets),stat=ierr)
        CALL LogMemAlloc('GraphRhoMat',NDets*NDets,8*HElementSize,this_routine,GraphRhoMatTag)
        CALL AZZERO(GraphRhoMat,NDets*NDets*HElementSize)
        ALLOCATE(HamElems(NDets),stat=ierr)
        CALL LogMemAlloc('HamElems',NDets,8*HElementSize,this_routine,HamElemsTag)
        CALL AZZERO(HamElems,NDets*HElementSize)

        HamElems(1)=Hii
        GraphRhoMat(1,1)=rhii
        do i=1,NEl
            GraphDets(1,i)=FDet(i)
        enddo
        MeanExcit=0.D0
        
        i=1
!Run through excitations
        do while(.true.)
            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,iExcit,0,nStore,exFlag)
            IF(nJ(1).eq.0) EXIT

!Since we already know that the excitations are double excitations of FDet, we can put those connections in seperatly (will be quicker)
            CALL CalcRho2(FDet,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,iExcit,ECore)

            IF(.not.(rh.agt.0.D0)) CYCLE
            i=i+1

            GraphRhoMat(1,i)=rh
            GraphRhoMat(i,1)=rh
            
!Store Path
            do j=1,NEl
                GraphDets(i,j)=nJ(j)
            enddo

!Find hamiltonian element coupling to FDet, and diagonal element of excitation
            HamElems(i)=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,iExcit,ECore)
            MeanExcit=MeanExcit+iExcit
            CALL CalcRho2(nJ,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,GraphRhoMat(i,i),nTay,0,ECore)

            IF(TDistrib) THEN
                Distribs(iExcit,1)=Distribs(iExcit,1)+1
            ENDIF
            IF(TNoSameExcit.or.TOneExcitConn) GraphExcitLevel(i)=iExcit

!Cycle through all excitations already generated to determine coupling to them
            IF((.not.TNoSameExcit).and.(.not.TOneExcitConn)) THEN
!If TNoSameExcit is on, then we are ignoring these crosslinks - the normal star graph should result
                do j=2,(i-1)
                    CALL CalcRho2(GraphDets(j,:),nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,-1,ECore)
                    GraphRhoMat(i,j)=rh
                    GraphRhoMat(j,i)=rh
                enddo
            ENDIF

        enddo

        MeanExcit=MeanExcit/(i-1)

        IF((ABS(MeanExcit-2.D0)).gt.1.D-08) THEN
            WRITE(6,*) "Error generating initial star graph"
            STOP 'Error generating initial star graph'
        ENDIF

        DEALLOCATE(nExcit)
        CALL LogMemDealloc(this_routine,nExcitTag)

        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in ConstructInitialStarGraph'
        CALL TIHALT('InitStarGraph',iSubInitStar)

    END SUBROUTINE ConstructInitialStarGraph

    SUBROUTINE ConstructInitialGraph()
        USE System , only : G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,Arr
        USE Calc , only : i_P,RhoEps,G_VMC_Pi,G_VMC_Seed
        USE Integrals , only : fck,nMax,nMsh,UMat,nTay
        USE Determinants , only : GetHElement2
        IMPLICIT NONE
        LOGICAL :: Attached,connected
        INTEGER :: ierr,nStore(6),nExcitMemLen,iMaxExcit,nJ(NEl),nExcitTag,iExcit
        INTEGER :: iPathTag,XijTag,RhoiiTag,RhoijTag,HijsTag,iSubInit,i,j,diff,k
        INTEGER :: iGetExcitLevel,IC
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
!        OldImport=G_VMC_Pi
!        G_VMC_Pi=1.D0
        
        IF(TDistrib) THEN
!If TDistribs is true, then we are recording the change in distributions of the graphs over the iterations
            ALLOCATE(Distribs(NEl,Iters+1),stat=ierr)
            CALL LogMemAlloc('Distribs',NEl*(Iters+1),4,this_routine,DistribsTag)
            CALL IAZZERO(Distribs,NEl*(Iters+1))
        ENDIF

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

        IF(TMoveDets) THEN
            connected=.false.
            k=1
            do while (.not.connected)
!Generate the initial graph...(Can speed this up as do not need to know prob, and are recalculating the rho matrix)
!Arr sent instead of NMax since getting values straight from modules (not passed through)
                CALL Fmcpr4d2GenGraph(FDet,NEl,Beta,i_P,iPath,NDets,Xij,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,Alat,UMat,nTay,RhoEps,Rhoii,Rhoij,ECore,Seed,Hijs,nExcit,0,ExcitGen)
!Need to check that all determinants in the graph are connected (if we are using movedets)
                do i=1,NDets
                    connected=.false.
                    do j=1,NDets
                        IF(i.eq.j) CYCLE
                        IF(Rhoij(i-1,j-1).agt.0.D0) THEN
                            connected=.true.
                            EXIT
                        ENDIF
                    enddo
                    IF(.not.connected) EXIT
                enddo
                IF(.not.connected) THEN
                    IF(k.gt.15) THEN
                        WRITE(6,*) "A fully connected initial graph could not be created - should now recreate the graph from scratch for one iteration..."
!                        STOP "A fully connected initial graph could not be created - exiting..."
                        EXIT
                    ENDIF
                    k=k+1
                    WRITE(6,"(A,I7,A)") "Determinant number ",i," not connected to rest of initial graph - trying again to create initial graph with new seed."
!change seed and repeat
                    seed=G_VMC_Seed-k
!Resetup arguments and excitation generator for fmcpr4d2gengraph
                    DEALLOCATE(nExcit)
                    CALL LogMemDealloc(this_routine,nExcitTag)
                    CALL IAZZERO(nStore,6)
                    CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
                    ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
                    CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag)
                    CALL IAZZERO(nExcit,nExcitMemLen)
                    CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,3)
                    CALL GenRandSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,nExcit,nJ,Seed,iExcit,0,UMat,nMax,PGen)
                    CALL IAZZERO(iPath,NEl*(NDets+1))
                    CALL AZZERO(Xij,NDets*NDets)
                    CALL AZZERO(Rhoii,NDets+1)
                    CALL AZZERO(Rhoij,(NDets+1)*(NDets+1)*HElementSize)
                    CALL AZZERO(Hijs,(NDets+1)*HElementSize)
                    do i=1,NEl
                        iPath(i,0)=FDet(i)
                        iPath(i,NDets)=FDet(i)
                    enddo
                    ExcitGen(:)=0
                ENDIF
            enddo
!            IF(.not.connected) THEN
!!Create a small artificial connection between disconnected determinants and root
!                do i=1,NDets
!                    connected=.false.
!                    do while(.not.connected)
!                        do j=1,NDets
!                            IF(i.eq.j) CYCLE
!                            IF(Rhoij(i-1,j-1).agt.0.D0) THEN
!                                connected=.true.
!                                EXIT
!                            ENDIF
!                        enddo
!                        IF(.not.connected) THEN
!                            Rhoij(i-1,1)=HElement(1.D-20)
!                            Rhoij(1,i-1)=HElement(1.D-20)
!                        ENDIF
!                    enddo
!                enddo
!            ENDIF
        ELSE
!Generate the initial graph...(Can speed this up as do not need to know prob, and are recalculating the rho matrix)
!Arr sent instead of NMax since getting values straight from modules (not passed through)
            CALL Fmcpr4d2GenGraph(FDet,NEl,Beta,i_P,iPath,NDets,Xij,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,Alat,UMat,nTay,RhoEps,Rhoii,Rhoij,ECore,Seed,Hijs,nExcit,0,ExcitGen)
        ENDIF

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
!                WRITE(6,"(E8.5)",advance='no') Rhoij(i-1,j-1)
                GraphRhoMat(i,j)=Rhoij(i-1,j-1)
            enddo
!            WRITE(6,*) ""
        enddo

!        do j=1,NDets
!            Attached=.false.
!            do i=1,NDets
!                IF(i.eq.j) CYCLE
!                IF((GraphRhoMat(i,j)%v).ne.0.D0) THEN
!                    Attached=.true.
!                    EXIT
!                ENDIF
!            enddo
!            IF(.not.attached) STOP 'determinant is not attached!'
!        enddo

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

!Find out manually the mean distance from the HF determinant...
        MeanExcit=0.D0
        do i=2,NDets
            IC=iGetExcitLevel(FDet(:),GraphDets(i,:),NEl)
            MeanExcit=MeanExcit+IC
            IF(TDistrib) THEN
                Distribs(IC,1)=Distribs(IC,1)+1
            ENDIF
        enddo
        MeanExcit=MeanExcit/(NDets-1.D0)

!Return G_VMC_Pi to original value (Just in case it wants to be used later)
!        G_VMC_Pi=OldImport
        
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

!This routine move determinants already in the graph, to excitations stocastically
!TO DO: since this only moves a few determinants, much of the excitation space is still the same
!this means that only a few of the excitations would need to be reproduced. However, bookkeeping
!may well be a nightmare
    SUBROUTINE MoveDetsGraph()
        USE System , only : G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,Arr
        USE Calc , only : i_P,RhoEps
        USE Integrals , only : fck,nMax,nMsh,UMat,nTay
        USE Determinants , only : GetHElement2
        IMPLICIT NONE
        INTEGER :: iSubMove,i,j,k,l,Success,Failure,iGetExcitLevel,IC,Excitation,IC1,IC2
        INTEGER , ALLOCATABLE :: MoveDetsFromPaths(:,:)
        INTEGER :: MoveDetsFromPathsTag=0
        INTEGER :: AttemptDet(NEl),IndexofDetsFrom(NoMoveDets),ierr,Tries,NoVerts
        LOGICAL :: Remove,Attach,SameDet
        TYPE(HElement) :: rh
        REAL*8 :: r,Ran2
        CHARACTER(len=*), PARAMETER :: this_routine='MoveDets'

        CALL TISET('MoveDets',iSubMove)

!Allocate space for the determinants to be moved, and to move to
        ALLOCATE(MoveDetsFromPaths(NoMoveDets,NEl),stat=ierr)
        CALL LogMemAlloc('MoveDetsFromPaths',NoMoveDets*NEl,4,this_routine,MoveDetsFromPathsTag)
        CALL IAZZERO(MoveDetsFromPaths,NoMoveDets*NEl)
        CALL IAZZERO(IndexofDetsFrom,NoMoveDets)

!Prepare for the change in MeanExcits
        MeanExcit=MeanExcit*(NDets-1)
        IF(TDistrib) THEN
!If we are recording the distributions, then we want to simply have the same distribution for the previous
!iteration, and then subtract the determinants that we remove
            do i=1,NEl
                Distribs(i,Iteration+1)=Distribs(i,Iteration)
            enddo
        ENDIF

!First need to pick NoMoveDets determinants stochastically to be moved from the original graph
        Tries=NoMoveDets*10000*NDets
        k=0
        Success=0
        Failure=0
!No vertices chosen to be removed yet
        NoVerts=0
        do while ((NoVerts.lt.NoMoveDets).and.(k.le.Tries))
            
            k=k+1
            r=RAN2(Seed)
!Set i=1, since we do not want to choose the first determinant of the original graph
            i=1

!Look through the normalised inverse eigenvector*eigenvalue (not squared now)
            do while ((r.gt.0.D0).and.(i.lt.NDets))
                i=i+1
                r=r-(Eigenvector(i)%v)
            enddo
            IF(r.gt.0.D0) STOP 'Problem in counting in MoveDets'
            
            IF(GraphDets(i,1).eq.0) THEN
!Graph previously been selected and accepted
                Remove=.false.
            ELSE
                Remove=.true.
                do j=1,NEl
                    AttemptDet(j)=GraphDets(i,j)
                enddo
            ENDIF

!Need to test whether graph is allowed to be removed - if not, they cycle around for another attempt at finding a valid determinant
!            Remove=.true.
!            DO j=1,NoVerts
!!Chosen determinant cannot already have been chosen, and must ensure that the other determinants in the graph are still attached.
!!If determinant has been chosen before, then its GraphDets-> 0 , so need another test
!                IF(AttemptDet(1).eq.0) THEN
!                    Remove=.false.
!                    EXIT
!                ENDIF
!                IF(SameDet(AttemptDet(:),MoveDetsFromPaths(j,:),NEl)) THEN
!                    Remove=.false.
!                    EXIT
!                ENDIF
!            ENDDO

            IF(Remove) THEN
!Need to look through the rho matrix to check that there is still a connection for the other determinants in the graph
                do j=1,NDets
!If we are at a diagonal element, ignore as it is not a connection
                    IF(j.eq.i) CYCLE
                    IF(CopyRhoMat(i,j).gt.0.D0) THEN
!Find that the determinant we are trying to move (i) is already connected to a different determinant (j)
!See if there are any other connections to determinant j
                        Remove=.false.
                        do l=1,NDets
!Ignore its diagonal element, and the connection to i
                            IF((l.eq.j).or.(l.eq.i)) CYCLE
                            IF(ABS(CopyRhoMat(l,j)).gt.0.D0) THEN
!We have found that j is connected to a different determinant, so it is ok to remove i
                                Remove=.true.
                                EXIT
                            ENDIF
                        enddo
                    ENDIF
!The determinant we want to remove would leave a disconnected graph - we cannot use it
                    IF(.not.Remove) EXIT
                enddo
            ENDIF

            IF(Remove) THEN
!We are able to successfully remove the chosen determinant
                Success=Success+1
                NoVerts=NoVerts+1
                do j=1,NEl
                    MoveDetsFromPaths(NoVerts,j)=AttemptDet(j)
                    GraphDets(i,j)=0
                enddo
!Store index of determinant to move
                IndexofDetsFrom(NoVerts)=i
!Remove Determinant from rho matrix
                do j=1,NDets
                    CopyRhoMat(j,i)=0.D0
                    CopyRhoMat(i,j)=0.D0
                enddo
                HamElems(i)=HElement(0.D0) 
!Calculate the change to the MeanExcits value...
                IC=iGetExcitLevel(FDet(:),AttemptDet(:),NEl)
                MeanExcit=MeanExcit-IC       
                IF(TDistrib) THEN
                    Distribs(IC,Iteration+1)=Distribs(IC,Iteration+1)-1
                ENDIF
            ELSE
                Failure=Failure+1
            ENDIF

!Cycle and try to find another determinant to remove
        enddo

!MC has failed - print out debugging info
        IF(k.gt.Tries) THEN
            WRITE(6,*) 'Error in trying to pick a determinant to move'
            WRITE(6,*) NoVerts, " removed from previous graph"
            WRITE(6,*) "Vertices which have been removed are: "
            do j=1,NoVerts
                WRITE(6,*) IndexofDetsFrom(j)
            enddo
            WRITE(6,*) Tries, " attempted moves made"
            WRITE(6,*) "Determinants in graph selected with probability: "
            do j=2,NDets
                WRITE(6,*) Eigenvector(j)
            enddo
            CALL FLUSH(6)
            STOP 'Error in trying to pick a determinant to move'
        ELSE
!            WRITE(6,*) "Determinants which have been removed are: "
!            do j=1,NoMoveDets
!                WRITE(6,*) IndexofDetsFrom(j)
!            enddo
        ENDIF
!        WRITE(6,*) "**************"
!        do j=1,NoMoveDets
!            WRITE(6,*) MoveDetsFromPaths(j,:)
!        enddo
!        CALL FLUSH(6)

        
!Now need to find a determinant from the excitations space to attach
        NoVerts=0
!NoVerts is the number of vertices already attached
        do while ((NoVerts.lt.NoMoveDets).and.(k.le.Tries))
            
            k=k+1
            r=RAN2(Seed)
            i=0

!Look through the normalised ExcitsVector array 
            do while ((r.gt.0.D0).and.(i.lt.TotExcits))
                i=i+1
                r=r-((ExcitsVector(i)%v)**2)
            enddo
            IF(r.gt.0.D0) STOP 'Problem in counting in MoveDets 2'
            
            do j=1,NEl
                AttemptDet(j)=ExcitsDets(i,j)
            enddo

!First we want to check whether the determinant we have selected has already been selected and is already in the graph
            Attach=.true.
            do j=1,NoVerts
                IF(SameDet(AttemptDet(:),GraphDets(IndexofDetsFrom(j),:),NEl)) THEN 
                    Attach=.false.
                    EXIT
                ENDIF
            enddo
!Allow the determinant to be one we have already got rid of; this is possible, since a determinant can be specified more than once
            IF(.not.Attach) THEN
                do l=(1+NoVerts),NoMoveDets
                    IF(j.eq.IndexofDetsFrom(l)) THEN
                        Attach=.true.
                        EXIT
                    ENDIF
                enddo
            ENDIF

            IF(Attach) THEN
!Before allowing it to be attached, we must check if the determinant chosen is attached to the graph at all
                Attach=.false.
                IF(TNoSameExcit) IC1=iGetExcitLevel(FDet(:),AttemptDet(:),NEl)

         loop2: do j=1,NDets
                    do l=(1+NoVerts),NoMoveDets
!Check we are not trying to attach to a determinant which we have chosen to remove! - (Unless it has already been replaced by a new det)
                        IF(j.eq.IndexofDetsFrom(l)) CYCLE loop2
                    enddo

                    IF(TNoSameExcit) THEN
!IC2 information for vertices already added to graph stored in GraphExcitLevel
                        IC2=GraphExcitLevel(j)
!                        IC2=iGetExcitLevel(FDet(:),GraphDets(j,:),NEl)
!If TNoSameDet is on, then do not allow connections between determinants from the same excitation level
                        IF(IC2.ne.IC1) THEN
                            IC=iGetExcitLevel(AttemptDet(:),GraphDets(j,:),NEl)
                            CALL CalcRho2(AttemptDet(:),GraphDets(j,:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
                            IF(rh.agt.0.D0) THEN
                                Attach=.true.
                                CopyRhoMat(IndexofDetsFrom(NoVerts+1),j)=rh%v
                                CopyRhoMat(j,IndexofDetsFrom(NoVerts+1))=rh%v
                            ENDIF
                        ENDIF
                    ELSE
!Allow all possible connections
                        IC=iGetExcitLevel(AttemptDet(:),GraphDets(j,:),NEl)
                        IF(IC.eq.0) THEN
!We have selected a determinant which is already in the graph
                            Attach=.false.
                            EXIT loop2
                        ENDIF
                        CALL CalcRho2(AttemptDet(:),GraphDets(j,:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
                        IF(rh.agt.0.D0) THEN
                            Attach=.true.
                            CopyRhoMat(IndexofDetsFrom(NoVerts+1),j)=rh%v
                            CopyRhoMat(j,IndexofDetsFrom(NoVerts+1))=rh%v
                        ENDIF
                    ENDIF

                enddo loop2
            ENDIF

!Attach the determinant
            IF(Attach) THEN
                Success=Success+1
                NoVerts=NoVerts+1
                do j=1,NEl
                    GraphDets(IndexofDetsFrom(NoVerts),j)=AttemptDet(j)
                enddo
!Find Diagonal rho matrix element, and hamiltonian element
                HamElems(IndexofDetsFrom(NoVerts))=GetHElement2(FDet,AttemptDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,IC1,ECore)
                CALL CalcRho2(AttemptDet(:),AttemptDet(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,0,ECore)
                CopyRhoMat(IndexofDetsFrom(NoVerts),IndexofDetsFrom(NoVerts))=rh%v
                IF(.not.TNoSameExcit) THEN
                    IC1=iGetExcitLevel(FDet(:),AttemptDet(:),NEl)
                ENDIF
!Calculate the change to the MeanExcits value...
                MeanExcit=MeanExcit+IC1          
!Store the excitation level of the new determinant we are adding
                IF(TNoSameExcit) GraphExcitLevel(IndexofDetsFrom(NoVerts+1))=IC1
                IF(TDistrib) THEN
                    Distribs(IC1,Iteration+1)=Distribs(IC1,Iteration+1)+1
                ENDIF

            ELSE
                Failure=Failure+1
            ENDIF
!Try to attach another determinant
        enddo
        
        IF(k.gt.Tries) THEN
            WRITE(6,*) 'Error in trying to pick a determinant to create'
            STOP 'Error in trying to pick a determinant to create'
        ENDIF

!Move RhoMatrix into the full one
        ALLOCATE(GraphRhoMat(NDets,NDets),stat=ierr)
        CALL LogMemAlloc('GraphRhoMat',NDets*NDets,8*HElementSize,this_routine,GraphRhoMatTag)
        do i=1,NDets
            do j=1,i
                GraphRhoMat(i,j)=HElement(CopyRhoMat(i,j))
                GraphRhoMat(j,i)=HElement(CopyRhoMat(j,i))
            enddo
        enddo

!Find the final value for MeanExcits, and the success ratio of the algorithm
        MeanExcit=MeanExcit/(NDets-1)
        SucRat=(Success+0.D0)/(Success+Failure+0.D0)
        Orig_Graph=NoMoveDets/(NDets+0.D0)

        DEALLOCATE(CopyRhoMat)
        CALL LogMemDealloc(this_routine,CopyRhoMatTag)
        DEALLOCATE(MoveDetsFromPaths)
        CALL LogMemDealloc(this_routine,MoveDetsFromPathsTag)
!In the future, these may not have to be completly changed, only modified slightly
        DEALLOCATE(ExcitsVector)
        CALL LogMemDealloc(this_routine,ExcitsVectorTag)
        DEALLOCATE(ExcitsDets)
        CALL LogMemDealloc(this_routine,ExcitsDetsTag)
        
        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in PickNewDets'
        CALL TIHALT('MoveDets',iSubMove)

    END SUBROUTINE MoveDetsGraph

!This routine stocastically picks NDets new determinants stochastically from the vector obtained.
    SUBROUTINE PickNewDets()
        USE System , only : G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,Arr
        USE Calc , only : i_P,RhoEps
        USE Integrals , only : fck,nMax,nMsh,UMat,nTay
        USE Determinants , only : GetHElement2
        IMPLICIT NONE
        INTEGER :: ierr,i,j,k,NoVerts,AttemptDet(NEl),GrowGraphTag,IC,dist
        INTEGER :: Success,Failure,OrigDets,ExcitDets,Tries,iGetExcitLevel,IC1,IC2
        INTEGER , ALLOCATABLE :: GrowGraph(:,:)
        REAL*8 :: r,Ran2
        LOGICAL :: Attach,OriginalPicked,SameDet
        TYPE(HElement) :: rh
        CHARACTER(len=*), PARAMETER :: this_routine='PickNewDets'

        GrowGraphTag=0
        IF(TNoSameExcit.or.TOneExcitConn) CALL IAZZERO(GraphExcitLevel,NDets)

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
                r=r-((Eigenvector(i)%v)**2)
            enddo

!If still not found, then look in Vector of excitations
            IF(r.gt.0.D0) THEN
                i=0
                do while ((r.ge.0.D0).and.(i.lt.TotExcits))
                    i=i+1
                    r=r-((ExcitsVector(i)%v)**2)
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
!IC1 is the excitation level of the attempted determinant to attach
            IC1=iGetExcitLevel(FDet,AttemptDet,NEl)
            DO i=1,NoVerts

!Check the number of excitations away from each other vertex in graph - or instead, just check if determinant already in graph
                IF(SameDet(AttemptDet(:),GrowGraph(i,:),NEl)) THEN
!Determinant is already in the graph - exit loop - determinant not valid
                    Attach=.false.
                    EXIT
                ENDIF

!Determine connectivity to other determinants in the graph, if no connection yet found
                IF(.not.Attach) THEN
                    IF(TNoSameExcit) THEN
!excitation level of determinants already in the graph stored in GraphExcitLevel
!                        IC2=iGetExcitLevel(FDet(:),GrowGraph(i,:),NEl)
                        IC2=GraphExcitLevel(i)
!Don't allow connections to the same excitation level
                        IF(IC1.ne.IC2) THEN
                            IC=iGetExcitLevel(GrowGraph(i,:),AttemptDet(:),NEl)
                            CALL CalcRho2(AttemptDet(:),GrowGraph(i,:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
                            IF(rh.agt.RhoEps) Attach=.true. 
                        ENDIF

                    ELSEIF(TOneExcitConn) THEN
!Only allow connections between excitation levels which differ by one
                        IC2=GraphExcitLevel(i)
                        IF((ABS(IC1-IC2)).eq.1) THEN
                            IC=iGetExcitLevel(GrowGraph(i,:),AttemptDet(:),NEl)
                            CALL CalcRho2(AttemptDet(:),GrowGraph(i,:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
                            IF(rh.agt.RhoEps) Attach=.true. 
                        ENDIF

                    ELSE
!Allow all connections
                        CALL CalcRho2(AttemptDet(:),GrowGraph(i,:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,-1,ECore)
                        IF(rh.agt.RhoEps) Attach=.true. 
                    ENDIF
                ENDIF

            ENDDO

            IF(Attach) THEN
!Attach new determinant to the list as it passes tests
!First check on attachment to HF, since wants to be stored in HamElems
!                IF(IC1.gt.2) WRITE(6,*) "Higher excitation attached ", IC1
     
                MeanExcit=MeanExcit+IC1
                IF(TDistrib) THEN
                    Distribs(IC1,Iteration+1)=Distribs(IC1,Iteration+1)+1
                ENDIF

!Only double excitations (or single? <- include for completness) of HF contribute
                IF(IC1.eq.2) THEN!.or.(IC1.eq.1)) THEN
                    HamElems(NoVerts+1)=GetHElement2(FDet,AttemptDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,IC1,ECore)
                    CALL CalcRho2(FDet,AttemptDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC1,ECore)
!Add rhoij contribution. No real need to set a rho epsilon here? We are diagonalising anyway...
                    GraphRhoMat(1,NoVerts+1)=rh
                    GraphRhoMat(NoVerts+1,1)=rh

                ELSE
!Added determinant not connected to root
                    HamElems(NoVerts+1)=HElement(0.D0)
                    GraphRhoMat(1,NoVerts+1)=HElement(0.D0)
                    GraphRhoMat(NoVerts+1,1)=HElement(0.D0)
                ENDIF

                IF(TNoSameExcit) THEN
!Don't allow connections to the same excitation level
                    do i=2,NoVerts
!Excitation level information stored in GraphExcitLevel
!                        IC2=iGetExcitLevel(FDet(:),GrowGraph(i,:),NEl)
                        IC2=GraphExcitLevel(i)
                        IF(IC1.ne.IC2) THEN
                            IC=iGetExcitLevel(AttemptDet(:),GrowGraph(i,:),NEl)
                            CALL CalcRho2(GrowGraph(i,:),AttemptDet(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
                            GraphRhoMat(i,NoVerts+1)=rh
                            GraphRhoMat(NoVerts+1,i)=rh
                        ENDIF
                    enddo
!Store the excitation level of the attached determinant
                    GraphExcitLevel(NoVerts+1)=IC1

                ELSEIF(TOneExcitConn) THEN
!Only allow connections between excitations which differ by one excitation level from HF
                    do i=2,NoVerts
!Excitation level information stored in GraphExcitLevel
                        IC2=GraphExcitLevel(i)
                        IF((ABS(IC1-IC2)).eq.1) THEN
                            IC=iGetExcitLevel(AttemptDet(:),GrowGraph(i,:),NEl)
                            CALL CalcRho2(GrowGraph(i,:),AttemptDet(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
                            GraphRhoMat(i,NoVerts+1)=rh
                            GraphRhoMat(NoVerts+1,i)=rh
                        ENDIF
                    enddo
!Store the excitation level of the attached determinant
                    GraphExcitLevel(NoVerts+1)=IC1

                ELSE

!Run through rest of excitations in the graph testing contributions and adding to rho matrix
                    do i=2,NoVerts
                        CALL CalcRho2(GrowGraph(i,:),AttemptDet(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,-1,ECore)
                        GraphRhoMat(i,NoVerts+1)=rh
                        GraphRhoMat(NoVerts+1,i)=rh
                    enddo

                ENDIF

!Include diagonal rho matrix element for chosen determinant
                CALL CalcRho2(AttemptDet,AttemptDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,0,ECore)
                GraphRhoMat(NoVerts+1,NoVerts+1)=rh

!Add determinant to growing graph
                do i=1,NEl
                    GrowGraph(NoVerts+1,i)=AttemptDet(i)
                enddo
                NoVerts=NoVerts+1
!                    WRITE(6,"A,I5") "Vertex Added - ",NoVerts
!                    CALL FLUSH(6)

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

!This is used to create normalised vectors for the MoveDets stochastic graph-growing algorithm
    SUBROUTINE NormaliseVectorSep()
        IMPLICIT NONE
        TYPE(HElement) :: Norm1,Norm2
        REAL*8 :: RootofNum
        INTEGER :: iSubNorm,i
        
        CALL TISET('NormVecSep',iSubNorm)
        
!Setup a normalised Inverse vector of determinants in the graph, so they can be chosen stochastically
!Since we no longer need the largest eigenvector, we can multiply its elements by its eigenvalue, then
!use it as the inverse vector
        do i=1,NDets
            Eigenvector(i)=Eigenvector(i)*HElement(Eigenvalue)
        enddo

        Norm1=HElement(0.D0)

        do i=2,NDets
!            IF(ABS(Eigenvector(i)%v).lt.1.D-17) THEN
!                WRITE(6,*) "For determinant ", i, ", connection is ", Eigenvector(i)
!                STOP 'Numerical errors likely to arise due to such small connection'
!            ENDIF
!Normalised as the reciprocal of the fourth root of the element. If root is changed, then numerical stability means that only one determinant is ever picked...
            Eigenvector(i)=HElement(1.D0/RootofNum(ABS(Eigenvector(i)%v),5))
!            WRITE(6,*) Eigenvector(i)
            Norm1=Norm1+Eigenvector(i)
        enddo

        do i=2,NDets
            Eigenvector(i)=Eigenvector(i)/Norm1
!            WRITE(6,*) Eigenvector(i)
        enddo

!The excitations need to be normalised separatly, according to their amplitude squared

        Norm2=HElement(0.D0)
        do i=1,TotExcits
            Norm2=Norm2+(ExcitsVector(i)*ExcitsVector(i))
        enddo
        Norm2=HElement(SQRT(Norm2%v))
        do i=1,TotExcits
            ExcitsVector(i)=ExcitsVector(i)/Norm2
        enddo

        WRITE(6,"(A,F12.5)") "Normalisation constant for the determinants in the graph is ",Norm1%v

        PStay=NoMoveDets/NDets

        CALL TIHALT('NormVecSep',iSubNorm)

    END SUBROUTINE NormaliseVectorSep

!This is used to normalise the vector of current and excited determinants from which we are going to pick the next graph.
!The vector is spread between the arrays for the original eigenvector (which needs to be multiplied by corresponding eigenvalue)
!and the ExcitsVector. The first element of the eigenvector should not be included in the normalisation, as
!it cannot be picked - the HF is always in each graph.
    SUBROUTINE NormaliseVector()
        IMPLICIT NONE
        INTEGER :: i
        REAL*8 :: Stay,Move
        TYPE(HElement) :: Norm,Norm1,Norm2

!The bias towards the determinants already in the graph is given by the largest eigenvector, multiplied by its eigenvalue.
!Since we no longer need the largest eigenvector, we can multiply its elements by its eigenvalue
        do i=2,NDets
            Eigenvector(i)=Eigenvector(i)*HElement(Eigenvalue)
        enddo

!An addititional bias against the determinants already in the graph can be designed.
        IF(TBiasing) THEN
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
!        CALL FLUSH(6)
!        WRITE(6,*) "Total Probability: ", Stay+Move
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
        REAL*8 :: Prob
        INTEGER , SAVE :: iSubConns
        INTEGER :: attempts,NoExcitsCurr,Noatt
        INTEGER :: ierr,i,j,DetCurr(NEl),nJ(NEl),nStore(6),iMaxExcit,nExcitMemLen
        INTEGER :: nExcitTag,iExcit,ExcitCurr,dist,iGetExcitLevel,IC,exFlag
        INTEGER , ALLOCATABLE :: nExcit(:)
        CHARACTER(len=*), PARAMETER :: this_routine='FindConnections'

        CALL TISET('FindConns',iSubConns)
        
        IF(TSinglesExcitSpace) THEN
!Only single excitations of the determinants in the graph are created
            exFlag=1
        ELSE
            exFlag=3
        ENDIF

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
            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
            ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
            CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag)
!            CALL IAZZERO(nExcit,nExcitMemLen)
            nExcit(1)=0
            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)

            IF(TMCExcitSpace) THEN
!Excitations are picked stocastically

                Attempts=0
                NoExcitsCurr=0

                do while(NoExcitsCurr.lt.NoMCExcits)

                    IF(Attempts.gt.(NoMCExcits*100)) THEN
                        WRITE(6,*) "Unable to find enough determinants attached to graph determinant ",DetCurr
                        WRITE(6,*) "Attempts: ", Attempts
                        STOP "Unable to find enough determinants attached to graph determinant "
                    ENDIF

                    CALL GenRandSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,nExcit,nJ,Seed,Noatt,0,UMat,Arr,Prob)
                    Attempts=Attempts+1
!Divide the coupling to a determinant by the probability of generating it, to account for possibility of bias in generation

                    IF(TMaxExcit) THEN

!A maximum excitation level is set - don't allow connection if its too high in excitation space
                        IC=iGetExcitLevel(FDet,nJ,NEl)
                        IF(IC.le.iMaxExcitLevel) THEN
!Throw the random excitation, and don't count it if it is above the excitation level threshold
                            ExcitCurr=ExcitCurr+1
                            NoExcitsCurr=NoExcitsCurr+1
                            CALL CalcRho2(DetCurr,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,-1,ECore)
!Store path of determinant in ExcitsDets, and the rho elements in ConnectionsToExcits
                            rh=rh/(HElement(Prob))
                            ConnectionsToExcits(ExcitCurr)=rh
                            do j=1,NEl
                                ExcitsDets(ExcitCurr,j)=nJ(j)
                            enddo

                        ENDIF

                    ELSE
!No restriction on excitation level 
        
                        ExcitCurr=ExcitCurr+1
                        NoExcitsCurr=NoExcitsCurr+1
                        CALL CalcRho2(DetCurr,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,-1,ECore)
!Store path of determinant in ExcitsDets, and the rho elements in ConnectionsToExcits
                        rh=rh/(HElement(Prob))
                        ConnectionsToExcits(ExcitCurr)=rh
                        do j=1,NEl
                            ExcitsDets(ExcitCurr,j)=nJ(j)
                        enddo

                    ENDIF

                enddo

            ELSE

                IF(i.eq.1) THEN
                    IF(iMaxExcit.ne.NoExcits(i)) STOP 'Error in counting in FindConnections'
                ELSE
                    IF(iMaxExcit.ne.(NoExcits(i)-NoExcits(i-1))) STOP 'Error in counting in FindConnections'
                ENDIF

!Cycle through all excitations of each determinant
            lp: do while(.true.)
                    CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,iExcit,0,nStore,exFlag)
                    IF(nJ(1).eq.0) exit lp
                    ExcitCurr=ExcitCurr+1

                    IF(TMaxExcit) THEN
!A maximum excitation level for the space of accessible determinants is imposed
                        IC=iGetExcitLevel(FDet,nJ,NEl)
                        IF(IC.gt.iMaxExcitLevel) THEN
!If excitation is further away than we want, then let connection to it = 0
                            ConnectionsToExcits(ExcitCurr)=HElement(0.D0)
                        ELSE
                            CALL CalcRho2(DetCurr,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,iExcit,ECore)
!Store path of determinant in ExcitsDets, and the rho elements in ConnectionsToExcits
                            ConnectionsToExcits(ExcitCurr)=rh
                            do j=1,NEl
                                ExcitsDets(ExcitCurr,j)=nJ(j)
                            enddo
                        ENDIF

                    ELSE
!No restriction on excitation level

                        CALL CalcRho2(DetCurr,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,iExcit,ECore)
!                    Dist=IGetExcitLevel(FDet,nJ,NEl)
!                    IF(Dist.gt.4) THEN
!                        WRITE(6,*) "Higher than double excitation found - attached to det:", DetCurr
!                        WRITE(6,*) "Determinant is: ", nJ
!                        WRITE(6,*) Dist," fold excitation"
!                    ENDIF

!Store path of determinant in ExcitsDets, and the rho elements in ConnectionsToExcits
                        ConnectionsToExcits(ExcitCurr)=rh
                        do j=1,NEl
                            ExcitsDets(ExcitCurr,j)=nJ(j)
                        enddo
                    ENDIF
                enddo lp

            ENDIF

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
        INTEGER :: nExcitTag,exFlag
        CHARACTER(len=*), PARAMETER :: this_routine='CountExcits'
        INTEGER , ALLOCATABLE :: nExcit(:)

        IF(TMCExcitSpace) THEN
!This searches for a given number of excitations stocastically per determinant - no need to count them
            
            TotExcits=NDets*NoMCExcits
            ALLOCATE(NoExcits(NDets),stat=ierr)
            CALL LogMemAlloc('NoExcits',NDets,4,this_routine,NoExcitsTag)
            CALL IAZZERO(NoExcits,NDets)

            do i=1,NDets
                IF(i.eq.1) THEN
                    NoExcits(i)=NoMCExcits
                ELSE
                    NoExcits(i)=NoExcits(i-1)+NoMCExcits
                ENDIF
            enddo

        ELSE
!Search through the entire space of single or single&double excitations of each determinant

            IF(TSinglesExcitSpace) THEN
!Only search through single excitations
                exFlag=1
            ELSE
                exFlag=3
            ENDIF
            
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
                CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
                ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
                CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag)
                CALL IAZZERO(nExcit,nExcitMemLen)
                CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
                
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

        ENDIF

!Perform check that all excitations accounted for
        IF(NoExcits(NDets).ne.TotExcits) THEN
            STOP 'Error in counting excits in CountExcits'
        ENDIF

        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in CountExcits'
    END SUBROUTINE CountExcits

!This routine uses a Lancsoz iterative diagonalisation technique to diagonalise the rho matrix.
    SUBROUTINE DiagGraphLanc()
        IMPLICIT NONE
        CHARACTER(len=*), PARAMETER :: this_routine='DiagGraphLanc'
        INTEGER , ALLOCATABLE :: Lab(:,:),NRow(:),ISCR(:),Index(:)
        TYPE(HElement) , ALLOCATABLE :: Mat(:,:),CK(:,:),CKN(:,:)
        REAL*8 , ALLOCATABLE :: A(:,:),V(:),AM(:),BM(:),T(:),WT(:),SCR(:)
        REAL*8 , ALLOCATABLE :: Work2(:),WH(:),V2(:,:),W(:)
        INTEGER :: LabTag,ATag,MatTag,NRowTag,VTag,WTTag
        INTEGER :: AMTag,BMTag,TTag,SCRTag,ISCRTag,IndexTag,WHTag
        INTEGER :: Work2Tag,V2Tag,WTag,CKTag,CKNTag
        INTEGER :: iSubLanc,ierr,LenMat,i,j,ICMax,RowElems
        REAL*8 :: B2L
        INTEGER :: NEval,NBlk,NKry,NCycle,NBlock,NKry1,LScr,LIScr

        CALL TISET('DiagGraphLanc',iSubLanc)
        
!Zero memory tags
        LabTag=0
        ATag=0
        MatTag=0
        NRowTag=0
        VTag=0
        WTTag=0
        AMTag=0
        BMTag=0
        TTag=0
        SCRTag=0
        ISCRTag=0
        IndexTag=0
        WHTag=0
        Work2Tag=0
        V2Tag=0
        WTag=0
        CKTag=0
        CKNTag=0

        IF(TMoveDets) THEN
!If using the MoveDets algorithm, need to keep the rho matrix for the previous graph
            ALLOCATE(CopyRhoMat(NDets,NDets),stat=ierr)
            CALL LogMemAlloc('CopyRhoMat',NDets**2,8,this_routine,CopyRhoMatTag)
            do i=1,NDets
                do j=i,NDets
                    CopyRhoMat(i,j)=GraphRhoMat(i,j)%v
                    CopyRhoMat(j,i)=GraphRhoMat(j,i)%v
                enddo
            enddo
        ENDIF

!Count the number of non-zero elements in the matrix
        LenMat=0
        ICMax=0
        do i=1,NDets
            RowElems=0
            do j=1,NDets
                IF(GraphRhoMat(i,j).agt.0.D0) THEN
                    LenMat=LenMat+1
                    RowElems=RowElems+1
                ENDIF
            enddo
            IF(RowElems.gt.ICMax) ICMax=RowElems
        enddo
        WRITE(6,"(A,I10,A,I10,A)") "Of ",NDets*NDets, " possible elements, only ", LenMat," are non-zero."

!This seems to be Alex's way of compressing the matrix, but not the way in the diagonaliser
!Allocate the rho matrix
!        ALLOCATE(Mat(LenMat),stat=ierr)
!        CALL LogMemAlloc('Mat',LenMat,8*HElementSize,this_routine,MatTag)
!        CALL AZZERO(Mat,LenMat*HElementSize)
!
!!Lab indicates the column that the 'i'th non-zero matrix element resides in the full matrix
!        ALLOCATE(Lab(LenMat),stat=ierr)
!        CALL LogMemAlloc('Lab',LenMat,4,this_routine,LabTag)
!        CALL IAZZERO(Lab,LenMat)

        ALLOCATE(Mat(NDets,ICMax),stat=ierr)
        CALL LogMemAlloc('Mat',NDets*ICMax,8*HElementSize,this_routine,MatTag)
        CALL AZZERO(Mat,NDets*ICMax*HElementSize)

!Lab now indicates the position of the non-zero element in the row specified
        ALLOCATE(Lab(NDets,ICMax),stat=ierr)
        CALL LogMemAlloc('Lab',NDets*ICMax,4,this_routine,LabTag)
        CALL IAZZERO(Lab,NDets*ICMax)

!NRow indicated the number of non-zero elements in the 'i'th row of the full matrix
        ALLOCATE(NRow(NDets),stat=ierr)
        CALL LogMemAlloc('NRow',NDets,4,this_routine,NRowTag)
        CALL IAZZERO(NRow,NDets)

!This compresses the rho matrix, so that only the non-zero elements are stored, in a suitable form for the Lanczos diagonaliser.
        CALL CompressMatrix(Mat,Lab,NRow,LenMat,ICMax)
        WRITE(6,"(A,I8)") "In compressed matrix, maximum number of non-zero elements in any one row is ", ICMax

!Set up parameters for the diagonaliser.
        B2L=1.D-13
        NBlk=4
        NKry=8
!NEval indicates the number of eigenvalues we want to calculate
        NEval=4
        NCycle=200
        NBlock=MIN(NEval,NBlk)
        NKry1=NKry+1
        LScr=MAX(NDets*NEval,8*NBlock*NKry)
        LIScr=6*NBlock*NKry

!Deallocate GraphRhoMat - no longer needed
        DEALLOCATE(GraphRhoMat)
        CALL LogMemDealloc(this_routine,GraphRhoMatTag)

!Allocate memory for diagonaliser
        ALLOCATE(A(NEval,NEval),stat=ierr)
        CALL LogMemAlloc('A',NEval*NEval,8,this_routine,ATag)
        CALL AZZERO(A,NEval*NEval)
        ALLOCATE(V(NDets*NBlock*NKry1),stat=ierr)
        CALL LogMemAlloc('V',NDets*NBlock*NKry1,8,this_routine,VTag)
        CALL AZZERO(V,NDets*NBlock*NKry1)
        ALLOCATE(AM(NBlock*NBlock*NKry1),stat=ierr)
        CALL LogMemAlloc('AM',NBlock*NBlock*NKry1,8,this_routine,AMTag)
        CALL AZZERO(AM,NBlock*NBlock*NKry1)
        ALLOCATE(BM(NBlock*NBlock*NKry),stat=ierr)
        CALL LogMemAlloc('BM',NBlock*NBlock*NKry,8,this_routine,BMTag)
        CALL AZZERO(BM,NBlock*NBlock*NKry)
        ALLOCATE(T(3*NBlock*NKry*NBlock*NKry),stat=ierr)
        CALL LogMemAlloc('T',NBlock*NKry*NBlock*NKry*3,8,this_routine,TTag)
        CALL AZZERO(T,NBlock*NKry*NBlock*NKry*3)
        ALLOCATE(WT(NBlock*NKry),stat=ierr)
        CALL LogMemAlloc('WT',NBlock*NKry,8,this_routine,WTTag)
        CALL AZZERO(WT,NBlock*NKry)
        ALLOCATE(SCR(LScr),stat=ierr)
        CALL LogMemAlloc('SCR',LScr,8,this_routine,SCRTag)
        CALL AZZERO(SCR,LScr)
        ALLOCATE(ISCR(LIScr),stat=ierr)
        CALL LogMemAlloc('ISCR',LIScr,4,this_routine,ISCRTag)
        CALL IAZZERO(ISCR,LIScr)
        ALLOCATE(Index(NEval),stat=ierr)
        CALL LogMemAlloc('Index',NEval,4,this_routine,IndexTag)
        CALL IAZZERO(Index,NEval)
        ALLOCATE(WH(NDets),stat=ierr)
        CALL LogMemAlloc('WH',NDets,8,this_routine,WHTag)
        CALL AZZERO(WH,NDets)
        ALLOCATE(Work2(3*NDets),stat=ierr)
        CALL LogMemAlloc('Work2',3*NDets,8,this_routine,Work2Tag)
        CALL AZZERO(Work2,3*NDets)
        ALLOCATE(V2(NDets,NEval),stat=ierr)
        CALL LogMemAlloc('V2',NDets*NEval,8,this_routine,V2Tag)
        CALL AZZERO(V2,NDets*NEval)

!W holds the eigenvalues 
        ALLOCATE(W(NEval),stat=ierr)
        CALL LogMemAlloc('W',NEval,8,this_routine,WTag)
        CALL AZZERO(W,NEval)

!CK holds the eigenvectors
        ALLOCATE(CK(NDets,NEval),stat=ierr)
        CALL LogMemAlloc('CK',NDets*NEval,8*HElementSize,this_routine,CKTag)
!The initial trial wavefuntion is set to zero        
        CALL AZZERO(CK,NDets*NEval*HElementSize)
        ALLOCATE(CKN(NDets,NEval),stat=ierr)
        CALL LogMemAlloc('CKN',NDets*NEval,8*HElementSize,this_routine,CKNTag)
        CALL AZZERO(CKN,NDets*NEval*HElementSize)

!Lanczos iterative diagonalisation routine
        CALL FRSBLKH(NDets,ICMax,NEval,Mat,Lab,CK,CKN,NKry,NKry1,NBlock,NRow,LScr,LIScr,A,W,V,AM,BM,T,WT,SCR,ISCR,Index,WH,Work2,V2,NCycle,B2L,.true.)
!Mulitply eigenvalues through by -1 to ensure they are positive
        CALL DSCAL(NEval,-1.D0,W,1)

!        do i=1,NDets
!            WRITE(6,*) NRow(i)
!        enddo
!        WRITE(6,*) "***********"

!        do i=1,NEval
!            WRITE(6,"(A,I4,A,F20.14)") "Eigenvalue ",i," is: ", W(i)
!        enddo

!Deallocate memory required by diagonaliser (including original matrix)
        DEALLOCATE(Mat)
        CALL LogMemDealloc(this_routine,MatTag)
        DEALLOCATE(Lab)
        CALL LogMemDealloc(this_routine,LabTag)
        DEALLOCATE(NRow)
        CALL LogMemDealloc(this_routine,NRowTag)
        DEALLOCATE(A)
        CALL LogMemDealloc(this_routine,ATag)
        DEALLOCATE(V)
        CALL LogMemDealloc(this_routine,VTag)
        DEALLOCATE(AM)
        CALL LogMemDealloc(this_routine,AMTag)
        DEALLOCATE(BM)
        CALL LogMemDealloc(this_routine,BMTag)
        DEALLOCATE(T)
        CALL LogMemDealloc(this_routine,TTag)
        DEALLOCATE(WT)
        CALL LogMemDealloc(this_routine,WTTag)
        DEALLOCATE(SCR)
        CALL LogMemDealloc(this_routine,ScrTag)
        DEALLOCATE(ISCR)
        CALL LogMemDealloc(this_routine,IScrTag)
        DEALLOCATE(Index)
        CALL LogMemDealloc(this_routine,IndexTag)
        DEALLOCATE(WH)
        CALL LogMemDealloc(this_routine,WHTag)
        DEALLOCATE(Work2)
        CALL LogMemDealloc(this_routine,Work2Tag)
        DEALLOCATE(V2)
        CALL LogMemDealloc(this_routine,V2Tag)

!Store largest eigenvector - ?first? column of CK (zero it)
        CALL AZZERO(Eigenvector,NDets*HElementSize)

        do i=1,NDets
            Eigenvector(i)=CK(i,1)
            WRITE(6,*) Eigenvector(i)
            IF(TMoveDets.and.(ABS(Eigenvector(i)%v)).eq.0.D0) THEN
!There are still the possibility of disconnected clusters - these will be removed by regrowing graph completly...
                IF(Iteration.eq.1) THEN
                    WRITE(6,*) "Disconnected cluster found in graph - performing one cycle of regrowing new graph from scratch..."
                    ReturntoTMoveDets=.true.
!Choose a low graph bias, which will allow the graph to be created easily
                    TBiasing=.true.
                    GraphBias=0.75
                ELSE
                    WRITE(6,*) "Disconnected clusters still found...exiting..."
                    STOP "Disconnected clusters still found...exiting..."
                ENDIF
            ENDIF
        enddo
        IF(ReturntoTMoveDets) THEN
            TMoveDets=.false.
        ENDIF

!Also need to save the largest Eigenvalue
        Eigenvalue=W(1)

!Find weight and energy of graph.
!There is no beta-dependance, so only largest eigenvector needed.
!Note - since we are not having a beta-dependance, 'weight' takes on a slightly different
!meaning. Here, the weight is simply the square of the first element of the largest eigenvector,
!i.e. the magnitude of the projection of the graph back onto the HF...
!        WRITE(6,*) "First element of eigenvector is: ", Eigenvector(1)%v
        SI=(Eigenvector(1)%v)*(Eigenvector(1)%v)
        DLWDB=0.D0
        do i=2,NDets
            DLWDB=DLWDB+(HamElems(i)%v)*(Eigenvector(i)%v)
        enddo
        DLWDB=(DLWDB/(Eigenvector(1)%v))+(HamElems(1)%v)

!Deallocate Eigenvectors and values
        DEALLOCATE(W)
        CALL LogMemDealloc(this_routine,WTag)
        DEALLOCATE(CK)
        CALL LogMemDealloc(this_routine,CKTag)
        DEALLOCATE(CKN)
        CALL LogMemDealloc(this_routine,CKNTag)
        WRITE(6,*) LabTag,ATag,MatTag,NRowTag,VTag,WTTag,AMTag,BMTag,TTag,SCRTag,ISCRTag,IndexTag,WHTag,Work2Tag,V2Tag,WTag,CKTag,CKNTag

        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in DiagGraphLanc'
        CALL TIHALT('DiagGraphLanc',iSubLanc)

    END SUBROUTINE DiagGraphLanc

!This compresses the rho matrix, so that only the non-zero elements are stored, in a suitable form for the Lanczos diagonaliser.
    SUBROUTINE CompressMatrix(Mat,Lab,NRow,LenMat,ICMax)
        IMPLICIT NONE
        INTEGER :: LenMat,i,j
        TYPE(HElement) :: Mat(NDets,ICMax)
        INTEGER :: Lab(NDets,ICMax),NRow(NDets),LabElem,RowElems,ICMax,SumElems

!This compression is how Alex uses the Lanczos diagonaliser - however, it seems to be used in a different way
!ICMax indicates the largest number of non-zero elements in any one row of the matrix
!        ICMax=0
!
!        LabElem=0
!        do i=1,NDets
!!Scan along each row of the matrix
!            RowElems=0
!            do j=1,NDets
!                IF(GraphRhoMat(i,j).agt.0.D0) THEN
!                    LabElem=LabElem+1
!                    RowElems=RowElems+1
!                    Mat(LabElem)=GraphRhoMat(i,j)
!                    Lab(LabElem)=j
!                ENDIF
!            enddo
!            IF(RowElems.gt.ICMax) ICMax=RowElems
!            NRow(i)=RowElems
!        enddo

        SumElems=0

        do i=1,NDets
            RowElems=0
            do j=1,NDets
                IF(GraphRhoMat(i,j).agt.0.D0) THEN
                    RowElems=RowElems+1
                    Mat(i,RowElems)=GraphRhoMat(i,j)
                    Lab(i,RowElems)=j
                ENDIF
            enddo
            IF(RowElems.gt.ICMax) THEN
                WRITE(6,*) "Error in compressing matrix 3"
                STOP 'Error in compressing matrix 3'
            ENDIF
            NRow(i)=RowElems
            SumElems=SumElems+RowElems
        enddo

        IF(SumElems.ne.LenMat) THEN
            WRITE(6,*) "Error in compressing matrix"
            STOP 'Error in compressing matrix'
        ENDIF

    END SUBROUTINE CompressMatrix

!This is a routine to find the energy of the graph by diagonalisation, and return as well its eigenvalues and largest eigenvector. Deallocate the RhoMatrix when done.
    SUBROUTINE DiagGraphMorph()
        IMPLICIT NONE
        INTEGER, SAVE :: iSubDiag
        INTEGER :: Info,ierr,i
        REAL*8 , ALLOCATABLE :: Work(:),Eigenvalues(:)
        INTEGER :: WorkTag,EigenvaluesTag,j
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

        IF(TMoveDets) THEN
!If using the MoveDets algorithm, need to keep the rho matrix for the previous graph
            ALLOCATE(CopyRhoMat(NDets,NDets),stat=ierr)
            CALL LogMemAlloc('CopyRhoMat',NDets**2,8,this_routine,CopyRhoMatTag)
            do i=1,NDets
                do j=i,NDets
!                    IF(i.eq.NDets) WRITE(6,*) j,GraphRhoMat(i,j)%v
!                    IF(j.eq.NDets) WRITE(6,*) i,GraphRhoMat(j,i)%v
                    CopyRhoMat(i,j)=GraphRhoMat(i,j)%v
                    CopyRhoMat(j,i)=GraphRhoMat(j,i)%v
                enddo
            enddo
        ENDIF

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
            IF(TMoveDets.and.(ABS(Eigenvector(i)%v)).eq.0.D0) THEN
!There are still the possibility of disconnected clusters - these will be removed by regrowing graph completly...
                IF(Iteration.eq.1) THEN
                    WRITE(6,*) "Disconnected cluster found in graph - performing one cycle of regrowing new graph from scratch..."
                    ReturntoTMoveDets=.true.
!Choose a low graph bias, which will allow the graph to be created easily
                    TBiasing=.true.
                    GraphBias=0.75
                ELSE
                    WRITE(6,*) "Disconnected clusters still found...exiting..."
                    STOP "Disconnected clusters still found...exiting..."
                ENDIF
            ENDIF
!            WRITE(6,*) i,Eigenvector(i)%v
        enddo
        IF(ReturntoTMoveDets) THEN
            TMoveDets=.false.
        ENDIF

!Also need to save the largest Eigenvalue
!        do i=1,5
!            WRITE(6,*) "Eigenvalue ", i," is: ",Eigenvalues(NDets-(i-1))
!        enddo
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
!        WRITE(6,*) "First element of eigenvector is: ", Eigenvector(1)%v
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

!LOGICAL FUNCTION ConnectGraph(Matrix,Dimen,NEl)
!    IMPLICIT NONE
!    INTEGER :: NotConnected(NEl,Dimen),NEl,Dimen
!    REAL*8 :: Matrix
!    LOGICAL :: ConnectGraph
!
!    j=0
!    do i=1,Dimen
!!Test if not connected to the excitation level above
!        IF(Matrix(1,i).eq.0) THEN
!            j=j+1
!            NotConnected(1,j)=i
!        ENDIF

REAL*8 FUNCTION RootofNum(Num,Root)
    INTEGER :: Root
    REAL*8 :: Num
    RootofNum=EXP(LOG(Num)/Root)
    RETURN
END FUNCTION RootofNum

!Tests determinants are the same - requires the same ordering of orbitals in them
LOGICAL FUNCTION SameDet(nI,nJ,NEl)
    IMPLICIT NONE
    INTEGER :: nI(NEl),nJ(NEl),NEl,i
    SameDet=.true.
    do i=1,NEl
        IF(nI(i).ne.nJ(i)) THEN
            SameDet=.false.
            EXIT
        ENDIF
    enddo
    RETURN
END FUNCTION
