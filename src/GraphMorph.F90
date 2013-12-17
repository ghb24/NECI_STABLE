!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!   TO DO   !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
!Fix MoveDets code
!Go to larger graphs, and look at possible convergence
!Look into MC sampling of rest of excitation space to aid convergence

!This code it designed to generate a random graph, with Ndets distinct 
!determinants in it (including HF), and improve it in an iterative fashion.
!Once a graph is found, all possible excitations from each determinant  in the 
!graph are found, and then the rho matrix from these additional
!excitations operate on the original eigenvector of the graph (null padded). 
!This produces a much larger vector (not eigenvector), from which
!excitations which have a strong connection can be chosen stochastically. This 
!process is then repeated with the new graph, 
!chosen from the vector (ensuring connectivity). This will improve the connections between the determinants, 
!which will hopefully lower the energy after only a few iterations.

MODULE GraphMorph
    
    use SystemData , only : NEl
    USE Determinants , only : FDet
!Iters is the number of interations of morphing the graph. Nd is the number of determinants in the graph.
    use CalcData , only : Iters,NDets,GraphBias,TBiasing,NoMoveDets,TMoveDets,TInitStar
    use CalcData , only : TNoSameExcit,TLanczos,TMaxExcit,iMaxExcitLevel,TOneExcitConn,THdiag
    use CalcData , only : TSinglesExcitSpace,TGrowInitGraph,GrowGraphsExpo
    USE LoggingData , only : TDistrib
    USE global_utilities
    use constants, only: dp
    IMPLICIT NONE
!    SAVE
!!This will hold all the determinants in the graph
!    INTEGER , ALLOCATABLE :: GraphDets(:,:)
!    INTEGER :: GraphDetsTag=0
!
!!This will hold all the Hamiltononian elements in the graph (zero if > double excit) (H_11) = HF
!    HElement_t , ALLOCATABLE :: HamElems(:)
!    INTEGER :: HamElemsTag=0
!    
!!This will hold all the determinants of the excitations of the graph
!    INTEGER , ALLOCATABLE :: ExcitsDets(:,:)
!    INTEGER :: ExcitsDetsTag=0
!
!!This is the number of possible excitations from each determinant in the graph, stored
!!in a cumulative way
!    INTEGER , ALLOCATABLE :: NoExcits(:)
!    INTEGER :: NoExcitsTag=0
!!The total number of excitations from the graph
!    INTEGER :: TotExcits
!
!!Holds the rho elements between each determinant and its excitations
!    HElement_t , ALLOCATABLE :: ConnectionsToExcits(:)
!    INTEGER :: ConnectionsToExcitsTag=0
!
!!The full Rho Matrix for the graph
!    HElement_t , ALLOCATABLE :: GraphRhoMat(:,:)
!    INTEGER :: GraphRhoMatTag=0
!!The Largest Eigenvector for the graph
!    HElement_t , ALLOCATABLE :: Eigenvector(:)
!    INTEGER :: EigenvectorTag=0
!!The largest Eigenvalue of the graph
!    real(dp) :: Eigenvalue
!
!!This is the vector of propensity to move towards the excited determinants of the graph
!    HElement_t , ALLOCATABLE :: ExcitsVector(:)
!    INTEGER :: ExcitsVectorTag=0
!
!!This is needed if we are using the MoveDets graph growing algorithm to store a copy of the rho matrix in
!    real(dp) , ALLOCATABLE :: CopyRhoMat(:,:)
!    INTEGER :: CopyRhoMatTag=0
!
!!If TDistrib is on, then this will show the distribution of determinants among the excitations
!    INTEGER , ALLOCATABLE :: Distribs(:,:)
!    INTEGER :: DistribsTag=0
!
!!If TNoSameExcit is on, then GraphExcitLevel stores the Excitation levels of the determinants in the graph
!!This is because we do not allow connections between determinants of the same excitation level.
!!We also require the array for TOneExcitConn for the same reasons
!    INTEGER , ALLOCATABLE :: GraphExcitLevel(:)
!    INTEGER :: GraphExcitLevelTag=0
!
!!The rho matrix and Hamiltonian element for the HF determinant
!    HElement_t :: rhii,Hii
!
!!These are the weight and energy of the current graph respectively
!    real(dp) :: SI,DLWDB
!
!!This is the seed for the random numbers needed in the routines
!    INTEGER :: Seed
!
!!This is the iteration of the GraphMorph that we are currently on...
!    INTEGER :: Iteration
!
!!This is needed in case a disconnected graph is generated - if so, a full new graph will be generated in the first MC
!    LOGICAL :: ReturntoTMoveDets
!
!!Various stats for printing
!    real(dp) :: PStay,Orig_Graph,SucRat,MeanExcit

    contains

    SUBROUTINE MorphGraph(Weight,Energyxw)
        use SystemData, only: Alat,Beta,Brr,ECore,G1,nBasis,nBasisMax,nMsh,Arr
        use CalcData , only : i_P,G_VMC_Seed
        use IntegralsData, only : fck,nMax,UMat,nTay
        USE Determinants , only : get_helement
        IMPLICIT NONE
        real(dp) :: Weight,Energyxw

        ! Avoid warnings
        weight = weight
        energyxw = energyxw 
!        real(dp) :: LowestE,BestSI
!        INTEGER :: ierr,i
!        CHARACTER(len=*), PARAMETER :: this_routine='MorphGraph'

!        OPEN(63,file='MorphStats',Status='unknown')
!        IF(TDistrib) THEN
!            OPEN(64,file='MorphDistrib',Status='unknown')
!        ENDIF
!        
!        IF(HElement_t_size.ne.1) STOP 'Only real orbitals allowed in GraphMorph so far'
!
!        IF(TMoveDets) THEN
!            WRITE(6,*) "***************************************"
!            WRITE(6,*) "*** WARNING - MOVEDETS MAY NOT WORK ***"
!            WRITE(6,*) "***************************************"
!        ENDIF
!        IF(TMoveDets.and.TOneExcitConn) THEN
!            WRITE(6,*) "TOneExcitConn is not yet implimented in TMoveDets"
!            STOP "TOneExcitConn is not yet implimented in TMoveDets"
!        ENDIF
!        IF(TMCExcits) THEN
!            WRITE(6,"(A,I10,A)") "Searching Excitations of each graph stochastically with ",NoMCExcits, " per determinant"
!        ENDIF
!
!        ReturntoTMoveDets=.false.
!!Initialise random number generator
!        Seed=G_VMC_Seed
!
!!Find rho_ii value, which all rho elements will be divided by, and Hii value
!        CALL CalcRho2(FDet,FDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhii,nTay,0,ECore)
!        Hii=GetHElement2(FDet,FDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)
!
!!It is first necessary to construct the original graph to work with.
!        IF(TInitStar) THEN
!!Let the initial graph be a complete connected star graph.
!            WRITE(6,"(A)") "Constructing random initial star graph to morph from."
!            CALL ConstructInitialStarGraph()
!        ELSEIF(TGrowInitGraph) THEN
!            WRITE(6,"(A,I8,A)") "Constructing initial graph in increasing excitation space with ",NDets, " vertices."
!            CALL ConstructExcitsInitGraph()
!        ELSE
!            WRITE(6,"(A,I8,A)") "Constructing random fully connected initial graph to morph from with ", NDets, " vertices."
!            CALL ConstructInitialGraph()
!        ENDIF
!
!        IF(TDistrib) THEN
!            WRITE(64,"(A8)",advance='no') "1"
!            do i=1,NEl
!                WRITE(64,"(I8)",advance='no') Distribs(i,1)
!            enddo
!            WRITE(64,*) ""
!            CALL neci_flush(64)
!        ENDIF
!        
!!Allocate space for largest Eigenvector for various graphs
!        ALLOCATE(Eigenvector(NDets),stat=ierr)
!        CALL LogMemAlloc('Eigenvector',NDets,8*HElement_t_size,this_routine,EigenvectorTag,ierr)
!
!        WRITE(63,*) "Iteration  Energy  P_stay  Orig_Graph  SucRat   MeanExcit"
!        IF(TBiasing) THEN
!            WRITE(6,"(A,F10.7)") "Graph growing biased towards original determinants with probability ", GraphBias
!        ELSEIF(TMoveDets) THEN
! WRITE(6,"(A,I3,A)") "Choosing new graph by moving ", NoMoveDets," determinants in previous graph to their excitations..."
!        ENDIF
!
!!Once the graph is found, the loop over the morphing iterations can begin
!        do Iteration=1,Iters
!
!            WRITE(6,*) ""
!            WRITE(6,"(A,I4,A)") "Starting Iteration ", Iteration, " ..."
!            WRITE(6,*) ""
!            CALL neci_flush(6)
!
!            IF((Iteration.eq.2).and.ReturntoTMoveDets) THEN
!!If one iteration of regrowing graphs is used, then return to the moving dets algorithm                
!                WRITE(6,*) "Returning to MoveDets algorithm"
!                TMoveDets=.true.
!                ReturntoTMoveDets=.false.
!                TBiasing=.false.
!            ENDIF
!
!!The graph is first diagonalised, and the energy of the graph found, along with the largest eigenvector and eigenvalues.
!            IF(TLanczos) THEN
!                CALL DiagGraphLanc()
!            ELSE
!                CALL DiagGraphMorph()
!            ENDIF
!
!            WRITE(6,"(A,2G20.12)") "Weight and Energy of current graph is: ", SI, DLWDB
!            IF(Iteration.eq.1) THEN
!                LowestE=DLWDB
!                BestSI=SI
!            ELSE
!                IF(DLWDB.lt.LowestE) THEN
!                    LowestE=DLWDB
!                    BestSI=SI
!                ENDIF
!            ENDIF
!            CALL neci_flush(6)
!
!!Write out stats
!            IF(Iteration.eq.1) THEN
!                WRITE(63,"(I10,G20.12,3A20,G20.12)") Iteration,DLWDB,"N/A","N/A","N/A",MeanExcit
!            ELSE
!                WRITE(63,"(I10,5G20.12)") Iteration,DLWDB,PStay,Orig_Graph,SucRat,MeanExcit
!            ENDIF
!            CALL neci_flush(63)
!
!!Excitation generators are initialised for each of the determinants in the graph, and the total number of possible
!!connected determinants calculated. Memory allocated for the ensuing calculation.
!            CALL CountExcits()
!!            WRITE(6,*) "Fraction of space which is space of excitations is: ", TotExcits/(TotExcits+NDets)
!            WRITE(6,*) "Total number of determinants available from current graph is: ",TotExcits
!!            CALL neci_flush(6)
!
!!Run through each determinant in the graph, calculating the excitations, storing them, and the rho elements to them.
!            CALL FindConnections()
!
!!Multiply the correct rho element, with the correct element of the eigenvector, to create the value for the 
!!tendancy to move to that excited determinant - store this in ExcitsVector
!            CALL CreateExcitsVector()
!
!            IF(TMoveDets) THEN
!!MoveDets means that a new graph-growing algorithm is used, where NoMoveDets determinants are 
!selected from the current graph, and replaced
!!by the same number of determinants from its excitations.
!
!!New arrays created - the inverse array from the original set of determinants, and the normalised array for the excitations
!                CALL NormaliseVectorSep()
!
!                CALL MoveDetsGraph()
!
!            ELSE
!
!!Add the original eigenvector*eigenvalue to the list of determinants - now have vector 
!(contained in Eigenvector*eigenvalue and ExcitsVector) with all
!!determinants in graph, and all possible determinants to excite to. This needs normalising.
!                CALL NormaliseVector()
!
!!Pick NDets new excitations stocastically from normalised list of determinants with probability |c|^2. Ensure connections,
!!allocate and create rho matrix for new graph. Deallocate info for old graph.
!                WRITE(6,*) "Choosing new graph stochastically from previous graph and its excitations..."
!                CALL neci_flush(6)
!                CALL PickNewDets()
!
!            ENDIF
!        
!            IF(TDistrib) THEN
!                WRITE(64,"(I8)",advance='no') Iteration+1
!                do i=1,NEl
!                    WRITE(64,"(I8)",advance='no') Distribs(i,Iteration+1)
!                enddo
!                WRITE(64,*) ""
!                CALL neci_flush(64)
!            ENDIF
!
!!Once graph is fully constructed, the next iteration can begin.
!        enddo
!        WRITE(6,*) ""
!
!!Diagonalise final graph, and find weight & energy. 
!!There is no beta-dependance, so only largest eigenvector and value needed.
!        IF(TLanczos) THEN
!            CALL DiagGraphLanc()
!        ELSE
!            CALL DiagGraphMorph()
!        ENDIF
!        WRITE(63,"(I10,5G20.12)") Iteration,DLWDB,PStay,Orig_Graph,SucRat,MeanExcit
!        
!        WRITE(6,"(A,2G20.12)") "Weight and Energy of final graph is: ", SI, DLWDB
!        IF(DLWDB.lt.LowestE) THEN
!            LowestE=DLWDB
!            BestSI=SI
!        ENDIF
!
!!Deallocate info...
!        DEALLOCATE(Eigenvector)
!        CALL LogMemDealloc(this_routine,EigenvectorTag)
!        DEALLOCATE(GraphDets)
!        CALL LogMemDealloc(this_routine,GraphDetsTag)
!        IF(.NOT.THDiag) THEN
!            DEALLOCATE(HamElems)
!            CALL LogMemDealloc(this_routine,HamElemsTag)
!        ENDIF
!        IF(TDistrib) THEN
!            DEALLOCATE(Distribs)
!            CALL LogMemDealloc(this_routine,DistribsTag)
!        ENDIF
!        IF(TNoSameExcit.or.TOneExcitConn) THEN
!            DEALLOCATE(GraphExcitLevel)
!            CALL LogMemDealloc(this_routine,GraphExcitLevelTag)
!        ENDIF
!
!        Weight=(BestSI-1.0_dp)
!        Energyxw=((LowestE*BestSI)-Hii)
!
!        Close(63)
!        IF(TDistrib) THEN
!            Close(64)
!        ENDIF
!        
!        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in ConstructInitialGraph'

    END SUBROUTINE MorphGraph

!!This routine constructs an initial graph in a non-stocastic manner by adding determinants from increasing excitation levels
!    SUBROUTINE ConstructExcitsInitGraph()
!        use SystemData , only : G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,nMsh,Arr
!        use CalcData , only : i_P,RhoEps
!        use IntegralsData , only : fck,nMax,UMat,nTay
!        USE Determinants , only : GetHElement2
!        IMPLICIT NONE
!        INTEGER :: nStore(6),exFlag,nExcitMemLen,iMaxExcit,nJ(NEl)
!        INTEGER , ALLOCATABLE :: nExcit(:)
!        INTEGER :: nExcitTag=0
!        HElement_t :: rh
!        INTEGER :: iExcit,excitcount,i,j,k,IC,ICRoot,Numberadded
!        CHARACTER(len=*), PARAMETER :: this_routine='ConstructExcitsInitGraph'
!        INTEGER :: ierr,Root,RootDet(NEl),iGetExcitLevel,NoNotAtt_Same,NoNotAtt_NoConn
!        type(timer), save :: proc_timerInitExcit
!        LOGICAL :: SameDet,Connection
!
!        proc_timerInitExcit%timer_name='InitExcitGraph'
!        call set_timer(proc_timerInitExcit)
!
!!Allow single and double excitations
!        exFlag=3
!
!        IF(TDistrib) THEN
!!We want to record the number of each excitation level for each determinant in the graph
!            ALLOCATE(Distribs(NEl,Iters+1),stat=ierr)
!            CALL LogMemAlloc('Distribs',NEl*(Iters+1),4,this_routine,DistribsTag,ierr)
!            Distribs(1:NEl,1:(Iters+1))=0
!        ENDIF
!
!!Allocate space for graph determinants and rho matrix
!        ALLOCATE(GraphDets(NDets,NEl),stat=ierr)
!        CALL LogMemAlloc('GraphDets',NDets*NEl,4,this_routine,GraphDetsTag,ierr)
!        GraphDets(1:NDets,1:NEl)=0
!        ALLOCATE(GraphRhoMat(NDets,NDets),stat=ierr)
!        CALL LogMemAlloc('GraphRhoMat',NDets*NDets,8*HElement_t_size,this_routine,GraphRhoMatTag,ierr)
!        GraphRhoMat=(0.0_dp)
!        IF(.NOT.THDiag) THEN
!            ALLOCATE(HamElems(NDets),stat=ierr)
!            CALL LogMemAlloc('HamElems',NDets,8*HElement_t_size,this_routine,HamElemsTag,ierr)
!            HamElems=(0.0_dp)
!        ENDIF
!        IF(TNoSameExcit.or.TOneExcitConn) THEN
!            ALLOCATE(GraphExcitLevel(NDets),stat=ierr)
!            CALL LogMemAlloc('GraphExcitLevel',NDets,4,this_routine,GraphExcitLevelTag,ierr)
!            GraphExcitLevel(1:NDets)=0
!        ENDIF
!
!!Put HF Determinant into first element of GraphDets
!        do i=1,NEl
!            GraphDets(1,i)=FDet(i)
!        enddo
!        IF(THDiag) THEN
!            GraphRhoMat(1,1)=Hii
!        ELSE
!            GraphRhoMat(1,1)=rhii
!            HamElems(1)=Hii
!        ENDIF
!        MeanExcit=0.0_dp
!
!        i=1
!        Root=1
!!        Numberadded=0
!!        NoNotAtt_Same=0
!!        NoNotAtt_NoConn=0
!        do while(i.lt.NDets)
!!Cycle through all excitations consecutivly, adding them where possible
!
!!            WRITE(6,*) "Number added for root ", Root-1," is = ", Numberadded
!!            WRITE(6,"(A,I5,A,I8)") "Number not added for root ", Root-1," due to already in graph = ", NoNotAtt_Same
!!            WRITE(6,"(A,I5,A,I8)") "Number not added for root ", Root-1," due to no connection = ", NoNotAtt_NoConn
!!            Numberadded=0
!!            NoNotAtt_Same=0
!!            NoNotAtt_NoConn=0
!            IF(Root.gt.i) THEN
!                WRITE(6,*) "Error - trying to make an unavailable determinant root"
!                WRITE(6,*) "Trying to select a root of ", Root
!                STOP "Error - trying to make an unavailable determinant root"
!            ENDIF
!
!!Let root of excitation generator be next available determinant
!            do j=1,NEl
!                RootDet(j)=GraphDets(Root,j)
!                IF(RootDet(j).eq.0) THEN
!                    WRITE(6,*) "Incorrect assignment of root determinant"
!                    STOP "Incorrect assignment of root determinant"
!                ENDIF
!            enddo
!
!!Setup excitation generator for new root
!            nStore(1:6)=0
!            CALL GenSymExcitIt2(RootDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
!            ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
!            CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag,ierr)
!            nExcit(1)=0
!            CALL GenSymExcitIt2(RootDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
!
!            do while(i.lt.NDets)
!
!!                WRITE(6,*) i
!!                CALL neci_flush(6)
!                CALL GenSymExcitIt2(RootDet,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,iExcit,0,nStore,exFlag)
!                IF(nJ(1).eq.0) THEN
!!In the next sweep, look at the next determinant being the root
!                    Root=Root+1
!!Deallocate excitation generators
!                    DEALLOCATE(nExcit)
!                    CALL LogMemDealloc(this_routine,nExcitTag)
!                    EXIT
!                ENDIF
!
!                Connection=.false.
!                do j=1,i
!!Look through determinants already in graph to check for double counting and connections
!                    SameDet=.true.
!                    do k=1,NEl
!                        IF(nJ(k).ne.GraphDets(j,k)) THEN
!                            SameDet=.false.
!                        ENDIF
!                    enddo
!!If we have found the determinant somewhere else in the graph, then discount it
!                    IF(SameDet) EXIT
!
!                    IF(.not.Connection) THEN
!                        IF(THDiag) THEN
!                            rh=GetHElement2(GraphDets(j,:),nJ(:),NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,-1,ECore)
!                        ELSE
!  CALL CalcRho2(GraphDets(j,:),nJ(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,-1,ECore)
!                        ENDIF
!                        IF(abs(rh).gt.0.0_dp) THEN
!!A connection has been found - we do not need to look for any others
!                            Connection=.true.
!                        ENDIF
!                    ENDIF
!                enddo
!
!!                IF(.not.Connection) NoNotAtt_NoConn=NoNotAtt_NoConn+1
!!                IF(SameDet) NoNotAtt_Same=NoNotAtt_Same+1
!                
!                IF(Connection.and.(.not.SameDet)) THEN
!!A valid determinant has been found - add it
!!                    Numberadded=Numberadded+1
!                    i=i+1
!                    do j=1,NEl
!                        GraphDets(i,j)=nJ(j)
!                    enddo
!
!!Find diagonal rho matrix element, and connection to HF for 
!                    IF(THDiag) THEN
!                        rh=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,0,ECore)
!                    ELSE
!                        CALL CalcRho2(nJ,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,0,ECore)
!                    ENDIF
!                    GraphRhoMat(i,i)=rh
!                    IF(Root.eq.1) THEN
!                        ICRoot=iExcit
!                    ELSE
!                        ICRoot=iGetExcitLevel(FDet,nJ,NEl)
!                    ENDIF
!                    MeanExcit=MeanExcit+(ICRoot+0.0_dp)
!                    IF(ICRoot.gt.2) THEN
!                        IF(.NOT.THDiag) HamElems(i)=(0.0_dp)
!                        GraphRhoMat(i,1)=(0.0_dp)
!                        GraphRhoMat(1,i)=(0.0_dp)
!                    ELSE
!                        IF(THDiag) THEN
!                            rh=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,ICRoot,ECore)
!                        ELSE
!                            HamElems(i)=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,ICRoot,ECore)
!                            CALL CalcRho2(FDet(:),nJ(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,
                                    !fck,Arr,ALat,UMat,rh,nTay,ICRoot,ECore)
!                        ENDIF
!                        GraphRhoMat(i,1)=rh
!                        GraphRhoMat(1,i)=rh
!                    ENDIF
!!Store the excitation level from HF if we are removing excitations from the same level
!                    IF(TNoSameExcit.or.TOneExcitConn) GraphExcitLevel(i)=ICRoot
!
!                    IF(TDistrib) THEN
!!Need to add excitation level of determinant. Also, store the excitation level of the determinants in the current graph
!                        Distribs(ICRoot,1)=Distribs(ICRoot,1)+1
!                    ENDIF
!
!                    do j=2,i-1
!!Need to find connection to all other determinants
!            
!                        IF(TNoSameExcit) THEN
!!Don't allow connections between excitations of the same level
!
!                            IF(GraphExcitLevel(j).ne.ICRoot) THEN
!                                IC=iGetExcitLevel(GraphDets(j,:),nJ(:),NEl)
!                                IF(THDiag) THEN
!                                    rh=GetHElement2(GraphDets(j,:),nJ(:),NEl,nBasisMax,G1,nBasis,Brr,
                                            !nMsh,fck,nMax,ALat,UMat,IC,ECore)
!                                ELSE
!                                    CALL CalcRho2(GraphDets(j,:),nJ(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,
                                        !fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
!                                ENDIF
!                                GraphRhoMat(i,j)=rh
!                                GraphRhoMat(j,i)=rh
!                            ELSE
!                                GraphRhoMat(i,j)=(0.0_dp)
!                                GraphRhoMat(j,i)=(0.0_dp)
!                            ENDIF
!
!                        ELSEIF(TOneExcitConn) THEN
!!Only allow connections between determinants which differ in excitation level by one.
!
!                            IF((ABS(GraphExcitLevel(j)-ICRoot)).eq.1) THEN
!                                IC=iGetExcitLevel(GraphDets(j,:),nJ(:),NEl)
!                                IF(THDiag) THEN
!                                    rh=GetHElement2(GraphDets(j,:),nJ(:),NEl,nBasisMax,G1,nBasis,Brr,nMsh,
                                                    !fck,nMax,ALat,UMat,IC,ECore)
!                                ELSE
!                                    CALL CalcRho2(GraphDets(j,:),nJ(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,
                                                    !nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
!                                ENDIF
!                                GraphRhoMat(i,j)=rh
!                                GraphRhoMat(j,i)=rh
!                            ELSE
!                                GraphRhoMat(i,j)=(0.0_dp)
!                                GraphRhoMat(j,i)=(0.0_dp)
!                            ENDIF
!
!                        ELSE
!                            
!                            IC=iGetExcitLevel(GraphDets(j,:),nJ(:),NEl)
!!Fully connect the graph
!                            IF(THDiag) THEN
!                                rh=GetHElement2(GraphDets(j,:),nJ(:),NEl,nBasisMax,G1,nBasis,Brr,
                                                !nMsh,fck,nMax,ALat,UMat,IC,ECore)
!                            ELSE
!                                CALL CalcRho2(GraphDets(j,:),nJ(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,
                                    !nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
!                            ENDIF
!                            GraphRhoMat(i,j)=rh
!                            GraphRhoMat(j,i)=rh
!
!                        ENDIF
!
!!End of searching for connections
!                    enddo
!!Whether we connect the determinant created to the graph
!                ENDIF
!!End of loop for a given root
!            enddo
!!End of searching for excitation - all found
!        enddo
!                    
!!Deallocate excitation generators
!        DEALLOCATE(nExcit)
!        CALL LogMemDealloc(this_routine,nExcitTag)
!
!        IF(i.ne.NDets) THEN
!            WRITE(6,*) "Graph Growing routine has not found all determinant requested"
!            STOP "Graph Growing routine has not found all determinant requested"
!        ENDIF
!
!        WRITE(6,"(A,I4)") "Number of roots needed to create initial graph = ", Root
!
!        MeanExcit=MeanExcit/(NDets-1)
!
!        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in ConstructExcitsInitGraph'
!        call halt_timer(proc_timerInitExcit)
!
!    END SUBROUTINE ConstructExcitsInitGraph
!
!
!!This routine constructs the complete connected star graph as the initial graph to begin morphing from
!    SUBROUTINE ConstructInitialStarGraph()
!        use SystemData , only : G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,nMsh,Arr
!        use CalcData , only : i_P,RhoEps
!        use IntegralsData , only : fck,nMax,UMat,nTay
!        USE Determinants , only : GetHElement2
!        IMPLICIT NONE
!        INTEGER :: ierr,nStore(6),exFlag,nExcitMemLen,iMaxExcit,nJ(NEl)
!        type(timer), save :: proc_timerInitStar
!        INTEGER , ALLOCATABLE :: nExcit(:)
!        INTEGER :: nExcitTag=0
!        HElement_t :: rh
!        INTEGER :: iExcit,excitcount,i,j
!        CHARACTER(len=*), PARAMETER :: this_routine='ConstructInitialStarGraph'
!        
!        proc_timerInitStar%timer_name='InitStarGraph'
!        call set_timer(proc_timerInitStar)
!        
!        nStore(1:6)=0
!!Having exFlag=2 means that only double excitations are generated
!        exFlag=2
!        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
!        ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
!        CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag,ierr)
!        nExcit(1:nExcitMemLen)=0
!        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
!!First run through all excitations, so we do not store more excitations than necessary
!        
!        excitcount=0
!        do while(.true.)
!            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,iExcit,0,nStore,exFlag)
!            IF(nJ(1).eq.0) exit
!            IF(THDiag) THEN
!                rh=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,iExcit,ECore)
!            ELSE
!                CALL CalcRho2(FDet,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,iExcit,ECore)
!            ENDIF
!            IF(abs(rh).gt.0.0_dp) excitcount=excitcount+1
!        enddo
!
!        DEALLOCATE(nExcit)
!        nStore(1:6)=0
!        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
!        ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
!        nExcit(1:nExcitMemLen)=0
!        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
!        iMaxExcit=excitcount
!!NDets is equal to the number of excitations plus 1 (for root)
!        NDets=excitcount+1
!
!        WRITE(6,"(A,I10)") "Total number of determinants in GraphMorph is ", NDets
!        
!        IF(TDistrib) THEN
!!If TDistribs is true, then we are recording the change in distributions of the graphs over the iterations
!            ALLOCATE(Distribs(NEl,Iters+1),stat=ierr)
!            CALL LogMemAlloc('Distribs',NEl*(Iters+1),4,this_routine,DistribsTag,ierr)
!            Distribs(1:NEl,1:(Iters+1))=0
!        ENDIF
!        IF(TNoSameExcit.or.TOneExcitConn) THEN
!!We want to store the excitation levels of the determinants in the graph
!            ALLOCATE(GraphExcitLevel(NDets),stat=ierr)
!            CALL LogMemAlloc('GraphExcitLevel',NDets,4,this_routine,GraphExcitLevelTag,ierr)
!!This will automatically set the first element (HF) to zero
!            GraphExcitLevel(1:NDets)=0
!        ENDIF
!
!!Allocate memory needed for calculations
!        ALLOCATE(GraphDets(NDets,NEl),stat=ierr)
!        CALL LogMemAlloc('GraphDets',NDets*NEl,4,this_routine,GraphDetsTag,ierr)
!        GraphDets(1:NDets,1:NEl)=0
!        ALLOCATE(GraphRhoMat(NDets,NDets),stat=ierr)
!        CALL LogMemAlloc('GraphRhoMat',NDets*NDets,8*HElement_t_size,this_routine,GraphRhoMatTag,ierr)
!        GraphRhoMat=(0.0_dp)
!        IF(.NOT.THDiag) THEN
!            ALLOCATE(HamElems(NDets),stat=ierr)
!            CALL LogMemAlloc('HamElems',NDets,8*HElement_t_size,this_routine,HamElemsTag,ierr)
!            HamElems=(0.0_dp)
!            HamElems(1)=Hii
!        ENDIF
!
!        GraphRhoMat(1,1)=rhii
!        do i=1,NEl
!            GraphDets(1,i)=FDet(i)
!        enddo
!        MeanExcit=0.0_dp
!        
!        i=1
!!Run through excitations
!        do while(.true.)
!            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,iExcit,0,nStore,exFlag)
!            IF(nJ(1).eq.0) EXIT
!
!!Since we already know that the excitations are double excitations of FDet, we can put 
                            !those connections in seperatly (will be quicker)
!            IF(.NOT.THDiag) THEN
!                CALL CalcRho2(FDet,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,iExcit,ECore)
!            ELSE
!                rh=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,iExcit,ECore)
!            ENDIF
!
!            IF(abs(.not.(rh).gt.0.0_dp)) CYCLE
!            i=i+1
!
!            GraphRhoMat(1,i)=rh
!            GraphRhoMat(i,1)=rh
!            
!!Store Path
!            do j=1,NEl
!                GraphDets(i,j)=nJ(j)
!            enddo
!
!!Find hamiltonian element coupling to FDet, and diagonal element of excitation
!            IF(.NOT.THDiag) THEN
!                HamElems(i)=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,iExcit,ECore)
!                CALL CalcRho2(nJ,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,GraphRhoMat(i,i),nTay,0,ECore)
!            ELSE
!                GraphRhoMat(i,i)=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,0,ECore)
!            ENDIF
!            
!            MeanExcit=MeanExcit+iExcit
!
!            IF(TDistrib) THEN
!                Distribs(iExcit,1)=Distribs(iExcit,1)+1
!            ENDIF
!            IF(TNoSameExcit.or.TOneExcitConn) GraphExcitLevel(i)=iExcit
!
!!Cycle through all excitations already generated to determine coupling to them
!            IF((.not.TNoSameExcit).and.(.not.TOneExcitConn)) THEN
!!If TNoSameExcit is on, then we are ignoring these crosslinks - the normal star graph should result
!                do j=2,(i-1)
!                    IF(THDiag) THEN
!                        rh=GetHElement2(GraphDets(j,:),nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,-1,ECore)
!                    ELSE
!                        CALL CalcRho2(GraphDets(j,:),nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,
                                        !fck,Arr,ALat,UMat,rh,nTay,-1,ECore)
!                    ENDIF
!                    GraphRhoMat(i,j)=rh
!                    GraphRhoMat(j,i)=rh
!                enddo
!            ENDIF
!
!        enddo
!
!        MeanExcit=MeanExcit/(i-1)
!
!        IF((ABS(MeanExcit-2.0_dp)).gt.1.0e-8_dp) THEN
!            WRITE(6,*) "Error generating initial star graph"
!            STOP 'Error generating initial star graph'
!        ENDIF
!
!        DEALLOCATE(nExcit)
!        CALL LogMemDealloc(this_routine,nExcitTag)
!
!        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in ConstructInitialStarGraph'
!        call halt_timer(proc_timerInitStar)
!
!    END SUBROUTINE ConstructInitialStarGraph
!
!    SUBROUTINE ConstructInitialGraph()
!        use SystemData , only : G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,nMsh,Arr
!        use CalcData , only : i_P,RhoEps,G_VMC_Pi,G_VMC_Seed
!        use IntegralsData , only : fck,nMax,UMat,nTay
!        USE Determinants , only : GetHElement2
!        IMPLICIT NONE
!        LOGICAL :: Attached,connected
!        INTEGER :: ierr,nStore(6),nExcitMemLen,iMaxExcit,nJ(NEl),nExcitTag,iExcit
!        INTEGER :: iPathTag,XijTag,RhoiiTag,RhoijTag,HijsTag,i,j,diff,k
!        type(timer), save :: proc_timerInit
!        INTEGER :: iGetExcitLevel,IC,MatFlag
!        INTEGER , ALLOCATABLE :: iPath(:,:)
!        real(dp) , ALLOCATABLE :: Xij(:,:)
!#if defined(POINTER8)
!        integer(int64) :: ExcitGen(0:NDets)
!#else
!        INTEGER :: ExcitGen(0:NDets)
!#endif
!        real(dp) , ALLOCATABLE :: Rhoii(:)
!        HElement_t , ALLOCATABLE :: Rhoij(:,:),Hijs(:)
!        HElement_t :: rh
!        INTEGER , ALLOCATABLE :: nExcit(:)
!        real(dp) :: PGen,OldImport
!        CHARACTER(len=*), PARAMETER :: this_routine='ConstructInitialGraph'
!        
!        proc_timerInit%timer_name='ConsInitGraph'
!        call set_timer(proc_timerInit)
!
!        IF(THDiag) THEN
!            MatFlag=-19
!        ELSE
!            MatFlag=0
!        ENDIF
!        
!!        WRITE(6,*) "FDET is ",FDet(:)
!
!!Set the importance parameter to be equal to 1 if we want random double excitation connected star graphs as initial graphs.
!!        OldImport=G_VMC_Pi
!!        G_VMC_Pi=1.0_dp
!        
!        IF(TDistrib) THEN
!!If TDistribs is true, then we are recording the change in distributions of the graphs over the iterations
!            ALLOCATE(Distribs(NEl,Iters+1),stat=ierr)
!            CALL LogMemAlloc('Distribs',NEl*(Iters+1),4,this_routine,DistribsTag,ierr)
!            Distribs(1:NEl,1:(Iters+1))=0
!        ENDIF
!
!!Set Tags to zero, so we know when they are allocated/deallocated
!        nExcitTag=0
!        iPathTag=0
!        XijTag=0
!        RhoiiTag=0
!        RhoijTag=0
!        HijsTag=0
!!Setup excitation generator
!        nStore(1:6)=0
!        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
!        ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
!        CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag,ierr)
!        nExcit(1:nExcitMemLen)=0
!        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,3)
!        CALL GenRandSymExcitIt2(FDet,NEl,nExcit,nJ,Seed,iExcit,PGen)
!
!!Allocate memory for graph
!        ALLOCATE(iPath(NEl,0:NDets),stat=ierr)
!        CALL LogMemAlloc('iPath',NEl*(NDets+1),4,this_routine,iPathTag,ierr)
!        iPath(1:NEl,0:NDets+1)=0
!        ALLOCATE(Xij(0:NDets-1,0:NDets-1),stat=ierr)
!        CALL LogMemAlloc('Xij',NDets*NDets,8,this_routine,XijTag,ierr)
!        Xij=0.0_dp
!        ALLOCATE(Rhoii(0:NDets),stat=ierr)
!        CALL LogMemAlloc('Rhoii',NDets+1,8,this_routine,RhoiiTag,ierr)
!        Rhoii=(0.0_dp)
!        ALLOCATE(Rhoij(0:NDets,0:NDets),stat=ierr)
!        CALL LogMemAlloc('Rhoij',(NDets+1)*(NDets+1),8*HElement_t_size,this_routine,RhoijTag,ierr)
!        Rhoij=(0.0_dp)
!        ALLOCATE(Hijs(0:NDets),stat=ierr)
!        CALL LogMemAlloc('Hijs',NDets+1,8*HElement_t_size,this_routine,HijsTag,ierr)
!        Hijs=(0.0_dp)
!
!!The first and last determinants of iPath want to be the FDet...
!        do i=1,NEl
!            iPath(i,0)=FDet(i)
!            iPath(i,NDets)=FDet(i)
!        enddo
!
!!ExcitGen is an array of pointers which points to the excitation generators for each vertex of the graph.
!!Is used, and the memory is deallocated in the graph generation algorithm, so do not need to worry about it.
!        ExcitGen(:)=0
!
!        IF(TMoveDets) THEN
!            connected=.false.
!            k=1
!            do while (.not.connected)
!!Generate the initial graph...(Can speed this up as do not need to know prob, and are recalculating the rho matrix)
!!Arr sent instead of NMax since getting values straight from modules (not passed through)
!                CALL Fmcpr4d2GenGraph(FDet,NEl,Beta,i_P,iPath,NDets,Xij,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,Alat,
                            !UMat,nTay,RhoEps,Rhoii,Rhoij,ECore,Seed,Hijs,nExcit,MatFlag,ExcitGen)
!!Need to check that all determinants in the graph are connected (if we are using movedets)
!                do i=1,NDets
!                    connected=.false.
!                    do j=1,NDets
!                        IF(i.eq.j) CYCLE
!                        IF(abs(Rhoij(i-1,j-1)).gt.0.0_dp) THEN
!                            connected=.true.
!                            EXIT
!                        ENDIF
!                    enddo
!                    IF(.not.connected) EXIT
!                enddo
!                IF(.not.connected) THEN
!                    IF(k.gt.15) THEN
!                        WRITE(6,*) "A fully connected initial graph could not be created - should now recreate 
                                                    !the graph from scratch for one iteration..."
!!                        STOP "A fully connected initial graph could not be created - exiting..."
!                        EXIT
!                    ENDIF
!                    k=k+1
!                    WRITE(6,"(A,I7,A)") "Determinant number ",i," not connected to rest of initial graph - 
                                                !trying again to create initial graph with new seed."
!!change seed and repeat
!                    seed=G_VMC_Seed-k
!!Resetup arguments and excitation generator for fmcpr4d2gengraph
!                    DEALLOCATE(nExcit)
!                    CALL LogMemDealloc(this_routine,nExcitTag)
!                    nStore(1:6)=0
!                    CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
!                    ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
!                    CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag,ierr)
!                    nExcit(1:nExcitMemLen)=0
!                    CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,3)
!                    CALL GenRandSymExcitIt2(FDet,NEl,nExcit,nJ,Seed,iExcit,PGen)
!                    iPath(1:NEl,0:NDets)=0
!                    Xij=0.0_dp
!                    Rhoii=(0.0_dp)
!                    Rhoij=(0.0_dp)
!                    Hijs=(0.0_dp)
!                    do i=1,NEl
!                        iPath(i,0)=FDet(i)
!                        iPath(i,NDets)=FDet(i)
!                    enddo
!                    ExcitGen(:)=0
!                ENDIF
!            enddo
!!            IF(.not.connected) THEN
!!!Create a small artificial connection between disconnected determinants and root
!!                do i=1,NDets
!!                    connected=.false.
!!                    do while(.not.connected)
!!                        do j=1,NDets
!!                            IF(i.eq.j) CYCLE
!!                            IF(abs(Rhoij(i-1,j-1)).gt.0.0_dp) THEN
!!                                connected=.true.
!!                                EXIT
!!                            ENDIF
!!                        enddo
!!                        IF(.not.connected) THEN
!!                            Rhoij(i-1,1)=(1.0e-20_dp)
!!                            Rhoij(1,i-1)=(1.0e-20_dp)
!!                        ENDIF
!!                    enddo
!!                enddo
!!            ENDIF
!        ELSE
!!Generate the initial graph...(Can speed this up as do not need to know prob, and are recalculating the rho matrix)
!!Arr sent instead of NMax since getting values straight from modules (not passed through)
!            CALL Fmcpr4d2GenGraph(FDet,NEl,Beta,i_P,iPath,NDets,Xij,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,
                        !Alat,UMat,nTay,RhoEps,Rhoii,Rhoij,ECore,Seed,Hijs,nExcit,MatFlag,ExcitGen)
!        ENDIF
!
!        IF(THDiag) THEN
!            DEALLOCATE(Hijs)
!            CALL LogMemDealloc(this_routine,HijsTag)
!        ENDIF
!
!!Store initial graph determinants in the GraphDets array
!        ALLOCATE(GraphDets(NDets,NEl),stat=ierr)
!        CALL LogMemAlloc('GraphDets',NDets*NEl,4,this_routine,GraphDetsTag,ierr)
!        GraphDets(1:NDets,1:NEl)=0
!        do i=1,NDets
!            do j=1,NEl
!                GraphDets(i,j)=iPath(j,i-1)
!                IF(GraphDets(i,j).eq.0) STOP 'Error in creating initial graph'
!            enddo
!        enddo
!        DEALLOCATE(iPath)
!        CALL LogMemDealloc(this_routine,iPathTag)
!
!!Xij is the probability matrix - we do not need this information...
!        DEALLOCATE(Xij)
!        CALL LogMemDealloc(this_routine,XijTag)
!
!!Rho matrix is stored as paths, and so do not need last column and row (should be same as first) - same with H elements
!        ALLOCATE(GraphRhoMat(NDets,NDets),stat=ierr)
!        CALL LogMemAlloc('GraphRhoMat',NDets*NDets,8*HElement_t_size,this_routine,GraphRhoMatTag,ierr)
!        GraphRhoMat=(0.0_dp)
!
!!Trust rho matrix from Fmcpr4d2GenGraph - no need to recalculate...
!        do i=1,NDets
!            do j=1,NDets
!!                WRITE(6,"(E8.5)",advance='no') Rhoij(i-1,j-1)
!                GraphRhoMat(i,j)=Rhoij(i-1,j-1)
!            enddo
!!            WRITE(6,*) ""
!        enddo
!
!!        do j=1,NDets
!!            Attached=.false.
!!            do i=1,NDets
!!                IF(i.eq.j) CYCLE
!!                IF((GraphRhoMat(i,j)).ne.0.0_dp) THEN
!!                    Attached=.true.
!!                    EXIT
!!                ENDIF
!!            enddo
!!            IF(.not.attached) STOP 'determinant is not attached!'
!!        enddo
!
!!To make sure, put in diagonal elements separatly
!        do i=1,NDets
!            GraphRhoMat(i,i)=Rhoii(i-1)
!        enddo
!
!!First element is not put in for you...
!        IF(THDiag) THEN
!            GraphRhoMat(1,1)=Hii
!        ELSE
!            GraphRhoMat(1,1)=rhii
!        ENDIF
!
!!        do i=1,NDets
!!            do j=i,NDets
!!                IF(i.eq.j) THEN
!!                    diff=0
!!                ELSE
!!                    diff=-1
!!                ENDIF
!!                CALL CalcRho2(GraphDets(i,:),GraphDets(j,:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,
                                        !ALat,UMat,rh,nTay,diff,ECore)
!!                GraphRhoMat(i,j)=rh
!!                GraphRhoMat(j,i)=rh
!!            enddo
!!        enddo
!
!!        IF(DREAL(GraphRhoMat(1,1)).ne.(rhii)) THEN
!!This could be because the elements haven't been divided by rhii - check value of rhii & then divide by it...
!!            STOP 'Rho matrix elements for initial graph incorrect'
!!        ENDIF
!
!!Write out rho matrix - debugging
!!        do i=1,NDets
!!            do j=1,NDets
!!                WRITE(6,"(E14.6)",advance='no') GraphRhoMat(i,j)
!!            enddo
!!            WRITE(6,*) ""
!!        enddo
!
!        IF(.NOT.THDiag) THEN
!            ALLOCATE(HamElems(NDets),stat=ierr)
!            CALL LogMemAlloc('HamElems',NDets,8*HElement_t_size,this_routine,HamElemsTag,ierr)
!            HamElems=(0.0_dp)
!
!!Trust the fmcpr4d2gengraph Hij values. Recalculated ones are identical
!            do i=1,NDets
!                HamElems(i)=Hijs(i-1)
!            enddo
!            HamElems(1)=Hii
!
!!        do i=1,NDets
!!            HamElems(i)=GetHElement2(FDet,GraphDets(i,:),NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,-1,ECore)
!!        enddo
!        
!            IF((HamElems(1)).ne.(Hii)) THEN
!                STOP 'H elements for initial graph incorrect'
!            ENDIF
!
!            DEALLOCATE(Hijs)
!            CALL LogMemDealloc(this_routine,HijsTag)
!        ENDIF
!
!!Find out manually the mean distance from the HF determinant...
!        MeanExcit=0.0_dp
!        do i=2,NDets
!            IC=iGetExcitLevel(FDet(:),GraphDets(i,:),NEl)
!            MeanExcit=MeanExcit+IC
!            IF(TDistrib) THEN
!                Distribs(IC,1)=Distribs(IC,1)+1
!            ENDIF
!        enddo
!        MeanExcit=MeanExcit/(NDets-1.0_dp)
!
!!Return G_VMC_Pi to original value (Just in case it wants to be used later)
!!        G_VMC_Pi=OldImport
!        
!        DEALLOCATE(Rhoii)
!        CALL LogMemDealloc(this_routine,RhoiiTag)
!        DEALLOCATE(Rhoij)
!        CALL LogMemDealloc(this_routine,RhoijTag)
!        DEALLOCATE(nExcit)
!        CALL LogMemDealloc(this_routine,nExcitTag)
!
!        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in ConstructInitialGraph'
!        call halt_timer(proc_timerInit)
!
!    END SUBROUTINE ConstructInitialGraph
!
!!This routine move determinants already in the graph, to excitations stocastically
!!TO DO: since this only moves a few determinants, much of the excitation space is still the same
!!this means that only a few of the excitations would need to be reproduced. However, bookkeeping
!!may well be a nightmare
!    SUBROUTINE MoveDetsGraph()
!        use SystemData , only : G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,nMsh,Arr
!        use CalcData , only : i_P,RhoEps
!        use IntegralsData , only : fck,nMax,UMat,nTay
!        USE Determinants , only : GetHElement2
!        IMPLICIT NONE
!        INTEGER :: i,j,k,l,Success,Failure,iGetExcitLevel,IC,Excitation,IC1,IC2
!        type(timer), save :: proc_timerMove 
!        INTEGER , ALLOCATABLE :: MoveDetsFromPaths(:,:)
!        INTEGER :: MoveDetsFromPathsTag=0
!        INTEGER :: AttemptDet(NEl),IndexofDetsFrom(NoMoveDets),ierr,Tries,NoVerts
!        LOGICAL :: Remove,Attach,SameDet
!        HElement_t :: rh
!        real(dp) :: r,Ran2
!        CHARACTER(len=*), PARAMETER :: this_routine='MoveDets'
!
!        proc_timerMove%timer_name='MoveDets'
!        call set_timer(proc_timerMove)
!
!!Allocate space for the determinants to be moved, and to move to
!        ALLOCATE(MoveDetsFromPaths(NoMoveDets,NEl),stat=ierr)
!        CALL LogMemAlloc('MoveDetsFromPaths',NoMoveDets*NEl,4,this_routine,MoveDetsFromPathsTag,ierr)
!        MoveDetsFromPaths(1:NoMoveDets,1:NEl)=0
!        IndexofDetsFrom(1:NoMoveDets)=0
!
!!Prepare for the change in MeanExcits
!        MeanExcit=MeanExcit*(NDets-1)
!        IF(TDistrib) THEN
!!If we are recording the distributions, then we want to simply have the same distribution for the previous
!!iteration, and then subtract the determinants that we remove
!            do i=1,NEl
!                Distribs(i,Iteration+1)=Distribs(i,Iteration)
!            enddo
!        ENDIF
!
!!First need to pick NoMoveDets determinants stochastically to be moved from the original graph
!        Tries=NoMoveDets*10000*NDets
!        k=0
!        Success=0
!        Failure=0
!!No vertices chosen to be removed yet
!        NoVerts=0
!        do while ((NoVerts.lt.NoMoveDets).and.(k.le.Tries))
!            
!            k=k+1
!            r=RAN2(Seed)
!!Set i=1, since we do not want to choose the first determinant of the original graph
!            i=1
!
!!Look through the normalised inverse eigenvector*eigenvalue (not squared now)
!            do while ((r.gt.0.0_dp).and.(i.lt.NDets))
!                i=i+1
!                r=r-(Eigenvector(i))
!            enddo
!            IF(r.gt.0.0_dp) STOP 'Problem in counting in MoveDets'
!            
!            IF(GraphDets(i,1).eq.0) THEN
!!Graph previously been selected and accepted
!                Remove=.false.
!            ELSE
!                Remove=.true.
!                do j=1,NEl
!                    AttemptDet(j)=GraphDets(i,j)
!                enddo
!            ENDIF
!
!!Need to test whether graph is allowed to be removed - if not, they cycle around for another attempt at 
    !finding a valid determinant
!!            Remove=.true.
!!            DO j=1,NoVerts
!!!Chosen determinant cannot already have been chosen, and must ensure that the other determinants 
        !in the graph are still attached.
!!!If determinant has been chosen before, then its GraphDets-> 0 , so need another test
!!                IF(AttemptDet(1).eq.0) THEN
!!                    Remove=.false.
!!                    EXIT
!!                ENDIF
!!                IF(SameDet(AttemptDet(:),MoveDetsFromPaths(j,:),NEl)) THEN
!!                    Remove=.false.
!!                    EXIT
!!                ENDIF
!!            ENDDO
!
!            IF(Remove) THEN
!!Need to look through the rho matrix to check that there is still a connection for the other determinants in the graph
!                do j=1,NDets
!!If we are at a diagonal element, ignore as it is not a connection
!                    IF(j.eq.i) CYCLE
!                    IF(CopyRhoMat(i,j).gt.0.0_dp) THEN
!!Find that the determinant we are trying to move (i) is already connected to a different determinant (j)
!!See if there are any other connections to determinant j
!                        Remove=.false.
!                        do l=1,NDets
!!Ignore its diagonal element, and the connection to i
!                            IF((l.eq.j).or.(l.eq.i)) CYCLE
!                            IF(ABS(CopyRhoMat(l,j)).gt.0.0_dp) THEN
!!We have found that j is connected to a different determinant, so it is ok to remove i
!                                Remove=.true.
!                                EXIT
!                            ENDIF
!                        enddo
!                    ENDIF
!!The determinant we want to remove would leave a disconnected graph - we cannot use it
!                    IF(.not.Remove) EXIT
!                enddo
!            ENDIF
!
!            IF(Remove) THEN
!!We are able to successfully remove the chosen determinant
!                Success=Success+1
!                NoVerts=NoVerts+1
!                do j=1,NEl
!                    MoveDetsFromPaths(NoVerts,j)=AttemptDet(j)
!                    GraphDets(i,j)=0
!                enddo
!!Store index of determinant to move
!                IndexofDetsFrom(NoVerts)=i
!!Remove Determinant from rho matrix
!                do j=1,NDets
!                    CopyRhoMat(j,i)=0.0_dp
!                    CopyRhoMat(i,j)=0.0_dp
!                enddo
!                IF(.NOT.THDiag) HamElems(i)=(0.0_dp) 
!!Calculate the change to the MeanExcits value...
!                IC=iGetExcitLevel(FDet(:),AttemptDet(:),NEl)
!                MeanExcit=MeanExcit-IC       
!                IF(TDistrib) THEN
!                    Distribs(IC,Iteration+1)=Distribs(IC,Iteration+1)-1
!                ENDIF
!            ELSE
!                Failure=Failure+1
!            ENDIF
!
!!Cycle and try to find another determinant to remove
!        enddo
!
!!MC has failed - print out debugging info
!        IF(k.gt.Tries) THEN
!            WRITE(6,*) 'Error in trying to pick a determinant to move'
!            WRITE(6,*) NoVerts, " removed from previous graph"
!            WRITE(6,*) "Vertices which have been removed are: "
!            do j=1,NoVerts
!                WRITE(6,*) IndexofDetsFrom(j)
!            enddo
!            WRITE(6,*) Tries, " attempted moves made"
!            WRITE(6,*) "Determinants in graph selected with probability: "
!            do j=2,NDets
!                WRITE(6,*) Eigenvector(j)
!            enddo
!            CALL neci_flush(6)
!            STOP 'Error in trying to pick a determinant to move'
!        ELSE
!!            WRITE(6,*) "Determinants which have been removed are: "
!!            do j=1,NoMoveDets
!!                WRITE(6,*) IndexofDetsFrom(j)
!!            enddo
!        ENDIF
!!        WRITE(6,*) "**************"
!!        do j=1,NoMoveDets
!!            WRITE(6,*) MoveDetsFromPaths(j,:)
!!        enddo
!!        CALL neci_flush(6)
!
!        
!!Now need to find a determinant from the excitations space to attach
!        NoVerts=0
!!NoVerts is the number of vertices already attached
!        do while ((NoVerts.lt.NoMoveDets).and.(k.le.Tries))
!            
!            k=k+1
!            r=RAN2(Seed)
!            i=0
!
!!Look through the normalised ExcitsVector array 
!            do while ((r.gt.0.0_dp).and.(i.lt.TotExcits))
!                i=i+1
!                r=r-((ExcitsVector(i))**2)
!            enddo
!            IF(r.gt.0.0_dp) STOP 'Problem in counting in MoveDets 2'
!            
!            do j=1,NEl
!                AttemptDet(j)=ExcitsDets(i,j)
!            enddo
!
!!First we want to check whether the determinant we have selected has already been selected and is already in the graph
!            Attach=.true.
!            do j=1,NoVerts
!                IF(SameDet(AttemptDet(:),GraphDets(IndexofDetsFrom(j),:),NEl)) THEN 
!                    Attach=.false.
!                    EXIT
!                ENDIF
!            enddo
!!Allow the determinant to be one we have already got rid of; this is possible, since a determinant 
!can be specified more than once
!            IF(.not.Attach) THEN
!                do l=(1+NoVerts),NoMoveDets
!                    IF(j.eq.IndexofDetsFrom(l)) THEN
!                        Attach=.true.
!                        EXIT
!                    ENDIF
!                enddo
!            ENDIF
!
!            IF(Attach) THEN
!!Before allowing it to be attached, we must check if the determinant chosen is attached to the graph at all
!                Attach=.false.
!                IF(TNoSameExcit) IC1=iGetExcitLevel(FDet(:),AttemptDet(:),NEl)
!
!         loop2: do j=1,NDets
!                    do l=(1+NoVerts),NoMoveDets
!!Check we are not trying to attach to a determinant which we have chosen to remove! - (Unless it has already 
!been replaced by a new det)
!                        IF(j.eq.IndexofDetsFrom(l)) CYCLE loop2
!                    enddo
!
!                    IF(TNoSameExcit) THEN
!!IC2 information for vertices already added to graph stored in GraphExcitLevel
!                        IC2=GraphExcitLevel(j)
!!                        IC2=iGetExcitLevel(FDet(:),GraphDets(j,:),NEl)
!!If TNoSameDet is on, then do not allow connections between determinants from the same excitation level
!                        IF(IC2.ne.IC1) THEN
!                            IC=iGetExcitLevel(AttemptDet(:),GraphDets(j,:),NEl)
!                            IF(THDiag) THEN
!                                rh=GetHElement2(AttemptDet(:),GraphDets(j,:),NEl,nBasisMax,G1,nBasis,Brr,nMsh,
                                                    !fck,nMax,ALat,UMat,IC,ECore)
!                            ELSE
!                                CALL CalcRho2(AttemptDet(:),GraphDets(j,:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,
                                        !nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
!                            ENDIF
!                            IF(abs(rh).gt.0.0_dp) THEN
!                                Attach=.true.
!                                CopyRhoMat(IndexofDetsFrom(NoVerts+1),j)=rh
!                                CopyRhoMat(j,IndexofDetsFrom(NoVerts+1))=rh
!                            ENDIF
!                        ENDIF
!                    ELSE
!!Allow all possible connections
!                        IC=iGetExcitLevel(AttemptDet(:),GraphDets(j,:),NEl)
!                        IF(IC.eq.0) THEN
!!We have selected a determinant which is already in the graph
!                            Attach=.false.
!                            EXIT loop2
!                        ENDIF
!                        IF(THDiag) THEN
!                            rh=GetHElement2(AttemptDet(:),GraphDets(j,:),NEl,nBasisMax,G1,nBasis,Brr,nMsh,
                                                !fck,nMax,ALat,UMat,IC,ECore)
!                        ELSE
!                            CALL CalcRho2(AttemptDet(:),GraphDets(j,:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,
                                            !nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
!                        ENDIF
!                        IF(abs(rh).gt.0.0_dp) THEN
!                            Attach=.true.
!                            CopyRhoMat(IndexofDetsFrom(NoVerts+1),j)=rh
!                            CopyRhoMat(j,IndexofDetsFrom(NoVerts+1))=rh
!                        ENDIF
!                    ENDIF
!
!                enddo loop2
!            ENDIF
!
!!Attach the determinant
!            IF(Attach) THEN
!                Success=Success+1
!                NoVerts=NoVerts+1
!                do j=1,NEl
!                    GraphDets(IndexofDetsFrom(NoVerts),j)=AttemptDet(j)
!                enddo
!                IF(.not.TNoSameExcit) THEN
!                    IC1=iGetExcitLevel(FDet(:),AttemptDet(:),NEl)
!                ENDIF
!!Find Diagonal rho matrix element, and hamiltonian element
!                IF(THDiag) THEN
!                    rh=GetHElement2(AttemptDet(:),AttemptDet(:),NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,0,ECore)
!                ELSE
!                    HamElems(IndexofDetsFrom(NoVerts))=GetHElement2(FDet,AttemptDet,NEl,nBasisMax,G1,
                                !nBasis,Brr,nMsh,fck,nMax,ALat,UMat,IC1,ECore)
!                    CALL CalcRho2(AttemptDet(:),AttemptDet(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,
                                    !fck,Arr,ALat,UMat,rh,nTay,0,ECore)
!                ENDIF
!                CopyRhoMat(IndexofDetsFrom(NoVerts),IndexofDetsFrom(NoVerts))=rh
!!Calculate the change to the MeanExcits value...
!                MeanExcit=MeanExcit+IC1          
!!Store the excitation level of the new determinant we are adding
!                IF(TNoSameExcit) GraphExcitLevel(IndexofDetsFrom(NoVerts+1))=IC1
!                IF(TDistrib) THEN
!                    Distribs(IC1,Iteration+1)=Distribs(IC1,Iteration+1)+1
!                ENDIF
!
!            ELSE
!                Failure=Failure+1
!            ENDIF
!!Try to attach another determinant
!        enddo
!        
!        IF(k.gt.Tries) THEN
!            WRITE(6,*) 'Error in trying to pick a determinant to create'
!            STOP 'Error in trying to pick a determinant to create'
!        ENDIF
!
!!Move RhoMatrix into the full one
!        ALLOCATE(GraphRhoMat(NDets,NDets),stat=ierr)
!        CALL LogMemAlloc('GraphRhoMat',NDets*NDets,8*HElement_t_size,this_routine,GraphRhoMatTag,ierr)
!        do i=1,NDets
!            do j=1,i
!                GraphRhoMat(i,j)=(CopyRhoMat(i,j))
!                GraphRhoMat(j,i)=(CopyRhoMat(j,i))
!            enddo
!        enddo
!
!!Find the final value for MeanExcits, and the success ratio of the algorithm
!        MeanExcit=MeanExcit/(NDets-1)
!        SucRat=(Success+0.0_dp)/(Success+Failure+0.0_dp)
!        Orig_Graph=NoMoveDets/(NDets+0.0_dp)
!
!        DEALLOCATE(CopyRhoMat)
!        CALL LogMemDealloc(this_routine,CopyRhoMatTag)
!        DEALLOCATE(MoveDetsFromPaths)
!        CALL LogMemDealloc(this_routine,MoveDetsFromPathsTag)
!!In the future, these may not have to be completly changed, only modified slightly
!        DEALLOCATE(ExcitsVector)
!        CALL LogMemDealloc(this_routine,ExcitsVectorTag)
!        DEALLOCATE(ExcitsDets)
!        CALL LogMemDealloc(this_routine,ExcitsDetsTag)
!        
!        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in MoveDetsGraph'
!        call halt_timer(proc_timerMove)
!
!    END SUBROUTINE MoveDetsGraph
!
!!This routine stocastically picks NDets new determinants stochastically from the vector obtained.
!    SUBROUTINE PickNewDets()
!        use SystemData , only : G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,nMsh,Arr
!        use CalcData , only : i_P,RhoEps
!        use IntegralsData , only : fck,nMax,UMat,nTay
!        USE Determinants , only : GetHElement2
!        IMPLICIT NONE
!        INTEGER :: ierr,i,j,k,NoVerts,AttemptDet(NEl),GrowGraphTag,IC,dist,Rej_SameDet
!        INTEGER :: Success,Failure,OrigDets,ExcitDets,Tries,iGetExcitLevel,IC1,IC2,IndexRemoved
!        INTEGER , ALLOCATABLE :: GrowGraph(:,:)
!        real(dp) :: r,Ran2,Sumdetsprob,RootofNum,NormFactor
!        LOGICAL :: Attach,OriginalPicked,SameDet
!        HElement_t :: rh
!        CHARACTER(len=*), PARAMETER :: this_routine='PickNewDets'
!
!        GrowGraphTag=0
!        IF(TNoSameExcit.or.TOneExcitConn) GraphExcitLevel(1:NDets)=0
!
!!Allocate Rho Matrix for growing graph
!        ALLOCATE(GraphRhoMat(NDets,NDets),stat=ierr)
!        CALL LogMemAlloc('GraphRhoMat',NDets*NDets,8*HElement_t_size,this_routine,GraphRhoMatTag,ierr)
!        GraphRhoMat=(0.0_dp)
!        IF(.NOT.THDiag) HamElems=(0.0_dp)
!
!!Initial determinant has to be HF determinant - ensure that this is stored in GraphDets(:,:)
!        ALLOCATE(GrowGraph(NDets,NEl),stat=ierr)
!        CALL LogMemAlloc('GrowGraph',NDets*NEl,4,this_routine,GrowGraphTag,ierr)
!        GrowGraph(1:NEl,1:NDets)=0
!        do i=1,NEl
!            GrowGraph(1,i)=FDet(i)
!        enddo
!        IF(THDiag) THEN
!            GraphRhoMat(1,1)=Hii
!        ELSE
!            GraphRhoMat(1,1)=rhii
!            HamElems(1)=Hii
!        ENDIF
!
!!NoVerts indicates the number of vertices in the graph so far
!        NoVerts=1
!
!!Sumdetsprob is the sum of the probabilities of the determinants in the graph already
!        Sumdetsprob=1.0_dp
!
!!For interesting statistics, keep a record of successful/unsuccessful attachments, and 
!!mean number of excitations away from HF
!        Success=0
!        Failure=0
!        Rej_SameDet=0
!        MeanExcit=0.0_dp
!!Also keep a ratio of determinants which were in the graph before, and new ones added.
!        OrigDets=0
!        ExcitDets=0
!!Allow a maximum average of 5000 attempts for every successfully attached determinant
!        Tries=NDets*5000
!        k=0
!
!!Continue trying to build graph until fully constructed
!        do while ((NoVerts.lt.NDets).and.(k.le.Tries))
!            
!            k=k+1
!!Stochastically pick determinant to add to the graph from vectors
!!Multiply the random number by 1-the sum of the probabilities of the determinants which are already in the graph
!!Since these probabilities have now been set to zero in the vectors, they will not be chosen again.
!            r=RAN2(Seed)*Sumdetsprob
!!            NormFactor=RootofNum(Sumdetsprob,GrowGraphsExpo)
!!            WRITE(6,*) Sumdetsprob,NormFactor
!!Set i=1, since we do not want to choose the first determinant of the original graph
!            i=1
!
!!First look through the normalised eigenvector*eigenvalue
!            do while ((r.gt.0.0_dp).and.(i.lt.NDets))
!                i=i+1
!!                r=r-(((Eigenvector(i))/NormFactor)**(GrowGraphsExpo))
!                r=r-((Eigenvector(i))**(GrowGraphsExpo))
!            enddo
!
!!If still not found, then look in Vector of excitations
!            IF(r.gt.0.0_dp) THEN
!                i=0
!                do while ((r.ge.0.0_dp).and.(i.lt.TotExcits))
!                    i=i+1
!!                    r=r-(((ExcitsVector(i))/NormFactor)**(GrowGraphsExpo))
!                    r=r-((ExcitsVector(i))**(GrowGraphsExpo))
!                enddo
!
!!Error - cannot find determinant to attach in original graph, or its excitations
!                IF(r.gt.0.0_dp) THEN
!                    WRITE(6,*) "k = ", k
!                    WRITE(6,*) "r = ", r
!                    STOP 'Error in stochastic sampling in PickNewDets'
!                ENDIF
!
!!Determinant picked is from excitations...
!                OriginalPicked=.false.
!                IndexRemoved=i
!                do j=1,NEl
!                    AttemptDet(j)=ExcitsDets(i,j)
!                enddo
!
!            ELSE
!!Determinant picked is from original graph...
!                OriginalPicked=.true.
!                IndexRemoved=i
!                do j=1,NEl
!                    AttemptDet(j)=GraphDets(i,j)
!                enddo
!
!            ENDIF
!
!!Need to test whether graph is accepted - if not, they cycle around for another attempt at finding a valid determinant
!            Attach=.false.
!!IC1 is the excitation level of the attempted determinant to attach
!            IC1=iGetExcitLevel(FDet,AttemptDet,NEl)
!            DO i=1,NoVerts
!
!!Check the number of excitations away from each other vertex in graph - or instead, just check if determinant already in graph
!                IF(SameDet(AttemptDet(:),GrowGraph(i,:),NEl)) THEN
!!Determinant is already in the graph - exit loop - determinant not valid
!                    Rej_SameDet=Rej_SameDet+1
!
!!The only way that the same determinant can be picked again is if it is multiply specified and already included from elsewhere
!!Therefore, we can remove this reference from the list in the same way
!                    IF(OriginalPicked) THEN
!                        Sumdetsprob=Sumdetsprob-((Eigenvector(IndexRemoved))**(GrowGraphsExpo))
!                        Eigenvector(IndexRemoved)=(0.0_dp)
!                    ELSE
!                        Sumdetsprob=Sumdetsprob-((ExcitsVector(IndexRemoved))**(GrowGraphsExpo))
!                        ExcitsVector(IndexRemoved)=(0.0_dp)
!                    ENDIF
!
!!Do not attach determinant, and no need to carry on testing connections                    
!                    Attach=.false.
!                    EXIT
!                ENDIF
!
!!Determine connectivity to other determinants in the graph, if no connection yet found
!                IF(.not.Attach) THEN
!                    IF(TNoSameExcit) THEN
!!excitation level of determinants already in the graph stored in GraphExcitLevel
!!                        IC2=iGetExcitLevel(FDet(:),GrowGraph(i,:),NEl)
!                        IC2=GraphExcitLevel(i)
!!Don't allow connections to the same excitation level
!                        IF(IC1.ne.IC2) THEN
!                            IC=iGetExcitLevel(GrowGraph(i,:),AttemptDet(:),NEl)
!                            IF(THDiag) THEN
!                                rh=GetHElement2(AttemptDet(:),GrowGraph(i,:),NEl,nBasisMax,G1,nBasis,Brr,
                                            !nMsh,fck,nMax,ALat,UMat,IC,ECore)
!                                IF(abs(rh).gt.0.0_dp) Attach=.true.
!                            ELSE
!                                CALL CalcRho2(AttemptDet(:),GrowGraph(i,:),Beta,i_P,NEl,nBasisMax,G1,nBasis,
                                        !Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
!                                IF(abs(rh).gt.RhoEps) Attach=.true. 
!                            ENDIF
!                        ENDIF
!
!                    ELSEIF(TOneExcitConn) THEN
!!Only allow connections between excitation levels which differ by one
!                        IC2=GraphExcitLevel(i)
!                        IF((ABS(IC1-IC2)).eq.1) THEN
!                            IC=iGetExcitLevel(GrowGraph(i,:),AttemptDet(:),NEl)
!                            IF(THDiag) THEN
!                                rh=GetHElement2(AttemptDet(:),GrowGraph(i,:),NEl,nBasisMax,G1,nBasis,Brr,
                                        !nMsh,fck,nMax,ALat,UMat,IC,ECore)
!                                IF(abs(rh).gt.0.0_dp) Attach=.true.
!                            ELSE
!                                CALL CalcRho2(AttemptDet(:),GrowGraph(i,:),Beta,i_P,NEl,nBasisMax,G1,nBasis,
                                            !Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
!                                IF(abs(rh).gt.RhoEps) Attach=.true. 
!                            ENDIF
!                        ENDIF
!
!                    ELSE
!!Allow all connections
!                        IF(THDiag) THEN
!                            rh=GetHElement2(AttemptDet(:),GrowGraph(i,:),NEl,nBasisMax,G1,nBasis,Brr,
                                        !nMsh,fck,nMax,ALat,UMat,-1,ECore)
!                            IF(abs(rh).gt.0.0_dp) Attach=.true.
!                        ELSE
!                            CALL CalcRho2(AttemptDet(:),GrowGraph(i,:),Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,
                                        !nMsh,fck,Arr,ALat,UMat,rh,nTay,-1,ECore)
!                            IF(abs(rh).gt.RhoEps) Attach=.true. 
!                        ENDIF
!                    ENDIF
!                ENDIF
!
!            ENDDO
!
!            IF(Attach) THEN
!!Attach new determinant to the list as it passes tests
!
!!Make the probability of picking the determinants again = 0.
!!Need to first sum the probability so that the random number can be adjusted accordingly
!                IF(OriginalPicked) THEN
!                    Sumdetsprob=Sumdetsprob-((Eigenvector(IndexRemoved))**(GrowGraphsExpo))
!                    Eigenvector(IndexRemoved)=(0.0_dp)
!                ELSE
!                    Sumdetsprob=Sumdetsprob-((ExcitsVector(IndexRemoved))**(GrowGraphsExpo))
!                    ExcitsVector(IndexRemoved)=(0.0_dp)
!                ENDIF
!
!!First check on attachment to HF, since wants to be stored in HamElems
!!                IF(IC1.gt.2) WRITE(6,*) "Higher excitation attached ", IC1
!     
!                MeanExcit=MeanExcit+IC1
!                IF(TDistrib) THEN
!                    Distribs(IC1,Iteration+1)=Distribs(IC1,Iteration+1)+1
!                ENDIF
!
!!Only double excitations (or single? <- include for completness) of HF contribute
!                IF(IC1.eq.2) THEN!.or.(IC1.eq.1)) THEN
!                    IF(THDiag) THEN
!                        rh=GetHElement2(FDet,AttemptDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,IC1,ECore)
!                    ELSE
!                        HamElems(NoVerts+1)=GetHElement2(FDet,AttemptDet,NEl,nBasisMax,G1,nBasis,Brr,
                                        !nMsh,fck,nMax,ALat,UMat,IC1,ECore)
!                        CALL CalcRho2(FDet,AttemptDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,
                                            !UMat,rh,nTay,IC1,ECore)
!                    ENDIF
!!Add rhoij contribution. No real need to set a rho epsilon here? We are diagonalising anyway...
!                    GraphRhoMat(1,NoVerts+1)=rh
!                    GraphRhoMat(NoVerts+1,1)=rh
!
!                ELSE
!!Added determinant not connected to root
!                    IF(.NOT.THDiag) HamElems(NoVerts+1)=(0.0_dp)
!                    GraphRhoMat(1,NoVerts+1)=(0.0_dp)
!                    GraphRhoMat(NoVerts+1,1)=(0.0_dp)
!                ENDIF
!
!                IF(TNoSameExcit) THEN
!!Don't allow connections to the same excitation level
!                    do i=2,NoVerts
!!Excitation level information stored in GraphExcitLevel
!!                        IC2=iGetExcitLevel(FDet(:),GrowGraph(i,:),NEl)
!                        IC2=GraphExcitLevel(i)
!                        IF(IC1.ne.IC2) THEN
!                            IC=iGetExcitLevel(AttemptDet(:),GrowGraph(i,:),NEl)
!                            IF(THDiag) THEN
!                                rh=GetHElement2(GrowGraph(i,:),AttemptDet(:),NEl,nBasisMax,G1,nBasis,Brr,
                                            !nMsh,fck,nMax,ALat,UMat,IC,ECore)
!                            ELSE
!                                CALL CalcRho2(GrowGraph(i,:),AttemptDet(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,
                                        !Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
!                            ENDIF
!                            GraphRhoMat(i,NoVerts+1)=rh
!                            GraphRhoMat(NoVerts+1,i)=rh
!                        ENDIF
!                    enddo
!!Store the excitation level of the attached determinant
!                    GraphExcitLevel(NoVerts+1)=IC1
!
!                ELSEIF(TOneExcitConn) THEN
!!Only allow connections between excitations which differ by one excitation level from HF
!                    do i=2,NoVerts
!!Excitation level information stored in GraphExcitLevel
!                        IC2=GraphExcitLevel(i)
!                        IF((ABS(IC1-IC2)).eq.1) THEN
!                            IC=iGetExcitLevel(AttemptDet(:),GrowGraph(i,:),NEl)
!                            IF(THDiag) THEN
!                                rh=GetHElement2(GrowGraph(i,:),AttemptDet(:),NEl,nBasisMax,G1,nBasis,Brr,
                                        !nMsh,fck,nMax,ALat,UMat,IC,ECore)
!                            ELSE
!                                CALL CalcRho2(GrowGraph(i,:),AttemptDet(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,
                                            !Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,IC,ECore)
!                            ENDIF
!                            GraphRhoMat(i,NoVerts+1)=rh
!                            GraphRhoMat(NoVerts+1,i)=rh
!                        ENDIF
!                    enddo
!!Store the excitation level of the attached determinant
!                    GraphExcitLevel(NoVerts+1)=IC1
!
!                ELSE
!
!!Run through rest of excitations in the graph testing contributions and adding to rho matrix
!                    do i=2,NoVerts
!                        IF(THDiag) THEN
!                            rh=GetHElement2(GrowGraph(i,:),AttemptDet(:),NEl,nBasisMax,G1,nBasis,Brr,
                                        !nMsh,fck,nMax,ALat,UMat,-1,ECore)
!                        ELSE
!                            CALL CalcRho2(GrowGraph(i,:),AttemptDet(:),Beta,i_P,NEl,nBasisMax,G1,nBasis,
                                        !Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,-1,ECore)
!                        ENDIF
!                        GraphRhoMat(i,NoVerts+1)=rh
!                        GraphRhoMat(NoVerts+1,i)=rh
!                    enddo
!
!                ENDIF
!
!!Include diagonal rho matrix element for chosen determinant
!                IF(THDiag) THEN
!                    rh=GetHElement2(AttemptDet,AttemptDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,0,ECore)
!                ELSE
!                    CALL CalcRho2(AttemptDet,AttemptDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,
                                            !ALat,UMat,rh,nTay,0,ECore)
!                ENDIF
!                GraphRhoMat(NoVerts+1,NoVerts+1)=rh
!
!!Add determinant to growing graph
!                do i=1,NEl
!                    GrowGraph(NoVerts+1,i)=AttemptDet(i)
!                enddo
!                NoVerts=NoVerts+1
!!                WRITE(6,"A,I5") "Vertex Added - ",NoVerts
!!                CALL neci_flush(6)
!
!                Success=Success+1
!                IF(OriginalPicked) THEN
!                    OrigDets=OrigDets+1
!                ELSE
!                    ExcitDets=ExcitDets+1
!                ENDIF
!
!            ELSE
!!Determinant chosen not attached to growing graph
!                Failure=Failure+1
!
!            ENDIF
!
!!Attempt to find another determinant to attach
!        enddo
!
!!Find average number of excitations away from HF determinant
!        MeanExcit=MeanExcit/(NDets-1)
!
!!Test if graph growing was successful. If not, possible reasons are that the space is too small for
!!such a large graph to be grown easily - reduce NDets, or increase time to search for determinants.
!!Also, check that Tries hasn't gone above the range allowed for an integer
!        IF(NoVerts.ne.NDets) THEN
!            WRITE(6,*) "Error in attaching determinants to new graph"
!            WRITE(6,*) "Out of ",Tries," attempts to attach determinants, ",Failure, " failed."
!            WRITE(6,*) "Of these failures, ", Rej_SameDet," of these were due to the same determinant being picked"
!            WRITE(6,*) "The rest were due to the determinants not being connected to the growing graph"
!            WRITE(6,*) "Eigenvector components are: "
!            do k=2,NDets
!                WRITE(6,*) Eigenvector(k)
!            enddo
!            STOP 'Error in attaching determinants to new graph'
!        ENDIF
!        IF(Tries.gt.NDets*5000) STOP 'Should never get here'
!        IF(Success.ne.(NDets-1)) STOP 'Error in attaching determinants to new graph 2'
!
!!Determine success ratio for attachment of determinants into new graph
!        SucRat=(Success+0.0_dp)/(Success+Failure+0.0_dp)
!        Orig_Graph=((OrigDets+0.0_dp)/(ExcitDets+OrigDets+0.0_dp))*100
!        WRITE(6,"(A,F9.5)") "New graph created. Success ratio for adding determinants: ",SucRat
!        WRITE(6,"(A,F9.4)") "Percentage of determinants found in original graph: ",Orig_Graph
!
!!Replace original list of determinants in graph with new list
!        GraphDets(1:NDets,1:NEl)=0
!        do i=1,NDets
!            do j=1,NEl
!                IF(GrowGraph(i,j).eq.0) STOP 'Error in allocating determinants to graph'
!                GraphDets(i,j)=GrowGraph(i,j)
!            enddo
!        enddo
!
!!Deallocate temporary list of determinants in graph, list of all excitations, and the ExcitsVector vector
!        DEALLOCATE(ExcitsVector)
!        CALL LogMemDealloc(this_routine,ExcitsVectorTag)
!        DEALLOCATE(ExcitsDets)
!        CALL LogMemDealloc(this_routine,ExcitsDetsTag)
!        DEALLOCATE(GrowGraph)
!        CALL LogMemDealloc(this_routine,GrowGraphTag)
!
!        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in PickNewDets'
!
!    END SUBROUTINE PickNewDets
!
!!This is used to create normalised vectors for the MoveDets stochastic graph-growing algorithm
!    SUBROUTINE NormaliseVectorSep()
!        IMPLICIT NONE
!        HElement_t :: Norm1,Norm2
!        real(dp) :: RootofNum
!        INTEGER :: i
!        type(timer), save :: proc_timerNorm 
!        
!        proc_timerNorm%timer_name='NormVecSep'
!        call set_timer(proc_timerNorm)
!        
!!Setup a normalised Inverse vector of determinants in the graph, so they can be chosen stochastically
!!Since we no longer need the largest eigenvector, we can multiply its elements by its eigenvalue, then
!!use it as the inverse vector
!        do i=1,NDets
!            Eigenvector(i)=Eigenvector(i)*(Eigenvalue)
!        enddo
!
!        Norm1=(0.0_dp)
!
!        do i=2,NDets
!!            IF(ABS(Eigenvector(i)).lt.1.0e-17_dp) THEN
!!                WRITE(6,*) "For determinant ", i, ", connection is ", Eigenvector(i)
!!                STOP 'Numerical errors likely to arise due to such small connection'
!!            ENDIF
!!Normalised as the reciprocal of the fourth root of the element. If root is changed, then 
    !numerical stability means that only one determinant is ever picked...
!            Eigenvector(i)=(1.0_dp/RootofNum(ABS(Eigenvector(i)),5.0_dp))
!!            WRITE(6,*) Eigenvector(i)
!            Norm1=Norm1+Eigenvector(i)
!        enddo
!
!        do i=2,NDets
!            Eigenvector(i)=Eigenvector(i)/Norm1
!!            WRITE(6,*) Eigenvector(i)
!        enddo
!
!!The excitations need to be normalised separatly, according to their amplitude squared
!
!        Norm2=(0.0_dp)
!        do i=1,TotExcits
!            Norm2=Norm2+(ExcitsVector(i)*ExcitsVector(i))
!        enddo
!        Norm2=(SQRT(Norm2))
!        do i=1,TotExcits
!            ExcitsVector(i)=ExcitsVector(i)/Norm2
!        enddo
!
!        WRITE(6,"(A,F12.5)") "Normalisation constant for the determinants in the graph is ",Norm1
!
!        PStay=NoMoveDets/NDets
!
!        call halt_timer(proc_timerNorm)
!
!    END SUBROUTINE NormaliseVectorSep
!
!!This is used to normalise the vector of current and excited determinants from which we are going to pick the next graph.
!!The vector is spread between the arrays for the original eigenvector (which needs to be multiplied by corresponding eigenvalue)
!!and the ExcitsVector. The first element of the eigenvector should not be included in the normalisation, as
!!it cannot be picked - the HF is always in each graph.
!    SUBROUTINE NormaliseVector()
!        IMPLICIT NONE
!        INTEGER :: i
!        real(dp) :: Stay,Move,RootofNum
!        real(dp) :: Norm,Norm1,Norm2
!
!!The bias towards the determinants already in the graph is given by the largest eigenvector, multiplied by its eigenvalue.
!!Since we no longer need the largest eigenvector, we can multiply its elements by its eigenvalue
!        do i=2,NDets
!            Eigenvector(i)=Eigenvector(i)*(Eigenvalue)
!        enddo
!
!!An addititional bias against the determinants already in the graph can be designed.
!        IF(TBiasing) THEN
!!GraphBias is the probability of picking a determinant which is already in the graph.
!!Remember that this is not actually a strict probability, since the way that the graph is grown,
!!means that some determinants will be automatically preferred over others.
!            
!            IF((GraphBias.gt.1.0_dp).or.(GraphBias.lt.0.0_dp)) THEN
!                STOP 'Value for Graphbias must be between 1 and 0'
!            ENDIF
!
!            Norm1=0.0_dp
!!The elements of the vectors can be turned into probabilities by raising to powers other than two
!            do i=2,NDets
!                Norm1=Norm1+(ABS(Eigenvector(i))**(GrowGraphsExpo))
!            enddo
!            Norm1=RootofNum((Norm1/GraphBias),GrowGraphsExpo)
!!            Norm1=SQRT(Norm1/GraphBias)
!
!!Divide elements of the eigenvector by the normalisation
!            Stay=0.0_dp
!            do i=2,NDets
!                Eigenvector(i)=(ABS(Eigenvector(i))/Norm1)
!                Stay=Stay+(ABS(Eigenvector(i))**(GrowGraphsExpo))
!            enddo
!
!!Now find normalisation for the excitations
!            Norm2=0.0_dp
!            do i=1,TotExcits
!                Norm2=Norm2+(ABS(ExcitsVector(i))**(GrowGraphsExpo))
!            enddo
!!            Norm2=SQRT(Norm2/(1.0_dp-GraphBias))
!            Norm2=RootofNum((Norm2/(1.0_dp-GraphBias)),GrowGraphsExpo)
!
!!Divide elements of ExcitsVector by new normalisation
!!            Move=0.0_dp
!            do i=1,TotExcits
!                ExcitsVector(i)=(ABS(ExcitsVector(i))/Norm2)
!!                IF(GraphBias.eq.1.0_dp) ExcitsVector(i)=(0.0_dp)
!!                Move=Move+ABS((ExcitsVector(i))**(GrowGraphsExpo))
!            enddo
!
!        ELSE
!!No biasing towards excitations
!
!!We need to find the normalisation constant, given by the sum of the squares of all the elements
!            Norm=0.0_dp
!!First, sum the squares of the original determinants in the graph
!            do i=2,NDets
!                Norm=Norm+(ABS(Eigenvector(i))**(GrowGraphsExpo))
!            enddo
!
!!Then, sum the squares of the vector for the excitations
!            do i=1,TotExcits
!                Norm=Norm+(ABS(ExcitsVector(i))**(GrowGraphsExpo))
!            enddo
!
!            Norm=RootofNum(Norm,GrowGraphsExpo)
!
!!Once the normalisation is found, all elements need to be divided by it.
!!Stay is the total probability of staying with original graph
!            Stay=0.0_dp
!!            Move=0.0_dp
!            do i=2,NDets
!                Eigenvector(i)=(ABS(Eigenvector(i))/Norm)
!                Stay=Stay+ABS((Eigenvector(i))**(GrowGraphsExpo))
!            enddo
!            do i=1,TotExcits
!                ExcitsVector(i)=(ABS(ExcitsVector(i))/Norm)
!!                Move=Move+ABS((ExcitsVector(i))**(GrowGraphsExpo))
!            enddo
!
!        ENDIF
!
!        PStay=Stay
!        WRITE(6,*) "Probability of staying at original determinants: ", Stay
!!        WRITE(6,*) "Probability of Moving: ", Move
!!        CALL neci_flush(6)
!!        WRITE(6,*) "Total Probability: ", Stay+Move
!!        WRITE(6,*) "Normalisation constant for propagation vector: ", Norm
!
!    END SUBROUTINE NormaliseVector
!
!!This routine creates a new vector (ExcitsVector), which is formed by the matrix product of the rho matrix for the 
!new excitations
!!(ConnectionsToExcits) and the original eigenvector.
!!This is the rho operator for the excitations acting on the original eigenstate to show the probabilities of being in the new
!!excitations. Since each excitation can only be attached to one vertex of the graph, there is only one element per
!!row of the operator, and only one multiplication is needed per excitation.
!    SUBROUTINE CreateExcitsVector()
!        IMPLICIT NONE
!        INTEGER :: ierr,NoExcitsAttached,Element,i,j
!        CHARACTER(len=*), PARAMETER :: this_routine='CreateExcitsVector'
!
!!Allocate space to hold this new vector - could get away without allocating any more memory by simply multiplying Connections
!!array by needed element of eigenvector...?
!        ALLOCATE(ExcitsVector(TotExcits),stat=ierr)
!        CALL LogMemAlloc('ExcitsVector',TotExcits,8*HElement_t_size,this_routine,ExcitsVectorTag,ierr)
!        ExcitsVector=(0.0_dp)
!
!!Separatly deal with first vertex in graph for clarity
!        do j=1,NoExcits(1)
!!Multiply rho element to excitation by eigenvector component corresponding to determinant being excited
!            ExcitsVector(j)=ConnectionsToExcits(j)*Eigenvector(1)
!        enddo
!
!!Cycle through all remaining vertices in graph
!        do i=2,NDets
!
!            NoExcitsAttached=NoExcits(i)-NoExcits(i-1)
!!Cycle through all determinants attached to vertex i
!            do j=1,NoExcitsAttached
!
!!Corresponding element in 'ConnectionsToExcits' and 'ExcitsVector' is calculated
!                Element=j+NoExcits(i-1)
!                ExcitsVector(Element)=ConnectionsToExcits(Element)*Eigenvector(i)
!
!            enddo
!
!        enddo
!
!!Deallocate NoExcits Array and ConnectionsToExcits array
!        DEALLOCATE(ConnectionsToExcits)
!        CALL LogMemDealloc(this_routine,ConnectionsToExcitsTag)
!        DEALLOCATE(NoExcits)
!        CALL LogMemDealloc(this_routine,NoExcitsTag)
!
!        IF(Element.ne.TotExcits) THEN
!            STOP 'Error in counting in CreateExcitsVector'
!        ENDIF
!
!        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in CreateExcitsVector'
!    END SUBROUTINE CreateExcitsVector
!
!!This subroutine will go through all determinants, find all excitations, and then calculate the rho element to
!!the determinant in the graph to which it is connected.
!    SUBROUTINE FindConnections()
!        use SystemData, only: G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,nMsh,Arr
!        use CalcData , only : i_P
!        use IntegralsData, only: fck,nMax,UMat,nTay
!        USE Determinants , only : GetHElement2
!        IMPLICIT NONE
!        HElement_t :: rh
!        real(dp) :: Prob
!        type(timer), save :: proc_timerConns
!        INTEGER :: attempts,NoExcitsCurr,Noatt
!        INTEGER :: ierr,i,j,DetCurr(NEl),nJ(NEl),nStore(6),iMaxExcit,nExcitMemLen
!        INTEGER :: nExcitTag,iExcit,ExcitCurr,dist,iGetExcitLevel,IC,exFlag
!        INTEGER , ALLOCATABLE :: nExcit(:)
!        CHARACTER(len=*), PARAMETER :: this_routine='FindConnections'
!
!        proc_timerConns%timer_name='FindConns'
!        call set_timer(proc_timerConns)
!        
!        IF(TSinglesExcitSpace) THEN
!!Only single excitations of the determinants in the graph are created
!            exFlag=1
!        ELSE
!            exFlag=3
!        ENDIF
!
!        nExcitTag=0
!
!!Allocate memory to hold connections, and form of the determinants of the excitations
!        ALLOCATE(ConnectionsToExcits(TotExcits),stat=ierr)
!        CALL LogMemAlloc('ConnectionsToExcits',TotExcits,8*HElement_t_size,this_routine,ConnectionsToExcitsTag,ierr)
!        ConnectionsToExcits=(0.0_dp)
!        ALLOCATE(ExcitsDets(TotExcits,NEl),stat=ierr)
!        CALL LogMemAlloc('ExcitsDets',TotExcits*NEl,4,this_routine,ExcitsDetsTag,ierr)
!        ExcitsDets(1:TotExcits,1:NEl)=0
!
!!Cycle over all vertices in graph
!        ExcitCurr=0
!        do i=1,NDets
!
!!Find determinant form for current vertex
!            do j=1,NEl
!                DetCurr(j)=GraphDets(i,j)
!            enddo
!
!!Setup excitation generators again for this determinant
!            iMaxExcit=0
!            nStore(1:6)=0
!            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
!            ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
!            CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag,ierr)
!!            nExcit(1:nExcitMemLen)=0
!            nExcit(1)=0
!            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
!
!            IF(TMCExcits) THEN
!!Excitations are picked stocastically
!
!                Attempts=0
!                NoExcitsCurr=0
!
!                do while(NoExcitsCurr.lt.NoMCExcits)
!
!                    IF(Attempts.gt.(NoMCExcits*100)) THEN
!                        WRITE(6,*) "Unable to find enough determinants attached to graph determinant ",DetCurr
!                        WRITE(6,*) "Attempts: ", Attempts
!                        STOP "Unable to find enough determinants attached to graph determinant "
!                    ENDIF
!
!                    CALL GenRandSymExcitIt2(DetCurr,NEl,nExcit,nJ,Seed,Noatt,Prob)
!                    Attempts=Attempts+1
!!Divide the coupling to a determinant by the probability of generating it, to account for possibility of bias in generation
!
!                    IF(TMaxExcit) THEN
!
!!A maximum excitation level is set - don't allow connection if its too high in excitation space
!                        IC=iGetExcitLevel(FDet,nJ,NEl)
!                        IF(IC.le.iMaxExcitLevel) THEN
!!Throw the random excitation, and don't count it if it is above the excitation level threshold
!                            ExcitCurr=ExcitCurr+1
!                            NoExcitsCurr=NoExcitsCurr+1
!                            IF(THDiag) THEN
!                                rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,-1,ECore)
!                            ELSE
!                                CALL CalcRho2(DetCurr,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,
                                        !ALat,UMat,rh,nTay,-1,ECore)
!                            ENDIF
!!Store path of determinant in ExcitsDets, and the rho elements in ConnectionsToExcits
!                            rh=rh/((Prob))
!                            ConnectionsToExcits(ExcitCurr)=rh
!                            do j=1,NEl
!                                ExcitsDets(ExcitCurr,j)=nJ(j)
!                            enddo
!
!                        ENDIF
!
!                    ELSE
!!No restriction on excitation level 
!        
!                        ExcitCurr=ExcitCurr+1
!                        NoExcitsCurr=NoExcitsCurr+1
!                        IF(THDiag) THEN
!                            rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,-1,ECore)
!                        ELSE
!                            CALL CalcRho2(DetCurr,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,
                                        !UMat,rh,nTay,-1,ECore)
!                        ENDIF
!!Store path of determinant in ExcitsDets, and the rho elements in ConnectionsToExcits
!                        rh=rh/((Prob))
!                        ConnectionsToExcits(ExcitCurr)=rh
!                        do j=1,NEl
!                            ExcitsDets(ExcitCurr,j)=nJ(j)
!                        enddo
!
!                    ENDIF
!
!                enddo
!
!            ELSE
!
!                IF(i.eq.1) THEN
!                    IF(iMaxExcit.ne.NoExcits(i)) STOP 'Error in counting in FindConnections'
!                ELSE
!                    IF(iMaxExcit.ne.(NoExcits(i)-NoExcits(i-1))) STOP 'Error in counting in FindConnections'
!                ENDIF
!
!!Cycle through all excitations of each determinant
!            lp: do while(.true.)
!                    CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,iExcit,0,nStore,exFlag)
!                    IF(nJ(1).eq.0) exit lp
!                    ExcitCurr=ExcitCurr+1
!
!                    IF(TMaxExcit) THEN
!!A maximum excitation level for the space of accessible determinants is imposed
!                        IC=iGetExcitLevel(FDet,nJ,NEl)
!                        IF(IC.gt.iMaxExcitLevel) THEN
!!If excitation is further away than we want, then let connection to it = 0
!                            ConnectionsToExcits(ExcitCurr)=(0.0_dp)
!                        ELSE
!                            IF(THDiag) THEN
!                                rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,iExcit,ECore)
!                            ELSE
!                                CALL CalcRho2(DetCurr,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,
                                                !ALat,UMat,rh,nTay,iExcit,ECore)
!                            ENDIF
!!Store path of determinant in ExcitsDets, and the rho elements in ConnectionsToExcits
!                            ConnectionsToExcits(ExcitCurr)=rh
!                            do j=1,NEl
!                                ExcitsDets(ExcitCurr,j)=nJ(j)
!                            enddo
!                        ENDIF
!
!                    ELSE
!!No restriction on excitation level
!                        IF(THDiag) THEN
!                            rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,iExcit,ECore)
!                        ELSE
!                            CALL CalcRho2(DetCurr,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,
                                            !UMat,rh,nTay,iExcit,ECore)
!                        ENDIF
!!                    Dist=IGetExcitLevel(FDet,nJ,NEl)
!!                    IF(Dist.gt.4) THEN
!!                        WRITE(6,*) "Higher than double excitation found - attached to det:", DetCurr
!!                        WRITE(6,*) "Determinant is: ", nJ
!!                        WRITE(6,*) Dist," fold excitation"
!!                    ENDIF
!
!!Store path of determinant in ExcitsDets, and the rho elements in ConnectionsToExcits
!                        ConnectionsToExcits(ExcitCurr)=rh
!                        do j=1,NEl
!                            ExcitsDets(ExcitCurr,j)=nJ(j)
!                        enddo
!                    ENDIF
!                enddo lp
!
!            ENDIF
!
!            IF(ExcitCurr.ne.NoExcits(i)) STOP 'Incorrect counting in FindConnections'
!
!!Deallocate excitation generators
!            DEALLOCATE(nExcit)
!            CALL LogMemDealloc(this_routine,nExcitTag)
!
!        enddo
!
!        IF(ExcitCurr.ne.TotExcits) STOP 'Incorrect counting in FondConnections'
!
!        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in FindConnections'
!        call halt_timer(proc_timerConns)
!
!    END SUBROUTINE FindConnections
!
!!This routine initialises and then destroys excitation generators, in order to calculate the space required for all excitations.
!!This could improve on space if all excitations were explicitly calculated, though would take longer.
!    SUBROUTINE CountExcits()
!
!        use SystemData, only : G1,nBasis,nBasisMax
!        IMPLICIT NONE
!        INTEGER :: ierr,i,j,nStore(6),DetCurr(NEl),nJ(NEl),iMaxExcit,nExcitMemLen
!        INTEGER :: nExcitTag,exFlag
!        CHARACTER(len=*), PARAMETER :: this_routine='CountExcits'
!        INTEGER , ALLOCATABLE :: nExcit(:)
!
!        IF(TMCExcits) THEN
!!This searches for a given number of excitations stocastically per determinant - no need to count them
!            
!            TotExcits=NDets*NoMCExcits
!            ALLOCATE(NoExcits(NDets),stat=ierr)
!            CALL LogMemAlloc('NoExcits',NDets,4,this_routine,NoExcitsTag,ierr)
!            NoExcits(1:NDets)=0
!
!            do i=1,NDets
!                IF(i.eq.1) THEN
!                    NoExcits(i)=NoMCExcits
!                ELSE
!                    NoExcits(i)=NoExcits(i-1)+NoMCExcits
!                ENDIF
!            enddo
!
!        ELSE
!!Search through the entire space of single or single&double excitations of each determinant
!
!            IF(TSinglesExcitSpace) THEN
!!Only search through single excitations
!                exFlag=1
!            ELSE
!                exFlag=3
!            ENDIF
!            
!            TotExcits=0
!            nExcitTag=0
!!Allocate Memory for NoExcits array
!            ALLOCATE(NoExcits(NDets),stat=ierr)
!            CALL LogMemAlloc('NoExcits',NDets,4,this_routine,NoExcitsTag,ierr)
!            NoExcits(1:NDets)=0
!
!!Cycle over all vertices in graph
!            do i=1,NDets
!                
!!Find determinant form for current vertex
!                do j=1,NEl
!                    DetCurr(j)=GraphDets(i,j)
!                enddo
!
!!Create excitation generator
!                iMaxExcit=0
!                nStore(1:6)=0
!                CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
!                ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
!                CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag,ierr)
!                nExcit(1:nExcitMemLen)=0
!                CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
!                
!!Store number of excitations from each determinant cumulativly.
!                IF(i.eq.1) THEN
!                    NoExcits(i)=iMaxExcit
!                ELSE
!                    NoExcits(i)=NoExcits(i-1)+iMaxExcit
!                ENDIF
!
!                TotExcits=TotExcits+iMaxExcit
!
!!Destroy excitation generator
!                DEALLOCATE(nExcit)
!                CALL LogMemDealloc(this_routine,nExcitTag)
!                nStore(1:6)=0
!            
!            enddo
!
!        ENDIF
!
!!Perform check that all excitations accounted for
!        IF(NoExcits(NDets).ne.TotExcits) THEN
!            STOP 'Error in counting excits in CountExcits'
!        ENDIF
!
!        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in CountExcits'
!    END SUBROUTINE CountExcits
!
!!This routine uses a Lancsoz iterative diagonalisation technique to diagonalise the rho matrix.
!    SUBROUTINE DiagGraphLanc()
!        USE DetCalc , only : B2L, nBlk, nKry, NEval
!        IMPLICIT NONE
!        CHARACTER(len=*), PARAMETER :: this_routine='DiagGraphLanc'
!        INTEGER , ALLOCATABLE :: Lab(:),NRow(:),ISCR(:),Index(:)
!        real(dp) , ALLOCATABLE :: Mat(:),CK(:,:),CKN(:,:), temp(:)
!        real(dp) , ALLOCATABLE :: A(:,:),V(:),AM(:),BM(:),T(:),WT(:),SCR(:)
!        real(dp) , ALLOCATABLE :: Work2(:),WH(:),V2(:,:),W(:)
!        INTEGER :: LabTag,ATag,MatTag,NRowTag,VTag,WTTag
!        INTEGER :: AMTag,BMTag,TTag,SCRTag,ISCRTag,IndexTag,WHTag
!        INTEGER :: Work2Tag,V2Tag,WTag,CKTag,CKNTag
!        INTEGER :: ierr,LenMat,i,j,ICMax,RowElems
!        type(timer), save :: proc_timerLanc 
!        real(dp) :: SumVec,LancVar
!        INTEGER :: NCycle,NBlock,NKry1,LScr,LIScr
!        LOGICAL :: TSeeded
!
!        proc_timerLanc%timer_name='DiagGraphLanc'
!        call set_timer(proc_timerLanc)
!        
!!Zero memory tags
!        LabTag=0
!        ATag=0
!        MatTag=0
!        NRowTag=0
!        VTag=0
!        WTTag=0
!        AMTag=0
!        BMTag=0
!        TTag=0
!        SCRTag=0
!        ISCRTag=0
!        IndexTag=0
!        WHTag=0
!        Work2Tag=0
!        V2Tag=0
!        WTag=0
!        CKTag=0
!        CKNTag=0
!
!        IF(TMoveDets) THEN
!!If using the MoveDets algorithm, need to keep the rho matrix for the previous graph
!            ALLOCATE(CopyRhoMat(NDets,NDets),stat=ierr)
!            CALL LogMemAlloc('CopyRhoMat',NDets**2,8,this_routine,CopyRhoMatTag,ierr)
!            do i=1,NDets
!                do j=i,NDets
!                    CopyRhoMat(i,j)=GraphRhoMat(i,j)
!                    CopyRhoMat(j,i)=GraphRhoMat(j,i)
!                enddo
!            enddo
!        ENDIF
!
!
!!This seems to be Alex's way of compressing the matrix, but not the way in the diagonaliser
!!In this method, only the non-zero elements are stored (and only the symmetric matrix)
!        LenMat=0
!        ICMax=0
!        do i=1,NDets
!            RowElems=0
!            do j=i,NDets
!                IF(abs(GraphRhoMat(i,j)).gt.0.0_dp) THEN
!!                    LenMat=LenMat+1
!                    RowElems=RowElems+1
!                ENDIF
!            enddo
!            LenMat=LenMat+RowElems
!            IF(RowElems.gt.ICMax) ICMax=RowElems
!        enddo
!
!!        WRITE(6,"(A,I10,A,I10,A)") "Of ",((NDets*NDets+NDets)/2), " possible elements, only ", LenMat," are non-zero."
!!Allocate the rho matrix
!        ALLOCATE(Mat(LenMat),stat=ierr)
!        CALL LogMemAlloc('Mat',LenMat,8*HElement_t_size,this_routine,MatTag,ierr)
!        Mat=0.0_dp
!!
!!!Lab indicates the column that the 'i'th non-zero matrix element resides in the full matrix
!        ALLOCATE(Lab(LenMat),stat=ierr)
!        CALL LogMemAlloc('Lab',LenMat,4,this_routine,LabTag,ierr)
!        Lab(1:LenMat)=0
!
!!Count the number of non-zero elements in the matrix - the way the diagonaliser seems to want to...
!!Also now only top triangle
!!        LenMat=0
!!        ICMax=0
!!        do i=1,NDets
!!            RowElems=0
!!            do j=i,NDets
!!                IF(abs(GraphRhoMat(i,j)).gt.0.0_dp) THEN
!!                    LenMat=LenMat+1
!!                    RowElems=RowElems+1
!!                ENDIF
!!            enddo
!!            IF(RowElems.gt.ICMax) ICMax=RowElems
!!        enddo
!!!        WRITE(6,"(A,I10,A,I10,A)") "Of ",((NDets*NDets+NDets)/2), " possible elements, only ", LenMat," are non-zero."
!!
!!        ALLOCATE(Mat(NDets,ICMax),stat=ierr)
!!        CALL LogMemAlloc('Mat',NDets*ICMax,8,this_routine,MatTag,ierr)
!!        Mat=0.0_dp
!!
!!!Lab now indicates the position of the non-zero element in the row specified
!!        ALLOCATE(Lab(NDets,ICMax),stat=ierr)
!!        CALL LogMemAlloc('Lab',NDets*ICMax,4,this_routine,LabTag,ierr)
!!        Lab(1:NDets,1:ICMax)=0
!
!!NRow indicated the number of non-zero elements in the 'i'th row of the full matrix
!        ALLOCATE(NRow(NDets),stat=ierr)
!        CALL LogMemAlloc('NRow',NDets,4,this_routine,NRowTag,ierr)
!        NRow(1:NDets)=0
!
!!This compresses the rho matrix, so that only the non-zero elements are stored, in a suitable form for the Lanczos diagonaliser.
!        CALL CompressMatrix(Mat,Lab,NRow,LenMat,ICMax)
!!        WRITE(6,"(A,I8)") "In compressed matrix, maximum number of non-zero elements in any one row is ", ICMax
!
!!Set up parameters for the diagonaliser.
!!NEval indicates the number of eigenvalues we want to calculate
!        IF(NEval.eq.0) THEN
!!NEval is set to 0 by default, computing all eigenvectors
!            IF(NDets.gt.25) THEN
!                WRITE(6,*) "Resetting NEval to 10 eigenvectors"
!                NEval=10
!            ELSE
!                WRITE(6,*) "Computing all eigenvectors"
!                NEval=NDets
!            ENDIF
!        ENDIF
!!        B2L=1.0e-25_dp
!!        NBlk=4
!!        NKry=8
!        NCycle=200
!        NKry1=NKry+1
!        NBlock=MIN(NEval,NBlk)
!        LScr=MAX(NDets*NEval,8*NBlock*NKry)
!        LIScr=6*NBlock*NKry
!
!!Check whether RAN2 has been initialised yet, or not
!        IF(Seed.eq.0) THEN
!            TSeeded=.true.
!        ELSE
!            TSeeded=.false.
!        ENDIF
!!        WRITE(6,*) "TSeeded: ", TSeeded
!
!!Deallocate GraphRhoMat - no longer needed
!        DEALLOCATE(GraphRhoMat)
!        CALL LogMemDealloc(this_routine,GraphRhoMatTag)
!
!!Allocate memory for diagonaliser
!        ALLOCATE(A(NEval,NEval),stat=ierr)
!        CALL LogMemAlloc('A',NEval*NEval,8,this_routine,ATag,ierr)
!        A=0.0_dp
!        ALLOCATE(V(NDets*NBlock*NKry1),stat=ierr)
!        CALL LogMemAlloc('V',NDets*NBlock*NKry1,8,this_routine,VTag,ierr)
!        V=0.0_dp
!        ALLOCATE(AM(NBlock*NBlock*NKry1),stat=ierr)
!        CALL LogMemAlloc('AM',NBlock*NBlock*NKry1,8,this_routine,AMTag,ierr)
!        AM=0.0_dp
!        ALLOCATE(BM(NBlock*NBlock*NKry),stat=ierr)
!        CALL LogMemAlloc('BM',NBlock*NBlock*NKry,8,this_routine,BMTag,ierr)
!        BM=0.0_dp
!        ALLOCATE(T(3*NBlock*NKry*NBlock*NKry),stat=ierr)
!        CALL LogMemAlloc('T',NBlock*NKry*NBlock*NKry*3,8,this_routine,TTag,ierr)
!        T=0.0_dp
!        ALLOCATE(WT(NBlock*NKry),stat=ierr)
!        CALL LogMemAlloc('WT',NBlock*NKry,8,this_routine,WTTag,ierr)
!        WT=0.0_dp
!        ALLOCATE(SCR(LScr),stat=ierr)
!        CALL LogMemAlloc('SCR',LScr,8,this_routine,SCRTag,ierr)
!        SCR=0.0_dp
!        ALLOCATE(ISCR(LIScr),stat=ierr)
!        CALL LogMemAlloc('ISCR',LIScr,4,this_routine,ISCRTag,ierr)
!        ISCR(1:LIScr)=0
!        ALLOCATE(Index(NEval),stat=ierr)
!        CALL LogMemAlloc('Index',NEval,4,this_routine,IndexTag,ierr)
!        Index(1:NEval)=0
!        ALLOCATE(WH(NDets),stat=ierr)
!        CALL LogMemAlloc('WH',NDets,8,this_routine,WHTag,ierr)
!        WH=0.0_dp
!        ALLOCATE(Work2(3*NDets),stat=ierr)
!        CALL LogMemAlloc('Work2',3*NDets,8,this_routine,Work2Tag,ierr)
!        Work2=0.0_dp
!        ALLOCATE(V2(NDets,NEval),stat=ierr)
!        CALL LogMemAlloc('V2',NDets*NEval,8,this_routine,V2Tag,ierr)
!        V2=0.0_dp
!
!!W holds the eigenvalues 
!        ALLOCATE(W(NEval),stat=ierr)
!        CALL LogMemAlloc('W',NEval,8,this_routine,WTag,ierr)
!        W=0.0_dp
!
!!CK holds the eigenvectors
!        ALLOCATE(CK(NDets,NEval),stat=ierr)
!        CALL LogMemAlloc('CK',NDets*NEval,8,this_routine,CKTag,ierr)
!!The initial trial wavefuntion is set to zero        
!        CK=0.0_dp
!        ALLOCATE(CKN(NDets,NEval),stat=ierr)
!        CALL LogMemAlloc('CKN',NDets*NEval,8,this_routine,CKNTag,ierr)
!        CKN=0.0_dp
!
!!Lanczos iterative diagonalisation routine
!        IF(THDiag) THEN
!!If using Hamiltonian matrix in the diagonliser, we want the smallest eigenvalues, not the largest
!            CALL NECI_FRSBLKH(NDets,ICMax,NEval,Mat,Lab,CK,CKN,NKry,NKry1,NBlock,NRow,LScr,LIScr,A,W,
                !V,AM,BM,T,WT,SCR,ISCR,Index,WH,Work2,V2,NCycle,B2L,.false.,.false.,TSeeded,.true.)
!        ELSE
!            CALL NECI_FRSBLKH(NDets,ICMax,NEval,Mat,Lab,CK,CKN,NKry,NKry1,NBlock,NRow,LScr,LIScr,A,
                !W,V,AM,BM,T,WT,SCR,ISCR,Index,WH,Work2,V2,NCycle,B2L,.false.,.true.,TSeeded,.true.)
!        ENDIF
!!Mulitply eigenvalues through by -1 to ensure they are positive - no longer needed
!!        CALL DSCAL(NEval,-1.0_dp,W,1)
!!If first element of eigenvector is positive, multiply whole lot through by -1
!!       IF(CK(1,1).gt.0.0_dp) THEN
!!           CALL DSCAL(NDets,-1.0_dp,CK(:,1),1)
!!       ENDIF
!
!!        do i=1,NDets
!!            WRITE(6,*) NRow(i)
!!        enddo
!!        WRITE(6,*) "***********"
!
!!        do i=1,NEval
!!            WRITE(6,"(A,I4,A,F20.14)") "Eigenvalue ",i," is: ", W(i)
!!        enddo
!
!!Deallocate memory required by diagonaliser (including original matrix)
!        DEALLOCATE(Mat)
!        CALL LogMemDealloc(this_routine,MatTag)
!        DEALLOCATE(Lab)
!        CALL LogMemDealloc(this_routine,LabTag)
!        DEALLOCATE(A)
!        CALL LogMemDealloc(this_routine,ATag)
!        DEALLOCATE(V)
!        CALL LogMemDealloc(this_routine,VTag)
!        DEALLOCATE(AM)
!        CALL LogMemDealloc(this_routine,AMTag)
!        DEALLOCATE(BM)
!        CALL LogMemDealloc(this_routine,BMTag)
!        DEALLOCATE(T)
!        CALL LogMemDealloc(this_routine,TTag)
!        DEALLOCATE(WT)
!        CALL LogMemDealloc(this_routine,WTTag)
!        DEALLOCATE(SCR)
!        CALL LogMemDealloc(this_routine,ScrTag)
!        DEALLOCATE(WH)
!        CALL LogMemDealloc(this_routine,WHTag)
!        DEALLOCATE(Work2)
!        CALL LogMemDealloc(this_routine,Work2Tag)
!        DEALLOCATE(V2)
!        CALL LogMemDealloc(this_routine,V2Tag)
!        
!        
!!        ALLOCATE(temp(NDets))
!!        temp=0.0_dp
!!        OPEN(24,FILE='fort.23',Status='old')
!!        do i=1,NDets
!!            READ(24,*) temp(i)
!!        enddo
!
!!Store largest eigenvector - ?first? column of CK (zero it)
!        Eigenvector=(0.0_dp)
!
!!        SumVec=0.0_dp
!!        LancVar=0.0_dp
!        do i=1,NDets
!!Need to record the largest/smallest eigenvector - depending on whether THDiag is on or off
!            Eigenvector(i)=(CK(i,1))
!!            LancVar=LancVar+ABS(temp(i)-Eigenvector(i))/ABS(Eigenvector(i))
!!            SumVec=SumVec+((Eigenvector(i))**2)
!!            WRITE(22,*) Eigenvector(i)
!!            IF(i.le.NEval) THEN
!!                WRITE(6,*) Eigenvector(i),W(i)
!!            ELSE
!!                WRITE(6,*) Eigenvector(i)
!!            ENDIF
!            IF(TMoveDets.and.(ABS(Eigenvector(i))).eq.0.0_dp) THEN
!!There are still the possibility of disconnected clusters - these will be removed by regrowing graph completly...
!                IF(Iteration.eq.1) THEN
!                    WRITE(6,*) "Disconnected cluster found in graph - performing one cycle of regrowing new graph 
                                        !from scratch..."
!                    ReturntoTMoveDets=.true.
!!Choose a low graph bias, which will allow the graph to be created easily
!                    TBiasing=.true.
!                    GraphBias=0.75
!                ELSE
!                    WRITE(6,*) "Disconnected clusters still found...exiting..."
!                    STOP "Disconnected clusters still found...exiting..."
!                ENDIF
!            ENDIF
!        enddo
!        
!!        DEALLOCATE(temp)
!!        WRITE(6,*) "The variance of the eigenvector components is: ", LancVar/NDets
!
!!        IF((ABS(SumVec-1.0_dp)).gt.1.0e-8_dp) THEN
!!            WRITE(6,*) "Eigenvector NOT NORMALISED!",SumVec
!!        ENDIF
!        IF(ReturntoTMoveDets) THEN
!            TMoveDets=.false.
!        ENDIF
!
!!Also need to save the largest/smallest Eigenvalue
!        Eigenvalue=W(1)
!
!!Find weight and energy of graph.
!!There is no beta-dependance, so only largest/smallest eigenvector needed.
!!Note - since we are not having a beta-dependance, 'weight' takes on a slightly different
!!meaning. Here, the weight is simply the square of the first element of the largest/smallest eigenvector,
!!i.e. the magnitude of the projection of the graph back onto the HF...
!!        WRITE(6,*) "First element of eigenvector is: ", Eigenvector(1)
!        SI=(Eigenvector(1))*(Eigenvector(1))
!        IF(THDiag) THEN
!!If diagonlising the Hamiltonian matrix, then the energy is given by the smallest eigenvalue
!            IF(Eigenvalue.lt.0.0_dp) THEN
!                DLWDB=Eigenvalue
!            ELSE
!                DLWDB=Eigenvalue*(-1.0_dp)
!            ENDIF
!        ELSE
!            DLWDB=0.0_dp
!            do i=2,NDets
!                DLWDB=DLWDB+(HamElems(i))*(Eigenvector(i))
!            enddo
!            DLWDB=(DLWDB/(Eigenvector(1)))+(HamElems(1))
!        ENDIF
!
!!Deallocate Eigenvectors and values
!        DEALLOCATE(W)
!        CALL LogMemDealloc(this_routine,WTag)
!        DEALLOCATE(CK)
!        CALL LogMemDealloc(this_routine,CKTag)
!        DEALLOCATE(CKN)
!        CALL LogMemDealloc(this_routine,CKNTag)
!        DEALLOCATE(Index)
!        CALL LogMemDealloc(this_routine,IndexTag)
!        DEALLOCATE(ISCR)
!        CALL LogMemDealloc(this_routine,IScrTag)
!        DEALLOCATE(NRow)
!        CALL LogMemDealloc(this_routine,NRowTag)
!!        WRITE(6,*) LabTag,ATag,MatTag,NRowTag,VTag,WTTag,AMTag,BMTag,TTag,SCRTag,ISCRTag,
                !IndexTag,WHTag,Work2Tag,V2Tag,WTag,CKTag,CKNTag
!
!        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in DiagGraphLanc'
!        call halt_timer(proc_timerLanc)
!
!    END SUBROUTINE DiagGraphLanc
!
!!This compresses the rho matrix, so that only the non-zero elements are stored, in a suitable form for the Lanczos diagonaliser.
!    SUBROUTINE CompressMatrix(Mat,Lab,NRow,LenMat,ICMax)
!        IMPLICIT NONE
!        INTEGER :: LenMat,i,j
!        real(dp) :: Mat(LenMat)
!        INTEGER :: Lab(LenMat),NRow(NDets),LabElem,RowElems,ICMax,SumElems,lenrow
!!        real(dp) :: Mat(NDets,ICMax)
!!        INTEGER :: Lab(NDets,ICMax),NRow(NDets),LabElem,RowElems,ICMax,SumElems
!
!!This compression is how Alex uses the Lanczos diagonaliser - however, it seems to be used in a different way
!!Here, only the non-zero elements are stored, and it is stored symmetrically
!!Lab now indicates the column that the ith non-zero matrix resides in
!!ICMax indicates the largest number of non-zero elements in any one row of the matrix
!
!        lenrow=ICMax
!        SumElems=0
!        ICMax=0
!        LabElem=0
!        do i=1,NDets
!!Scan along each row of the matrix
!            RowElems=0
!!Only scan the upper triangle of the matrix
!            do j=i,NDets
!                IF(abs(GraphRhoMat(i,j)).gt.0.0_dp) THEN
!                    LabElem=LabElem+1
!                    RowElems=RowElems+1
!                    Mat(LabElem)=GraphRhoMat(i,j)
!                    Lab(LabElem)=j
!                ENDIF
!            enddo
!            IF(RowElems.gt.ICMax) ICMax=RowElems
!            NRow(i)=RowElems
!            SumElems=SumElems+RowElems
!        enddo
!        IF(ICMax.ne.lenrow) THEN
!            WRITE(6,*) "Stop - error in compressing matrix"
!            STOP "Stop - error in compressing matrix"
!        ENDIF
!
!!This is the way to compress the matrix, as it seems in the actual diagonaliser...
!!Here, the matrix is only partially compressed, with it being of size (NDets,ICMax)
!!        SumElems=0
!!!
!!        do i=1,NDets
!!            RowElems=0
!!            do j=i,NDets
!!                IF(abs(GraphRhoMat(i,j)).gt.0.0_dp) THEN
!!                    RowElems=RowElems+1
!!                    Mat(i,RowElems)=GraphRhoMat(i,j)
!!                    Lab(i,RowElems)=j
!!                ENDIF
!!            enddo
!!            IF(RowElems.gt.ICMax) THEN
!!                WRITE(6,*) "Error in compressing matrix 3"
!!                STOP 'Error in compressing matrix 3'
!!            ENDIF
!!            NRow(i)=RowElems
!!            SumElems=SumElems+RowElems
!!        enddo
!
!        IF(SumElems.ne.LenMat) THEN
!            WRITE(6,*) "Error in compressing matrix"
!            STOP 'Error in compressing matrix'
!        ENDIF
!
!    END SUBROUTINE CompressMatrix
!
!!This is a routine to find the energy of the graph by diagonalisation, and return as well its 
!eigenvalues and largest eigenvector. Deallocate the RhoMatrix when done.
!    SUBROUTINE DiagGraphMorph()
!        IMPLICIT NONE
!        type(timer), save :: proc_timerDiag
!        INTEGER :: Info,ierr,i
!        real(dp) , ALLOCATABLE :: Work(:),Eigenvalues(:)
!        INTEGER :: WorkTag,EigenvaluesTag,j
!        CHARACTER(len=*), PARAMETER :: this_routine='DiagGraphMorph'
!
!        proc_timerDiag%timer_name='DiagGraphMorph'
!        call set_timer(proc_timerDiag)
!
!        WorkTag=0
!        EigenvaluesTag=0
!
!        ALLOCATE(Eigenvalues(NDets),stat=ierr)
!        CALL LogMemAlloc('Eigenvalues',NDets,8,this_routine,EigenvaluesTag,ierr)
!        Eigenvalues=0.0_dp
!
!!Workspace needed for diagonaliser
!        ALLOCATE(Work(3*NDets),stat=ierr)
!        CALL LogMemAlloc('WORK',NDets*3,8,this_routine,WorkTag,ierr)
!
!        IF(TMoveDets) THEN
!!If using the MoveDets algorithm, need to keep the rho matrix for the previous graph
!            ALLOCATE(CopyRhoMat(NDets,NDets),stat=ierr)
!            CALL LogMemAlloc('CopyRhoMat',NDets**2,8,this_routine,CopyRhoMatTag,ierr)
!            do i=1,NDets
!                do j=i,NDets
!!                    IF(i.eq.NDets) WRITE(6,*) j,GraphRhoMat(i,j)
!!                    IF(j.eq.NDets) WRITE(6,*) i,GraphRhoMat(j,i)
!                    CopyRhoMat(i,j)=GraphRhoMat(i,j)
!                    CopyRhoMat(j,i)=GraphRhoMat(j,i)
!                enddo
!            enddo
!        ENDIF
!
!!Diagonalise...
!!Eigenvalues now stored in ascending order
!        CALL DSYEV('V','U',NDets,GraphRhoMat,NDets,Eigenvalues,Work,3*NDets,Info)
!        IF(Info.ne.0) THEN
!            WRITE(6,*) "DYSEV error in DiagGraphMorph: ",Info
!            STOP
!        ENDIF
!
!        DEALLOCATE(Work)
!        CALL LogMemDealloc(this_routine,WorkTag)
!
!!Store largest eigenvector - last column of GraphRhoMat (zero it)
!        Eigenvector=(0.0_dp)
!
!        do i=1,NDets
!            IF(THDiag) THEN
!!If we are diagonalising the hamiltonian matrix, rather than the rho-matrix, then we want the 
!eigenvector corresponding to the smallest, not the largest eigenvalue
!                Eigenvector(i)=GraphRhoMat(i,1)
!            ELSE
!                Eigenvector(i)=GraphRhoMat(i,NDets)
!            ENDIF
!!            WRITE(6,*) Eigenvector(i),Eigenvalues(NDets-(i-1))
!!            WRITE(23,*) Eigenvector(i)
!            IF(TMoveDets.and.(ABS(Eigenvector(i))).eq.0.0_dp) THEN
!!There are still the possibility of disconnected clusters - these will be removed by regrowing graph completly...
!                IF(Iteration.eq.1) THEN
!                    WRITE(6,*) "Disconnected cluster found in graph - performing one cycle of regrowing new graph 
                                    !from scratch..."
!                    ReturntoTMoveDets=.true.
!!Choose a low graph bias, which will allow the graph to be created easily
!                    TBiasing=.true.
!                    GraphBias=0.75
!                ELSE
!                    WRITE(6,*) "Disconnected clusters still found...exiting..."
!                    STOP "Disconnected clusters still found...exiting..."
!                ENDIF
!            ENDIF
!!            WRITE(6,*) i,Eigenvector(i)
!        enddo
!
!        IF(ReturntoTMoveDets) THEN
!            TMoveDets=.false.
!        ENDIF
!
!!Also need to save the largest Eigenvalue
!!        do i=1,5
!!            WRITE(6,*) "Eigenvalue ", i," is: ",Eigenvalues(NDets-(i-1))
!!        enddo
!        IF(THDiag) THEN
!!If we are diagonlising the hamiltonian matrix, then we want to take the smallest eigenvalue, not the largest
!            Eigenvalue=Eigenvalues(1)
!        ELSE
!            Eigenvalue=Eigenvalues(NDets)
!        ENDIF
!
!!Deallocate RhoMatrix and Eigenvalues (not needed)
!        DEALLOCATE(Eigenvalues)
!        CALL LogMemDealloc(this_routine,EigenvaluesTag)
!        DEALLOCATE(GraphRhoMat)
!        CALL LogMemDealloc(this_routine,GraphRhoMatTag)
!
!!Find weight and energy of graph.
!!There is no beta-dependance, so only largest eigenvector needed.
!!Note - since we are not having a beta-dependance, 'weight' takes on a slightly different
!!meaning. Here, the weight is simply the square of the first element of the largest eigenvector,
!!i.e. the magnitude of the projection of the graph back onto the HF...
!!        WRITE(6,*) "First element of eigenvector is: ", Eigenvector(1)
!        SI=(Eigenvector(1))*(Eigenvector(1))
!        IF(THDiag) THEN
!            IF(Eigenvalue.lt.0.0_dp) THEN
!                DLWDB=Eigenvalue
!            ELSE
!                DLWDB=Eigenvalue*(-1.0_dp)
!            ENDIF
!        ELSE
!            DLWDB=0.0_dp
!            do i=2,NDets
!                DLWDB=DLWDB+(HamElems(i))*(Eigenvector(i))
!            enddo
!            DLWDB=(DLWDB/(Eigenvector(1)))+(HamElems(1))
!        ENDIF
!
!        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in DiagGraphMorph'
!        call halt_timer(proc_timerDiag)
!
!    END SUBROUTINE DiagGraphMorph

END MODULE GraphMorph

!LOGICAL FUNCTION ConnectGraph(Matrix,Dimen,NEl)
!    IMPLICIT NONE
!    INTEGER :: NotConnected(NEl,Dimen),NEl,Dimen
!    real(dp) :: Matrix
!    LOGICAL :: ConnectGraph
!
!    j=0
!    do i=1,Dimen
!!Test if not connected to the excitation level above
!        IF(Matrix(1,i).eq.0) THEN
!            j=j+1
!            NotConnected(1,j)=i
!        ENDIF

FUNCTION RootofNum(Num,Root)
    use constants, only: dp
    real(dp) :: Root
    real(dp) :: Num,RootofNum
    IF(Num.lt.1.0e-16_dp) THEN
        RootofNum=0.0_dp
    ELSE
        RootofNum=EXP(LOG(Num)/Root)
    ENDIF
    RETURN
END FUNCTION RootofNum

!Tests determinants are the same - requires the same ordering of orbitals in them
LOGICAL FUNCTION SameDet(nI,nJ,NEl)
    IMPLICIT NONE
    INTEGER :: NEl,nI(NEl),nJ(NEl),i
    SameDet=.true.
    do i=1,NEl
        IF(nI(i).ne.nJ(i)) THEN
            SameDet=.false.
            EXIT
        ENDIF
    enddo
    RETURN
END FUNCTION
