#include "macros.h"
module GreensFuncUtils
    use GreensFuncData
    use SystemData, only: nel, nBasis,tHPHF,tFixLz, G1,tNoBrillouin,tUseBrillouin
    use SystemData, only: nBasisMax,tKPntSym,ElecPairs,nOccAlpha,nOccBeta
    use global_utilities
    use constants, only: dp, int64, n_int
    use FciMCData, only: SpawnedParts,SpawnedParts2,AllNoatHF,OldAllHFCyc,NoatHF
    use FciMCData, only: HFDet,HFDetTag,ProjEDet,iLutRef,HashIndex,iLutHF
    use FciMCData, only: CurrentH,HolesInList,AllTotWalkers,AllTotWalkersOld
    use FciMCData, only: AllTotParts,AllTotPartsOld,TotParts,TotPartsOld,nClashMax
    use FciMCData, only: tHashWalkerList,TotWalkers,TotWalkersOld,CurrentDets,MaxSpawned
    use FciMCData, only: nWalkerHashes,OldAllAvWalkersCyc,OldAllNoatHF
    use FciMCData, only: WalkVecDetsTag,WalkVecHTag,SpawnVecTag,SpawnVec2Tag
    use FciMCData, only: ValidSpawnedList,InitialSpawnedSlots,Hii
    use FciMCParMod, only: BinSearchParts3
    use SymExcitDataMod, only: SpinOrbSymLabel
    use GenRandSymExcitNUMod , only: RandExcitSymLabelProd
    use PopsfileMod, only: readpopsheadv4
    use Logging, only: tIncrementPops
    use util_mod, only: get_unique_filename,get_free_unit
    use sort_mod
    use sym_mod, only: GetSym,GetLz
    use bit_reps, only: NIfD,NIfDBO,NIfTot,extract_sign,decode_bit_det,decode_bit_det_nel
    use AnnihilationMod, only: IsUnoccDet,DetermineDetNode,FindWalkerHash,DetermineDetNode_nel
    use DetBitOps, only: DetBitEQ
    use Parallel_neci

    contains

    !Function to check whether the given j,i orbital pair make a desired transisition matrix
    !element of the greens function to calculate, based on symmetry, spin, whether we want
    !diagonal contributions, and the specific transistions requested
    logical function CheckAllowed_i(j,i)
        implicit none
        integer, intent(in) :: j,i
        integer :: PairIndex,k
        logical :: tFoundPair

        CheckAllowed_i=.true.
        if(tNoDiag.and.(i.eq.j)) then
            CheckAllowed_i = .false. 
            return
        endif
        if(G1(i)%Ms.ne.G1(j)%Ms) then
            CheckAllowed_i = .false.  !Check spin symmetry
            return
        endif
        !Check spation/kpoint symmetry
        if(RandExcitSymLabelProd(SpinOrbSymLabel(i),SpinOrbSymLabel(j)).ne.0) then
            CheckAllowed_i = .false.
            return
        endif
        PairIndex=FindPairIndex(j,i)
        if(tSpecificMatEls) then
            tFoundPair=.false.
            do k=1,iSpecificMatEls
                if(SpecificMatEls(k).eq.PairIndex) then
                    tFoundPair=.true.
                    exit
                endif
            enddo
            if(.not.tFoundPair) then   !Not a desired i transistion
                CheckAllowed_i=.false.
            endif
        endif

    end function CheckAllowed_i

    subroutine SetupSpecificMatEls()
        implicit none
        integer :: FoundPair,i,j,k,l,PairIndex,ierr
        logical :: jAlreadyCalc,tFoundPair
        character(len=*), parameter :: t_r="SetupSpecificMatEls"

        !First, turn the list of desired transistions into triangular indexed array
        allocate(SpecificMatEls(iSpecificMatEls),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc Error")
        do i=1,iSpecificMatEls
            PairIndex=FindPairIndex(SpecificMatEls_temp(1,i),SpecificMatEls_temp(2,i))
            SpecificMatEls(i)=PairIndex
        enddo
        deallocate(SpecificMatEls_temp)

        write(6,"(A)") "Only certain elements of the greens function will be considered."
        write(6,"(A)") "These are: "
        FoundPair = 0
        Contributing_Js = 0

        tNoDiag=.true. !This will be turned off if any of the transistions are diagonal

        do j=1,nBasis
            jAlreadyCalc=.false.
            do i=1,nBasis
!                if(tNoDiag.and.(i.eq.j)) cycle
                if(G1(i)%Ms.ne.G1(j)%Ms) cycle  !Check spin symmetry
                !Check spation/kpoint symmetry
                if(RandExcitSymLabelProd(SpinOrbSymLabel(i),SpinOrbSymLabel(j)).ne.0) cycle
                PairIndex=FindPairIndex(j,i)
!                write(6,*) j,i,PairIndex

                tFoundPair=.false.
                do k=1,iSpecificMatEls
!                    write(6,*) k,PairIndex,SpecificMatEls(k)
                    if(PairIndex.eq.SpecificMatEls(k)) then
                        tFoundPair=.true.
                        exit
                    endif
                enddo
                if(tFoundPair) then
                    !We want to calculate this pair
                    FoundPair=FoundPair+1
                    write(6,"(A,I5,A,I5,A)") "(",j,")->-(",i,")"

                    if(i.eq.j) tNoDiag=.false.

                    if(.not.jAlreadyCalc) then
                        jAlreadyCalc = .true.
                        Contributing_Js = Contributing_Js + 1
                    endif
                endif
            enddo
        enddo
        if(FoundPair.ne.iSpecificMatEls) then
            write(6,*) FoundPair, iSpecificMatEls
            call stop_all(t_r,"Cannot find the same number of pairs as in the input file")
        endif
        if(.not.tJustSpectrum) then
            if(Contributing_Js.gt.1) then
                write(6,"(I6,A)") Contributing_Js," j orbitals will need to be run over to &
                                    &calculate desired matrix elements"
            else
                write(6,"(I6,A)") Contributing_Js," j orbital will need to be run over to &
                                    &calculate desired matrix elements"
            endif
        endif
        allocate(Allowed_j(Contributing_Js),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Error allocating")
        Allowed_j(:) = 0

        !Run over allowed orbitals again, to now store the allowed j's
        l=0
        do j=1,nBasis
            do i=1,nBasis
                if(tNoDiag.and.(i.eq.j)) cycle
                if(G1(i)%Ms.ne.G1(j)%Ms) cycle  !Check spin symmetry
                !Check spation/kpoint symmetry
                if(RandExcitSymLabelProd(SpinOrbSymLabel(i),SpinOrbSymLabel(j)).ne.0) cycle
                PairIndex=FindPairIndex(j,i)

                tFoundPair=.false.
                do k=1,iSpecificMatEls
                    if(PairIndex.eq.SpecificMatEls(k)) then
                        tFoundPair=.true.
                        exit
                    endif
                enddo
                if(tFoundPair) then
                    !We want to calculate this pair
                    l=l+1
                    Allowed_j(l)=j
                    exit    !No need to run through the other i's
                endif
            enddo
        enddo

        if(l.ne.Contributing_Js) then
            call stop_all(t_r,"Wrong number of j orbitals calculated to run over")
        endif
        !Should already be sorted in increasing j
        if(.not.tJustSpectrum) then
            if(Contributing_Js.gt.1) then
                write(6,"(A)",advance='no') "These orbitals are: "
            else
                write(6,"(A)",advance='no') "This orbital is: "
            endif
            write(6,"(I6)") Allowed_j(1:Contributing_Js)
        endif

    end subroutine SetupSpecificMatEls

    !Return nel to its original value here at the end of a GF calculation with given GFType
    subroutine ReturnElectronNumber(GFType)
        implicit none
        integer, intent(in) :: GFType

        nel=nel_orig
        elecpairs=(nel*(nel-1))/2
        nOccAlpha=nOccAlpha_orig
        nOccBeta=nOccBeta_orig
    end subroutine ReturnElectronNumber


    function FindPairIndex(j,i) 
        implicit none
        integer , intent(in) :: j,i
        integer :: FindPairIndex

        FindPairIndex=(j-1)*nBasis+i
    end function FindPairIndex

    subroutine FindDetPos(k,MinInd,PartInd,tSuccess)
        implicit none
        integer, intent(in) :: k,MinInd
        integer, intent(out) :: PartInd
        logical, intent(out) :: tSuccess
        integer :: clash,FinalVal,nI(nel),DetHash
        character(len=*) , parameter :: this_routine="FindDetPos"

        if(tHashWalkerList) then
            tSuccess=.false.
            call decode_bit_det (nI, SpawnedParts(:,k))              
            DetHash=FindWalkerHash(nI)
            FinalVal=HashIndexSaved(0,DetHash)-1
            do clash=1,FinalVal
                ASSERT(HashIndexSaved(clash,DetHash).le.TotWalkersSaved)
                if(DetBitEQ(SpawnedParts(:,k),CurrentDetsSaved(:,HashIndexSaved(clash,DetHash)),NIfDBO)) then
                    tSuccess=.true.
                    PartInd=HashIndexSaved(clash,DetHash)
                    exit
                endif
            enddo
        else
            CALL BinSearchParts3(SpawnedParts(:,k),CurrentDetsSaved,TotWalkersSaved,MinInd,TotWalkersSaved,PartInd,tSuccess)
        endif
    end subroutine FindDetPos
                
    subroutine SendDetsToProcessors(MaxIndex)
        implicit none
        integer(MPIArg), dimension(nProcessors) :: sendcounts, disps 
        integer(MPIArg), dimension(nProcessors) :: recvcounts, recvdisps
        integer(kind=n_int) , pointer :: PointTemp(:,:)
        integer, intent(out) :: MaxIndex
        integer :: error,MaxSendIndex,i
        character(*), parameter :: t_r="SendDetsToProcessors"

        if(nProcessors.eq.1) then
            sendcounts(1)=ValidSpawnedList(0)-1
            disps(1)=0
        else
            do i=0,nProcessors-1
                sendcounts(i+1)=ValidSpawnedList(ProcNode(i))-InitialSpawnedSlots(ProcNode(i))
                disps(i+1)=InitialSpawnedSlots(ProcNode(i))-1
            enddo
        endif
        MaxSendIndex=ValidSpawnedList(nNodes-1)-1
        recvcounts(:)=0
        call MPIBarrier(error)
        call MPIAlltoAll(sendcounts,1,recvcounts,1,error)
        recvdisps(1)=0
        do i=2,nProcessors
            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
        enddo
        MaxIndex=recvdisps(nProcessors)+recvcounts(nProcessors)
        do i=1,nProcessors
            recvdisps(i)=recvdisps(i)*(NIfTot+1)
            recvcounts(i)=recvcounts(i)*(NIfTot+1)
            sendcounts(i)=sendcounts(i)*(NIfTot+1)
            disps(i)=disps(i)*(NIfTot+1)
        enddo

        if(MaxIndex.gt.(0.9*MaxSpawned)) then
            write(6,*) MaxIndex,MaxSpawned
            call Warning_neci(t_r,"Maximum index of spawned array is close " &
                & //"to maximum length after send. Increase MemoryFacSpawn")
        endif
        call MPIAlltoAllv(SpawnedParts,sendcounts,disps,SpawnedParts2,recvcounts,recvdisps,error)
        PointTemp => SpawnedParts2
        SpawnedParts2 => SpawnedParts
        SpawnedParts => PointTemp
        
    end subroutine SendDetsToProcessors

    subroutine SetupForGF(Orb,GFType)
        use SystemData, only : BasisFN
        use CalcData, only: InitWalkers,MemoryFacPart,tRegenDiagHEls,iReadWalkersRoot
        use FciMCParMod, only: SetupValidSpawned,init_excit_gen_store,CalcApproxpDoubles
        use FciMCParMod, only: CountExcitsOld
        use FciMCData, only: WalkVecDets,SpawnVec,SpawnVec2,HashIndexArr1,MaxWalkersPart
        use FciMCData, only: freeslot,WalkVecH,fcimc_excit_gen_store,HighestPopDet
        use SymExcit3, only: CountExcitations3
        implicit none
        integer, intent(in) :: GFType,Orb
        type(BasisFN) :: HFSym
        integer :: ierr,WalkerListSize,ReadBatch,HFLz,i,Ms
        integer :: HFConn,nSingles,nDoubles,error,exflag
        character(*), parameter :: t_r="SetupForGF"

        if(GFType.eq.1) then
            allocate(HFDet(NEl-1),stat=ierr)
            call LogMemAlloc('HFDet',NEl-1,8,t_r,HFDetTag,ierr)
            allocate(ProjEDet(NEl-1),stat=ierr)
        else
            allocate(HFDet(NEl+1),stat=ierr)
            call LogMemAlloc('HFDet',NEl+1,8,t_r,HFDetTag,ierr)
            allocate(ProjEDet(NEl+1),stat=ierr)
        endif
        allocate(ilutHF(0:NIfTot),stat=ierr)
        allocate(iLutRef(0:NIfTot),stat=ierr)
        allocate(HighestPopDet(0:NIfTot),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc Err")
        !These will all be set later

        exflag=3
        tNoBrillouin=.true.
        tUseBrillouin=.false.
        WalkerListSize=InitWalkers
        MaxWalkersPart=NINT(MemoryFacPart*WalkerListSize)
        Call SetupValidSpawned(int(WalkerListSize,int64))  
        ALLOCATE(WalkVecDets(0:NIfTot,MaxWalkersPart),stat=ierr)
        call LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfTot+1),size_n_int,t_r,WalkVecDetsTag,ierr)
        WalkVecDets(0:NIfTot,1:MaxWalkersPart)=0
        if(tRegenDiagHEls) call stop_all(t_r,"Diagonal HElements should be stored")
        ALLOCATE(WalkVecH(MaxWalkersPart),stat=ierr)
        call LogMemAlloc('WalkVecH',MaxWalkersPart,8,t_r,WalkVecHTag,ierr)
        WalkVecH(:)=0.d0
        ALLOCATE(SpawnVec(0:NIftot,MaxSpawned),stat=ierr)
        call LogMemAlloc('SpawnVec',MaxSpawned*(NIfTot+1),size_n_int,t_r,SpawnVecTag,ierr)
        SpawnVec(:,:)=0
        ALLOCATE(SpawnVec2(0:NIfTot,MaxSpawned),stat=ierr)
        call LogMemAlloc('SpawnVec2',MaxSpawned*(NIfTot+1),size_n_int,t_r,SpawnVec2Tag,ierr)
        SpawnVec2(:,:)=0
        SpawnedParts=>SpawnVec
        SpawnedParts2=>SpawnVec2
        if(tHashWalkerList) then
            nClashMax=int(real(WalkerListSize,dp)/real(nWalkerHashes,dp))+1
            !HashIndex: (0,:) is first free slot in the hash list.
            allocate(HashIndexArr1(0:nClashMax,nWalkerHashes),stat=ierr)
            HashIndex => HashIndexArr1
            HashIndex(:,:)=0
            HashIndex(0,:)=1    !First free slot is first slot
            !Also need to allocate memory for the freeslot array
            allocate(FreeSlot(MaxWalkersPart),stat=ierr)
            freeslot(:)=0
        endif
        CurrentDets=>WalkVecDets
        CurrentH=>WalkVecH
        TotParts(:)=0
        TotPartsOld(:)=0
        NoatHF=0
        if(iReadWalkersRoot.eq.0) then
            ReadBatch=MaxSpawned  
        else
            ReadBatch = iReadWalkersRoot
        endif
        ! Initialise excitation generation storage
        call init_excit_gen_store (fcimc_excit_gen_store)
        
        !nel changes in here
        call ReadInForGF(Orb,GFType,ReadBatch,CurrentDets,MaxWalkersPart)

        !Now calculate these variables - shouldn't actually need to change them
        !since they are never explicitly used - only conserved.
        if(tFixLz) then
            call GetLz(HFDet,NEl,HFLz)
            WRITE(6,"(A,I5)") "Ml value of reference determinant is: ",HFLz
        endif
        CALL GetSym(HFDet,NEl,G1,NBasisMax,HFSym)
        WRITE(6,"(A,I10)") "Symmetry of reference determinant is: ",INT(HFSym%Sym%S,sizeof_int)
        MS=0
        do i=1,nel
            MS=MS+G1(HFDet(i))%Ms
        enddo
        write(6,*) "2 x Spin polarization of reference determinant is: ",MS
        call CalcApproxpDoubles()

        IF(.not.tKPntSym) THEN
            !Count all possible excitations - put into HFConn
            CALL CountExcitations3(HFDet,exflag,nSingles,nDoubles)
        ELSE
            !use Alex's old excitation generators to enumerate all excitations.
            CALL CountExcitsOld(HFDet,exflag,nSingles,nDoubles)
        ENDIF
        HFConn=nSingles+nDoubles
        !Setup global variables
        TotWalkersOld=TotWalkers
        TotPartsOld = TotParts
        call MPISumAll(TotWalkers,AllTotWalkers)
        AllTotWalkersOld = AllTotWalkers
        call MPISumAll(TotParts,AllTotParts)
        AllTotPartsOld=AllTotParts
        call MPISumAll(NoatHF,AllNoatHF)
        OldAllNoatHF=AllNoatHF
        OldAllAvWalkersCyc=real(sum(AllTotParts),dp)
        if(lenof_sign.eq.1) then
            OldAllHFCyc = real(AllNoatHF(1),dp)
        else
            OldAllHFCyc = cmplx(real(AllNoatHF(1),dp),real(AllNoatHF(lenof_sign),dp), dp)
        endif

        CALL MPIBarrier(error)

    end subroutine SetupforGF
    
    !EndPopsList = no. dets in file
    !ReadBatch = how many to read in one go
    subroutine ReadInForGF(Orb,GFType,ReadBatch,Dets,DetsLen)
        use PopsfileMod, only: open_pops_head,FindPopsfileVersion,read_popsfile_det
        use CalcData, only: iPopsFileNoRead
        use FciMCData, only: RandomHash
        use AnnihilationMod, only: EnlargeHashTable
        use Determinants, only : get_helement
        implicit none
        integer , intent(in) :: Orb,GFType,ReadBatch,DetsLen
        integer(int64) :: EndPopsList
        integer(n_int) , intent(out) :: Dets(0:NIfTot,DetsLen) 
        integer :: iunit,PopsVersion,iPopLenof_sign,iPopNel,i,j,ierr,err
        integer :: iPopIter,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot
        integer(int64) :: iPopAllTotWalkers,acc
        integer(int64) , dimension(lenof_sign) :: PopSumNoatHF
        integer(n_int) , allocatable :: BatchRead(:,:)
        HElement_t :: PopAllSumENum,HElemTemp
        real(dp) :: PopDiagSft,read_tau,BatchSize,LowestEnergy
        real(dp) :: LowDet(2),LowDetOut(2)
        integer :: LowestIndex,PopBlockingIter
        logical :: formpops,binpops,tPop64Bit,tPopHPHF,tPopLz
        character(255) :: popsfile
        integer , dimension(lenof_sign) :: CurrHF,SignTemp
        integer(int64), dimension(lenof_sign) :: CurrParts
        integer :: CurrWalkers,nBatches,DetRemoved,WalkersRemoved,tempnel
        integer :: AllWalkersRemoved,AllDetRemoved,DetHash
        integer :: PopsSendList(0:nNodes-1),proc,TempNi(NEl),newelec
        integer :: PopsInitialSlots(0:nNodes-1),MaxSendIndex,Slot
        integer(int64) :: Det,TempCurrWalkers,AllCurrWalkers
        integer(n_int) :: WalkerTemp(0:NIfTot)
        logical :: tReadAllPops
        integer , allocatable :: SampledGFnI(:)
        integer(MPIArg) :: sendcounts(nNodes),disps(nNodes),recvcount
        character(*), parameter :: t_r="ReadInForGF"

        if(allocated(SampledGFnI)) deallocate(SampledGFnI)
        if(GFType.eq.1) then
            allocate(SampledGFnI(nel-1),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err")
        else
            allocate(SampledGFnI(nel+1),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err")
        endif

        sendcounts=0
        disps=0
        MaxSendIndex=1

        call open_pops_head(iunit,formpops,binpops)
        IF(FormPops) THEN
            !determine version number
            PopsVersion=FindPopsfileVersion(iunit)
            if(PopsVersion.ne.4) call stop_all(t_r,"Wrong POPS verion")
            call ReadPopsHeadv4(iunit,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot,read_tau,PopBlockingIter)
        ELSE
            call ReadPopsHeadv4(iunit,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot,read_tau,PopBlockingIter)
            if(iProcIndex.eq.root) then
                close(iunit)
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                OPEN(iunit,FILE=popsfile,Status='old',form='unformatted')
            endif
        ENDIF
        EndPopsList=iPopAllTotWalkers
        BatchSize=REAL(ReadBatch,dp)/REAL(nNodes,dp)
        if(iProcIndex.eq.Root) then
            !Create PopsInitialSlots
            do i=0,nNodes-1
                PopsInitialSlots(i)=NINT(BatchSize*i)+1
            enddo
            !Allocate array to store particle to distribute
            allocate(BatchRead(0:NIfTot,1:ReadBatch),stat=ierr)
        else
            allocate(BatchRead(0:NIfTot,1:MaxSendIndex),stat=ierr)
        endif

        CurrHF=0        !Number of HF walkers on each node.
        CurrParts=0     !Number of walkers on each node.
        CurrWalkers=0   !Number of determinants on each node.
        nBatches=0      !Number of batches of walkers it takes to distribute popsfile.
        Det=1
        tReadAllPops=.false.
        DetRemoved=0
        WalkersRemoved=0
        HolesInList=0
        if(GFType.eq.1) then
            tempnel=nel-1
        else
            tempnel=nel+1
        endif
        do while(.not.tReadAllPops)

            if(iProcIndex.eq.Root) then

                !Get ready for reading in the next batch of walkers
                nBatches=nBatches+1
                BatchRead(:,:)=0
                PopsSendList(:)=PopsInitialSlots(:)
                do while(Det.le.EndPopsList)

                    !Read the next entry, and store the walker in WalkerTemp and TempnI
                    call read_popsfile_det(iunit,Det,BinPops,WalkerTemp,TempnI)

!                    write(6,*) Det,WalkerTemp,TempnI(:),Orb
                    !Remove/add the orbital
                    !GFType 1 = ADVANCED GF: simulate N-1 electron system and propagate forwards in IT
                    !GFType 2 = RETARDED GF: simulate N+1 electron system and propagate backwards in IT
                    if(GFType.eq.1) then
                        !Remove Orb if present, and keep walkers
                        !If not present, then remove all walkers
                        if(IsOcc(WalkerTemp,Orb)) then
                            clr_orb(WalkerTemp,Orb)
                            newelec=1
                            do i=1,nel  !(original nel)
                                if(TempnI(i).ne.Orb) then
                                    SampledGFnI(newelec)=TempnI(i)
                                    newelec=newelec+1
                                endif
                            enddo
                        else
                            DetRemoved=DetRemoved+1
                            call extract_sign(WalkerTemp(:),SignTemp)
                            WalkersRemoved=WalkersRemoved+sum(abs(SignTemp))
                            cycle   !Ignore all walkers here
                        endif
                    else
                        !Add orb if not present, and keep walkers
                        !If present, then remove all walkers
                        if(IsNotOcc(WalkerTemp,Orb)) then
                            set_orb(WalkerTemp,Orb)
                            SampledGFnI(1:nel)=TempnI(:)
                            SampledGFnI(nel+1)=Orb
                            call sort(SampledGFnI)
                        else
                            DetRemoved=DetRemoved+1
                            call extract_sign(WalkerTemp(:),SignTemp)
                            WalkersRemoved=WalkersRemoved+sum(abs(SignTemp))
                            cycle
                        endif
                    endif

!                    proc = DetermineDetNode (TempnI,0)
                    !Do this explicitly just here, since nel is currently incorrect
                    do i = 1, tempnel
!                        write(6,*) "***",tempnel,i,SampledGFnI(i)
                        acc = (1099511628211_int64 * acc) + &
                                (RandomHash(mod(SampledGFnI(i)-1,int(nBasis,int64))+1) * i)
                    enddo
                    proc = abs(mod(acc, int(nNodes,int64)))

                    BatchRead(:,PopsSendList(proc)) = WalkerTemp(:)
                    PopsSendList(proc) = PopsSendList(proc) + 1
                    if(proc.ne.(nNodes-1)) then
                        if(PopsInitialSlots(proc+1)-PopsSendList(proc).lt.2) then
                            exit  !Now distribute the particles
                        endif
                    else
                        if(ReadBatch-PopsSendList(proc).lt.2) then
                            exit  !Now distribute the particles
                        endif
                    endif

                enddo

                if(Det.gt.EndPopsList) tReadAllPops=.true.

                do j=0,nNodes-1
!                    sendcounts(j+1)=(PopsSendList(j)-(NINT(BatchSize*j)+1))*(NIfTot+1)
                    sendcounts(j+1)=(PopsSendList(j)-PopsInitialSlots(j))*(NIfTot+1)
!                    disps(j+1)=(NINT(BatchSize*j))*(NIfTot+1)
                    disps(j+1)=(PopsInitialSlots(j)-1)*(NIfTot+1)
                enddo
                MaxSendIndex=(disps(nNodes)+sendcounts(nNodes))/(nIfTot+1)

            endif

            !Now scatter the particles read in to their correct processors.
            if(bNodeRoot) call MPIScatter(sendcounts,recvcount,err,Roots)
            if(err.ne.0) call stop_all(t_r,"MPI scatter error")
            if(bNodeRoot) then
!                call MPIScatterV(BatchRead(:,1:MaxSendIndex),sendcounts,disps,Dets(:,CurrWalkers+1:DetsLen),recvcount,err,Roots)
                call MPIScatterV(BatchRead(:,1:MaxSendIndex),sendcounts,disps,  &
                    Dets(:,CurrWalkers+1:(recvcount/NIfTot+1)),recvcount,err,Roots)
            endif
            if(err.ne.0) call stop_all(t_r,"MPI error")
            if(bNodeRoot) CurrWalkers=CurrWalkers+recvcount/(NIfTot+1)
            call MPIBCast(tReadAllPops)

        enddo

        if(iProcIndex.eq.Root) close(iunit)
        
        if(GFType.eq.1) then
            !Advanced GF
            nel=nel-1
            if(is_alpha(Orb)) then
                nOccAlpha=nOccAlpha-1
            else
                nOccBeta=nOccBeta-1
            endif
        else
            !Retarded GF
            nel=nel+1
            if(is_alpha(Orb)) then
                nOccAlpha=nOccAlpha+1
            else
                nOccBeta=nOccBeta+1
            endif
        endif
        elecpairs=(nel*(nel-1))/2
        write(6,*) "Changing NEL to :",nel

        !Test we have still got all determinants
        TotWalkers=int(CurrWalkers,int64)
        call MPISum(TotWalkers,1,AllTotWalkers)
        call MPISum(WalkersRemoved,1,AllWalkersRemoved)
        call MPISum(DetRemoved,1,AllDetRemoved)
        if(iProcIndex.eq.Root) then
            if(AllTotWalkers.ne.(EndPopsList-DetRemoved)) then
                write(6,*) "AllTotWalkers: ",AllTotWalkers
                write(6,*) "EndPopsList: ",EndPopsList
                call Stop_All(t_r,"Not all walkers accounted for when reading in")
            endif
            write(6,*) "Determinants removed from simulation: ",AllDetRemoved
            write(6,*) "Walkers removed from simulation: ",AllWalkersRemoved
        endif

        deallocate(BatchRead)
        
        write(6,*) "Number of configurations read in to this process: ",TotWalkers 
        if(tHashWalkerList) then
            do i=1,TotWalkers
                call decode_bit_det (SampledGFnI, dets(:,i))              
                DetHash=FindWalkerHash(SampledGFnI)
                Slot=HashIndex(0,DetHash)
                HashIndex(Slot,DetHash)=i
                HashIndex(0,DetHash)=HashIndex(0,DetHash)+1
                if(HashIndex(0,DetHash).gt.nClashMax) then
                    call EnlargeHashTable()
                endif
            enddo
        else
            !Order the determinants on all the lists.
            call sort (dets(:,1:TotWalkers))
        endif

        LowestEnergy=1.0e10_dp
        LowestIndex=0
        !Run through all determinants on each node, and calculate the total number of walkers, and noathf
        do i=1,TotWalkers
!            WRITE(6,*) i,Dets(:,i)

            !Find diagonal element
            call decode_bit_det (SampledGFnI, Dets(:,i))
            HElemTemp = get_helement (SampledGFnI, SampledGFnI, 0)
            CurrentH(i)=REAL(HElemTemp,dp)

            if(real(HElemTemp,dp).lt.LowestEnergy) then
                LowestEnergy=real(HElemTemp,dp)
                LowestIndex=i
            endif

            call extract_sign(Dets(:,i),SignTemp)
!            write(6,*) i,SampledGFnI(:),SignTemp
            CurrParts=CurrParts+abs(SignTemp)
        enddo
        TotParts=CurrParts
        call MPISum(TotParts,AllTotParts)

        !Find the *global* lowest energy determinant
        LowDet(1)=LowestEnergy
        LowDet(2)=real(iProcIndex,dp)
        call MPIAllReduceDatatype(LowDet,1,MPI_MINLOC,MPI_2DOUBLE_PRECISION,LowDetOut)
        LowestEnergy = LowDetOut(1)

        if(iProcIndex.eq.nint(LowDetOut(2))) then
            !This has the reference det
            ilutHF(:)=Dets(:,LowestIndex)
            call extract_sign(Dets(:,LowestIndex),SignTemp)
            CurrHF=CurrHF+SignTemp
        endif
        NoatHF=CurrHF
        call MPISum(NoatHF,AllNoatHF)

        call MPIBCast(ilutHF,NIfTot+1,nint(LowDetOut(2)))
        !Broadcast the lowest energy determinant to all processes

        !Now, we know what we want our reference to be.
        do i=1,TotWalkers
            CurrentH(i)=CurrentH(i)-LowestEnergy
        enddo
        
        ilutRef(:)=ilutHF(:)
        call decode_bit_det (HFDet, ilutHF)
        ProjEDet(:)=HFDet(:)
        write(6,*) "Lowest energy reference found to be: ",HFDet(:)
        Hii = LowestEnergy
        write(6,*) "Reference energy: ",Hii           

        call mpibarrier(err)

    end subroutine ReadInForGF

    subroutine SetGFDefaults()
        implicit none

        tFreqGF = .true.    !By default, sample frequency domain GF
        tGreensFuncs=.false.
        tJustAdvGF=.false.
        tJustRetGF=.false.
        nIterEquilGF=1000
        GFSamples=100
        GreensImTime(:)=1.0_dp
        nTimePnts_ret=-1
        nTimePnts_adv=-1
        tNoDiag=.false.
        iSpecificMatEls=0
        tSpecificMatEls=.false.
        tJustSpectrum=.false.
        tCalcSpectrum=.true.
        StartFreq_adv = -1.0_dp
        FinalFreq_adv = 1.0_dp
        dFreqStep_adv = 0.1_dp
        GreensIters(:) = 100
        GreenTimeslice(:) = -1
        tFinalOverlapOnly = .false.

    end subroutine SetGFDefaults 

    subroutine GFReadInput()
        use input_neci
        implicit none
        logical eof
        integer :: i_orbital,j_orbital,ierr
        integer, allocatable :: temp(:,:)
        character(len=*), parameter :: t_r="GFReadInput"
        character(len=100) w

        !The presence of a GF block is enough to activate the GF code
        tGreensFuncs=.true.

        GF: do
            call read_line(eof)
            if(eof) then
                exit
            endif
            call readu(w)
            select case(w)
            case("SAMPLES")
                call readi(GFSamples)
            case("IMAGTIMEGF")
                tFreqGF=.false. !Go into the imaginary time GF calculation
            case("ADVFREQRANGE")
                !The advanced frequency GF range
                call readf(StartFreq_adv)
                call readf(FinalFreq_adv)
                call readf(dFreqStep_adv)
            case("JUSTADVANCED")
                tJustAdvGF=.true.
            case("JUSTRETARDED")
                tJustRetGF=.true.
            case("ADVIMTIME")
                call readf(GreensImTime(1))
            case("RETIMTIME")
                call readf(GreensImTime(2))
            case("ITERATIONS_ADV")
                call readi(GreensIters(1))
            case("ITERATIONS_RET")
                call readi(GreensIters(2))
            case("ADVINTEGRATIONSTEP")
                call readi(GreenTimeslice(1))
            case("RETINTEGRATIONSTEP")
                call readi(GreenTimeslice(2))
            case("ADVTIMEPOINTS")
                call readi(nTimePnts_adv)
            case("RETTIMEPOINTS")
                call readi(nTimePnts_ret)
            case("EQUILIBRATION")
                call readi(nIterEquilGF)
            case("NODIAGGF")
                tNoDiag=.true.  !Do not calculate diagonal elements of GFs
            case("TRANSISTION")
                !This may be called multiple times.
                !If specified, it will only calculate these matrix elements
                call readi(j_orbital)
                call readi(i_orbital)
                iSpecificMatEls=iSpecificMatEls+1
                if(.not.tSpecificMatEls) then
                    tSpecificMatEls=.true.
                    allocate(SpecificMatEls_temp(2,1),stat=ierr)
                    if(ierr.ne.0) call stop_all(t_r,"Alloc Err")
                    SpecificMatEls_temp(1,1) = j_orbital
                    SpecificMatEls_temp(2,1) = i_orbital
                else
                    allocate(temp(2,iSpecificMatEls-1),stat=ierr)
                    if(ierr.ne.0) call stop_all(t_r,"Alloc Err")
                    temp(:,:)=SpecificMatEls_temp(:,:)
                    deallocate(SpecificMatEls_temp)
                    allocate(SpecificMatEls_temp(2,iSpecificMatEls),stat=ierr)
                    if(ierr.ne.0) call stop_all(t_r,"Alloc Err")
                    SpecificMatEls_temp(:,1:iSpecificMatEls-1)=temp(:,:)
                    SpecificMatEls_temp(1,iSpecificMatEls) = i_orbital 
                    SpecificMatEls_temp(2,iSpecificMatEls) = j_orbital
                    deallocate(temp)
                endif
            case("JUSTSPECTRUM")
                tJustSpectrum=.true.
            case("NOSPECTRUM")
                tCalcSpectrum=.false.
            case("FINALOVERLAPONLY")
                tFinalOverlapOnly=.true.
            case("END")
                exit
            case default
                write(6,"(A)") "ALLOWED KEYWORDS IN GREENSFUNCS BLOCK: "
                write(6,"(A)") "ADVFREQRANGE"
                write(6,"(A)") "IMAGTIMEGF"
                write(6,"(A)") "SAMPLES"
                write(6,"(A)") "JUSTADVANCED"
                write(6,"(A)") "JUSTRETARDED"
                write(6,"(A)") "ADVIMTIME"
                write(6,"(A)") "RETIMTIME"
                write(6,"(A)") "ADVINTEGRATIONSTEP"
                write(6,"(A)") "RETINTEGRATIONSTEP"
                write(6,"(A)") "ADVTIMEPOINTS"
                write(6,"(A)") "RETTIMEPOINTS"
                write(6,"(A)") "EQUILIBRATION"
                write(6,"(A)") "NODIAGGF"
                write(6,"(A)") "TRANSISTION"
                write(6,"(A)") "JUSTSPECTRUM"
                write(6,"(A)") "NOSPECTRUM"
                write(6,"(A)") "FINALOVERLAPONLY"
                call report("Keyword "//trim(w)//" not recognized",.true.)
            end select
        enddo GF

    end subroutine GFReadInput

end module GreensFuncUtils
    
!Wrap greens function routines since we want to use the FCIMC module routines in it
subroutine CalcGreensFunc_wrap
    use GreensFuncData, only: tFreqGF
    implicit none

    if(tFreqGF) then
        call CalcGreensFuncFreq_wrap()
    else
        call CalcGreensFuncIT_wrap()
    endif

end subroutine CalcGreensFunc_wrap


