#include "macros.h"

MODULE PopsfileMod

    use SystemData, only: nel, tHPHF, tFixLz, tCSF, nBasis, tNoBrillouin,tMomInv
    use CalcData, only: tTruncInitiator, DiagSft, tWalkContGrow, nEquilSteps, &
                        ScaleWalkers, tReadPopsRestart, tRegenDiagHEls, &
                        InitWalkers, tReadPopsChangeRef, nShiftEquilSteps, &
                        iWeightPopRead, iPopsFileNoRead, tPopsMapping, Tau, &
                        InitiatorWalkNo
    use DetBitOps, only: DetBitLT,FindBitExcitLevel,DetBitEQ
    use hash , only : DetermineDetNode
    use Determinants, only : get_helement,write_det
    use hphf_integrals, only: hphf_diag_helement
    use MI_integrals, only: MI_diag_helement
    USE dSFMT_interface , only : genrand_real2_dSFMT
    use FciMCData
    use bit_reps
    use Parallel_neci
!    use AnnihilationMod, only: DetermineDetNode,FindWalkerHash,EnlargeHashTable,IsUnoccDet
    use AnnihilationMod, only: FindWalkerHash,EnlargeHashTable
    USE Logging , only : iWritePopsEvery,tPopsFile,iPopsPartEvery,tBinPops
    USE Logging , only : tPrintPopsDefault,tIncrementPops, tPrintInitiators
    use sort_mod
    use util_mod, only: get_free_unit,get_unique_filename

    implicit none

    contains

    !   V.3/4 POPSFILE ROUTINES   !
!This routine reads in particle configurations from a POPSFILE v.3-4.
!EndPopsList is the number of entries in the POPSFILE to read, and ReadBatch is the number of determinants
!which can be read in in a single batch.
    SUBROUTINE ReadFromPopsfile(EndPopsList,ReadBatch,CurrWalkers64,CurrParts,CurrHF,Dets,DetsLen)
        use MemoryManager, only: TagIntType
        integer(int64) , intent(in) :: EndPopsList  !Number of entries in the POPSFILE.
        integer , intent(in) :: ReadBatch       !Size of the batch of determinants to read in in one go.
        integer(int64) , intent(out) :: CurrWalkers64    !Number of determinants which end up on a given processor.
        real(dp), intent(out) :: CurrHF(lenof_sign)
        real(dp) :: CurrParts(lenof_sign)
        integer :: CurrWalkers,Slot,nJ(nel)
        integer :: iunit,i,j,ierr,PopsInitialSlots(0:nNodes-1)
        INTEGER(TagIntType) :: BatchReadTag=0
        real(dp) :: BatchSize
        integer :: PopsSendList(0:nNodes-1),proc
        integer(MPIArg) :: sendcounts(nNodes), disps(nNodes), recvcount
        integer :: MaxSendIndex,err,DetHash
        integer(n_int) , allocatable :: BatchRead(:,:)
        integer(n_int) :: WalkerTemp(0:NIfTot)
        integer(int64) :: Det,AllCurrWalkers,TempCurrWalkers
        logical :: FormPops,BinPops,tReadAllPops,tStoreDet
        real(dp) , dimension(lenof_sign) :: SignTemp
        integer :: TempNI(NEl),nBatches,PopsVersion
        character(len=*) , parameter :: this_routine='ReadFromPopsfile'
        HElement_t :: HElemTemp
        character(255) :: popsfile
        !variables from header file
        logical :: tPop64Bit,tPopHPHF,tPopLz
        integer :: iPopLenof_sign,iPopNEl,iPopIter,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot
        integer :: PopBlockingIter
        integer(int64) :: iPopAllTotWalkers
        real(dp) :: PopDiagSft,read_tau
        real(dp) , dimension(lenof_sign) :: PopSumNoatHF
        integer, intent(in) :: DetsLen
        INTEGER(kind=n_int), intent(out) :: Dets(0:nIfTot,DetsLen)
        HElement_t :: PopAllSumENum

        sendcounts=0
        disps=0
        MaxSendIndex=1
      
        call open_pops_head(iunit,formpops,binpops)
        IF(FormPops) THEN
            !determine version number
            PopsVersion=FindPopsfileVersion(iunit)
            if(PopsVersion.eq.3) then
                call ReadPopsHeadv3(iunit,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                    iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                    PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot)
            else
                call ReadPopsHeadv4(iunit,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                    iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                    PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot,read_tau,PopBlockingIter)
            endif

            if(EndPopsList.ne.iPopAllTotWalkers) then
                call stop_all(this_routine,"Error in assessing number of entries in POPSFILE")
            endif

        ELSEIF(BinPops) THEN
            if(iProcIndex.eq.root) then
                close(iunit)    !iunit here refers to the header file.
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                OPEN(iunit,FILE=popsfile,Status='old',form='unformatted')
            endif
        ENDIF

        call mpibarrier(err)

        IF(iProcIndex.eq.Root) THEN
            IF(iWeightPopRead.ne.0) THEN
                WRITE(6,"(A,I15,A,I4,A)") "Although ",EndPopsList, &
                " configurations will be read in, only determinants with a weight of over ",iWeightPopRead," will be stored."
            else
                write(6,"(A,I15,A)") "Reading in a total of ",EndPopsList, " configurations from POPSFILE."
            ENDIF
            if(ScaleWalkers.ne.1) then
                call warning_neci(this_routine,"ScaleWalkers parameter found, but not implemented in POPSFILE v3 - ignoring.")
            endif
            call neci_flush(6)
        ENDIF

        if(tPopsMapping) then
            write(6,"(A)") "Reading in popsfile of smaller basis, and distributing on larger basis determinants"
            call init_popsfile_mapping()
        endif

        BatchSize=REAL(ReadBatch,dp)/REAL(nNodes,dp)
        if(iProcIndex.eq.Root) then
            !Create PopsInitialSlots
            do i=0,nNodes-1
                PopsInitialSlots(i)=NINT(BatchSize*i)+1
            enddo
            !Allocate array to store particle to distribute
            allocate(BatchRead(0:NIfTot,1:ReadBatch),stat=ierr)
            CALL LogMemAlloc('BatchRead',ReadBatch*(NIfTot+1),size_n_int,this_routine,BatchReadTag,ierr)
            write(6,"(A,I12,A)") "Reading in a maximum of ",ReadBatch," determinants at a time from POPSFILE."
            call neci_flush(6)
        else
            allocate(BatchRead(0:NIfTot,1:MaxSendIndex),stat=ierr)
            CALL LogMemAlloc('BatchRead',MaxSendIndex*(NIfTot+1),size_n_int,this_routine,BatchReadTag,ierr)
        endif

        CurrHF=0.0_dp        !Number of HF walkers on each node.
        CurrParts=0.0     !Number of walkers on each node.
        write(6,*) "ReadBatch: ",ReadBatch
        write(6,*) "MaxSendIndex: ",MaxSendIndex

        CurrWalkers=0   !Number of determinants on each node.
        nBatches=0      !Number of batches of walkers it takes to distribute popsfile.
        Det=1
        tReadAllPops=.false.
        do while(.not.tReadAllPops)

            if(iProcIndex.eq.Root) then

                !Get ready for reading in the next batch of walkers
                nBatches=nBatches+1
                BatchRead(:,:)=0
                PopsSendList(:)=PopsInitialSlots(:)
                do while(Det.le.EndPopsList)

!                    write(6,*) Det,EndPopsList
                    !Read the next entry, and store the walker in WalkerTemp and TempnI
                    call read_popsfile_det(iunit,Det,BinPops,WalkerTemp,TempnI)

                    proc = DetermineDetNode (TempnI,0)
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
                    sendcounts(j+1)=int((PopsSendList(j)-PopsInitialSlots(j))*(NIfTot+1),MPIArg)
!                    disps(j+1)=(NINT(BatchSize*j))*(NIfTot+1)
                    disps(j+1)=int((PopsInitialSlots(j)-1)*(NIfTot+1),MPIArg)
                enddo
                MaxSendIndex=(disps(nNodes)+sendcounts(nNodes))/(nIfTot+1)

            endif

            !Now scatter the particles read in to their correct processors.
            if(bNodeRoot) call MPIScatter(sendcounts,recvcount,err,Roots)
            if(err.ne.0) call stop_all(this_routine,"MPI scatter error")
            if(bNodeRoot) then
!                call MPIScatterV(BatchRead(:,1:MaxSendIndex),sendcounts,disps,Dets(:,CurrWalkers+1:DetsLen),recvcount,err,Roots)
                call MPIScatterV(BatchRead(:,1:MaxSendIndex),sendcounts,disps,  &
                    Dets(:,CurrWalkers+1:(recvcount/NIfTot+1)),recvcount,err,Roots)
            endif
            if(err.ne.0) call stop_all(this_routine,"MPI error")
            if(bNodeRoot) CurrWalkers=CurrWalkers+recvcount/(NIfTot+1)
            call MPIBCast(tReadAllPops)

        enddo

        if(iProcIndex.eq.Root) close(iunit)

        !Test we have still got all determinants
        write(6,*) "CurrWalkers: ",CurrWalkers
        TempCurrWalkers=int(CurrWalkers,int64)
        call MPISum(TempCurrWalkers,1,AllCurrWalkers)
        if(iProcIndex.eq.Root) then
            if((iWeightPopRead.eq.0).and.(AllCurrWalkers.ne.EndPopsList)) then
                write(6,*) "AllCurrWalkers: ",AllCurrWalkers
                write(6,*) "EndPopsList: ",EndPopsList
                call Stop_All(this_routine,"Not all walkers accounted for when reading in")
            endif
        endif

        deallocate(BatchRead)
        CALL LogMemDealloc(this_routine,BatchReadTag)

        write(6,"(A,I8)") "Number of batches required to distribute all determinants in POPSFILE: ",nBatches
        write(6,*) "Number of configurations read in to this process: ",CurrWalkers 

        if(tHashWalkerList) then
            do i=1,CurrWalkers
                call decode_bit_det (nJ, dets(:,i))              
                DetHash=FindWalkerHash(nJ)
                Slot=HashIndex(0,DetHash)
                HashIndex(Slot,DetHash)=i
                HashIndex(0,DetHash)=HashIndex(0,DetHash)+1
                if(HashIndex(0,DetHash).gt.nClashMax) then
                    call EnlargeHashTable()
                endif
            enddo
        else
            !Order the determinants on all the lists.
            call sort (dets(:,1:CurrWalkers))
        endif

        !Run through all determinants on each node, and calculate the total number of walkers, and noathf
        do i=1,CurrWalkers
!            WRITE(6,*) i,Dets(:,i)
            call extract_sign(Dets(:,i),SignTemp)
            CurrParts=CurrParts+abs(SignTemp)
            if(DetBitEQ(Dets(:,i),iLutRef,NIfDBO)) then
                if(CurrHF(1).ne.0) then
                    call stop_all(this_routine,"HF already found, but shouldn't have")
                endif
                CurrHF=CurrHF+SignTemp 
                IF(.not.tRegenDiagHEls) CurrentH(1,i)=0.0_dp
            else
                if(.not.tRegenDiagHEls) THEN
                !Calculate diagonal matrix element
                    call decode_bit_det (TempnI, currentDets(:,i))
                    if (tHPHF) then
                        HElemTemp = hphf_diag_helement (TempnI,Dets(:,i))
                    elseif(tMomInv) then
                        HElemTemp = MI_diag_helement(TempnI,Dets(:,i))
                    else
                        HElemTemp = get_helement (TempnI, TempnI, 0)
                    endif
                    CurrentH(1,i)=REAL(HElemTemp,dp)-Hii
                endif
            endif
        enddo

        CurrWalkers64=int(CurrWalkers,int64)    !Since this variable is eventually going to be
        call mpibarrier(err)
                                                !Totwalkers, it wants to be a 64 bit int.

        if(allocated(PopsMapping)) deallocate(PopsMapping)
    
    end subroutine ReadFromPopsfile

    !This routine reads the next determinant entry from a popsfile and stores it in WalkerTemp,
    !ready to be distributed.
    !This should only be called from root node.
    !iunit = unit
    !Det = Current determinant entry in popsfile list
    !BinPops = Binary popsfile or formatted
    !WalkerTemp = Determinant entry returned (in new basis if using mapping)
    subroutine read_popsfile_det(iunit,Det,BinPops,WalkerTemp,nI)
        use DetBitOps , only : EncodeBitDet
        integer , intent(in) :: iunit
        integer(int64) , intent(inout) :: Det
        logical , intent(in) :: BinPops
        integer(n_int), intent(out) :: WalkerTemp(0:NIfTot)
        integer , intent(out) :: nI(NEl)
        integer(n_int) :: WalkerToMap(0:MappingNIfD),WalkerTemp2(0:NIfTot)
        integer :: elec, flg, i, j
        real(dp) :: sgn(lenof_sign)
        logical :: tStoreDet

        tStoreDet=.false.
        do while(.not.tStoreDet)

            if (.not. tPopsMapping) then
                ! All basis parameters match --> Read in directly.
                if (BinPops) then
                    read(iunit) WalkerTemp(0:NIfD), sgn
                    if (tUseFlags) read(iunit) flg
                else
                    if (tUseFlags) then
                        read(iunit,*) WalkerTemp(0:NIfD), sgn, flg
                    else
                        read(iunit,*) WalkerTemp(0:NIfD), sgn
                    end if
                end if
            else
                ! we are mapping from a smaller to larger basis.
                if (BinPops) then
                    read(iunit) WalkerToMap(0:MappingNIfD), sgn
                    if (tUseFlags) read(iunit) flg
                else
                    if (tUseFlags) then
                        read(iunit,*) WalkerToMap(0:MappingNIfD), sgn, flg
                    else
                        read(iunit,*) WalkerToMap(0:MappingNIfD), sgn
                    end if
                end if

                ! Decode to natural ordered integers, and then re-encode to
                ! the new representation.
                elec = 0
outer_map:      do i = 0, MappingNIfD
                    do j = 0, end_n_int
                        if (btest(WalkerToMap(i), j)) then
                            elec = elec + 1
                            nI(elec) = PopsMapping((i * bits_n_int) + (j + 1))
                            if (elec == nel) exit outer_map
                        end if
                    end do
                end do outer_map

                ! Encode this new determinant.
                call EncodeBitDet(nI, WalkerTemp)
            end if

            ! Store the sign and flag information in the determinant.
            call encode_sign (WalkerTemp, sgn)
            if (tUseFlags) call encode_flags (WalkerTemp, flg)

            ! Increment the determinant counter
            det = det + 1

            ! Test if we actually want to store this walker...
            if (iWeightPopRead /= 0) then
                do i = 1, lenof_sign
                    if (sgn(i) >= iWeightPopRead) then
                        tStoreDet = .true.
                        exit
                    end if
                end do
            else
                tStoreDet = .true.
            end if
        enddo

        ! Decode the determinant as required if not using mapping
        if (.not. tPopsMapping) then
            call decode_bit_det (nI, WalkerTemp)
        endif

    end subroutine read_popsfile_det

    
    !Read in mapping file, and store the mapping between orbitals in the old and new bases.
    subroutine init_popsfile_mapping()
        integer :: OldBasisSize,NewBasisSize,iunit,NewBasis,kkx,kky,kkz,dumint,OldBasis,i
        real(dp) :: dumEn
        logical :: exists
        character(len=*), parameter :: t_r='init_popsfile_mapping'

        if(iProcIndex.eq.Root) then
            inquire(file='mapping',exist=exists)
            if(.not.exists) then
                call stop_all(t_r,"No mapping file found")
            endif

            iunit=get_free_unit()
            open(iunit,file='mapping',status='old',action='read',form='formatted')

            OldBasisSize=0
            NewBasisSize=0
            do i=1,nBasis
                read(iunit,*) NewBasis,kkx,kky,kkz,dumint,dumEn,OldBasis
                if(NewBasis.gt.NewBasisSize) NewBasisSize=NewBasis
                if(OldBasis.gt.OldBasisSize) OldBasisSize=OldBasis
            enddo
            if(NewBasisSize.ne.nBasis) then
                call stop_all(t_r,"Size of new bases not consistent in mapping")
            endif
            if(int(OldBasisSize/bits_n_int).ne.MappingNIfD) then
                call stop_all(t_r,"Size of old basis not consistent with its POPSFILEHEAD NIfD")
            endif

            allocate(PopsMapping(OldBasisSize))   !Mapping from old basis to new basis
            PopsMapping=0
            rewind(iunit)
            do i=1,nBasis 
                read(iunit,*) NewBasis,kkx,kky,kkz,dumint,dumEn,OldBasis
                if(OldBasis.ne.0) then
                    PopsMapping(OldBasis)=NewBasis
                endif
            enddo
            
            do i=1,OldBasis
                if(PopsMapping(i).eq.0) call stop_all(t_r,"Not all orbital mappings accounted for")
            enddo

            close(iunit)

            write(6,"(A,I4,A,I4)") "Original NIfTot: ",MappingNIfTot, " New NIfTot: ",NIfTot

        endif

    end subroutine init_popsfile_mapping
                

    subroutine CheckPopsParams(tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                    iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                    PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot,WalkerListSize,read_tau,PopBlockingIter)
        use Logging , only : tZeroProjE
        logical , intent(in) :: tPop64Bit,tPopHPHF,tPopLz
        integer , intent(in) :: iPopLenof_sign,iPopNel,iPopIter,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot
        integer , intent(in) :: PopBlockingIter
        integer(int64) , intent(in) :: iPopAllTotWalkers
        real(dp) , intent(in) :: PopDiagSft,read_tau
        real(dp) , dimension(lenof_sign) , intent(in) :: PopSumNoatHF
        HElement_t , intent(in) :: PopAllSumENum
        integer , intent(out) :: WalkerListSize
        character(len=*) , parameter :: this_routine='CheckPopsParams'

        !Ensure all NIF and symmetry options the same as when popsfile was written out.
#ifdef __INT64
        if(.not.tPop64Bit) call stop_all(this_routine,"Popsfile created with 32 bit walkers, but now using 64 bit.")
#else
        if(tPop64Bit) call stop_all(this_routine,"Popsfile created with 64 bit walkers, but now using 32 bit.")
#endif
        if(tPopHPHF.neqv.tHPHF) call stop_all(this_routine,"Popsfile HPHF and input HPHF not same")
        if(tPopLz.neqv.tFixLz) call stop_all(this_routine,"Popsfile Lz and input Lz not same")
        if(iPopLenof_sign.ne.lenof_sign) call stop_all(this_routine,"Popsfile lenof_sign and input lenof_sign not same")
        if(iPopNEl.ne.NEl) call stop_all(this_routine,"Popsfile NEl and input NEl not same")
        if(.not.tPopsMapping) then
            if(PopNIfD.ne.NIfD) call stop_all(this_routine,"Popsfile NIfD and calculated NIfD not same")
        else
            MappingNIfD=PopNIfD
        endif
        if(PopNIfY.ne.NIfY) call stop_all(this_routine,"Popsfile NIfY and calculated NIfY not same")
        if(PopNIfSgn.ne.NIfSgn) call stop_all(this_routine,"Popsfile NIfSgn and calculated NIfSgn not same")
        if(PopNIfFlag.ne.NIfFlag) call stop_all(this_routine,"Popsfile NIfFlag and calculated NIfFlag not same")
        if(.not.tPopsMapping) then
            if(PopNIfTot.ne.NIfTot) call stop_all(this_routine,"Popsfile NIfTot and calculated NIfTot not same")
        else
            MappingNIfTot=PopNIfTot
        endif


        IF(.not.tWalkContGrow) THEN
!If we want the walker number to be stable, take the shift from the POPSFILE, otherwise, keep the input value.
            DiagSft=PopDiagSft
        ENDIF

        if(PopDiagSft.eq.0.0_dp) then
            !If the popsfile has a shift of zero, continue letting the population grow
            tWalkContGrow=.true.
            DiagSft=PopDiagSft
        endif

        if(tWalkContGrow) then
            tSinglePartPhase=.true.
            !If continuing to grow, ensure we can allocate enough memory for what we hope to get the walker population to,
            !rather than the average number of determinants in the popsfile.
            WalkerListSize=int(max(int(initwalkers,int64),NINT(real(iPopAllTotWalkers,dp)/real(nNodes,dp),int64)),sizeof_int)
        else
            tSinglePartPhase=.false.
            WalkerListSize=NINT(real(iPopAllTotWalkers,dp)/real(nNodes,dp))
        endif

        AllSumNoatHF=PopSumNoatHF
        AllSumENum=PopAllSumENum
        PreviousCycles=iPopIter
        WRITE(6,*) "Number of iterations in previous simulation: ",PreviousCycles
        IF(NEquilSteps.gt.0) THEN
            WRITE(6,*) "Removing equilibration steps since reading in from POPSFILE."
            NEquilSteps=0
        ENDIF
        IF(TZeroProjE) THEN
            !Reset energy estimator
            WRITE(6,*) "Resetting projected energy counters to zero..."
            AllSumENum=0.0_dp
            AllSumNoatHF = 0
        ENDIF
        if(read_tau.eq.0.0_dp) then
            !Using popsfile v.3, where tau is not written out.
            !Exit if trying to dynamically search for timestep
            if(tSearchTau) then
                call stop_all(this_routine,"Cannot dynamically search for timestep if reading &
                    &in POPSFILE v.3. Manually specify timestep.")
            endif
            write(6,*) "Old popsfile detected."
            write(6,*) "Therefore automatic blocking will only start from current run"
            iBlockingIter = PreviousCycles
        else
            !Using popsfile v.4, where tau is written out and read in
            if(tSearchTau) then
                if(.not.tSinglePartPhase) then
                    tSearchTau=.false.
                endif
                Tau=read_tau
                write(6,"(A)") "Using timestep specified in POPSFILE, although continuing to dynamically adjust to optimise this"
            else
                !Tau specified. if it is different, write this here.
                if(abs(read_tau-Tau).gt.1.0e-5_dp) then
                    call warning_neci(this_routine,"Timestep specified in input file is different to that in the popsfile.")
                    write(6,"(A,F12.8)") "Old timestep: ",read_tau
                    write(6,"(A,F12.8)") "New timestep: ",tau
                endif
            endif
            iBlockingIter = PopBlockingIter
        endif
    
    end subroutine CheckPopsParams


!Routine for reading in from iunit the header information from a popsile v3 file.
    subroutine ReadPopsHeadv3(iunithead,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot)
        integer , intent(in) :: iunithead
        logical , intent(out) :: tPop64Bit,tPopHPHF,tPopLz
        integer , intent(out) :: iPopLenof_sign,iPopNel,iPopIter,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot
        integer(int64) , intent(out) :: iPopAllTotWalkers
        real(dp) , intent(out) :: PopDiagSft
        real(dp) , dimension(lenof_sign) , intent(out) :: PopSumNoatHF
        HElement_t , intent(out) :: PopAllSumENum
        character(len=24) :: junk,junk2,junk3,junk4,junk5
        integer :: PopsVersion

        PopsVersion=FindPopsfileVersion(iunithead)
        if(PopsVersion.ne.3) call stop_all("ReadPopsfileHeadv3","Wrong popsfile version for this routine.")
            
        if(iProcIndex.eq.root) then
            read(iunithead,'(A12,L5,A8,L5,A8,L5,A13,I5,A7,I6)') junk,tPop64Bit,junk2,tPopHPHF,junk3, &
                tPopLz,junk4,iPopLenof_sign,junk5,iPopNEl
            read(iunithead,*) iPopAllTotWalkers
            read(iunithead,*) PopDiagSft 
            read(iunithead,*) PopSumNoatHF 
            read(iunithead,*) PopAllSumENum 
            read(iunithead,*) iPopIter 
            read(iunithead,*) PopNIfD 
            read(iunithead,*) PopNIfY 
            read(iunithead,*) PopNIfSgn 
            read(iunithead,*) PopNIfFlag 
            read(iunithead,*) PopNIfTot
        endif
        !Broadcast the read in values from the header to all nodes.
        call MPIBCast(tPop64Bit)
        call MPIBCast(tPopHPHF)
        call MPIBCast(tPopLz)
        call MPIBCast(iPopLenof_sign)
        call MPIBCast(iPopNEl)
        call MPIBCast(iPopAllTotWalkers)
        call MPIBCast(PopDiagSft)
        call MPIBCast(PopSumNoatHF)
        call MPIBCast(PopAllSumENum)
        call MPIBCast(iPopIter)
        call MPIBCast(PopNIfD)
        call MPIBCast(PopNIfY)
        call MPIBCast(PopNIfSgn)
        call MPIBCast(PopNIfFlag)
        call MPIBCast(PopNIfTot)

    end subroutine ReadPopsHeadv3

    !Routine for reading in from iunit the header information from a popsile v4 file.
    subroutine ReadPopsHeadv4(iunithead,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot,read_tau,PopBlockingIter)
        integer , intent(in) :: iunithead
        logical , intent(out) :: tPop64Bit,tPopHPHF,tPopLz
        integer , intent(out) :: iPopLenof_sign,iPopNel,iPopIter,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot
        integer , intent(out) :: PopBlockingIter
        integer(int64) , intent(out) :: iPopAllTotWalkers
        real(dp) , intent(out) :: PopDiagSft,read_tau
        real(dp) , dimension(lenof_sign) , intent(out) :: PopSumNoatHF
        HElement_t , intent(out) :: PopAllSumENum
        integer :: PopsVersion
        !Variables for the namelist
        logical :: Pop64Bit,PopHPHF,PopLz
        integer :: PopLensign,PopNEl,PopCyc,PopiBlockingIter
        integer(int64) :: PopTotwalk
        real(dp) :: PopSft,PopTau
        HElement_t :: PopSumENum
        namelist /POPSHEAD/ Pop64Bit,PopHPHF,PopLz,PopLensign,PopNEl,PopTotwalk,PopSft,PopSumNoatHF,PopSumENum, &
                    PopCyc,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot,PopTau,PopiBlockingIter

        PopsVersion=FindPopsfileVersion(iunithead)
        if(PopsVersion.ne.4) call stop_all("ReadPopsfileHeadv4","Wrong popsfile version for this routine.")
            
        if(iProcIndex.eq.root) then
            read(iunithead,POPSHEAD)
        endif
        !Broadcast the read in values from the header to all nodes.
        call MPIBCast(Pop64Bit)
        call MPIBCast(PopHPHF)
        call MPIBCast(PopLz)
        call MPIBCast(PopLensign)
        call MPIBCast(PopNEl)
        call MPIBCast(PopTotwalk)
        call MPIBCast(PopSft)
        call MPIBCast(PopSumNoatHF)
        call MPIBCast(PopSumENum)
        call MPIBCast(PopCyc)
        call MPIBCast(PopNIfD)
        call MPIBCast(PopNIfY)
        call MPIBCast(PopNIfSgn)
        call MPIBCast(PopNIfFlag)
        call MPIBCast(PopNIfTot)
        call MPIBCast(PopTau)
        call MPIBCast(PopiBlockingIter)
        tPop64Bit=Pop64Bit
        tPopHPHF=PopHPHF
        tPopLz=PopLz
        iPopLenof_sign=PopLensign
        iPopNel=PopNel
        iPopAllTotWalkers=PopTotwalk
        PopDiagSft=PopSft
        PopAllSumENum=PopSumENum
        iPopIter=PopCyc
        read_tau=PopTau 
        PopBlockingIter=PopiBlockingIter

    end subroutine ReadPopsHeadv4
    
    !NOTE: This should only be used for the v3 POPSFILEs, since we only open the POPSFILE on the head node.
    subroutine open_pops_head(iunithead,formpops,binpops)
        integer , intent(out) :: iunithead
        logical , intent(out) :: formpops,binpops
        character(255) :: popsfile

        if(iProcIndex.eq.root) then
            iunithead=get_free_unit()
            call get_unique_filename('POPSFILE',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            inquire(file=popsfile,exist=formpops)
            if(formpops) then
                open(iunithead,file=popsfile,status='old')
                binpops=.false.
            else
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                inquire(file=popsfile,exist=binpops)
                if(binpops) then
                    call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                    open(iunithead,file=popsfile,status='old')
                else 
                    call stop_all("open_pops_head","No POPSFILEs detected...")
                endif
            endif
            rewind(iunithead)
        endif
        call MPIBCast(binpops) 
        call MPIBCast(formpops) 

    end subroutine open_pops_head

!Return the version number of the popsfile
    integer function FindPopsfileVersion(iunithead)
        integer, intent(in) :: iunithead
        logical :: formpops,binpops
        character(255) :: FirstLine

        if(iProcIndex.eq.root) then
            rewind(iunithead)
            read(iunithead,'(a255)') FirstLine

            if(index(FirstLine,'VERSION').eq.0) then
                FindPopsfileVersion=1
            else
                rewind(iunithead)
                read(iunithead,*) FirstLine,FirstLine,FirstLine,FindPopsfileVersion
            endif
        endif
        call MPIBCast(FindPopsfileVersion)

    end function FindPopsfileVersion


!This routine is the same as WriteToPopsfilePar, but does not require two 
! main arrays to hold the data.
!The root processors data will be stored in a temporary array while it 
! recieves the data from the other processors.
!This routine will write out to a popsfile. It transfers all walkers to the 
! head node sequentially, so does not want to be called too often
    SUBROUTINE WriteToPopsfileParOneArr(Dets,nDets)
        use CalcData, only: iPopsFileNoWrite
        use constants, only: size_n_int,n_int
        use MemoryManager, only: TagIntType
        integer(int64),intent(in) :: nDets !The number of occupied entries in Dets
        integer(kind=n_int),intent(in) :: Dets(0:nIfTot,1:nDets)
        INTEGER :: error
        integer(int64) :: WalkersonNodes(0:nNodes-1),writeoutdet
        INTEGER :: Tag
        INTEGER :: Total,i,j,k
        INTEGER(KIND=n_int), ALLOCATABLE :: Parts(:,:)
        INTEGER(TagIntType) :: PartsTag=0
        INTEGER :: nMaxDets, TempDet(0:NIfTot), TempFlags
        integer :: iunit, iunit_2, Initiator_Count
        CHARACTER(len=*) , PARAMETER :: this_routine='WriteToPopsfileParOneArr'
        character(255) :: popsfile
        real(dp) :: TempSign(lenof_sign)
        integer :: excit_lev

        CALL MPIBarrier(error)  !sync
!        WRITE(6,*) "Get Here",nDets
!        CALL neci_flush(6)

!First, make sure we have up-to-date information - again collect AllTotWalkers
! ,AllSumNoatHF and AllSumENum...
!Calculate the energy by summing all on HF and doubles - convert number at HF
!  to a real since no int*8 MPI data type
        CALL MPISum(SumNoatHF,1,AllSumNoatHF)
        CALL MPISum(SumENum,1,AllSumENum)

!We also need to tell the root processor how many particles to expect from 
!each node - these are gathered into WalkersonNodes
        if(bNodeRoot) CALL MPIAllGather(nDets,WalkersonNodes,error,Roots)

        Tag=125

!We have to make the distinction here between the number of entries to expect,
!and the number of determinants we are writing out. Since the list is not
!necessarily contiguous any more, we have to calculate Alltotwalkers seperately.
        if(tHashwalkerlist) then
            Writeoutdet=0
            do i=1,int(nDets,sizeof_int)
                call extract_sign(Dets(:,i),TempSign)
                if(.not.IsUnoccDet(TempSign)) then
                    !Count this det in AllTotWalkers
                    Writeoutdet=Writeoutdet+1
                endif
            enddo
            writeoutdet=int(writeoutdet/iPopsPartEvery)
            call mpisum(writeoutdet,1,AllTotWalkers)
        else
            if(iProcIndex.eq.Root) then
                AllTotWalkers=0
                do i=0,nNodes-1
                    AllTotWalkers=AllTotWalkers+INT(WalkersonNodes(i)/iPopsPartEvery)
                enddo
            endif
        endif


        IF(iProcIndex.eq.root) THEN
!First, check that we are going to receive the correct number of particles...

!            IF(Total.ne.AllTotWalkers) THEN
!                CALL Stop_All("WriteToPopsfilePar","Not all walkers accounted for...")
!            ENDIF

!Write header information
#ifdef __INT64
            IF(iPopsPartEvery.ne.1) THEN
                IF(tBinPops) THEN
                    WRITE(6,"(A,I12,A)") "Writing a 64-bit binary reduced POPSFILEBIN, printing a total of ", &
                        AllTotWalkers, " determinants."
                ELSE
                    WRITE(6,"(A,I12,A)") "Writing a 64-bit reduced POPSFILE, printing a total of ", &
                        AllTotWalkers, " determinants."
                ENDIF
            ELSE
                IF(tBinPops) THEN
                    WRITE(6,*) "Writing to 64-bit binary POPSFILEBIN..."
                ELSE
                    WRITE(6,*) "Writing to 64-bit POPSFILE..."
                ENDIF
            ENDIF
#else
            IF(iPopsPartEvery.ne.1) THEN
                IF(tBinPops) THEN
                    WRITE(6,"(A,I12,A)") "Writing a binary reduced POPSFILEBIN, printing a total of ", &
                        AllTotWalkers, " determinants."
                ELSE
                    WRITE(6,"(A,I12,A)") "Writing a reduced POPSFILE, printing a total of ",AllTotWalkers, " determinants."
                ENDIF
            ELSE
                IF(tBinPops) THEN
                    WRITE(6,*) "Writing to binary POPSFILEBIN..."
                ELSE
                    WRITE(6,*) "Writing to POPSFILE..."
                ENDIF
            ENDIF
#endif
            IF(tBinPops) THEN
                call get_unique_filename('POPSFILEHEAD',tIncrementPops,.true.,iPopsFileNoWrite,popsfile)
            ELSE
                call get_unique_filename('POPSFILE',tIncrementPops,.true.,iPopsFileNoWrite,popsfile)
            ENDIF
            iunit = get_free_unit()
            OPEN(iunit,FILE=popsfile,Status='replace')
            WRITE(iunit,"(A)") "# POPSFILE VERSION 4"

!v.3 POPSFILE HEADER - now depreciated.
!#ifdef __INT64
!            WRITE(iunit,'(A12,L5,A8,L5,A8,L5,A13,I5,A7,I6)') '64BitDets=',.TRUE.,'HPHF=',tHPHF,'Lz=', &
!                tFixLz,'Lenof_sign=',lenof_sign,'NEl=',NEl
!#else
!            WRITE(iunit,'(A12,L5,A8,L5,A8,L5,A13,I5,A7,I6)') '64BitDets=',.FALSE.,'HPHF=',tHPHF,'Lz=', &
!                tFixLz,'Lenof_sign=',lenof_sign,'NEl=',NEl
!#endif
!            WRITE(iunit,*) AllTotWalkers,"   TOTWALKERS (all nodes)"
!            WRITE(iunit,*) DiagSft,"   DIAG SHIFT"
!            WRITE(iunit,*) AllSumNoatHF,"   SUMNOATHF (all nodes)"
!            WRITE(iunit,*) AllSumENum,"   SUMENUM ( \sum<D0|H|Psi> - all nodes)"
!            WRITE(iunit,*) Iter+PreviousCycles,"   PREVIOUS CYCLES"
!            WRITE(iunit,*) NIfD,"    NIfD"
!            WRITE(iunit,*) NIfY,"    NIfY"
!            WRITE(iunit,*) NIfSgn,"    NIfSgn"
!            WRITE(iunit,*) NIfFlag,"    NIfFlag"
!            WRITE(iunit,*) NIfTot,"    NIfTot"

#ifdef __INT64
            write(iunit,'(A)') '&POPSHEAD Pop64Bit=.TRUE.,'
#else
            write(iunit,'(A)') '&POPSHEAD Pop64Bit=.FALSE.,'
#endif
            write(iunit,'(A,L1,A,L1,A,I2,A,I3,A)') 'PopHPHF=',tHPHF,',PopLz=',tFixLz,',PopLensign=',lenof_sign,',PopNEl=',NEl,','
            write(iunit,'(A,I15,A,F18.12,A)') 'PopTotwalk=',AllTotWalkers,',PopSft=',DiagSft,','
            write(iunit,*) 'PopSumNoatHF=',AllSumNoatHF,','
            write(iunit,*) 'PopSumENum=',AllSumENum,','
            write(iunit,'(A,I16,A,I2,A,I2,A,I2,A)') 'PopCyc=',Iter+PreviousCycles,',PopNIfD=',  &
                                NIfD,',PopNIfY=',NIfY,',PopNIfSgn=',NIfSgn,','
            write(iunit,'(A,I2,A,I2,A,F18.12,A)') 'PopNIfFlag=',NIfFlag,',PopNIfTot=',NIfTot,',PopTau=',Tau,','
            write(iunit,'(A,I16)') 'PopiBlockingIter=',iBlockingIter
            write(iunit,'(A5)') '&END'

            IF(tBinPops) THEN
                CLOSE(iunit)
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.true.,iPopsFileNoWrite,popsfile)
                OPEN(iunit,FILE=popsfile,Status='replace',form='unformatted')
            ENDIF
            if(tPrintInitiators) then
                iunit_2 = get_free_unit()
                OPEN(iunit_2,FILE='INITIATORS',Status='UNKNOWN')
                Initiator_Count = 0
            endif

            ! Write out the dets from this node.
            do j = 1, int(ndets, sizeof_int)
                call write_pops_det (iunit, iunit_2, Dets(:,j), j)
            end do

!            WRITE(6,*) "Written out own walkers..."
!            write(6,*) WalkersOnNodes
!            CALL neci_flush(6)

!Now, we copy the head nodes data to a new array...
            nMaxDets=int(maxval(WalkersOnNodes),sizeof_int)
            ALLOCATE(Parts(0:NIfTot,nMaxDets),stat=error)
            CALL LogMemAlloc('Parts',int(nMaxDets,int32)*(NIfTot+1),size_n_int,this_routine,PartsTag,error)

!Now we need to receive the data from each other processor sequentially
!We can overwrite the head nodes information, since we have now stored it elsewhere.
            do i=1,nNodes-1
!Run through all other processors...receive the data...
                j=int(WalkersonNodes(i),sizeof_int)*(NIfTot+1)
                CALL MPIRecv(Parts(0:NIfTot,1:WalkersonNodes(i)),j,NodeRoots(i),Tag,error)
!                WRITE(6,*) "Recieved walkers for processor ",i,WalkersOnNodes(i)
!                CALL neci_flush(6)
                
                ! Then write it out...
                do j = 1, int(WalkersonNodes(i), sizeof_int)
                    call write_pops_det(iunit, iunit_2, Parts(:,j), j)
                end do

!                WRITE(6,*) "Writted out walkers for processor ",i
!                CALL neci_flush(6)

            enddo

            CLOSE(iunit)
            if (tPrintInitiators) CLOSE(iunit_2)

!Deallocate memory for temporary storage of information.
            DEALLOCATE(Parts)
            CALL LogMemDealloc(this_routine,PartsTag)

        ELSEIF(bNodeRoot) THEN
!All other processors need to send their data to root...
            j=int(nDets,sizeof_int)*(NIfTot+1)
            CALL MPISend(Dets(0:NIfTot,1:nDets),j,root,Tag,error)
!            WRITE(6,*) "Have sent info to head node..."
!            CALL neci_flush(6)
        ENDIF

!Reset the values of the global variables
        AllSumNoatHF = 0
        AllSumENum=0.0_dp
        AllTotWalkers=0
        RETURN

    END SUBROUTINE WriteToPopsfileParOneArr


    subroutine write_pops_det (iunit, iunit_2, det, j)

        ! Output a particle to a popsfile in format acceptable for popsfile v4

        integer, intent(in) :: iunit, iunit_2
        integer(n_int), intent(in) :: det(0:NIfTot)
        integer :: flg, j, k, tmp_flag, ex_level
        real(dp) :: real_sgn(lenof_sign)

        ! Note, we don't want to bother outputting empty particles
        call extract_sign(det, real_sgn)
        if (.not. IsUnoccDet(real_sgn)) then

            ! If we are using a reduced
            if (mod(j, iPopsPartEvery) == 0) then

                ! Write output in the desired format. If __INT64, we are 
                ! including the flag information with the signs in storage in
                ! memory --> need to extract these before outputting them.
                tmp_flag = extract_flags(det)
                if (tBinPops) then
                    write(iunit) det(0:NIfD), real_sgn
                    if (tUseFlags) write(iunit) int(flg, n_int)
                else
                    do k = 0, NIfD
                        write(iunit, '(i24)', advance='no') det(k)
                    end do
                    do k = 1, lenof_sign
                        write(iunit, '(f30.8)', advance='no') real_sgn(k)
                    end do
                    if (tUseFlags) write(iunit, '(i24)', advance='no') flg
                    write(iunit, *)
                end if

                if (tPrintInitiators .and. &
                    abs(real_sgn(1)) > InitiatorWalkNo) then
                    ! Testing using the sign now, because after annihilation
                    ! the current flag will not necessarily be correct.
                    ex_level = FindBitExcitLevel(ilutRef, det, nel)
                    write(iunit_2, '(f30.8,a20)', advance='no') &
                        abs(real_sgn(1)), ''
                    call writebitdet (iunit_2, det, .false.)
                    write(iunit_2, '(i30)') ex_level

                end if
            end if
        end if

    end subroutine




!This routine reads in particle configurations from a POPSFILE.
    SUBROUTINE ReadFromPopsfilePar()
        use CalcData , only : MemoryFacPart,MemoryFacAnnihil,MemoryFacSpawn,iWeightPopRead
        use Logging, only: tZeroProjE, tRDMonFly, tExplicitAllRDM, tHF_Ref_Explicit 
        use constants, only: size_n_int,bits_n_int
        LOGICAL :: exists,tBinRead
        Real(dp) :: NodeSumNoatHF(nProcessors)
        INTEGER :: WalkerstoReceive(nProcessors)
        integer(int64) :: AvWalkers
        integer(int64) :: TempTotParts(lenof_sign),TempCurrWalkers
        INTEGER :: TempInitWalkers,error,i,j,l,total,ierr,MemoryAlloc,Tag,Proc,CurrWalkers,ii
        INTEGER , DIMENSION(lenof_sign) :: TempSign
        real(dp) :: RealTempSign(lenof_sign)
        integer(int64) :: iLutTemp64(0:nBasis/64+1)
        INTEGER :: iLutTemp32(0:nBasis/32+1)
        INTEGER(KIND=n_int) :: iLutTemp(0:NIfTot)
        INTEGER :: AvSumNoatHF,IntegerPart,TempnI(NEl),ExcitLevel
        INTEGER :: NIfWriteOut,pos,orb,PopsVersion, iunit
        real(dp) :: r, FracPart, Gap, DiagSftTemp, tmp_dp
        HElement_t :: HElemTemp
        CHARACTER(len=*), PARAMETER :: this_routine='ReadFromPopsfilePar'
        character(255) :: popsfile,FirstLine
        character(len=24) :: junk,junk2,junk3,junk4
        LOGICAL :: tPop64BitDets,tPopHPHF,tPopLz,tPopInitiator
        integer(n_int) :: ilut_largest(0:NIfTot)
        real(dp) :: sign_largest

        IF(lenof_sign.ne.1) CALL Stop_All("ReadFromPopsfilePar","Popsfile V.2 does not work with complex walkers")
        
        PreviousCycles=0    !Zero previous cycles
        SumENum=0.0_dp
        TotParts=0
        SumNoatHF=0
        DiagSft=0.0_dp
        Tag=124             !Set Tag
        
        call get_unique_filename('POPSFILE',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
        iunit = get_free_unit()
        INQUIRE(FILE=popsfile,EXIST=exists)
        IF(exists) THEN
            OPEN(iunit,FILE=popsfile,Status='old')
            tBinRead=.false.
        ELSE
            tBinRead=.true.
            call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            INQUIRE(FILE=popsfile,EXIST=exists)
            IF(.not.exists) THEN
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"No POPSFILEs of any kind found.")
                ELSE
                    CALL Stop_All(this_routine,"POPSFILEBIN(.x) found, but POPSFILEHEAD(.x) also needed for header information")
                ENDIF
            ELSE
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"POPSFILEHEAD(.x) found, but no POPSFILEBIN(.x) " &
                    & //"for particle information - this is also needed")
                ELSE
                    call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                    OPEN(iunit,FILE=popsfile,Status='old')
                ENDIF
            ENDIF
        ENDIF

        READ(iunit,'(a255)') FirstLine

        IF(INDEX(FirstLine,'VERSION').eq.0) THEN
!No version number to be found
            PopsVersion=1
            REWIND(iunit)
        ELSE
            !Found version - which number is it?
            REWIND(iunit)
            READ(iunit,*) FirstLine,FirstLine,FirstLine,PopsVersion
        ENDIF
        WRITE(6,"(A,I5,A)") "Version",PopsVersion," POPSFILE detected"

!Read in initial data on processors which have a popsfile
        IF(PopsVersion.eq.2) THEN
            READ(iunit,'(A12,L5,A8,L5,A8,L5,A12,L5)') junk,tPop64BitDets,junk2,tPopHPHF,junk3,tPopLz,junk4,tPopInitiator
        ELSE
            WRITE(6,'(A)') "Reading in from depreciated POPSFILE - assuming that parameters " &
            & //"are the same as when POPSFILE was written"
        ENDIF
        READ(iunit,*) tmp_dp
        AllTotWalkers = int(tmp_dp,int64)
        READ(iunit,*) DiagSftTemp
        READ(iunit,*) AllSumNoatHF
        READ(iunit,*) AllSumENum
        READ(iunit,*) PreviousCycles

        IF(iProcIndex.eq.Root) THEN
            IF(iWeightPopRead.ne.0) THEN
                WRITE(6,"(A,I15,A,I4,A)") "Although ",AllTotWalkers, &
                 " configurations will be read in, only determinants with a weight of over ",iWeightPopRead," will be stored."
            ENDIF
        ENDIF

        IF(.not.tWalkContGrow) THEN
!If we want the walker number to continue growing, then take the diagonal 
! shift from the input, rather than the POPSFILE.
            DiagSft=DiagSftTemp
        ENDIF

        IF(DiagSftTemp.eq.0.0_dp) THEN
            tWalkContGrow=.true.
            DiagSft=DiagSftTemp
        ENDIF

        IF(tBinRead) THEN
!Test for the end of the file.
!If this is not the end of the file, there is one more keyword that tells us 
! the calculation had not entered variable shift mode yet.
!Want to put this test at the end of the non-binary file too.
            CLOSE(iunit)
            call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            OPEN(iunit,FILE=popsfile,Status='old',form='unformatted')
        ENDIF

        IF(iProcIndex.eq.Root) THEN

            WRITE(6,*) "Number of cycles in previous simulation: ",PreviousCycles
            IF(NEquilSteps.gt.0) THEN
                WRITE(6,*) "Removing equilibration steps since reading in from POPSFILE."
                NEquilSteps=0
            ENDIF
            IF(TZeroProjE) THEN
!Reset energy estimator
                WRITE(6,*) "Resetting projected energy counters to zero..."
                AllSumENum=0.0_dp
                AllSumNoatHF = 0
            ENDIF

!Need to calculate the number of walkers each node will receive...
            AvWalkers=NINT(AllTotWalkers/real(nProcessors,dp),int64)

!Divide up the walkers to receive for each node
            do i=1,nProcessors-1
                WalkerstoReceive(i)=int(AvWalkers,sizeof_int)
            enddo
!The last processor takes the 'remainder'
            WalkerstoReceive(nProcessors)=int(AllTotWalkers-(AvWalkers*(nProcessors-1)),sizeof_int)

!Quick check to ensure we have all walkers accounted for
            total=0
            do i=1,nProcessors
                total=total+WalkerstoReceive(i)
            enddo
            if (total /= AllTotWalkers) then
                CALL Stop_All("ReadFromPopsfilePar","All Walkers not accounted for when reading in from POPSFILE")
            endif
            
!InitWalkers needs to be reset for the culling criteria
            IF(.not.tWalkContGrow) THEN
!Now, let the total space allocated for storing walkers which have been read in to be equal to the initwalkers from the input file.
!                InitWalkers=AvWalkers
            ELSE
                TSinglePartPhase=.true.
            ENDIF
            SumENum=AllSumENum/REAL(nProcessors,dp)     !Divide up the SumENum over all processors
            AvSumNoatHF = int(AllSumNoatHF(1)/nProcessors,sizeof_int) !This is the average Sumnoathf
            do i=1,nProcessors-1
                NodeSumNoatHF(i)=real(INT(AvSumNoatHF,int64),dp)
            enddo
            NodeSumNoatHF(nProcessors)=real(AllSumNoatHF(1)-INT((AvSumNoatHF*(nProcessors-1)),int64),dp)

            ProjectionE=AllSumENum/real(AllSumNoatHF(1),dp)
                
!Reset the global variables
            AllSumENum=0.0_dp
            AllSumNoatHF = 0

        ENDIF

        CALL MPIBarrier(error)  !Sync

!Now we need to scatter the WalkerstoReceive to each node, and allocate the desired memory to each node...
!Broadcast info which needs to go to all processors
        CALL MPIBCast(DiagSft)
        CALL MPIBCast(SumENum)
        CALL MPIBCast(InitWalkers)
        CALL MPIBCast(NEquilSteps)
        CALL MPIBCast(NShiftEquilSteps)
        CALL MPIBCast(TSinglePartPhase)
!        CALL MPI_BCast(tChangenProcessors,1,MPI_LOGICAL,root,MPI_COMM_WORLD,error)
!Scatter the number of walkers each node will receive to TempInitWalkers, and 
! the SumNoatHF for each node which is distributed approximatly equally
        CALL MPIScatter(WalkerstoReceive,TempInitWalkers,error)
        CALL MPIScatter(NodeSumNoatHF,SumNoatHF(1),error)

        IF(MemoryFacPart.le.1.0_dp) THEN
            WRITE(6,*) 'MemoryFacPart must be larger than 1.0 when reading in a POPSFILE - increasing it to 1.50.'
            MemoryFacPart=1.50
        ENDIF
        
!Now we want to allocate memory on all nodes.
!InitWalkers here is simply the average number of walkers per node, not actual
        MaxWalkersPart=NINT(MemoryFacPart*(NINT(InitWalkers*ScaleWalkers)))
        MaxSpawned=NINT(MemoryFacSpawn*(NINT(InitWalkers*ScaleWalkers)))

        Gap=REAL(MaxSpawned)/REAL(nProcessors)
        do i=0,nProcessors-1
            InitialSpawnedSlots(i)=NINT(Gap*i)+1
        enddo
!ValidSpawndList now holds the next free position in the newly-spawned list, 
! but for each processor.
        ValidSpawnedList(:)=InitialSpawnedSlots(:)

        CALL MPIBarrier(error)  !Sync

!Allocate memory to hold walkers at least temporarily
        ALLOCATE(WalkVecDets(0:NIfTot,MaxWalkersPart),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating WalkVecDets array.')
        CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfTot+1),size_n_int,this_routine,WalkVecDetsTag,ierr)
        WalkVecDets(:,:)=0
        MemoryAlloc=(NIfTot+1)*MaxWalkersPart*size_n_int    !Memory Allocated in bytes

!Just allocating this here, so that the SpawnParts arrays can be used for sorting the determinants when using direct annihilation.
        WRITE(6,"(A,I12,A)") " Spawning vectors allowing for a total of ",MaxSpawned, &
         &" particles to be spawned in any one iteration."
        ALLOCATE(SpawnVec(0:NIfTot,MaxSpawned),stat=ierr)
        CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfTot+1),size_n_int,this_routine,SpawnVecTag,ierr)
        SpawnVec(:,:)=0
        ALLOCATE(SpawnVec2(0:NIfTot,MaxSpawned),stat=ierr)
        CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfTot+1),size_n_int,this_routine,SpawnVec2Tag,ierr)
        SpawnVec2(:,:)=0
        MemoryAlloc=MemoryAlloc+MaxSpawned*(1+NIfTot)*2*size_n_int

!Point at correct spawning arrays
        SpawnedParts=>SpawnVec
        SpawnedParts2=>SpawnVec2
!Allocate pointer to the correct walker array...
        CurrentDets=>WalkVecDets

!Need to now allocate other arrays
        IF(.not.tRegenDiagHEls) THEN
            if(tRDMonFly.and.(.not.tExplicitAllRDM)) then
                ALLOCATE(WalkVecH(3,MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVecH',3*MaxWalkersPart,8,this_routine,WalkVecHTag,ierr)
                WalkVecH(:,:)=0.0_dp
                MemoryAlloc=MemoryAlloc+8*MaxWalkersPart*3
                WRITE(6,"(A)") " The current signs before death will be store for use in the RDMs."
                WRITE(6,"(A,F14.6,A)") " This requires ", REAL(MaxWalkersPart*8*3,dp)/1048576.0_dp," Mb/Processor"
                NCurrH = 3
            else
                ALLOCATE(WalkVecH(1,MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVecH',MaxWalkersPart,8,this_routine,WalkVecHTag,ierr)
                WalkVecH(:,:)=0.0_dp
                MemoryAlloc=MemoryAlloc+8*MaxWalkersPart
                NCurrH = 1
            endif
        ELSE
            WRITE(6,"(A,F14.6,A)") " Diagonal H-Elements will not be stored. This will *save* ", &
                & REAL(MaxWalkersPart*8,dp)/1048576.0_dp," Mb/Processor"
        ENDIF

        if(tRDMonFly.and.(.not.tExplicitAllRDM).and.(.not.tHF_Ref_Explicit)) then
!Allocate memory to hold walkers spawned from one determinant at a time.
!Walkers are temporarily stored here, so we can check if we're spawning onto the same Dj multiple times.
            ALLOCATE(TempSpawnedParts(0:NIfDBO,20000),stat=ierr)
            CALL LogMemAlloc('TempSpawnedParts',20000*(NIfDBO+1),size_n_int,this_routine,TempSpawnedPartsTag,ierr)
            TempSpawnedParts(0:NIfDBO,1:20000)=0
            MemoryAlloc=MemoryAlloc + (NIfDBO+1)*20000*size_n_int    !Memory Allocated in bytes
            WRITE(6,"(A)") " Allocating temporary array for walkers spawned from a particular Di."
            WRITE(6,"(A,F14.6,A)") " This requires ", REAL(((NIfDBO+1)*20000*size_n_int),dp)/1048576.0_dp," Mb/Processor"
        endif

        IF(.not.tRegenDiagHEls) THEN
            CurrentH=>WalkVecH
        ENDIF

! The hashing will be different in the new calculation from the one where the
!  POPSFILE was produced, this means we must recalculate the processor each 
! determinant wants to go to.
! This is done by reading in all walkers to the root and then distributing 
! them in the same way as the spawning steps are done - by finding the
!  determinant and sending it there.
        IF((PopsVersion.ne.1).and.tHPHF.and.(.not.tPopHPHF)) THEN
            CALL Stop_All(this_routine,"HPHF on, but HPHF was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.tFixLz.and.(.not.tPopLz)) THEN
            CALL Stop_All(this_routine,"Lz on, but Lz was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tHPHF).and.tPopHPHF) THEN
            CALL Stop_All(this_routine,"HPHF off, but HPHF was used for creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tFixLz).and.tPopLz) THEN
            CALL Stop_All(this_routine,"Lz off, but Lz was used for creation of the POPSFILE")
        ENDIF
        ! TODO: Add tests for CSFs here.
        IF(PopsVersion.eq.1) THEN
            tPop64BitDets=.false.
            NIfWriteOut=nBasis/32
            IF(tCSF) NIfWriteOut=NIfWriteOut+1
        ELSE
            IF(.not.tPop64BitDets) THEN
                NIfWriteOut=nBasis/32
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ELSE
                NIfWriteOut=nBasis/64
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ENDIF
        ENDIF

        CurrWalkers=0
        sign_largest = 0
        ilut_largest = 0
        do i=1,int(AllTotWalkers,sizeof_int)
            iLutTemp(:)=0
            IF(PopsVersion.ne.1) THEN
                IF(tBinRead) THEN
                    IF(tPop64BitDets) THEN
                        READ(iunit) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ELSE
                    IF(tPop64BitDets) THEN
                        READ(iunit,*) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ENDIF
            ELSE
                !POPSFILE v. 1 only printed out 32 bit determinant strings.
                IF(tBinRead) THEN
                    READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                ELSE
                    READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                ENDIF
            ENDIF
            RealTempSign = transfer(TempSign, RealTempSign)

#ifdef __INT64
            if (.not.tPop64BitDets) then
                ! If we are using 64 bit integers, but have read in 32 bit 
                ! integers, then we need to convert them.
                do ii=0,nBasis/32
                    do j=0,31
                        if(btest(iLutTemp32(ii),j)) then
                            orb=(ii*32)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp64(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp64(0:NIfDBO)
            endif

#else
            ! If we are using 32 bit integers, but have read in 64 bit 
            ! integers, then we need to convert them.
            if (tPop64BitDets) then
                do ii=0,nBasis/64
                    do j=0,63
                        if(btest(iLutTemp64(ii),j)) then
                            orb=(ii*64)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp32(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp32(0:NIfDBO)
            endif
        
#endif
            call decode_bit_det (TempnI, iLutTemp)
            Proc = DetermineDetNode(TempnI,0)
            IF((Proc.eq.iNodeIndex).and.(abs(RealTempSign(1)).ge.iWeightPopRead)) THEN
                CurrWalkers=CurrWalkers+1
                !Do not need to send a flag here...
                call encode_bit_rep(CurrentDets(:,CurrWalkers),iLutTemp(0:NIfDBO),RealTempSign,0) 
                !TODO: Add flag for complex walkers to read in both
            ENDIF

            ! Keep track of what the most highly weighted determinant is
            if (abs(RealTempSign(1)) > sign_largest) then
                sign_largest = abs(RealTempSign(1))
                ilut_largest = iLutTemp
            endif
        enddo
        CLOSE(iunit)
        TempCurrWalkers=int(CurrWalkers,int64)

        ! Sort the lists so that they are in order if we change the number
        ! of processors.
        call sort (currentdets(:,1:CurrWalkers))

        ! Check that the bit-det comparisons agree that it is in order.
        do i=2,currwalkers
            if(DetBitLT(CurrentDets(:,i),CurrentDets(:,i-1),NIfDBO) == 1) then
                print*, 'Walkers: ', i-1, i
                print*, 'bit reps: '
                print*, currentdets(:, i-1)
                print*, currentdets(:, i)
                call stop_all (this_routine, 'Main list out of order')
            endif
        enddo

        CALL MPIBarrier(error)  !Sync
        CALL MPIAllReduce(TempCurrWalkers,MPI_SUM,AllTotWalkers)

        IF(iProcIndex.eq.root) WRITE(6,'(I10,A)') INT(AllTotWalkers,int64)," configurations read in from POPSFILE and distributed."

        IF(ScaleWalkers.ne.1) THEN

            WRITE(6,*) "Rescaling walkers by a factor of: ",ScaleWalkers

! CurrWalkers is the number of determinants on a particular node, AllTotWalkers is the total over all nodes.
            IntegerPart=INT(ScaleWalkers)
            FracPart=ScaleWalkers-REAL(IntegerPart)

            do l=1,CurrWalkers
                call extract_sign(CurrentDets(:,l),RealTempSign)
                RealTempSign=RealTempSign*IntegerPart
                r = genrand_real2_dSFMT() 
                IF(r.lt.FracPart) THEN
!Stochastically create another particle
                    IF(RealTempSign(1).lt.0) THEN
                        RealTempSign(1)=RealTempSign(1)-1
                    ELSE
                        RealTempSign(1)=RealTempSign(1)+1
                    ENDIF
                ENDIF
                call encode_sign(CurrentDets(:,l),RealTempSign)
            enddo

            InitWalkers=NINT(InitWalkers*ScaleWalkers)  !New (average) number of initial particles for culling criteria
!Other parameters don't change (I think) because the number of determinants isn't changing.                
            TotWalkers=CurrWalkers
            TotWalkersOld=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                AllTotWalkersOld=AllTotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ELSE
!We are not scaling the number of walkers...

            TotWalkers=CurrWalkers
            TotWalkersOld=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                AllTotWalkersOld=AllTotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ENDIF
            
        WRITE(6,*) "Initial Diagonal Shift (ECorr guess) is now: ",DiagSft
        WRITE(6,"(A,F14.6,A)") " Initial memory (without excitgens) consists of : ",REAL(MemoryAlloc,dp)/1048576.0_dp," Mb"
        WRITE(6,*) "Initial memory allocation successful..."
        WRITE(6,*) "Excitgens will be regenerated when they are needed..."
        CALL neci_flush(6)

        ! If we are changing the reference determinant to the largest
        ! weighted one in the file, do it here
        if (tReadPopsChangeRef .or. tReadPopsRestart) then
            if (.not. DetBitEq(ilut_largest, iLutRef, NIfDBO)) then
                
                ! Set new reference
                iLutRef = ilut_largest
                call decode_bit_det (ProjEDet, iLutRef)
                tNoBrillouin = .true.

                ! Recalculate the reference E
                if (tHPHF) then
                    HElemTemp = hphf_diag_helement (ProjEDet, iLutRef)
                elseif(tMomInv) then
                    HElemTemp = MI_diag_helement(ProjEDet,iLutRef)
                else
                    HElemTemp = get_helement (ProjEDet, ProjEDet, 0)
                endif
                Hii = real(HElemTemp, dp)

                ! Output info on root node.
                if (iProcIndex == root) then
                    write(6, '(a)', advance='no') &
                        "Changing projected energy reference determinant to: "
                    call write_det (6, ProjEDet, .true.)
                    write (6, '(a)') &
                        "Ensuring that Brillouin's theorem is no longer used."
                    write (6, '(a,g25.15)') &
                        "Reference energy now set to: ", Hii
                endif
            endif
        endif

!Now find out the data needed for the particles which have been read in...
        TotParts=0
        do j=1,int(TotWalkers,sizeof_int)
            call decode_bit_det (TempnI, currentDets(:,j))
            Excitlevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j), 2)
            IF(Excitlevel.eq.0) THEN
                IF(.not.tRegenDiagHEls) CurrentH(1,j)=0.0_dp
            ELSE
                IF(.not.tRegenDiagHEls) THEN
                    if (tHPHF) then
                        HElemTemp = hphf_diag_helement (TempnI, &
                                                        CurrentDets(:,j))
                    elseif(tMomInv) then
                        HElemTemp = MI_diag_helement(TempnI,CurrentDets(:,j))
                    else
                        HElemTemp = get_helement (TempnI, TempnI, 0)
                    endif
                    CurrentH(1,j)=REAL(HElemTemp,dp)-Hii
                ENDIF

            ENDIF
            call extract_sign(CurrentDets(:,j),RealTempSign)
            TotParts=TotParts+abs(RealTempSign(1))

        enddo

        CALL MPIBarrier(error)  !Sync
        CALL MPIReduce(TotParts,MPI_SUM,AllTotParts)

        IF(iProcIndex.eq.root) then
            AllTotPartsOld=AllTotParts
            iter_data_fciqmc%tot_parts_old = AllTotParts
        endif
        write(6,'(A,i20)') ' The total number of particles read from the POPSFILE is: ',AllTotParts(1)

        if (tReadPopsRestart) then
            tPopsAlreadyRead = .true.
            call ChangeRefDet (ProjEDet)
            tPopsAlreadyRead = .false.
        endif

    END SUBROUTINE ReadFromPopsfilePar

!This routine reads in particle configurations from a POPSFILE.
! It's a bastardisation of ReadFromPopsfilePar but only reads a popsfile into an array and leaves the rest until later.
! nDets comes in with the max number of amplitudes possible in Dets.
! Dets is then filled with all the amplitudes in the POPSFILE and nDets is returned with that number.

! This will very likely still be tweaking the inner heart-strings of FciMCPar.  Caveatis stulti!  AJWT
! GHB says he will incorporate this functionality into a rewrite of ReadFromPopsfilePar. 19/8/2010
    SUBROUTINE ReadFromPopsfileOnly(Dets,nDets)
        use CalcData , only : MemoryFacPart,MemoryFacAnnihil,MemoryFacSpawn,iWeightPopRead
        use Logging, only: tZeroProjE
        use constants, only: size_n_int,bits_n_int
        integer(int64),intent(inout) :: nDets !The number of occupied entries in Dets
        integer(kind=n_int),intent(out) :: Dets(0:nIfTot,1:nDets)
        LOGICAL :: exists,tBinRead
        integer(int64) :: TempCurrWalkers
        REAL(dp) :: TempTotParts(lenof_sign)
        INTEGER :: TempInitWalkers,error,i,j,l,total,ierr,MemoryAlloc,Tag,Proc,CurrWalkers,ii
        INTEGER , DIMENSION(lenof_sign) :: TempSign
        REAL(dp) , DIMENSION(lenof_sign) :: RealTempSign
        integer(int64) :: iLutTemp64(0:nBasis/64+1)
        INTEGER :: iLutTemp32(0:nBasis/32+1)
        INTEGER(KIND=n_int) :: iLutTemp(0:NIfTot)
        INTEGER :: AvSumNoatHF,IntegerPart,TempnI(NEl),ExcitLevel
        INTEGER :: NIfWriteOut,pos,orb,PopsVersion, iunit
        real(dp) :: r,FracPart,Gap,DiagSftTemp
        HElement_t :: HElemTemp
        CHARACTER(len=*), PARAMETER :: this_routine='ReadFromPopsfilePar'
        character(255) :: popsfile,FirstLine
        character(len=24) :: junk,junk2,junk3,junk4
        LOGICAL :: tPop64BitDets,tPopHPHF,tPopLz,tPopInitiator
        integer(n_int) :: ilut_largest(0:NIfTot)
        real(dp) :: sign_largest

        MemoryAlloc = 0

        IF(lenof_sign.ne.1) CALL Stop_All("ReadFromPopsfilePar","Popsfile does not work with complex walkers")
        
        PreviousCycles=0    !Zero previous cycles
        SumENum=0.0_dp
        TotParts=0
        SumNoatHF=0
        DiagSft=0.0_dp
        Tag=124             !Set Tag
        
        call get_unique_filename('POPSFILE',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
        iunit = get_free_unit()
        INQUIRE(FILE=popsfile,EXIST=exists)
        IF(exists) THEN
            OPEN(iunit,FILE=popsfile,Status='old')
            tBinRead=.false.
        ELSE
            tBinRead=.true.
            call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            INQUIRE(FILE=popsfile,EXIST=exists)
            IF(.not.exists) THEN
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"No POPSFILEs of any kind found.")
                ELSE
                    CALL Stop_All(this_routine,"POPSFILEBIN(.x) found, but POPSFILEHEAD(.x) also needed for header information")
                ENDIF
            ELSE
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"POPSFILEHEAD(.x) found, but no POPSFILEBIN(.x) for " &
                    & //"particle information - this is also needed")
                ELSE
                    call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                    OPEN(iunit,FILE=popsfile,Status='old')
                ENDIF
            ENDIF
        ENDIF

        READ(iunit,'(a255)') FirstLine

        IF(INDEX(FirstLine,'VERSION').eq.0) THEN
!No version number to be found
            PopsVersion=1
            REWIND(iunit)
        ELSE
            !Found version - which number is it?
            REWIND(iunit)
            READ(iunit,*) FirstLine,FirstLine,FirstLine,PopsVersion
        ENDIF
        WRITE(6,"(A,I5,A)") "Version",PopsVersion," POPSFILE detected"



!Read in initial data on processors which have a popsfile
        IF(PopsVersion.eq.2) THEN
            READ(iunit,'(A12,L5,A8,L5,A8,L5,A12,L5)') junk,tPop64BitDets,junk2,tPopHPHF,junk3,tPopLz,junk4,tPopInitiator
        ELSE
            WRITE(6,'(A)') "Reading in from depreciated POPSFILE - assuming that parameters are same as when POPSFILE was written"
        ENDIF
        READ(iunit,*) AllTotWalkers
        READ(iunit,*) DiagSftTemp
        READ(iunit,*) AllSumNoatHF
        READ(iunit,*) AllSumENum
        READ(iunit,*) PreviousCycles

        IF(iProcIndex.eq.Root) THEN
            IF(iWeightPopRead.ne.0) THEN
                WRITE(6,"(A,I15,A,I4,A)") "Although ",AllTotWalkers,    &
                " configurations will be read in, only determinants with a weight of over ",iWeightPopRead," will be stored."
            ENDIF
        ENDIF

        IF(.not.tWalkContGrow) THEN
!If we want the walker number to continue growing, then take the diagonal shift from the input, rather than the POPSFILE.
            DiagSft=DiagSftTemp
        ENDIF

        IF(DiagSftTemp.eq.0.0_dp) THEN
            tWalkContGrow=.true.
            DiagSft=DiagSftTemp
        ENDIF

        IF(tBinRead) THEN
!Test for the end of the file.
!If this is not the end of the file, there is one more keyword that tells us 
! the calculation had not entered variable shift mode yet.
!Want to put this test at the end of the non-binary file too.
            CLOSE(iunit)
            call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            OPEN(iunit,FILE=popsfile,Status='old',form='unformatted')
        ENDIF

        IF(iProcIndex.eq.Root) THEN

            WRITE(6,*) "Number of cycles in previous simulation: ",PreviousCycles
            IF(NEquilSteps.gt.0) THEN
                WRITE(6,*) "Removing equilibration steps since reading in from POPSFILE."
                NEquilSteps=0
            ENDIF
            IF(TZeroProjE) THEN
!Reset energy estimator
                WRITE(6,*) "Resetting projected energy counters to zero..."
                AllSumENum=0.0_dp
                AllSumNoatHF = 0
            ENDIF

!InitWalkers needs to be reset for the culling criteria
            IF(.not.tWalkContGrow) THEN
!Now, let the total space allocated for storing walkers which have been read in to be equal to the initwalkers from the input file.
!                InitWalkers=AvWalkers
            ELSE
                TSinglePartPhase=.true.
            ENDIF
            SumENum=AllSumENum/REAL(nProcessors,dp)     !Divide up the SumENum over all processors
            AvSumNoatHF = int(AllSumNoatHF(1),sizeof_int)/nProcessors !This is the average Sumnoathf

            ProjectionE=AllSumENum/real(AllSumNoatHF(1),dp)
                
!Reset the global variables
            AllSumENum=0.0_dp
            AllSumNoatHF = 0

        ENDIF

        CALL MPIBarrier(error)  !Sync

!Now we need to scatter the WalkerstoReceive to each node, and allocate the desired memory to each node...
!Broadcast info which needs to go to all processors
        CALL MPIBCast(DiagSft)
        CALL MPIBCast(SumENum)
        CALL MPIBCast(InitWalkers)
        CALL MPIBCast(NEquilSteps)
        CALL MPIBCast(NShiftEquilSteps)
        CALL MPIBCast(TSinglePartPhase)
!        CALL MPI_BCast(tChangenProcessors,1,MPI_LOGICAL,root,MPI_COMM_WORLD,error)
!Scatter the number of walkers each node will receive to TempInitWalkers, 
! and the SumNoatHF for each node which is distributed approximatly equally

        IF(MemoryFacPart.le.1.0_dp) THEN
            WRITE(6,*) 'MemoryFacPart must be larger than 1.0 when reading in a POPSFILE - increasing it to 1.50.'
            MemoryFacPart=1.50
        ENDIF
        
        CALL MPIBarrier(error)  !Sync

        if(AllTotWalkers>nDets) CALL Stop_All(this_routine,'Not enough memory to read in POPSFILE.')

! The hashing will be different in the new calculation from the one where the
! POPSFILE was produced, this means we must recalculate the processor each 
! determinant wants to go to.                
! This is done by reading in all walkers to the root and then distributing 
! them in the same way as the spawning steps are done - by finding the 
! determinant and sending it there.
        IF((PopsVersion.ne.1).and.tHPHF.and.(.not.tPopHPHF)) THEN
            CALL Stop_All(this_routine,"HPHF on, but HPHF was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.tFixLz.and.(.not.tPopLz)) THEN
            CALL Stop_All(this_routine,"Lz on, but Lz was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tHPHF).and.tPopHPHF) THEN
            CALL Stop_All(this_routine,"HPHF off, but HPHF was used for creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tFixLz).and.tPopLz) THEN
            CALL Stop_All(this_routine,"Lz off, but Lz was used for creation of the POPSFILE")
        ENDIF
        ! TODO: Add tests for CSFs here.
        IF(PopsVersion.eq.1) THEN
            tPop64BitDets=.false.
            NIfWriteOut=nBasis/32
            IF(tCSF) NIfWriteOut=NIfWriteOut+1
        ELSE
            IF(.not.tPop64BitDets) THEN
                NIfWriteOut=nBasis/32
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ELSE
                NIfWriteOut=nBasis/64
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ENDIF
        ENDIF

        CurrWalkers=0
        sign_largest = 0
        ilut_largest = 0
        write(6,*) "Reading in ", AllTotWalkers, " walkers"
        do i=1,int(AllTotWalkers,sizeof_int)
            iLutTemp(:)=0
            IF(PopsVersion.ne.1) THEN
                IF(tBinRead) THEN
                    IF(tPop64BitDets) THEN
                        READ(iunit) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ELSE
                    IF(tPop64BitDets) THEN
                        READ(iunit,*) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ENDIF
            ELSE
                !POPSFILE v. 1 only printed out 32 bit determinant strings.
                IF(tBinRead) THEN
                    READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                ELSE
                    READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                ENDIF
            ENDIF
            RealTempSign = transfer(TempSign, RealTempSign)

#ifdef __INT64
            if (.not.tPop64BitDets) then
                ! If we are using 64 bit integers, but have read in 32 bit 
                ! integers, then we need to convert them.
                do ii=0,nBasis/32
                    do j=0,31
                        if(btest(iLutTemp32(ii),j)) then
                            orb=(ii*32)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp64(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp64(0:NIfDBO)
            endif

#else
            ! If we are using 32 bit integers, but have read in 64 bit 
            ! integers, then we need to convert them.
            if (tPop64BitDets) then
                do ii=0,nBasis/64
                    do j=0,63
                        if(btest(iLutTemp64(ii),j)) then
                            orb=(ii*64)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp32(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp32(0:NIfDBO)
            endif
        
#endif
            call decode_bit_det (TempnI, iLutTemp)
            Proc=0  !DetermineDetNode(TempnI)   !This wants to return a value between 0 -> nProcessors-1
            IF((Proc.eq.iProcIndex).and.(abs(RealTempSign(1)).ge.iWeightPopRead)) THEN
                CurrWalkers=CurrWalkers+1
                !Do not need to send a flag here...
                call encode_bit_rep(Dets(:,CurrWalkers),iLutTemp(0:NIfDBO),RealTempSign,0)
                !TODO: Add flag for complex walkers to read in both
            ENDIF

            ! Keep track of what the most highly weighted determinant is
            if (abs(RealTempSign(1)) > sign_largest) then
                sign_largest = abs(RealTempSign(1))
                ilut_largest = iLutTemp
            endif
        enddo
        CLOSE(iunit)
        TempCurrWalkers=int(CurrWalkers,int64)

        ! Sort the lists so that they are in order if we change the number
        ! of processors.
        call sort(Dets(:,1:CurrWalkers))

        ! Check that the bit-det comparisons agree that it is in order.
        do i=2,currwalkers
            if(DetBitLT(Dets(:,i),Dets(:,i-1),NIfDBO) == 1) then
                print*, 'Walkers: ', i-1, i
                print*, 'bit reps: '
                print*, dets(:, i-1)
                print*, dets(:, i)
                call stop_all (this_routine, 'Main list out of order')
            endif
        enddo

        CALL MPIBarrier(error)  !Sync
        CALL MPIAllReduce(TempCurrWalkers,MPI_SUM,AllTotWalkers)

        IF(iProcIndex.eq.root) WRITE(6,'(I10,A)') INT(AllTotWalkers,int64)," configurations read in from POPSFILE and distributed."

        IF(ScaleWalkers.ne.1) THEN

            WRITE(6,*) "Rescaling walkers by a factor of: ",ScaleWalkers

! CurrWalkers is the number of determinants on a particular node, AllTotWalkers is the total over all nodes.
            IntegerPart=INT(ScaleWalkers)
            FracPart=ScaleWalkers-REAL(IntegerPart)

            do l=1,CurrWalkers
                call extract_sign(Dets(:,l),RealTempSign)
                RealTempSign=TempSign*IntegerPart
                r = genrand_real2_dSFMT() 
                IF(r.lt.FracPart) THEN
!Stochastically create another particle
                    IF(RealTempSign(1).lt.0) THEN
                        RealTempSign(1)=RealTempSign(1)-1
                    ELSE
                        RealTempSign(1)=RealTempSign(1)+1
                    ENDIF
                ENDIF
                call encode_sign(Dets(:,l),RealTempSign)
            enddo

            InitWalkers=NINT(InitWalkers*ScaleWalkers)  !New (average) number of initial particles for culling criteria
!Other parameters don't change (I think) because the number of determinants isn't changing.                
            nDets=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ELSE
!We are not scaling the number of walkers...

            nDets=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ENDIF
            
        WRITE(6,*) "Initial Diagonal Shift (ECorr guess) is now: ",DiagSft
        WRITE(6,"(A,F14.6,A)") " Initial memory (without excitgens) consists of : ",REAL(MemoryAlloc,dp)/1048576.0_dp," Mb"
        WRITE(6,*) "Initial memory allocation successful..."
        WRITE(6,*) "Excitgens will be regenerated when they are needed..."
        CALL neci_flush(6)

!Now find out the data needed for the particles which have been read in...
        TotParts=0
        do j=1,int(nDets,sizeof_int)
            call extract_sign(Dets(:,j),RealTempSign)
            TotParts=TotParts+abs(RealTempSign(1))
        enddo

        TempTotParts=TotParts

        CALL MPIBarrier(error)  !Sync
        CALL MPIReduce(TempTotParts,MPI_SUM,AllTotParts)

        IF(iProcIndex.eq.root) AllTotPartsOld=AllTotParts
        write(6,'(A,i20)') ' The total number of particles read from the POPSFILE is: ',AllTotParts(1)

        if (tReadPopsRestart) then
            tPopsAlreadyRead = .true.
            call ChangeRefDet (ProjEDet)
            tPopsAlreadyRead = .false.
        endif

    END SUBROUTINE ReadFromPopsfileOnly

END MODULE PopsfileMod
