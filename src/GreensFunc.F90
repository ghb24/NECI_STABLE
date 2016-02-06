#include "macros.h"
!TODO:  Symmetry in storage (more difficult if not inverse)
!       Even better would be a mapping function, so we can choose which orbitals to consider

!       Writing out walker distridutions to disk should be done with a different filename, and
!       deleted at the end

!       If all particles die in any of the runs, then the overlap is 0, and come right out

!       Remove reading a writing to disk - we need to store both lists anyway

module GreensFuncMod
    use SystemData, only: nel, nBasis,tHPHF,tFixLz, G1,tNoBrillouin,tUseBrillouin
    use SystemData, only: nBasisMax,tKPntSym,ElecPairs,nOccAlpha,nOccBeta
    use global_utilities
    use constants, only: dp, int64, n_int
    use FciMCParMod, only: DeallocFCIMCMemPar,SetupParameters,InitFCIMCCalcPar
    use FciMCParMod, only: PerformFciMCycPar,calculate_new_shift_wrapper,BinSearchParts3
    use FciMCData, only: tRunningGF,tBackwardsTime,tSearchTau,CurrentDets,TotWalkers
    use FciMCData, only: ProjEDet,iLutRef,Hii,iter_data_fciqmc,TotParts
    use FciMCData, only: tRestart,AllTotParts,AllTotWalkers,tHashWalkerlist,HashIndex
    use FciMCData, only: ProjectionE,tSaveHashList,SpawnedParts,SpawnedParts2,TotPartsOld
    use FciMCData, only: ValidSpawnedList,InitialSpawnedSlots,HFDet,ProjEDet,ilutHF,iLutRef
    use FciMCData, only: TotWalkersOld,AllTotWalkersOld,NoatHF,AllNoatHF,AllTotPartsOld
    use FciMCData, only: OldAllNoatHF,OldAllAvWalkersCyc,OldAllHFCyc,CurrentH,nClashMax
    use FciMCData, only: nWalkerHashes,MaxSpawned,HolesInList
    use FciMCData, only: HFDetTag,WalkVecDetsTag,WalkVecHTag,SpawnVecTag,SpawnVec2Tag
    use FciMCData, only: ptr_excit_generator,ptr_attempt_create,ptr_get_spawn_helement
    use FciMCData, only: ptr_new_child_stats,ptr_encode_child,ptr_attempt_die,ptr_iter_data
    use CalcData, only: tReadPops,tChangeProjEDet,tRestartHighPop
    use CalcData, only: tau,DiagSft,SftDamp,StepsSft
    use bit_reps, only: NIfD,NIfDBO,NIfTot,extract_sign,decode_bit_det,decode_bit_det_nel
    use DetBitOps, only: DetBitEQ
    use sym_mod, only: GetSym,GetLz
    use Logging, only: tPopsfile,tIncrementPops,GFDebug
    use AnnihilationMod, only: IsUnoccDet,DetermineDetNode,FindWalkerHash,DetermineDetNode_nel
    use DeterminantData, only: FDet
    use PopsfileMod, only: WriteToPopsfileParOneArr,readpopsheadv4
    use GenRandSymExcitNUMod , only: RandExcitSymLabelProd
    use SymExcitDataMod, only: SpinOrbSymLabel
    use util_mod, only: get_unique_filename,get_free_unit
    use Parallel_neci
    use sort_mod
    use GreensFuncData
    use MEMmod, only : MEM

    implicit none

    contains

    !Call this after appropriate equilibration time
    subroutine CalcGreensFunc()
    implicit none
    integer :: ierr,Ex(2,1),sample,EquilIter,i,j,k,GFType,MinInd
    integer(int64) :: norm, allnorm
    integer :: MaxIndex,PartInd,neltemp,StartInd
    integer, dimension(2) :: GreenTimeSlice
    real(dp) :: OrigRefEnergy,rat,SftDampInput
    integer, allocatable :: OrigRefDet(:),nI_GF(:)
    integer(n_int) :: OrigRefiLut(0:NIfTot),TempDet(0:NIfTot)
    integer :: TimeNumber,iter,walk,proc,nI(nel),PairIndex,Pairs,jval
    integer, dimension(lenof_sign) :: SignCurr,SignRight,SignLeft
    logical :: tCommunicate,tAllCommunicate,tSuccess,tSign,tSkipDet,tSendDet
    logical :: toverride
    integer(int64) :: int64_tmp(1+lenof_sign),TotWalkersTemp,MaxWalkersProc
    real(dp) , allocatable :: G_ret_loc(:,:,:),G_adv_loc(:,:,:)
    real(dp) , allocatable :: G_adv_sq_loc(:,:,:),G_ret_sq_loc(:,:,:)
    character(len=35) :: abstr
    character(len=*), parameter :: this_routine='CalcGreensFunc'
    integer, allocatable :: UnitNos(:,:)

    !This block is just if we want to calculate the spectrum, without 
    !worrying about calculating the greens function
    if(tJustSpectrum) then
        write(6,"(A)") "Skipping calculation of greens function."
        if(tSpecificMatEls) then
            !Work out which j orbitals will actually be run over if tSpecificMatEls is set
            call SetupSpecificMatEls()
        endif
        jval=1
        lp: do j=1,nBasis
            if(tSpecificMatEls) then
                if(Allowed_j(jval).ne.j) cycle  !No desired i's here
                jval=jval+1
            endif
            do i=1,nBasis

                if(.not.CheckAllowed_i(j,i)) cycle   !Not an allowed i orbital with this j
                do GFType=1,2
                    if(tJustAdvGF.and.(GFType.eq.2)) cycle
                    if(tJustRetGF.and.(GFType.eq.1)) cycle

                    abstr=''
                    write(abstr,'(I3)') j
                    if(GFType.eq.1) then
                        abstr='Adv_Greens-'//adjustl(abstr)
                    else
                        abstr='Ret_Greens-'//adjustl(abstr)
                    endif
                    do k=1,30
                        if(abstr(k:k).eq.' ') then
                            abstr(k:k)='_'
                            if(j.lt.10) then
                                write(abstr(k+1:k+1),'(I1)') i
                            elseif(j.lt.100) then
                                write(abstr(k+1:k+2),'(I2)') i
                            elseif(j.lt.1000) then
                                write(abstr(k+1:k+3),'(I3)') i
                            endif
                            exit
                        endif
                    enddo
                    abstr=adjustl(abstr)

                    call MEM(abstr)
                    inquire(file='expdata',exist=toverride)
                    if(toverride) exit lp
                enddo
            enddo
            if(tSpecificMatEls.and.(jval.gt.Contributing_Js)) exit  !We have run through all allowed j's
        enddo lp
        if(tSpecificMatEls) deallocate(SpecificMatEls,Allowed_j)
        return
    endif

    write(6,*) 
    write(6,"(A)") "==========================================="
    write(6,"(A)") " ENTERING GREENS FUNCTION CALCULATION ... "
    write(6,"(A)") "==========================================="
    write(6,*) 

    Greens_Time%timer_name='GreensFunc'
    call set_timer(Greens_Time,30)

    if(tHPHF) then
        call stop_all(this_routine,"Cannot calculate GF with HPHF - bug ghb if needed")
    endif
    if(lenof_sign.ne.1) then
        call stop_all(this_routine,"Cannot yet calculate greens functions with complex walkers")
    endif

    tRunningGF=.true.   !Indicate that we are in the GF routines
    tSearchTau=.false.  !Cannot have tau changing
    tIncrementPops=.false.  !Walkers always written out simply to "POPSFILE"
    tReadPops=.true.    !Ensure that each time we reinitialise we read walkers from disk

    !Initialization
    tSkipDet=.false.
    tAllCommunicate=.false.

    !Allocate Gij(tau) : [i,j,time]
    if(.not.tJustAdvGF) then
        allocate(G_ret(nBasis,nBasis,0:nTimePnts_ret),stat=ierr)
        allocate(G_ret_loc(nBasis,nBasis,0:nTimePnts_ret),stat=ierr)
        allocate(G_ret_sq_loc(nBasis,nBasis,0:nTimePnts_ret),stat=ierr)
        allocate(G_ret_sq(nBasis,nBasis,0:nTimePnts_ret),stat=ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Alloc error")
        G_ret(:,:,:)=0.0_dp
        G_ret_sq(:,:,:)=0.0_dp
    endif
    if(.not.tJustRetGF) then
        allocate(G_adv(nBasis,nBasis,0:nTimePnts_adv),stat=ierr)
        allocate(G_adv_loc(nBasis,nBasis,0:nTimePnts_adv),stat=ierr)
        allocate(G_adv_sq_loc(nBasis,nBasis,0:nTimePnts_adv),stat=ierr)
        allocate(G_adv_sq(nBasis,nBasis,0:nTimePnts_adv),stat=ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Alloc error")
        G_adv(:,:,:)=0.0_dp
        G_adv_sq(:,:,:)=0.0_dp
    endif

    !Write out popsfile
    !Ensure the initial walker distribution is saved to disk
    !Idealy, this would actually be written out to a file with a different name, and then
    !overwritten, and deleted at the end.
    if(.not.tPopsfile) then
        call WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
    endif

    GreenTimeslice(1)=nint((GreensImTime(1)/real(nTimePnts_adv,dp))/tau,sizeof_int)
    GreensIters(1)=GreenTimeslice(1)*nTimePnts_adv
    GreenTimeslice(2)=nint((GreensImTime(2)/real(nTimePnts_ret,dp))/tau,sizeof_int)
    GreensIters(2)=GreenTimeslice(2)*nTimePnts_ret

    IFDEBUGTHEN(GFDebug,2)
        if(.not.tJustRetGF) then
            write(6,*) "Advanced Imaginary time: ",GreensImTime(1)
        endif
        if(.not.tJustAdvGF) then
            write(6,*) "Retarded Imaginary time: ",GreensImTime(2)
        endif
    ENDIFDEBUG

    if(.not.tJustRetGF) then
        write(6,*) "Number of iterations for sweep through advanced imaginary time: ",GreensIters(1)
    endif
    if(.not.tJustAdvGF) then
        write(6,*) "Number of iterations for sweep through retarded imaginary time: ",GreensIters(2)
    endif

    if(.not.tSpecificMatEls) then
        if(tNoDiag) then
            write(6,"(A)") "No diagonal matrix elements of the greens functions will be calculated"
        else
            write(6,"(A)") "All matrix elements of green function will be calculated"
        endif
    else
        !Work out which j orbitals will actually be run over if tSpecificMatEls is set
        call SetupSpecificMatEls()
    endif

    !Save original N-electron reference determinant and energy
    allocate(OrigRefDet(nel))
    OrigRefDet(:)=ProjEDet(:)
    OrigRefiLut(:)=iLutRef(:)
    OrigRefEnergy=Hii

    !Also, save the initial value of SftDamp, and some other variables which we will want to return to
    SftDampInput=SftDamp
    nel_orig=nel
    nOccAlpha_orig=nOccAlpha
    nOccBeta_orig=nOccBeta

    IFDEBUG(GFDebug,2) write(6,*) "Saved Det: ",OrigRefDet(:),OrigRefEnergy

    tChangeProjEDet=.false. !Do not allow reference to change any more
    tRestartHighPop=.false.

    if(iProcIndex.eq.Root) then
        !Work out how many i,j's we want to use
        allocate(UnitNos(2,nBasis*nBasis),stat=ierr)
        UnitNos(:,:)=-1
        pairs=0
        jval = 1
        do j=1,nBasis
            if(tSpecificMatEls) then
                if(Allowed_j(jval).ne.j) cycle  !No desired i's here
                jval=jval+1
            endif
            do i=1,nBasis

                if(.not.CheckAllowed_i(j,i)) cycle  !Not a desired i transistion
                pairs=pairs+1
                PairIndex=FindPairIndex(j,i)

                do GFType=1,2
                    if(tJustAdvGF.and.(GFType.eq.2)) cycle
                    if(tJustRetGF.and.(GFType.eq.1)) cycle

                    if(GFType.eq.1) then
                        UnitNos(1,PairIndex) = get_free_unit()
                    else
                        UnitNos(2,PairIndex) = get_free_unit()
                    endif
                    abstr=''
                    write(abstr,'(I3)') j
                    if(GFType.eq.1) then
                        abstr='Adv_Greens-'//adjustl(abstr)
                    else
                        abstr='Ret_Greens-'//adjustl(abstr)
                    endif
                    do k=1,30
                        if(abstr(k:k).eq.' ') then
                            abstr(k:k)='_'
                            if(j.lt.10) then
                                write(abstr(k+1:k+1),'(I1)') i
                            elseif(j.lt.100) then
                                write(abstr(k+1:k+2),'(I2)') i
                            elseif(j.lt.1000) then
                                write(abstr(k+1:k+3),'(I3)') i
                            endif
                            exit
                        endif
                    enddo
                    abstr=adjustl(abstr)
                    if(GFType.eq.1) then
                        open(unit=UnitNos(1,PairIndex),file=abstr,status='unknown')
                    else
                        open(unit=UnitNos(2,PairIndex),file=abstr,status='unknown')
                    endif
                enddo
            enddo
            if(tSpecificMatEls.and.(jval.gt.Contributing_Js)) exit  !We have run through all allowed j's
        enddo
        write(6,"(I6,A)") pairs," pairs of orbitals will be used for the calculation of the 1-body GFs."
    endif

    !Deallocate all arrays and knowledge of previous simulation
    IFDEBUG(GFDebug,2) write(6,*) "Deallocating memory"
    call DeallocFCIMCMemPar()

    !run over the number of times that we want to sample the greens function for statistics
    do sample=1,GFSamples
        
        if(.not.tJustRetGF) then
            G_adv_loc(:,:,:)=0.0_dp
            G_adv_sq_loc(:,:,:)=0.0_dp
        endif
        if(.not.tJustAdvGF) then
            G_ret_loc(:,:,:)=0.0_dp
            G_ret_sq_loc(:,:,:)=0.0_dp
        endif

        !Initialise by reading in popsfile (do not delete any orbitals/walkers)
        root_write (6,"(A,I5)") "Reading walkers in from disk to start sample ",sample
        FDet(:)=OrigRefDet(:)   !Return the reference determinant to the original one
        call SetupParameters()
        call InitFCIMCCalcPar()
        IFDEBUG(GFDebug,2) write(6,*) "Finished read from disk"

        if(sample.ne.1) then
            !continue running simulation for time b', so that memory of previous config forgotton

            root_write (6,"(A,I7,A,F9.5,A)") "Running for ",nIterEquilGF," iterations, or ",    &
                    nIterEquilGF*Tau," IT before next sampling"
            SftDamp=SftDampInput    !Allow the shift to vary again

            do EquilIter=1,nIterEquilGF

                call sub_dispatcher_7(PerformFciMCycPar, &
                                      ptr_excit_generator, &
                                      ptr_attempt_create, &
                                      ptr_get_spawn_helement, &
                                      ptr_encode_child, &
                                      ptr_new_child_stats, &
                                      ptr_attempt_die, &
                                      ptr_iter_data)

                if(mod(EquilIter,StepsSft).eq.0) then
                    call calculate_new_shift_wrapper(iter_data_fciqmc,TotParts)
                    if(tRestart) call stop_all(this_routine,"All walkers died during equilibration time")
                    root_write (6,"(A,I6,F12.7,2I12)") "Equilibration iter: ",   &
                                EquilIter,DiagSft,sum(AllTotParts),AllTotWalkers
!                    call neci_flush(6)
                endif
            enddo

            root_write (6,"(A)") "Equilibration time finished"
            !Save initial configuration of walkers by writing out to POPSFILE
            root_write (6,"(A)") "Saving final walker configuration by writing to disk"
            call WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
        endif

        allocate(CurrentDetsSaved(0:NIfTot,1:TotWalkers),stat=ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error allocating")
        CurrentDetsSaved(:,:)=CurrentDets(:,1:TotWalkers)
        if(tHashWalkerList) then
            allocate(HashIndexSaved(0:nClashMax,nWalkerHashes),stat=ierr)
            if(ierr.ne.0) call stop_all(this_routine,"Error allocating")
            HashIndexSaved(:,:)=HashIndex(:,:)
        endif
        TotWalkersSaved=TotWalkers

        !Calculate normalization of the wavefunction
        norm=0
        do k=1,TotWalkers
            call extract_sign(CurrentDets(:,k),SignCurr)
            norm = norm + sum(int(SignCurr,int64)**2)
        enddo
        call MPISumAll(norm,Allnorm)

        !Fix shift at current averaged energy value (apply level shifting?)
        root_write (6,"(A,F20.10)") "Fixing shift value at: ",ProjectionE
        DiagSft=ProjectionE
        SftDamp=0.0_dp  !Stop shift changing
        tSaveHashList=.true.    !Stop the has lists getting deallocated

        jval = 1    !The allowed j we are on

        !Loop over orbital j (rightmost annihilation/creation operator acting on FCIQMC wavefuncion)
        do j=1,nBasis

            if(tSpecificMatEls) then
                if(Allowed_j(jval).ne.j) cycle  !No desired i's here
                jval=jval+1
            endif

            !GFType 1 = ADVANCED GF: simulate N-1 electron system and propagate forwards in IT
            !GFType 2 = RETARDED GF: simulate N+1 electron system and propagate backwards in IT
            do GFType=1,2

                if(tJustAdvGF.and.(GFType.eq.2)) cycle
                if(tJustRetGF.and.(GFType.eq.1)) cycle

                if(GFType.eq.1) then
                    root_write (6,"(A,I4,A)") "Deleting orbital :",     &
                        j," to calculate advanced GF, and running N-1 system forwards."
                else
                    root_write (6,"(A,I4,A)") "Creating orbital :",     &
                        j," to calculate retarded GF, and running N+1 system backwards."
                    tBackwardsTime=.true.   !Imaginary time runs backwards in simulations
                endif

                !Scrub all previous information 
                call DeallocFCIMCMemPar()
                
                !Pick new reference det, deallocate all previous FCIQMC arrays, and restart, with N-/+1 electrons and a reduced number of determinants
                !return to saved initial configuration - read in walkers from POPSFILE
                !As loading in, remove orbital j from all determinants (annihilate all walkers if not occupied, otherwise create N-1 determinant)
                !Huge reduction in number of walkers
                IFDEBUG(GFDebug,2) write(6,*) "Setting up for GF calculation"
                call SetupForGF(j,GFType)

                !Temporary storage for a normal ordered list of the new length
                allocate(nI_GF(nel),stat=ierr)
                if(ierr.ne.0) call stop_all(this_routine,"Alloc err")
                
                root_write(6,"(A,I10,A,I10)") "Setup finished, running for total iterations: ",     &
                    GreensIters(GFType)," and calculating GF every ",GreenTimeslice(GFType)

                TimeNumber=0
                do iter=0,GreensIters(GFType)
                    !Call FCIQMC for a certain amout of imaginary time to loop over time range for GF

                    if(iter.ne.0) then
                
                        call sub_dispatcher_7(PerformFciMCycPar, &
                                          ptr_excit_generator, &
                                          ptr_attempt_create, &
                                          ptr_get_spawn_helement, &
                                          ptr_encode_child, &
                                          ptr_new_child_stats, &
                                          ptr_attempt_die, &
                                          ptr_iter_data)
                    endif

                    !Collate some stats
                    TotWalkersTemp=TotWalkers-HolesInList
                    call MPISumAll((/TotWalkersTemp,TotParts/),int64_tmp(1:1+lenof_sign))
                    AllTotWalkers = int64_tmp(1)
                    AllTotParts = int64_tmp(2:1+lenof_sign)

                    if(mod(iter,StepsSft).eq.0) then
                        if(GFType.eq.1) then
                            root_write (6,"(A,I8,2I12)") "N-1 Adv GF Iter: ",iter,sum(AllTotParts),AllTotWalkers
                        else
                            root_write (6,"(A,I8,2I12)") "N+1 Adv GF Iter: ",iter,sum(AllTotParts),AllTotWalkers
                        endif
                    endif
                    
                    if(mod(iter,GreenTimeslice(GFType)).eq.0) then
                        !Every t timeslices, run over all i 
                        if(iter.ne.0) TimeNumber=TimeNumber+1

                        if(GFType.eq.1) then
                            root_write (6,"(A,I6,A,I6)") "Calculating advanced GF point",       &
                                    TimeNumber,' / ',nTimePnts_adv
                        else
                            root_write (6,"(A,I6,A,I6)") "Calculating retarded GF point",       &
                                    TimeNumber,' / ',nTimePnts_ret
                        endif

                        do i=1,nBasis
                
                            if(.not.CheckAllowed_i(j,i)) cycle  !Not a desired i transistion 
                            if((i.ne.j).and.(GFType.eq.2).and.(.not.tJustRetGF).and.(iter.eq.0)) cycle
                                !We don't need to sample the t=0 bit for both advance and retarded
                                !(Except for diagonal matrix elements, where they'll be different

                            IFDEBUGTHEN(GFDebug,2)
                                if(GFType.eq.1) then
                                    write(6,*) "Creating symmetry-allowed orbital: ",i
                                else
                                    write(6,*) "Annihilating symmetry-allowed orbital: ",i
                                endif
                            ENDIFDEBUG

!We must all loop for the same number of times, so that processors are synchronized and can communicate when necessary
                            call MPIAllReduce(TotWalkers,MPI_MAX,MaxWalkersProc)
!                            write(6,*) "Max Determinant number is: ",MaxWalkersProc

                            ValidSpawnedList = InitialSpawnedSlots
                            do walk=1,MaxWalkersProc
                                !Run through determinant list on this core

                                if(walk.le.TotWalkers) then
                                    IFDEBUGTHEN(GFDebug,6)
                                        call extract_sign(CurrentDets(:,walk),SignCurr)
                                        call decode_bit_det(nI_GF,CurrentDets(:,walk))
                                        write(6,*) walk,TotWalkers,SignCurr,nI_GF(:)
                                    ENDIFDEBUG

                                    if(tHashWalkerList) then
                                        call extract_sign(CurrentDets(:,walk),SignCurr)
                                        if(IsUnoccDet(SignCurr)) then
                                            tSkipDet=.true.
                                        else
                                            tSkipDet=.false.
                                        endif
                                    endif

                                    if(.not.tSkipDet) then
                                        if(GFType.eq.1) then
                                            !if i is occupied in determinant, then cycle (it is destroyed)
                                            !Check if its occupied
                                            if(btest(CurrentDets((i-1)/bits_n_int,walk),mod(i-1,bits_n_int))) then
                                                tSendDet=.false.
                                            else
                                                tSendDet=.true.
                                                !i is not occupied - occupy it
                                                TempDet(:)=CurrentDets(:,walk)
                                                set_orb(TempDet,i)
                                            endif
                                        else
                                            !Check it is not occupied
                                            if(.not.btest(CurrentDets((i-1)/bits_n_int,walk),mod(i-1,bits_n_int))) then
                                                tSendDet=.false.
                                            else
                                                tSendDet=.true.
                                                !i is occupied - distroy it
                                                TempDet(:)=CurrentDets(:,walk)
                                                clr_orb(TempDet,i)
                                            endif
                                        endif
                                        if(tSendDet) then
                                            !TempDet now contains the N-electron determinant
                                            !We now need to find the weight of this determinant in the original list
                                            !However, this will involve communication since it is on proc
                                            call decode_bit_det_nel(nI,TempDet(:),nel_orig)
                                            proc = DetermineDetNode_nel(nI,0,nel_orig)
                                            
                                            SpawnedParts(:,ValidSpawnedList(proc))=TempDet(:)
                                            ValidSpawnedList(proc) = ValidSpawnedList(proc) + 1
                                            IFDEBUGTHEN(GFDebug,6)
                                                write(6,*) "Found allowed weight on determinant: ", nI(:)
                                                write(6,*) "Sending to processor: ",proc
                                                call flush(6)
                                            ENDIFDEBUG

                                            !Do we need to comminicate now?
                                            if(nNodes.gt.1) then
                                                rat = real(ValidSpawnedList(proc)-InitialSpawnedSlots(proc),dp)/&
                                                        real(InitialSpawnedSlots(1),dp)
                                            else
                                                rat = real(ValidSpawnedList(0),dp) / real(MaxSpawned,dp)
                                            endif
                                            if(rat.gt.0.95) then
                                                tCommunicate=.true.
                                            else
                                                tCommunicate=.false.
                                            endif
                                        endif   !Whether the determinant is actually allowed or not
                                    endif   !Whether to skip the det (for non-contiguous tHashWalkerList storage)
                                else
                                    tCommunicate=.false.
                                endif   !Whether we have actually run through all our determinants already
                                call MPIAllReduce(tCommunicate,MPI_LOR,tAllCommunicate)

                                if(tAllCommunicate) then

                                    IFDEBUG(GFDebug,3) write(6,*) "Communicating lists, even though we have not run through all"
                                    !Communicate the spawned lists
                                    call SendDetsToProcessors(MaxIndex)
                                    if(.not.tHashWalkerList) then
                                        !Sort for fasterbinary searching later on
                                        call sort(SpawnedParts(:,1:MaxIndex))
                                    endif
                                    !Now we have the determinants on the correct processors
                                    !There will only be a maximum of one instance of each determinant
                                    MinInd=1
                                    do k=1,MaxIndex
                                        !look up weight of determinant (with i) in original
                                        !list, and multiply it by the weight in the new set.
                                        call FindDetPos(k,MinInd,PartInd,tSuccess)
                                        MinInd=PartInd  !We can search a smaller list next time
                                    
                                        !We now want to sum in Ci(t)*Cj(t')
                                        if(tSuccess) then
                                            !Find the two weights
                                            call extract_sign(SpawnedParts(:,k),SignRight)
                                            call extract_sign(CurrentDetsSaved(:,PartInd),SignLeft)

                                            Ex(1,1)=1
                                            call GetBitExcitation(SpawnedParts(0:NIfD,k),CurrentDetsSaved(0:NIfD,PartInd),Ex,tSign)
                                            if(tSign) SignRight=-SignRight  !Flip the sign of one of them

                                            !Sum into this timeslice in the greensfunction
                                            if(GFType.eq.1) then
                                                !Advanced
                                                G_adv_loc(i,j,TimeNumber)=G_adv_loc(i,j,TimeNumber)+    &
                                                    real(SignRight(1)*SignLeft(1),dp)/AllNorm
                                            else
                                                G_ret_loc(i,j,TimeNumber)=G_ret_loc(i,j,TimeNumber)+    &
                                                    real(SignRight(1)*SignLeft(1),dp)/AllNorm
                                            endif
                                        endif

                                    enddo
                                    ValidSpawnedList = InitialSpawnedSlots
                                    tAllCommunicate=.false.
                                endif

                            enddo   !Enddo sum over walkers 
                            IFDEBUG(GFDebug,2) write(6,*) "Finished running over all walkers to check can make allowed determinants"

                            !Final communication and summing in
                            IFDEBUGTHEN(GFDebug,2) 
                                write(6,*) "Sending final walkers to correct processors"
                                IFDEBUG(GFDebug,4) write(6,*) "VSL: ",ValidSpawnedList(:)
                            ENDIFDEBUG
                            !Communicate the spawned lists
                            call SendDetsToProcessors(MaxIndex)
                            IFDEBUG(GFDebug,4) write(6,*) "MaxIndex: ",MaxIndex
                            if(.not.tHashWalkerList) then
                                !Sort for fasterbinary searching later on
                                call sort(SpawnedParts(:,1:MaxIndex))
                            endif
!                            do k=1,MaxIndex
!                                write(6,*) SpawnedParts(:,k)
!                            enddo
                            !Now we have the determinants on the correct processors
                            !** WORK OUT SIGN ISSUES WHICH WILL BE PRESENT HERE **
                            !There will only be a maximum of one instance of each determinant
                            MinInd=1
                            do k=1,MaxIndex
                                call FindDetPos(k,MinInd,PartInd,tSuccess)
                                MinInd=PartInd  !We can search a smaller list next time
                            
                                !We now want to sum in Ci(t)*Cj(t')
                                if(tSuccess) then
                                    !Find the two weights
                                    call extract_sign(SpawnedParts(:,k),SignRight)
                                    call extract_sign(CurrentDetsSaved(:,PartInd),SignLeft)

                                    Ex(1,1)=1
                                    call GetBitExcitation(SpawnedParts(0:NIfD,k),CurrentDetsSaved(0:NIfD,PartInd),Ex,tSign)
                                    if(tSign) SignRight=-SignRight  !Flip the sign of one of them
                                    
                                    IFDEBUG(GFDebug,5) write(6,*) "Found matching det: ",SignRight,SignLeft

                                    !Normalize the contribution to the greens function from this sample
                                    !by the normalization of the original wavefunction
                                    if(GFType.eq.1) then
                                        !Advanced
                                        G_adv_loc(i,j,TimeNumber)=G_adv_loc(i,j,TimeNumber)+    &
                                            real(SignRight(1)*SignLeft(1),dp)/AllNorm
                                    else
                                        G_ret_loc(i,j,TimeNumber)=G_ret_loc(i,j,TimeNumber)+    &
                                            real(SignRight(1)*SignLeft(1),dp)/AllNorm
                                    endif
                                endif

                            enddo
                            ValidSpawnedList = InitialSpawnedSlots

                        enddo   !Enddo loop over i

                    endif   !Endif time datapoint

                enddo   !Enddo GF time range

                !Return electron number to full N system
                call ReturnElectronNumber(GFType)
                    
                if(GFType.eq.1) then
                    if(TimeNumber.ne.nTimePnts_adv) call stop_all(this_routine,"Incorrect time points")
                else
                    if(TimeNumber.ne.nTimePnts_ret) call stop_all(this_routine,"Incorrect time points")
                endif
                deallocate(nI_GF)

            enddo   !End retarded/advanced
            tBackwardsTime=.false.
            
            if(tSpecificMatEls.and.(jval.gt.Contributing_Js)) exit  !We have run through all allowed j's

        enddo   !Enddo orbital j

        if(.not.tJustRetGF) then
 !For Squared Values (cannot square each contribution as not linear)
            call MPISumAll(G_adv_loc,G_adv_sq_loc) 
        endif
        if(.not.tJustAdvGF) then
            call MPISumAll(G_ret_loc,G_ret_sq_loc)  
        endif
        
        do i=1,nBasis
            do j=1,nBasis
                if(.not.tJustRetGF) then
                    do k=0,nTimePnts_adv
                        G_adv(i,j,k)=G_adv(i,j,k)+G_adv_loc(i,j,k)  !Sum into averaged value
                        G_adv_sq(i,j,k)=G_adv_sq(i,j,k)+G_adv_sq_loc(i,j,k)**2.0_dp  !Square contribution
                    enddo
                endif
                if(.not.tJustAdvGF) then
                    do k=0,nTimePnts_ret
                        G_ret(i,j,k)=G_ret(i,j,k)+G_ret_loc(i,j,k)  !Sum into averaged value
                        G_ret_sq(i,j,k)=G_ret_sq(i,j,k)+G_ret_sq_loc(i,j,k)**2.0_dp  !Square contribution
                    enddo
                endif
            enddo
        enddo

        tSaveHashList=.false.   !Regenerate the lists to avoid correlation?
        call DeallocFCIMCMemPar()
        deallocate(CurrentDetsSaved)
        if(tHashWalkerList) deallocate(HashIndexSaved)

    enddo   !Enddo samples

    !Sum in parallel
    if(.not.tJustRetGF) then
        deallocate(G_adv_loc,G_adv_sq_loc)
        allocate(G_adv_all(nbasis,nbasis,0:nTimePnts_adv),stat=ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Alloc err")
        G_adv_all(:,:,:)=0.0_dp
        call MPISumAll(G_adv,G_adv_all)
    endif
    if(.not.tJustAdvGF) then
        deallocate(G_ret_loc,G_ret_sq_loc)
        allocate(G_ret_all(nbasis,nbasis,0:nTimePnts_ret),stat=ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Alloc err")
        G_ret_all(:,:,:)=0.0_dp
        call MPISumAll(G_ret,G_ret_all)
    endif
    
    !Find average now by dividing by GFSamples
    do i=1,nBasis
        do j=1,nBasis
            if(.not.tJustRetGF) then
                do k=0,nTimePnts_adv
                    G_adv_all(i,j,k)=G_adv_all(i,j,k)/real(GFSamples,dp)
                    G_adv_sq(i,j,k)=G_adv_sq(i,j,k)/real(GFSamples,dp)
                enddo
            endif
            if(.not.tJustAdvGF) then
                do k=0,nTimePnts_ret
                    G_ret_all(i,j,k)=G_ret_all(i,j,k)/real(GFSamples,dp)
                    G_ret_sq(i,j,k)=G_ret_sq(i,j,k)/real(GFSamples,dp)
                enddo
            endif
        enddo
    enddo

    !Symmetrize greens function?
    !(Work out symmetrization error)

    !Write out greens function
    if(iProcIndex.eq.Root) then
        jval=1
        do j=1,nBasis
            if(tSpecificMatEls) then
                if(Allowed_j(jval).ne.j) cycle  !No desired i's here
                jval=jval+1
            endif
            do i=1,nBasis

                if(.not.CheckAllowed_i(j,i)) cycle   !Not an allowed i orbital with this j
        
                PairIndex=FindPairIndex(j,i)

                if(.not.tJustRetGF) then
                    if(UnitNos(1,PairIndex).eq.-1) call stop_all(this_routine,"Unit not open for writing")
                    write(UnitNos(1,PairIndex),"(A)") "# -Im.Time  Av.imaginary_GF  Av.Im_GF_sq    Error    "
!                    do k=nTimePnts_adv,0,-1
                    do k=0,nTimePnts_adv
                        write(UnitNos(1,PairIndex),"(4G25.15)") k*tau,-G_adv_all(i,j,k),G_adv_sq(i,j,k),   &
                            sqrt((G_adv_sq(i,j,k)-(G_adv_all(i,j,k)**2.0_dp))/real(GFSamples,dp))
                    enddo
                    close(UnitNos(1,PairIndex))
                endif
                if(.not.tJustAdvGF) then
                    if((.not.tJustRetGF).and.(i.ne.j)) then
                        !We don't want to write out the t=0 bit twice, apart from diagonal
                        !elements where they'll be different
                        startInd=1
                    else
                        startInd=0
                    endif
                    if(UnitNos(2,PairIndex).eq.-1) call stop_all(this_routine,"Unit not open for writing")
                    write(UnitNos(2,PairIndex),"(A)") "# Im.Time  Av.imaginary_GF  Av.Im_GF_sq    Error    "
                    do k=startInd,nTimePnts_ret
                        write(UnitNos(2,PairIndex),"(4G25.15)") k*tau,-G_ret_all(i,j,k),G_ret_sq(i,j,k),    &
                            sqrt((G_ret_sq(i,j,k)-(G_ret_all(i,j,k)**2.0_dp))/real(GFSamples,dp))
                    enddo
                    close(UnitNos(2,PairIndex))
                endif
            enddo
            if(tSpecificMatEls.and.(jval.gt.Contributing_Js)) exit  !We have run through all allowed j's
        enddo

        deallocate(UnitNos)
    endif
        
    if(.not.tJustAdvGF) then
        deallocate(G_ret,G_ret_all,G_ret_sq)
    endif
    if(.not.tJustRetGF) then
        deallocate(G_adv,G_adv_all,G_adv_sq)
    endif
    deallocate(OrigRefDet)
    
    call halt_timer(Greens_Time)

    if(tCalcSpectrum) then
        !Perform maximum entropy analysis to infer the spectrum from the greens functions
        jval=1
        do j=1,nBasis
            if(tSpecificMatEls) then
                if(Allowed_j(jval).ne.j) cycle  !No desired i's here
                jval=jval+1
            endif
            do i=1,nBasis

                if(.not.CheckAllowed_i(j,i)) cycle   !Not an allowed i orbital with this j
                do GFType=1,2
                    if(tJustAdvGF.and.(GFType.eq.2)) cycle
                    if(tJustRetGF.and.(GFType.eq.1)) cycle

                    abstr=''
                    write(abstr,'(I3)') j
                    if(GFType.eq.1) then
                        abstr='Adv_Greens-'//adjustl(abstr)
                    else
                        abstr='Ret_Greens-'//adjustl(abstr)
                    endif
                    do k=1,30
                        if(abstr(k:k).eq.' ') then
                            abstr(k:k)='_'
                            if(j.lt.10) then
                                write(abstr(k+1:k+1),'(I1)') i
                            elseif(j.lt.100) then
                                write(abstr(k+1:k+2),'(I2)') i
                            elseif(j.lt.1000) then
                                write(abstr(k+1:k+3),'(I3)') i
                            endif
                            exit
                        endif
                    enddo
                    abstr=adjustl(abstr)

                    call MEM(abstr)
                enddo
            enddo
            if(tSpecificMatEls.and.(jval.gt.Contributing_Js)) exit  !We have run through all allowed j's
        enddo
    endif

    if(tSpecificMatEls) deallocate(SpecificMatEls,Allowed_j)

    end subroutine CalcGreensFunc

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

        tNoDiag=.true.  !This will be turned off if any of the transistions are diagonal

        do j=1,nBasis
            jAlreadyCalc=.false.
            do i=1,nBasis
                if(tNoDiag.and.(i.eq.j)) cycle
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
        use FciMCData, only: freeslot,WalkVecH,fcimc_excit_gen_store
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
        use AnnihilationMod, only: FindWalkerHash,EnlargeHashTable,IsUnoccDet
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
            tempnel=tempnel-1
        else
            tempnel=tempnel+1
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

                    !Remove/add the orbital
                    !GFType 1 = ADVANCED GF: simulate N-1 electron system and propagate forwards in IT
                    !GFType 2 = RETARDED GF: simulate N+1 electron system and propagate backwards in IT
                    if(GFType.eq.1) then
                        !Remove Orb if present, and keep walkers
                        !If not present, then remove all walkers
                        if(IsOcc(WalkerTemp,Orb)) then
                            clr_orb(WalkerTemp,Orb)
                            do i=1,nel  !(original nel)
                                newelec=1
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
        CurrParts=TotParts
        call MPISum(TotParts,1,AllTotParts)

        !Find the *global* lowest energy determinant
        LowDet(1)=LowestEnergy
        LowDet(2)=real(iProcIndex,dp)
        call MPIAllReduceDatatype(LowDet,1,MPI_MINLOC,MPI_2DOUBLE_PRECISION,LowDetOut)
        call MPIBCast(LowDetOut(1),1,nint(LowDetOut(2)))

        if(iProcIndex.eq.nint(LowDetOut(2))) then
            !This has the reference det
            ilutHF(:)=Dets(:,LowestIndex)
            call extract_sign(Dets(:,LowestIndex),SignTemp)
            CurrHF=CurrHF+SignTemp
        endif
        NoatHF=CurrHF
        call MPISum(NoatHF,1,AllNoatHF)

        call MPIBCast(ilutHF,NIfTot+1,nint(LowDetOut(2)))

        !Now, we know what we want our reference to be.
        do i=1,TotWalkers
            CurrentH(i)=CurrentH(i)-LowestEnergy
        enddo
        
        ilutRef(:)=ilutHF(:)
        call decode_bit_det (HFDet, ilutHF)
        ProjEDet(:)=HFDet(:)
        write(6,*) "Lowest energy reference found to be: ",HFDet(:)

        call mpibarrier(err)

    end subroutine ReadInForGF

    subroutine SetGFDefaults()
        implicit none

        tGreensFuncs=.false.
        tJustAdvGF=.false.
        tJustRetGF=.false.
        nIterEquilGF=1000
        GFSamples=100
        GreensImTime(:)=1.0_dp
        nTimePnts_ret=10
        nTimePnts_adv=10
        tNoDiag=.false.
        iSpecificMatEls=0
        tSpecificMatEls=.false.
        tJustSpectrum=.false.
        tCalcSpectrum=.true.

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
            case("JUSTADVANCED")
                tJustAdvGF=.true.
            case("JUSTRETARDED")
                tJustRetGF=.true.
            case("ADVIMTIME")
                call readf(GreensImTime(1))
            case("RETIMTIME")
                call readf(GreensImTime(2))
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
            case("END")
                exit
            case default
                write(6,"(A)") "ALLOWED KEYWORDS IN GREENSFUNCS BLOCK: "
                write(6,"(A)") "SAMPLES"
                write(6,"(A)") "JUSTADVANCED"
                write(6,"(A)") "JUSTRETARDED"
                write(6,"(A)") "ADVIMTIME"
                write(6,"(A)") "RETIMTIME"
                write(6,"(A)") "ADVTIMEPOINTS"
                write(6,"(A)") "RETTIMEPOINTS"
                write(6,"(A)") "EQUILIBRATION"
                write(6,"(A)") "NODIAGGF"
                write(6,"(A)") "TRANSISTION"
                write(6,"(A)") "JUSTSPECTRUM"
                write(6,"(A)") "NOSPECTRUM"
                call report("Keyword "//trim(w)//" not recognized",.true.)
            end select
        enddo GF

    end subroutine GFReadInput


end module GreensFuncMod

!Wrap greens function routines since we want to use the FCIMC module routines in it
subroutine CalcGreensFunc_wrap
    use GreensFuncMod, only: CalcGreensFunc
    implicit none

    call CalcGreensFunc()

end subroutine CalcGreensFunc_wrap

