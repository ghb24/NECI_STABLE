#include "macros.h"

module GreensFuncFreqMod
    use SystemData, only: nel, nBasis,tHPHF,tFixLz
    use SystemData, only: tKPntSym,nOccAlpha,nOccBeta
    use global_utilities
    use constants, only: dp, int64, n_int
    use FciMCParMod, only: DeallocFCIMCMemPar,SetupParameters,InitFCIMCCalcPar
    use FciMCParMod, only: PerformFciMCycPar,calculate_new_shift_wrapper
    use FciMCData, only: tRunningGF,tBackwardsTime,tSearchTau,CurrentDets,TotWalkers
    use FciMCData, only: ProjEDet,iLutRef,Hii,iter_data_fciqmc,TotParts
    use FciMCData, only: tRestart,AllTotParts,AllTotWalkers,tHashWalkerlist,HashIndex
    use FciMCData, only: ProjectionE,tSaveHashList,SpawnedParts
    use FciMCData, only: ValidSpawnedList,InitialSpawnedSlots
    use FciMCData, only: nClashMax,tMultipleH,DriveFactor
    use FciMCData, only: nWalkerHashes,MaxSpawned,HolesInList
    use FciMCData, only: ptr_excit_generator,ptr_attempt_create,ptr_get_spawn_helement
    use FciMCData, only: ptr_new_child_stats,ptr_encode_child,ptr_attempt_die,ptr_iter_data
    use FciMCData, only: tDrivePropagation, SftDampDrive,HashIndexArr1,HashIndexArr2,FreeSlot 
    use CalcData, only: tReadPops,tChangeProjEDet,tRestartHighPop,tRegenDiagHEls
    use CalcData, only: tau,DiagSft,SftDamp,StepsSft, tImagEnergyShift,tCheckHighestPop
    use bit_reps, only: NIfD,NIfDBO,NIfTot,extract_sign,decode_bit_det,decode_bit_det_nel
    use FciMCData, only: iMaxH_App
    use Logging, only: tPopsfile,tIncrementPops,GFDebug
    use AnnihilationMod, only: IsUnoccDet,DetermineDetNode_nel
    use DeterminantData, only: FDet
    use PopsfileMod, only: WriteToPopsfileParOneArr
    use util_mod, only: get_unique_filename,get_free_unit
    use Parallel_neci
    use sort_mod
    use GreensFuncData
    use GreensFuncUtils

    implicit none

    contains

    !Call this after appropriate equilibration time
    subroutine CalcGreensFuncFreq()
    implicit none
    integer :: ierr,Ex(2,1),sample,EquilIter,i,j,k,GFType,MinInd
    integer(int64) :: norm, allnorm
    integer :: MaxIndex,PartInd,neltemp,StartInd,nFreq,nFreqSteps
    real(dp) :: OrigRefEnergy,rat,SftDampInput,TauSq,SumBetaSq
    real(dp), dimension(lenof_sign) :: OverlapNew,AllOverlapNew
    real(dp) :: SavedDiagSft,Freq
    integer, allocatable :: OrigRefDet(:),nI_GF(:)
    integer(n_int) :: OrigRefiLut(0:NIfTot),TempDet(0:NIfTot)
    integer :: TimeNumber,iter,walk,proc,nI(nel),PairIndex,Pairs,jval
    integer, dimension(lenof_sign) :: SignCurr,SignRight,SignLeft
    logical :: tCommunicate,tAllCommunicate,tSuccess,tSign,tSkipDet,tSendDet
    logical :: toverride,exists
    integer(int64) :: int64_tmp(1+lenof_sign),TotWalkersTemp,MaxWalkersProc
    real(dp) , allocatable :: G_ret_loc(:,:,:,:),G_adv_loc(:,:,:,:)
    real(dp), allocatable :: OldOverlap(:,:,:)
    character(len=35) :: abstr
    character(len=*), parameter :: t_r='CalcGreensFuncFreq'
    integer, allocatable :: UnitNos(:,:)
    integer :: tempunit,ios,allsamples
    real(dp) :: G_sum_re,G_sum_im,G_sum_sq_re,G_sum_sq_im

    write(6,*) 
    write(6,"(A)") "============================================================"
    write(6,"(A)") " ENTERING FREQUENCY-DOMAIN GREENS FUNCTION CALCULATION ... "
    write(6,"(A)") "============================================================"
    write(6,*) 

    Greens_Time%timer_name='GreensFuncFreq'
    call set_timer(Greens_Time,30)

    !Both these restrictions NEED to go
    if(tHPHF) then
        call stop_all(t_r,"Cannot calculate GF with HPHF - bug ghb if needed")
    endif
    if(tMultipleH) then
        call stop_all(t_r,"Should not be using gaussian propagator yet - need GS distribution")
    endif
    if(tImagEnergyShift) then
        call stop_all(t_r,"Should not have imaginary component in diagonals to start")
    endif
    if(.not.tJustAdvGF) then
        write(6,*) "Currently, can only calculate advanced Greens functions. Hard coding this."
        tJustAdvGF = .true.
    endif
    if(SftDampDrive.ne.0.0_dp) then
        write(6,*) "Turning off changing of drive factor..."
        SftDampDrive = 0.0_dp
    endif

    tRunningGF=.true.   !Indicate that we are in the GF routines
    tSearchTau=.false.  !Cannot have tau changing
    tIncrementPops=.false.  !Walkers always written out simply to "POPSFILE"
    tReadPops=.true.    !Ensure that each time we reinitialise we read walkers from disk
    TauSq=Tau*Tau

    !Initialization
    tSkipDet=.false.
    tAllCommunicate=.false.

    allsamples=0
                
    !Allocate Gij(frequency) : [i,j,Frequency]. Now also potentially complex implementation.
    if(.not.tJustAdvGF) then
        nFreqPnts_ret = nint((FinalFreq_ret-StartFreq_ret)/dFreqStep_ret)
        root_write(6,"(A,I7,A)") "Looping over ",nFreqPnts_ret," different frequencies for retarded GFs"
        allocate(G_ret(lenof_sign,nBasis,nBasis,0:nFreqPnts_ret),stat=ierr)
        allocate(G_ret_loc(lenof_sign,nBasis,nBasis,0:nFreqPnts_ret),stat=ierr)
!        allocate(G_ret_sq_loc(nBasis,nBasis,0:nFreqPnts_ret),stat=ierr)
        allocate(G_ret_sq(lenof_sign,nBasis,nBasis,0:nFreqPnts_ret),stat=ierr)
        allocate(G_ret_all(lenof_sign,nbasis,nbasis,0:nFreqPnts_ret),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc error")
        G_ret_loc(:,:,:,:) = 0.0_dp
        G_ret(:,:,:,:) = 0.0_dp
        G_ret_sq(:,:,:,:) = 0.0_dp
        if(iProcIndex.eq.Root) then
            inquire(file='GFRestart_ret',exist=exists)
            if(exists) then
                !Read in from previous calculation
                write(6,"(A)") "Reading in retarded GF from previous calculation"
                tempunit=get_free_unit()
                open(unit=tempunit,file='GFRestart_ret',status='old')
                read(tempunit,*) allsamples
                do while(.true.)
                    if(lenof_sign.eq.1) then
                        read(tempunit,"(3I6,2G25.15)",iostat=ios) i,j,k,G_sum_re,G_sum_sq_re
                    else
                        read(tempunit,"(3I6,4G25.15)",iostat=ios) i,j,k,G_sum_re,G_sum_im,G_sum_sq_re,G_sum_sq_im
                    endif
                    if(ios.gt.0) call stop_all(t_r,'Error reading restart adv')
                    if(ios.lt.0) exit   !eof
                    G_ret(1,i,j,k)=G_sum_re
                    G_ret_sq(1,i,j,k)=G_sum_sq_re
                    if(lenof_sign.eq.2) then
                        G_ret(lenof_sign,i,j,k)=G_sum_im
                        G_ret_sq(lenof_sign,i,j,k)=G_sum_sq_im
                    endif
                enddo
            endif
            close(tempunit)
        endif
        call MPIBCast(G_ret,lenof_sign*nBasis*nBasis*(nFreqPnts_ret+1))
        call MPIBCast(allsamples)
    endif
    if(.not.tJustRetGF) then
        nFreqPnts_adv = nint((FinalFreq_adv-StartFreq_adv)/dFreqStep_adv)
        root_write(6,"(A,I7,A)") "Looping over ",nFreqPnts_adv," different frequencies for advanced GFs"
        allocate(G_adv(lenof_sign,nBasis,nBasis,0:nFreqPnts_adv),stat=ierr)
        allocate(G_adv_loc(lenof_sign,nBasis,nBasis,0:nFreqPnts_adv),stat=ierr)
!        allocate(G_adv_sq_loc(nBasis,nBasis,0:nFreqPnts_adv),stat=ierr)
        allocate(G_adv_sq(lenof_sign,nBasis,nBasis,0:nFreqPnts_adv),stat=ierr)
        allocate(G_adv_all(lenof_sign,nbasis,nbasis,0:nFreqPnts_adv),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc error")
        G_adv(:,:,:,:) = 0.0_dp
        G_adv_sq(:,:,:,:) = 0.0_dp
        G_adv_loc(:,:,:,:) = 0.0_dp
        if(iProcIndex.eq.Root) then
            inquire(file='GFRestart_adv',exist=exists)
            if(exists) then
                !Read in from previous calculation
                write(6,"(A)") "Reading in advanced GF from previous calculation"
                tempunit=get_free_unit()
                open(unit=tempunit,file='GFRestart_adv',status='old')
                read(tempunit,*) allsamples
                do while(.true.)
                    if(lenof_sign.eq.1) then
                        read(tempunit,"(3I6,2G25.15)",iostat=ios) i,j,k,G_sum_re,G_sum_sq_re
                    else
                        read(tempunit,"(3I6,4G25.15)",iostat=ios) i,j,k,G_sum_re,G_sum_im,G_sum_sq_re,G_sum_sq_im
                    endif
                    if(ios.gt.0) call stop_all(t_r,'Error reading restart adv')
                    if(ios.lt.0) exit   !eof
                    G_adv(1,i,j,k)=G_sum_re
                    G_adv_sq(1,i,j,k)=G_sum_sq_re
                    if(lenof_sign.eq.2) then
                        G_adv(lenof_sign,i,j,k)=G_sum_im
                        G_adv_sq(lenof_sign,i,j,k)=G_sum_sq_im
                    endif
                enddo
            endif
            close(tempunit)
        endif
        call MPIBCast(G_adv,lenof_sign*nBasis*nBasis*(nFreqPnts_adv+1))
        call MPIBCast(allsamples)
    endif
    allocate(OldOverlap(lenof_sign,nBasis,nBasis),stat=ierr)
    if(ierr.ne.0) call stop_all(t_r,"Alloc error")
    OldOverlap(:,:,:) = 0.0_dp

    !This indicates how long to run for in the integration, and the fineness of the integration points
    if(.not.tJustRetGF) then
        if(nTimePnts_adv.eq.-1) then
            if(GreenTimeslice(1).eq.-1) then
                call stop_all(t_r,"Neither ADVTIMEPOINTS nor ADVINTEGRATIONSTEP specified. Specify one.")
            endif
            nTimePnts_adv = nint(real(GreensIters(1),dp)/real(GreenTimeslice(1),dp))
        else
            GreenTimeslice(1)=nint(real(GreensIters(1),dp)/real(nTimePnts_adv,dp))
        endif
    endif
    if(.not.tJustAdvGF) then
        if(nTimePnts_ret.eq.-1) then
            if(GreenTimeslice(2).eq.-1) then
                call stop_all(t_r,"Neither RETTIMEPOINTS nor RETINTEGRATIONSTEP specified. Specify one.")
            endif
            nTimePnts_adv = nint(real(GreensIters(2),dp)/real(GreenTimeslice(2),dp))
        else
            GreenTimeslice(2)=nint(real(GreensIters(2),dp)/real(nTimePnts_ret,dp))
        endif
    endif
    
    !Write out popsfile
    !Ensure the initial walker distribution is saved to disk
    !Idealy, this would actually be written out to a file with a different name, and then
    !overwritten, and deleted at the end.
    if(.not.tPopsfile) then
        call WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
    endif

    !Restrict which GF elements will actually be calculated
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
    SavedDiagSft=DiagSft
    nel_orig=nel
    nOccAlpha_orig=nOccAlpha
    nOccBeta_orig=nOccBeta

    IFDEBUG(GFDebug,2) write(6,*) "Saved Det: ",OrigRefDet(:),OrigRefEnergy

    tChangeProjEDet=.false. !Do not allow reference to change any more
    tRestartHighPop=.false.
    tCheckHighestPop=.false.

    if(iProcIndex.eq.Root) then
        !Work out how many i,j's we want to use
        allocate(UnitNos(2,nBasis*nBasis),stat=ierr)
        UnitNos(:,:)=-1
        pairs=0
        jval = 1
        do j=1,nBasis   !This cycles through all j's, but generally only picks on the Allowed_j ones, otherwise it will cycle
            if(tSpecificMatEls) then
                if(jval.gt.iSpecificMatEls) cycle   !This will now just skip to the end
                if(Allowed_j(jval).ne.j) cycle  !No desired i's here
                jval=jval+1 !jval denotes the 'next' allowed index
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
                        abstr='Adv_FreqGreens-'//adjustl(abstr)
                    else
                        abstr='Ret_FreqGreens-'//adjustl(abstr)
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
        allsamples=allsamples+1 !Add another sample
        
        !Ensure standard exponential propagation dynamics here
        iMaxH_App = 1   !Exponential propagator
        tImagEnergyShift = .false.
        tDrivePropagation = .false.
        tMultipleH=.false.
        
        if(.not.tJustRetGF) then
            G_adv_loc(:,:,:,:)=0.0_dp !Save all the different frequencies, to average between different runs later
        endif
        if(.not.tJustAdvGF) then
            G_ret_loc(:,:,:,:)=0.0_dp
        endif

        !Initialise by reading in popsfile (do not delete any orbitals/walkers)
        root_write (6,"(A,I5)") "Reading walkers in from disk to start sample ",sample
        FDet(1:NEl)=OrigRefDet(1:NEl)   !Return the reference determinant to the original one

        call SetupParameters()
        call InitFCIMCCalcPar()
        IFDEBUG(GFDebug,2) write(6,*) "Finished read from disk"

        if(sample.ne.1) then
            !continue running simulation for time b', so that memory of previous config forgotton

            root_write (6,"(A,I7,A)") "Running for ",nIterEquilGF," iterations before next sampling"
            SftDamp=SftDampInput    !Allow the shift to vary again
            DiagSft=SavedDiagSft

            root_write (6,"(A,I6,F12.7,2I12)") "Equilibration iter: ",   &
                        0,DiagSft,sum(AllTotParts),AllTotWalkers

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
                    if(tRestart) call stop_all(t_r,"All walkers died during equilibration time")
                    root_write (6,"(A,I6,F12.7,2I12)") "Equilibration iter: ",   &
                                EquilIter,DiagSft,sum(AllTotParts),AllTotWalkers
!                    call neci_flush(6)
                endif
            enddo

            root_write (6,"(A)") "Equilibration time finished"
            !Save initial configuration of walkers by writing out to POPSFILE
            root_write (6,"(A)") "Saving final walker configuration by writing to disk"
            call WriteToPopsfileParOneArr(CurrentDets,TotWalkers)
            SavedDiagSft=DiagSft
        endif

        !CurrentDetsSaved saves the ground state
        allocate(CurrentDetsSaved(0:NIfTot,1:TotWalkers),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Error allocating")
        CurrentDetsSaved(:,:)=CurrentDets(:,1:TotWalkers)
        if(tHashWalkerList) then
            allocate(HashIndexSaved(0:nClashMax,nWalkerHashes),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Error allocating")
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
        SftDamp=0.0_dp  !Stop shift changing
        tSaveHashList=.true.    !Stop the has lists getting deallocated

        jval = 1    !The allowed j we are on
        !Loop over orbital j (rightmost annihilation/creation operator acting on FCIQMC wavefuncion)
        do j=1,nBasis

            if(tSpecificMatEls) then
                if(jval.gt.iSpecificMatEls) cycle   !This will now just skip to the end
                if(Allowed_j(jval).ne.j) cycle   !No desired i's here
                jval=jval+1
            endif

            !GFType 1 = ADVANCED GF: simulate N-1 electron system and propagate forwards in IT
            !GFType 2 = RETARDED GF: simulate N+1 electron system and propagate backwards in IT
            do GFType=1,2

                if(tJustAdvGF.and.(GFType.eq.2)) cycle
                if(tJustRetGF.and.(GFType.eq.1)) cycle

                if(GFType.eq.1) then
                    root_write (6,"(A,I4,A)") "Deleting orbital :",     &
                        j," to calculate advanced GF, and integrating N-1 system forwards."
                else
                    root_write (6,"(A,I4,A)") "Creating orbital :",     &
                        j," to calculate retarded GF, and integrating N+1 system backwards."
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
                call SetupForGaussProp()    !Turn on driving, imagdiagshift & MultipleH (as well as allocating space for it)
                call SaveStartingVec()      !Save the a_j|0> vector, so we can return to in when we sweep over difference frequencies

                !Temporary storage for a normal ordered list of the new length
                allocate(nI_GF(nel),stat=ierr)
                if(ierr.ne.0) call stop_all(t_r,"Alloc err")
                    
                root_write(6,"(A)") "Setup finished"
                
                !Now, we want to loop over the different frequencies that we want to calc
                if(GFType.eq.1) then
                    Freq = StartFreq_adv - dFreqStep_adv
                    nFreqSteps = nFreqPnts_adv
                    write(6,*) "Calculating freq-domain GF between values ",StartFreq_adv," and ",FinalFreq_adv
                else
                    Freq = StartFreq_ret - dFreqStep_ret
                    nFreqSteps = nFreqPnts_ret
                    write(6,*) "Calculating freq-domain GF between values ",StartFreq_ret," and ",FinalFreq_ret
                endif
                do nFreq=0,nFreqSteps
                    if(GFType.eq.1) then
                        Freq=Freq + dFreqStep_adv
                    else
                        Freq=Freq + dFreqStep_ret
                    endif
                    SumBetaSq = -(TauSq)/2  !We need to bias all overlaps by the amount of time elapsed (beta). This wants to be midway between integration points
                    root_write (6,"(A,I8,A,I8)") "Taking frequency point: ",nFreq," / ",nFreqSteps
                    root_write (6,"(A,F20.10)") "Fixing total shift value at: ",Freq
                    DiagSft=Freq-Hii    !This actually needs to be offset by the reference energy... 
                    root_write(6,"(A,F20.10,A,F20.10)") "With a reference energy of ",Hii,  &
                        " this is a shift of ",DiagSft
                    call ReloadStartingVec()    !Reload the N +/- 1 electron wavefunction data

                    root_write(6,"(A,I10,A,I10)") "Setup finished, running for total iterations: ",     &
                        GreensIters(GFType)," and calculating integral every ",GreenTimeslice(GFType)

                    TimeNumber=0
                    integrate: do iter=0,GreensIters(GFType)
                        !Call FCIQMC for a certain amout of imaginary time to loop over time range for GF
                        SumBetaSq = SumBetaSq + TauSq   !Calculate the amount of beta^2 we have gone through

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
                                root_write(6,"(A,I8,2I12,G20.10,I8)") "N-1 Adv Iter ",iter,sum(AllTotParts),   &
                                    AllTotWalkers,Freq,Sample
                            else
                                root_write (6,"(A,I8,2I12,G20.10,I8)") "N+1 Adv Iter ",iter,sum(AllTotParts),   &
                                    AllTotWalkers,Freq,Sample
                            endif
                        endif
                    
                        if(mod(iter,GreenTimeslice(GFType)).eq.0) then
                            !Every t timeslices, run over all i's and get overlap 
                            if(iter.ne.0) TimeNumber=TimeNumber+1
                            if(iter.eq.0.and.tFinalOverlapOnly) cycle

                            IFDEBUGTHEN(GFDebug,1)
                                if(GFType.eq.1) then
                                    root_write (6,"(A,I6,A,I6)") "Taking overlap for advanced integral: ",       &
                                            TimeNumber,' / ',nTimePnts_adv
                                else
                                    root_write (6,"(A,I6,A,I6)") "Taking overlap for retarded integral: ",       &
                                            TimeNumber,' / ',nTimePnts_ret
                                endif
                            ENDIFDEBUG

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
                                
                                !Now find overlap, and sum into the integral
                                OverlapNew(:) = 0.0_dp

    !We must all loop for the same number of times, so that processors are synchronized and can communicate when necessary
                                call MPIAllReduce(TotWalkers,MPI_MAX,MaxWalkersProc)
    !                            write(6,*) "Max Determinant number is: ",MaxWalkersProc

                                ValidSpawnedList = InitialSpawnedSlots
                                do walk=1,MaxWalkersProc
                                    !Run through determinant list on this core
                                    !Here, we apply the annihilation/creation operator on orbital i to the sampled wavefunction

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
                                                    TempDet((i-1)/bits_n_int)=ibset(TempDet((i-1)/bits_n_int),  &
                                                        mod(i-1,bits_n_int))
!                                                    set_orb(TempDet,i)
                                                endif
                                            else
                                                !Check it is not occupied
                                                if(.not.btest(CurrentDets((i-1)/bits_n_int,walk),mod(i-1,bits_n_int))) then
                                                    tSendDet=.false.
                                                else
                                                    tSendDet=.true.
                                                    !i is occupied - distroy it
                                                    TempDet(:)=CurrentDets(:,walk)
                                                    TempDet((i-1)/bits_n_int)=ibclr(TempDet((i-1)/bits_n_int),  &
                                                        mod(i-1,bits_n_int))
!                                                    clr_orb(TempDet,i)
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
                                                call GetBitExcitation(SpawnedParts(0:NIfD,k),       &
                                                    CurrentDetsSaved(0:NIfD,PartInd),Ex,tSign)
                                                if(tSign) SignRight(:)=-SignRight(:)  !Flip the sign of one of them

                                                !Sum this determinants contribution into the new overlap
                                                if(lenof_sign.eq.1) then
                                                    OverlapNew(1) = OverlapNew(1) + real(SignRight(1)*SignLeft(1),dp)
                                                else
                                                    !Overlap of complex weights
                                                    OverlapNew(1) = OverlapNew(1) + real(SignRight(1)*SignLeft(1) - &
                                                        SignRight(lenof_sign)*SignLeft(lenof_sign),dp)
                                                    OverlapNew(lenof_sign) = OverlapNew(lenof_sign) +   &
                                                        real(SignRight(1)*SignLeft(lenof_sign) + &
                                                        SignRight(lenof_sign)*SignLeft(1),dp)
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
                                !Now we have the determinants on the correct processors
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
                                        !Sum this determinants contribution into the new overlap
                                        if(lenof_sign.eq.1) then
                                            OverlapNew(1) = OverlapNew(1) + real(SignRight(1)*SignLeft(1),dp)
                                        else
                                            !Overlap of complex weights
                                            OverlapNew(1) = OverlapNew(1) + real(SignRight(1)*SignLeft(1) - &
                                                SignRight(lenof_sign)*SignLeft(lenof_sign),dp)
                                            OverlapNew(lenof_sign) = OverlapNew(lenof_sign) +   &
                                                real(SignRight(1)*SignLeft(lenof_sign) + &
                                                SignRight(lenof_sign)*SignLeft(1),dp)
                                        endif

                                    endif

                                enddo
                                ValidSpawnedList = InitialSpawnedSlots
                                AllOverlapNew(:) = 0.0_dp
                                call MPISumAll(OverlapNew,AllOverlapNew)
                                AllOverlapNew(:) = AllOverlapNew(:)/AllNorm !Because we never actually take the sqrt of the norm, this is actually the norm^2
                                if(TimeNumber.eq.0) then
                                    !This is simply the initial point, having not actually gone through any iterations yet.
                                    !This should just be the overlap <0|a_ia_j|0> (an element of the 1RDM)
                                    !Don't actually sum this in, but it'll be used for the next point.
                                    IFDEBUGTHEN(GFDebug,2) 
                                        write(6,*) "Initial overlap of <0|a_ia_j|0> is: ",AllOverlapNew
                                    ENDIFDEBUG
                                else
                                    if(GFType.eq.1) then
                                        if(tFinalOverlapOnly) then
                                            if(iter.ne.GreensIters(GFType)) call stop_all(t_r,"Not taking final overlap")
                                            G_adv_loc(:,i,j,nFreq) = AllOverlapNew(:)*sqrt(SumBetaSq)
                                        else
                                            G_adv_loc(:,i,j,nFreq) = G_adv_loc(:,i,j,nFreq) +   &
                                                (((AllOverlapNew(:)+OldOverlap(:,i,j))/2.0_dp)*sqrt(SumBetaSq))
                                        endif
                                    else
                                        if(tFinalOverlapOnly) then
                                            if(iter.ne.GreensIters(GFType)) call stop_all(t_r,"Not taking final overlap")
                                            G_ret_loc(:,i,j,nFreq) = AllOverlapNew(:)*sqrt(SumBetaSq)
                                        else
                                            G_ret_loc(:,i,j,nFreq) = G_ret_loc(:,i,j,nFreq) +   &
                                                (((AllOverlapNew(:)+OldOverlap(:,i,j))/2.0_dp)*sqrt(SumBetaSq))
                                        endif
                                    endif
                                    IFDEBUGTHEN(GFDebug,2)
                                        if(lenof_sign.eq.1) then
                                            write(6,"(A,2I4,A,I6,A,F18.10)") "Overlap for transistion: ",i,j,   &
                                                " and frequency point: ",nFreq," is: ",AllOverlapNew(:)
                                        else
                                            write(6,"(A,2I4,A,I6,A,2F18.10)") "Overlap for transistion: ",i,j,   &
                                                " and frequency point: ",nFreq," is: ",AllOverlapNew(:)
                                        endif
                                    ENDIFDEBUG

                                endif
                                OldOverlap(:,i,j) = AllOverlapNew(:)

                                if(AllTotWalkers.eq.0) then
                                    write(6,*) "All walkers died. Exiting this frequency and transistion loop"
                                    exit integrate !No point in carrying on if there are no particles left.
                                endif

                            enddo   !Enddo loop over i

                        endif   !Endif integration datapoint

                    enddo integrate   !Enddo GF integration range

                    if(AllTotWalkers.ne.0) then
                        write(6,*) "Did not manage to integrate this frequency point to infinity"
                        write(6,*) "Remaining walkers: ",AllTotParts
                    endif
                    !This should be roughly equivalent to 1/delta
                    write(6,*) "Finished integration, which was over a quantity beta of: ",sqrt(SumBetaSq)

                enddo   !Enddo over frequency

                !Return electron number to full N system.
                call ReturnElectronNumber(GFType)
                    
                deallocate(nI_GF)

            enddo   !End retarded/advanced
            tBackwardsTime=.false.
            
            if(tSpecificMatEls.and.(jval.gt.Contributing_Js)) exit  !We have run through all allowed j's

        enddo   !Enddo orbital j

        do i=1,nBasis
            do j=1,nBasis
                if(.not.tJustRetGF) then
                    do k=0,nFreqPnts_adv
                        if(tFinalOverlapOnly) then
                            G_adv_loc(:,i,j,k) = G_adv_loc(:,i,j,k) /1.77245385091_dp
                        else
                            G_adv_loc(:,i,j,k) = G_adv_loc(:,i,j,k) * (tauSq*GreenTimeslice(1)/1.77245385091_dp)   !Number is sqrt(pi)
                        endif
                        G_adv(:,i,j,k)=G_adv(:,i,j,k)+G_adv_loc(:,i,j,k)  !Sum into averaged value
                        G_adv_sq(:,i,j,k)=G_adv_sq(:,i,j,k)+G_adv_loc(:,i,j,k)**2.0_dp  !Square contribution
                    enddo
                endif
                if(.not.tJustAdvGF) then
                    do k=0,nFreqPnts_ret
                        if(tFinalOverlapOnly) then
                            G_ret_loc(:,i,j,k) = G_ret_loc(:,i,j,k) /1.77245385091_dp   !Number is sqrt(pi)
                        else
                            G_ret_loc(:,i,j,k) = G_ret_loc(:,i,j,k) * (tauSq*GreenTimeslice(2)/1.77245385091_dp)   !Number is sqrt(pi)
                        endif
                        G_ret(:,i,j,k)=G_ret(:,i,j,k)+G_ret_loc(:,i,j,k)  !Sum into averaged value
                        G_ret_sq(:,i,j,k)=G_ret_sq(:,i,j,k)+G_ret_loc(:,i,j,k)**2.0_dp  !Square contribution
                    enddo
                endif
            enddo
        enddo

        tSaveHashList=.false.   !Regenerate the lists to avoid correlation?
        call DeallocFCIMCMemPar()
        deallocate(CurrentDetsSaved)
        if(tHashWalkerList) deallocate(HashIndexSaved)
    
        
        !Sum in parallel
        if(.not.tJustRetGF) then
            G_adv_all(:,:,:,:)=G_adv(:,:,:,:)
        endif
        if(.not.tJustAdvGF) then
            G_ret_all(:,:,:,:)=G_ret(:,:,:,:)
        endif
        
        !Find average now by dividing by GFSamples
        do i=1,nBasis
            do j=1,nBasis
                if(.not.tJustRetGF) then
                    do k=0,nFreqPnts_adv
                        G_adv_all(:,i,j,k)=G_adv_all(:,i,j,k)/real(allsamples,dp)
!                        G_adv_sq(:,i,j,k)=G_adv_sq(:,i,j,k)/real(GFSamples,dp)
                    enddo
                endif
                if(.not.tJustAdvGF) then
                    do k=0,nFreqPnts_ret
                        G_ret_all(:,i,j,k)=G_ret_all(:,i,j,k)/real(allsamples,dp)
!                        G_ret_sq(:,i,j,k)=G_ret_sq(:,i,j,k)/real(GFSamples,dp)
                    enddo
                endif
            enddo
        enddo
        !Write out greens function
        if(iProcIndex.eq.Root) then
            jval=1
            do j=1,nBasis
                if(tSpecificMatEls) then
                    if(jval.gt.iSpecificMatEls) cycle   !This will now just skip to the end
                    if(Allowed_j(jval).ne.j) cycle  !No desired i's here
                    jval=jval+1
                endif
                do i=1,nBasis

                    if(.not.CheckAllowed_i(j,i)) cycle   !Not an allowed i orbital with this j
            
                    PairIndex=FindPairIndex(j,i)

                    if(.not.tJustRetGF) then
                        if(UnitNos(1,PairIndex).eq.-1) call stop_all(t_r,"Unit not open for writing")
                        rewind(unit=UnitNos(1,PairIndex))
                        if(lenof_sign.eq.1) then
                            write(UnitNos(1,PairIndex),"(A)") "# Freq  Av.GF  Av.GF_sq    Error    Samples"
                        else
                            write(UnitNos(1,PairIndex),"(A)") "# Freq  Av.GF_Re  Av.GF_Im  Av.GF_sq_Re  ",  &
                                " Av.GF_sq_Im   Error_Re   Error_Im   Samples"
                        endif
                        Freq = StartFreq_adv-dFreqStep_adv
                        do k=0,nFreqPnts_adv
                            Freq = Freq + dFreqStep_adv
                            if(lenof_sign.eq.1) then
                                write(UnitNos(1,PairIndex),"(4G25.15,I7)") Freq,G_adv_all(1,i,j,k),    &
                                    G_adv_sq(1,i,j,k)/real(allsamples,dp),   &
                                    sqrt(((G_adv_sq(1,i,j,k)/allsamples)-(G_adv_all(1,i,j,k)**2.0_dp))/real(allsamples,dp)),    &
                                    allsamples
                            else
                                write(UnitNos(1,PairIndex),"(7G25.15,I7)") Freq,G_adv_all(1,i,j,k),G_adv_all(2,i,j,k), &
                                    G_adv_sq(1,i,j,k)/real(allsamples), &
                                    (G_adv_sq(2,i,j,k)/real(allsamples,dp)),     &
                                    sqrt(((G_adv_sq(1,i,j,k)/real(allsamples,dp))-(G_adv_all(1,i,j,k)**2.0_dp))/real(allsamples,dp)),    &
                                    sqrt(((G_adv_sq(2,i,j,k)/real(allsamples,dp))-(G_adv_all(2,i,j,k)**2.0_dp))/real(allsamples,dp)),   &
                                    allsamples
                            endif
                        enddo
                        call neci_flush(UnitNos(1,PairIndex))
                    endif
                    if(.not.tJustAdvGF) then
                        if(UnitNos(2,PairIndex).eq.-1) call stop_all(t_r,"Unit not open for writing")
                        rewind(unit=UnitNos(2,PairIndex))
                        if(lenof_sign.eq.1) then
                            write(UnitNos(2,PairIndex),"(A)") "# Freq  Av.GF  Av.GF_sq    Error    Samples"
                        else
                            write(UnitNos(1,PairIndex),"(A)") "# Freq  Av.GF_Re  Av.GF_Im  Av.GF_sq_Re  ",  &
                                " Av.GF_sq_Im   Error_Re   Error_Im   Samples"
                        endif
                        Freq = StartFreq_adv-dFreqStep_adv
                        do k=1,nFreqPnts_ret
                            Freq = Freq + dFreqStep_adv
                            if(lenof_sign.eq.1) then
                                write(UnitNos(2,PairIndex),"(4G25.15,I7)") Freq,G_ret_all(1,i,j,k),    &
                                    (G_ret_sq(1,i,j,k)/real(allsamples,dp)),   &
                                    sqrt(((G_ret_sq(1,i,j,k)/real(allsamples,dp))-(G_ret_all(1,i,j,k)**2.0_dp))/real(allsamples,dp)), &
                                    allsamples
                            else
                                write(UnitNos(2,PairIndex),"(7G25.15,I7)") Freq,G_ret_all(1,i,j,k),G_ret_all(2,i,j,k), &
                                    (G_ret_sq(1,i,j,k)/real(allsamples,dp)), &
                                    (G_ret_sq(2,i,j,k)/real(allsamples,dp)),    &
                                    sqrt(((G_ret_sq(1,i,j,k)/real(allsamples,dp))-(G_ret_all(1,i,j,k)**2.0_dp))/real(allsamples,dp)),    &
                                    sqrt(((G_ret_sq(2,i,j,k)/real(allsamples,dp))-(G_ret_all(2,i,j,k)**2.0_dp))/real(allsamples,dp)),   &
                                    allsamples
                            endif
                        enddo
                        call neci_flush(UnitNos(2,PairIndex))
                    endif
                enddo
                if(tSpecificMatEls.and.(jval.gt.Contributing_Js)) exit  !We have run through all allowed j's
            enddo

            !Write out current stats so that they can be read in again later if necessary
            !We need to write out the current value of allsamples, as well as all the values of G_adv and G_adv_sq
            tempunit = get_free_unit()
            if(.not.tJustRetGF) then
                open(unit=tempunit,file='GFRestart_adv',status='unknown')
                write(tempunit,*) allsamples," allsamples"
                jval=1
                do j=1,nBasis
                    if(tSpecificMatEls) then
                        if(jval.gt.iSpecificMatEls) cycle   !This will now just skip to the end
                        if(Allowed_j(jval).ne.j) cycle  !No desired i's here
                        jval=jval+1
                    endif
                    do i=1,nBasis
                        if(.not.CheckAllowed_i(j,i)) cycle   !Not an allowed i orbital with this j
                        do k=0,nFreqPnts_adv
                            if(lenof_sign.eq.1) then
                                write(tempunit,"(3I6,2G25.15)") i,j,k,G_adv(1,i,j,k),G_adv_sq(1,i,j,k)
                            else
                                write(tempunit,"(3I6,4G25.15)") i,j,k,G_adv(1,i,j,k),G_adv(lenof_sign,i,j,k), &
                                    G_adv_sq(1,i,j,k),G_adv_sq(lenof_sign,i,j,k)
                            endif
                        enddo
                    enddo
                enddo
                close(tempunit)
            endif
            if(.not.tJustAdvGF) then
                open(unit=tempunit,file='GFRestart_ret',status='unknown')
                write(tempunit) allsamples," allsamples"
                jval=1
                do j=1,nBasis
                    if(tSpecificMatEls) then
                        if(jval.gt.iSpecificMatEls) cycle   !This will now just skip to the end
                        if(Allowed_j(jval).ne.j) cycle  !No desired i's here
                        jval=jval+1
                    endif
                    do i=1,nBasis
                        if(.not.CheckAllowed_i(j,i)) cycle   !Not an allowed i orbital with this j
                        do k=0,nFreqPnts_adv
                            if(lenof_sign.eq.1) then
                                write(tempunit,"(3I6,2G25.15)") i,j,k,G_ret(1,i,j,k),G_ret_sq(1,i,j,k)
                            else
                                write(tempunit,"(3I6,4G25.15)") i,j,k,G_ret(1,i,j,k),G_ret(lenof_sign,i,j,k), &
                                    G_ret_sq(1,i,j,k),G_ret_sq(lenof_sign,i,j,k)
                            endif
                        enddo
                    enddo
                enddo
                close(tempunit)
            endif

        endif

    enddo   !Enddo samples

    !Sum in parallel
    if(.not.tJustRetGF) then
        G_adv_all(:,:,:,:)=G_adv(:,:,:,:)
    endif
    if(.not.tJustAdvGF) then
        G_ret_all(:,:,:,:)=G_ret(:,:,:,:)
    endif
    
    !Find average now by dividing by GFSamples
    do i=1,nBasis
        do j=1,nBasis
            if(.not.tJustRetGF) then
                do k=0,nFreqPnts_adv
                    G_adv_all(:,i,j,k)=G_adv_all(:,i,j,k)/real(allsamples,dp)
                    G_adv_sq(:,i,j,k)=G_adv_sq(:,i,j,k)/real(allsamples,dp)
                enddo
            endif
            if(.not.tJustAdvGF) then
                do k=0,nFreqPnts_ret
                    G_ret_all(:,i,j,k)=G_ret_all(:,i,j,k)/real(allsamples,dp)
                    G_ret_sq(:,i,j,k)=G_ret_sq(:,i,j,k)/real(allsamples,dp)
                enddo
            endif
        enddo
    enddo

    !Write out greens function
    if(iProcIndex.eq.Root) then
        jval=1
        do j=1,nBasis
            if(tSpecificMatEls) then
                if(jval.gt.iSpecificMatEls) cycle   !This will now just skip to the end
                if(Allowed_j(jval).ne.j) cycle  !No desired i's here
                jval=jval+1
            endif
            do i=1,nBasis

                if(.not.CheckAllowed_i(j,i)) cycle   !Not an allowed i orbital with this j
        
                PairIndex=FindPairIndex(j,i)

                if(.not.tJustRetGF) then
                    if(UnitNos(1,PairIndex).eq.-1) call stop_all(t_r,"Unit not open for writing")
                    rewind(unit=UnitNos(1,PairIndex))
                    if(lenof_sign.eq.1) then
                        write(UnitNos(1,PairIndex),"(A)") "# Freq  Av.GF  Av.GF_sq    Error   Samples "
                    else
                        write(UnitNos(1,PairIndex),"(A)") "# Freq Av.GF_Re Av.GF_Im Av.GF_sq_Re Av.GF_Im Error_Re Error_Im Samples"
                    endif
                    Freq = StartFreq_adv-dFreqStep_adv
                    do k=0,nFreqPnts_adv
                        Freq = Freq + dFreqStep_adv
                        if(lenof_sign.eq.1) then
                            write(UnitNos(1,PairIndex),"(4G25.15,I7)") Freq,G_adv_all(1,i,j,k),G_adv_sq(1,i,j,k),   &
                                sqrt((G_adv_sq(1,i,j,k)-(G_adv_all(1,i,j,k)**2.0_dp))/real(allsamples,dp)), &
                                allsamples
                        else
                            write(UnitNos(1,PairIndex),"(7G25.15,I7)") Freq,G_adv_all(1,i,j,k),G_adv_all(lenof_sign,i,j,k),    &
                                G_adv_sq(1,i,j,k), &
                                G_adv_sq(lenof_sign,i,j,k), &
                                sqrt((G_adv_sq(1,i,j,k)-(G_adv_all(1,i,j,k)**2.0_dp))/real(allsamples,dp)),    &
                                sqrt((G_adv_sq(lenof_sign,i,j,k)-(G_adv_all(lenof_sign,i,j,k)**2.0_dp))/real(allsamples,dp)),   &
                                allsamples
                        endif
                    enddo
                    close(UnitNos(1,PairIndex))
                endif
                if(.not.tJustAdvGF) then
                    if(UnitNos(2,PairIndex).eq.-1) call stop_all(t_r,"Unit not open for writing")
                    rewind(unit=UnitNos(2,PairIndex))
                    write(UnitNos(2,PairIndex),"(A)") "# Freq  Av.GF  Av.GF_sq    Error    samples"
                    Freq = StartFreq_adv-dFreqStep_adv
                    do k=1,nFreqPnts_ret
                        Freq = Freq + dFreqStep_adv
                        if(lenof_sign.eq.1) then
                            write(UnitNos(2,PairIndex),"(4G25.15,I7)") Freq,G_ret_all(1,i,j,k),G_ret_sq(1,i,j,k),   &
                                sqrt((G_ret_sq(1,i,j,k)-(G_ret_all(1,i,j,k)**2.0_dp))/real(allsamples,dp)), &
                                allsamples
                        else
                            write(UnitNos(2,PairIndex),"(7G25.15,I7)") Freq,G_ret_all(1,i,j,k),G_ret_all(2,i,j,k),G_ret_sq(1,i,j,k), &
                                G_ret_sq(2,i,j,k),sqrt((G_ret_sq(1,i,j,k)-(G_ret_all(1,i,j,k)**2.0_dp))/real(allsamples,dp)),    &
                                sqrt((G_ret_sq(2,i,j,k)-(G_ret_all(2,i,j,k)**2.0_dp))/real(allsamples,dp)), &
                                allsamples
                        endif
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

    if(tSpecificMatEls) deallocate(SpecificMatEls,Allowed_j)

    end subroutine CalcGreensFuncFreq
                
    !Go from exponential propagator, to gaussian
    subroutine SetupForGaussProp()
        use FciMCData, only: SingleHIterH,SingleHAppH,SingleHIterHTag,MaxSpawned
        integer :: ierr
        character(len=*), parameter :: t_r="SetupForGaussProp"

        write(6,*) "Setting up for gaussian propagation"

        tMultipleH = .true.
        if(iMaxH_App.ne.1) call stop_all(t_r,"Should be exponential propagation")
        iMaxH_App = 2   !Gaussian prop
        if(.not.tRegenDiagHEls) then
            if(.not.allocated(SingleHIterH)) then
                !Allocate space for intermediate propagation 
                allocate(SingleHIterH(MaxSpawned),stat=ierr)
                call LogMemAlloc('SingleHIterH',MaxSpawned,8,t_r,SingleHIterHTag,ierr)
                SingleHIterH(:)=0.0_dp
            endif
            SingleHAppH => SingleHIterH
        endif
        if(lenof_sign.eq.2) then
            tImagEnergyShift = .true.   !Turn on imaginary energy shift TODO
        endif
        if(DriveFactor.ne.1.0_dp) then
            tDrivePropagation = .true.
            write(6,*) "Setting drivefactor to :",DriveFactor
        endif
        SftDampDrive = 0.0_dp   !Keep driving fixed

    end subroutine SetupForGaussProp

    !Save the a_i|0> starting vector to make it easier to loop over frequency range
    subroutine SaveStartingVec()
        implicit none
        integer :: ierr
        character(len=*), parameter :: t_r='SaveStartingVec'

        !Walkers are initially stored contiguously, so we can just store the walkers, even with
        !tHashWalkerList.
        if(allocated(SavedInitialVec)) deallocate(SavedInitialVec)
        allocate(SavedInitialVec(0:NIfTot,TotWalkers),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"alloc err")
        SavedInitialVec(:,1:TotWalkers) = CurrentDets(:,1:TotWalkers) 

        if(allocated(SavedCurrentH)) deallocate(SavedCurrentH)
        allocate(SavedCurrentH(TotWalkers),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"alloc err")
        SavedCurrentH(1:TotWalkers) = CurrentH(1:TotWalkers)

        if(tHashWalkerList) then
            if(allocated(SavedHashIndex)) deallocate(SavedHashIndex)
            allocate(SavedHashIndex(0:nClashMax,nWalkerHashes),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"alloc err")
            SavedHashIndex(0:nClashMax,:) = HashIndex(0:nClashMax,:)
            SavednClashMax = nClashMax
        endif

        SavedTotWalkers = TotWalkers
        SavedTotParts = TotParts
        SavedNoatHF = NoatHF

    end subroutine SaveStartingVec

    !Routine to reload the starting vector
    subroutine ReloadStartingVec()
        implicit none
        integer :: ierr
        character(len=*), parameter :: t_r='ReloadStartingVec'

        TotWalkers = SavedTotWalkers
        TotWalkersOld = TotWalkers           
        TotParts = SavedTotParts
        TotPartsOld = TotParts
        NoatHF = SavedNoatHF

        call MPISum(TotParts,AllTotParts)
        call MPISum(NoatHF,AllNoatHF)
        call MPISum(TotWalkers,AllTotWalkers)
        AllTotWalkersOld = AllTotWalkers
        AllTotPartsOld = AllTotParts
        OldAllNoatHF = AllNoatHF
        OldAllAvWalkersCyc = real(sum(AllTotParts),dp)
        if(lenof_sign.eq.1) then
            OldAllHFCyc = real(AllNoatHF(1),dp)
        else
            OldAllHFCyc = cmplx(real(AllNoatHF(1),dp),real(AllNoatHF(lenof_sign),dp), dp)
        endif

        if(tHashWalkerList) then
            nClashMax = SavednClashMax
            if(allocated(HashIndexArr1)) then
                deallocate(HashIndexArr1)
                nullify(HashIndex)
                allocate(HashIndexArr1(0:nClashMax,nWalkerHashes),stat=ierr)
                if(ierr.ne.0) call stop_all(t_r,"alloc err")
                HashIndex => HashIndexArr1
            else
                deallocate(HashIndexArr2)
                nullify(HashIndex)
                allocate(HashIndexArr2(0:nClashMax,nWalkerHashes),stat=ierr)
                if(ierr.ne.0) call stop_all(t_r,"alloc err")
                HashIndex => HashIndexArr2
            endif
            HashIndex(:,:) = SavedHashIndex(:,:)
            FreeSlot(:) = 0
        endif

        CurrentDets(:,:) =0
        CurrentDets(0:NIfTot,1:TotWalkers) = SavedInitialVec(0:NIfTot,1:TotWalkers)

        CurrentH(1:TotWalkers) = SavedCurrentH(1:TotWalkers)

    end subroutine ReloadStartingVec

end module GreensFuncFreqMod

subroutine CalcGreensFuncFreq_wrap()
    use GreensFuncFreqMod, only: CalcGreensFuncFreq
    implicit none

    call CalcGreensFuncFreq()

end subroutine CalcGreensFuncFreq_wrap

