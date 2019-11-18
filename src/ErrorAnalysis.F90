#include "macros.h"
module errors
    use constants
    use SystemData, only: BasisFN, nBasisMax, G1, nel
    use FciMCData, only: tSinglePartPhase, ProjectionE, iBlockingIter
    use Determinants, only: get_helement
    use DeterminantData, only: fdet
    use Parallel_neci
    use LoggingData, only: iBlockEquilProjE, iBlockEquilShift
    use util_mod, only: get_free_unit
    use CalcData, only: SftDamp
    implicit none
    real(dp), allocatable :: numerator_data(:)
    real(dp), allocatable :: imnumerator_data(:)
    real(dp), allocatable :: pophf_data(:)
    real(dp), allocatable :: shift_data(:)
    integer :: Errordebug
    logical :: tGivenEquilibrium
    integer :: Shift_Start, ProjE_Start
    integer :: relaxation_time_shift, relaxation_time_proje

    contains

    !Perform automatic FCIMC blocking analysis, reading in whether we have started varying the shift or not,
    !and which iteration we did so, and
    !returning the mean and error for all energy estimates.
    subroutine error_analysis(tSingPartPhase,iShiftVary,mean_ProjE_re,ProjE_Err_re, &
        mean_ProjE_im,ProjE_Err_im,mean_Shift,Shift_Err,tNoProjEValue,tNoShiftValue, &
        equilshift,equilproje)
        implicit none
        logical, intent(in) :: tSingPartPhase
        integer , intent(inout) :: iShiftVary
        integer , intent(in) , optional :: equilshift,equilproje
        integer :: corrlength_denom,corrlength_num,corrlength_shift,corrlength_imnum
        integer :: equil_time_denom,equil_time_num,equil_time_shift,equil_time_imnum
        real(dp), allocatable :: temp(:)
        integer :: final_corrlength_proj,ierr
        real(dp) :: mean_denom, error_denom, eie_denom
        real(dp) :: mean_num, error_num, eie_num
        real(dp) :: eie_shift
        real(dp) :: mean_imnum, error_imnum, eie_imnum
        real(dp) :: covariance_re, correction_re
        real(dp) :: covariance_im, correction_im
        logical :: tFail_num,tFail_denom,tFail_imnum,tFail_shift
        logical :: tPrintIntermediateBlocking
        logical, intent(out) :: tNoProjEValue,tNoShiftValue
        real(dp), intent(out) :: mean_ProjE_re,ProjE_Err_re,mean_ProjE_im,ProjE_Err_im,mean_Shift,Shift_Err
        character(len=*), parameter :: t_r='Error_analysis'
        logical :: tFailRead

        if(iProcIndex.ne.Root) return   !Just do this in serial

        write(6,"(A)")
        write(6,"(A)") "Calculating approximate errorbar for energy estimates..."
        write(6,"(A)")

        if(present(equilshift).or.present(equilproje)) then
            tGivenEquilibrium = .true.
            if(present(equilshift)) then
                Shift_Start = equilshift
            else
                Shift_Start = 0
            endif
            if(present(equilproje)) then
                ProjE_Start = equilproje
            else
                ProjE_Start = 0
            endif
            iShiftVary = min(ProjE_Start,Shift_Start)
        else
            tGivenEquilibrium = .false.
        endif

        tNoProjEValue = .false.
        tNoShiftValue = .false.
        if((tSingPartPhase.and.(.not.tGivenEquilibrium)).or.(SftDamp.lt.1.0e-7_dp)) then
            write(6,"(A)") "Calculation has not entered variable shift phase. Error analysis therefore not performed."
            write(6,"(A)") "Direct reblocking of instantaneous energy possible, but this would contain finite sampling bias."
            tNoProjEValue = .true.
            tNoShiftValue = .true.
            return
        else
            write(6,"(A,I12)") "Attempting automatic reblocking analysis on data from iteration: ",iShiftVary
        endif

        !Read and store numerator_data, pophf_data, shift_data, and if complex, also imnumerator_data
        call read_fcimcstats(iShiftVary,tFailRead)
        if(tFailRead) then
            tNoProjEValue = .true.
            tNoShiftValue = .true.
            return
        endif

        ! STEP 2) Performs an preliminary blocking analysis, with a fixed blocklength of 2
        ! This is predominantly to yield the correlation length
        ! Files for plotting (for manual reblocking) are outputted in PLOT files
        if(Errordebug.gt.0) then
            tPrintIntermediateBlocking=.true.
        else
            tPrintIntermediateBlocking=.false.
        endif
        call automatic_reblocking_analysis(pophf_data,2,corrlength_denom,   &
            'PLOT_denom',tPrintIntermediateBlocking,1)
        call automatic_reblocking_analysis(numerator_data,2,corrlength_num, &
            'PLOT_num',tPrintIntermediateBlocking,2)
        call automatic_reblocking_analysis(shift_data,2,corrlength_shift,   &
            'PLOT_shift',tPrintIntermediateBlocking,3)
        if(lenof_sign > inum_runs) then
            call automatic_reblocking_analysis(imnumerator_data,2,  &
                corrlength_imnum,'PLOT_imnum',tPrintIntermediateBlocking,4)
        endif

        if(.not.tGivenEquilibrium) then

            ! STEP 3) Attempts an automatic removal of equilibration time
            ! From the data that's had serial correlation removed
            ! (See subroutine for how this is achieved)
            allocate(temp(size(pophf_data)),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,'Alloc error')
            temp = pophf_data
            call reblock_data(temp,corrlength_denom)
            if(errordebug.gt.0) call print_vector(temp,filename='PLOT_blocked_pophf')
            call attempt_detect_equil_time(temp,equil_time_denom,tFail_denom)
            if(tFail_denom) then
                write(6,"(A)")
                write(6,"(A)") "*** ERROR *** Failure to automatically detect equilibration time for " &
                    & //"projected energy denominator"
                write(6,"(A)") "Skipping blocking analysis of projected energy, and energy estimate " &
                    & //"will be simple average over "
                write(6,"(A)") "all iterations (including growth phase), which may contain correlated " &
                    & //"sampling bias. Use with caution."
                write(6,"(A)") "Manual reblocking or continued running suggested for accurate " &
                    & //"projected energy estimate."
                tNoProjEValue = .true.
            endif
            deallocate(temp)

            allocate(temp(size(numerator_data)),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,'Alloc error')
            temp = numerator_data
            call reblock_data(temp,corrlength_num)
            if(errordebug.gt.0) call print_vector(temp,filename='PLOT_blocked_num')
            call attempt_detect_equil_time(temp,equil_time_num,tFail_num)
            if(tFail_num) then
                write(6,"(A)")
                write(6,"(A)") "*** ERROR *** Failure to automatically detect equilibration time for " &
                    & //"projected energy numerator"
                write(6,"(A)") "Skipping blocking analysis of projected energy, and energy estimate " &
                    & //"will be simple average over "
                write(6,"(A)") "all iterations (including growth phase), which may contain correlated " &
                    & //"sampling bias. Use with caution."
                write(6,"(A)") "Manual reblocking or continued running suggested for accurate " &
                    & //"projected energy estimate."
                tNoProjEValue = .true.
            endif
            deallocate(temp)

            allocate(temp(size(shift_data)),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,'Alloc error')
            temp = shift_data
            call reblock_data(temp,corrlength_shift)
            if(errordebug.gt.0) call print_vector(temp,filename='PLOT_blocked_shift')
            call attempt_detect_equil_time(temp,equil_time_shift,tFail_shift)
            if(tFail_shift) then
                write(6,"(A)")
                write(6,"(A)") "*** ERROR *** Failure to automatically detect equilibration time for shift value."
                write(6,"(A)") "Skipping blocking analysis and calculation of average shift."
                write(6,"(A)") "Continued running suggested if accurate shift estimate required."
                tNoShiftValue = .true.
            endif
            deallocate(temp)

            if(lenof_sign > inum_runs) then
                allocate(temp(size(imnumerator_data)),stat=ierr)
                if(ierr.ne.0) call stop_all(t_r,'Alloc error')
                temp = imnumerator_data
                call reblock_data(temp,corrlength_imnum)
                if(errordebug.gt.0) call print_vector(temp,filename='PLOT_blocked_imnum')
                call attempt_detect_equil_time(temp,equil_time_imnum,tFail_imnum)
                if(tFail_imnum) then
                    write(6,"(A)")
                    write(6,"(A)") "*** ERROR *** Failure to automatically detect equilibration time for " &
                        & //"imaginary projected energy numerator"
                    write(6,"(A)") "Skipping blocking analysis of projected energy, and energy estimate " &
                        & //"will be simple average over "
                    write(6,"(A)") "all iterations (including growth phase), which may contain correlated " &
                        & //"sampling bias. Use with caution."
                    write(6,"(A)") "Continued running suggested for accurate projected energy estimate."
                    tNoProjEValue = .true.
                endif
                deallocate(temp)
            endif

            ! STEP 4) Removes the equilibration time from the start of the data set
            ! Allow shift and projected energy to have seperate equilibration times
            ! If complex, the equilibration time for projected energy is the longest of all three quantities
            if(.not.tNoProjEValue) then
                relaxation_time_proje=max(equil_time_denom*corrlength_denom,equil_time_num*corrlength_num)
                if(lenof_sign > inum_runs) then
                    relaxation_time_proje=max(relaxation_time_proje,equil_time_imnum*corrlength_imnum)
                endif
                write(6,"(A,I8,A)") "Relaxation time for projected energy estimated to be ", relaxation_time_proje, &
                    " update cycles."
            endif
            if(.not.tNoShiftValue) then
                relaxation_time_shift = equil_time_shift*corrlength_shift
                write(6,"(A,I8,A)") "Relaxation time for shift estimated to be ", relaxation_time_shift,    &
                    " update cycles."
            endif

        else
            write(6,"(A,I8,A)") "Relaxation time for projected energy estimated to be ", relaxation_time_proje, &
                " update cycles."
            write(6,"(A,I8,A)") "Relaxation time for shift estimated to be ", relaxation_time_shift," update cycles."
        endif   !Endif automatic relaxation time calculation

        if(.not.tNoProjEValue) then
            call resize(pophf_data,relaxation_time_proje)
            call resize(numerator_data,relaxation_time_proje)
            if(lenof_sign > inum_runs) call resize(imnumerator_data,relaxation_time_proje)

            !Now find all block lengths, and write them to file, so that blocklengths can be checked.
            !Call routine here to write out a file (Blocks_proje) with the projected energy mean and error for all block sizes.
            if(lenof_sign > inum_runs) then
                call print_proje_blocks(pophf_data,numerator_data,2,'Blocks_proje_re')
#ifdef __CMPLX
                call print_proje_blocks(pophf_data,imnumerator_data,2,'Blocks_proje_im')
#endif
            else
                call print_proje_blocks(pophf_data,numerator_data,2,'Blocks_proje')
            endif

            !Now perform automatic reblocking again, to get expected blocklength
            call automatic_reblocking_analysis(pophf_data,2,corrlength_denom,'Blocks_denom',.true.,1)
            call automatic_reblocking_analysis(numerator_data,2,corrlength_num,'Blocks_num',.true.,2)
            if(lenof_sign > inum_runs) then
                call automatic_reblocking_analysis(imnumerator_data,2,corrlength_imnum,'Blocks_imnum',.true.,4)
            endif
        endif

        if(.not.tNoShiftValue) then
            call resize(shift_data,relaxation_time_shift)
            call automatic_reblocking_analysis(shift_data,2,corrlength_shift,'Blocks_shift',.true.,3)
        endif

        ! STEP 5) Now gathers together the properly reblocked data and find statistics
        if(.not.tNoProjEValue) then
            final_corrlength_proj = max(corrlength_num,corrlength_denom)
            if(lenof_sign > inum_runs) then
                final_corrlength_proj = max(final_corrlength_proj,corrlength_imnum)
            endif
            if(Errordebug.gt.0) write(6,"(A,I10)") "Final projected energy correlation length",final_corrlength_proj
            call reblock_data(pophf_data,final_corrlength_proj)
            call reblock_data(numerator_data,final_corrlength_proj)
            call analyze_data(pophf_data,mean_denom,error_denom,eie_denom)
            call analyze_data(numerator_data,mean_num,error_num,eie_num)
            write(6,"(A,F22.10,A,G20.8,A,G22.10)") "ProjE_denominator:", mean_denom, " +/- ", error_denom, &
                " Relative error: ", abs(error_denom/mean_denom)
            if(lenof_sign > inum_runs) then
                call reblock_data(imnumerator_data,final_corrlength_proj)
                call analyze_data(imnumerator_data,mean_imnum,error_imnum,eie_imnum)
                write(6,"(A,F22.10,A,G20.8,A,G22.10)") "ProjE_numerator (Re):", mean_num, " +/- ", error_num, &
                    " Relative error: ", abs(error_num/mean_num)
                write(6,"(A,F22.10,A,G20.8,A,G22.10)") "ProjE_numerator (Im):", mean_imnum, " +/- ", error_imnum, &
                    " Relative error: ", abs(error_imnum/mean_imnum)
            else
                write(6,"(A,F22.10,A,G20.8,A,G22.10)") "ProjE_numerator:  ", mean_num, " +/- ", error_num, &
                    " Relative error: ", abs(error_num/mean_num)
            endif
        endif
        ! write(6,*) "OVERALL", mean2/mean1, "+-",
        !abs(mean2/mean1)*((abs(error1/mean1))**2.0_dp+(abs(error2/mean2))**2.0_dp)**0.5_dp
        ! this is before the covariance, so don't print it now

        if(.not.tNoShiftValue) then
            !Now, also do this for the shift estimate
            call reblock_data(shift_data,corrlength_shift)
            call analyze_data(shift_data,mean_shift,shift_Err,eie_shift)
            if(Errordebug.gt.0) then
                write(6,"(A,F20.10,A,G20.8,A,G22.10)") "Shift:            ", &
                    mean_shift, " +/- ", shift_Err, " Relative error: ", abs(shift_Err/mean_shift)
            endif
        endif

        ! STEP 6) Refine statistics using covariance
        if(.not.tNoProjEValue) then
            covariance_re=calc_covariance(pophf_data,numerator_data)
            correction_re=2.0_dp*covariance_re/((size(pophf_data))*mean_denom*mean_num)
            if(lenof_sign > inum_runs) then
                covariance_im=calc_covariance(pophf_data,imnumerator_data)
                correction_im=2.0_dp*covariance_im/((size(pophf_data))*mean_denom*mean_imnum)
            endif

            mean_ProjE_re = mean_num/mean_denom
            ProjE_Err_re = abs(mean_ProjE_re)* sqrt( (abs(error_denom/mean_denom))**2.0_dp &
                + (abs(error_num/mean_num))**2.0_dp - correction_re )
            if(lenof_sign/inum_runs .eq. 2) then
                mean_ProjE_im = mean_imnum/mean_denom

                ProjE_Err_im = abs(mean_ProjE_im)* sqrt( (abs(error_denom/mean_denom))**2.0_dp &
                    + (abs(error_imnum/mean_imnum))**2.0_dp - correction_im )
                if(ErrorDebug.gt.0) then
                    write(6,"(A,F20.8)") "Covariance correction (Re):", correction_re
                    write(6,"(A,F20.8)") "Covariance correction (Im):", correction_im
                    write(6,"(A,F20.10,A,G20.8)") "Final projected energy (Re): ", mean_ProjE_re, " +/- ", ProjE_Err_re
                    write(6,"(A,F20.10,A,G20.8)") "Final projected energy (Im): ", mean_ProjE_im, " +/- ", ProjE_Err_im
                endif
            else
                if(ErrorDebug.gt.0) then
                    write(6,"(A,F20.8)") "Covariance correction:", correction_re
                    write(6,"(A,F20.10,A,G20.8)") "Final projected energy: ", mean_ProjE_re, " +/- ", ProjE_Err_re
                endif
            endif
        endif

        deallocate(pophf_data,numerator_data,shift_data)
        if(lenof_sign > inum_runs) deallocate(imnumerator_data)

    end subroutine error_analysis

    subroutine print_proje_blocks(hf_data,num_data,blocklength,filename)
        implicit none
        integer :: iunit,length,blocking_events,i
        real(dp), intent(in) :: hf_data(:),num_data(:)
        character(len=*), intent(in) :: filename
        integer,intent(in) :: blocklength ! this is the minimum step size (e.g. 2)
        real(dp) :: mean1, error1, eie1,mean2,error2,eie2
        real(dp) :: covariance,mean_proje,final_error,final_eie
        real(dp), allocatable :: this(:),that(:) ! stores this array

        iunit = get_free_unit()
        open(iunit,file=filename)
        write(iunit,"(A)") "# Block_length   Mean       Error        Error_in_Error"

        length=size(hf_data,1)
        allocate(this(length))
        allocate(that(length))
        that = hf_data  !Denominator data
        this = num_data !Numerator data

        call analyze_data(this,mean1,error1,eie1)   !means & error in num
        call analyze_data(that,mean2,error2,eie2)   !means & error in denom
        covariance=calc_covariance(that,this)
        mean_proje = mean1/mean2
        final_error=abs(mean_proje) * &
            sqrt(   (error2/mean2)**2.0_dp + (error1/mean1)**2.0_dp &
                    - 2.0_dp*covariance/((size(this))*mean1*mean2)  )
        final_eie = final_error/sqrt(2.0_dp*(size(this)-1))
        write(iunit,*) size(that), mean_proje, final_error, final_eie

        blocking_events=int(log(real(length,dp))/log(2.0_dp)-1.0_dp)
        do i=1,blocking_events
            call reblock_data(this,blocklength)
            call analyze_data(this,mean1,error1,eie1)   !means & error in num
            call reblock_data(that,blocklength)
            call analyze_data(that,mean2,error2,eie2)   !means & error in denom
            ! now have the correct means, but incorrect errors
            covariance=calc_covariance(that,this)
            mean_proje = mean1/mean2
            final_error=abs(mean_proje) * &
                sqrt( (error2/mean2)**2.0_dp + (error1/mean1)**2.0_dp &
                        - 2.0_dp*covariance/((size(this))*mean1*mean2)  )
            final_eie = final_error/sqrt(2.0_dp*(size(this)-1))
            write(iunit,*) size(that), mean_proje, final_error, final_eie
        enddo
        close(iunit)

    end subroutine print_proje_blocks

    subroutine automatic_reblocking_analysis(this,blocklength,corrlength,filename,tPrint,iValue)
    ! Performs automatic reblocking (Flyvbjerg and Petersen)
    ! plus some crude selection scheme
    ! General routine, does not require global data
    ! iValue indicates the value which is currently being blocked:
    ! 1 = denominantor
    ! 2 = real numerator
    ! 3 = shift
    ! 4 = imaginary numerator
        integer, intent(in) :: iValue
        real(dp), intent(in) :: this(:)
        character(len=*), intent(in) :: filename
        logical, intent(in) :: tPrint
        integer :: length
        integer,intent(in) :: blocklength ! this is the minimum step size (e.g. 2)
        integer,intent(out) :: corrlength ! this is the correlation length in iterations (e.g. 2**N)
        integer :: i
        real(dp) :: mean, error, eie
        real(dp), allocatable :: that(:) ! stores this array
        integer :: blocking_events
        real(dp), allocatable :: mean_array(:), error_array(:), eie_array(:)
        integer :: which_element,iunit
        real(dp) :: final_error
        character(len=*), parameter :: t_r="automatic_reblocking_analysis"

        length=size(this,1)
        allocate(that(length))
        that=this
        blocking_events=int(log(real(size(that),dp))/log(2.0_dp)-1.0_dp)
                ! Yuck! This justs finds the
                ! expected number of reblocking analyses automatically
        allocate(mean_array(blocking_events+1))
        allocate(error_array(blocking_events+1))
        allocate(eie_array(blocking_events+1))

        iunit = 0
        if(tPrint) then
            iunit = get_free_unit()
            open(iunit,file=filename)
            write(iunit,"(A)") "# No.Blocks      Mean       Error        Error_in_Error"
        endif

        call analyze_data(that,mean,error,eie)
        if(tPrint) write(iunit,*) size(that), mean, error, eie
        mean_array(1)=mean
        error_array(1)=error
        eie_array(1)=eie
        do i=1,blocking_events
            call reblock_data(that,blocklength)
            call analyze_data(that,mean,error,eie)
            mean_array(i+1)=mean
            error_array(i+1)=error
            eie_array(i+1)=eie
            if(tPrint) write(iunit,*) size(that), mean, error, eie
        enddo
        call check_reblock_monotonic_inc(error_array,tPrint,iValue)
        call find_max_error(error_array,final_error,which_element)
        corrlength=blocklength**(which_element-1)
        if(tPrint) then
            if(iValue.eq.1) then
                write(6,"(A,I7)") "Number of blocks assumed for calculation of error in projected energy denominator: ", &
                    length/corrlength
            elseif(iValue.eq.2) then
                write(6,"(A,I7)") "Number of blocks assumed for calculation of error in projected energy numerator: ", &
                    length/corrlength
            elseif(iValue.eq.3) then
                write(6,"(A,I7)") "Number of blocks assumed for calculation of error in shift: ",   &
                    length/corrlength
            elseif(iValue.eq.4) then
                write(6,"(A,I7)") "Number of blocks assumed for calculation of error in Im projected energy numerator: ", &
                    length/corrlength
            else
                call stop_all(t_r,"Error in iValue")
            endif
        endif
        if(errordebug.gt.0) then
            write(6,*) "Mean", mean_array(1)
            write(6,*) "Final error", final_error, "number of blocks", length/blocklength**(which_element-1)
            write(6,*) "Correlation length", corrlength
        endif
        if(tPrint) close(iunit)

        deallocate(mean_array,error_array,eie_array)

    end subroutine automatic_reblocking_analysis

    subroutine read_fcimcstats(iShiftVary,tFailRead)
        use SystemData, only: tMolpro, MolproID, tMolproMimic
        integer, intent(inout) :: iShiftVary
        logical, intent(out) :: tFailRead
        character(len=1) :: readline
        character(len=*), parameter :: t_r="read_fcimcstats"
        character(len=24) :: filename
        logical :: exists,tRefToZero
        integer :: eof,comments,i,ierr
        integer :: iunit,WalkersDiffProc
        real(dp) :: doubs,change,Ann,Died,Born, rewalkers, imwalkers
        real(dp) :: shift,rate,reproje,improje,reinstproje,iminstproje
        real(dp) :: AccRat,IterTime,FracFromSing,TotImagTime,HFShift,InstShift
        real(dp) :: denom,renum,imnum,normhf,norm_psi,curr_S2,Avshift,dud,tote
        real(dp) :: curr_S2_init,AbsProjE
        integer(int64) :: iters,validdata,datapoints,totdets
        real(dp), dimension(lenof_sign) :: insthf

        !Open file (FCIMCStats or FCIQMCStats)
        iunit = get_free_unit()
        if(tMolpro .and. .not. tMolproMimic) then
            filename = 'FCIQMCStats_' // adjustl(MolproID)
            inquire(file=filename,exist=exists)
            if (.not. exists) then
                write(iout,*) 'No FCIQMCStats file found for error analysis'
                tFailRead = .true.
                return
            end if
            OPEN(iunit,file=filename,status='old',action='read',position='rewind')
            write(6,"(A)") "Reading back in FCIQMCStats datafile..."
        else
            filename = 'FCIMCStats'
            inquire(file='FCIMCStats',exist=exists)
            if (.not. exists) then
                write(iout, *) 'No FCIMCStats file found for error analysis'
                tFailRead = .true.
                return
            end if
            OPEN(iunit,file='FCIMCStats',status='old',action='read',position='rewind')
            write(6,"(A)") "Reading back in FCIMCStats datafile..."
            if (inum_runs == 2 .and. exists) then
                write(6,"(A)") " This is a DNECI run! But we just analyse the first FCIMCStats!"
            end if
        endif


        comments=0
        datapoints=0
        validdata=0
        tRefToZero = .false.
        do while(.true.)

            read(iunit,"(A)",advance='no',iostat=eof) readline
            if(eof.lt.0) then
                exit
            elseif(eof.gt.0) then
                call stop_all(t_r,"Error reading FCIMCStats file")
            endif

            if(readline.ne.'#') then
                !Valid data on line

                if(lenof_sign > inum_runs) then
                    !complex fcimcstats
                    read(iunit, *, iostat=eof) &
                        iters, &                !1.
                        shift, &                              !2.
                        change, &   !3.
                        rate, &   !4.
                        rewalkers, &                       !5.
                        imwalkers, &                       !6.
                        reproje, &                !7.     real \sum[ nj H0j / n0 ]
                        improje, &                   !8.     Im   \sum[ nj H0j / n0 ]
                        reinstproje, &                 !9.
                        iminstproje, &                    !10.
                        tote, &       ! Tot.ProjE.iter (Re)
                        insthf(1), &                         !11.
                        insthf(lenof_sign), &                         !12.
                        doubs, &                         !13.
                        AccRat, &                               !14.
                        TotDets, &                        !15.
                        IterTime, &                             !16.
                        FracFromSing, &     !17.
                        WalkersDiffProc, &                           !18.
                        TotImagTime, &                               !19.
                        HFShift, &                                   !20.
                        InstShift, &                                 !21.
                        denom     !24     |n0|^2  This is the denominator for both calcs
                else
                    read(iunit, *, iostat=eof) &
                        Iters, &
                        shift, &
                        change, &
                        rate, &
                        rewalkers, &
                        Ann, &
                        Died, &
                        Born, &
                        reproje, &
                        Avshift, &
                        reinstproje, &
                        insthf(1), &
                        doubs, &
                        AccRat, &
                        totdets, &
                        IterTime, &
                        FracFromSing, &
                        WalkersDiffProc, &
                        TotImagTime, &
                        dud, &
                        HFShift, &
                        InstShift, &
                        tote, &
                        denom
                endif

                if(eof.lt.0) then
                    call stop_all(t_r,"Should not be at end of file")
                elseif(eof.gt.0) then
                    ! This is normally due to a difficulty reading NaN or
                    ! Infinity. Assume that this occurs when NoAtRef --> 0.
                    ! Therefore we can safely wipe the stats.
                    if (iters > iShiftVary) then
                        tRefToZero = .true.
                        validdata=0
                        iShiftVary = int(Iters,sizeof_int) + 1
                    endif
                else
                    if(iters.gt.iShiftVary) then
                        if(abs(denom).lt.1.0e-5_dp) then
                            !Denominator gone to zero - wipe the stats
                            tRefToZero = .true.
                            validdata=0
                            iShiftVary = int(Iters,sizeof_int) + 1
                        else
                            validdata=validdata+1
                        endif
                    end if
                endif
                datapoints=datapoints+1

            else
                !Just read it again without the advance=no to move onto the next line
                read(iunit,"(A)",iostat=eof) readline
                if(eof.lt.0) then
                    call stop_all(t_r,"Should not be at end of file")
                elseif(eof.gt.0) then
                    call stop_all(t_r,"Error reading FCIMCStats file")
                endif
                comments=comments+1
            endif
        enddo

        if(tRefToZero) then
            write(6,"(A)") "After shift varies, reference population goes to zero"
            write(6,"(A,I14)") "Blocking will only start after non-zero population, at iteration ",iShiftVary-1
        endif

        write(6,"(A,I12)") "Number of comment lines found in file: ",comments
        write(6,"(A,I12)") "Number of data points found in file: ",datapoints
        write(6,"(A,I12)") "Number of useable data points: ",validdata

        if(validdata.le.1) then
            write(6,"(A)") "No valid datapoints found in file. Exiting error analysis."
            tFailRead = .true.
            return
        else
            tFailRead = .false.
        endif

        !Allocate arrays
        allocate(numerator_data(validdata),stat=ierr)
        if(lenof_sign > inum_runs) then
            allocate(imnumerator_data(validdata),stat=ierr)
        endif
        allocate(pophf_data(validdata),stat=ierr)
        allocate(shift_data(validdata),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Memory allocation error")
        numerator_data(:) = 0.0_dp
        pophf_data(:) = 0.0_dp
        shift_data(:) = 0.0_dp
        rewind(iunit)

        !Now read all the data in again
        i=0
        do while(.true.)

            read(iunit,"(A)",advance='no',iostat=eof) readline
            if(eof.lt.0) then
                exit
            elseif(eof.gt.0) then
                call stop_all(t_r,"Error reading FCIMCStats file")
            endif

            if(readline.ne.'#') then
                !Valid data on line

                if(lenof_sign > inum_runs) then
                    ! complex fcimcstats
                    read(iunit, *, iostat=eof) &
                        iters, &                !1.
                        shift, &                              !2.
                        change, &   !3.
                        rate, &   !4.
                        rewalkers, &                       !5.
                        imwalkers, &                       !6.
                        reproje, &                !7.     real \sum[ nj H0j / n0 ]
                        improje, &                   !8.     Im   \sum[ nj H0j / n0 ]
                        reinstproje, &                 !9.
                        iminstproje, &                    !10.
                        tote,  &     ! Tot.PorjE.iter (Rm)
                        insthf(1), &                         !11.
                        insthf(lenof_sign), &                         !12.
                        doubs, &                         !13.
                        AccRat, &                               !14.
                        TotDets, &                        !15.
                        IterTime, &                             !16.
                        FracFromSing, &     !17.
                        WalkersDiffProc, &                           !18.
                        TotImagTime, &                               !19.
                        HFShift, &                                   !20.
                        InstShift, &                                 !21.
                        denom, &     !24     |n0|^2  This is the denominator for both calcs
                        renum, &   !22.    Re[\sum njH0j] x Re[n0] + Im[\sum njH0j] x Im[n0]   No div by StepsSft
                        imnum, &       !23.    Im[\sum njH0j] x Re[n0] - Re[\sum njH0j] x Im[n0]   since no physicality
                        normhf, &
                        norm_psi, &
                        curr_S2
                else
                    read(iunit, *, iostat=eof) &
                        Iters, &
                        shift, &
                        change, &
                        rate, &
                        rewalkers, &
                        Ann, &
                        Died, &
                        Born, &
                        reproje, &
                        Avshift, &
                        reinstproje, &
                        insthf(1), &
                        doubs, &
                        AccRat, &
                        totdets, &
                        IterTime, &
                        FracFromSing, &
                        WalkersDiffProc, &
                        TotImagTime, &
                        dud, &
                        HFShift, &
                        InstShift, &
                        tote, &
                        denom, &
                        renum, &
                        normhf, &
                        norm_psi, &
                        curr_S2, curr_S2_init, &
                        AbsProjE
                endif

                if(eof.lt.0) then
                    call stop_all(t_r,"Should not be at end of file")
                elseif(eof.gt.0) then
                    ! This is normally due to a difficulty reading NaN or
                    ! Infinity. Assume that this occurs when NoAtRef --> 0.
                    ! Therefore, we should not be trying to store the stats,
                    ! and can ignore the error.
                    if (iters > iShiftVary) &
                        call stop_all (t_r, "This is not a valid data line - &
                                            &should not be trying to store")
                endif
                if(iters.gt.iShiftVary) then
                    i=i+1
                    if(tGivenEquilibrium) then
                        !Count the update cycles we are to ignore
                        if(iters.lt.Shift_Start) then
                            relaxation_time_shift = relaxation_time_shift + 1
                        endif
                        if(iters.lt.ProjE_Start) then
                            relaxation_time_proje = relaxation_time_proje + 1
                        endif
                    endif
                    numerator_data(i) = renum
                    if(lenof_sign > inum_runs) then
                        imnumerator_data(i) = imnum
                    endif
                    pophf_data(i) = denom
                    shift_data(i) = shift
                endif
            else
                !Just read it again without the advance=no to move onto the next line
                read(iunit,"(A)",iostat=eof) readline
                if(eof.lt.0) then
                    call stop_all(t_r,"Should not be at end of file")
                elseif(eof.gt.0) then
                    call stop_all(t_r,"Error reading FCIMCStats file")
                endif
            endif
        enddo

        if(i.ne.validdata) then
            call stop_all(t_r,"Data arrays not filled correctly 1")
        endif

        if(abs(shift_data(validdata)) < 1.0e-10_dp .or. abs(pophf_data(validdata)) < 1.0e-10_dp .or. &
            abs(numerator_data(validdata)) < 1.0e-10_dp) then
            write(6,*) shift_data(validdata), pophf_data(validdata), numerator_data(validdata)
            call stop_all(t_r,"Data arrays not filled correctly 2")
        endif

        close(iunit)

    end subroutine read_fcimcstats

    subroutine resize(this,remove)
    ! Takes a chunk of size remove off the beginning of array
    ! General routine, does not require global data

        real(dp), allocatable :: this(:)
        real(dp), allocatable :: new(:)
        integer :: remove
        integer :: length
        integer :: i
        character(len=*), parameter :: t_r="resize"

        if(.not. allocated(this)) call stop_all(t_r,"Error, array not allocated on entry into resize")
        length=size(this,1)
        allocate(new(length-remove))
        new=0.0_dp
        do i=1,length
            if (i.gt.remove) then
                new(i-remove)=this(i)
            endif
        enddo
        if (abs(new(length-remove)) < 1.0e-10_dp) call stop_all(t_r,"Resize failed")
        deallocate(this)
        allocate(this(length-remove))
        this=new
        deallocate(new)

    end subroutine resize

    subroutine attempt_detect_equil_time(this,equil_time,tFail)
    ! General routine, does not require global data
    ! Tries to remove the equilibration time from the start of the array

        real(dp), allocatable :: this(:)
        integer :: length
        real(dp), allocatable :: tmp(:)
        integer, allocatable :: i_tmp(:)
        logical :: swapped
        logical, intent(out) :: tFail
        real(dp) :: store
        integer :: i,i_store
        integer :: equil_time
        character(len=*), parameter :: t_r='attempt_detect_equil_time'

        equil_time=0
        if(Errordebug.gt.0) write(6,*) "Attempting to detect equilibration time"

        if(.not. allocated(this)) call stop_all(t_r,"Error, array not allocated on entry into " &
            & //"attempt_remove_equil_time")
        length=size(this,1)
        allocate(tmp(length))
        allocate(i_tmp(length))
        tmp=this

        do i=1,length
            i_tmp(i)=i
        enddo

        ! bubble sort - yay
        do
            swapped=.false.
            do i=2,length
                if (tmp(i).lt.tmp(i-1)) then
                    store=tmp(i-1)
                    tmp(i-1)=tmp(i)
                    tmp(i)=store
                    i_store=i_tmp(i-1)
                    i_tmp(i-1)=i_tmp(i)
                    i_tmp(i)=i_store
                    swapped=.true.
                endif
            enddo
            if (.not.swapped) exit
        enddo

!        do i=1,length
            !write(6,*) tmp(i),i_tmp(i) ! this just prints the sorted data
                                        ! useful to check routine does what is expected
!        enddo

        ! A crude way to remove equilibration time automatically
        ! is to remove the first value (labelled 1st in i_tmp) in the list
        ! if it's ranked first or last in the sorted list
        ! Then to repeat this process.
        ! This can be done without resorting the list because now the 2nd element
        ! must lie at the (n-1)th position or the 2nd position (depending on where
        ! the first element was).

        tFail=.false.

        if (i_tmp(1).eq.1) then
            if(Errordebug.gt.0) write(6,*) "Equilib time detected at lowest value"
            do i=1,length
                if (i_tmp(i).ne.i) then
                    if(Errordebug.gt.0) write(6,*) "Equilib time detected to be", i-1, "values"
                    equil_time=i-1
                    exit
                endif
                if (i.eq.length) then
                    if(Errordebug.gt.0) then
                        write(6,*) "ERROR: the whole data file would be removed by automatic equilibrium detection"
                    endif
                    tFail = .true.
                    return
                endif
            enddo
        endif
        if (i_tmp(length).eq.1) then
            if(Errordebug.gt.0) write(6,*) "Equilib time detected at the highest value"
            do i=1,length
                if (i_tmp(length-i+1).ne.i) then ! Yuck, integer adding, prone to error
                    if(Errordebug.gt.0) write(6,*) "Equilib time detected to be", i-1, "values"
                    equil_time=i-1
                    exit
                endif
                if (i.eq.length) then
                    if(Errordebug.gt.0) then
                        write(6,*) "ERROR: the whole data file would be removed by automatic equilibrium detection"
                    endif
                    tFail = .true.
                    return
                endif
            enddo
        endif

    end subroutine attempt_detect_equil_time

    subroutine analyze_data(this,mean,error,eie)
    ! Finds the mean, error and error in error for Flyvbjerg and Petersen
    ! blocking analysis
    ! General routine, does not require global data

        real(dp), intent(in) :: this(:)
        real(dp), intent(out) :: mean,error,eie
        real(dp) :: mean2 ! squared value mean
        integer :: length
        integer :: i
        real(dp) :: s,s2 ! sum x_i and sum x_i^2

        s=0.0_dp
        s2=0.0_dp
        mean=0.0_dp
        mean2=0.0_dp
        error=0.0_dp
        eie=0.0_dp
        length=size(this,1)
        do i=1,length
            s=s+this(i)
            s2=s2+this(i)**2
        enddo
        mean=s/length
        mean2=s2/length
        error=sqrt(mean2-mean**2.0_dp)/sqrt(real(length-1,dp))

        eie=error/sqrt(2.0_dp*(length-1))

    end subroutine analyze_data

    subroutine reblock_data(this,blocklength)
    ! Reduces the data vector for Flyvbjerg and Petersen blocking analysis
    ! by a factor of blocklength (commonly 2, only tested for 2)
    ! General routine, does not require global data

        real(dp), allocatable :: this(:)
        integer :: length,new_length,ind_end,i,j
        real(dp), allocatable :: tmp(:)
        integer,intent(in) :: blocklength
        integer :: ierr
        character(len=*), parameter :: t_r='reblock_data'

        if(.not. allocated(this)) call stop_all(t_r,"Error, array not allocated on entry into "&
            & //"reblock_data")
        length=size(this,1)
        new_length=length/blocklength ! truncates towards zero deliberate, loses data
        !write(6,*) "Length is", length
        !write(6,*) "Will be reduced to", new_length
        !write(6,*) "Loss of data", length-new_length*blocklength
        allocate(tmp(new_length),stat=ierr)
        if (ierr < 0) &
            call stop_all(t_r, 'Bad allocation')
        tmp=0.0_dp
        j=1 ! lazy but foolproof - counting elements
        do i=1,length,blocklength
            ind_end=i+blocklength-1 ! integer addition is disgusting
            if (ind_end.le.length) then
                tmp(j)=average_vector(this(i:ind_end))
            endif
            j=j+1
        enddo
#ifdef __DEBUG
!        if (abs(tmp(new_length)) < 1.0e-10_dp) call stop_all(t_r,"Whole length of new vector not properly used")
#endif
        deallocate(this)
        allocate(this(new_length))
        this=0.0_dp
        this=tmp
        deallocate(tmp)

    end subroutine reblock_data

    function average_vector(this)
    ! Returns average of a vector
    ! General routine, does not require global data

        real(dp) :: this(:)
        integer :: length
        integer :: i
        real(dp) :: s ! sum of array elements
        real(dp) :: average_vector

        s=0
        length=size(this,1)
        do i=1,length
            s=s+this(i)
        enddo
        if (length == 0) then
            average_vector = 0
        else
            average_vector = s / length
        endif

    end function average_vector

    subroutine print_vector(this,filename)
    ! Just prints a vector (useful for debug)
    ! General routine, does not require global data

        real(dp) :: this(:)
        character(len=*), intent(in), optional :: filename
        integer :: length
        integer :: i,iunit

        length=size(this,1)
        iunit = get_free_unit()
        if (present(filename)) open(iunit,file=filename)
        do i=1,length
            if (present(filename)) then
                write(iunit,*) i,this(i)
            else
                write(6,*) i,this(i)
            endif
        enddo
        if (present(filename)) close(iunit)

    end subroutine print_vector

    subroutine check_reblock_monotonic_inc(these_errors,tPrint,iValue)
    ! One of the simplest checks on the errors for F&P blocking analysis
    ! is to look for a monotonic increase in errors
    ! which indicates no tail-off/plateauing.
    ! General routine, does not require global data
    ! iValue indicates the value which is currently being blocked:
    ! 1 = denominantor
    ! 2 = real numerator
    ! 3 = shift
    ! 4 = imaginary numerator

        real(dp), intent(in) :: these_errors(:)
        logical , intent(in) :: tPrint
        integer , intent(in) :: iValue
        integer :: length
        integer :: i
        logical :: monotonic
        character(len=*), parameter :: t_r='check_reblock_monotonic_inc'

        monotonic=.true.
        length=size(these_errors)
        do i=2,length
            if (these_errors(i).lt.these_errors(i-1)) monotonic=.false.
        enddo
        if (monotonic.and.tPrint) then
            if(iValue.eq.1) then
                write(6,"(A)") "WARNING: Error increases monotonically on the blocking graph for " &
                    & //"*denominator of projected energy*"
            elseif(iValue.eq.2) then
                write(6,"(A)") "WARNING: Error increases monotonically on the blocking graph for " &
                    & //"*numerator of projected energy*"
            elseif(iValue.eq.3) then
                write(6,"(A)") "WARNING: Error increases monotonically on the blocking graph for *shift*"
            elseif(iValue.eq.4) then
                write(6,"(A)") "WARNING: Error increases monotonically on the blocking graph for " &
                    & //"*imaginary numerator of projected energy*"
            else
                call stop_all(t_r,"Unknown iValue passed in")
            endif
            write(6,"(A)") "         whilst performing Flyvbjerg and Petersen blocking analysis."
            write(6,"(A)") "         Inspect BLOCKING files carefully. Manual reblocking may be necessary."
        elseif(monotonic.and.(ErrorDebug.gt.0)) then
            if(iValue.eq.1) then
                write(6,"(A)") "WARNING: Error increases monotonically on the blocking graph for " &
                    & //"*denominator of projected energy*"
            elseif(iValue.eq.2) then
                write(6,"(A)") "WARNING: Error increases monotonically on the blocking graph for " &
                    & //"*numerator of projected energy*"
            elseif(iValue.eq.3) then
                write(6,"(A)") "WARNING: Error increases monotonically on the blocking graph for *shift*"
            elseif(iValue.eq.4) then
                write(6,"(A)") "WARNING: Error increases monotonically on the blocking graph for " &
                    & //"*imaginary numerator of projected energy*"
            else
                call stop_all(t_r,"Unknown iValue passed in")
            endif
            write(6,"(A)") "         whilst performing Flyvbjerg and Petersen blocking analysis. If this warning"
            write(6,"(A)") "         appears after equilibration time has been removed, then inspect"
            write(6,"(A)") "         BLOCKING files carefully"
        endif

    end subroutine check_reblock_monotonic_inc

    subroutine find_max_error(these_errors,error,which_element)
    ! One of the simplest ways to choose the error in F+P reblocking
    ! is to choose the largest error that's not the one of the last two points
    ! (which guarentees 8 bits of data)
    ! General routine, does not require global data

        real(dp) :: these_errors(:), error
        integer :: length,which_element,i

        length=size(these_errors,1)
        error=these_errors(1)
        which_element=1
        do i=2,length-2
            if (these_errors(i).gt.error) then
                error=these_errors(i)
                which_element=i
            endif
        enddo

    end subroutine find_max_error

    function calc_covariance(this,that)
    ! Covariance between two arrays
    ! General routine, no global data

        real(dp), intent(in) :: this(:),that(:)
        integer :: length1,length2
        integer :: i
        real(dp) :: sx,sy,sxy ! sums
        real(dp) :: meanx,meany
        real(dp) :: calc_covariance
        character(len=*), parameter :: t_r='calc_covariance'

        length1=size(this)
        length2=size(that)
        if (length1.ne.length2) call stop_all(t_r,"ERROR: Something has gone very wrong indeed...")
        sx=0.0_dp
        sy=0.0_dp
        do i=1,length1
            sx=sx+this(i)
            sy=sy+that(i)
            if(Errordebug.gt.0) write (6,*) this(i), that(i)
        enddo
        meanx=sx/length1
        meany=sy/length1
        sxy=0.0_dp
        do i=1,length1
            sxy=sxy+(this(i)-meanx)*(that(i)-meany)
        enddo
        calc_covariance=sxy/(length1)

    end function calc_covariance


    !Routine to just calculate errors from FCIMCStats file
    subroutine Standalone_Errors()
        use sym_mod, only: getsym
        USE MolproPlugin, only : MolproPluginResult
#ifdef MOLPRO
        use outputResult
        integer :: nv,ityp(1)
#endif
        real(dp) :: mean_ProjE_re,mean_ProjE_im,mean_Shift
        real(dp) :: ProjE_Err_re,ProjE_Err_im,Shift_Err
        logical :: tNoProjEValue,tNoShiftValue
        integer :: iroot,isymh
        TYPE(BasisFn) RefSym
        HElement_t(dp) :: h_tmp
        real(dp) :: Hii,BestEnergy,EnergyDiff
#ifdef MOLPRO
        real(dp) :: get_scalar
        include "common/molen"
#endif

        !Automatic error analysis
        call error_analysis(tSinglePartPhase(1),iBlockingIter(1),mean_ProjE_re,ProjE_Err_re,  &
            mean_ProjE_im,ProjE_Err_im,mean_Shift,Shift_Err,tNoProjEValue,tNoShiftValue, &
            equilshift=iBlockEquilShift,equilproje=iBlockEquilProjE)
        call MPIBCast(ProjectionE)
        call MPIBCast(mean_ProjE_re)
        call MPIBCast(ProjE_Err_re)
        call MPIBCast(mean_ProjE_im)
        call MPIBCast(ProjE_Err_im)
        call MPIBCast(mean_Shift)
        call MPIBCast(Shift_Err)
        call MPIBCast(tNoProjEValue)
        call MPIBCast(tNoShiftValue)

        h_tmp = get_helement (FDet, FDet, 0)
        Hii = real(h_tmp, dp)

        iroot=1
        CALL GetSym(FDet,NEl,G1,NBasisMax,RefSym)
        isymh=int(RefSym%Sym%S,sizeof_int)+1
        write (iout,10101) iroot,isymh
10101   format(//'RESULTS FOR STATE',i2,'.',i1/'====================='/)
        write (iout,'('' Current reference energy'',T52,F19.12)') Hii
        if(tNoProjEValue) then
            write(iout,'('' No projected energy value could be obtained'',T52)')
        else
            write (iout,'('' Projected correlation energy'',T52,F19.12)') mean_ProjE_re
            write (iout,'('' Estimated error in Projected correlation energy'',T52,F19.12)') ProjE_Err_re
        endif
        if(lenof_sign > inum_runs) then
            write (iout,'('' Projected imaginary energy'',T52,F19.12)') mean_ProjE_im
            write (iout,'('' Estimated error in Projected imaginary energy'',T52,F19.12)') ProjE_Err_im
        endif
        if(tNoShiftValue) then
            write(iout,'('' No shift energy value could be obtained'',T52)')
        else
            write (iout,'('' Shift correlation energy'',T52,F19.12)') mean_Shift
            write (iout,'('' Estimated error in shift correlation energy'',T52,F19.12)') shift_err
        endif

        !Do shift and projected energy agree?
        write(iout,"(A)")
        if(tNoProjEValue.and.tNoShiftValue) return
        EnergyDiff = abs(mean_Shift-mean_ProjE_re)
        if(EnergyDiff.le.sqrt(shift_err**2+ProjE_Err_re**2)) then
            write(iout,"(A,F15.8)") " Projected and shift energy estimates agree " &
                & //"within errorbars: EDiff = ",EnergyDiff
        elseif(EnergyDiff.le.sqrt((max(shift_err,ProjE_Err_re)*2)**2+min(shift_err,ProjE_Err_re)**2)) then
            write(iout,"(A,F15.8)") " Projected and shift energy estimates agree to within two sigma " &
                & //"of largest error: EDiff = ",EnergyDiff
        else
            write(iout,"(A,F15.8)") " Projected and shift energy estimates do not agree to " &
                & //"within approximate errorbars: EDiff = ",EnergyDiff
        endif
        if(ProjE_Err_re.lt.shift_err) then
            BestEnergy = mean_ProjE_re + Hii
        else
            BestEnergy = mean_shift + Hii
        endif
        write(iout,"(A)")
        write(iout,"(A,F20.8,A,G15.6)") " Total projected energy ", &
            mean_ProjE_re+Hii," +/- ",ProjE_Err_re
        write(iout,"(A,F20.8,A,G15.6)") " Total shift energy     ", &
            mean_shift+Hii," +/- ",shift_err

#ifdef MOLPRO
        call output_result('FCIQMC','ENERGY',BestEnergy,iroot,isymh)
        if (iroot.eq.1) call clearvar('ENERGY')
        ityp(1)=1
        call setvar('ENERGY',BestEnergy,'AU',ityp,1,nv,iroot)
        do i=10,2,-1
            gesnam(i)=gesnam(i-1)
            energ(i)=energ(i-1)
        enddo
        gesnam(i) = 'FCIQMC'
        energ(i) = get_scalar("ENERGY")
        call output_result('FCIQMC','FCIQMC_ERR',min(ProjE_Err_re,shift_err),iroot,isymh)
        if (iroot.eq.1) call clearvar('FCIQMC_ERR')
        call setvar('FCIQMC_ERR',min(ProjE_Err_re,shift_err),'AU',ityp,1,nv,iroot)
#endif
        CALL MolproPluginResult('ENERGY',[BestEnergy])
        CALL MolproPluginResult('FCIQMC_ERR',[min(ProjE_Err_re,shift_err)])
        write(iout,"(/)")

    end subroutine Standalone_Errors


end module
