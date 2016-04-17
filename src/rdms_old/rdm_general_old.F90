#include "macros.h"

module rdm_general_old

    ! This module contains general routines related to RDM calculation in
    ! FCIQMC. This includes the initialisation and end routines, and also a
    ! few routines (such as calc_rdmbiasfac and store_parent_with_spawned) used
    ! during the main simulation.

    use bit_rep_data, only: NIfTot, NIfDBO
    use SystemData, only: NEl, nBasis
    use util_mod
    use constants

    implicit none

contains

    ! Initialisation routines.

    subroutine InitRDMs_old(nrdms)

        ! This routine initialises any of the arrays needed to calculate the
        ! reduced density matrix. It is used for both the explicit and
        ! stochastic RDMs.

        use FciMCData, only: tSinglePartPhase
        use LoggingData, only: tNo_RDMs_To_Read, tWrite_RDMs_to_read
        use LoggingData, only: tRDMInstEnergy, RDMExcitLevel, tPrint1RDM
        use LoggingData, only: tDiagRDM, tReadRDMs, tPopsfile, tDumpForcesInfo, tDipoles
        use Parallel_neci, only: iProcIndex
        use rdm_data, only: tOpenShell, tCalc_RDMEnergy
        use rdm_data_old, only: rdms, rdm_write_unit_old, rdm_estimates_old
        use RotateOrbsData, only: NoOrbs, SpatOrbs
        use util_mod, only: get_free_unit, LogMemAlloc

        integer, intent(in) :: nrdms

        integer :: ierr, i, iproc, rdm_size_1, rdm_size_2
        integer :: MemoryAlloc, MemoryAlloc_Root
        character(len=*), parameter :: t_r = 'InitRDMs_old'

        if (.not. allocated(rdms)) allocate(rdms(nrdms))
        if (.not. allocated(rdm_estimates_old)) allocate(rdm_estimates_old(nrdms))

        ! There are two different size arrays to allocate, as given by these
        ! array sizes (the RDMs are there sizes squared).
        rdm_size_1 = (SpatOrbs*(SpatOrbs-1))/2
        rdm_size_2 = (SpatOrbs*(SpatOrbs+1))/2

        MemoryAlloc = 0
        ! Memory allocated in bytes.
        MemoryAlloc_Root = 0

        ! First for the storage of the actual 1- or 2-RDM.
        if (RDMExcitLevel .eq. 1) then

            do i = 1, nrdms
                allocate(rdms(i)%matrix(NoOrbs, NoOrbs), stat=ierr)
                if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating 1-RDM array,')
                call LogMemAlloc('rdms(i)%matrix', NoOrbs**2, 8, t_r, rdms(i)%matrix_tag, ierr)
                rdms(i)%matrix(:,:) = 0.0_dp

                MemoryAlloc = MemoryAlloc + ( NoOrbs * NoOrbs * 8 )
                MemoryAlloc_Root = MemoryAlloc_Root + ( NoOrbs * NoOrbs * 8 )
            end do
        else
            ! If we're calculating the 2-RDM, the 1-RDM does not need to be
            ! calculated as well because all its info is in the 2-RDM anyway.

            ! The 2-RDM of the type alpha alpha alpha alpha (= beta beta beta beta).
            ! These *do not* include any 2-RDM(i,j,a,b) terms where i=j or a=b (if
            ! they're the same spin this can't happen).

            do i = 1, nrdms
                if (tRDMInstEnergy) then
                    ! We will be filling up rdms(i)%aaaa_inst as we go along, which need
                    ! to be allocated per core.
                    ! When calculating the energy, these will be summed over cores
                    ! using an _inplace type command.
                    ! To calculate the full energy of the RDM (i.e. over full accum.
                    ! period), we need to allocate rdms(i)%aaaa_full on the head nodes.

                    allocate(rdms(i)%aaaa_inst(rdm_size_1, rdm_size_1), stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r,'Problem allocating aaaa_inst RDM array,')
                    call LogMemAlloc('rdms(i)%aaaa_inst', (rdm_size_1**2), 8, t_r, rdms(i)%aaaa_instTag, ierr)
                    rdms(i)%aaaa_inst(:,:) = 0.0_dp

                    ! The 2-RDM of the type alpha beta beta alpha (= beta alpha alpha beta).
                    ! These also *do not* also include 2-RDM(i,j,a,b) terms where i=j or a=b
                    ! (these are the same as the abab elements).
                    allocate(rdms(i)%abba_inst(rdm_size_1, rdm_size_1), stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating abba_inst RDM array,')
                    call LogMemAlloc('rdms(i)%abba_inst', (rdm_size_1**2), 8, t_r, rdms(i)%abba_instTag, ierr)
                    rdms(i)%abba_inst(:,:) = 0.0_dp

                    MemoryAlloc = MemoryAlloc + (rdm_size_1**2)*2*8
                    MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_1**2)*2*8

                    ! The 2-RDM of the type alpha beta alpha beta ( = beta alpha beta alpha).
                    ! These *do* include 2-RDM(i,j,a,b) terms where i=j or a=b, if they're
                    ! different spin this is possible - hence the slightly different size to
                    ! the aaaa array.
                    allocate(rdms(i)%abab_inst(rdm_size_2, rdm_size_2), stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating abab_inst RDM array,')
                    call LogMemAlloc('rdms(i)%abab_inst', (rdm_size_2**2), 8, t_r, rdms(i)%abab_instTag, ierr)
                    rdms(i)%abab_inst(:,:) = 0.0_dp

                    MemoryAlloc = MemoryAlloc + (rdm_size_2**2)*8
                    MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_2**2)*8 

                    ! Extra arrays for open shell systems.
                    if (tOpenShell) then
                        allocate(rdms(i)%bbbb_inst(rdm_size_1, rdm_size_1), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating bbbb_inst RDM array,')
                        call LogMemAlloc('rdms(i)%bbbb_inst', (rdm_size_1**2), 8, t_r, rdms(i)%bbbb_instTag, ierr)
                        rdms(i)%bbbb_inst(:,:) = 0.0_dp

                        allocate(rdms(i)%baab_inst(rdm_size_1, rdm_size_1), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating baab_inst RDM array,')
                        call LogMemAlloc('rdms(i)%baab_inst', (rdm_size_1**2), 8, t_r, rdms(i)%baab_instTag, ierr)
                        rdms(i)%baab_inst(:,:) = 0.0_dp
                        MemoryAlloc = MemoryAlloc + (rdm_size_1**2)*2*8
                        MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_1**2)*2*8

                        allocate(rdms(i)%baba_inst(rdm_size_2, rdm_size_2), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating baba_inst RDM array,')
                        call LogMemAlloc('rdms(i)%baba_inst', (rdm_size_2**2), 8, t_r, rdms(i)%baba_instTag, ierr)
                        rdms(i)%baba_inst(:,:) = 0.0_dp
                        MemoryAlloc = MemoryAlloc + (rdm_size_2**2)*8
                        MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_2**2)*8
                    end if

                    if (iProcIndex .eq. 0) then
                        allocate(rdms(i)%aaaa_full(rdm_size_1, rdm_size_1), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating aaaa_full RDM array,')
                        call LogMemAlloc('rdms(i)%aaaa_full', (rdm_size_1**2), 8, t_r, rdms(i)%aaaa_fullTag, ierr)
                        rdms(i)%aaaa_full(:,:) = 0.0_dp

                        allocate(rdms(i)%abab_full(rdm_size_2, rdm_size_2), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating abab_full RDM array,')
                        call LogMemAlloc('rdms(i)%abab_full', (rdm_size_2**2), 8, t_r, rdms(i)%abab_fullTag, ierr)
                        rdms(i)%abab_full(:,:) = 0.0_dp

                        allocate(rdms(i)%abba_full(rdm_size_1, rdm_size_1), stat = ierr)
                        if (ierr .ne. 0) call stop_all(t_r,'Problem allocating abba_full RDM array,')
                        call LogMemAlloc('rdms(i)%abba_full', (rdm_size_1**2), 8, t_r, rdms(i)%abba_fullTag, ierr)
                        rdms(i)%abba_full(:,:) = 0.0_dp

                        MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_1**2)*2*8
                        MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_2**2)*8
                        MemoryAlloc = MemoryAlloc + (rdm_size_1**2)*2*8
                        MemoryAlloc = MemoryAlloc + (rdm_size_2**2)*8

                        if (tOpenShell) then
                            allocate(rdms(i)%bbbb_full(rdm_size_1, rdm_size_1), stat=ierr)
                            if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating bbbb_full RDM array,')
                            call LogMemAlloc('rdms(i)%bbbb_full', (rdm_size_1**2), 8, t_r,rdms(i)%bbbb_fullTag,ierr)
                            rdms(i)%bbbb_full(:,:) = 0.0_dp

                            allocate(rdms(i)%baba_full(rdm_size_2,rdm_size_2),stat=ierr)
                            if (ierr .ne. 0) call stop_all(t_r,'Problem allocating baba_full RDM array,')
                            call LogMemAlloc('rdms(i)%baba_full', (rdm_size_2**2), 8,t_r, rdms(i)%baba_fullTag, ierr)
                            rdms(i)%baba_full(:,:) = 0.0_dp

                            allocate(rdms(i)%baab_full(rdm_size_1, rdm_size_1), stat=ierr)
                            if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating baab_full array,')
                            call LogMemAlloc('rdms(i)%baab_full', (rdm_size_1**2), 8,t_r, rdms(i)%baab_fullTag, ierr)
                            rdms(i)%baab_full(:,:) = 0.0_dp

                            MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_1**2)*2*8
                            MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_2**2)*8
                            MemoryAlloc = MemoryAlloc + (rdm_size_1**2)*2*8
                            MemoryAlloc = MemoryAlloc + (rdm_size_2**2)*8

                        end if

                    end if
                    
                    rdms(i)%aaaa => rdms(i)%aaaa_inst
                    rdms(i)%abba => rdms(i)%abba_inst
                    rdms(i)%abab => rdms(i)%abab_inst

                    if (tOpenShell) then
                        rdms(i)%bbbb => rdms(i)%bbbb_inst
                        rdms(i)%baab => rdms(i)%baab_inst
                        rdms(i)%baba => rdms(i)%baba_inst
                    end if

                else
                    ! We're not calculating an instantaneous RDM energy.
                    ! Put RDM contributions directly into 'full' arrays, which are
                    ! now allocated every core
                    allocate(rdms(i)%aaaa_full(rdm_size_1,rdm_size_1), stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating aaaa_full RDM array,')
                    call LogMemAlloc('rdms(i)%aaaa_full', (rdm_size_1**2), 8, t_r,rdms(i)%aaaa_fullTag,ierr)
                    rdms(i)%aaaa_full(:,:) = 0.0_dp

                    ! The 2-RDM of the type alpha beta beta alpha ( = beta alpha alpha beta).
                    ! These also *do not* also include 2-RDM(i,j,a,b) terms where i=j or a=b
                    ! (these are the same as the abab elements).
                    allocate(rdms(i)%abba_full(rdm_size_1,rdm_size_1),stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r,'Problem allocating abba_full RDM array,')
                    call LogMemAlloc('rdms(i)%abba_full', (rdm_size_1**2), 8, t_r, rdms(i)%abba_fullTag, ierr)
                    rdms(i)%abba_full(:,:) = 0.0_dp

                    MemoryAlloc = MemoryAlloc + (rdm_size_1**2)*2*8
                    MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_1**2)*2*8

                    ! The 2-RDM of the type alpha beta alpha beta ( = beta alpha beta alpha).
                    ! These *do* include 2-RDM(i,j,a,b) terms where i=j or a=b, if they're
                    ! different spin this is possible - hence the slightly different size
                    ! to the aaaa array.
                    allocate(rdms(i)%abab_full(rdm_size_2, rdm_size_2), stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating abab_full RDM array,')
                    call LogMemAlloc('rdms(i)%abab_full', (rdm_size_2**2), 8, t_r, rdms(i)%abab_fullTag, ierr)
                    rdms(i)%abab_full(:,:) = 0.0_dp

                    MemoryAlloc = MemoryAlloc + (rdm_size_2**2)*8
                    MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_2**2)*8
                    
                    rdms(i)%aaaa => rdms(i)%aaaa_full
                    rdms(i)%abba => rdms(i)%abba_full
                    rdms(i)%abab => rdms(i)%abab_full

                    if (tOpenShell) then
                        allocate(rdms(i)%bbbb_full(rdm_size_1, rdm_size_1), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating bbbb_full RDM array,')
                        call LogMemAlloc('rdms(i)%bbbb_full', (rdm_size_1**2), 8, t_r, rdms(i)%bbbb_fullTag, ierr)
                        rdms(i)%bbbb_full(:,:) = 0.0_dp

                        allocate(rdms(i)%baab_full(rdm_size_1,rdm_size_1),stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating baab_full array,')
                        call LogMemAlloc('rdms(i)%baab_full', (rdm_size_1**2), 8, t_r, rdms(i)%baab_fullTag, ierr)
                        rdms(i)%baab_full(:,:) = 0.0_dp
                        MemoryAlloc = MemoryAlloc + (rdm_size_1**2)*2*8
                        MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_1**2)*2*8

                        allocate(rdms(i)%baba_full(rdm_size_2, rdm_size_2), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating baba_full RDM array,')
                        call LogMemAlloc('rdms(i)%baba_full', (rdm_size_2**2), 8, t_r, rdms(i)%baba_fullTag, ierr)
                        rdms(i)%baba_full(:,:) = 0.0_dp
                        MemoryAlloc = MemoryAlloc + (rdm_size_2**2)*8
                        MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_2**2)*8
          
                        rdms(i)%bbbb => rdms(i)%bbbb_full
                        rdms(i)%baab => rdms(i)%baab_full
                        rdms(i)%baba => rdms(i)%baba_full
                    end if
                end if ! Not instantaneous

                if (iProcindex .eq. 0) then
                    if (tDiagRDM .or. tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) then
                        ! Still need to allocate 1-RDM to get nat orb occupation numbers.
                        allocate(rdms(i)%matrix(NoOrbs, NoOrbs), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating 1-RDM array,')
                        call LogMemAlloc('rdms(i)%matrix', NoOrbs**2,8, t_r, rdms(i)%matrix_tag, ierr)
                        rdms(i)%matrix(:,:) = 0.0_dp
                        MemoryAlloc_Root = MemoryAlloc_Root + ( NoOrbs * NoOrbs * 8 )
                    end if
                end if
            end do ! Looping over all RDMs.

        end if

        if ((RDMExcitLevel .eq. 1) .or. tDiagRDM .or. tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) then

            if ((iProcIndex .eq. 0) .and. tDiagRDM) then
                do i = 1, nrdms
                    allocate(rdms(i)%Evalues(NoOrbs), stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating Evalues array,')
                    call LogMemAlloc('rdms(i)%Evalues', NoOrbs, 8, t_r, rdms(i)%EvaluesTag, ierr)
                    rdms(i)%Evalues(:) = 0.0_dp

                    allocate(rdms(i)%Rho_ii(NoOrbs), stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating Rho_ii array,')
                    call LogMemAlloc('Rho_ii', NoOrbs, 8, t_r, rdms(i)%Rho_iiTag, ierr)
                    rdms(i)%Rho_ii(:) = 0.0_dp
                end do
            end if

        end if

        ! Open file to keep track of RDM estimates. 
        if ((iProcIndex .eq. 0) .and. tCalc_RDMEnergy) then
            rdm_write_unit_old = get_free_unit()
            open(rdm_write_unit_old, file='RDMEstimates_old', status='unknown', position='append')

            write(rdm_write_unit_old, '("#", 4X, "Iteration")', advance='no')

            do i = 1, nrdms
                write(rdm_write_unit_old, '(4x,"Energy numerator",1x,i2)', advance='no') i
                write(rdm_write_unit_old, '(4x,"Spin^2 numerator",1x,i2)', advance='no') i
                write(rdm_write_unit_old, '(7x,"Normalisation",1x,i2)', advance='no') i
            end do

            write(rdm_write_unit_old,'()')
        end if

        ! Reads in the RDMs from a previous calculation, sets the accumulating
        ! normalisations, writes out the starting energy.
        if (tReadRDMs) then
            if (nrdms > 1) call stop_all(t_r, "Reading in multiple RDMs is not yet supported.")
            if (tSinglePartPhase(1) .or. tSinglePartPhase(inum_runs)) then
                write(6,'(A)') 'WARNING - Asking to read in the RDMs, but not varying shift from &
                                & the beginning of the calculation.'
                write(6,'(A)') 'Ignoring the request to read in the RDMs and starting again.'
                tReadRDMs = .false.
            else
                call Read_In_RDMs_old(rdms(1), rdm_estimates_old(1))
            end if
        end if

        ! By default, if we're writing out a popsfile (and doing an RDM
        ! calculation), we also write out the unnormalised RDMs that can be
        ! read in when restarting a calculation. If the NORDMSTOREAD option
        ! is on, these wont be printed.  
        if (tPopsfile .and. (.not. tno_RDMs_to_read)) twrite_RDMs_to_read = .true.

    end subroutine InitRDMs_old

    subroutine Read_In_RDMs_old(rdm, est_old)

        ! Reads in the arrays to restart the RDM calculation (and continue
        ! accumulating). These arrays are not normalised, so the trace is
        ! also calculated. The energy is then calculated (if required) from
        ! the RDMs read in only.

        use LoggingData, only: IterRDMonFly
        use LoggingData, only: RDMExcitLevel
        use Parallel_neci, only: iProcIndex
        use rdm_data, only: tOpenShell, tCalc_RDMEnergy
        use rdm_data_old, only: rdm_t, rdm_estimates_old_t, rdm_estimates_old
        use rdm_estimators_old, only: rdm_output_wrapper_old, write_rdm_estimates_old
        use RotateOrbsData, only: SymLabelListInv_rot
        use util_mod, only: get_free_unit

        type(rdm_t), intent(inout) :: rdm
        type(rdm_estimates_old_t), intent(inout) :: est_old

        logical :: exists_one
        logical :: exists_aaaa, exists_abab, exists_abba
        logical :: exists_bbbb, exists_baba, exists_baab
        integer :: RDM_unit, FileEnd
        integer :: i, j, a, b, Ind1, Ind2
        real(dp) :: Temp_RDM_Element
        character(len=*), parameter :: t_r = 'Read_In_RDMs_old'

        if (iProcIndex .eq. 0) then 

            if (RDMExcitLevel .eq. 1) then

                write(6,'(1X,"Reading in the 1-RDM")')

                ! The OneRDM will have been printed exactly as is. Without
                ! having been made hermitian, without being normalised, and in
                ! spatial orbitals if tStoreSpinOrbs is false.

                inquire(file='OneRDM_POPS', exist=exists_one)
                if (exists_one) then
                    RDM_unit = get_free_unit()
                    open(RDM_unit, file='OneRDM_POPS', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i, j, Temp_RDM_Element
                        if (FileEnd .gt. 0) call stop_all(t_r, "Error reading OneRDM_POPS")
                        if (FileEnd .lt. 0) exit

                        rdm%matrix(SymLabelListInv_rot(i), SymLabelListInv_rot(j)) = Temp_RDM_Element
                    end do
                    close(RDM_unit)
                else
                    call stop_all(t_r, "Attempting to read in the OneRDM, but the OneRDM_POPS file does not exist.")
                end if                                    

            else

                write(6,'(1X,"Reading in the 2-RDMs")')

                ! The TwoRDMs will have been printed exactly as they were.
                ! Without having been made hermitian, without being
                ! normalised, and in spatial orbitals. 

                ! Only read in the 2-RDMs (the 1-RDM becomes redundant).
                inquire(file='TwoRDM_POPS_aaaa', exist=exists_aaaa)
                inquire(file='TwoRDM_POPS_abab', exist=exists_abab)
                inquire(file='TwoRDM_POPS_abba', exist=exists_abba)

                if (tOpenShell)then
                    inquire(file='TwoRDM_POPS_bbbb', exist=exists_bbbb)
                    inquire(file='TwoRDM_POPS_baba', exist=exists_baba)
                    inquire(file='TwoRDM_POPS_baab', exist=exists_baab)
                end if

                if (exists_aaaa .and. exists_abab .and. exists_abba) then
                    ! All TOREAD RDM files are present - read in.
                    RDM_unit = get_free_unit()
                    open(RDM_unit, file='TwoRDM_POPS_aaaa', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all(t_r, "Error reading TwoRDM_POPS_aaaa")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        rdm%aaaa_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                    open(RDM_unit, file='TwoRDM_POPS_abab', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit, iostat=FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all(t_r, "Error reading TwoRDM_POPS_abab")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-1) * j ) / 2 ) + i
                        Ind2 = ( ( (b-1) * b ) / 2 ) + a
                        rdm%abab_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                    open(RDM_unit, file='TwoRDM_POPS_abba', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all(t_r, "Error reading TwoRDM_POPS_abba")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        rdm%abba_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                else
                    write(6,*) 'exists_aaaa', exists_aaaa
                    write(6,*) 'exists_abab', exists_abab
                    write(6,*) 'exists_abba', exists_abba
                    call neci_flush(6)
                    call stop_all(t_r,"Attempting to read in the RDMs, &
                                    &but at least one of the TwoRDM_a***_TOREAD files are missing.")
                end if

                if (tOpenShell .and. exists_bbbb .and. exists_baba .and. exists_baab) then
                    ! All TOREAD RDM files for open shell RDMs are present - read in.
                    RDM_unit = get_free_unit()
                    open(RDM_unit, file='TwoRDM_POPS_bbbb', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit, iostat=FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all(t_r, "Error reading TwoRDM_POPS_bbbb")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        rdm%bbbb_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                    open(RDM_unit, file='TwoRDM_POPS_baba', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit, iostat = FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd.gt.0) call stop_all(t_r, "Error reading TwoRDM_POPS_baba")
                        if (FileEnd.lt.0) exit

                        Ind1 = ( ( (j-1) * j ) / 2 ) + i
                        Ind2 = ( ( (b-1) * b ) / 2 ) + a
                        rdm%baba_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                    open(RDM_unit,file='TwoRDM_POPS_baab',status='old',form='unformatted')
                    do while (.true.)
                        read(RDM_unit, iostat = FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all(t_r, "Error reading TwoRDM_POPS_baab")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        rdm%baab_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                else if (tOpenShell) then
                    write(6,*) 'exists_bbbb', exists_bbbb
                    write(6,*) 'exists_baba', exists_baba
                    write(6,*) 'exists_baab', exists_baab
                    call neci_flush(6)
                    call stop_all(t_r, "Attempting to read in the RDMs, &
                                  &but at least one of the TwoRDM_b***_TOREAD files are missing.")
                end if

            end if
        end if

        ! Calculate the energy for the matrices read in (if we're calculating more
        ! than the 1-RDM).
        if (tCalc_RDMEnergy) then
            call rdm_output_wrapper_old(rdm, 1, est_old)
            if (iProcIndex == 0) call write_rdm_estimates_old(rdm_estimates_old)
        end if

        ! Continue calculating the RDMs from the first iteration when the popsfiles
        ! (and RDMs) are read in. This overwrites the iteration number put in the input.
        IterRDMonFly = 0

    end subroutine Read_In_RDMs_old

    ! Routines called at the end of a simulation.

    subroutine FinaliseRDMs_old(rdms, rdm_estimates_old)

        ! This routine performs some finalisation, including summing each of
        ! the individual matrices from each processor, and calling the
        ! diagonalisation routines if we want to get the occupation numbers.

#ifdef _MOLCAS_
        USE EN2MOLCAS, only : NECI_E
#endif
        use FciMCData, only: tFinalRDMEnergy
        use LoggingData, only: tBrokenSymNOs, occ_numb_diff, RDMExcitLevel, tExplicitAllRDM
        use LoggingData, only: tPrint1RDM, tDiagRDM, tDumpForcesInfo, tDipoles
        use Parallel_neci, only: iProcIndex, MPIBarrier, MPIBCast, MPISumAll
        use rdm_data, only: rdm_estimates_t, tRotatedNos
        use rdm_data, only: rdm_main, one_rdms, two_rdm_spawn, two_rdm_recv, tOpenShell
        use rdm_data_old, only: rdm_t, rdm_estimates_old_t
        use rdm_estimators_old, only: Calc_Lagrangian_from_RDM, convert_mats_Molpforces
        use rdm_estimators_old, only: rdm_output_wrapper_old, CalcDipoles
        use rdm_estimators, only: rdm_output_wrapper, write_rdm_estimates
        use rdm_general, only: Finalise_1e_RDM, calc_1e_norms
        use rdm_nat_orbs, only: find_nat_orb_occ_numbers, BrokenSymNo
        use rdm_parallel, only: calc_rdm_trace, calc_1rdms_from_2rdms
        use rdm_parallel, only: create_spinfree_2rdm, calc_1rdms_from_spinfree_2rdms

        type(rdm_t), intent(inout) :: rdms(:)
        type(rdm_estimates_old_t), intent(inout) :: rdm_estimates_old(:)

        integer :: i, ierr
        real(dp) :: Norm_1RDM, Trace_1RDM, SumN_Rho_ii

        ! Combine the 1- or 2-RDM from all processors, etc.

        do i = 1, size(rdms)

            if (RDMExcitLevel .eq. 1) then
                call Finalise_1e_RDM(rdms(i)%matrix, rdms(i)%Rho_ii, i, Norm_1RDM, .true.)
            else
                ! We always want to calculate one final RDM energy, whether or not we're 
                ! calculating the energy throughout the calculation.
                ! Unless of course, only the 1-RDM is being calculated.

                ! Calculate the energy one last time - and write out everything we need.
                tFinalRDMEnergy = .true.

                ! 1-RDM is constructed here (in calc_1RDM_and_1RDM_energy).
                call rdm_output_wrapper_old(rdms(i), i, rdm_estimates_old(i))

                if (tPrint1RDM) then
                    call Finalise_1e_RDM(rdms(i)%matrix, rdms(i)%Rho_ii, i, Norm_1RDM, .true.)
                else if (tDiagRDM .and. (iProcIndex .eq. 0)) then
                    call calc_1e_norms(rdms(i)%matrix, rdms(i)%Rho_ii, Trace_1RDM, Norm_1RDM, SumN_Rho_ii)
                    write(6,'(/,1X,"SUM OF 1-RDM(i,i) FOR THE N LOWEST ENERGY HF ORBITALS:",1X,F20.13)') SumN_Rho_ii
                end if

                if (tDumpForcesInfo) then
                    if (.not. tPrint1RDM) call Finalise_1e_RDM(rdms(i)%matrix, rdms(i)%Rho_ii, i, Norm_1RDM, .true.)
                    call Calc_Lagrangian_from_RDM(rdms(i), Norm_1RDM, rdm_estimates_old(i)%Norm_2RDM)
                    call convert_mats_Molpforces(rdms(i), Norm_1RDM, rdm_estimates_old(i)%Norm_2RDM)
                end if

            end if

            call MPIBarrier(ierr)

            ! Call the routines from NatOrbs that diagonalise the one electron
            ! reduced density matrix.
            tRotatedNOs = .false. ! Needed for BrokenSymNo routine
            if (tDiagRDM) call find_nat_orb_occ_numbers(rdms(i), i)

            ! This is where we would likely call any further calculations of
            ! forces, etc.
            if (tDipoles) then
                if (.not. tPrint1RDM) call Finalise_1e_RDM(rdms(i)%matrix, rdms(i)%Rho_ii, i, Norm_1RDM, .true.)
                call CalcDipoles(rdms(i), Norm_1RDM)
            end if

            ! After all the NO calculations are finished we'd like to do another
            ! rotation to obtain symmetry-broken natural orbitals
            if (tBrokenSymNOs) then
                call BrokenSymNO(rdms(i)%Evalues, occ_numb_diff)
            end if

        end do

    end subroutine FinaliseRDMs_old

    subroutine DeallocateRDMs_old()

        ! This routine just deallocates the arrays allocated in InitRDMs.
        ! If the NECI calculation softexits before the RDMs start to fill,
        ! this is all that is called at the end.

        use FciMCData, only: Spawned_Parents, Spawned_Parents_Index
        use FciMCData, only: Spawned_ParentsTag, Spawned_Parents_IndexTag
        use LoggingData, only: RDMExcitLevel, tExplicitAllRDM
        use rdm_data, only: Sing_ExcDjs, Doub_ExcDjs
        use rdm_data, only: Sing_ExcDjs2, Doub_ExcDjs2, Sing_ExcDjsTag, Doub_ExcDjsTag
        use rdm_data, only: Sing_ExcDjs2Tag, Doub_ExcDjs2Tag
        use rdm_data, only: Sing_InitExcSlots, Doub_InitExcSlots, Sing_ExcList, Doub_ExcList
        use rdm_data_old, only: rdms
        use RotateOrbsData, only: SymLabelCounts2_rot,SymLabelList2_rot, SymLabelListInv_rot
        use RotateOrbsData, only: SymLabelCounts2_rotTag, SymLabelList2_rotTag
        use RotateOrbsData, only: SymLabelListInv_rotTag
        use RotateOrbsMod, only: FourIndInts, FourIndIntsTag
        use util_mod, only: LogMemDealloc

        integer :: i
        character(len=*), parameter :: t_r = 'DeallocateRDMs_old'

        do i = 1, size(rdms)
            if (allocated(rdms(i)%matrix)) then
                deallocate(rdms(i)%matrix)
                call LogMemDeAlloc(t_r,rdms(i)%matrix_tag)
            end if

            if (allocated(rdms(i)%Evalues)) then
                deallocate(rdms(i)%Evalues)
                call LogMemDeAlloc(t_r,rdms(i)%EvaluesTag)
            end if

            if (allocated(rdms(i)%Rho_ii)) then
                deallocate(rdms(i)%Rho_ii)
                call LogMemDeAlloc(t_r,rdms(i)%Rho_iiTag)
            end if

            if (associated(rdms(i)%aaaa_inst)) then
                deallocate(rdms(i)%aaaa_inst)
                nullify(rdms(i)%aaaa_inst)
                call LogMemDeAlloc(t_r,rdms(i)%aaaa_instTag)
                nullify(rdms(i)%aaaa_inst)
            end if

            if (associated(rdms(i)%abab_inst)) then
                deallocate(rdms(i)%abab_inst)
                nullify(rdms(i)%abab_inst)
                call LogMemDeAlloc(t_r,rdms(i)%abab_instTag)
                nullify(rdms(i)%abab_inst)
            end if

            if (associated(rdms(i)%abba_inst)) then
                deallocate(rdms(i)%abba_inst)
                nullify(rdms(i)%abba_inst)
                call LogMemDeAlloc(t_r,rdms(i)%abba_instTag)
                nullify(rdms(i)%abba_inst)
            end if

            if (associated(rdms(i)%bbbb_inst)) then
                deallocate(rdms(i)%bbbb_inst)
                nullify(rdms(i)%bbbb_inst)
                call LogMemDeAlloc(t_r,rdms(i)%bbbb_instTag)
                nullify(rdms(i)%bbbb_inst)
            end if

            if (associated(rdms(i)%baba_inst)) then
                deallocate(rdms(i)%baba_inst)
                nullify(rdms(i)%baba_inst)
                call LogMemDeAlloc(t_r,rdms(i)%baba_instTag)
                nullify(rdms(i)%baba_inst)
            end if

            if (associated(rdms(i)%baab_inst)) then
                deallocate(rdms(i)%baab_inst)
                nullify(rdms(i)%baab_inst)
                call LogMemDeAlloc(t_r,rdms(i)%baab_instTag)
                nullify(rdms(i)%baab_inst)
            end if

            if (associated(rdms(i)%aaaa_full)) then
                deallocate(rdms(i)%aaaa_full)
                nullify(rdms(i)%aaaa_full)
                call LogMemDeAlloc(t_r,rdms(i)%aaaa_fullTag)
                nullify(rdms(i)%aaaa_full)
            end if

            if (associated(rdms(i)%abab_full)) then
                deallocate(rdms(i)%abab_full)
                nullify(rdms(i)%abab_full)
                call LogMemDeAlloc(t_r,rdms(i)%abab_fullTag)
                nullify(rdms(i)%abab_full)
            end if

            if (associated(rdms(i)%abba_full)) then
                deallocate(rdms(i)%abba_full)
                nullify(rdms(i)%abab_full)
                call LogMemDeAlloc(t_r,rdms(i)%abba_fullTag)
                nullify(rdms(i)%abba_full)
            end if

            if (associated(rdms(i)%bbbb_full)) then
                deallocate(rdms(i)%bbbb_full)
                nullify(rdms(i)%bbbb_full)
                call LogMemDeAlloc(t_r,rdms(i)%bbbb_fullTag)
                nullify(rdms(i)%bbbb_full)
            end if

            if (associated(rdms(i)%baba_full)) then
                deallocate(rdms(i)%baba_full)
                nullify(rdms(i)%baba_full)
                call LogMemDeAlloc(t_r,rdms(i)%baba_fullTag)
                nullify(rdms(i)%baba_full)
            end if

            if (associated(rdms(i)%baab_full)) then
                deallocate(rdms(i)%baab_full)
                nullify(rdms(i)%baba_full)
                call LogMemDeAlloc(t_r,rdms(i)%baab_fullTag)
                nullify(rdms(i)%baab_full)
            end if

            if (allocated(rdms(i)%sym_list_no)) deallocate(rdms(i)%sym_list_no)
            if (allocated(rdms(i)%sym_list_inv_no)) deallocate(rdms(i)%sym_list_inv_no)

            if (associated(rdms(i)%aaaa)) nullify(rdms(i)%aaaa)
            if (associated(rdms(i)%abab)) nullify(rdms(i)%abab)
            if (associated(rdms(i)%abba)) nullify(rdms(i)%abba)
            if (associated(rdms(i)%bbbb)) nullify(rdms(i)%bbbb)
            if (associated(rdms(i)%baba)) nullify(rdms(i)%baba)
            if (associated(rdms(i)%baab)) nullify(rdms(i)%baab)
        end do

        deallocate(rdms)

    end subroutine DeallocateRDMs_old

end module rdm_general_old
