#include "macros.h"

module hist

    use DeterminantData, only: get_lexicographic, calculated_ms
    use MemoryManager
    use SystemData, only: tHistSpinDist, ilut_spindist, nbasis, nel, LMS, &
                          hist_spin_dist_iter, nI_spindist, LMS, tHPHF, &
                          tOddS_HPHF, G1
    use DetBitOps, only: count_open_orbs, EncodeBitDet, spatial_bit_det, &
                         DetBitEq, count_open_orbs, TestClosedShellDet, &
                         CalcOpenOrbs, IsAllowedHPHF
    use hash , only : DetermineDetNode                     
    use CalcData, only: tFCIMC, tTruncInitiator
    use DetCalcData, only: FCIDetIndex, det
    use FciMCData, only: tFlippedSign, TotWalkers, CurrentDets, iter, &
                         all_norm_psi_squared
    use util_mod, only: choose, get_free_unit, binary_search
    use HPHFRandExcitMod, only: FindExcitBitDetSym
    use hphf_integrals, only: hphf_sign
    use constants, only: n_int, bits_n_int, size_n_int, lenof_sign
    use bit_rep_data, only: NIfTot, NIfD, extract_sign
    use bit_reps, only: encode_sign, extract_bit_rep, NOffSgn, &
                        decode_bit_det, flag_is_initiator, test_flag
    use parallel_neci
    use csf, only: get_num_csfs, csf_coeff, csf_get_yamas, write_yama, &
                   extract_dorder
!    use AnnihilationMod, only: DetermineDetNode
    use hist_data
    use timing_neci
    use Determinants, only: write_det

    implicit none

contains

    subroutine init_hist_spin_dist ()

        ! Initialise the spin histogramming. This precalculates nopen and
        ! ncsf, so that dynamic local arrays can be used without using
        ! heap allocation.

        integer :: nopen, ncsf, S
        integer(n_int) :: ilut_tmp(0:NIfTot)

        ! Encode spin distribution
        if (.not. allocated(ilut_spindist)) &
            allocate(ilut_spindist(0:NIfTot))
        call EncodeBitDet(nI_spindist(1:nel), ilut_tmp)
        ilut_spindist = spatial_bit_det(ilut_tmp)

        ! How many unpaired electrons
        nopen = count_open_orbs (ilut_spindist)

        ! How many csfs are there?
        ncsf = 0
        do S = LMS, nopen, 2
            ncsf = ncsf + get_num_csfs (nopen, S)
        enddo

        call init_hist_spin_dist_local (nopen, ncsf)

    end subroutine


    subroutine init_hist_spin_dist_local (nopen, ncsf)
    
        integer, intent(in) :: ncsf, nopen
        integer :: open_orbs(nel), dorder(nel), ierr, orb, i, S, tmp
        integer :: nfound, nup, ndets, fd
        integer(n_int), pointer :: p_ilut(:)
        integer :: yamas(ncsf, nopen), abo
        character(*), parameter :: this_routine = 'init_hist_spin_dist'
        character(20) :: fmt_str

        ! Initialise the spin distribution histograms
        if (.not. allocated(ilut_spindist)) &
            call stop_all (this_routine, 'Spin distribution must be specified&
                                & and initialised to use spin distribution &
                                &histogramming.')

        ! How many spin alpha electrons are there?
        nup = (nopen + LMS) / 2

        ! How many dets with this config?
        ndets = int(choose(nopen, nup))

        ! Allocate and initialise storage
        allocate(hist_spin_dist(0:NIfTot, ndets), stat=ierr)
        LogAlloc(ierr, 'hist_spin_dist', (niftot+1)*ndets, size_n_int, tag_spindist)
        hist_spin_dist = 0

        ! Obtain a list of the open spatial orbitals.
        ! N.B. This is a ZERO BASED list
        nfound = 0
        do i = 1, nbasis-1, 2
            if (IsOcc(ilut_spindist, i) .and. &
             IsNotOcc(ilut_spindist, i+1)) then
                nfound = nfound + 1
                open_orbs(nfound) = i
                if (nfound == nopen) exit
            endif
        enddo

        ! Obtain a list of the CSFs required, and open a file to output their
        ! coefficients
        call csf_get_yamas (nopen, LMS, yamas, ncsf)
        fd = get_free_unit ()
        open(unit=fd, file='det-csf-coeffs', status='replace')
        write(fd,'("# Print ilut(0:NIfD), then coeffs of each avail csf.")')
        write(fd,'("# NIfD = ",i2)') NIfD
        write(fd,'("# 2*Ms = ",i2)') LMS
        write(fd,'("# Total num csfs = ", i4)') ncsf
        nfound = 0
        do S = LMS, nopen, 2
            tmp = get_num_csfs (nopen, S)
            write(fd,'("# S = ",i2,", Num csfs = ",i4)') S, tmp
            call csf_get_yamas (nopen, S, yamas(nfound+1:nfound+tmp,:), tmp)
            do i = 1, tmp
                write(fd,'("# ",i4,": ")', advance='no') nfound + i
                call write_yama(fd, yamas(nfound + i,:), .true.)
            enddo
            nfound = nfound + tmp
        enddo
        if (nfound /= ncsf) &
            call stop_all (this_routine,"Incorrect number of CSFs found")

        ! Allocate list of csf coefficients
        if (allocated(hist_csf_coeffs)) deallocate(hist_csf_coeffs)
        allocate(hist_csf_coeffs(ndets, ncsf), stat=ierr)
        LogAlloc(ierr, 'hist_csf_coeffs', ndets*ncsf, sizeof_dp, tag_histcsfs)

        ! Generate all possible determinants with the given spatial structure
        ! and value of Ms
        write(6,*) 'Initialising Spin distribution histogramming'
        write(fmt_str,'("(",i2,"i20,",i4,"f20.15)")') NIfD+1, ncsf
        call neci_flush(6)
        nfound = 0
        dorder(1) = -1
        call get_lexicographic (dorder, nopen, nup)
        do while (dorder(1) /= -1)

            nfound = nfound + 1
            call ilut_nifd_pointer_assign(p_ilut, &
                                          hist_spin_dist(0:NIfD, nfound))
            p_ilut = ilut_spindist(0:NIfD)

            ! Obtain bit representation by shifting things
            do i = 1, nopen
                if (dorder(i) == 0) then
                    orb = get_alpha(open_orbs(i))
                else
                    orb = get_beta(open_orbs(i))
                endif
                set_orb(p_ilut, orb)
                abo=ab_pair(orb)
                clr_orb(p_ilut, abo)
            enddo
            write(6,'(i5,": ")', advance='no')
            call writebitdet(6, p_ilut, .false.)

            ! Look at the csf coeffs.
            ! This could be a forall, but PGI doesn't like it.
            do i = 1, ncsf
                hist_csf_coeffs(nfound, i) = &
                    csf_coeff(yamas(i,:), dorder, nopen) 
            enddo
            write(fd,fmt_str) p_ilut(0:NIfD), hist_csf_coeffs(nfound, :)

            ! Continue the loop
            call get_lexicographic (dorder, nopen, nup)
        enddo
        p_ilut => null()

        if (nfound /= ndets) &
            call stop_all (this_routine, &
                           "Generated incorrect number of determinants")

    end subroutine

    subroutine write_clear_hist_spin_dist (iter, nsteps)

        ! Output a spin histogram to the file: spin-hist-%iter
        ! The values contain the values in hist_spin_dist / niter
        !
        ! --> Only need to print out every nsteps, and get averaged value
        !     over those steps.

        integer, intent(in) :: iter, nsteps
        integer :: fd, i, j, ierr
        integer(n_int) :: all_hist(0:NIfTot, ubound(hist_spin_dist, 2))
        real(dp) :: csf_contrib(lenof_sign, ubound(hist_csf_coeffs, 2))
        real(dp) :: sgn(lenof_sign)
        character(22) :: fname, iterstr

        ! Collate the data on the head node.
        call MPIReduce (hist_spin_dist, MPI_SUM, all_hist)

        if (iProcIndex == Root) then

            ! Open the file for writing
            fd = get_free_unit ()
            write(iterstr, '(i12)') iter
            write(fname, '("spin-hist-",a)') trim(adjustl(iterstr))
            open(unit=fd, file=trim(fname), status='replace')

            ! Output file header
            write(fd, '("# 1. Determinant(niftot+1)\t2. sign(lenof_sign)")')

            csf_contrib = 0
            do i = 1, ubound(all_hist, 2)

                ! Extract sign
                call extract_sign (all_hist(:,i), sgn)

                ! Output to file
                write(fd, *) all_hist(0:NIfD, i), dble(sgn)/nsteps
                
                ! Add csf contribs
                do j = 1, ubound(csf_contrib, 2)
                    csf_contrib(:,j) = csf_contrib(:,j) + &
                        hist_csf_coeffs(i, j) * sgn
                enddo

            enddo

            ! Close the file
            close(fd)

            ! Open csf histogram file
            write(fname, '("csf-hist-",a)') trim(adjustl(iterstr))
            open(unit=fd, file=trim(fname), status='replace')
            
            ! Output csf components to hist file
            do i = 1, ubound(csf_contrib, 2)
                write(fd,*) csf_contrib(:,i)/nsteps
            enddo

            ! Close the file
            close(fd)

        endif


        ! Zero all of the elements
        sgn = 0
        do i = 1, ubound(hist_spin_dist, 2)
            call encode_sign(hist_spin_dist(:,i), sgn)
        enddo

        call MPIBarrier (ierr)

    end subroutine


    subroutine clean_hist_spin_dist ()

        character(*), parameter :: this_routine = 'clean_hist_spin_dist'

        if (allocated(ilut_spindist)) &
            deallocate(ilut_spindist)

        if (allocated(nI_spindist)) &
            deallocate(nI_spindist)

        if (allocated(hist_spin_dist)) then
            deallocate(hist_spin_dist)
            LogDealloc(tag_spindist)
        endif

    end subroutine


    subroutine ilut_nifd_pointer_assign (ptr, ilut)

        ! This is a useful helper function to get around some irritating
        ! behaviour (namely that pointer array slices are indexed to
        ! begin at 0, not 1 --> wrong for iluts).

        integer(n_int), intent(in), target :: ilut(0:NIfTot)
        integer(n_int), intent(out), pointer :: ptr(:)

        ptr => ilut

    end subroutine

    subroutine test_add_hist_spin_dist_det (ilut, Sign)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer(n_int) :: ilut_tmp(0:NIfTot)
        real(dp), intent(in) :: Sign(lenof_sign)
        real(dp) :: sgn_old(lenof_sign), new_sgn(lenof_sign)
        integer :: pos, i

        ! Should we add this ilut to the histogram?
        ilut_tmp = spatial_bit_det (ilut)
        
        if (DetBitEq(ilut_tmp, ilut_spindist)) then
            pos = binary_search (hist_spin_dist, ilut, NIfD+1)
            if (pos < 0) then
                call writebitdet(6, ilut, .false.)
                write(6,*) ilut
                write(6,*) '================'
                do i=1,ubound(hist_spin_dist, 2)
                    write(6,*) hist_spin_dist(:,i)
                enddo
                call stop_all ("test_add_hist_spin_dist_det", &
                               "Determinant not found in spin histogram list")
            endif
            call extract_sign (hist_spin_dist(:,pos), sgn_old)
            new_sgn = sgn_old + sign
            call encode_sign(hist_spin_dist(:,pos), new_sgn)
        endif

    end subroutine

    subroutine add_hist_spawn (ilut, sign, ExcitLevel, dProbFin)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: ExcitLevel
        real(dp), intent(in) :: dProbFin
        real(dp), intent(in) :: sign(lenof_sign)

        integer :: PartInd, open_orbs
        integer(n_int) :: ilut_sym(0:NIfTot)
        real(dp) :: delta(lenof_sign)
        logical :: tSuccess
        character(*), parameter :: t_r = 'add_hist_spawn'
        
        if (ExcitLevel == nel) then
            call BinSearchParts2 (ilut, HistMinInd(ExcitLevel), det, PartInd,&
                                  tSuccess)
            ! CCMC doesn't sum particle contributions in order, so we must 
            ! search the whole space again!
            if (tFCIMC) HistMinInd(ExcitLevel) = PartInd
        elseif (ExcitLevel == 0) then
            PartInd = 1
            tSuccess = .true.
        else
            call BinSearchParts2 (ilut, HistMinInd(ExcitLevel), &
                                  FCIDetIndex(ExcitLevel+1)-1, PartInd, &
                                  tSuccess)
            ! CCMC doesn't sum particle contributions in order, so we must
            ! search the whole space again!
            if (tFCIMC) HistMinInd(ExcitLevel) = PartInd
        endif

        if (tSuccess) then
            delta = sign / dProbFin

            if (tHPHF) then
                call FindExcitBitDetSym (ilut, ilut_sym)
                if (.not. DetBitEq(ilut, ilut_sym)) delta = delta / sqrt(2.0_dp)
            endif

            if (tFlippedSign) delta = -delta

            Histogram(:, PartInd) = Histogram(:, PartInd) + delta
            if (tHistSpawn) &
                InstHist(:, PartInd) = InstHist(:, PartInd) + delta

            ! In HPHF we also include the spin-coupled determinant, which will
            ! have the same amplitude as the original determinant, unless it
            ! is antisymmetric.
            if (tHPHF .and. .not. DetBitEQ(ilut, ilut_sym)) then
                if (ExcitLevel == nel) then
                    call BinSearchParts2 (ilut_sym, FCIDetIndex(ExcitLevel), &
                                          Det, PartInd, tSuccess)
                elseif (ExcitLevel == 0) then
                    PartInd = 1
                    tSuccess = .true.
                else
                    call BinSearchParts2 (ilut_sym, FCIDetIndex(ExcitLevel), &
                                          FCIDetIndex(ExcitLevel+1)-1, &
                                          PartInd, tSuccess)
                endif
                if (tSuccess) then
                    delta = (sign / sqrt(2.0_dp)) / dProbFin

                    call CalcOpenOrbs(ilut_sym, open_orbs)
                    if ((mod(open_orbs, 2) == 1) .neqv. tOddS_HPHF) &
                        delta = -delta

                    if (tFlippedSign .neqv. tOddS_HPHF) delta = -delta

                    Histogram(:, PartInd) = Histogram(:, PartInd) + delta
                    if (tHistSpawn) &
                        InstHist(:, PartInd) = InstHist(:, PartInd) + delta
                endif
            endif
        else
            call writebitdet (6, ilut, .true.)
            write(6,*) '***', ilut
            write(6,*) '***', ExcitLevel, HistMinInd(ExcitLevel), Det
            call stop_all (t_r, "Cannot find corresponding FCI determinant &
                                &when histogramming")
        endif

    end subroutine

    subroutine find_hist_coeff_explicit (ilut, ExcitLevel, PartInd, tSuccess)

        implicit none
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: ExcitLevel
        integer, intent(out) :: PartInd
        logical , intent(out) :: tSuccess

        integer :: open_orbs
        integer(n_int) :: ilut_sym(0:NIfTot)
        real(dp) :: delta(lenof_sign)
        character(*), parameter :: t_r = 'add_hist_spawn'

        tSuccess = .false.
        if (ExcitLevel == nel) then
            call BinSearchParts2 (ilut,  FCIDetIndex(ExcitLevel), det, PartInd,&
                                  tSuccess)
        elseif (ExcitLevel == 0) then
            PartInd = 1
            tSuccess = .true.
        else
            call BinSearchParts2 (ilut, FCIDetIndex(ExcitLevel), &
                                  FCIDetIndex(ExcitLevel+1)-1, PartInd, &
                                  tSuccess)
        endif

    end subroutine

    subroutine add_hist_energies (ilut, Sign, HDiag, ExcitLevel)

        ! This will histogram the energies of the particles, rather than the
        ! determinants themselves.

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp), dimension(lenof_sign), intent(in) :: Sign
        real(dp), intent(in) :: HDiag
        integer, intent(in) :: ExcitLevel
        integer :: bin
        character(*), parameter :: t_r = "add_hist_energies"
        
        bin = int(HDiag / BinRange) + 1
        if (bin > iNoBins) &
            call stop_all (t_r, "Histogramming energies higher than the &
                         &arrays can cope with. Increase iNoBins or BinRange")
        HistogramEnergy(bin) = HistogramEnergy(bin) + sum(abs(Sign))

    end subroutine

    subroutine project_spin_csfs ()

        type nopen_type
            integer, allocatable :: offsets (:)
            integer, allocatable :: nyama (:)
            integer, allocatable :: yamas(:,:)
            real(dp), allocatable :: coeff (:)
        end type
        type(nopen_type), target :: y_storage(LMS:nel)

        integer :: nopen, ntot, S, ncsf, off, flg, j, k, epos
        integer :: dorder(nel), nI(nel)
        real(dp) :: coeff, S_coeffs(LMS:nel), norm, S2, S22
        real(dp) :: sgn(lenof_sign)
        real(dp) :: AllNode_S_coeffs(LMS:nel)

        ! TODO: This bit could be done just once, couldn't it...
        !       We could store all the intermediates.
        do nopen = LMS, nel, 2
            ! How many possibililities are there?
            allocate(y_storage(nopen)%offsets(LMS:nopen), &
                     y_storage(nopen)%nyama(LMS:nopen))

            ! Find how many csfs there are for each spin value, and obtain
            ! their offsets into the array.
            ntot = 0
            do S = LMS, nopen, 2
                ncsf = get_num_csfs (nopen, S)
                y_storage(nopen)%offsets(S) = ntot + 1
                y_storage(nopen)%nyama(S) = ncsf
                ntot = ntot + ncsf
            enddo

            ! Allocate storage space for the Yamanouchi Symbols
            ! TODO: We really need to swap these indices
            allocate(y_storage(nopen)%yamas(ntot, nopen), &
                     y_storage(nopen)%coeff(ntot))
            y_storage(nopen)%yamas = 0
            y_storage(nopen)%coeff = 0

            ! Populate the Yamanouchi symbol arrays.
            do S = LMS, nopen, 2
                ncsf = y_storage(nopen)%nyama(S)
                off = y_storage(nopen)%offsets(S)
                if (ncsf > 0 .and. nopen > 0) then
                    call csf_get_yamas (nopen, S, &
                                    y_storage(nopen)%yamas(off:off+ncsf-1,:),&
                                    ncsf)
                endif
            enddo
        enddo

        do j = 1, int(TotWalkers,sizeof_int)

            ! Extract the current walker
            call extract_bit_rep (CurrentDets(:,j), nI, sgn, flg)
            call extract_dorder (nI, dorder, nopen)

            ! TODO: Be careful when we adjust the order of indices
            do k = 1, ubound(y_storage(nopen)%yamas, 1)
                coeff = csf_coeff (y_storage(nopen)%yamas(k,:), dorder, nopen)
                coeff = coeff * sgn(1)
                y_storage(nopen)%coeff(k) = y_storage(nopen)%coeff(k) + coeff
            enddo
            
        enddo

        ! Sum all of the S components --> get values.
        S_coeffs = 0
        do nopen = LMS, nel, 2
            j = 1
            do S = LMS, nopen, 2
                epos = j + y_storage(nopen)%nyama(S) - 1
                do while (j <= epos)
                    S_coeffs(S) = S_coeffs(S) &
                                + (y_storage(nopen)%coeff(j)**2)
                    j = j + 1
                enddo
            enddo
        enddo

        call MPISum(S_coeffs, AllNode_S_coeffs)
        S_coeffs=AllNode_S_coeffs

        if (iProcIndex == Root) then
            norm = sum(S_coeffs)
            S_coeffs = sqrt(S_coeffs / norm)

            S2 = 0
            do S = LMS, nel, 2
                S2 = S2 + real(S * (S + 2) * S_coeffs(S),dp) / 4
            enddo
            
            write(6,*) 'Scoeffs', iter, S_coeffs
            write(6,*) 'Scsf2', S2
            write(6,*) 'norm compare', norm, all_norm_psi_squared
        endif

        ! Deallocate stuff
        do nopen = LMS, nel, 2
            if (allocated(y_storage(nopen)%offsets)) &
                deallocate(y_storage(nopen)%offsets)
            if (allocated(y_storage(nopen)%yamas)) &
                deallocate(y_storage(nopen)%yamas)
            if (allocated(y_storage(nopen)%nyama)) &
                deallocate(y_storage(nopen)%nyama)
            if (allocated(y_storage(nopen)%coeff)) &
                deallocate(y_storage(nopen)%coeff)
        enddo


    end subroutine


    subroutine project_spins ()

        ! Sum the squares of all of the components.

        real(dp) :: spin_cpts (0:nel), norm
        real(dp) :: AllNode_spin_cpts (0:nel)
        integer :: j, nup, S, flg, ncsf, nopen
        integer :: nI(nel), dorder(nel)
        real(dp) :: sgn(lenof_sign)

        ! Initially, have no component on any of the spins
        spin_cpts = 0

        do j = 1, int(TotWalkers,sizeof_int)

            ! Extract the current walker
            call extract_bit_rep (CurrentDets(:,j), nI, sgn, flg)
            call extract_dorder (nI, dorder, nopen)
            nup = (nopen + LMS) / 2

            do S = LMS, nopen

                ! How many csfs are there?
                ncsf = get_num_csfs (nopen, S)
                call sum_cpts_S (spin_cpts, dorder(1:nopen), ncsf, nopen, &
                                 sgn, S)
            enddo

        enddo

        call MPISum(spin_cpts,AllNode_spin_cpts)
        spin_cpts=AllNode_spin_cpts

        if (iProcIndex == Root) then
            norm = sum(spin_cpts)
            spin_cpts = sqrt(spin_cpts / norm)

            write(6,'(a,i8,13f15.10)') 'spn: ', iter, spin_cpts
        endif


    end subroutine
    
    subroutine sum_cpts_S (spin_cpts, dorder, ncsf, nopen, sgn, S)

        real(dp), intent(inout) :: spin_cpts(0:nel)
        integer, intent(in) :: ncsf, nopen, S
        integer, intent(in) :: dorder(1:nopen)
        real(dp), intent(in) :: sgn(lenof_sign)
        integer :: yamas(ncsf, nopen), y
        real(dp) :: coeff

        call csf_get_yamas (nopen, S, yamas, ncsf)

        ! Iterate through the yamas, and collect coeffs.
        ! TODO: Reorder (yamas(ncsf, nopen) --> yamas(nopen, ncsf)
        do y = 1, ncsf

            coeff = csf_coeff (yamas(y,:), dorder, nopen)
            coeff = (coeff * sgn(1)) ** 2
            spin_cpts (S) = spin_cpts(S) + coeff

        enddo
    
    end subroutine

    function calc_s_squared (only_init) result (ssq)

        ! Calculate the instantaneous value of S^2 for the walkers stored
        ! in CurrentDets
        !
        ! --> This could be generalised to an arbitrary list of iluts. We 
        !     would also then need to calculate the value of psi_squared

        real(dp) :: ssq, tmp
        integer :: i
        logical, intent(in) :: only_init
        type(timer), save :: s2_timer

        ! TODO: Deal with HPHF. Should be fairly easy.
!        write(6,*) 'totwalkers', totwalkers
!        write(6,*) 'dets', currentdets(noffsgn, 1:totwalkers)
        s2_timer%timer_name = 'S2 local'
        call set_timer (s2_timer)

        ssq = 0
        do i = 1, int(TotWalkers,sizeof_int)
            if ((test_flag(CurrentDets(:,i), flag_is_initiator(1)) .or. &
                 test_flag(CurrentDets(:,i), flag_is_initiator(lenof_sign)))&
                 .and. .not. TestClosedShellDet(CurrentDets(:,i))) then
                ssq = ssq + ssquared_contrib (CurrentDets(:,i), only_init)
            end if
        enddo

        ! Sum over all processors and normalise
        call MPISum (ssq, tmp)
        ssq = tmp

        if (all_norm_psi_squared == 0) then
            ssq = 0.0_dp
        else

            ssq = ssq / all_norm_psi_squared

            ! TODO: n.b. This is a hack. LMS appears to contain -2*Ms of the
            !            system I am somewhat astounded I haven't noticed this
            !            before...
            ssq = ssq + real(calculated_ms * (calculated_ms + 2), dp) / 4
        end if

        call halt_timer (s2_timer)

    end function

    function calc_s_squared_multi () result (ssq)
    
        integer :: max_linked, max_per_proc, max_spawned
        real(dp) :: ssq
        type(timer), save :: s2_timer

        s2_timer%timer_name = 'S^2'
        call set_timer (s2_timer)

        max_linked = int(choose(nel, (nel + LMS) / 2))
        max_per_proc = 2 * (max_linked / nProcessors) + 1
        max_spawned = max_per_proc * nProcessors

        ssq = calc_s_squared_multi_worker (max_per_proc, max_spawned)

        call halt_timer (s2_timer)

    end function

    function calc_s_squared_multi_worker (max_per_proc, max_spawned) &
             result(ssq)

        integer :: i, j, k, orb2, orbtmp, pos, ierr
        integer :: nI(nel), nJ(nel), proc
        integer(n_int), pointer :: detcurr(:)
        integer(n_int) :: splus(0:NIfTot), sminus(0:NIfTot)
        logical :: running, any_running
        real(dp) :: ssq, Allssq
        integer :: max_per_proc, max_spawned
        real(dp) :: sgn1(lenof_sign), sgn2(lenof_sign), tmp


        ! Could we pre-initialise all of these data structures
        !integer :: max_linked = int(choose(nel, (nel + LMS)/2))
        integer(n_int) :: det_list(0:NIfTot, max_spawned)
        integer(n_int) :: recv_dets(0:NIfTot, max_spawned)
        integer :: proc_pos (nProcessors), proc_pos_init(nProcessors)
        integer :: send_count(nProcessors), recv_count(nProcessors)
        integer(MPIArg) :: send_data(nProcessors), recv_data(nProcessors)
        integer(MPIArg) :: send_off(nProcessors), recv_off(nProcessors)

        running = (TotWalkers > 0)
        any_running = .true.
        j = 1
        ssq = 0
        forall (i=1:nProcessors) proc_pos_init(i) = (i-1)*max_per_proc + 1

        do while (any_running)

            ! Clear transmission lists
            proc_pos = proc_pos_init
            
            if (running) then

                ! Generate items, add to list (and use the sgn of initial
                ! walker, so we send it to the target processor)
                call ilut_nifd_pointer_assign(detcurr, CurrentDets(0:NIfTot,j))
                call decode_bit_det (nI, detcurr)

                do i = 1, nel
                    orbtmp = get_alpha(nI(i))
                    if (is_beta(nI(i)) &
                        .and. IsNotOcc(detcurr, orbtmp)) then
                        splus = detcurr
                        clr_orb(splus, nI(i))
                        set_orb(splus, orbtmp)

                        do k = 1, nel
                            orb2 = nI(k)
                            if (k == i) orb2 = get_alpha(orb2)

                            orbtmp = get_beta(orb2)
                            if (is_alpha(orb2) &
                                .and. IsNotOcc(splus, orbtmp)) then
                                sminus = splus
                                clr_orb(sminus, orb2)
                                set_orb(sminus, orbtmp)

                                ! Store this det (n.b. contains original sgn)
                                call decode_bit_det(nJ, sminus)
                                proc = DetermineDetNode(nJ, 0) + 1
                                det_list(:,proc_pos(proc)) = sminus
                                proc_pos(proc) = proc_pos(proc) + 1
                            endif
                        enddo
                    endif
                enddo

                ! Walk through the list. Stop when we get to the end.
                j = j + 1
                if (j > TotWalkers) running = .false.

            endif
            
            ! How many elements are there in each list?
            send_count = proc_pos - proc_pos_init
            if (any(send_count > max_per_proc)) &
                send_count = max_per_proc + 1

            call MPIAlltoAll (send_count, 1, recv_count, 1, ierr)

            send_off = int((proc_pos_init - 1) * (NIfTot + 1),MPIArg)
            recv_off(1) = 0
            do i = 2, nProcessors
                recv_off(i) = recv_off(i - 1) + int(recv_count(i - 1),MPIArg)
            enddo
            recv_off = recv_off * int(NIfTot + 1,MPIArg)
            send_data = int(send_count * (NIfTot + 1),MPIArg)
            recv_data = int(recv_count * (NIfTot + 1),MPIArg)

            call MPIAlltoAllv (det_list, send_data, send_off, &
                               recv_dets, recv_data, recv_off, ierr)

            ! Find the det in list, and sum in its term.
            do i = 1, sum(recv_count)

                ! The sign of the source particle
                call extract_sign (recv_dets(:,i), sgn1)

                ! And the generated, connected particle
                pos = binary_search(CurrentDets(:,1:TotWalkers), &
                                    recv_dets(:,i), NIfD+1)
                if (pos > 0) then
                    call extract_sign (CurrentDets(:,pos), sgn2)
                    ssq = ssq + (sgn1(1) * sgn2(1))
                endif

            enddo

            ! Is there anything left to do on any process?
            call MPIAllReduce (running, MPI_LOR, any_running)

        enddo

        call MPISum (ssq, tmp)
        ssq = tmp

        if (all_norm_psi_squared == 0) then
            ssq = 0.0_dp
        else
            ssq = ssq / all_norm_psi_squared

            ! TODO: n.b. This is a hack. LMS appears to contain -2Ms of the
            !            system. I am somewhat astounded I haven't noticed
            !            this before...
            ssq = ssq + real(calculated_ms * (calculated_ms + 2), dp) / 4
        end if

    end function


    function calc_s_squared_star (only_init) result (ssq)

        real(dp) :: ssq
        integer, parameter :: max_per_proc = 1000
        integer(n_int) :: recv_dets(0:NIfTot,max_per_proc)
        integer :: proc_dets, start_pos, nsend, i, p
        integer :: bcast_tmp(2)
        real(dp) :: sgn_tmp(lenof_sign)
        type(timer), save :: s2_timer, s2_timer_init
        real(dp) :: ssq_sum, psi_squared
        real(dp) :: All_ssq_sum, All_psi_squared
        logical, intent(in) :: only_init


        s2_timer%timer_name = 'S^2 star'
        s2_timer_init%timer_name = 'S^2 init'
        if (only_init) then
            call set_timer(s2_timer_init)
        else
            call set_timer (s2_timer)
        endif

        ssq_sum = 0
        psi_squared = 0
        do p = 0, nProcessors-1

            ! How many dets are on processor p
            proc_dets = int(TotWalkers,sizeof_int)
            call MPIBcast (proc_dets, iProcIndex == p)

            ! Send the dets around bit by bit
            start_pos = 1
            do while(start_pos <= proc_dets)

                if (tTruncInitiator .and. only_init) then
                    ! Loop over walkers and only add initiators to bcast list
                    nsend = 0
                    if (p == iProcIndex) then
                        do i = start_pos, int(TotWalkers,sizeof_int)
                            ! Break up the list into correctly sized chunks
                            if (nsend == max_per_proc) exit

                            if (test_flag(CurrentDets(:,i), &
                                          flag_is_initiator(1)) .or. &
                                test_flag(CurrentDets(:,i), &
                                          flag_is_initiator(lenof_sign))) then
                                nsend = nsend + 1
                                recv_dets(:,nsend) = CurrentDets(:,i)

                                ! If we are using initiators only, keep track
                                ! of the overall magnitude of the init-only
                                ! wavefunction.
                                call extract_sign (CurrentDets(:,i), sgn_tmp)
                                psi_squared = psi_squared + sum(sgn_tmp ** 2)
                            endif
                        enddo
                        start_pos = i
                    endif
                    bcast_tmp = (/nsend, start_pos/)
                    call MPIBcast(bcast_tmp, p == iProcIndex)
                    nsend = bcast_tmp(1)
                    start_pos = bcast_tmp(2)

                else
                    ! How many dets to send in this batch?
                    if (start_pos + max_per_proc - 1 > proc_dets) then
                        nsend = proc_dets - start_pos + 1
                    else
                        nsend = max_per_proc
                    endif

                    ! Update list of received dets
                    if (p == iProcIndex) recv_dets(:,1:nsend) = &
                                    CurrentDets(:,start_pos:start_pos+nsend-1)

                    ! Increment position
                    start_pos = start_pos + nsend

                endif
                
                ! Broadcast to all processors
                call MPIBcast (recv_dets(:,1:nsend), iProcIndex == p)

                ! All processors loop over these dets, and calculate their
                ! contribution to S^2
                do i = 1, nsend
                    if (.not. TestClosedShellDet(recv_dets(:,i))) then
                        ! If nopen == 2, and tHPHF, then this can be
                        ! optimised further
                        ssq_sum = ssq_sum + ssquared_contrib (recv_dets(:,i),&
                                                              only_init)
                    endif
                enddo

            enddo

        enddo ! Loop over processors


        ! Sum all of the s squared terms
        if (tTruncInitiator .and. only_init) then
            call MPISum(psi_squared, All_psi_squared)
            psi_squared=All_psi_squared
        else
            psi_squared = all_norm_psi_squared
        end if
        call MPISum(ssq_sum, All_ssq_sum)
        ssq_sum=All_ssq_sum

        if (psi_squared == 0) then
            ssq = 0.0_dp
        else
            ssq = real(ssq_sum,dp) / psi_squared
         
            ! TODO: n.b. This is a hack. LMS appears to contain -2Ms of the
            !            system. I am somewhat astounded I haven't noticed
            !            this before...
            ssq = ssq + real(calculated_ms * (calculated_ms + 2), dp) / 4

            if (only_init) then
                call halt_timer(s2_timer_init)
            else
                call halt_timer (s2_timer)
            endif
        end if

    end function

    function ssquared_contrib (ilut, only_init_) result(ssq)

        ! Calculate the contribution to s-squared from the determinant
        ! provided (from walkers on this processor).
        !
        ! This applies the operator S-S+, returning the result:
        !
        !      <Psi(iProcIndex) | S-S+ | D_i>

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        logical, intent(in), optional :: only_init_
        integer(n_int) :: splus(0:NIfD), sminus(0:NIfD)
        integer(n_int) :: ilut_srch(0:NIfD), ilut_sym(0:NIfD)
        real(dp) :: sgn(lenof_sign), sgn2(lenof_sign), sgn_hphf
        integer :: flg, nI(nel), j, k, orb2, pos, orb_tmp
        integer(int64) :: ssq
        logical :: only_init, inc

        if (present(only_init_)) then
            only_init = only_init_
        else
            only_init = .false.
        end if

        ! Extract details of determinant
        call extract_bit_rep (ilut, nI, sgn, flg)

        ssq = 0
        do j = 1, nel
            if (is_beta(nI(j)) &
                .and. IsNotOcc(ilut, get_alpha(nI(j)))) then
                splus = ilut(0:NIfD)
                orb_tmp = get_alpha(nI(j))
                clr_orb(splus, nI(j))
                set_orb(splus, orb_tmp)

                do k = 1, nel
                    orb2 = nI(k)
                    if (k == j) orb2 = get_alpha(orb2)

                    orb_tmp = get_beta(orb2)
                    if (is_alpha(orb2) &
                        .and. IsNotOcc(splus, orb_tmp)) then
                        sminus = splus
                        clr_orb(sminus, orb2)
                        set_orb(sminus, orb_tmp)

                        ! Adjust for the sign of the paired det in HPHF.
                        sgn_hphf = 1.0
                        if (tHPHF) then
                            if (IsAllowedHPHF(sminus, ilut_sym)) then
                                ilut_srch = sminus
                            else
                                ilut_srch = ilut_sym
                                sgn_hphf = hphf_sign (ilut_srch)
                            endif
                        else
                            ilut_srch = sminus
                        endif

                        ! --> sminus is an allowed result of applying S-S+
                        pos = binary_search (CurrentDets(:,1:TotWalkers), &
                                             ilut_srch, NIfD+1)
                        if (pos > 0) then

                            ! If we are looking for the spin of the initiator
                            ! only wavefunction, we need to ensure that we
                            ! are projecting onto an initiator...
                            inc = .true.
                            if (tTruncInitiator .and. only_init) then
                                if (test_flag(CurrentDets(:,pos), &
                                              flag_is_initiator(1)) .or. &
                                    test_flag(CurrentDets(:,pos), &
                                              flag_is_initiator(lenof_sign)))&
                                                                    then
                                    inc = .true.
                                else
                                    inc = .false.
                                end if
                            end if

                            call extract_sign (CurrentDets(:,pos), sgn2)
                            ssq = ssq + int(sgn(1) * sgn2(1) * sgn_hphf,int64) 
                        endif
                    endif
                enddo
            endif
        enddo

    end function

end module
