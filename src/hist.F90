#include "macros.h"

module hist

    use DeterminantData, only: get_lexicographic
    use MemoryManager
    use SystemData, only: tHistSpinDist, ilut_spindist, nbasis, nel, LMS, &
                          hist_spin_dist_iter, nI_spindist, LMS, tHPHF, &
                          tOddS_HPHF, G1
    use DetBitOps, only: count_open_orbs, EncodeBitDet, spatial_bit_det, &
                         DetBitEq
    use CalcData, only: tFCIMC
    use DetCalcData, only: FCIDetIndex, det
    use FciMCData, only: tFlippedSign, TotWalkers, CurrentDets, iter, &
                         norm_psi_squared
    use util_mod, only: choose, get_free_unit, binary_search
    use HPHFRandExcitMod, only: FindExcitBitDetSym
    use constants, only: n_int, bits_n_int, size_n_int, lenof_sign
    use bit_rep_data, only: NIfTot, NIfD
    use bit_reps, only: extract_sign, encode_sign, extract_bit_rep
    use parallel
    use csf, only: get_num_csfs, csf_coeff, csf_get_yamas, write_yama, &
                   extract_dorder
    use hist_data

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
            call stop_all ("Incorrect number of CSFs found")

        ! Allocate list of csf coefficients
        if (allocated(hist_csf_coeffs)) deallocate(hist_csf_coeffs)
        allocate(hist_csf_coeffs(ndets, ncsf), stat=ierr)
        LogAlloc(ierr, 'hist_csf_coeffs', ndets*ncsf, sizeof_dp, tag_histcsfs)

        ! Generate all possible determinants with the given spatial structure
        ! and value of Ms
        write(6,*) 'Initialising Spin distribution histogramming'
        write(fmt_str,'("(",i2,"i20,",i4,"f20.15)")') NIfD+1, ncsf
        call flush(6)
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
        integer :: fd, i, j, sgn(lenof_sign), ierr
        integer(n_int) :: all_hist(0:NIfTot, ubound(hist_spin_dist, 2))
        real(dp) :: csf_contrib(lenof_sign, ubound(hist_csf_coeffs, 2))
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
                write(fd, *) all_hist(0:NIfD, i), float(sgn)/nsteps
                
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

        integer(n_int), intent(in), target :: ilut(0:NIfD)
        integer(n_int), intent(out), pointer :: ptr(:)

        ptr => ilut

    end subroutine

    subroutine test_add_hist_spin_dist_det (ilut, sgn)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer(n_int) :: ilut_tmp(0:NIfTot)
        integer, dimension(lenof_sign), intent(in) :: sgn
        integer, dimension(lenof_sign) :: sgn_old
        integer :: pos, i

        ! Should we add this ilut to the histogram?
        ilut_tmp = spatial_bit_det (ilut)
        
        if (DetBitEq(ilut_tmp, ilut_spindist)) then
            pos = binary_search (hist_spin_dist, ilut, NIfD+1)
            if (pos < 0) then
                call writebitdet(6, ilut, .false.)
                write(6,*) ilut
                write(6,*) '----------------'
                do i=1,ubound(hist_spin_dist, 2)
                    write(6,*) hist_spin_dist(:,i)
                enddo
                call stop_all ("test_add_hist_spin_dist_det", &
                               "Determinant not found in spin histogram list")
            endif
            call extract_sign (hist_spin_dist(:,pos), sgn_old)
            call encode_sign(hist_spin_dist(:,pos), sgn + sgn_old)
        endif

    end subroutine

    subroutine add_hist_spawn (ilut, sgn, ExcitLevel, dProbFin)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: sgn(lenof_sign), ExcitLevel
        real(dp), intent(in) :: dProbFin

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
            delta = real(sgn(:), dp) / dProbFin

            if (tHPHF) then
                call FindExcitBitDetSym (ilut, ilut_sym)
                if (.not. DetBitEq(ilut, ilut_sym)) delta = delta / sqrt(2.0)
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
                    delta = (real(sgn(:), dp) / sqrt(2.0)) / dProbFin

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

    subroutine add_hist_energies (ilut, sgn, HDiag)

        ! This will histogram the energies of the particles, rather than the
        ! determinants themselves.

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, dimension(lenof_sign), intent(in) :: sgn
        real(dp), intent(in) :: HDiag
        integer :: bin
        character(*), parameter :: t_r = "add_hist_energies"

        bin = int(HDiag / BinRange) + 1
        if (bin > iNoBins) &
            call stop_all (t_r, "Histogramming energies higher than the &
                         &arrays can cope with. Increase iNoBins or BinRange")
        HistogramEnergy(bin) = HistogramEnergy(bin) + real(sum(abs(sgn)),dp)

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
        integer :: sgn(lenof_sign), dorder(nel), nI(nel)
        real(dp) :: coeff, S_coeffs(LMS:nel), norm, S2, S22

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

        do j = 1, TotWalkers

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

        call MPISum_inplace (S_coeffs)
        if (iProcIndex == Root) then
            norm = sum(S_coeffs)
            S_coeffs = sqrt(S_coeffs / norm)

            S2 = 0
            do S = LMS, nel, 2
                S2 = S2 + real(S * (S + 2) * S_coeffs(S)) / 4
            enddo
            
            write(6,*) 'Scoeffs', iter, S_coeffs
            write(6,*) 'Scsf2', S2
            write(6,*) 'norm compare', norm, norm_psi_squared
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
        integer :: j, nup, S, flg, ncsf, nopen
        integer :: nI(nel), sgn(lenof_sign), dorder(nel)

        ! Initially, have no component on any of the spins
        spin_cpts = 0

        do j = 1, TotWalkers

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

        call MPISum_inplace (spin_cpts)

        if (iProcIndex == Root) then
            norm = sum(spin_cpts)
            spin_cpts = sqrt(spin_cpts / norm)

            write(6,'(a,i8,13f15.10)') 'spn: ', iter, spin_cpts
        endif


    end subroutine
    
    subroutine sum_cpts_S (spin_cpts, dorder, ncsf, nopen, sgn, S)

        real(dp), intent(inout) :: spin_cpts(0:nel)
        integer, intent(in) :: ncsf, nopen, S
        integer, intent(in) :: dorder(1:nopen), sgn(lenof_sign)
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

            ! TODO: And not an open shell det.
!            if (tHPHF) then
!                !
!                ! Deal with HPHF
!                ! TODO: Does this work with both odd and even HPHF?
!                !
!                call FindExcitBitDetSym(detcurr, detsym)
!                call CalcOpenOrbs (detcurr, nop_pairs)
!                if (mod(OpenOrbs, 2) == 0) then
!                    pair_sgn = 1
!                else
!                    pair_sgn = -1
!                endif
!
!                do j = 1, nel
!                    orb = nI(j)
!                    sgn_carry = 1
!                    if (.not. IsDoub(detcurr, orb)) then
!                        if (is_beta(orb)) then
!                            splus = detcurr
!                            clr_orb(splus, orb)
!                            set_orb(splus, get_alpha(orb))
!                        else
!                            sgn_carry = sgn_carry * pair_sgn
!                            splus = detsym
!                            clr_orb(splus, get_alpah(orb))
!                            set_orb(splus, orb)
!                        endif
!                    endif
!
!
!                enddo
!
!            else
    
    function calc_s_squared () result (ssq)

        real(dp) :: ssq
        integer :: i, j, k, l, orb, orb2, pos, pair_sgn, nop_pairs
        integer :: sgn_carry
        integer(n_int) :: splus(0:nifd), sminus(0:nifd), detsym(0:nifd)
        integer(n_int), pointer :: detcurr(:)
        integer :: nI(nel), flg, sgn(lenof_sign), sgn2(lenof_sign)
        integer :: lms_tmp



        ! Loop over beta electrons, and consider promoting them to alpha
        ssq = 0
        do i = 1, TotWalkers

            call ilut_nifd_pointer_assign(detcurr, CurrentDets(0:NIfD, i))
            call extract_bit_rep (CurrentDets(0:NifD, i), nI, sgn, flg)

            do j = 1, nel
                if (is_beta(nI(j)) &
                    .and. IsNotOcc(detcurr, get_alpha(nI(j)))) then
                    splus = detcurr
                    clr_orb(splus, nI(j))
                    set_orb(splus, get_alpha(nI(j)))

                    do k = 1, nel
                        orb2 = nI(k)
                        if (k == j) orb2 = get_alpha(orb2)

                        if (is_alpha(orb2) &
                            .and. IsNotOcc(splus, get_beta(orb2))) then
                            sminus = splus
                            clr_orb(sminus, orb2)
                            set_orb(sminus, get_beta(orb2))

                            ! --> sminus is an allowed result of applying S-S+
                            pos = binary_search(CurrentDets(:,1:TotWalkers), &
                                            sminus, NIfD+1)
                            if (pos > 0) then
                                call extract_sign (CurrentDets(:,pos), sgn2)
                                ssq = ssq + (sgn(1) * sgn2(1))
                            endif
                        endif
                    enddo
                endif
            enddo
        enddo

        ! Sum over all processors and normalise
        call MPISum_inplace (ssq)
        ssq = ssq / norm_psi_squared

        ! TODO: n.b. This is a hack. LMS appears to contain -2*Ms of the system
        !            I am somewhat astounded I haven't noticed this before...
        lms_tmp = -LMS
        ssq = ssq + real(lms_tmp * (lms_tmp + 2), dp) / 4


    end function

end module
