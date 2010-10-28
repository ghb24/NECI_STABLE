#include "macros.h"

module hist

    use DeterminantData, only: get_lexicographic
    use MemoryManager
    use SystemData, only: tHistSpinDist, ilut_spindist, nbasis, nel, LMS, &
                          hist_spin_dist_iter, nI_spindist
    use DetBitOps, only: count_open_orbs, EncodeBitDet, spatial_bit_det, &
                         DetBitEq
    use util_mod, only: choose, get_free_unit, binary_search
    use constants, only: n_int, bits_n_int, size_n_int, lenof_sign
    use bit_rep_data, only: NIfTot, NIfD
    use bit_reps, only: extract_sign, encode_sign

    implicit none

    ! Should we histogram the distribution of spin dets within a given
    ! spatial structure --> Analyse spin development
    integer(n_int), allocatable, target :: hist_spin_dist(:,:)
    integer :: tag_spindist=0

contains

    subroutine init_hist_spin_dist ()

        integer :: open_orbs(nel), dorder(nel), ierr, orb, i
        integer :: nfound, nopen, nup, ndets
        integer(n_int), pointer :: p_ilut(:)
        integer(n_int) :: ilut_tmp(0:NIfTot)
        character(*), parameter :: this_routine = 'init_hist_spin_dist'

        ! Encode spin distribution
        if (.not. allocated(ilut_spindist)) &
            allocate(ilut_spindist(0:NIfTot))
        call EncodeBitDet(nI_spindist(1:nel), ilut_tmp)
        print*, 'encoded', ilut_tmp
        call writebitdet(6, ilut_tmp, .true.)
        call flush(6)
        ilut_spindist = spatial_bit_det(ilut_tmp)
        print*, 'encoded', ilut_spindist
        call writebitdet(6, ilut_spindist, .true.)
        call flush(6)

        ! Initialise the spin distribution histograms
        if (.not. allocated(ilut_spindist)) &
            call stop_all (this_routine, 'Spin distribution must be specified&
                                & and initialised to use spin distribution &
                                &histogramming.')

        ! How many open shell electrons are there in this spatial config
        nopen = count_open_orbs(ilut_spindist)
        nup = (nopen + LMS) / 2

        ! How many dets with this config?
        ndets = int(choose(nopen, nup))
        print*, 'nop, up, dets', nopen, nup, ndets
        call flush(6)

        ! Allocate and initialise storage
        allocate(hist_spin_dist(0:NIfTot, ndets), stat=ierr)
        LogAlloc(ierr, 'hist_spin_dist', (niftot+1)*ndets, size_n_int, tag_spindist)
        hist_spin_dist = 0

        ! Obtain a list of the open spatial orbitals.
        ! N.B. This is a ZERO BASED list
        nfound = 0
        do i = 1, nbasis-1, 2
            if (IsOcc(ilut_spindist, i) .and. IsNotOcc(ilut_spindist, i+1)) then
                nfound = nfound + 1
                open_orbs(nfound) = i
                if (nfound == nopen) exit
            endif
        enddo
        print*, 'open orbs', open_orbs(1:nopen)

        ! Generate all possible determinants with the given spatial structure
        ! and value of Ms
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
                clr_orb(p_ilut, ab_pair(orb))
            enddo
            call writebitdet(6, p_ilut, .false.)
            print*, p_ilut

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
        integer :: fd, i, sgn(lenof_sign)
        character(22) :: fname, iterstr

        ! Open the file for writing
        fd = get_free_unit ()
        write(iterstr, '(i12)') iter
        write(fname, '("spin-hist-",a)') trim(adjustl(iterstr))
        open(unit=fd, file=trim(fname), status='replace')

        ! Output file header
        write(fd, '("# 1. Determinant(niftot+1)\t2. sign(lenof_sign)")')

        do i = 1, ubound(hist_spin_dist, 2)

            ! Extract sign
            call extract_sign (hist_spin_dist(:,i), sgn)

            ! Output to file
            write(fd, *) hist_spin_dist(0:NIfD, i), float(sgn)/nsteps

            ! Zero the element for the next time around
            sgn = 0
            call encode_sign(hist_spin_dist(:,i), sgn)

        enddo

        ! Close the file
        close(fd)

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
            call encode_sign(hist_spin_dist(:,pos), sgn)
        endif

    end subroutine

end module
