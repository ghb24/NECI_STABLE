#include "macros.h"

module hist

    use DeterminantData, only: get_lexicographic
    use MemoryManager
    use SystemData, only: tHistSpinDist, ilut_spindist, nbasis, nel, LMS
    use DetBitOps, only: count_open_orbs
    use util_mod, only: choose
    use constants, only: n_int, bits_n_int, size_n_int
    use bit_rep_data, only: NIfTot, NIfD

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
        character(*), parameter :: this_routine = 'init_hist_spin_dist'

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

        ! Allocate and initialise storage
        allocate(hist_spin_dist(0:NIfTot, ndets), stat=ierr)
        LogAlloc(ierr, 'hist_spin_dist', (niftot+1)*ndets, size_n_int, &
                 tag_spindist)
        hist_spin_dist = 0

        ! Obtain a list of the open spatial orbitals.
        ! N.B. This is a ZERO BASED list
        nfound = 0
        do i = 0, nbasis-2, 2
            if (IsOcc(ilut_spindist, i)) then
                nfound = nfound + 1
                open_orbs(nfound) = i
                if (nfound == nopen) exit
            endif
        enddo

        ! Generate all possible determinants with the given spatial structure
        ! and value of Ms
        nfound = 0
        dorder(1) = -1
        call get_lexicographic (dorder, nopen, nup)
        do while (dorder(1) /= -1)

            nfound = nfound + 1
            p_ilut => hist_spin_dist(0:NIfD, nfound)
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

            call get_lexicographic (dorder, nopen, nup)
        enddo

        if (nfound /= ndets) &
            call stop_all (this_routine, &
                           "Generated incorrect number of determinants")

    end subroutine

    subroutine clean_hist_spin_dist ()

        character(*), parameter :: this_routine = 'clean_hist_spin_dist'

        if (allocated(ilut_spindist)) &
            deallocate(ilut_spindist)

        if (allocated(hist_spin_dist)) then
            deallocate(hist_spin_dist)
            LogDealloc(tag_spindist)
        endif

    end subroutine

end module
