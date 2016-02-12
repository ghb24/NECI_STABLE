#include "macros.h"

module rdm_parallel

    ! This module contains all routines used for the calculations of reduced
    ! density matrices, when distributed across all MPI processes.

    use bit_rep_data, only: NIfTot, NIfDBO
    use SystemData, only: NEl, nBasis
    use util_mod
    use constants

    implicit none

    subroutine init_parallel_rdms(nrdms)

        ! This routine initialises all arrays and data used for the parallel
        ! implementation of RDMs.

        use rdm_data, only: rdm_arr , rdm_arr_ht, rdm_spawn, rdm_spawn_ht
        use rdm_data, only: rdm_spawn_recv, rdm_arr_size

        integer, intent(in) :: nrdms

        ! For now, create RDM arrays big enough so that *all* RDM elements on
        ! a particular processor can be stored, using the usual approximations
        ! to take symmetry into account. Include a factor of 1.2 to account for
        ! factors such as imperfect load balancing (which affects the spawned
        ! array).
        rdm_arr_size = 1.2*nbasis*(nbasis-1)/(16*nProcessors)
        nhashes_rdm = 0.8*rdm_arr_size

        allocate(rdm_spawn(nrdms+2, rdm_arr_size))
        allocate(rdm_spawn_recv(nrdms+2, rdm_arr_size))
        allocate(rdm_arr(nrdms+2, rdm_arr_size))

        ! Allocate and initialise hash tables.
        allocate(rdm_spawn_ht(nhashes_rdm), stat=ierr)
        call init_hash_table(rdm_spawn_ht)

        allocate(rdm_spawn_ht(nhashes_rdm), stat=ierr)
        call init_hash_table(rdm_arr_ht)

    end subroutine init_parallel_rdms

end module rdm_parallel
