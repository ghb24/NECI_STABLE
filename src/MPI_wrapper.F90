#include "macros.h"

module MPI_wrapper
    use constants
!All use of mpi routines come from this module
#if defined(USE_MPI)
    use mpi
#endif
    use timing_neci, only: timer, set_timer, halt_timer
    implicit none

    type(timer), save :: Sync_Time

    !
    ! If we are using C-bindings, certain things need to be defined
    !

#ifdef USE_MPI
    ! MpiDetInt needs to be defined here, so that it can make use of the
    ! above
#ifdef INT64_
    integer(MPIArg), parameter :: MpiDetInt = MPI_INTEGER8
#else
    integer(MPIArg), parameter :: MpiDetInt = MPI_INTEGER4
#endif

#else
    ! In serial, set this to a nonsense value
    integer(MPIArg), parameter :: MpiDetInt = -1
#endif

#ifndef USE_MPI
    ! These don't exist in serial, so fudge them
    integer(MPIArg), parameter :: MPI_2INTEGER = 0
    integer(MPIArg), parameter :: MPI_2DOUBLE_PRECISION = 0
    integer(MPIArg), parameter :: MPI_MIN = 0
    integer(MPIArg), parameter :: MPI_MAX = 0
    integer(MPIArg), parameter :: MPI_SUM = 0
    integer(MPIArg), parameter :: MPI_LOR = 0
    integer(MPIArg), parameter :: MPI_MAXLOC = 0
    integer(MPIArg), parameter :: MPI_MINLOC = 0
    integer(MPIArg), parameter :: MPI_MAX_ERROR_STRING = 255
#endif

    Type :: CommI
        Integer n
    End Type
    ! Rank of the root processor
    integer, parameter :: root = 0
    integer              :: iProcIndex
    integer              :: nNodes            !The total number of nodes
    integer              :: iIndexInNode      !The index (zero-based) of this processor in its node
    integer iNodeIndex  ! Set from ParallelHelper.  Use this if an integer rather than a CommI object is needed.
    type(CommI)          :: Node              !The index of this node - this is a type to allow overloading
    logical              :: bNodeRoot         !Set if this processor is root of its node
    integer, allocatable :: CommNodes(:)      !Each node has a separate communicator
    integer(MPIArg), allocatable :: GroupNodes(:)     !Each node has a separate communicator
    integer(MPIArg), allocatable :: GroupNodesDum(:), CommNodesDum(:)
    type(CommI), allocatable :: Nodes(:)      !The node for each processor
    integer, allocatable :: ProcNode(:)   !The node for each processor (as a zero-based integer)
    integer, allocatable :: NodeRoots(:)      !The root for each node (zero-based)
    integer, allocatable :: NodeLengths(:)    !The number of procs in each node

    ! A communicator to all processors
    integer(MPIArg)      :: CommGlobal

    ! A group with all node roots in it
    integer(MPIArg)      :: GroupRoots

    ! A communicator between the roots on each node
    integer(MPIArg)      :: CommRoot

    ! Communicator/indices for MPI3 version of shared memory communication.
    ! Probably this can eventually be merged with the variables above
    integer(MPIArg):: mpi_comm_inter, mpi_comm_intra
    integer(MPIArg):: iProcIndex_inter, iProcIndex_intra

    ! A null-info structure
    integer(MPIArg)      :: mpiInfoNull

    ! A 'node' which communicates between roots on each node
    type(CommI)          :: Roots

contains
    Subroutine GetComm(Comm, Node, rt, tMe)
        implicit none
        type(CommI), intent(in), optional :: Node
        integer(MPIArg), intent(out) :: Comm
        integer(MPIArg), intent(out), optional :: rt
        logical, intent(in), optional :: tMe
        logical tMe2
        if (present(tMe)) then
            tMe2 = tMe
        else
            tMe2 = .false.
        end if
        if (nNodes == 0) then
            Comm = CommGlobal
            if (present(rt)) then
                if (tMe2) then
                    rt = int(iProcIndex, MPIArg)
                else
                    rt = Root
                end if
            end if
            return
        end if

        if (present(Node)) then
            if (Node%n == Roots%n) then
                Comm = CommRoot
                if (present(rt)) then
                    if (tMe2) then
                        rt = int(iNodeIndex, MPIArg)
                    else
                        rt = Root
                    end if
                end if
            else
                Comm = int(CommNodes(Node%n), MPIArg)
                if (present(rt)) then
                    if (tMe2) then
                        rt = int(iIndexInNode, MPIArg)
                    else
                        rt = 0 !NodeRoots(Node%n) is the procindex of the root, but not the index within the communicator
                    end if
                end if
            end if
        else
            Comm = CommGlobal
            if (present(rt)) then
                if (tMe2) then
                    rt = int(iProcIndex, MPIArg)
                else
                    rt = Root
                end if
            end if
        end if
    end subroutine

    subroutine MPIErr(iunit, err)
        integer, intent(in) :: err, iunit
        integer(MPIArg) :: l, e
#ifdef USE_MPI
        character(len=MPI_MAX_ERROR_STRING) :: s

        l = 0
        e = 0
        call MPI_Error_string(int(err, MPIArg), s, l, e)

        write(iunit, *) s
#endif

    end subroutine

    subroutine MPIBarrier(err, Node, tTimeIn)

        integer, intent(out) :: err
        type(CommI), intent(in), optional :: Node
        logical, intent(in), optional :: tTimeIn
        integer(MPIArg) :: comm, ierr
        logical :: tTime

        ! By default, do time the call.
        if (.not. present(tTimeIn)) then
            tTime = .true.
        else
            tTime = tTimeIn
        end if

        if (tTime) call set_timer(Sync_Time)

#ifdef USE_MPI
        call GetComm(comm, node)

        call MPI_Barrier(comm, ierr)
        err = ierr
#else
        err = 0
#endif

        if (tTime) call halt_timer(Sync_Time)

    end subroutine

    subroutine MPIGroupIncl(grp, n, rnks, ogrp, ierr)

        integer, intent(in) :: grp, n
        integer, intent(in) :: rnks(:)
        integer, intent(out) :: ierr
        integer(MPIArg), intent(out) :: ogrp
        integer(MPIArg) :: err

#ifdef USE_MPI
        call MPI_Group_incl(int(grp, MPIArg), int(n, MPIArg), &
                            int(rnks, MPIArg), ogrp, err)
        ierr = err
#else
        ogrp = 0
        ierr = 0
#endif

    end subroutine

    subroutine MPICommcreate(comm, group, ncomm, ierr)

        integer(MPIArg), intent(in) :: comm
        integer(MPIArg), intent(in) :: group
        integer(MPIArg), intent(out) :: ncomm
        integer, intent(out) :: ierr
        integer(MPIArg) :: err

#ifdef USE_MPI
        call MPI_Comm_create(int(comm, MPIArg), int(group, MPIArg), &
                             ncomm, err)
        ierr = err
#else
        ncomm = 0
        ierr = 0
#endif

    end subroutine

    subroutine MPICommGroup(comm, grp, ierr)

        integer(MPIArg), intent(in) :: comm
        integer(MPIArg), intent(out) :: grp
        integer, intent(out) :: ierr
        integer(MPIArg) :: err, gout

#ifdef USE_MPI
        call MPI_Comm_Group(comm, gout, err)
        ierr = err
        grp = gout
#else
        grp = 0
        ierr = 0
#endif

    end subroutine

    subroutine MPIGather_hack(v, ret, nchar, nprocs, ierr, Node)

        integer, intent(in) :: nchar, nprocs
        character(len=nchar), target :: v
        character(len=nchar), target :: ret(nprocs)
        integer, intent(out) :: ierr
        type(CommI), intent(in), optional :: Node
        integer(MPIArg) :: Comm, rt, err

#ifdef USE_MPI
        call GetComm(Comm, Node, rt)

        call MPI_Gather(v, int(nchar, MPIArg), MPI_CHARACTER, &
                        Ret, int(nchar, MPIArg), MPI_CHARACTER, &
                        rt, comm, err)

        ierr = err
#else
        ret(1) = v
        ierr = 0
#endif
    end subroutine

    subroutine MPIAllreduceRt(rt, nrt, comm, ierr)
        integer(MPIArg), intent(in) :: rt, comm
        integer(MPIArg), intent(out) :: nrt, ierr
#ifdef USE_MPI
        call MPI_Allreduce(rt, nrt, 1_MPIArg, MPI_INTEGER, MPI_MAX, &
                           comm, ierr)
#else
        ierr = 0
        nrt = rt
#endif

    end subroutine

end module

subroutine mpibarrier_c(error) bind(c)
    use MPI_wrapper, only: MPIBarrier
    use constants
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    integer(c_int), intent(inout) :: error
    integer :: ierr

#ifdef USE_MPI
    call MPIBarrier(ierr)
    error = int(ierr, kind=kind(error))
#else
    error = 0
#endif
end subroutine
