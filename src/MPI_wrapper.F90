#include "macros.h"

module MPI_wrapper
    use constants, only: dp, MPIArg, int32
!All use of mpi routines come from this module
#if defined(USE_MPI)
    use mpi_f08
#endif
    use timing_neci, only: timer, set_timer, halt_timer
    better_implicit_none

    type(timer), save :: Sync_Time

    !
    ! If we are using C-bindings, certain things need to be defined
    !

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
    type(MPI_Comm), allocatable :: CommNodes(:)      !Each node has a separate communicator
    type(MPI_Group), allocatable :: GroupNodes(:)     !Each node has a separate communicator
    type(MPI_Group), allocatable :: GroupNodesDum(:)
    type(MPI_Comm), allocatable :: CommNodesDum(:)
    type(CommI), allocatable :: Nodes(:)      !The node for each processor
    integer, allocatable :: ProcNode(:)   !The node for each processor (as a zero-based integer)
    integer, allocatable :: NodeRoots(:)      !The root for each node (zero-based)
    integer, allocatable :: NodeLengths(:)    !The number of procs in each node

    ! A communicator to all processors
    type(MPI_Comm)      :: CommGlobal

    ! A group with all node roots in it
    type(MPI_Group)     :: GroupRoots

    ! A communicator between the roots on each node
    type(MPI_Comm)      :: CommRoot

    ! Communicator/indices for MPI3 version of shared memory communication.
    ! Probably this can eventually be merged with the variables above
    type(MPI_Comm):: mpi_comm_inter, mpi_comm_intra
    integer(MPIArg):: iProcIndex_inter, iProcIndex_intra

    ! A null-info structure
    type(MPI_Info)      :: mpiInfoNull

    ! A 'node' which communicates between roots on each node
    type(CommI)          :: Roots

contains
    Subroutine GetComm(Comm, Node, rt, tMe)
        type(CommI), intent(in), optional :: Node
        type(MPI_Comm), intent(out) :: Comm
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
                Comm = CommNodes(Node%n) !int(CommNodes(Node%n), MPIArg)
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
        call MPI_Error_string(err, s, l, e)

        write(iunit, *) s
#endif

    end subroutine

    subroutine MPIBarrier(err, Node, tTimeIn)

        integer, intent(out) :: err
        type(CommI), intent(in), optional :: Node
        logical, intent(in), optional :: tTimeIn
        integer(MPIArg) :: ierr
        type(MPI_Comm) :: comm
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

        type(MPI_Group), intent(in) :: grp 
        integer, intent(in) :: n
        integer, intent(in) :: rnks(:)
        integer, intent(out) :: ierr
        type(MPI_Group), intent(out) :: ogrp
        integer(MPIArg) :: err

#ifdef USE_MPI
        call MPI_Group_incl(grp, int(n, MPIArg), &
                            int(rnks, MPIArg), ogrp, err)
        ierr = err
#else
        ogrp = 0
        ierr = 0
#endif

    end subroutine

    subroutine MPICommcreate(comm, group, ncomm, ierr)

        type(MPI_Comm), intent(in) :: comm
        type(MPI_Group), intent(in) :: group
        type(MPI_Comm), intent(out) :: ncomm
        integer, intent(out) :: ierr
        integer(MPIArg) :: err

#ifdef USE_MPI
        call MPI_Comm_create(comm, group, &
                             ncomm, err)
        ierr = err
#else
        ncomm = 0
        ierr = 0
#endif

    end subroutine

    subroutine MPICommGroup(comm, grp, ierr)

        type(MPI_Comm), intent(in) :: comm
        type(MPI_Group), intent(out) :: grp
        integer, intent(out) :: ierr
        type(MPI_Group) :: gout
        integer(MPIArg) :: err

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
        type(MPI_Comm) :: Comm
        integer(MPIArg) :: rt, err

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
        integer(MPIArg), intent(in) :: rt
        type(MPI_Comm), intent(in) :: comm
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
