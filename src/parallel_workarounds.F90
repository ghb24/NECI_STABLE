#include "macros.h"

module Parallel_neci

    ! This is a wrapper module around Parallel_neci_internal, located in
    ! Parallel.F90.template.
    !
    ! As a result of the aggressive templating, and relatively inselective
    ! use statements (designed to pull all the relevant stuff through the
    ! templating structure), this can cause a certain number of conflicts.
    ! Especially using the *.h config files
    !
    ! --> This module acts as a filter. All the things we want to expose are
    !     explicitly used. And then the routines that were tempramental with
    !     broad imports are expressed.

    ! n.b. be explicit about uses from Parallel, to avoid multiple
    !      imports of mpi_ functions which annoy ifort
    use MPI_wrapper
    use par_internal, only: MPIBCast, MPIAllReduce, MpiAllGather, &
                            MPIReduce, MPISum, MPISumAll, MPIScatter, MPIAllGather, &
                            MPIAllGatherV, MPIGather, MPIGatherV, MPIScatterV, &
                            MPIAllReduceDatatype, MPIAllToAll, MPIAllToAllV, &
                            MPIStopAll, MPINodes, MPIInit, MPIEnd, clean_parallel, &
                            MPISend, MPIRecv, GetProcElectrons, nProcessors, &
                            neci_MPIInit_called, neci_MPINodes_called

    use constants
    use error_handling_neci, only: stop_all
    implicit none

    interface MPIBcast
        module procedure MPIBcast_character

        module procedure MPIBCastLogical
        module procedure MPIBCastLogical_logic
        module procedure MPIBCastLogicalArr
        module procedure MPIBCastLogicalArr_logic
    end interface

    interface MPIAllGather
        module procedure MPIAllGatherLogical
    end interface

contains

    ! These MPI functions are here to avoid ifort getting overly upset about
    ! re-use of mpi_functions if we put it in the parallel supermodule
    ! (the mpi_ functions are imported via too many routes).

    subroutine MPIBcast_character(v, Node)

        ! A special case for broadcasting string data, which doesn't apply
        ! to other datatypes.

        ! n.b. we may need to explicitly broadcast the length of the strig,
        !      as that is not part of the normal data.

        character(*), intent(inout), target :: v
        integer(MPIArg) :: ierr, comm, rt
        type(CommI), intent(in), optional :: Node
        character(*), parameter :: t_r = 'MPIBcast_character'
        integer(MPIArg) :: length

#ifdef USE_MPI
        call GetComm(Comm, Node, rt)

        ! Broadcast the length
        length = int(len_trim(v), MPIArg)
        call MPIBCast(length, node)

        call MPI_BCast(v, length, MPI_CHARACTER, rt, comm, ierr)

        ! Use v, not v, as v is explicitly a character array.
        v = v(1:length)

        if (ierr /= MPI_SUCCESS) &
            call stop_all(t_r, 'MPI Error. Terminating')
#endif
    end subroutine

    subroutine MPIAllLORLogical(param_in, param_out, node)

        logical, intent(in) :: param_in
        logical, intent(out) :: param_out

        type(CommI), intent(in), optional :: Node
        integer(MPIArg) :: comm
        integer :: v, ret, ierr

#ifdef USE_MPI

        call GetComm(Comm, Node)

        ! We cast the logical to an integer. It is a bit yucky, but avoids
        ! oddities in behaviour with some MPI implementations
        if (param_in) then
            v = 1
        else
            v = 0
        end if

        call MPIAllReduce(v, MPI_SUM, ret, node)

        param_out = (ret > 0)
#else
        param_out = param_in
#endif

    end subroutine

    subroutine MPIBCastLogical(param_inout, node)

        logical, intent(inout) :: param_inout
        type(CommI), intent(in), optional :: node
        character(*), parameter :: t_r = 'MPIBCastLogical'

        integer(MPIArg) :: ierr, comm, rt
        logical(int32), target :: v

#ifdef USE_MPI
        call GetComm(comm, node, rt)
        v = param_inout

        ! Use MPI_INTEGER instead of MPI_LOGICAL. Makes no difference for bcast
        call MPI_Bcast(v, 1_MPIArg, MPI_INTEGER4, rt, comm, ierr)

        if (Ierr /= MPI_SUCCESS) &
            call stop_all(t_r, 'MPIError. Terminating')

        param_inout = v
#endif

    end subroutine

    subroutine MPIBCastLogical_logic(param_inout, tMe, node)

        logical, intent(inout) :: param_inout
        type(CommI), intent(in), optional :: Node
        logical, intent(in) :: tMe
        integer(MPIArg) :: ierr, comm, rt, nrt
        character(*), parameter :: t_r = 'MPIBCastLogical_logic'
        logical(int32), target :: v

#ifdef USE_MPI
        call GetComm(Comm, Node, rt, tMe)

        v = param_inout

        ! Which processor is root?
        call MPIAllReducert(rt, nrt, comm, ierr)

        if (ierr == MPI_SUCCESS) then
            ! Use MPI_INTEGER instead of MPI_LOGICAL. Makes no difference for
            ! bcast
            call MPI_Bcast(v, 1_MPIArg, MPI_INTEGER4, nrt, comm, ierr)
        end if

        if (ierr /= MPI_SUCCESS) &
            call stop_all(t_r, "MPI Error. Terminating")

        param_inout = v
#endif

    end subroutine

    subroutine MPIBCastLogicalArr(param_inout, node)

        logical, intent(inout) :: param_inout(:)
        type(CommI), intent(in), optional :: node
        character(*), parameter :: t_r = 'MPIBCastLogical'

        integer(MPIArg) :: ierr, comm, rt, length
        logical(int32), target :: v(size(param_inout))

#ifdef USE_MPI
        call GetComm(comm, node, rt)

        v = param_inout

        length = int(size(v), MPIArg)
        ! Use MPI_INTEGER instead of MPI_LOGICAL. Makes no difference for bcast
        call MPI_Bcast(v, length, MPI_INTEGER4, rt, comm, ierr)

        if (Ierr /= MPI_SUCCESS) &
            call stop_all(t_r, 'MPIError. Terminating')

        param_inout = v
#endif

    end subroutine

    subroutine MPIBCastLogicalArr_logic(param_inout, tMe, node)

        logical, intent(inout) :: param_inout(:)
        type(CommI), intent(in), optional :: Node
        logical, intent(in) :: tMe
        integer(MPIArg) :: ierr, comm, rt, nrt, length
        character(*), parameter :: t_r = 'MPIBCastLogical'
        logical(int32), target :: v(size(param_inout))

#ifdef USE_MPI
        call GetComm(Comm, Node, rt, tMe)
        v = param_inout

        ! Which processor is root?
        call MPIAllReducert(rt, nrt, comm, ierr)

        if (ierr == MPI_SUCCESS) then
            length = int(size(v), MPIArg)
            ! Use MPI_INTEGER instead of MPI_LOGICAL. Makes no difference for
            ! bcast
            call MPI_Bcast(v, length, MPI_INTEGER4, nrt, comm, ierr)
        end if

        if (ierr /= MPI_SUCCESS) &
            call stop_all(t_r, "MPI Error. Terminating")
        param_inout = v
#endif

    end subroutine

    Subroutine MPIAllGatherLogical(param_in, param_out, ierr, Node)

        logical, intent(in), target :: param_in
        logical, intent(inout), target :: param_out(:)
        integer, intent(out) :: ierr
        type(CommI), intent(in), optional :: Node

        integer :: v, ret(lbound(param_out, 1):ubound(param_out, 1)), j

        ! Cast this gather statement via integers. A bit of a faff, but works
        ! around some nasties between MPI libraries - and it isn't performance
        ! limiting, so so what!

        if (param_in) then
            v = 1
        else
            v = 0
        end if

        call MPIAllGather(v, ret, ierr, node)

        do j = lbound(param_out, 1), ubound(param_out, 1)
            param_out(j) = ret(j) > 0
        end do

    end subroutine

end module
