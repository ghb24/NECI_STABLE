! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
#include "macros.h"

module ParallelHelper
   use constants
   use iso_c_hack
   use timing_neci, only: timer, set_timer, halt_timer
    implicit none

    type(timer), save :: Sync_Time

    !
    ! If we are using C-bindings, certain things need to be defined
    !
#ifdef CBINDMPI

    ! These are not defined, if using MPI in C
    integer(MPIArg), parameter :: MPI_SUCCESS = 0
    integer(MPIArg) :: MPI_COMM_WORLD
    integer, parameter :: MPI_STATUS_SIZE = 1

! ****** HACK ********
! We would like to define these consts as here, but this breaks gfortran 4.5.1
! --> See macros.h
! ********************
!#if defined(__PATHSCALE__) || defined(__ISO_C_HACK) || defined (__OPEN64__)
!    c_ptr_t, parameter :: MPI_IN_PLACE = 0
!#else
!    c_ptr_t, parameter :: MPI_IN_PLACE = C_NULL_PTR
!#endif

    ! Define values so our C-wrapper can work nicely
    integer(MPIArg), parameter :: MPI_INTEGER4 = 0, &
                                  MPI_INTEGER8 = 1, &
                                  MPI_DOUBLE_PRECISION = 2, &
                                  MPI_DOUBLE_COMPLEX = 3, &
                                  MPI_2INTEGER = 4, &
                                  MPI_INTEGER = 5, &
                                  MPI_CHARACTER = 6, &
                                  MPI_2DOUBLE_PRECISION = 7

    ! Similarly for operations
    integer(MPIArg), parameter :: MPI_SUM = 0, &
                                  MPI_MAX = 1, &
                                  MPI_MIN = 2, &
                                  MPI_LOR = 3, &
                                  MPI_MINLOC = 4, &
                                  MPI_MAXLOC = 5

    ! Useful lengths
    integer, parameter :: MPI_MAX_ERROR_STRING = 500

    interface
        subroutine MPI_Error_string (err, s, l, ierr) &
            bind(c, name='mpi_error_string_wrap')
            use iso_c_hack
            integer(c_int), intent(in), value :: err, l
            integer(c_int), intent(out) :: ierr
            character(c_char), intent(out) :: s(*)
        end subroutine
        subroutine MPI_Barrier (comm, ierr) bind(c, name='mpi_barrier_wrap')
            use iso_c_hack
            integer(c_int), intent(in), value :: comm
            integer(c_int), intent(out) :: ierr
        end subroutine
        subroutine MPI_Init (ierr) bind(c, name='mpi_init_wrap')
            use iso_c_hack
            integer(c_int), intent(out) :: ierr
        end subroutine
        subroutine MPI_Finalize (ierr) bind(c, name='mpi_finalize_wrap')
            use iso_c_hack
            integer(c_int), intent(out) :: ierr
        end subroutine
        subroutine MPI_Abort (comm, err, ierr) bind(c, name='mpi_abort_wrap')
            use iso_c_hack
            integer(c_int), intent(in), value :: comm, err
            integer(c_int), intent(out) :: ierr
        end subroutine
        subroutine MPI_Comm_rank (comm, rank, ierr) &
            bind(c, name='mpi_comm_rank_wrap')
            use iso_c_hack
            integer(c_int), intent(in), value :: comm
            integer(c_int), intent(out) :: rank, ierr
        end subroutine
        subroutine MPI_Comm_size (comm, sz, ierr) &
            bind(c, name='mpi_comm_size_wrap')
            use iso_c_hack
            integer(c_int), intent(in), value :: comm
            integer(c_int), intent(out) :: sz, ierr
        end subroutine
        subroutine MPI_Comm_create (comm, group, ncomm, ierr) &
            bind(c, name='mpi_comm_create_wrap')
            use iso_c_hack
            integer(c_int), intent(in), value :: comm, group
            integer(c_int), intent(out) :: ncomm, ierr
        end subroutine
        subroutine MPI_Group_incl (grp, n, rnks, ogrp, ierr) &
            bind(c, name='mpi_group_incl_wrap')
            use iso_c_hack
            integer(c_int), intent(in), value :: grp, n
            integer(c_int), intent(in) :: rnks(*)
            integer(c_int), intent(out) :: ierr, ogrp
        end subroutine
        subroutine MPI_Comm_group (comm, group, ierr) &
            bind(c, name='mpi_comm_group_wrap')
            use iso_c_hack
            integer(c_int), intent(in), value :: comm
            integer(c_int), intent(out) :: group, ierr
        end subroutine
        subroutine MPI_Gather (sbuf, scnt, stype, rbuf, rcnt, rtype, &
                                 rt, comm, ierr) &
                                 bind(c, name='mpi_gather_wrap')
            use iso_c_hack
            use constants
            c_ptr_t, intent(in), value :: sbuf, rbuf
            integer(c_int), intent(in), value :: scnt, stype, rcnt, rtype
            integer(c_int), intent(in), value :: comm, rt
            integer(c_int), intent(out) :: ierr
        end subroutine
        function mpicommworld_c2f () result(cw) &
            bind(c, name='mpicommworld_c2f')
            use constants
            integer(MPIArg) :: cw
        end function
    end interface
#endif

#ifdef PARALLEL
    ! MpiDetInt needs to be defined here, so that it can make use of the
    ! above
#ifdef __INT64
    integer(MPIArg), parameter :: MpiDetInt = MPI_INTEGER8
#else
    integer(MPIArg), parameter :: MpiDetInt = MPI_INTEGER4
#endif

! SDS: We just no longer use MPI logical variables, and work around with
!      integers instead
!! This is a hack to work around disagreement between compilers on what
!! datatype is acceptable for logical variables in MPI routines.
!#ifdef __MPILOGTYPE
!    integer(MPIArg), parameter :: MPI_LOGTYPE4 = MPI_LOGICAL4
!    integer(MPIArg), parameter :: MPI_LOGTYPE8 = MPI_LOGICAL8
!#else
!    integer(MPIArg), parameter :: MPI_LOGTYPE4 = MPI_INTEGER4
!    integer(MPIArg), parameter :: MPI_LOGTYPE8 = MPI_INTEGER8
!#endif
#else
    ! In serial, set this to a nonsense value
    integer(MPIArg), parameter :: MpiDetInt = -1
#endif

#ifndef PARALLEL
    ! These don't exist in serial, so fudge them
    integer(MPIArg), parameter :: MPI_2INTEGER=0
    integer(MPIArg), parameter :: MPI_2DOUBLE_PRECISION=0
    integer(MPIArg), parameter :: MPI_MIN=0
    integer(MPIArg), parameter :: MPI_MAX=0
    integer(MPIArg), parameter :: MPI_SUM=0
    integer(MPIArg), parameter :: MPI_LOR=0
    integer(MPIArg), parameter :: MPI_MAXLOC=0
    integer(MPIArg), parameter :: MPI_MINLOC=0
    integer(MPIArg), parameter :: MPI_MAX_ERROR_STRING=255
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

   ! A null-info structure
   integer(MPIArg)      :: mpiInfoNull

   ! A 'node' which communicates between roots on each node
   type(CommI)          :: Roots


contains
   Subroutine GetComm(Comm,Node,rt,tMe)
       implicit none
       type(CommI), intent(in),optional :: Node
       integer(MPIArg), intent(out) :: Comm
       integer(MPIArg), intent(out), optional :: rt
       logical, intent(in), optional :: tMe
       logical tMe2
       if(present(tMe)) then
         tMe2=tMe
       else
         tMe2=.false.
       endif
       if(nNodes==0) then
         Comm=CommGlobal
         if(present(rt)) then
            if(tMe2) then
               rt=int(iProcIndex,MPIArg)
            else
               rt=Root
            endif
         endif
         return
       endif
            
       if (present(Node)) then
         if(Node%n==Roots%n) then
            Comm=CommRoot
            if (present(rt)) then
               if(tMe2) then
                  rt=int(iNodeIndex,MPIArg)
               else
                  rt=Root
               endif
            endif
         else
            Comm=int(CommNodes(Node%n),MPIArg)
            if (present(rt)) then
               if(tMe2) then
                  rt=int(iIndexInNode,MPIArg)
               else
                  rt=0 !NodeRoots(Node%n) is the procindex of the root, but not the index within the communicator
               endif
            endif
         endif
       else
         Comm=CommGlobal
         if (present(rt)) then
            if(tMe2) then
               rt=int(iProcIndex,MPIArg) 
            else
               rt=Root
            endif
         endif
       endif
   end subroutine



    subroutine MPIErr (iunit, err)
        integer, intent(in) :: err, iunit
        integer(MPIArg) :: l, e
#ifdef PARALLEL
        character(len=MPI_MAX_ERROR_STRING) :: s

        l=0
        e=0
        call MPI_Error_string (int(err, MPIArg), s, l, e)

        write(iunit,*) s
#endif

    end subroutine



    subroutine MPIBarrier (err, Node, tTimeIn)

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

#ifdef PARALLEL
        call GetComm (comm, node)

        call MPI_Barrier (comm, ierr)
        err = ierr
#else
        err = 0
#endif

        if (tTime) call halt_timer(Sync_Time)

    end subroutine

    subroutine MPIGroupIncl (grp, n, rnks, ogrp, ierr) 

        integer, intent(in) :: grp, n
        integer, intent(in) :: rnks(:)
        integer, intent(out) :: ierr
        integer(MPIArg), intent(out) :: ogrp
        integer(MPIArg) :: err

#ifdef PARALLEL
        call MPI_Group_incl (int(grp, MPIArg), int(n, MPIArg), &
                int(rnks, MPIArg), ogrp, err)
        ierr = err
#else
        ogrp = 0
        ierr = 0
#endif

    end subroutine

    subroutine MPICommcreate (comm, group, ncomm, ierr)

        integer(MPIArg), intent(in) :: comm
        integer(MPIArg), intent(in) :: group
        integer(MPIArg), intent(out) :: ncomm
        integer, intent(out) :: ierr
        integer(MPIArg) :: err

#ifdef PARALLEL
        call MPI_Comm_create (int(comm, MPIArg), int(group, MPIArg), &
                              ncomm, err)
        ierr = err
#else
        ncomm = 0
        ierr = 0
#endif

    end subroutine

    subroutine MPICommGroup (comm, grp, ierr)

        integer(MPIArg), intent(in) :: comm
        integer(MPIArg), intent(out) :: grp
        integer, intent(out) :: ierr
        integer(MPIArg) :: err, gout

#ifdef PARALLEL
        call MPI_Comm_Group (comm, gout, err)
        ierr = err
        grp = gout
#else
        grp = 0
        ierr = 0
#endif

    end subroutine
!
! n.b HACK
! We need to be able to do a bit of hackery when using C-based MPI
!
! --> We relabel things a bit...
#ifdef CBINDMPI
#define val_in vptr
#define val_out rptr
#else
#define val_in v
#define val_out Ret
#endif

    subroutine MPIGather_hack (v, ret, nchar, nprocs, ierr, Node)

        integer, intent(in) :: nchar, nprocs
        character(len=nchar), target :: v
        character(len=nchar), target :: ret(nprocs)
        integer, intent(out) :: ierr
        type(CommI), intent(in), optional :: Node
        integer(MPIArg) :: Comm, rt, err

#ifdef CBINDMPI
        character(c_char) :: in_tmp(nchar)
        character(len=nchar*nprocs), target :: out_tmp
        integer :: i, st, fn
#endif


#ifdef PARALLEL
#ifdef CBINDMPI
#ifdef __GFORTRAN__
        type(c_ptr) :: g_loc
#endif
        c_ptr_t :: vptr, rptr
        vptr = loc_neci(v)
        rptr = loc_neci(ret)
#endif

        call GetComm (Comm, Node, rt)

        call MPI_Gather (val_in, int(nchar, MPIArg), MPI_CHARACTER, &
                         val_out, int(nchar, MPIArg), MPI_CHARACTER, &
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
#ifdef PARALLEL
#ifdef CBINDMPI
        interface
            ! Put this here to avoid polluting the global namespace
            subroutine MPI_Allreduce_rt(val, ret, cnt, dtype, op, comm, ierr)&
                bind(c, name='mpi_allreduce_wrap')
                use iso_c_hack
                use constants
                integer(MPIArg), intent(in) :: val
                integer(MPIArg), intent(out) :: ret
                integer(c_int), intent(in), value :: cnt, dtype, op, comm
                integer(c_int), intent(out) :: ierr
            end subroutine
        end interface
        call MPI_Allreduce_rt (rt, nrt, 1_MPIArg, MPI_INTEGER, MPI_MAX, &
                               comm, ierr)
#else
        call MPI_Allreduce (rt, nrt, 1_MPIArg, MPI_INTEGER, MPI_MAX, &
                            comm, ierr)
#endif
#else
        ierr=0  !Avoid compiler warnings
        nrt=rt
#endif

    end subroutine

end module

subroutine mpibarrier_c (error) bind(c)
    use ParallelHelper, only: MPIBarrier
    use constants
    use iso_c_hack
    implicit none
    integer(c_int), intent(inout) :: error
    integer :: ierr

#ifdef PARALLEL
    call MPIBarrier (ierr)
    error = int(ierr,kind=kind(error))
#else
    error = 0
#endif
end subroutine
