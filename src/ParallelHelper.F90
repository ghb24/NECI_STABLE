module ParallelHelper
   use constants
   use iso_c_hack
    implicit none



    !
    ! If we are using C-bindings, certain things need to be defined
    !
#ifdef CBINDMPI

    ! These are not defined, if using MPI in C
    integer(MPIArg), parameter :: MPI_SUCCESS = 0
    integer(MPIArg), parameter :: MPI_COMM_WORLD = 0
    integer, parameter :: MPI_STATUS_SIZE = 1

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
        subroutine MPI_Allreduce_rt (val, ret, cnt, dtype, op, comm, ierr) &
            bind(c, name='mpi_allreduce_wrap')
            use iso_c_hack
            use constants
            integer(MPIArg), intent(in) :: val
            integer(MPIArg), intent(out) :: ret
            integer(c_int), intent(in), value :: cnt, dtype, op, comm
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
    end interface
#endif

    ! MpiDetInt needs to be defined here, so that it can make use of the
    ! above
#ifdef PARALLEL
#ifdef __INT64
    integer(MPIArg), parameter :: MpiDetInt = MPI_INTEGER8
#else
    integer(MPIArg), parameter :: MpiDetInt = MPI_INTEGER4
#endif
#else
    ! In serial, set this to a nonsense value
    integer(MPIArg), parameter :: MpiDetInt = -1
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
               rt=iProcIndex
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
                  rt=iNodeIndex
               else
                  rt=Root
               endif
            endif
         else
            Comm=CommNodes(Node%n)
            if (present(rt)) then
               if(tMe2) then
                  rt=iIndexInNode 
               else
                  rt=0 !NodeRoots(Node%n) is the procindex of the root, but not the index within the communicator
               endif
            endif
         endif
       else
         Comm=CommGlobal
         if (present(rt)) then
            if(tMe2) then
               rt=iProcIndex 
            else
               rt=Root
            endif
         endif
       endif
   end subroutine



    subroutine MPIErr (err)
        integer, intent(in) :: err
        integer(MPIArg) :: l, e
#ifdef PARALLEL
        character(len=MPI_MAX_ERROR_STRING) :: s

        call MPI_Error_string (int(err, MPIArg), s, l, e)

        write(6,*) s
#endif

    end subroutine



    subroutine MPIBarrier (err, Node)

        integer, intent(out) :: err
        type(CommI), intent(in), optional :: Node
        integer(MPIArg) :: comm, ierr

#ifdef PARALLEL
        call GetComm (comm, node)

        call MPI_Barrier (comm, ierr)
        err = ierr
#endif

    end subroutine

    subroutine MPIGroupIncl (grp, n, rnks, ogrp, ierr) 

        integer, intent(in) :: grp, n
        integer, intent(in) :: rnks(:)
        integer, intent(out) :: ierr
        integer(MPIArg), intent(out) :: ogrp
        integer(MPIArg) :: err, out_grp

        call MPI_Group_incl (int(grp, MPIArg), int(n, MPIArg), &
                int(rnks, MPIArg), ogrp, err)
        ierr = err

    end subroutine

    subroutine MPICommcreate (comm, group, ncomm, ierr)

        integer(MPIArg), intent(in) :: comm
        integer(MPIArg), intent(in) :: group
        integer(MPIArg), intent(out) :: ncomm
        integer, intent(out) :: ierr
        integer(MPIArg) :: err

        call MPI_Comm_create (int(comm, MPIArg), int(group, MPIArg), &
                              ncomm, err)
        ierr = err

    end subroutine

    subroutine MPICommGroup (comm, grp, ierr)

        integer(MPIArg), intent(in) :: comm
        integer(MPIArg), intent(out) :: grp
        integer, intent(out) :: ierr
        integer(MPIArg) :: err, gout

        call MPI_Comm_Group (comm, gout, err)
        ierr = err
        grp = gout

    end subroutine

end module

subroutine mpibarrier_c (error) bind(c)
    use ParallelHelper, only: MPIBarrier
    use constants
    use iso_c_hack
    implicit none
    integer(c_int), intent(inout) :: error
    integer :: ierr

    call MPIBarrier (ierr)
    error = ierr
end subroutine
