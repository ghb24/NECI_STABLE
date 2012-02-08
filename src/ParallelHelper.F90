module ParallelHelper
   use constants
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
   integer, allocatable :: GroupNodes(:)     !Each node has a separate communicator
   type(CommI), allocatable :: Nodes(:)      !The node for each processor
   integer, allocatable :: ProcNode(:)   !The node for each processor (as a zero-based integer)
   integer, allocatable :: NodeRoots(:)      !The root for each node (zero-based)
   integer, allocatable :: NodeLengths(:)    !The number of procs in each node
   integer              :: CommGlobal        !A Communicator to all processors
   integer              :: GroupRoots        ! A group with all node roots in it
   integer              :: CommRoot          !A Communicator between the Roots on each nodes
   type(CommI)          :: Roots             !A 'node' which communiccates between roots on each node

contains
   Subroutine GetComm(Comm,Node,rt,tMe)
       implicit none
       type(CommI), intent(in),optional :: Node
       integer, intent(out) :: Comm
       integer, intent(out), optional :: rt
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

   subroutine MPIErr(error)
      INTEGER(MPIArg) :: error,l,e
#ifdef PARALLEL
      character(len=MPI_MAX_ERROR_STRING) :: s
#ifdef CBINDMPI
      e = MPI_ERROR_STRING(error,s,l)
#else
      call MPI_ERROR_STRING(error,s,l,e)
#endif
      write(6,*) s
#endif
   end subroutine

    Subroutine MPIBarrier(error,Node)
        INTEGER :: error
        INTEGER(MPIArg) :: err
        type(CommI), intent(in),optional :: Node
        integer Comm
        integer(MPIArg) :: CommArg
#ifdef CBINDMPI
        interface
            function MPI_Barrier (comm) result(ret) bind(c)
                integer, intent(in), value :: comm
                integer :: ret
            end function
        end interface
#endif
#ifdef PARALLEL
        call GetComm(Comm,Node)
        CommArg=Comm
#ifdef CBINDMPI
        err = MPI_Barrier(CommArg)
#else
        CALL MPI_Barrier(CommArg,err)
#endif
        error=err
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

    call MPIBarrier (ierr)
    error = ierr
end subroutine
