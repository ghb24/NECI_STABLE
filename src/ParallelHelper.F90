module ParallelHelper
   Type :: CommI
      Integer n
   End Type
   ! Rank of the root processor
   integer, parameter :: root = 0
   integer              :: nNodes            !The total number of nodes
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
   Subroutine GetComm(Comm,Node,rt)
       implicit none
       type(CommI), intent(in),optional :: Node
       integer, intent(out) :: Comm
       integer, intent(out), optional :: rt
       if (present(Node)) then
         if(Node%n==Roots%n) then
            Comm=CommRoot
            if (present(rt)) rt=Root
         else
            Comm=CommNodes(Node%n)
            if (present(rt)) rt=0 !NodeRoots(Node%n) is the procindex of the root, but not the index within the communicator
         endif
       else
         Comm=CommGlobal
         if (present(rt)) rt=Root
       endif
   end subroutine
end module
