! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
module default_sets

! Flags for setting up sets of defaults.

! By default, the default set of defaults is used, unless one of the flags
! below is set.

! The defaults are set in the relevant input parsing routine (e.g.
! SysReadInput,IntReadInput etc).

! Feb08 set:
!   * Use the Fock-Partition-Lowdiag partitioning scheme for evaluating
!     rho-matrix elements
!   * RhoEpsilon=10^-8, the threshold below which an element of the rho matrix
!     is taken to be zero.
!   * MCPATHS logging option is turned on.
logical :: Feb08 = .false.
logical :: Nov11 = .false.

end module default_sets
