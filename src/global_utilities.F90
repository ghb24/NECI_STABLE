! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
module global_utilities

!= A convenient wrapper module for collecting together the common
!= variables and procedures that need to be accessible from utility
!= modules.  Rather than doing, say,::
!=
!=    use MemoryManager, only: LogMemAlloc
!=    use timing_neci, only: set_timer, halt_timer
!=
!= we can just USE global_utilities without too much namespace pollution::
!=
!=    use global_utilities
!=
!= Note that global_utilities doesn't contain all the functionality of all 
!= the utility modules, but rather access to the most common routines 
!= required (e.g. it doesn't include initialisation and termination routines,
!= which are typically only required once).

use MemoryManager, only: LogMemAlloc,LogMemDealloc
use timing_neci, only: set_timer,halt_timer,get_total_time, timer

end module global_utilities
