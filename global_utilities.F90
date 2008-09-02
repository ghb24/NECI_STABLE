module global_utilities

!= A convenient wrapper module for collecting together the common
!= variables and procedures that need to be accessible from utility
!= modules.  Rather than doing, say,::
!=
!=    use MemoryManager, only: LogMemAlloc
!=    use timing, only: set_timer, halt_timer
!=
!= we can just USE global_utilities without too much namespace pollution::
!=
!=    use global_utilities

use MemoryManager, only: LogMemAlloc,LogMemDealloc
use timing, only: set_timer,halt_timer,get_total_time, timer

end module global_utilities
