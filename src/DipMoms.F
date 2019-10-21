      subroutine null_wrapper
          implicit none
          character(*), parameter :: this_routine = 'null_wrapper'
          call stop_all(this_routine, 'Should not be here')
      end subroutine null_wrapper
