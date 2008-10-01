module soft_exit

! Based almost entirely on work by GHB.

!= During the calculation, test for the existence of the file SOFTEXIT.  If
!= test_SOFTEXIT is true, then SOFTEXIT has been created (by, we hope, the
!= user) and the programmer's code should finish up loose ends (e.g. and 
!= post-processing of data required) and exit cleanly.
!=
!= The existence of this file can be tested via the function TestSoftExit,
!= which allows a clean exit to be performed, e.g. to enable termination
!= of a MC calculation before the number of steps specified in the input
!= file have been performed.
!=
!= If you wish to exit immediately, rather than trapping the value returned by
!= test_SOFTEXIT and then tracing back through the calling stack and hence avoid
!= any further unwanted calculations, use test_SOFTEXIT with the termination
!= routines provided in NECICore.F90 and in error_handling.F90::
!=     if (test_SOFTEXIT) then
!=         call NECICalcEnd(0)
!=         call NECICodeEnd(tCPMD,tVASP)
!=         call quiet_stop()
!=     end if
!= Alternatively, you can call the SOFTEXIT_exit subroutine, which does exactly
!= the above.

implicit none

contains

#if PARALLEL

    logical function test_SOFTEXIT()
       != Test if the file SOFTEXIT exists.
       != Return True means that the file is there and you should start to exit cleanly.
       != False means that the calculation should proceed as normal. 

       use Parallel
       implicit none
       logical :: exists
       logical :: AnyExist
       integer :: error

       inquire(file='SOFTEXIT',exist=exists)
       !This collective will check the exists logical on all nodes, and perform a logical or operation,
       !before broadcasting the result back to all nodes.
       CALL MPI_AllReduce(exists,AnyExist,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,error)
       test_SOFTEXIT=AnyExist
!       if(exists) THEN
           !At the moment, SOFTEXIT isn't removed after it is detected.
           !This is because it is likely to be detected on multiple processors and they can't all try to close the file
           !at the same time. Anyone have a good solution to this problem? Anyone?!
!           open(13,file='SOFTEXIT')
!           close(13,status='delete')
!       endif
       if(test_SOFTEXIT) then
           write (6,'(X,a)') 'Request for SOFTEXIT detected on a node.'
       endif

    end function test_SOFTEXIT
    
#else
    
    logical function test_SOFTEXIT()
       != Test if the file SOFTEXIT exists.
       != Return True means that the file is there and you should start to exit cleanly.
       != False means that the calculation should proceed as normal. 

       implicit none
       logical :: exists

       inquire(file='SOFTEXIT',exist=exists)
       test_SOFTEXIT=exists
       ! We'll also do some house-keeping whilst we're here.
       if (test_SOFTEXIT) then
           ! Remove it so it doesn't catch us out the next time the calculation
           ! is run!
           open(13,file='SOFTEXIT')
           close(13,status='delete')
           write (6,'(X,a)') 'Request for SOFTEXIT detected.'
       end if

    end function test_SOFTEXIT
    
#endif
    
    subroutine SOFTEXIT_exit()
       != If SOFTEXIT has been created, then exit cleanly and immediately, with
       != no further calculations.
       != If SOFTEXIT has not been created, do nothing.
    
        use SystemData, only: tCPMD,tVASP
        implicit none

        if (test_SOFTEXIT()) then
            call NECICalcEnd(0)
            call NECICodeEnd(tCPMD,tVASP)
            call quiet_stop()
        end if

    end subroutine SOFTEXIT_exit

end module soft_exit
