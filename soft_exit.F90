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

!= To do:
!=  * Parallel version.  The file SOFTEXIT is created on all nodes.  However,
!=    test_SOFTEXIT should be false if SOFTEXIT is removed from *any* node.
!=    This involves MPI communication (and JSS doesn't know all the MPI routines by
!=    heart yet).

implicit none

contains

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
