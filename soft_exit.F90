module soft_exit

! Based almost entirely on work by GHB.

!= We create a temporary file called SOFTEXIT, which is removed at the end of
!= the calculation.  The existence of this file can be tested via the function 
!= TestSoftExit, which allows a clean exit to be performed, e.g. to enable
!= termination of a MC calculation before the number of steps specified in the
!= input file have been performed.

!= 1. At the start of a calculation, create the file SOFTEXIT::
!=        call init_soft_exit()
!= 2. During the calculation, test for the existence of the file SOFTEXIT.  If
!=    test_SOFTEXIT is false, then SOFTEXIT has been removed and the
!=    programmer's code should finish up loose ends (e.g. and post-processing of data
!=    required) and exit cleanly.
!=    If you wish to exit immediately, rather than trapping the value returned by
!=    test_SOFTEXIT and then tracing back through the calling stack and hence avoid
!=    any further unwanted calculations, use test_SOFTEXIT with the termination
!=    routines provided in NECICore.F90 and in error_handling.F90::
!=        if (.not.test_SOFTEXIT) then
!=            call NECICalcEnd(0)
!=            call NECICodeEnd()
!=            call quiet_stop()
!=        end if
!= 3. At the end of the calculation, if SOFTEXIT is still around, we remove it
!=    because it's polite to tidy up after ourselves::
!=        call end_soft_exit()

!= To do:
!=  * Parallel version.  The file SOFTEXIT is created on all nodes.  However,
!=    test_SOFTEXIT should be false if SOFTEXIT is removed from *any* node.
!=    This involves MPI communication (and JSS doesn't know all the MPI routines by
!=    heart yet).

implicit none

contains

    subroutine init_soft_exit()
       != Create the SOFTEXIT file.

       implicit none

       open (22,file='SOFTEXIT',status='unknown')
       write (22,*) ""
       close (22)

    end subroutine init_soft_exit



    subroutine end_soft_exit()
       != Remove the SOFTEXIT file to avoid polluting the user's disk.

       implicit none

       open(22,file='SOFTEXIT',status='unknown')
       close(22,status='delete')

    end subroutine end_soft_exit



    logical function test_SOFTEXIT()
       != Test if the file SOFTEXIT still exists.
       != Return True means that the file is still there and the calculation should 
       != proceed as normal. False means that you should start to exit cleanly.

       implicit none
       logical :: exists

       inquire(file='SOFTEXIT',exist=exists)
       test_SOFTEXIT=exists

    end function test_SOFTEXIT

end module soft_exit
