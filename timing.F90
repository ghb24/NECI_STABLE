module timing
! JSS.  Routines for timing code blocks.

! To do:
!  * convert NECI to new timers.
!  * save tag rather than have to search for name.

! ========================================================================
!= Usage:
!= 1. Start the global timer at the start of the calculation:
!=       call init_timing()
!= 
!= 2. For each procedure to be timed, declare:
!=       integer, save :: isub=0
!=    And then time the procedure:
!=       call set_timer('procedure name',isub)
!=       [...]
!=       call halt_timer(isub)
!=    This will time all the statements between the two timing calls.
!=    If a procedure is called between the timing calls which is itself 
!=    timed, then that is only included in one set of timings (in its own
!=    timer)..
!=
!= 3. At the end of the calculation, stop the global timer and print out
!=    the timing report:
!=       call end_timing()
!=       call print_timing_report()
!=
!= set_timer and print_timing_report take optional arguments.
!= See the individual routines for more information.
! ========================================================================

implicit none
save

integer, parameter :: ntimer=200
integer :: itimer=0

type timer_object
    character(25) :: timer_name
    integer :: ncalls=0
    real(4) :: time_cpu=0.0     ! For timing of the current
    real(4) :: time_system=0.0  ! call to the procedure.
    real(4) :: sum_time_cpu=0.0    ! Sum of time spent in the 
    real(4) :: sum_time_system=0.0 ! procedure.
    logical :: timing_on=.false.
end type timer_object

type(timer_object) :: timers(ntimer)

! For total calculation time.
real(4) :: global_time_cpu=0.d0
real(4) :: global_time_system=0.d0
! If global_timing_on is true, then handle the total time differently in the timing output,
! as then have requested timing output without halting the global timer.
logical :: global_timing_on=.false. 

contains



   subroutine init_timing()
      != Start global timer for timing the total calculation time.

      implicit none
      real(4) :: etime,s,t(2)

      s=etime(t)
      global_time_cpu=t(1)
      global_time_system=t(2)
      global_timing_on=.true.

   end subroutine init_timing



   subroutine end_timing()
      != Stop global timer for timing the total calculation time.

      implicit none
      real(4) :: etime,s,t(2)

      s=etime(t)
      global_time_cpu=t(1)-global_time_cpu
      global_time_system=t(2)-global_time_system
      global_timing_on=.false.

   end subroutine end_timing



   subroutine set_timer(obj_name,obj_id,obj_level)
      != Start the timer for the specified object.
      != In:
      !=   obj_name: Name of procedure to be timed.
      !=   obj_level (optional): timing level of the procedure.  Procedures with
      !=           a timing level above the iGlobalTimerLevel (specified in the
      !=           LOGGING block) are not timed.  The default timing level is 30.
      != In/Out:
      !=   obj_id: The timer id of the procedure being timed.  It should be
      !=           initialised to 0 in the calling routine and SAVEd.  set_timer 
      !=           returns the timer id of the procedure.  Subsequent calls to 
      !=           the same procedure then re-use the same timer id, to obtain
      !=           the total time spent in a given procudure.
      !=           If the timer_id is returned to be -127, then the procedure is
      !=           not timed.

      use Logging, only: iGlobalTimerLevel
      implicit none
      character(*), intent(in) :: obj_name
      integer, intent(inout) :: obj_id
      integer, optional, intent(in) :: obj_level
      real(4) :: etime,s,t(2)
      integer :: timer_level
      integer :: i
      
      if (present(obj_level)) then
          timer_level=obj_level
      else
          timer_level=30
      end if
      if (timer_level.gt.iGlobalTimerLevel) then
          ! This object is too low-level to be timed.
          ! Return a id value of -127 to indicate the object is not timed.
          obj_id=-127
      else
          if (obj_id.eq.0) then
              ! Have a new object.
              itimer=itimer+1
              obj_id=itimer
          end if
          if (itimer.gt.ntimer) then
              call stop_all('set_timer','ntimer parameter too small for the number of objects to be timed.')
          else if (obj_id.le.0) then
              call stop_all('set_timer','invalid timer id from'//obj_name)
          end if
          timers(obj_id)%timer_name=obj_name
          timers(obj_id)%ncalls=timers(obj_id)%ncalls+1
          if (.not.timers(obj_id)%timing_on) then
              ! Not in the middle of a recursive function.
              ! Start the clock.
              s=etime(t)
              timers(obj_id)%time_cpu=t(1)
              timers(obj_id)%time_system=t(2)
              timers(obj_id)%timing_on=.true.
          else
              obj_id=-127
          end if
      end if

   end subroutine set_timer



   subroutine halt_timer(obj_id)
      != Halt the timer for the specified object.
      != In:
      !=   obj_id: The timer id of the procedure being timed, as returned in
      !=           set_timer. 
      
      implicit none
      integer, intent(in) :: obj_id
      integer :: i
      real(4) :: etime,s,t(2)
      real(4) :: time_cpu,time_system

      if (obj_id.eq.-127) then
          ! Not timing this object: its level is below that of the
          ! iGLobalTimerLevel given via the logging option TIMING or
          ! in the middle of a recursive function.
      else if (obj_id.le.0.or.obj_id.gt.ntimer) then
          call stop_all('halt_timer','Invalid obj_id given to time object.')
      else
          s=etime(t)
          time_cpu=t(1)-timers(obj_id)%time_cpu
          time_system=t(2)-timers(obj_id)%time_system
          timers(obj_id)%sum_time_cpu=timers(obj_id)%sum_time_cpu+time_cpu
          timers(obj_id)%sum_time_system=timers(obj_id)%sum_time_system+time_system
          ! Have to remove the time spent in this routine from the other
          ! timers, so that the currently active timers exclude time spent
          ! in other timed procedures.
          do i=1,itimer
              timers(i)%time_cpu=timers(i)%time_cpu+time_cpu
              timers(i)%time_system=timers(i)%time_system+time_system
          end do
          timers(obj_id)%timing_on=.false.
      end if

   end subroutine halt_timer



   real(4) function get_total_time(obj_id)
      != Return the (current) total time for a given timed procedure.
      != In:
      !=   obj_id: The timer id of the procedure, as returned in set_timer. 
      implicit none
      integer :: obj_id
      get_total_time=timers(obj_id)%sum_time_cpu+timers(obj_id)%sum_time_system
   end function get_total_time



   subroutine print_timing_report(ntimer_objects,iunit)
      != Output a timing report.
      != In:
      !=    ntimer_objects (optional): the timing report prints out the objects 
      !=    took the largest amount of time in total.  ntimer_objects gives the 
      !=    number of objects to print out, in descending order of total time.
      !=    Default value: 10, as set in Logging module.
      !=    iunit (optional): file unit to which the timing  report is printed.
      !=    Default value: 6 (stdout).

      use Logging, only: nPrintTimer,iGlobalTimerLevel
      implicit none
      integer, optional, intent(in) :: ntimer_objects
      integer, optional, intent(in) :: iunit
      integer :: io=6
      integer :: nobjs
      integer :: i,it,id(1)
      real(4) :: etime,s,t(2)
      real(4) :: sum_times(ntimer),total_cpu,total_system

      sum_times=timers(:)%sum_time_system+timers(:)%sum_time_cpu

      if (present(iunit)) io=iunit
      if (present(ntimer_objects)) then
          nobjs=ntimer_objects
      else
          nobjs=nPrintTimer
      end if

      write (io,'(/a65)') '================================================================'
      write (io,'(a15/)') 'Timing report.'
      if (min(itimer,nobjs).gt.0) then
          write (io,'(a37/)') 'Timing of most expensive procedures.'
          write (io,'(a65)') 'Procedure                    Calls       CPU    system     total'
          write (io,'(a65)') '----------------------------------------------------------------'
          
          do i=1,min(itimer,nobjs)
              ! Find i-th most expensive procedure.
              id=maxloc(sum_times)
              it=id(1)
              sum_times(it)=0.d0 ! Don't find this object again.
              write (io,'(X,a25,i9,3f10.2)') trim(timers(it)%timer_name),timers(it)%ncalls,timers(it)%sum_time_cpu,timers(it)%sum_time_system,timers(it)%sum_time_cpu+timers(it)%sum_time_system
              total_cpu=total_cpu+timers(it)%sum_time_cpu
              total_system=total_system+timers(it)%sum_time_system
          end do
          write (io,'(a65)') '----------------------------------------------------------------'
          write (io,'(a35,3f10.2/)') 'Total                             ',total_cpu,total_system,total_cpu+total_system
      end if
      if (.not.global_timing_on) then
          write (io,'(a20,f10.2)') 'Global CPU time    ',global_time_cpu
          write (io,'(a20,f10.2)')  'Global system time ',global_time_system
          write (io,'(a20,f10.2)')  'Global total time  ',global_time_cpu+global_time_system
      else
          s=etime(t)
          write (io,'(/a20,f10.2)') 'Global CPU time    ',t(1)-global_time_cpu
          write (io,'(a20,f10.2)') 'Global system time',t(2)-global_time_system
      end if
      write (io,'(a65)') '================================================================'

   end subroutine print_timing_report



end module timing
