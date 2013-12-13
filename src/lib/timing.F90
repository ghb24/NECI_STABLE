module timing_neci
!= JSS.  Routines for timing code blocks.

!= To do: 
!=   *  parallel timing.  To be honest, I'm not sure how useful this will be,
!=      nor what form it should take.  Some kind of collation?  Average time
!=      over processors (mean/std deviation?)  Questions that can only be
!=      answered with some exerience, I suspect!

! ========================================================================
!= Usage:
!= 1. Start the global timer at the start of the calculation:
!=       call init_timing()
!= 
!= 2. For each procedure to be timed, declare a timer object and pass the
!=    procedure name:
!=       type(timer), save :: proc_timer
!=       proc_timer%timer_name='procedure name'
!=    And then time the procedure:
!=       call set_timer(proc_timer)
!=       [...]
!=       call halt_timer(proc_timer)
!=    This will time all the statements between the two timing calls.
!=    If a procedure is called between the timing calls which is itself 
!=    timed, then that is only included in one set of timings (in its own
!=    timer).
!=    The save is necessary for timing all calls to the same routine in one
!=    group.  The timer object is, however, very lightweight.
!=
!= 3. At the end of the calculation, stop the global timer and print out
!=    the timing report:
!=       call end_timing()
!=       call print_timing_report()
!=
!= 4. The time spent in a given timer can be obtained at runtime using
!=    the get_total_time function.  See this function for more details.
!=
!= set_timer and print_timing_report take optional arguments.
!= See the individual routines for more information.
! ========================================================================

use constants
implicit none
save

integer, parameter :: ntimer=200
integer :: itimer=0

type timer
    character(25) :: timer_name=''
    type(timer_object), pointer :: store=>null()
    logical :: time=.true. ! False if object is too low-level to be timed.
end type

type timer_object
    character(25) :: timer_name
    integer :: ncalls=0
    real(sp) :: time_cpu=0.0     ! For timing of the current
    real(sp) :: time_system=0.0  ! call to the procedure.
    real(sp) :: sum_time_cpu=0.0    ! Sum of time spent in the 
    real(sp) :: sum_time_system=0.0 ! procedure.
    logical :: timing_on=.false.   ! true whilst the timer is active.
end type timer_object

type(timer_object),allocatable,target :: timers(:)

! For total calculation time.
real(sp) :: global_time_cpu=0.0
real(sp) :: global_time_system=0.0
! If global_timing_on is true, then handle the total time differently in the timing output,
! as then have requested timing output without halting the global timer.
logical :: global_timing_on=.false. 

! To see if errors were encountered.
logical :: timer_error=.false.

contains



   subroutine init_timing()
      != Start global timer for timing the total calculation time.

      implicit none
      integer :: i

      call cpu_time(global_time_cpu)
      global_time_system=real(0.0,sp)
      global_timing_on=.true.

      if (.not.allocated(timers)) allocate(timers(ntimer))

      do i=1,itimer
          ! Have already done one run if itimer>0.  Clear existing timing info.
          timers(i)%ncalls=0
          timers(i)%time_cpu=0.0
          timers(i)%time_system=0.0
          timers(i)%sum_time_cpu=0.0
          timers(i)%sum_time_system=0.0
          timers(i)%timing_on=.false.
      end do

   end subroutine init_timing



   subroutine end_timing()
      != Stop global timer for timing the total calculation time.

      implicit none
      real(sp) :: t(2)

      if (global_timing_on) then
          call cpu_time(t(1))
          t(2)=real(0.0,sp)
          
          global_time_cpu=t(1)-global_time_cpu
          global_time_system=t(2)-global_time_system
          global_timing_on=.false.
      else
          call warning_neci('end_timing','Global timing never initialised via call to init_timing.')
      end if

   end subroutine end_timing



   subroutine set_timer(proc_timer,obj_level)
      != Start the timer for the specified object.
      != In:
      !=   obj_level (optional): timing level of the procedure.  Procedures with
      !=           a timing level above the iGlobalTimerLevel (specified in the
      !=           LOGGING block) are not timed.  The default timing level is 30.
      != In/Out:
      !=   proc_timer: The procedure timer.  Should contain the name of the
      !=           procedure and be SAVEd.  On exit, proc_timer%store points to the
      !=           appropriate entry in the timers array, which contains the
      !=           timing information for this object.  If the procedure is
      !=           called multiple times, the timer is not reinitialised, but
      !=           rather updated with new timing information (i.e. the current
      !=           timer is set).

      Use LoggingData, only: iGlobalTimerLevel
      implicit none
      type(timer) :: proc_timer
      integer, optional, intent(in) :: obj_level
      real(sp) :: t(2)
      integer :: timer_level

      if (.not.global_timing_on) then
          ! Initialise global timer.
          call init_timing()
      end if
      
      if (present(obj_level)) then
          timer_level=obj_level
      else
          timer_level=30
      end if

      if (timer_level.gt.iGlobalTimerLevel) then
          ! This object is too low-level to be timed.
          proc_timer%time=.false.
      else
          proc_timer%time=.true.
          if (.not.associated(proc_timer%store)) then
              ! Have a new object.
              itimer=itimer+1
              if (itimer.gt.ntimer) then
                  call warning_neci('set_timer','ntimer parameter too small for the number of objects to be timed.')
                  proc_timer%time=.false.
                  timer_error=.true.
                  return
              end if
              proc_timer%store=>timers(itimer)
          end if
          proc_timer%store%timer_name=proc_timer%timer_name
          proc_timer%store%ncalls=proc_timer%store%ncalls+1
          if (.not.proc_timer%store%timing_on) then
              ! Not in the middle of a recursive function.
              ! A recursive function will have the recursive section between the
              ! set_timer and halt_timer calls.  If we avoid refreshing the
              ! start time for the timer of the recursive procedure, then the
              ! correct timings are obtained.
              ! Start the clock.
              call cpu_time(t(1))
              t(2)=real(0.0,sp)
              proc_timer%store%time_cpu=t(1)
              proc_timer%store%time_system=t(2)
              proc_timer%store%timing_on=.true.
          end if
      end if

   end subroutine set_timer



   subroutine halt_timer(proc_timer)
      != Halt the timer for the specified object.
      != In/Out:
      !=   proc_timer: the timer object of the procedure.  Must be intialised by
      !=               set_timer.  The timer is stopped and the total cpu and
      !=               system time spent in the procedure is updated with the time
      !=               spent for the current call.
      
      implicit none
      type(timer), intent(inout) :: proc_timer
      integer :: i
      real(sp) :: t(2)
      real(sp) :: time_cpu,time_system

      if (.not.proc_timer%time) then
          ! Not timing this object: its level is below that of the
          ! iGLobalTimerLevel given via the logging option TIMING.
      else if (.not.associated(proc_timer%store)) then
          call warning_neci('halt_timer','proc_timer not intialised: '//proc_timer%timer_name)
          timer_error=.true.
      else
          call cpu_time(t(1))
          t(2)=real(0.0,sp)
          time_cpu=t(1)-proc_timer%store%time_cpu
          time_system=t(2)-proc_timer%store%time_system
          proc_timer%store%sum_time_cpu=proc_timer%store%sum_time_cpu+time_cpu
          proc_timer%store%sum_time_system=proc_timer%store%sum_time_system+time_system
          ! Have to remove the time spent in this routine from the other
          ! timers, so that the currently active timers exclude time spent
          ! in other timed procedures.
          do i=1,itimer
              timers(i)%time_cpu=timers(i)%time_cpu+time_cpu
              timers(i)%time_system=timers(i)%time_system+time_system
          end do
          ! Unset timer behaviour flags.
          proc_timer%store%timing_on=.false.
          proc_timer%time=.true.
      end if

   end subroutine halt_timer



   real(sp) function get_total_time(proc_timer,t_elapsed)
      != Return the (current) total time for a given timed procedure.
      != By default this does not include the elapsed time of the current
      != run of proc_timer's routine, so if proc_timer is active then
      != the default call to get_total_time returns the time spent in 
      != proc_timer%timer_name up to the most recent call.
      != In:
      !=   proc_timer: the timer object of the procedure.  Must be intialised by
      !=               set_timer.
      !=   t_elapsed(optional): include the elapsed time.  Warning: involves an
      !=               additional call to etime, so will affect performance if
      !=               called large numbers (10s of millions) of times.

      implicit none
      type(timer) :: proc_timer
      logical,optional :: t_elapsed
      real(sp) :: t(2)

      if (.not.associated(proc_timer%store)) then
          call warning_neci('get_total_time.','proc_timer not intialised: '//adjustl(proc_timer%timer_name))
          get_total_time=-1000.0 ! Helpfully return insane value, so it is obvious something went wrong. ;-)
      else
          get_total_time=proc_timer%store%sum_time_cpu+proc_timer%store%sum_time_system
          if (present(t_elapsed)) then
              if (t_elapsed) then
                  call cpu_time(t(1))
                  t(2)=real(0.0,sp)
                  get_total_time=get_total_time+t(1)+t(2)-proc_timer%store%time_cpu-proc_timer%store%time_system
              end if
          end if
      end if

   end function get_total_time



   subroutine print_timing_report(ntimer_objects,iunit)
      != Output a timing report.
      != If the global timer has been turned off (ie end_timing has been
      != called), then also deallocate the timers array).
      != In:
      !=    ntimer_objects (optional): the timing report prints out the objects 
      !=    took the largest amount of time in total.  ntimer_objects gives the 
      !=    number of objects to print out, in descending order of total time.
      !=    Default value: 10, as set in Logging module.
      !=    iunit (optional): file unit to which the timing  report is printed.
      !=    Default value: 6 (stdout).

      Use LoggingData, only: nPrintTimer,iGlobalTimerLevel
      implicit none
      integer, optional, intent(in) :: ntimer_objects
      integer, optional, intent(in) :: iunit
      integer :: io=6
      integer :: nobjs
      integer :: i,it,id(1)
      real(sp) :: t(2)
      real(sp) :: sum_times(ntimer),total_cpu,total_system
      integer :: date_values(8)

      ! Add on a small perturbation for the cases where the total time is 
      ! zero to single-precision.  This forces the procedure times to be printed
      ! out, if required, even if they are 0.0000, by avoiding issues with
      ! maxloc as the elements of the sum_times array are set to zero.
      sum_times=timers(:)%sum_time_system+timers(:)%sum_time_cpu+real(1.e-4,sp)

      if (present(iunit)) io=iunit
      if (present(ntimer_objects)) then
          nobjs=ntimer_objects
      else
          nobjs=nPrintTimer
      end if

      write (io,'(/a65)') '================================================================'
      
      write (io,'(a15/)') 'Timing report.'
      if (timer_error) write (io,'(a61/)') 'Timer encountered errors.  The following might be incorrect.'
      if (min(itimer,nobjs).gt.0) then
          write (io,'(a37/)') 'Timing of most expensive procedures.'
          write (io,'(a65)') 'Procedure                    Calls       CPU    system     total'
          write (io,'(a65)') '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
          
          total_cpu=0.0_dp
          total_system=0.0_dp
          do i=1,min(itimer,nobjs)
              ! Find i-th most expensive procedure.
              id=maxloc(sum_times)
              it=id(1)
              sum_times(it)=0.0_dp ! Don't find this object again.
              if(timers(it)%ncalls.gt.0) then
                  write (io,'(1X,a25,i9,3f10.2)') adjustl(timers(it)%timer_name),timers(it)%ncalls,  &
                                                 timers(it)%sum_time_cpu,timers(it)%sum_time_system,&
                                                 timers(it)%sum_time_cpu+timers(it)%sum_time_system
              endif
              total_cpu=total_cpu+timers(it)%sum_time_cpu
              total_system=total_system+timers(it)%sum_time_system
          end do
          write (io,'(a65)') '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
          write (io,'(a35,3f10.2/)') 'Total                             ',total_cpu,total_system,total_cpu+total_system
      end if
      if (.not.global_timing_on) then
          write (io,'(a20,f10.2)') 'Global CPU time    ',global_time_cpu
          write (io,'(a20,f10.2)')  'Global system time ',global_time_system
          write (io,'(a20,f10.2)')  'Global total time  ',global_time_cpu+global_time_system
      else
          call cpu_time(t(1))
          t(2)=real(0.0,sp)
          write (io,'(/a20,f10.2)') 'Global CPU time    ',t(1)-global_time_cpu
          write (io,'(a20,f10.2)') 'Global system time',t(2)-global_time_system
      end if
      write (io,'(a65)') '================================================================'

      if (.not.global_timing_on) deallocate(timers) ! Assume we're done as end_timing has been called.

      call date_and_time(VALUES=date_values)
      write(io,'(1X,"Calculation ended",1X,i2.2,"/",i2.2,"/",i4.4,1X,"at",1X,i2.2,2(":",i2.2))') &
            date_values(3:1:-1),date_values(5:7)
      write (io,'(a65)') '================================================================'


   end subroutine print_timing_report



end module timing_neci
