! JSS.  Memory book-keeping routines.  Contains some elements of the initialisation, 
! output and structure of the memory_manager module from CamCASP (formerly SITUS), 
! written by Alston Misquitta, with permission.

! To do:
!   TICK   Peak memory table.
!   TICK   Large memory table.
!   TICK   Have log as a cache (psuedo-LIFO): only store active allocations.
!   TICK   Debug option: store as much as the cache can hold.
!   TICK   Optional argument: calls from a routine to be counted.

! Log memory usage in one of two ways:
!   1. Store everything.  The size of the log (MaxLen) had best be suitably large.
!   2. Store active allocations.  When the top most slot in use in the log is
!   deallocated, free up all the slots at the top of the log that have been
!   deallocated for logging later actions (i.e. psuedo-LIFO [last in, first
!   out]).  This does allow a fractured log to occur, but memory
!   allocation-deallocation is often LIFO and it allows us to be efficient in
!   storing the actions without wasting much effort searching the log.
! This is controlled by the CachingMemLog (true==approach 2) flag.

module MemoryManager

implicit none

private
public :: MemoryLeft, MemoryUsed,MaxMemory,li,LookupPointer
! Allow users to do the potentially dangerous thing of changing how the log is run.
! We'll hope they'll only use this for good...
public :: CachingMemLog 
! Routines that need to be accessible.
public :: InitMemoryManager,LogMemAlloc,LogMemDealloc,LeaveMemoryManager

integer, parameter :: li = selected_int_kind(18) !ints between +-10^18

type MemLogEl
    character(len=25) :: ObjectName=''
    character(len=25) :: AllocRoutine=''
    character(len=25) :: DeallocRoutine='not deallocated'
    integer(li) :: ObjectSize=0
end type MemLogEl

! All memory variables stored in bytes.
integer(li), save :: MaxMemory
integer(li), save :: MemoryUsed
integer(li), save :: MemoryLeft
integer(li), save :: MaxMemoryUsed

! Warnings, debug flags, output parameters.
logical, save :: initialised = .false.
logical, save :: warned = .false.
integer, parameter :: MaxWarn = 10 ! maximum number of low memory warning messages to be printed.
integer, parameter :: nLargeObjects = 10 ! maximum number of the largest memory allocations remember.
integer, save :: nWarn = 0
logical, save :: debug = .false.
logical, save :: CachingMemLog = .true. ! See above for how MemLog is used.

! Log of memory allocations.
integer, parameter :: MaxLen = 5000
type(MemLogEl), save :: MemLog(MaxLen)
integer, save :: ipos=1  ! Next available empty slot in the log.

! Capture the state of the MemLog at peak usage.  Currently not outputted, but
! useful for diagnostics.
type(MemLogEl), save :: PeakMemLog(MaxLen)

type(MemLogEl), save :: LargeObjLog(nLargeObjects) ! Store the largest allocations.
integer, save :: ismall=1 ! The smallest large object (remember to avoid repeating minloc again and again...)

! For backwards compatibility with the existing scheme, where the IP address is
! stored as the tag. Use long integer (li) so can handle POINTER8.
integer(li), save :: LookupPointer(MaxLen)

contains

    subroutine InitMemoryManager(MaxMemBytes)
    ! Initialise memory manager.

    use common_routines, only: internal_error
    implicit none

    integer(li), intent(in) :: MaxMemBytes
    character(len=*), parameter :: ThisRoutine = 'InitMemoryManager'

    if (MaxMemBytes.le.0) then
        write (6,*) 'Illegal maximum memory value passed to memorymanager.'
        write (6,*) 'MaxMemgbytes = ',dfloat(MaxMemBytes)/(1024**2)
        call internal_error(ThisRoutine,__LINE__,&
           & 'Illegal maximum memory. Check MEMORY in your input file.')
    endif

    MaxMemory = MaxMemBytes
    MemoryUsed = 0
    MemoryLeft = MaxMemory
    MaxMemoryUsed = 0
    initialised = .true.
    nWarn = 0
!   Deal with debug options at a later date.
!   debug = gmemdebug

    write (6,'(a33,f8.1,a3)') ' Memory Manager initialised with ',dfloat(MaxMemBytes)/(1024**2),' MB'

    return
    end subroutine InitMemoryManager



    subroutine LogMemAlloc(ObjectName,ObjectSize,ElementSize,AllocRoutine,tag,nCalls)
    ! Log a memory allocation.
    ! INPUT:
    !       ObjectName - Name of object.
    !       ObjectSize - Number of elements in object.
    !       ElementSize - Number of bytes per element.
    !       AllocRoutine - routine in which object is allocated.
    ! OUTPUT:
    !       tag - position in memory log the object is stored at.
    !             If -1, then the log is full and it's not been stored.
    !       nCalls - optional.  Increments nCalls: counts the number of times
    !       a routine has called the LogMemAlloc routine.

    implicit none

    character(len=*),intent(in) :: ObjectName,AllocRoutine
    integer, intent(in) :: ObjectSize
    integer, intent(in) :: ElementSize
    integer, intent(out) :: tag
    integer, intent(inout), optional :: nCalls
 
    integer(li), parameter :: DefaultMem=1024**3
    integer :: ObjectSizeBytes,ismallloc(1)

    if (present(nCalls)) nCalls=nCalls+1

    if (.not.initialised) then
        write (6,*) 'Memory manager not initialised. Doing so now with 1GB limit.'
        call InitMemoryManager(DefaultMem)
    end if

    ObjectSizeBytes=ObjectSize*ElementSize

    MemoryUsed=MemoryUsed+ObjectSizeBytes
    MaxMemoryUsed=max(MemoryUsed,MaxMemoryUsed)
    MemoryLeft=MaxMemory-MemoryUsed

    if (MemoryLeft.lt.0.and.nWarn.lt.MaxWarn) then
        write (6,*) 'WARNING: Memory used exceeds maximum memory set',MemoryLeft
        nWarn=nWarn+1
    end if

    if (ipos.gt.MaxLen) then
        if (.not.warned) then
            warned=.true.
            write (6,*) 'Warning: Array capacity of memory manager exceeded.'
            write (6,*) 'Required array length is ',ipos
            write (6,*) 'Max memory used is likely to be incorrect.'
        end if
        tag=-1
        ! If we're not putting it in the log, test if it's a huge array:
        ! it's always the biggest fishes that get away!
        if (ObjectSizeBytes.gt.LargeObjLog(ismall)%ObjectSize) then
            LargeObjLog(ismall)%ObjectName=ObjectName
            LargeObjLog(ismall)%AllocRoutine=AllocRoutine
            LargeObjLog(ismall)%ObjectSize=ObjectSizeBytes
            ismallloc=minloc(LargeObjLog(:)%ObjectSize)
            ismall=ismallloc(1)
        end if
    else
        MemLog(ipos)%ObjectName=ObjectName
        MemLog(ipos)%AllocRoutine=AllocRoutine
        MemLog(ipos)%DeallocRoutine='not deallocated' ! In case this slot has already been used and abandoned.
        MemLog(ipos)%ObjectSize=ObjectSizeBytes
        tag=ipos
        ipos=ipos+1
    end if

    if (debug) then
        write (6,*) 'Allocating memory: ',tag,ObjectSizeBytes,ObjectName,AllocRoutine,MemoryUsed
    end if

    return
    end subroutine LogMemAlloc



    subroutine LogMemDealloc(DeallocRoutine,tag)
    ! Log a memory deallocation.
    ! INPUT:
    !       DeallocRoutine - routine in which object is deallocated.
    !       tag - position in memory log the object is stored at.
    ! OUTPUT:
    !       tag - 0 if successfully logged (or noted that it wasn't stored in
    !             the log in the first place).

    implicit none

    character(len=*), intent(in) :: DeallocRoutine
    integer, intent(inout) :: tag
    integer :: i,ismallloc(1)

    if (tag.eq.0) then
        write (6,*) 'Warning: attempting to log deallocation but never logged allocation.'
        tag=-1
    else if(tag.gt.MaxLen.or.tag.lt.-1) then
        write (6,*) 'Warning: attempting to log deallocation but tag does not exist: ',tag
        tag=-1
    else 

        if (MemoryUsed.eq.MaxMemoryUsed) then
            ! Are at peak memory usage.  Copy the memory log.
            ! Useful to see what's around when memory usage is at a maximum.
            PeakMemLog(:)=MemLog(:)
        end if

        MemoryUsed=MemoryUsed-MemLog(tag)%ObjectSize
        MemoryLeft=MaxMemory-MemoryUsed

        if (tag.eq.-1) then
            ! No record of it in the log: can only print out a debug message.
            if (debug) then
                write (6,*) 'Deallocating memory in: ',DeallocRoutine,tag
            end if
        else
            ! Object was stored in the cache.
            MemLog(tag)%DeallocRoutine=DeallocRoutine

            ! Check to see if object is larger than the smallest of the large
            ! objects: if so, keep a record of it.
            if (MemLog(tag)%ObjectSize.gt.LargeObjLog(ismall)%ObjectSize) then
                LargeObjLog(ismall)=MemLog(tag)
                ismallloc=minloc(LargeObjLog(:)%ObjectSize)
                ismall=ismallloc(1)
            end if
            
            if (CachingMemLog) then
                ! Then we free up this slot and slots of all objects directly below it in
                ! the log that have also been deallocated.  This is not the most
                ! efficient storage (we still can have a fractured log) but
                ! works well for LIFO approaches, which are most common for us.
                MemLog(tag)%ObjectSize=0 ! Nothing to see here now.
                if (tag.eq.ipos-1) then
                    do i=tag,1,-1
                        if (MemLog(i)%ObjectSize.eq.0) then
                            ipos=ipos-1
                        else
                            exit
                        end if
                    end do
                end if
            end if
            if (debug) then
                write (6,*) 'Deallocating memory: ',tag,MemLog(tag)
            end if
        end if

        ! Set tag to zero: there was no problem with the logging deallocation
        ! (apart from maybe a too small cache).
        tag=0

    end if

    return
    end subroutine LogMemDealloc



    subroutine LeaveMemoryManager
    ! Call this to print out the largest memory allocations.
    ! If debug flag is on, then the full memory log is dumped to file.

    use common_routines, only: getunit
    implicit none

    integer :: iunit,iobjloc(1),iobj,i
    integer(li) :: ObjectSizes(nLargeObjects+MaxLen)
    type(MemLogEl) :: AllMemEl(nLargeObjects+MaxLen)
    character(len=*), parameter :: memoryfile = 'TMPMemoryusage.dat'
    character(len=*), parameter :: fmt1='(3a19,f8.1)'

    if (MemoryUsed.eq.MaxMemoryUsed) then
        ! Peak memory usage is now.
        PeakMemLog(:)=MemLog(:)
    end if

    call WriteMemLogHeader(6)

    if (CachingMemLog) then
        ! Large objects might be residing in the MemLog, but not deallocated
        ! (and so haven't been moved to the large object store).
        AllMemEl(:)=(/ MemLog(:),LargeObjLog(:) /)
    else
        ! Everything really ought to be held in just the MemLog: if not, then
        ! this is a "feature".
        AllMemEl(1:MaxLen)=MemLog(:)
    end if

    ! Copy the sizes to an integer array: we use maxloc on the copy.  This
    ! allows us to check for arrays of the same size without writing over 
    ! information in our log.
    ObjectSizes(:)=AllMemEl(:)%ObjectSize
    iobjloc(:)=maxloc(ObjectSizes)
    iobj=iobjloc(1)
    ObjectSizes(iobj)=ObjectSizes(iobj)+1
    do i=2,min(nLargeObjects+1,ipos)
        write (6,fmt1) ' '//AllMemEl(iobj)%ObjectName,AllMemEl(iobj)%AllocRoutine,AllMemEl(iobj)%DeallocRoutine,dfloat(AllMemEl(iobj)%ObjectSize)/1024**2
        iobjloc=maxloc(ObjectSizes,mask=ObjectSizes.lt.ObjectSizes(iobj))
        iobj=iobjloc(1)
        ObjectSizes(iobj)=ObjectSizes(iobj)+1
    end do
    if (warned) then
        write (6,*) '== NOTE: Length of logging arrays exceeded. Length needed is ',ipos
    endif
    write (6,*) '================================================================'

    if (debug) then
        ! Dump entire memory log to file.
        call getunit(iunit)
        open(unit=iunit,file=memoryfile,form='formatted',status='unknown')
        call WriteMemLogHeader(iunit)
        do iobj = 1, min(ipos,MaxLen)
            write (6,fmt1) ' '//MemLog(iobj)%ObjectName,MemLog(iobj)%AllocRoutine,MemLog(iobj)%DeallocRoutine,dfloat(MemLog(iobj)%ObjectSize)/1024**2
        enddo
        if (warned) then
            write (iunit,*) '== NOTE: Length of logging arrays exceeded. Length needed is ',ipos
        endif
        write (iunit,*) '================================================================'
        close(iunit)
    end if

    return
    end subroutine LeaveMemoryManager

    

    subroutine WriteMemLogHeader(iunit)
    implicit none
    integer :: iunit
    write(iunit,*)
    write (iunit,*) '================================================================'
    write (iunit,*) 'Memory usage'
    write (iunit,'(a34,f9.1)') ' Maximum memory defined is (MB) : ',dfloat(MaxMemory)/1024**2
    write (iunit,'(a34,f9.1)') ' Maximum memory used is    (MB) : ',dfloat(MaxMemoryUsed)/1024**2
    if (nWarn.gt.0) then
        write (iunit,*)'Maximum memory exceeded ',nWarn,' times.'
    endif
    write (iunit,*) ''
    if (iunit.eq.6) then
        write (iunit,*) 'Large memory allocations:'
        write (iunit,*) ''
    end if
    write (iunit,*) 'Name              Allocated in       Deallocated in     Size(MB)'
    write (iunit,*) '----------------------------------------------------------------'
    return
    end subroutine WriteMemLogHeader



end module MemoryManager
