module aliasSampling
  ! This module contains the utility to use alias table lookup on lists,
  ! requiring to precompute biases but making the lookup O(1)
  use constants
  use shared_memory_mpi
  use ParallelHelper, only: iProcIndex_intra
  use dSFMT_interface , only : genrand_real2_dSFMT  
  implicit none

  private
  public :: aliasSampler_t, aliasTable_t, aliasSamplerArray_t, clear_sampler_array
  
  ! type for tables: contains everything you need to get a random number
  ! with given biases
  type aliasTable_t
     private
     ! WARNING: DO NOT MANUALLY RE-ASSIGN THESE POINTERS, THIS WILL MOST LIKELY BREAK STUFF     
     ! this is the table of bias
     real(dp), pointer :: biasTable(:) => null()
     ! this is the lookup table for the resulting random number
     integer, pointer :: aliasTable(:) => null()

     ! shared memory windows
     integer(MPIArg) :: biasTableShmw
     integer(MPIArg) :: aliasTableShmw
   contains
     ! constructor
     procedure :: setupTable
     ! only load the data, without allocation
     procedure :: initTable
     ! destructor - final would be suited better, but is not supported by all compilers
     procedure :: tableDestructor
     ! get a random value from the alias table
     procedure :: getRand

     ! this hurts, but the array class has to have a mean of setting the pointers
     ! WHY ARE THERE NO FRIENDS IN FORTRAN !!!! THIS IS NOT OOP
     procedure :: setBiasTablePtr
     procedure :: setAliasTablePtr
  end type aliasTable_t

  !------------------------------------------------------------------------------------------!

  ! type for alias samplers: given an (un-normalized) array of weights, draws
  ! random numbers from this distribution
  type aliasSampler_t
     private
     ! alias table used for sampling
     type(aliasTable_t) :: table
     ! WARNING: DO NOT MANUALLY RE-ASSIGN THIS POINTER, THIS WILL MOST LIKELY BREAK STUFF          
     ! the probabilities     
     real(dp), pointer :: probs(:) => null()
     ! the shm window for the probabilities
     integer(MPIArg) :: probsShmw

   contains
     ! constructor
     procedure :: setupSampler
     ! only load the data, without allocation
     procedure :: initSampler
     ! initialize the probabilities
     procedure :: initProbs
     ! destructor - see above re final
     procedure :: samplerDestructor
     ! get a random element and the generation probability
     procedure :: sample
     ! get the probability to produce a given value
     procedure :: getProb
     ! I'm so sorry, but this is required for the array class to access the probs pointer
     procedure :: setProbsPtr
     procedure :: setBiasPtr
     procedure :: setAliasPtr
  end type aliasSampler_t

  !------------------------------------------------------------------------------------------!
  ! sampler array class - required for technical reasons: we need multiple samplers to share
  ! the same shared memory windows
  !------------------------------------------------------------------------------------------!  

  type aliasSamplerArray_t
     private
     ! this is an array of aliasSamplers, in the end
     type(aliasSampler_t), allocatable :: samplerArray(:)

     ! shared resources of the array entries
     real(dp), pointer :: allProbs(:) => null()
     integer, pointer :: allAliasTable(:) => null()
     real(dp), pointer :: allBiasTable(:) => null()
     integer(MPIArg) :: allProbsShmw, allAliasShmw, allBiasShmw
   contains
     ! constructor
     procedure :: setupSamplerArray
     procedure :: setupEntry     
     ! destructor
     procedure :: samplerArrayDestructor
     ! get a random element and the generation probability from one of the samplers
     procedure :: aSample
     procedure :: aGetProb
  end type aliasSamplerArray_t

contains

  !------------------------------------------------------------------------------------------!
  ! Initialization / Finalization routines of the aliasTable
  !------------------------------------------------------------------------------------------!

  subroutine setupTable(this, arr)
    ! pseudo-constructor for alias tables      
    ! Input: arr - array containing the (not necessarily normalized) probabilities we
    !              want to use for sampling
    implicit none
    class(aliasTable_t) :: this
    real(dp), intent(in) :: arr(:)

    character(*), parameter :: t_r = "setupTable"
    integer(int64) :: arrSize
    if(sum(arr) < eps) call stop_all(t_r,&
         "Trying to setup empty alias table")

    ! allocate the shared memory segment for the alias table
    arrSize = size(arr)

    call safe_shared_memory_alloc(this%biasTableShmw, this%biasTable, arrSize)
    if(associated(this%aliasTable)) &
         call shared_deallocate_mpi(this%aliasTableShmw, this%aliasTable)    
    call shared_allocate_mpi(this%aliasTableShmw, this%aliasTable, (/arrSize/))

    call this%initTable(arr)

  end subroutine setupTable

  !------------------------------------------------------------------------------------------!

  subroutine initTable(this,arr)
    ! Input: arr - array containing the (not necessarily normalized) probabilities we
    !              want to use for sampling    
    implicit none
    class(aliasTable_t) :: this
    real(dp), intent(in) :: arr(:)
    
    integer :: i,j,cV, cU
    integer(int64) :: arrSize
    integer, allocatable :: overfull(:), underfull(:)

    arrSize = size(arr)

    ! as this is shared memory, only node-root has to do this
    if(iProcIndex_intra .eq. 0) then
       ! initialize the probabilities
       this%biasTable = arr/sum(arr)*arrSize

       ! indices of subarrays
       allocate(overfull(arrSize))
       allocate(underfull(arrSize))

       cV = 0
       cU = 0
       do i = 1, arrSize
          call assignLabel(i)
       end do
       ! we now labeled each entry

       ! it is more efficient to start with the largest biases
       ! -> reverse overfull
       overfull(1:cV) = overfull(cV:1:-1)
       do 
          if((cV .eq. 0) .or. (cU .eq. 0)) then            
             exit
          end if
          ! pick one overfull and one underfull index
          i = overfull(cV)
          j = underfull(cU)
          ! set the alias of the underfull to be the other
          this%aliasTable(j) = i
          ! correct the bias
          this%biasTable(i) = this%biasTable(i) + this%biasTable(j) - 1

          ! unmark j
          cU = cU - 1
          ! reassign i based on the new bias
          cV = cV - 1
          call assignLabel(i)
       end do

       ! make sure we do not leave anything unfilled
       call roundTo1(overfull,cV)
       call roundTo1(underfull,cU)

    endif

  contains 

    subroutine assignLabel(i)
      integer, intent(in) :: i

      if(this%biasTable(i) > 1) then
         cV = cV + 1
         overfull(cV) = i
      else
         cU = cU + 1
         underfull(cU) = i
      end if
    end subroutine assignLabel

    subroutine roundTo1(labels,cI)
      integer, intent(in) :: labels(:)
      integer, intent(in) :: cI

      ! if, due to floating point errors, one of the categories is not empty, empty it
      ! (error is negligible then)
      if(cI > 0) then
         do i = 1, cI
            this%biasTable(labels(i)) = 1.0_dp
            this%aliasTable(labels(i)) = labels(i)
         end do
      endif

    end subroutine roundTo1
        
  end subroutine initTable

  !------------------------------------------------------------------------------------------!

  subroutine tableDestructor(this)
    ! clear the memory used by the alias table
    implicit none
    class(aliasTable_t) :: this

    if(associated(this%aliasTable)) &
         call shared_deallocate_mpi(this%aliasTableShmw, this%aliasTable)
    call safe_shared_memory_dealloc(this%biasTableShmw, this%biasTable)
  end subroutine tableDestructor

  !------------------------------------------------------------------------------------------!
  ! Sampling function
  !------------------------------------------------------------------------------------------!

  function getRand(this) result(ind)
    ! Draw a random number from an alias table created with the corresponding probabilities
    ! A bit tricky: this cannot be pure because a random number has to be drawn from
    ! the rng
    ! (maybe can be resolved, but current implementation requires side-effects for random numbers)
    ! output: ind - random number between 1 and the size of the array used to create the
    !               aliasTable object
    implicit none    
    class(aliasTable_t) :: this
    integer :: ind
    real(dp) :: r, bias
    integer :: sizeArr, pos

    sizeArr = size(this%biasTable)
    ! random number between 0 and 1
    r = genrand_real2_dSFMT()
    ! random position in arr
    pos = int(sizeArr*r)+1
    ! remainder of the integer conversion
    bias = sizeArr*r + 1 - pos

    if(bias < this%biasTable(pos)) then
       ind = pos
    else
       ind = this%aliasTable(pos)
    endif
  end function getRand

  !------------------------------------------------------------------------------------------!
  ! Functions to .... set the pointers manually (do NOT use them unless you know what you do)
  !------------------------------------------------------------------------------------------!

  subroutine setBiasTablePtr(this,ptr)
    ! Input: ptr - new target for the biasTable
    implicit none
    class(aliasTable_t) :: this
    real(dp), pointer :: ptr(:)

    this%biasTable => ptr
  end subroutine setBiasTablePtr

  !------------------------------------------------------------------------------------------!
  
  subroutine setAliasTablePtr(this,ptr)
    ! Input: ptr - new target for the aliasTable
    implicit none
    class(aliasTable_t) :: this
    integer, pointer :: ptr(:)

    this%aliasTable => ptr
  end subroutine setAliasTablePtr
  
  !------------------------------------------------------------------------------------------!
  ! Initialization / Finalization routines of the sampler
  !------------------------------------------------------------------------------------------!

  subroutine setupSampler(this, arr)
    ! allocate the resources of this and load the probability distribution from arr into this
    ! Input: arr - array containing the (not necessarily normalized) probabilities we
    !              want to use for sampling
    implicit none
    class(aliasSampler_t) :: this
    real(dp), intent(in) :: arr(:)

    integer(int64) :: arrSize
    ! if all weights are 0, throw an error
    if(sum(arr) < eps) then
       write(iout,*) &
            "Warning: trying to initialize sampler with empty probability distribution"
       ! probs defaults to null(), so it is not associated at this point (i.e. in a well-defined state)
       return
    endif
    
    ! initialize the alias table
    call this%table%setupTable(arr)
    arrSize = size(arr)

    ! allocate the probabilities
    call safe_shared_memory_alloc(this%probsShmw, this%probs, arrSize)

    ! set the probabilities
    call this%initProbs(arr)
  end subroutine setupSampler

  !------------------------------------------------------------------------------------------!
  
  subroutine initSampler(this, arr)
    ! load the probability distribution from arr into this
    ! we only use this in the sampler array, but fortran has no friend classes, so its public
    ! Input: arr - array containing the (not necessarily normalized) probabilities we
    !              want to use for sampling    
    implicit none
    class(aliasSampler_t) :: this
    real(dp), intent(in) :: arr(:)

    ! if all weights are 0, throw an error
    if(sum(arr) < eps) then
       write(iout,*) &
            "Warning: trying to initialize sampler with empty probability distribution"
       ! probs defaults to null(), so it is not associated at this point (i.e. in a well-defined state)
       return
    endif

    ! load the data - assume this is pre-allocated
    call this%table%initTable(arr)
    call this%initProbs(arr)
  end subroutine initSampler

  !------------------------------------------------------------------------------------------!
  
  subroutine initProbs(this, arr)
    ! load the probability distribution from arr into this%probs
    ! Input: arr - array containing the (not necessarily normalized) probabilities we
    !              want to use for sampling    
    implicit none
    class(aliasSampler_t) :: this
    real(dp), intent(in) :: arr(:)

    ! the array is shared memory, so only node-root has to do this
    if(iProcIndex_intra.eq.0) then
       ! the probabilities are taken from input and normalized
       this%probs = arr / sum(arr)
    end if    
  end subroutine initProbs
  
  !------------------------------------------------------------------------------------------!

  subroutine samplerDestructor(this)
    implicit none
    class(aliasSampler_t) :: this

    ! free the stored probabilities
    call safe_shared_memory_dealloc(this%probsShmw, this%probs)
    ! free the corresponding alias table
    call this%table%tableDestructor()
  end subroutine samplerDestructor

  !------------------------------------------------------------------------------------------!
  
  subroutine sample(this, tgt, prob)
    ! draw a random element from 1:size(this%probs) with the probabilities listed in prob
    ! Input: tgt - on return, this is a random number in the sampling range of this
    !        prob - on return, the probability of picking tgt
    implicit none
    class(aliasSampler_t), intent(in) :: this
    integer, intent(out) :: tgt
    real(dp), intent(out) :: prob

    ! empty samplers don't return anything - since probs defaults to null(), this check is safe
    if(.not.associated(this%probs)) then
       tgt = 0
       prob = 1.0
       return
    end if
    ! get the drawn number from the alias table
    tgt = this%table%getRand()
    ! and its probability
    prob = this%probs(tgt)
  end subroutine sample

    !------------------------------------------------------------------------------------------!

  pure function getProb(this, tgt) result(prob)
    ! Returns the probability to draw tgt from this sampler
    ! Input: tgt - the number for which we request the probability of sampling
    ! Output: prob - the probability of drawing tgt with the sample routine
    implicit none
    class(aliasSampler_t), intent(in) :: this
    integer, intent(in) :: tgt
    real(dp) :: prob

    ! the probability of drawing anything from an empty sampler is 0
    if(.not.associated(this%probs)) then
       prob = 0.0
    else
       prob = this%probs(tgt)
    endif
  end function getProb

  !------------------------------------------------------------------------------------------!
  ! Manually set the member ptr - do not use unless you know what you are doing
  !------------------------------------------------------------------------------------------!

  subroutine setProbsPtr(this,ptr)
    ! Input: ptr - new target for the probs
    implicit none
    class(aliasSampler_t) :: this
    real(dp), pointer :: ptr(:)

    this%probs => ptr
  end subroutine setProbsPtr

  !------------------------------------------------------------------------------------------!
  
  subroutine setBiasPtr(this,ptr)
    ! Input: ptr - new target for the bias table
    implicit none
    class(aliasSampler_t) :: this
    real(dp), pointer :: ptr(:)

    call this%table%setBiasTablePtr(ptr)
  end subroutine setBiasPtr

  !------------------------------------------------------------------------------------------!
  
  subroutine setAliasPtr(this,ptr)
    ! Input: ptr - new target for the alias table
    implicit none
    class(aliasSampler_t) :: this
    integer, pointer :: ptr(:)

    call this%table%setAliasTablePtr(ptr)
  end subroutine setAliasPtr

  !------------------------------------------------------------------------------------------!
  ! This is where stuff gets technical...
  ! A sampler array class is required since intel mpi cannot have more than 16381 shared
  ! memory windows (i.e. we could not handle more than ~5500 samplers, which is easily
  ! required for larger systems)
  !------------------------------------------------------------------------------------------!

  subroutine setupSamplerArray(this, nEntries, entrySize)
    implicit none
    class(aliasSamplerArray_t) :: this
    integer(int64), intent(in) :: nEntries, entrySize
    integer(int64) :: totalSize
    integer :: iEntry, windowStart, windowEnd
    real(dp), pointer :: probPtr(:)
    real(dp), pointer :: biasPtr(:)
    integer, pointer :: aliasPtr(:)

    allocate(this%samplerArray(nEntries))

    ! all entries in the array use the same shared memory window, just different
    ! portions of it
    totalSize = nEntries * entrySize
    call safe_shared_memory_alloc(this%allProbsShmw, this%allProbs, totalSize)
    call safe_shared_memory_alloc(this%allBiasShmw, this%allBiasTable, totalSize)
    if(associated(this%allAliasTable)) &
         call shared_deallocate_mpi(this%allAliasShmw, this%allAliasTable)
    call shared_allocate_mpi(this%allAliasShmw, this%allAliasTable, (/totalSize/))
    
    do iEntry = 1, nEntries
       ! from where to where this entry has memory access in the shared resources
       windowStart = (iEntry-1)*entrySize + 1
       windowEnd = windowStart + entrySize - 1

       ! set this entry's pointers
       probPtr => this%allProbs(windowStart:windowEnd)
       call this%samplerArray(iEntry)%setProbsPtr(probPtr)

       aliasPtr => this%allAliasTable(windowStart:windowEnd)
       call this%samplerArray(iEntry)%setAliasPtr(aliasPtr)

       biasPtr => this%allBiasTable(windowStart:windowEnd)
       call this%samplerArray(iEntry)%setBiasPtr(biasPtr)
    end do
    
  end subroutine setupSamplerArray

  !------------------------------------------------------------------------------------------!

  subroutine setupEntry(this, iEntry, arr)
    ! Input: iEntry - index of the entry to initialize
    !        arr - data to be loaded by that entry
    implicit none
    class(aliasSamplerArray_t) :: this
    integer, intent(in) :: iEntry
    real(dp), intent(in) :: arr(:)

    call this%samplerArray(iEntry)%initSampler(arr)
  end subroutine setupEntry

  !------------------------------------------------------------------------------------------!

  subroutine samplerArrayDestructor(this)
    implicit none
    class(aliasSamplerArray_t) :: this

    ! free the collective resources
    call safe_shared_memory_dealloc(this%allBiasShmw, this%allBiasTable)
    call safe_shared_memory_dealloc(this%allProbsShmw, this%allProbs)
    if(associated(this%allAliasTable)) &
         call shared_deallocate_mpi(this%allAliasShmw, this%allAliasTable)
    
    ! TODO: null all entry-pointers
  end subroutine samplerArrayDestructor

  !------------------------------------------------------------------------------------------!
  ! Array access functions
  !------------------------------------------------------------------------------------------!

  subroutine aSample(this, iEntry, tgt, prob)
    ! draw a random element from 1:entrySize with the probabilities listed in this entry's prob
    ! Input: iEntry - index of the sampler to use
    !        tgt - on return, this is a random number in the sampling range of entrySize
    !        prob - on return, the probability of picking tgt
    implicit none
    class(aliasSamplerArray_t), intent(in) :: this
    integer, intent(in) :: iEntry
    integer, intent(out) :: tgt
    real(dp), intent(out) :: prob

    call this%samplerArray(iEntry)%sample(tgt,prob)
  end subroutine aSample

  pure function aGetProb(this, iEntry, tgt) result(prob)
    ! Returns the probability to draw tgt from the sampler with index iEntry
    ! Input: iEntry - index of the sampler to use
    !        tgt - the number for which we request the probability of sampling
    ! Output: prob - the probability of drawing tgt with the sample routine
    implicit none
    class(aliasSamplerArray_t), intent(in) :: this
    integer, intent(in) :: iEntry
    integer, intent(in) :: tgt
    real(dp) :: prob

    prob = this%samplerArray(iEntry)%getProb(tgt)
  end function aGetProb
  
  !------------------------------------------------------------------------------------------!
  ! Auxiliary functions to prevent code duplication
  !------------------------------------------------------------------------------------------!

  subroutine safe_shared_memory_alloc(win,ptr,size)
    ! wrapper for shared_allocate_mpi that tests if the pointer is associated          
    ! Input: win - MPI shared memory window for internal MPI usage
    !        ptr - pointer to be allocated, on return points to a shared memory segment of given size
    !        size - size of the memory segment to be allocated
    implicit none
    integer(MPIArg) :: win
    real(dp), pointer :: ptr(:)
    integer(int64) :: size

    ! if pointer was allocated prior, re-allocate the probabilities
    ! WARNING: DO NOT MANUALLY RE-ASSIGN ptr, THIS WILL MOST LIKELY BREAK STUFF    
    call safe_shared_memory_dealloc(win,ptr)
    call shared_allocate_mpi(win, ptr, (/size/))
  end subroutine safe_shared_memory_alloc

  !------------------------------------------------------------------------------------------!

  subroutine safe_shared_memory_dealloc(win,ptr)
    ! wrapper for shared_deallocate_mpi that tests if the pointer is associated
    ! Input: win - MPI shared memory window for internal MPI usage
    !        ptr - pointer to be deallocated (if associated)
    ! WARNING: THIS ASSUMES THAT IF ptr IS ASSOCIATED, IT POINTS TO AN MPI SHARED MEMORY
    !          WINDOW win
    implicit none
    integer(MPIArg) :: win
    real(dp), pointer :: ptr(:)

    ! assume that if ptr is associated, it points to mpi shared memory
    if(associated(ptr)) call shared_deallocate_mpi(win, ptr)    
  end subroutine safe_shared_memory_dealloc
  
  !------------------------------------------------------------------------------------------!
  ! Public non-member function to deallocate 1d-arrays of samplers (common task)
  !------------------------------------------------------------------------------------------!
  
    subroutine clear_sampler_array(arr)
      ! call the destructor on all elements of an array, then deallocate it
      type(aliasSampler_t), allocatable :: arr(:)

      integer :: i

      do i = 1, size(arr)
         call arr(i)%samplerDestructor()
      end do
      deallocate(arr)
    end subroutine clear_sampler_array
  
end module aliasSampling
