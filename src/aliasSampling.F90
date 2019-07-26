module aliasSampling
  ! This module contains the utility to use alias table lookup on lists,
  ! requiring to precompute biases but making the lookup O(1)
  use constants
  use shared_memory_mpi
  use ParallelHelper, only: iProcIndex_intra
  use dSFMT_interface , only : genrand_real2_dSFMT  
  implicit none

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
     ! destructor - final would be suited better, but is not supported by all compilers
     procedure :: tableDestructor
     ! get a random value from the alias table
     procedure :: getRand
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
     ! destructor - see above re final
     procedure :: samplerDestructor
     ! get a random element and the generation probability
     procedure :: sample
  end type aliasSampler_t
 
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

    integer :: i,j,cV, cU
    integer(int64) :: arrSize
    integer, allocatable :: overfull(:), underfull(:)

    ! allocate the shared memory segment for the alias table
    arrSize = size(arr)

    call safe_shared_memory_alloc(this%biasTableShmw, this%biasTable, arrSize)
    if(associated(this%aliasTable)) &
         call shared_deallocate_mpi(this%aliasTableShmw, this%aliasTable)    
    call shared_allocate_mpi(this%aliasTableShmw, this%aliasTable, (/arrSize/))

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
    
  end subroutine setupTable

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
  ! Initialization / Finalization routines of the sampler
  !------------------------------------------------------------------------------------------!

  subroutine setupSampler(this, arr)
    ! load the probability distribution from arr into this
    implicit none
    class(aliasSampler_t) :: this
    real(dp), intent(in) :: arr(:)

    integer(int64) :: arrSize
    character(*), parameter :: t_r = "setupSampler"
    ! if all weights are 0, throw an error
    if(sum(arr) < eps) then
       call stop_all(t_r, "Trying to initialize sampler with empty probability distribution")
    endif
    
    ! initialize the alias table
    call this%table%setupTable(arr)
    arrSize = size(arr)

    ! allocate the probabilities
    call safe_shared_memory_alloc(this%probsShmw, this%probs, arrSize)
    ! the array is shared memory, so only node-root has to do this
    if(iProcIndex_intra.eq.0) then
       ! the probabilities are taken from input and normalized
       this%probs = arr / sum(arr)
    end if
  end subroutine setupSampler

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
    ! draw a random element from 1:size(prob) with the probabilities listed in prob
    implicit none
    class(aliasSampler_t) :: this
    integer, intent(out) :: tgt
    real(dp), intent(out) :: prob

    ! in debug, do a sanity check: is this initialized?
    if(.not.associated(this%probs)) then
       call stop_all("Sample","Trying to call routines of an unitialized sampler object")
    end if
    ! get the drawn number from the alias table
    tgt = this%table%getRand()
    ! and its probability
    prob = this%probs(tgt)
  end subroutine sample
  
  !------------------------------------------------------------------------------------------!
  ! Auxiliary functions to prevent code duplication
  !------------------------------------------------------------------------------------------!

  ! wrapper for shared_allocate_mpi that tests if the pointer is associated
  subroutine safe_shared_memory_alloc(win,ptr,size)
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

  ! wrapper for shared_deallocate_mpi that tests if the pointer is associated
  subroutine safe_shared_memory_dealloc(win,ptr)
    implicit none
    integer(MPIArg) :: win
    real(dp), pointer :: ptr(:)

    if(associated(ptr)) call shared_deallocate_mpi(win, ptr)    
  end subroutine safe_shared_memory_dealloc
  
end module aliasSampling
