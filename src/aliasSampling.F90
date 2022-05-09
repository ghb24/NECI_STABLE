#include "macros.h"
module aliasSampling
    ! This module contains the utility to use alias table lookup on lists,
    ! requiring to precompute biases but making the lookup O(1)
    use constants
    use util_mod, only: near_zero, stop_all, binary_search_first_ge
    use shared_memory_mpi
    use MPI_wrapper, only: iProcIndex_intra
    use dSFMT_interface, only: genrand_real2_dSFMT
    use shared_array
    implicit none

    private
    public :: aliasSampler_t, AliasSampler_1D_t, AliasSampler_2D_t, AliasSampler_3D_t, &
        clear_sampler_array, aliasTable_t

    ! type for tables: contains everything you need to get a random number
    ! with given biases
    type aliasTable_t
        private
        ! WARNING: DO NOT MANUALLY RE-ASSIGN THESE POINTERS, THIS WILL MOST LIKELY BREAK STUFF
        ! this is the table of bias
        type(shared_array_real_t) :: biasTable
        ! this is the lookup table for the resulting random number
        type(shared_array_int32_t) :: aliasTable

    contains
        private
        ! constructor
        procedure, public:: setupTable
        ! only load the data, without allocation
        procedure :: initTable
        ! destructor - final would be suited better, but is not supported by all compilers
        procedure, public :: tableDestructor
        ! get a random value from the alias table
        procedure, public :: getRand
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
        type(shared_array_real_t) :: probs

    contains
        ! constructor
        procedure :: setupSampler
        ! only load the data, without allocation
        procedure, private :: initSampler
        ! initialize the probabilities
        procedure, private :: initProbs
        ! destructor - see above re final
        procedure :: samplerDestructor
        ! get a random element and the generation probability
        procedure :: sample
        ! get a random element from a constrained set and its normalized! generation probability
        procedure :: constrained_sample
        ! get the probability to produce a given value
        procedure :: getProb
        ! get the probability to draw a given value from a constrained set
        procedure :: constrained_getProb
    end type aliasSampler_t

    !------------------------------------------------------------------------------------------!
    ! sampler array class - required for technical reasons: we need multiple samplers to share
    ! the same shared memory windows, because the number of shared memory windows
    ! is bounded by the number of communicators which cannot exceed 16381 on most
    ! implementations. https://community.intel.com/t5/Intel-oneAPI-HPC-Toolkit/INTEL-MPI-5-0-Bug-in-MPI-3-shared-memory-allocation-MPI-WIN/td-p/1016993
    !------------------------------------------------------------------------------------------!
    type AliasSampler_3D_t
        private
        ! this is an array of aliasSamplers, in the end
        type(aliasSampler_t), allocatable :: samplerArray(:, :, :)

        ! shared resources of the array entries
        type(shared_array_real_t) :: allProbs
        type(shared_array_real_t) :: allBiasTable
        type(shared_array_int32_t) :: allAliasTable
    contains
        ! constructor
        procedure :: shared_alloc => setupSamplerArray_3D
        procedure :: setup_entry => setupEntry_3D
        ! destructor
        procedure :: finalize => samplerArrayDestructor_3D
        ! get a random element and the generation probability from one of the samplers
        procedure :: sample => aSample_3D
        procedure :: constrained_sample => constrained_sample_3D
        procedure :: get_prob => aGetProb_3D
        procedure :: constrained_getProb => constrained_get_prob_3D
    end type AliasSampler_3D_t


    type AliasSampler_2D_t
        private
        type(AliasSampler_3D_t) :: alias_sampler
    contains
        ! constructor
        procedure :: shared_alloc => setupSamplerArray_2D
        procedure :: setup_entry => setupEntry_2D
        ! destructor
        procedure :: finalize => samplerArrayDestructor_2D
        ! get a random element and the generation probability from one of the samplers
        procedure :: sample => aSample_2D
        procedure :: constrained_sample => constrained_sample_2D
        procedure :: get_prob => aGetProb_2D
        procedure :: constrained_getProb => constrained_get_prob_2D
    end type AliasSampler_2D_t


    type AliasSampler_1D_t
        private
        type(AliasSampler_3D_t) :: alias_sampler
    contains
        ! constructor
        procedure :: shared_alloc => setupSamplerArray_1D
        procedure :: setup_entry => setupEntry_1D
        ! destructor
        procedure :: finalize => samplerArrayDestructor_1D
        ! get a random element and the generation probability from one of the samplers
        procedure :: sample => aSample_1D
        procedure :: constrained_sample => constrained_sample_1D
        procedure :: get_prob => aGetProb_1D
        procedure :: constrained_getProb => constrained_get_prob_1D
    end type AliasSampler_1D_t



contains

    !------------------------------------------------------------------------------------------!
    ! Initialization / Finalization routines of the aliasTable
    !------------------------------------------------------------------------------------------!

    !> pseudo-constructor for alias tables
    !> @param[in] arr  array containing the (not necessarily normalized) probabilities we
    !!              want to use for sampling
    subroutine setupTable(this, arr)
        class(aliasTable_t), intent(inout) :: this
        real(dp), intent(in) :: arr(:)

        character(*), parameter :: t_r = "setupTable"
        integer(int64) :: arrSize
        if (sum(arr) < eps) call stop_all(t_r, &
                                          "Trying to setup empty alias table")

        ! allocate the shared memory segment for the alias table
        arrSize = size(arr, kind=int64)

        call this%biasTable%shared_alloc(arrSize)
        call this%aliasTable%shared_alloc(arrSize)

        call this%initTable(arr)

        ! Sync the shared resource between tasks
        call this%biasTable%sync()
        call this%aliasTable%sync()

    end subroutine setupTable

    !------------------------------------------------------------------------------------------!

    !> Set the bias and alias values for each value in range
    !> @param[in] arr - array containing the (not necessarily normalized) probabilities we
    !>              want to use for sampling
    subroutine initTable(this, arr)
        class(aliasTable_t), intent(inout) :: this
        real(dp), intent(in) :: arr(:)

        integer :: i, j, cV, cU
        integer(int64) :: arrSize
        integer, allocatable :: overfull(:), underfull(:)

        arrSize = size(arr, kind=int64)

        ! as this is shared memory, only node-root has to do this
        if (iProcIndex_intra == 0) then
            ! initialize the probabilities
            this%biasTable%ptr = arr / sum(arr) * arrSize

            ! indices of subarrays
            allocate(overfull(arrSize))
            allocate(underfull(arrSize))

            cV = 0
            cU = 0
            do i = 1, int(arrSize)
                call assignLabel(i)
            end do
            ! we now labeled each entry

            ! it is more efficient to start with the largest biases
            ! -> reverse overfull
            overfull(1:cV) = overfull(cV:1:-1)
            do
                if ((cV == 0) .or. (cU == 0)) then
                    exit
                end if
                ! pick one overfull and one underfull index
                i = overfull(cV)
                j = underfull(cU)
                ! set the alias of the underfull to be the other
                this%aliasTable%ptr(j) = i
                ! correct the bias
                this%biasTable%ptr(i) = this%biasTable%ptr(i) + this%biasTable%ptr(j) - 1.0_dp

                ! unmark j
                cU = cU - 1
                ! reassign i based on the new bias
                cV = cV - 1
                call assignLabel(i)
            end do

            ! make sure we do not leave anything unfilled
            call roundTo1(overfull, cV)
            call roundTo1(underfull, cU)

        end if

    contains

        subroutine assignLabel(i)
            integer, intent(in) :: i

            if (this%biasTable%ptr(i) > 1) then
                cV = cV + 1
                overfull(cV) = i
            else
                cU = cU + 1
                underfull(cU) = i
            end if
        end subroutine assignLabel

        subroutine roundTo1(labels, cI)
            integer, intent(in) :: labels(:)
            integer, intent(in) :: cI

            ! if, due to floating point errors, one of the categories is not empty, empty it
            ! (error is negligible then)
            if (cI > 0) then
                do i = 1, cI
                    this%biasTable%ptr(labels(i)) = 1.0_dp
                    this%aliasTable%ptr(labels(i)) = labels(i)
                end do
            end if

        end subroutine roundTo1

    end subroutine initTable
    !------------------------------------------------------------------------------------------!

    !> clear the memory used by the alias table
    subroutine tableDestructor(this)
        class(aliasTable_t), intent(inout) :: this

        call this%aliasTable%shared_dealloc()
        call this%biasTable%shared_dealloc()
    end subroutine tableDestructor

    !------------------------------------------------------------------------------------------!
    ! Sampling function
    !------------------------------------------------------------------------------------------!

    !> Draw a random number from an alias table created with the corresponding probabilities
    !> @return ind  random number between 1 and the size of the array used to create the
    !!               aliasTable object
    function getRand(this) result(ind)
        class(aliasTable_t), intent(in) :: this
        integer :: ind
        real(dp) :: r, bias
        integer :: sizeArr, pos

        sizeArr = size(this%biasTable%ptr)
        ! random number between 0 and 1
        r = genrand_real2_dSFMT()
        ! random position in arr
        pos = int(sizeArr * r) + 1
        ! remainder of the integer conversion
        ! floating point errors can lead to very small negative values of bias here
        ! this would allow for picking elements which have probability 0 (-> biasTable entry 0)
        ! -> ensure that bias>=0
        bias = max(sizeArr * r + 1.0_dp - real(pos, dp), 0.0_dp)

        if (bias < this%biasTable%ptr(pos)) then
            ind = pos
        else
            ind = this%aliasTable%ptr(pos)
        end if
    end function getRand

    !------------------------------------------------------------------------------------------!
    ! Initialization / Finalization routines of the sampler
    !------------------------------------------------------------------------------------------!

    !> allocate the resources of this and load the probability distribution from arr into this
    !> @param[in] arr  array containing the (not necessarily normalized) probabilities we
    !!              want to use for sampling
    subroutine setupSampler(this, arr)
        class(aliasSampler_t), intent(inout) :: this
        real(dp), intent(in) :: arr(:)

        integer(int64) :: arrSize
        ! if all weights are 0, throw an error
        if (sum(arr) < eps) then
            ! probs defaults to null(), so it is not associated at this point (i.e. in a well-defined state)
            return
        end if

        ! initialize the alias table
        call this%table%setupTable(arr)
        arrSize = size(arr, kind=int64)

        ! allocate the probabilities
        call this%probs%shared_alloc(arrSize)

        ! set the probabilities
        call this%initProbs(arr)

        call this%probs%sync()
    end subroutine setupSampler

    !------------------------------------------------------------------------------------------!

    !> load the probability distribution from arr into this
    !! we only use this in the sampler array, but fortran has no friend classes, so its public
    !> @param[in] arr  array containing the (not necessarily normalized) probabilities we
    !!              want to use for sampling
    subroutine initSampler(this, arr)
        class(aliasSampler_t), intent(inout) :: this
        real(dp), intent(in) :: arr(:)

        ! if all weights are 0, throw an error
        if (sum(arr) < eps) then
            ! if we reach this point, probs is uninitialized -> null it
            this%probs%ptr => null()
            return
        end if

        ! load the data - assume this is pre-allocated
        call this%table%initTable(arr)
        call this%initProbs(arr)
    end subroutine initSampler

    !------------------------------------------------------------------------------------------!

    !> load the probability distribution from arr into this%probs
    !> @param[in] arr  array containing the (not necessarily normalized) probabilities we
    !!              want to use for sampling
    subroutine initProbs(this, arr)
        class(aliasSampler_t), intent(inout) :: this
        real(dp), intent(in) :: arr(:)

        ! the array is shared memory, so only node-root has to do this
        if (iProcIndex_intra == 0) then
            ! the probabilities are taken from input and normalized
            this%probs%ptr = arr / sum(arr)
        end if
    end subroutine initProbs

    !------------------------------------------------------------------------------------------!

    subroutine samplerDestructor(this)
        class(aliasSampler_t), intent(inout) :: this

        ! free the stored probabilities
        call this%probs%shared_dealloc()
        ! free the corresponding alias table
        call this%table%tableDestructor()
    end subroutine samplerDestructor

    !------------------------------------------------------------------------------------------!

    !> draw a random element from 1:size(this%probs) with the probabilities listed in prob
    !> @param[in] tgt  on return, this is a random number in the sampling range of this
    !> @param[out] prob  on return, the probability of picking tgt
    subroutine sample(this, tgt, prob)
        class(aliasSampler_t), intent(in) :: this
        integer, intent(out) :: tgt
        real(dp), intent(out) :: prob

        ! empty samplers don't return anything - since probs defaults to null(), this check is safe
        if (.not. associated(this%probs%ptr)) then
            tgt = 0
            prob = 1.0
            return
        end if
        ! get the drawn number from the alias table
        tgt = this%table%getRand()
        ! and its probability
        prob = this%probs%ptr(tgt)
    end subroutine sample

    !------------------------------------------------------------------------------------------!

    !> draw a random element from 1:size(this%probs) with the probabilities listed in prob
    !> @param[in] constraint pick only elements from constraint
    !> @param[out] tgt  on return, this is a random number in the sampling range of this
    !> @param[out] prob  on return, the probability of picking tgt from constraint
    subroutine constrained_sample(this, constraint, tgt, prob)
        class(aliasSampler_t), intent(in) :: this
        integer, intent(in) :: constraint(:)
        integer, intent(out) :: tgt
        real(dp), intent(out) :: prob
        character(*), parameter :: this_routine = 'constrained_sample'

        real(dp) :: renormalization
        ASSERT(1 <= minval(constraint) .and. maxval(constraint) <= size(this%probs%ptr))
        ! TODO Assert that it is unique

        if (.not. associated(this%probs%ptr)) then
            renormalization = 0.0
        else
            renormalization = sum(this%probs%ptr(constraint))
        end if

        if (near_zero(renormalization)) then
            tgt = 0
            prob = 1.0
        else
            tgt = this%table%getRand()
            do while (binary_search_first_ge(constraint, tgt) == -1)
                tgt = this%table%getRand()
            end do
            prob = this%probs%ptr(tgt) / renormalization
        end if
    end subroutine constrained_sample

    !------------------------------------------------------------------------------------------!

    !> Returns the probability to draw tgt from this sampler
    !> @param[in] tgt  the number for which we request the probability of sampling
    !> @param[out] prob  the probability of drawing tgt with the sample routine
    elemental function getProb(this, tgt) result(prob)
        class(aliasSampler_t), intent(in) :: this
        integer, intent(in) :: tgt
        real(dp) :: prob

        ! the probability of drawing anything from an empty sampler is 0
        if (.not. associated(this%probs%ptr)) then
            prob = 0.0
        else
            prob = this%probs%ptr(tgt)
        end if
    end function getProb

    !------------------------------------------------------------------------------------------!

    !> Returns the probability to draw tgt from this sampler
    !> @param[in] tgt  the number for which we request the probability of sampling
    !> @param[in] constraint pick only elements from constraint
    !> @param[out] prob  the probability of picking tgt from constraint
    pure function constrained_getProb(this, constraint, tgt) result(prob)
        class(aliasSampler_t), intent(in) :: this
        integer, intent(in) :: constraint(:)
        integer, intent(in) :: tgt
        character(*), parameter :: this_routine = 'constrained_getProb'
        real(dp) :: prob

        real(dp) :: renormalization
        ASSERT(1 <= minval(constraint) .and. maxval(constraint) <= size(this%probs%ptr))

        renormalization = sum(this%getProb(constraint))
        ! the probability of drawing anything from an empty sampler is 0
        if (near_zero(renormalization)) then
            prob = 0.0
        else
            prob = this%probs%ptr(tgt) / renormalization
        end if
    end function constrained_getProb

    !------------------------------------------------------------------------------------------!

    !------------------------------------------------------------------------------------------!
    ! Public non-member function to deallocate 1d-arrays of samplers (common task)
    !------------------------------------------------------------------------------------------!

    !> call the destructor on all elements of an array, then deallocate it. This is for
    !! intrinsic arrays, the sampler array class has its own deallocate routine.
    !> @param[in, out] arr  array to deallocate
    subroutine clear_sampler_array(arr)
        type(aliasSampler_t), allocatable, intent(inout) :: arr(:)

        integer :: i

        do i = 1, size(arr)
            call arr(i)%samplerDestructor()
        end do
        deallocate(arr)
    end subroutine clear_sampler_array


    !------------------------------------------------------------------------------------------!
    ! This is where stuff gets technical...
    ! A sampler array class is required since intel mpi cannot have more than 16381 shared
    ! memory windows (i.e. we could not handle more than ~5500 samplers, which is easily
    ! required for larger systems)
    !------------------------------------------------------------------------------------------!

    !> Setup an array of samplers using a single shared resource (split into parts associated
    !! with one of them each). This only does the allocation.
    !> @param[in] nEntries  number of samplers to initialise
    !> @param[in] entrySize  number of values per sampler
    subroutine setupSamplerArray_1D(this, nEntries, entrySize, name)
        class(AliasSampler_1D_t) :: this
        integer, intent(in) :: nEntries, entrySize
        character(*), intent(in) :: name
        call this%alias_sampler%shared_alloc([nEntries, 1, 1], entrySize, name)
    end subroutine setupSamplerArray_1D

    !------------------------------------------------------------------------------------------!

    !> Initialise one sampler of an array
    !> @param[in] iEntry  index of the entry to initialize
    !> @param[in] arr  data to be loaded by that entry
    subroutine setupEntry_1D(this, iEntry, arr)
        class(AliasSampler_1D_t), intent(inout) :: this
        integer, intent(in) :: iEntry
        real(dp), intent(in) :: arr(:)
        call this%alias_sampler%setup_entry(iEntry, 1, 1, arr)
    end subroutine setupEntry_1D

    !------------------------------------------------------------------------------------------!

    !> Deallocate an array of samplers
    subroutine samplerArrayDestructor_1D(this)
        class(AliasSampler_1D_t), intent(inout) :: this
        call this%alias_sampler%finalize()
    end subroutine samplerArrayDestructor_1D

    !------------------------------------------------------------------------------------------!
    ! Array access functions
    !------------------------------------------------------------------------------------------!

    !> draw a random element from 1:entrySize with the probabilities listed in this entry's prob
    !> @param[in] iEntry  index of the sampler to use
    !> @param[out] tgt  on return, this is a random number in the sampling range of entrySize
    !> @param[out] prob  on return, the probability of picking tgt
    subroutine aSample_1D(this, iEntry, tgt, prob)
        class(AliasSampler_1D_t), intent(in) :: this
        integer, intent(in) :: iEntry
        integer, intent(out) :: tgt
        real(dp), intent(out) :: prob
        call this%alias_sampler%sample(iEntry, 1, 1, tgt, prob)
    end subroutine aSample_1D

    !> draw a random element from 1:entrySize with the probabilities listed in this entry's prob
    !> @param[in] i Index of the sampler to use
    !> @param[in] constraint pick only elements from constraint
    !> @param[out] tgt  on return, this is a random number in the sampling range of entrySize and in constraint
    !> @param[out] prob  on return, the probability of picking tgt from constraint
    subroutine constrained_sample_1D(this, i, constraint, tgt, prob)
        class(AliasSampler_1D_t), intent(in) :: this
        integer, intent(in) :: i
        integer, intent(in) :: constraint(:)
        integer, intent(out) :: tgt
        real(dp), intent(out) :: prob

        call this%alias_sampler%constrained_sample(i, 1, 1, constraint, tgt, prob)
    end subroutine

    !> Returns the probability to draw tgt from the sampler with index iEntry
    !> @param[in] iEntry  index of the sampler to use
    !> @param[in] tgt  the number for which we request the probability of sampling
    !> @return prob  the probability of drawing tgt with the sample routine
    pure function aGetProb_1D(this, iEntry, tgt) result(prob)
        class(AliasSampler_1D_t), intent(in) :: this
        integer, intent(in) :: iEntry
        integer, intent(in) :: tgt
        real(dp) :: prob
        prob = this%alias_sampler%get_prob(iEntry, 1, 1, tgt)
    end function aGetProb_1D

    !> Returns the probability to draw tgt from the sampler with index iEntry
    !> @param[in] i Index of the sampler to use
    !> @param[in] constraint pick only elements from constraint
    !> @param[in] tgt  the number for which we request the probability of sampling
    !> @return prob  the probability of drawing tgt with the sample routine from constraint
    pure function constrained_get_prob_1D(this, i, constraint, tgt) result(prob)
        class(AliasSampler_1D_t), intent(in) :: this
        integer, intent(in) :: i
        integer, intent(in) :: constraint(:)
        integer, intent(in) :: tgt
        real(dp) :: prob

        prob = this%alias_sampler%constrained_getProb(i, 1, 1, constraint, tgt)
    end function

    !> Setup an array of samplers using a single shared resource (split into parts associated
    !! with one of them each). This only does the allocation.
    !> @param[in] nEntries  number of samplers to initialise
    !> @param[in] entrySize  number of values per sampler
    subroutine setupSamplerArray_2D(this, dims, entrySize, name)
        class(AliasSampler_2D_t) :: this
        integer, intent(in) :: dims(2), entrySize
        character(*), intent(in) :: name
        call this%alias_sampler%shared_alloc([dims(1), dims(2), 1], entrySize, name)
    end subroutine setupSamplerArray_2D

    !------------------------------------------------------------------------------------------!

    !> Initialise one sampler of an array
    !> @param[in] iEntry  index of the entry to initialize
    !> @param[in] arr  data to be loaded by that entry
    subroutine setupEntry_2D(this, i, j, arr)
        class(AliasSampler_2D_t), intent(inout) :: this
        integer, intent(in) :: i, j
        real(dp), intent(in) :: arr(:)
        call this%alias_sampler%setup_entry(i, j, 1, arr)
    end subroutine setupEntry_2D

    !------------------------------------------------------------------------------------------!

    !> Deallocate an array of samplers
    subroutine samplerArrayDestructor_2D(this)
        class(AliasSampler_2D_t), intent(inout) :: this
        call this%alias_sampler%finalize()
    end subroutine samplerArrayDestructor_2D

    !------------------------------------------------------------------------------------------!
    ! Array access functions
    !------------------------------------------------------------------------------------------!

    !> draw a random element from 1:entrySize with the probabilities listed in this entry's prob
    !> @param[in] iEntry  index of the sampler to use
    !> @param[out] tgt  on return, this is a random number in the sampling range of entrySize
    !> @param[out] prob  on return, the probability of picking tgt
    subroutine aSample_2D(this, i, j, tgt, prob)
        class(AliasSampler_2D_t), intent(in) :: this
        integer, intent(in) :: i, j
        integer, intent(out) :: tgt
        real(dp), intent(out) :: prob
        call this%alias_sampler%sample(i, j, 1, tgt, prob)
    end subroutine aSample_2D

    !> draw a random element from 1:entrySize with the probabilities listed in this entry's prob
    !> @param[in] i Index of the sampler to use
    !> @param[in] j Index of the sampler to use
    !> @param[in] constraint pick only elements from constraint
    !> @param[out] tgt  on return, this is a random number in the sampling range of entrySize and in constraint
    !> @param[out] prob  on return, the probability of picking tgt from constraint
    subroutine constrained_sample_2D(this, i, j, constraint, tgt, prob)
        class(AliasSampler_2D_t), intent(in) :: this
        integer, intent(in) :: i, j
        integer, intent(in) :: constraint(:)
        integer, intent(out) :: tgt
        real(dp), intent(out) :: prob

        call this%alias_sampler%constrained_sample(i, j, 1, constraint, tgt, prob)
    end subroutine

    !> Returns the probability to draw tgt from the sampler with index iEntry
    !> @param[in] iEntry  index of the sampler to use
    !> @param[in] tgt  the number for which we request the probability of sampling
    !> @return prob  the probability of drawing tgt with the sample routine
    pure function aGetProb_2D(this, i, j, tgt) result(prob)
        class(AliasSampler_2D_t), intent(in) :: this
        integer, intent(in) :: i, j
        integer, intent(in) :: tgt
        real(dp) :: prob
        prob = this%alias_sampler%get_prob(i, j, 1, tgt)
    end function aGetProb_2D

    !> Returns the probability to draw tgt from the sampler with index iEntry
    !> @param[in] i Index of the sampler to use
    !> @param[in] j Index of the sampler to use
    !> @param[in] constraint pick only elements from constraint
    !> @param[in] tgt  the number for which we request the probability of sampling
    !> @return prob  the probability of drawing tgt with the sample routine from constraint
    pure function constrained_get_prob_2D(this, i, j, constraint, tgt) result(prob)
        class(AliasSampler_2D_t), intent(in) :: this
        integer, intent(in) :: i, j
        integer, intent(in) :: constraint(:)
        integer, intent(in) :: tgt
        real(dp) :: prob

        prob = this%alias_sampler%constrained_getProb(i, j, 1, constraint, tgt)
    end function


    !> Setup an array of samplers using a single shared resource (split into parts associated
    !! with one of them each). This only does the allocation.
    !> @param[in] dims Dimension of the three-dimensional array of samplers.
    !> @param[in] entry_size number of values per sampler
    subroutine setupSamplerArray_3D(this, dims, entry_size, name)
        class(AliasSampler_3D_t), intent(inout) :: this
        ! NOTE: We might have to change dims and entry_size to int64 in the near future... :-(
        integer, intent(in) :: dims(3), entry_size
        character(*), intent(in) :: name

        integer :: i, j, k
        integer(int64) :: window_start, window_end, total_size

        allocate(this%samplerArray(dims(1), dims(2), dims(3)))

        ! all entries in the array use the same shared memory window, just different
        ! portions of it
        total_size = entry_size * product(int(dims, kind=int64))
        call this%allProbs%shared_alloc(total_size, name//'_Probs')
        call this%allBiasTable%shared_alloc(total_size, name//'_Bias')
        call this%allAliasTable%shared_alloc(total_size, name//'_Alias')

        window_start = 1
        do k = 1, dims(3)
            do j = 1, dims(2)
                do i = 1, dims(1)
                    ! from where to where this entry has memory access in the shared resources
                    window_end = window_start + entry_size - 1

                    ! set this entry's pointers
                    this%samplerArray(i, j, k)%probs%ptr => this%allProbs%ptr(window_start:window_end)
                    this%samplerArray(i, j, k)%table%aliasTable%ptr => this%allAliasTable%ptr(window_start:window_end)
                    this%samplerArray(i, j, k)%table%biasTable%ptr => this%allBiasTable%ptr(window_start:window_end)
                    window_start = window_end + 1
                end do
            end do
        end do
    end subroutine

    !> @brief
    !> Initialise one sampler of an array
    !>
    !> @param[in] i Index of the entry to initialize
    !> @param[in] j Index of the entry to initialize
    !> @param[in] k Index of the entry to initialize
    !> @param[in] arr  data to be loaded by that entry
    subroutine setupEntry_3D(this, i, j, k, arr)
        class(AliasSampler_3D_t), intent(inout) :: this
        integer, intent(in) :: i, j, k
        real(dp), intent(in) :: arr(:)

        call this%samplerArray(i, j, k)%initSampler(arr)

        ! Sync the shared resources
        call this%allBiasTable%sync()
        call this%allAliasTable%sync()
        call this%allProbs%sync()
    end subroutine


    !> @brief
    !> Deallocate an array of samplers
    subroutine samplerArrayDestructor_3D(this)
        class(AliasSampler_3D_t), intent(inout) :: this

        ! free the collective resources
        call this%allAliasTable%shared_dealloc()
        call this%allProbs%shared_dealloc()
        call this%allBiasTable%shared_dealloc()

        if (allocated(this%samplerArray)) deallocate(this%samplerArray)
    end subroutine


    !------------------------------------------------------------------------------------------!
    ! Array access functions
    !------------------------------------------------------------------------------------------!

    !> draw a random element from 1:entrySize with the probabilities listed in this entry's prob
    !> @param[in] i Index of the sampler to use
    !> @param[in] j Index of the sampler to use
    !> @param[in] k Index of the sampler to use
    !> @param[out] tgt  on return, this is a random number in the sampling range of entrySize
    !> @param[out] prob  on return, the probability of picking tgt
    subroutine aSample_3D(this, i, j, k, tgt, prob)
        class(AliasSampler_3D_t), intent(in) :: this
        integer, intent(in) :: i, j, k
        integer, intent(out) :: tgt
        real(dp), intent(out) :: prob

        call this%samplerArray(i, j, k)%sample(tgt, prob)
    end subroutine


    !> draw a random element from 1:entrySize with the probabilities listed in this entry's prob
    !> @param[in] i Index of the sampler to use
    !> @param[in] j Index of the sampler to use
    !> @param[in] k Index of the sampler to use
    !> @param[in] constraint pick only elements from constraint
    !> @param[out] tgt  on return, this is a random number in the sampling range of entrySize and in constraint
    !> @param[out] prob  on return, the probability of picking tgt from constraint
    subroutine constrained_sample_3D(this, i, j, k, constraint, tgt, prob)
        class(AliasSampler_3D_t), intent(in) :: this
        integer, intent(in) :: i, j, k
        integer, intent(in) :: constraint(:)
        integer, intent(out) :: tgt
        real(dp), intent(out) :: prob

        call this%samplerArray(i, j, k)%constrained_sample(constraint, tgt, prob)
    end subroutine

    !> Returns the probability to draw tgt from the sampler with index iEntry
    !> @param[in] i Index of the sampler to use
    !> @param[in] j Index of the sampler to use
    !> @param[in] k Index of the sampler to use
    !> @param[in] tgt  the number for which we request the probability of sampling
    !> @return prob  the probability of drawing tgt with the sample routine
    pure function aGetProb_3D(this, i, j, k, tgt) result(prob)
        class(AliasSampler_3D_t), intent(in) :: this
        integer, intent(in) :: i, j, k
        integer, intent(in) :: tgt
        real(dp) :: prob

        prob = this%samplerArray(i, j, k)%getProb(tgt)
    end function


    !> Returns the probability to draw tgt from the sampler with index iEntry
    !> @param[in] i Index of the sampler to use
    !> @param[in] j Index of the sampler to use
    !> @param[in] k Index of the sampler to use
    !> @param[in] constraint pick only elements from constraint
    !> @param[in] tgt  the number for which we request the probability of sampling
    !> @return prob  the probability of drawing tgt with the sample routine from constraint
    pure function constrained_get_prob_3D(this, i, j, k, constraint, tgt) result(prob)
        class(AliasSampler_3D_t), intent(in) :: this
        integer, intent(in) :: i, j, k
        integer, intent(in) :: constraint(:)
        integer, intent(in) :: tgt
        real(dp) :: prob

        prob = this%samplerArray(i, j, k)%constrained_getProb(constraint, tgt)
    end function

end module aliasSampling
