#include "macros.h"

module LMat_class
    use LMat_indexing, only: lMatIndSym
    use shared_array
    use SystemData, only: nBasis
    use ParallelHelper
    use Parallel_neci
    use procedure_pointers, only: lMatInd_t
    use constants
    use shared_rhash, only: shared_rhash_t
    use mpi
    use util_mod, only: get_free_unit, operator(.div.)
    use tc_three_body_data, only: lMatEps, tHDF5LMat
    use HElem, only: HElement_t_sizeB
    use UMatCache, only: numBasisIndices
    use LoggingData, only: tHistLMat    
#ifdef USE_HDF_
    use hdf5
    use hdf5_util
#endif
    implicit none
    private
    public :: lMat_t, sparse_lMat_t, dense_lMat_t

    !> Abstract base class for lMat_t objects (6-index integrals)
    type, abstract :: lMat_t
        private
        
        ! Generic indexing routine
        procedure(lMatInd_t), nopass, pointer, public :: indexFunc => lMatIndSym        
    contains
        ! Element getters/setters
        procedure(get_elem_t), deferred :: get_elem
        procedure(set_elem_t), deferred, private :: set_elem

        ! Allocation routines
        procedure :: alloc => void_one
        procedure :: dealloc => void_zero

        ! I/O routines
        ! Interfaced read routine: Delegates to the read_kernel
        procedure :: read
        ! Routine to read in any format (ascii/hdf5) to this object
        procedure(read_t), deferred, private :: read_kernel
        ! When reading hdf5, the read is done by a lMat_hdf5_read object. This object
        ! calls a read_op_t from the respective lMat to load the data to memory
        ! This has to be threadsafe
        procedure(read_op_t), deferred, private :: read_op_hdf5

        procedure :: lMat_size
        procedure :: histogram_lMat
    end type lMat_t

    !------------------------------------------------------------------------------------------!

    !> Implementation for densely stored 6-index objects
    type, extends(lMat_t) :: dense_lMat_t
        private
        ! The values of the integrals
#ifdef CMPLX_
        type(shared_array_cmplx_t) :: lMat_vals
#else
        type(shared_array_real_t) :: lMat_vals
#endif
    contains
        ! Element getters/setters
        procedure :: get_elem => get_elem_dense
        procedure, private :: set_elem => set_elem_dense

        ! Allocation routines 
        procedure :: alloc => alloc_dense
        procedure :: dealloc => dealloc_dense

        ! I/O routines
        procedure, private :: read_kernel => read_dense
        procedure, private :: read_hdf5_dense
        procedure, private :: read_op_hdf5 => read_op_dense_hdf5
    end type dense_lMat_t

    !------------------------------------------------------------------------------------------!    

    !> Implementation for sparsely stored 6-index objects
    type, extends(lMat_t) :: sparse_lMat_t
        private
        ! The values of the nonzero integrals
#ifdef CMPLX_
        type(shared_array_cmplx_t) :: nonzero_vals
#else
        type(shared_array_real_t) :: nonzero_vals
#endif
        ! read-only shared memory hash table
        type(shared_rhash_t) :: htable
    contains
        ! Element getters/setters
        procedure :: get_elem => get_elem_sparse
        procedure, private :: set_elem => set_elem_sparse

        ! Allocation routines
        procedure :: alloc => alloc_sparse
        procedure :: dealloc => dealloc_sparse

        ! I/O routines
        procedure, private :: read_kernel => read_sparse
        procedure, private :: read_op_hdf5 => read_op_sparse

        ! These are auxiliary internal I/O routines
        procedure, private :: read_data
        procedure, private :: count_conflicts        
    end type sparse_lMat_t

    !------------------------------------------------------------------------------------------!

    !> Handler for reading hdf5 tcdump files. Calls the read_op_hdf5 of the calling lMat_t
    ! This cannot go into a separate module since this would introduce a cyclical dependency.
    ! Every way of lifting this cyclical dependency would couple different modules tightly
    ! and introduce logic belonging to the lMat objects and their components to other, unrelated
    ! objects. Hence, this is contained here.
#ifdef USE_HDF_
    type :: lMat_hdf5_read_t
        private
        integer(hid_t) :: err, file_id, plist_id, grp_id, ds_vals, ds_inds
        integer(hsize_t), allocatable :: offsets(:)
        integer(hsize_t) :: countsEnd
    contains
        procedure :: open
        procedure :: close
        procedure :: loop_file
    end type lMat_hdf5_read_t
#endif
    
    ! Names of the LMat hdf5 datasets
    character(*), parameter :: nm_grp = "tcdump", nm_nInts = "nInts", nm_vals = "values", nm_indices = "indices"    

    !------------------------------------------------------------------------------------------!    

    ! Interfaces for deferred functions
    abstract interface
        !> Set an element of a lMat
        subroutine set_elem_t(this, index, element)
            use constants            
            import :: lMat_t
            class(lMat_t), intent(inout) :: this
            integer(int64), intent(in) :: index
            HElement_t(dp), intent(in) :: element
        end subroutine set_elem_t

        !> Get an element of a lMat. This replaces the old lMatAccess function pointer
        function get_elem_t(this, index) result(element)
            use constants            
            import :: lMat_t
            class(lMat_t), intent(in) :: this
            integer(int64), intent(in) :: index
            HElement_t(dp) :: element
        end function get_elem_t

        !> Read a lMat from a file
        subroutine read_t(this, filename)
            import :: lMat_t
            class(lMat_t), intent(inout) :: this
            character(*), intent(in) :: filename
        end subroutine read_t

        !> Read operation on a single block of data read from an hdf5 file
        subroutine read_op_t(this, indices, entries)
            use constants            
            import :: lMat_t
            class(lMat_t), intent(inout) :: this
            ! The read operation is allowed to deallocate the input to make
            ! memory available
            integer(int64), allocatable, intent(inout) :: indices(:,:), entries(:,:)
        end subroutine read_op_t
    end interface

contains

    !------------------------------------------------------------------------------------------!
    ! Generic routines
    !------------------------------------------------------------------------------------------!    

    !> Empty routine: This is the default operation without arguments
    subroutine void_zero(this)
        class(lMat_t), intent(inout) :: this

        unused_var(this)
    end subroutine void_zero

    !> Empty routine: This is the default operation with one int64 argument
    !> @param[in] size input argument (no operations performed)
    subroutine void_one(this, size)
        class(lMat_t), intent(inout) :: this
        integer(int64), intent(in) :: size

        unused_var(this)
        unused_var(size)
    end subroutine void_one

    !------------------------------------------------------------------------------------------!    

    !> Return the max. index appearing in this lMat_t (i.e. the number of 6-index integrals)
    !> @return size The number of 6-index integrals of this object, depending on the symmetry.
    function lMat_size(this) result(size)
        class(lMat_t), intent(in) :: this
        integer(int64) :: size

        integer(int64) :: nBI

        nBI = int(numBasisIndices(nBasis),int64)
        ! The size is given by the largest index (LMatInd is monotonous in all arguments)        
        size = this%indexFunc(nBI,nBI,nBI,nBI,nBI,nBI)
    end function lMat_size

    !------------------------------------------------------------------------------------------!

    !> Read in the 6-index integrals from disk and histogram the integrals. The file itself
    !! only has to store the nonzero integrals, either in ASCII or in HDF5 format. For conventions
    !! in the HDF5 format, please refer to the developer's guide.
    !> @param[in] filename name of the integrals file
    subroutine read(this, filename)
        class(lMat_t), intent(inout) :: this
        character(*), intent(in) :: filename

        call this%read_kernel(filename)

        if(tHistLMat) call this%histogram_lMat()

    end subroutine read

    !------------------------------------------------------------------------------------------!
    ! Dense lMat routines
    !------------------------------------------------------------------------------------------!   

    !> Allocate the 6-index integrals for the dense storage
    !> @param[in] size size of the integral container to be allocated
    subroutine alloc_dense(this, size)
        class(dense_lMat_t), intent(inout) :: this
        integer(int64), intent(in) :: size

        write(iout,*) "Six-index integrals require", real(size)*real(HElement_t_SizeB)/(2.0**30), "GB"
        call this%lmat_vals%shared_alloc(size, "LMat")
        if(iProcIndex_intra == 0) then
            this%lmat_vals%ptr = 0.0_dp
        end if        
    end subroutine alloc_dense

    !------------------------------------------------------------------------------------------!    

    !> Deallocate the 6-index integrals (dense)
    subroutine dealloc_dense(this)
        class(dense_lMat_t), intent(inout) :: this

        call this%lMat_vals%shared_dealloc()
    end subroutine dealloc_dense

    !------------------------------------------------------------------------------------------!    

    !> Set an element in the dense 6-index integrals to a new value
    !> @param[in] index  position of the element
    !> @param[in] element  new value of the element
    subroutine set_elem_dense(this, index, element)
        class(dense_lMat_t), intent(inout) :: this
        integer(int64), intent(in) :: index
        HElement_t(dp), intent(in) :: element

        this%lMat_vals%ptr(index) = element
    end subroutine set_elem_dense

    !------------------------------------------------------------------------------------------!   

    !> Get an element of the 6-index integrals from the densely stored container
    !> @param[in] index  position of the element
    !> @return element  value of the element
    function get_elem_dense(this, index) result(element)
        class(dense_lMat_t), intent(in) :: this
        integer(int64), intent(in) :: index
        HElement_t(dp) :: element

        element = this%lMat_vals%ptr(index)
    end function get_elem_dense

    !------------------------------------------------------------------------------------------!

    !> Read the 6-index integrals from a file to dense format
    !> @param[in] filename  name of the file to read from
    subroutine read_dense(this, filename)
        class(dense_lMat_t), intent(inout) :: this
        character(*), intent(in) :: filename    
        integer :: iunit, ierr
        integer(int64) :: a,b,c,i,j,k
        HElement_t(dp) :: matel
        character(*), parameter :: t_r = "readLMat"
        integer(int64) :: counter

        call this%alloc(this%lMat_size())

        if(tHDF5LMat) then
#ifdef USE_HDF_
            call this%read_hdf5_dense(filename)
#else
            call stop_all(t_r, "HDF5 integral files disabled at compile time")
#endif
        else
            if(iProcIndex_intra .eq. 0) then
                iunit = get_free_unit()
                open(iunit,file = filename, status = 'old')
                counter = 0
                do
                    read(iunit,*,iostat = ierr) matel, a,b,c,i,j,k
                    ! end of file reached?
                    if(ierr < 0) then
                        exit
                    else if(ierr > 0) then
                        ! error while reading?
                        call stop_all(t_r,"Error reading TCDUMP file")
                    else
                        ! else assign the matrix element
                        if(this%indexFunc(a,b,c,i,j,k) > this%lMat_size()) then
                            counter = this%indexFunc(a,b,c,i,j,k)
                            write(iout,*) "Warning, exceeding size"
                        endif
                        if(abs(3.0_dp*matel) > LMatEps) &
                            this%lMat_vals%ptr(this%indexFunc(a,b,c,i,j,k)) = 3.0_dp * matel
                        if(abs(matel)> 0.0_dp) counter = counter + 1
                    endif

                end do

                counter = counter / 12

                write(iout, *) "Sparsity of LMat", real(counter)/real(this%lMat_size())
                write(iout, *) "Nonzero elements in LMat", counter
                write(iout, *) "Allocated size of LMat", this%lMat_size()
            endif
            call MPIBcast(counter)
        end if
    end subroutine read_dense

    !------------------------------------------------------------------------------------------!    

    !> Read the integrals from an hdf5 file to dense format
    !> @param[in] filename  name of the file to read from
    subroutine read_hdf5_dense(this, filename)
        class(dense_lMat_t), intent(inout) :: this
        character(*), intent(in) :: filename
#ifdef USE_HDF_        
        type(lMat_hdf5_read_t) :: reader
        integer(hsize_t) :: nInts

        call reader%open(filename, nInts)
        call reader%loop_file(this)
        call reader%close()
#else
        character(*), parameter :: t_r = "read_hdf5_dense"
        unused_var(this)
        if(len(filename) /= 0) continue
        call stop_all(t_r, "hdf5 support disabled at compile time")
#endif
    end subroutine read_hdf5_dense

    !------------------------------------------------------------------------------------------!
    ! Read/count operations
    !------------------------------------------------------------------------------------------!

    !> This is the operation to be performed on each block of data read from an hdf5 file
    !! both arguments may or may not be still allocated upon return
    !> @param[in,out] indices  chunk of indices read in from the file
    !> @param[in,out] entries  chunk of corresponding values
    subroutine read_op_dense_hdf5(this, indices, entries)
        class(dense_lMat_t), intent(inout) :: this
        integer(int64), allocatable, intent(inout) :: indices(:,:), entries(:,:)
        integer(int64) :: i, this_blocksize
        HElement_t(dp) :: rVal

        this_blocksize = size(entries, dim=2)        

        do i = 1, this_blocksize
            ! truncate down to lMatEps
            rVal = 3.0_dp * transfer(entries(1,i),rVal)
            if(abs(rVal)>lMatEps) then
                call this%set_elem(this%indexFunc(int(indices(1,i),int64),int(indices(2,i),int64),&
                    int(indices(3,i),int64),&
                    int(indices(4,i),int64),int(indices(5,i),int64),int(indices(6,i),int64)) &
                    , rVal)
            endif
        end do
    end subroutine read_op_dense_hdf5

    !------------------------------------------------------------------------------------------!
    ! Sparse lMat functions
    !------------------------------------------------------------------------------------------!    

    !> Allocate memory for the sparse storage of the 6-index integrals
    !> @param[in] size  number of non-zero integrals
    subroutine alloc_sparse(this, size)
        class(sparse_lMat_t), intent(inout) :: this
        integer(int64), intent(in) :: size
        character(*), parameter :: t_r = "alloc_sparse"

        write(iout,*) "Six-index integrals require", real(size)*real(HElement_t_SizeB)/(2.0**30), "GB"
        call this%nonzero_vals%shared_alloc(size, "LMat")
        ! For now, have the htable of the same size as the integrals
        call this%htable%alloc(size,size)
        write(iout,*) "Sparse format overhead is", 2*real(size)*real(sizeof_int64)/(2.0**30), "GB"
    end subroutine alloc_sparse

    !------------------------------------------------------------------------------------------!    

    !> Deallocate memory used for the sparse storage of the 6-index integrals
    subroutine dealloc_sparse(this)
        class(sparse_lMat_t), intent(inout) :: this
        character(*), parameter :: t_r = "dealloc_sparse"

        ! Requires deallocation of the values and the hash table for the indices
        call this%nonzero_vals%shared_dealloc()
        call this%htable%dealloc()
    end subroutine dealloc_sparse

    !------------------------------------------------------------------------------------------!    

    !> Set an element to the sparsely stored 6-index integrals. This requires the hash
    !! table to be set up and CANNOT be done once htable%finalize_setup has been called
    !> @param[in] index  contiguous index of the element (not the one in the sparse array)
    !> @param[in] element  new value of the element
    subroutine set_elem_sparse(this, index, element)
        class(sparse_lMat_t), intent(inout) :: this
        integer(int64), intent(in) :: index
        HElement_t(dp), intent(in) :: element

        integer(int64) :: pos
        character(*), parameter :: t_r = "set_elem_sparse"

        ! Add the new entry to the hashtable
        call this%htable%add_index(index, pos)
        ! And the entry to the values
        this%nonzero_vals%ptr(pos) = element
    end subroutine set_elem_sparse

    !------------------------------------------------------------------------------------------!

    !> Retrieve an element from the 6-index integrals stored in sparse format
    !> @param[in] index  contiguous index of the element to be retrieved
    !> @return element  value of the element with the given contiguous index
    function get_elem_sparse(this, index) result(element)
        class(sparse_lMat_t), intent(in) :: this
        integer(int64), intent(in) :: index
        HElement_t(dp) :: element

        integer(int64) :: lower, upper, i, pos
        logical :: t_found

        ! Lookup the element
        call this%htable%lookup(index, pos, t_found)
        ! If it is there, return the entry
        if(t_found) then
            element = this%nonzero_vals%ptr(pos)
        else
            ! Else, return 0 (zero matrix element)
            element = 0.0_dp
        end if

    end function get_elem_sparse

    !------------------------------------------------------------------------------------------!    

    !> Read the 6-index integrals from a file to sparse format
    !> @param[in] filename  name of the file to read from
    subroutine read_sparse(this, filename)
        class(sparse_lMat_t), intent(inout) :: this
        character(*), intent(in) :: filename
        character(*), parameter :: t_r = "read_sparse"        
#ifdef USE_HDF_        
        type(lMat_hdf5_read_t) :: reader
        integer(hsize_t) :: nInts
        ! There is no sparse ascii reader yet, so filename is never used
        if(.not. tHDF5LMat) call stop_all(t_r, "Sparse 6-index integrals require hdf5 format")

        call reader%open(filename, nInts)
        call this%alloc(nInts)
        call reader%loop_file(this)
        call this%htable%setup_offsets()
        ! The core energy has already been updated, no need to do so again
        call reader%loop_file(this)
        call this%htable%finalize_setup()
        call reader%close()
#else        

        unused_var(this)
        ! unused_var on strings is not supported by some older compilers
        if(len(filename) /= 0) continue
        call stop_all(t_r, "Sparse 6-index integrals are only available for hdf5 format")
#endif
    end subroutine read_sparse

    !------------------------------------------------------------------------------------------!

    !> This is the operation to be performed for sparse storage on each block of data read from an hdf5 file
    !! both arguments may or may not be still allocated upon return.
    !> @param[in,out] indices  chunk of indices read in from the file
    !> @param[in,out] entries  chunk of corresponding values
    subroutine read_op_sparse(this, indices, entries) 
        class(sparse_lMat_t), intent(inout) :: this
        ! We allow deallocation of indices/entries
        integer(int64), allocatable, intent(inout) :: indices(:,:), entries(:,:)

        integer(int64) :: block_size, i
        integer(int64), allocatable :: combined_inds(:)        

        block_size = size(indices, dim=2)
        allocate(combined_inds(block_size))
        ! Transfer the 6 orbitals to one contiguous index
        do i = 1, block_size
            combined_inds(i) = this%indexFunc(indices(1,i), indices(2, i), indices(3, i), &
                indices(4, i), indices(5, i), indices(6,i))
        end do
        ! We might need this memory - all these operations can be memory critical
        deallocate(indices)

        if(this%htable%known_conflicts()) then
            call this%read_data(combined_inds, entries(1,:))
        else
            call this%count_conflicts(combined_inds)
        endif

        deallocate(combined_inds)        
    end subroutine read_op_sparse

    !------------------------------------------------------------------------------------------!

    !> Loop through a chunk of indices and count the number of hash conflicts. This is required
    !! for setting up the hash table
    !> @param[in] indices  chunk combined 6-index values for the 6-index integrals
    subroutine count_conflicts(this, indices)
        class(sparse_lMat_t), intent(inout) :: this
        integer(int64), intent(in) :: indices(:)

        integer(int64) :: i, total_size
        integer(int64), allocatable :: tmp(:)

        ! Gather all read indices on node-root
        call gather_block(indices, tmp)
        total_size = size(tmp)

        if(iProcIndex_intra == 0) then
            do i = 1, total_size
                ! count_index is not threadsafe => only do it on node-root
                call this%htable%count_index(tmp(i))
            end do
        end if
        deallocate(tmp)
    end subroutine count_conflicts

    !------------------------------------------------------------------------------------------!

    !> Add the (combined) indices and the corresponding integral values to the sparse storage
    !> @param[in] indices  chunk of combined 6-index values
    !> @param[in] entries  corresponding values of the 6-index integrals
    subroutine read_data(this, indices, entries)
        class(sparse_lMat_t), intent(inout) :: this
        integer(int64), intent(in) :: indices(:), entries(:)

        integer(int64), allocatable :: tmp_inds(:), tmp_entries(:)
        integer(int64) :: i, total_size, pos
        ! For typecasts
        HElement_t(dp), parameter :: rVal = 0.0_dp
        ! Gather all data on node root
        call gather_block(indices, tmp_inds)
        call gather_block(entries, tmp_entries)
        total_size = size(tmp_inds)

        ! Then write there
        if(iProcIndex_intra == 0) then
            do i = 1, total_size
                call this%set_elem(tmp_inds(i), 3.0_dp * transfer(tmp_entries(i),rVal))
            end do
        endif
        deallocate(tmp_inds)
        deallocate(tmp_entries)
    end subroutine read_data

    !------------------------------------------------------------------------------------------!
    ! Generic auxiliary routine
    !------------------------------------------------------------------------------------------!

    !> Gather a chunk of data on node-root.
    !> @param[in] data_block  on each proc, the data from this proc to be gathered
    !> @param[out] tmp  on return, on node-root the gathered data from all procs on this node, empty on all other procs. Guaranteed to be allocated on return (of size 0 on other than node-root).
    subroutine gather_block(data_block, tmp)
        integer(int64), intent(in) :: data_block(:)
        integer(int64), allocatable, intent(out) :: tmp(:)

        integer(MPIArg) :: procs_per_node
        integer(MPIArg) :: this_block_size, total_size
        integer(MPIArg), allocatable :: block_sizes(:), displs(:)
        integer(MPIArg) :: ierr
        integer :: i

        call MPI_Comm_Size(mpi_comm_intra, procs_per_node, ierr)

        allocate(block_sizes(0:procs_per_node-1))
        this_block_size = int(size(data_block), MPIArg)
        ! Check how much data was read in in total (needs to be known on node root)
        call MPI_Gather(this_block_size, 1, MPI_INT, &
            block_sizes, 1, MPI_INT, 0, mpi_comm_intra, ierr)

        allocate(displs(0:procs_per_node-1))
        displs(0) = 0
        do i = 1, procs_per_node-1
            displs(i) = displs(i-1) + block_sizes(i-1)
        end do

        total_size = sum(block_sizes)
        if(iProcIndex_intra == 0) then
            allocate(tmp(total_size))
        else
            allocate(tmp(0))
        end if

        ! Gather all data on node-root
        call MPI_GatherV(data_block, this_block_size, MPI_INTEGER8, &
            tmp, block_sizes, displs, MPI_INTEGER8, 0, mpi_comm_intra, ierr)

        deallocate(block_sizes)
        deallocate(displs)
    end subroutine gather_block

    !------------------------------------------------------------------------------------------!
    ! Histogramming
    !------------------------------------------------------------------------------------------!

    !> Generate a histogram of the 6-index integrals and write it to stdout
    subroutine histogram_lMat(this)
        class(lMat_t), intent(in) :: this
        integer(int64) :: i
        integer :: thresh
        integer, parameter :: minExp = 10
        integer :: histogram(0:minExp)
        real :: ratios(0:minExp)

        histogram = 0
        do i = 1, this%lMat_size()
            do thresh = minExp,1,-1
                ! in each step, count all matrix elements that are below the threshold and
                ! have not been counted yet
                if(abs(this%get_elem(i))<0.1**(thresh)) then
                    histogram(thresh) = histogram(thresh) + 1
                    ! do not count this one again
                    exit
                endif
            end do
            ! the last check has a different form: everything that is bigger than 0.1 counts here
            if(abs(this%get_elem(i)) > 0.1) histogram(0) = histogram(0) + 1
        end do

        ratios(:) = real(histogram(:)) / real(this%lMat_size())
        ! print the ratios
        write(iout,*) "Matrix elements below", 0.1**(minExp), ":", ratios(minExp)
        do i = minExp-1,1,-1
            write(iout,*) "Matrix elements from", 0.1**(i+1),"to",0.1**(i),":",ratios(i)
        end do
        write(iout,*) "Matrix elements above", 0.1,":",ratios(0)
        write(iout,*) "Total number of logged matrix elements", sum(histogram)

    end subroutine histogram_lMat

    !------------------------------------------------------------------------------------------!
    ! HDF5 I/O class methods
    !------------------------------------------------------------------------------------------!
    
#ifdef USE_HDF_
    !> Open an hdf5 file containing 6-index integrals
    !> @param[in] filename  name of the file
    !> @param[out] nInts  number of integrals stored in the file (normally only nonzeros)
    subroutine open(this, filename, nInts)
        class(lMat_hdf5_read_t) :: this
        character(*), intent(in) :: filename
        integer(hsize_t), intent(out) :: nInts

        integer :: proc, i
        integer(hid_t) :: err
        integer(hsize_t) :: rest
        integer(hsize_t), allocatable :: counts(:)
        integer(MPIArg) :: procs_per_node, ierr
        character(*), parameter :: t_r = "lMat_hdf5_read_t%open"

        call h5open_f(err)
        call h5pcreate_f(H5P_FILE_ACCESS_F, this%plist_id, err)
        call h5pset_fapl_mpio_f(this%plist_id, mpi_comm_intra, mpiInfoNull, err)

        ! open the file
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, this%file_id, err, access_prp=this%plist_id)

        call h5gopen_f(this%file_id, nm_grp, this%grp_id, err)

        ! get the number of integrals
        call read_int64_attribute(this%grp_id, nm_nInts, nInts, required=.true.)
        write(iout,*) "Reading", nInts, "integrals"

        ! how many entries does each proc get?
        call MPI_Comm_Size(mpi_comm_intra, procs_per_node, ierr)
        allocate(counts(0:procs_per_node-1))
        allocate(this%offsets(0:procs_per_node-1))
        counts = nInts / int(procs_per_node, hsize_t)
        rest = mod(nInts, procs_per_node)
        if(rest>0) counts(0:rest-1) = counts(0:rest-1) + 1

        this%offsets(0) = 0
        do proc = 1, procs_per_node - 1
            this%offsets(proc) = this%offsets(proc-1) + counts(proc-1)
        end do
        ! the last element to read on each proc
        if(iProcIndex_intra.eq.procs_per_node - 1) then
            this%countsEnd = nInts - 1
        else
            this%countsEnd = this%offsets(iProcIndex_intra + 1) - 1
        end if

        call h5dopen_f(this%grp_id, nm_vals, this%ds_vals, err)
        call h5dopen_f(this%grp_id, nm_indices, this%ds_inds, err)        

    end subroutine open

    !------------------------------------------------------------------------------------------!

    !> Apply the read_op_hdf5 of an lMat to the data in the currently opened file
    !! The file will be read chunkwise and the read_op_hdf5 operation applied per chunk
    !> @param[in] lMat  the lMat object to read the data to
    subroutine loop_file(this, lMat)
        class(lMat_hdf5_read_t), intent(inout) :: this
        class(lMat_t), intent(inout) :: lMat

        real(dp) :: rVal
        logical :: running, any_running
        integer(hsize_t) :: blocksize, blockstart, blockend, this_blocksize
        integer(hsize_t), allocatable :: indices(:,:), entries(:,:)
        integer(MPIArg) :: ierr
        rVal = 0.0_dp        

        ! reserve max. 128MB buffer size for dumpfile I/O
        blocksize = 2_hsize_t**27 .div. (7*sizeof(0_int64))
        blockstart = this%offsets(iProcIndex_intra)

        blockend = min(blockstart + blocksize - 1, this%countsEnd)
        any_running = .true.
        running = .true.
        do while(any_running)
            if(running) then
                ! the number of elements to read in this block
                this_blocksize = blockend - blockstart + 1
            else
                this_blocksize = 0
            end if

            allocate(indices(6,this_blocksize), source = 0_int64)
            allocate(entries(1,this_blocksize), source = 0_int64)

            ! read in the data
            call read_2d_multi_chunk(&
                this%ds_vals, entries, H5T_NATIVE_REAL_8, &
                [1_hsize_t, this_blocksize],&
                [0_hsize_t, blockstart],&
                [0_hsize_t, 0_hsize_t])

            call read_2d_multi_chunk(&
                this%ds_inds, indices, H5T_NATIVE_INTEGER_8, &
                [6_hsize_t, this_blocksize], &
                [0_hsize_t, blockstart], &
                [0_hsize_t, 0_hsize_t])

            ! Do something with the read-in values
            ! This has to be threadsafe !!!
            call lMat%read_op_hdf5(indices, entries)

            ! the read_op is allowed to deallocate if memory has to be made available
            if(allocated(entries)) deallocate(entries)
            if(allocated(indices)) deallocate(indices)            

            ! set the size/offset for the next block
            if(running) then
                blockstart = blockend + 1
                blockend = min(blockstart + blocksize - 1, this%countsEnd)
                if(blockstart > this%countsEnd) running = .false.
            end if

            ! once all procs on this node are done reading, we can exit
            call MPI_ALLREDUCE(running, any_running, 1, MPI_LOGICAL, MPI_LOR, mpi_comm_intra, ierr)
        end do

    end subroutine loop_file

    !------------------------------------------------------------------------------------------!    

    !> Close the currently opened hdf5 file - requires a previous call to open()
    subroutine close(this)
        class(lMat_hdf5_read_t) :: this
        integer(hid_t) :: err
        call h5dclose_f(this%ds_vals, err)
        call h5dclose_f(this%ds_inds, err)
        ! close the file, finalize hdf5
        call h5gclose_f(this%grp_id, err)
        call h5pclose_f(this%plist_id, err)
        call h5fclose_f(this%file_id, err)
        call h5close_f(err)

    end subroutine close
#endif

end module LMat_class
