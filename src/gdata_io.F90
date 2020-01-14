module gdata_io
    ! This module manages if and how much global determinant data is read/write from/to
    ! the popsfile. 
    use LoggingData, only: tAccumPopsActive
    use global_det_data, only: writeFVals, readFVals, write_max_ratio, &
        set_all_max_ratios, writeAPVals, readAPVals
    use constants
    use CalcData, only: tAutoAdaptiveShift, tScaleBlooms
    use util_mod, only: operator(.div.)
#ifdef USE_HDF_
    use global_det_data, only: writeFValsAsInt, write_max_ratio_as_int, &
        readFValsAsInt, set_max_ratio_hdf5Int
    use hdf5, only: hsize_t
#endif
    implicit none

    private
    public :: gdata_io_t, resize_attribute

#ifdef USE_HDF_
    public :: clone_signs
#endif

    type gdata_io_t
        ! The gdata_io_t type does the transfer of global determinant data (gdata)
        ! to and from a popsfile. It contains the information on how the data is
        ! arranged in the popsfile and is the interface between the popsfile modules
        ! and the global_determinant_data module
        ! The read/write buffers themselves are not part of this class and are at
        ! the responsibility of the popsfile modules
        private
        ! size of the buffer per entry used for global data i/o
        integer :: gdata_size = 0
        ! ranges to read/write the data to in the buffer
        integer :: fvals_start, fvals_end
        integer :: max_ratio_start, max_ratio_end
        integer :: apvals_start, apvals_end
    contains
        ! initialization routine
        procedure :: init_gdata_io
        ! size per entry
        procedure :: entry_size, t_io
        ! write the read data to the global det data
        procedure :: write_gdata, read_gdata
#ifdef USE_HDF_
        ! hdf5 currently uses another output format
        procedure :: write_gdata_hdf5, read_gdata_hdf5
        ! isolate the dynamically sized part (i.e. the part depending on the number of runs
        procedure :: clone_gdata
#endif
    end type gdata_io_t

contains

    !------------------------------------------------------------------------------------------!    

    subroutine init_gdata_io(this, t_aas, t_ms, t_ap, fvals_size_in, max_ratio_size_in, apvals_size_in)
        ! Initialize a gdata_io_t object. This sets the read/write ranges for buffers
        ! and the size of required buffers
        ! Input: t_aas - Is auto adaptive shift data read/written (acc. rates)?
        !        t_ms - Is scale bloom data read/written (hij/pgen ratios)?
        !        t_ap - Is accumplated populations data read/written?
        !        fvals_size - Size of the aas data (only referenced when t_aas is true)
        !        max_ratio_size - size of the ms data (only referenced when t_ms is true)
        !        apvals_size - size of the accumlated populations data (only referenced when t_apvals is true)
        class(gdata_io_t) :: this
        logical, intent(in) :: t_aas, t_ms, t_ap
        integer, intent(in) :: fvals_size_in, max_ratio_size_in, apvals_size_in
        integer :: fvals_size, max_ratio_size, apvals_size

        ! how much data is to be read?        
        this%gdata_size = 0

        if(.not. t_aas) then
            fvals_size = 0
        else
            fvals_size = fvals_size_in
        endif
        if(.not. t_ms) then
            max_ratio_size = 0
        else
            max_ratio_size = max_ratio_size_in
        endif
        if(.not. t_ap) then
            apvals_size = 0
        else
            apvals_size = apvals_size_in
        endif

        ! set the range of each section 
        this%fvals_start = 1
        this%fvals_end = this%fvals_start + fvals_size - 1
        this%gdata_size = this%gdata_size + fvals_size

        this%max_ratio_start = this%fvals_start + fvals_size
        this%max_ratio_end = this%max_ratio_start + max_ratio_size - 1
        this%gdata_size = this%gdata_size + max_ratio_size

        this%apvals_start = this%max_ratio_end + max_ratio_size
        this%apvals_end = this%apvals_start + apvals_size - 1
        this%gdata_size = this%gdata_size + apvals_size

    end subroutine init_gdata_io

    !------------------------------------------------------------------------------------------!    

    function entry_size(this) result(e_size)
        ! Return the size of one entry in the gdata buffer (first dimension of buffers)
        ! Output: e_size - first dimension of a buffer to be used with this gdata_io_t object
        class(gdata_io_t), intent(in) :: this
        integer :: e_size

        e_size = this%gdata_size
    end function entry_size

    !------------------------------------------------------------------------------------------!    

    function t_io(this) result(t_do_io)
        ! Returns if any gdata read/write is happening
        ! Output: t_do_io - true if gdata shall be read/written
        class(gdata_io_t), intent(in) :: this
        logical :: t_do_io

        t_do_io = this%entry_size() > 0
    end function t_io
    
    !------------------------------------------------------------------------------------------!    
    
    subroutine read_gdata(this, gdata_buf, ndets, initial)
        ! Reads the gdata from a buffer to the global_det_data array
        ! Input: gdata_buf - Buffer containing the read data, first dimension has
        !                    to be this%entry_size()
        !        ndets - number of entries to read (can be less than the size of gdata_buf)
        !        initial - optionally, an offset where to start writing
        class(gdata_io_t), intent(in) :: this
        real(dp), intent(in) :: gdata_buf(:,:)
        integer, intent(in) :: ndets
        integer, intent(in), optional :: initial

        ! do a sanity check: there has to be some data
        if(this%t_io()) then
            if(this%entry_size() <= size(gdata_buf, dim=1)) then
                if(tAutoAdaptiveShift) then
                    ! set the global det data for auto adaptive shift
                    call readFVals(&
                        gdata_buf(this%fvals_start:this%fvals_end,:), ndets, initial)
                endif
                if(tScaleBlooms) then
                    ! set the global det data for bloom scaling
                    call set_all_max_ratios(&
                        gdata_buf(this%max_ratio_start:this%max_ratio_end,:), ndets, initial)
                end if
                if(tAccumPopsActive) then
                    ! set the global det data for accumlated population
                    call readAPVals(&
                        gdata_buf(this%apvals_start:this%apvals_end,:), ndets, initial)
                endif
            else
                write(iout,*) "WARNING: Dimension mismatch in read_gdata, ignoring all read data"
            endif
        endif
    end subroutine read_gdata

    !------------------------------------------------------------------------------------------!    

    subroutine write_gdata(this, gdata_buf, ndets)
        ! Write this calculations gdata to a buffer
        ! Input: gdata_buf - on return, contains the gdata of the first ndets determinants
        !                    has to have a first dimension of this%entry_size()
        !        ndets - number of determinants to write
        class(gdata_io_t), intent(inout) :: this
        real(dp), intent(out) :: gdata_buf(:,:)
        integer, intent(in) :: ndets

        ! sanity check
        if(this%t_io()) then
            if(this%entry_size() <= size(gdata_buf, dim=1)) then
                ! write the global det data to the buffer
                if(tAutoAdaptiveShift) call writeFVals(&
                    gdata_buf(this%fvals_start:this%fvals_end,:), ndets)
                if(tScaleBlooms) call write_max_ratio(&
                    gdata_buf(this%max_ratio_start:this%max_ratio_end,:), ndets)
                if(tAccumPopsActive) call writeAPVals(&
                    gdata_buf(this%apvals_start:this%apvals_end,:), ndets)
            else
                write(iout,*) "WARNING: Dimension mismatch in write_gdata, writing 0"
                gdata_buf = 0.0_dp
            endif
        end if
    end subroutine write_gdata

    !------------------------------------------------------------------------------------------!
    ! HDF popsfile funcitonalities - require building with hdf5 enabled
    !------------------------------------------------------------------------------------------!    
#ifdef USE_HDF_
    
    subroutine read_gdata_hdf5(this, gdata_buf, pos)
        ! Read the gdata of a single determinant from an hdf5 file
        ! Input: gdata_buf - gdata for the determinant at pos, has to be of size this%entry_size
        !        pos - position to put read data into 
        class(gdata_io_t), intent(in) :: this
        integer(hsize_t), intent(in) :: gdata_buf(:)
        integer, intent(in) :: pos
        ! only print one warning
        logical, save :: t_warn = .true.

        if(this%t_io()) then
            if(this%entry_size() <= size(gdata_buf)) then
                if(tAutoAdaptiveShift) then
                    call readFValsAsInt(gdata_buf(this%fvals_start:this%fvals_end), pos)
                endif
                if(tScaleBlooms) then
                    ! set the global det data for bloom scaling
                    call set_max_ratio_hdf5Int(&
                        gdata_buf(this%max_ratio_start:this%max_ratio_end), pos)
                end if
                if(tAccumPopsActive) then
                    call readAPValsAsInt(gdata_buf(this%apvals_start:this%apvals_end), pos)
                endif
            else
                if(t_warn) then
                    write(iout,*) "WARNING: Dimension mismatch in read_gdata_hdf5, ignoring read data"
                    t_warn = .false.
                endif
            endif
        end if
    end subroutine read_gdata_hdf5

    !------------------------------------------------------------------------------------------!

    subroutine write_gdata_hdf5(this, gdata_buf, pos)
        ! Write the gdata of a single determinant to a buffer usable in hdf5 popsfiles
        ! Input: gdata_buf - gdata for the determinant at pos, has to be of size this%entry_size
        !        pos - position to get the written data from
        class(gdata_io_t), intent(in) :: this
        integer(hsize_t), intent(out) :: gdata_buf(:)
        integer, intent(in) :: pos

        logical :: t_aas, t_sb, t_ap

        ! if these are above 0, the option has been set and memory is reserved
        t_aas = this%fvals_end - this%fvals_start + 1> 0
        t_sb = this%max_ratio_end - this%max_ratio_start + 1> 0
        t_ap = this%apvals_end - this%apvals_start + 1> 0

        if(this%entry_size() <= size(gdata_buf, dim=1)) then
            if(t_aas) then
                ! write the fvals to the buffer at the respective position
                call writeFValsAsInt(gdata_buf(this%fvals_start:this%fvals_end), pos)
            endif

            if(t_sb) then
                ! write the ratios to the buffer at the respective position
                call write_max_ratio_as_int(gdata_buf(this%max_ratio_start:this%max_ratio_end), pos)
            endif

            if(t_ap) then
                ! write the apvals to the buffer at the respective position
                call writeAPValsAsInt(gdata_buf(this%apvals_start:this%apvals_end), pos)
            endif
        else
            write(iout,*) "WARNING: Dimension mismatch in write_gdata_hdf5, writing 0"
            gdata_buf = 0_hsize_t
        endif
    end subroutine write_gdata_hdf5

    !------------------------------------------------------------------------------------------!
    ! Generic HDF functionality
    !------------------------------------------------------------------------------------------!    
    
    subroutine clone_gdata(this, gdata_buf, tmp_fvals_size, fvals_size, tmp_apvals_size, & 
                           apvals_size, nsigns)
        ! expand the global det data:
        ! clone the fvals from tmp_fvals_size to fvals_size, leaving the rest of
        ! the data as it is
        ! Input: gdata_buf - On input, the read-in gdata sized with tmp_fvals_size
        !                    on return, the re-sized gdata sized with fvals_size
        !        tmp_fvals_size - size of the acc. rate part of gdata on input
        !        fvals_size - size of the acc. rate part of gdata on return
        !        nsigns - number of determinants affected
        class(gdata_io_t), intent(inout) :: this
        integer(hsize_t), allocatable, intent(inout) :: gdata_buf(:,:)
        integer, intent(in) :: tmp_fvals_size, fvals_size, tmp_apvals_size, apvals_size
        integer(int64), intent(in) :: nsigns

        integer(hsize_t), allocatable :: tmp_fvals_acc(:,:), tmp_fvals_tot(:,:), tmp_mr(:,:)
        integer(hsize_t), allocatable :: tmp_apvals_sum(:,:), tmp_apvals_iter(:,:)
        integer :: max_ratio_size, gdata_size
        logical :: t_aas, t_sb, t_ap
        integer :: tmp_tot_start, tot_start, tmp_acc_size, acc_size 
        integer :: tmp_iter_start, iter_start, tmp_sum_size, sum_size 

        ! if size of data is above 0, the option has been set and memory is reserved
        t_aas = this%fvals_end - this%fvals_start + 1> 0
        max_ratio_size = this%max_ratio_end - this%max_ratio_start + 1
        t_sb = max_ratio_size> 0
        t_ap = this%apvals_end - this%apvals_start + 1> 0

        
        ! Only the sizes fvals and apvals depend on lenof_sign
        if(t_aas .or. t_ap) then
            ! can only resize buffers with correct size
            if(this%entry_size() == size(gdata_buf, dim=1)) then
                if(t_aas) then
                    ! copy the fvals to a temporary - one for tot. and one for acc. spawns
                    ! each half the total size
                    ! we do this to clone each of them independently
                    tmp_acc_size = tmp_fvals_size .div. 2
                    tmp_tot_start = this%fvals_start + tmp_acc_size                
                    allocate(tmp_fvals_acc(tmp_acc_size, nsigns))
                    allocate(tmp_fvals_tot(tmp_acc_size, nsigns))
                    tmp_fvals_acc(:,:) = gdata_buf(this%fvals_start:tmp_tot_start - 1,:)
                    tmp_fvals_tot(:,:) = gdata_buf(tmp_tot_start:this%fvals_end,:)

                    acc_size = fvals_size .div. 2
                    ! clone the content of the temporary - clone the first
                    ! and second half seperately
                    call clone_signs(tmp_fvals_acc, tmp_acc_size, acc_size, nsigns)
                    ! tot and acc have to have the same size
                    call clone_signs(tmp_fvals_tot, tmp_acc_size, acc_size, nsigns)
                endif
                if(t_sb) then
                    allocate(tmp_mr(max_ratio_size, nsigns))
                    tmp_mr(:,:) = gdata_buf(this%max_ratio_start:this%max_ratio_end,:)
                endif
                if(t_ap) then
                    ! copy the apvals (pops sum and pops iter) to a temporary - 
                    ! only the size of pops sum dependes on lenof_sign
                    tmp_sum_size = tmp_apvals_size - 1
                    tmp_iter_start = this%apvals_start + tmp_sum_size
                    allocate(tmp_apvals_sum(tmp_sum_size, nsigns))
                    allocate(tmp_apvals_iter(1, nsigns))
                    tmp_apvals_sum(:,:) = gdata_buf(this%apvals_start:tmp_iter_start - 1,:)
                    tmp_apvals_iter(:,:) = gdata_buf(tmp_iter_start:this%apvals_end,:)

                    sum_size = apvals_size -1
                    ! clone the content of the temporary pops sum
                    call clone_signs(tmp_apvals_sum, tmp_sum_size, sum_size, nsigns)
                endif
                deallocate(gdata_buf)
                ! adjust the gdata offsets of this io handler
                call this%init_gdata_io(t_aas, t_sb, t_ap, fvals_size, max_ratio_size, apvals_size)

                ! resize the buffer - with the new gdata_size
                gdata_size = this%entry_size()            
                allocate(gdata_buf(gdata_size,nsigns))

                if(t_aas)then
                    tot_start = this%fvals_start + acc_size
                    ! fill in the resized data                
                    gdata_buf(this%fvals_start:tot_start-1,:) = tmp_fvals_acc(:,:)
                    gdata_buf(tot_start:this%fvals_end,:) = tmp_fvals_tot(:,:)
                    deallocate(tmp_fvals_tot)
                    deallocate(tmp_fvals_acc)
                endif
                if(t_sb) then
                    gdata_buf(this%max_ratio_start:this%max_ratio_end,:) = tmp_mr(:,:)
                    deallocate(tmp_mr)
                end if
                if(t_ap)then
                    iter_start = this%apvals_start + sum_size
                    ! fill in the resized data                
                    gdata_buf(this%apvals_start:iter_start-1,:) = tmp_apvals_sum(:,:)
                    gdata_buf(iter_start:this%apvals_end,:) = tmp_apvals_iter(:,:)
                    deallocate(tmp_apvals_sum)
                    deallocate(tmp_apvals_iter)

                endif
            else
                write(iout,*) "WARNING: Dimension mismatch in clone_gdata. No data read"
            endif
        endif
    end subroutine clone_gdata

    !------------------------------------------------------------------------------------------!
    
    subroutine clone_signs(tmp_sgns, tmp_lenof_sign, lenof_sign, num_signs)
        ! Resize a 2-D array from one first dimension to another by either deleting
        ! or copying entries
        implicit none
        ! Input: tmp_sgns - temporary storing the signs to be adapted to this runs number
        !                   of replicas
        !        tmp_lenof_sign - first dimension of tmp_sgns on input
        !        lenof_sign - first dimension of  tmp_sgns on return
        !        num_signs - number of entries in tmp_sgns to copy
        integer(hsize_t), allocatable, intent(inout) :: tmp_sgns(:,:)
        integer(hsize_t), intent(in) :: num_signs
        integer, intent(in) :: tmp_lenof_sign, lenof_sign

        ! a temporary buffer to store the old signs while reallocating tmp_sgns
        integer(hsize_t), allocatable :: sgn_store(:,:)
        integer :: ierr, i

        if(allocated(tmp_sgns)) then
            ! copy the signs to a temporary
            allocate(sgn_store(tmp_lenof_sign,num_signs),stat=ierr)
            sgn_store(:,:) = tmp_sgns(:,:)

            ! now, resize tmp_sgns
            deallocate(tmp_sgns)
            allocate(tmp_sgns(lenof_sign,num_signs),stat=ierr)

            ! and clone the signs to match lenof_sign numbers per entry
            do i = 1, int(num_signs)
                ! depending on if we want to remove or add replicas,
                ! shrink or expand the signs
                call resize_sign(tmp_sgns(:,i),sgn_store(:,i))
            end do

            deallocate(sgn_store)
        else
            write(6,*) "WARNING: Attempted to adjust lenof_sign for an empty input"
            ! throw a warning
        endif

    end subroutine clone_signs

    !------------------------------------------------------------------------------------------!

    subroutine resize_sign(out_sgn, in_sgn)
        ! Copy data between two arrays of different size
        ! Input: in_sgn - array to copy from
        !        out_sgn - on return, contains the data from in_sgn
        !                  if out_sgn is larger than in_sgn, the last entry is
        !                  multiplied, if it is smaller, the last entries are left out
        implicit none

        integer(hsize_t), intent(out) :: out_sgn(:)
        integer(hsize_t), intent(in) :: in_sgn(:)

        integer :: out_size, in_size

        out_size = size(out_sgn)
        in_size = size(in_sgn)

        if(out_size < in_size) then
            ! remove the last entries from the input
            out_sgn(1:out_size) = in_sgn(1:out_size)
        else
            ! copy the last replicas to fill up to the desired number
            out_sgn(1:in_size) = in_sgn(1:in_size)
            out_sgn(in_size+1:out_size) = in_sgn(in_size)
        endif
    end subroutine resize_sign

#endif

    !------------------------------------------------------------------------------------------!

    subroutine resize_attribute(attribute, new_size)
        ! take an array and expand/shrink it to a new size
        ! Input: attribute - array to resize
        !        new_size - new size of the array. If larger than the current one,
        !                   data will be duplicated, if smaller, data will be deleted
        implicit none
        integer, intent(in) :: new_size
        real(dp), allocatable, intent(inout) :: attribute(:)

        real(dp), allocatable :: tmp(:)
        integer :: old_size
        integer :: ierr

        old_size = size(attribute)
        !store the old entries
        allocate(tmp(old_size), stat = ierr)
        tmp(:) = attribute(:)

        deallocate(attribute)
        allocate(attribute(new_size), stat = ierr)

        ! resize
        if(old_size < new_size) then
            attribute(1:old_size) = tmp(1:old_size)
            attribute(old_size+1:new_size) = tmp(old_size)
        else
            attribute(1:new_size) = tmp(1:new_size)
        end if

        deallocate(tmp)
    end subroutine resize_attribute

!------------------------------------------------------------------------------------------!
   
end module gdata_io
