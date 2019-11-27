module gdata_io
    ! This module manages if and how much global determinant data is read/write from/to
    ! the popsfile
    use LoggingData, only: tPopAutoAdaptiveShift, tPopScaleBlooms
    use global_det_data, only: writeFFunc, write_max_ratio, set_all_max_ratios, &
        set_tot_acc_spawns, writEFFuncAsInt, write_max_ratio_as_int, &
        set_tot_acc_spawns_hdf5Int, set_max_ratio_hdf5Int
    use constants
    use CalcData, only: tAutoAdaptiveShift, tScaleBlooms
#ifdef USE_HDF_
    use hdf5, only: hsize_t
#endif
    implicit none

    private
    public :: gdata_io_t, resize_attribute
#ifdef USE_HDF_
    public :: clone_signs
#endif

    type gdata_io_t
        private
        ! size of the buffer per entry used for global data i/o
        integer :: gdata_size = 0
        ! ranges to read/write the data to in the buffer
        integer :: fvals_start, fvals_end
        integer :: max_ratio_start, max_ratio_end
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

    subroutine init_gdata_io(this, t_aas, t_ms, fvals_size, max_ratio_size)
        class(gdata_io_t) :: this
        logical, intent(in) :: t_aas, t_ms
        integer, intent(in) :: fvals_size, max_ratio_size
        ! how much data is to be read?        
        this%gdata_size = 0
        ! set the offset for reading to the gdata buffer
        this%fvals_start = 1
        this%fvals_end = 0
        this%max_ratio_start = 1
        this%max_ratio_end = 0
        if(t_aas) then
            this%gdata_size = this%gdata_size + fvals_size
            ! increase the offset for all following data
            this%fvals_end = this%fvals_start + fvals_size - 1
            this%max_ratio_start = this%fvals_end + 1
        endif
        if(t_ms) then
            this%gdata_size = this%gdata_size + max_ratio_size
            this%max_ratio_end = this%max_ratio_start + max_ratio_size - 1
        endif        
    end subroutine init_gdata_io

    !------------------------------------------------------------------------------------------!    

    function entry_size(this) result(e_size)       
        class(gdata_io_t), intent(in) :: this
        integer :: e_size

        e_size = this%gdata_size
    end function entry_size

    !------------------------------------------------------------------------------------------!    

    function t_io(this) result(t_do_io)
        class(gdata_io_t), intent(in) :: this
        logical :: t_do_io

        t_do_io = this%entry_size() > 0
    end function t_io
    
    !------------------------------------------------------------------------------------------!    
    
    subroutine read_gdata(this, gdata_buf, ndets, initial)
        class(gdata_io_t), intent(in) :: this
        real(dp), intent(in) :: gdata_buf(:,:)
        integer, intent(in) :: ndets
        integer, intent(in), optional :: initial

        ! do a sanity check: there has to be some data
        if(this%t_io()) then
            if(tPopAutoAdaptiveShift .and. tAutoAdaptiveShift) then
                ! set the global det data for auto adaptive shift
                call set_tot_acc_spawns(&
                    gdata_buf(this%fvals_start:this%fvals_end,:), ndets, initial)
            endif
            if(tPopScaleBlooms .and. tScaleBlooms) then
                ! set the global det data for bloom scaling
                call set_all_max_ratios(&
                    gdata_buf(this%max_ratio_start:this%max_ratio_end,:), ndets, initial)
            end if
        endif
    end subroutine read_gdata

    !------------------------------------------------------------------------------------------!    

    subroutine write_gdata(this, gdata_buf, ndets)
        class(gdata_io_t), intent(inout) :: this
        real(dp), intent(out) :: gdata_buf(:,:)
        integer, intent(in) :: ndets

        ! sanity check
        if(this%t_io()) then
            ! write the global det data to the buffer
            if(tAutoAdaptiveShift) call writeFFunc(&
                gdata_buf(this%fvals_start:this%fvals_end,:), ndets)
            if(tScaleBlooms) call write_max_ratio(&
                gdata_buf(this%max_ratio_start:this%max_ratio_end,:), ndets)
        end if
    end subroutine write_gdata

    !------------------------------------------------------------------------------------------!
    ! HDF popsfile funcitonalities - require building with hdf5 enabled
    !------------------------------------------------------------------------------------------!    
#ifdef USE_HDF_
    
    subroutine read_gdata_hdf5(this, gdata_buf, pos)
        class(gdata_io_t), intent(in) :: this
        integer(hsize_t), intent(in) :: gdata_buf(:)
        integer, intent(in) :: pos

        if(this%t_io()) then
            if(tPopAutoAdaptiveShift .and. tAutoAdaptiveShift) then
                call set_tot_acc_spawns_hdf5Int(gdata_buf(this%fvals_start:this%fvals_end), pos)
            endif
            if(tPopScaleBlooms .and. tScaleBlooms) then
                ! set the global det data for bloom scaling
                call set_max_ratio_hdf5Int(&
                    gdata_buf(this%max_ratio_start:this%max_ratio_end), pos)
            end if
        end if
    end subroutine read_gdata_hdf5

    !------------------------------------------------------------------------------------------!

    subroutine write_gdata_hdf5(this, gdata_buf, ndets, max_ex)
        class(gdata_io_t), intent(in) :: this
        integer(hsize_t), intent(out) :: gdata_buf(:,:)
        integer, intent(in) :: ndets
        integer, intent(in), optional :: max_ex

        logical :: t_aas, t_sb

        ! if these are above 0, the option has been set and memory is reserved
        t_aas = this%fvals_end > 0
        t_sb = this%max_ratio_end > 0

        if(t_aas) then
            ! write the fvals to the buffer at the respective position
            call writeFFuncAsInt(gdata_buf(this%fvals_start:this%fvals_end,:), ndets, &
                max_ex)
        endif
        
        if(t_sb) then
            ! write the ratios to the buffer at the respective position
            call write_max_ratio_as_int(gdata_buf(this%max_ratio_start:this%max_ratio_end,:), &
                ndets, max_ex)
        endif
    end subroutine write_gdata_hdf5

    !------------------------------------------------------------------------------------------!
    ! Generic HDF functionality
    !------------------------------------------------------------------------------------------!    
    
    subroutine clone_gdata(this, gdata_buf, tmp_fvals_size, fvals_size)
        ! expand the global det data:
        ! clone the fvals from tmp_fvals_size to fvals_size, leaving the rest of
        ! the data as it is
        class(gdata_io_t), intent(inout) :: this
        integer(hsize_t), allocatable, intent(inout) :: gdata_buf(:,:)
        integer, intent(in) :: tmp_fvals_size, fvals_size

        integer(hsize_t), allocatable :: tmp_fvals(:,:), tmp_mr(:,:)
        integer :: nsigns, max_ratio_size, gdata_size, tmp_gdata_size
        logical :: t_aas, t_sb

        ! get the number of determinants
        nsigns = size(gdata_buf, dim=2)
        ! is there any aas data? If not, nothing to be done
        t_aas = this%fvals_end > 0
        ! this is the size of the rest of the data
        max_ratio_size = this%max_ratio_end - this%max_ratio_start + 1
        ! is there any max ratio data?
        t_sb = max_ratio_size > 0
        
        if(t_aas) then
            ! copy the fvals to a temporary
            allocate(tmp_fvals(tmp_fvals_size, nsigns))        
            tmp_fvals(:,:) = gdata_buf(this%fvals_start:this%fvals_end,:)

            ! clone the content of the temporary
            call clone_signs(tmp_fvals, tmp_fvals_size, fvals_size, int(nsigns,int64))

            ! resize the full buffer to fit the cloned data -> deallocate and then allocate
            tmp_gdata_size = this%entry_size()
            ! if there is data remaining, copy it along
            if(t_sb) then
                allocate(tmp_mr(tmp_gdata_size - tmp_fvals_size, nsigns))
                tmp_mr(:,:) = gdata_buf(this%max_ratio_start:this%max_ratio_end,:)
            endif
            deallocate(gdata_buf)
            ! adjust the gdata offsets of this io handler
            call this%init_gdata_io(t_aas, t_sb, fvals_size, max_ratio_size)

            ! resize the buffer - with the new gdata_size
            gdata_size = this%entry_size()            
            allocate(gdata_buf(gdata_size,nsigns))
            
            ! fill in the resized data
            gdata_buf(this%fvals_start:this%fvals_end,:) = tmp_fvals(:,:)
            if(t_sb) then
                gdata_buf(this%max_ratio_start:this%max_ratio_end,:) = tmp_mr(:,:)
                deallocate(tmp_mr)
            end if
            deallocate(tmp_fvals)
        endif
    end subroutine clone_gdata

    !------------------------------------------------------------------------------------------!
    
    subroutine clone_signs(tmp_sgns, tmp_lenof_sign, lenof_sign, num_signs)
      implicit none
      ! expand/shrink the sign to the target lenof_sign
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
            if(tmp_lenof_sign > lenof_sign) then
               call shrink_sign(tmp_sgns(:,i),lenof_sign,sgn_store(:,i),tmp_lenof_sign)
            else
               call expand_sign(tmp_sgns(:,i),lenof_sign,sgn_store(:,i),tmp_lenof_sign)
            endif
         end do

         deallocate(sgn_store)
      else
         write(6,*) "WARNING: Attempted to adjust lenof_sign for an empty input"
         ! throw a warning
      endif

    end subroutine clone_signs

!------------------------------------------------------------------------------------------!

    subroutine shrink_sign(out_sgn, out_size, in_sgn, in_size)
      implicit none

      integer, intent(in) :: out_size, in_size
      integer(hsize_t), intent(out) :: out_sgn(out_size)
      integer(hsize_t), intent(in) :: in_sgn(in_size)

      ! remove the last entries from the input
      out_sgn(1:out_size) = in_sgn(1:out_size)
    end subroutine shrink_sign

!------------------------------------------------------------------------------------------!

    subroutine expand_sign(out_sgn, out_size, in_sgn, in_size)
      implicit none

      integer, intent(in) :: out_size, in_size
      integer(hsize_t), intent(out) :: out_sgn(out_size)
      integer(hsize_t), intent(in) :: in_sgn(in_size)

      ! copy the last replica to fill up to the desired number
      out_sgn(1:in_size) = in_sgn(1:in_size)
      out_sgn(in_size+1:out_size) = in_sgn(in_size)
    end subroutine expand_sign        
#endif

        !------------------------------------------------------------------------------------------!

    subroutine resize_attribute(attribute, new_size, old_size)
      ! take an array and expand/shrink it to a new size
      implicit none
      integer, intent(in) :: new_size, old_size
      real(dp), allocatable, intent(inout) :: attribute(:)

      real(dp), allocatable :: tmp(:)
      integer :: ierr

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
