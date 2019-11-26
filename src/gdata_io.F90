module gdata_io
    ! This module manages if and how much global determinant data is read/write from/to
    ! the popsfile
    use LoggingData, only: tPopAutoAdaptiveShift, tPopScaleBlooms
    use global_det_data, only: writeFFunc, write_max_ratio, set_all_max_ratios, &
        set_tot_acc_spawns, fvals_size, max_ratio_size
    use constants
    use CalcData, only: tAutoAdaptiveShift, tScaleBlooms
    implicit none

    private
    public :: gdata_io_t

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
        procedure :: entry_size

        ! write the read data to the global det data
        procedure :: write_gdata, read_gdata
    end type gdata_io_t

contains

    !------------------------------------------------------------------------------------------!    

    subroutine init_gdata_io(this, t_aas, t_ms)
        class(gdata_io_t) :: this
        logical, intent(in) :: t_aas, t_ms
        ! how much data is to be read?        
        this%gdata_size = 0
        ! set the offset for reading to the gdata buffer
        this%fvals_start = 1
        this%max_ratio_start = 1
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
    
    subroutine read_gdata(this, gdata_buf, ndets, initial)
        class(gdata_io_t), intent(in) :: this
        real(dp), intent(in) :: gdata_buf(:,:)
        integer, intent(in) :: ndets
        integer, intent(in), optional :: initial

        ! do a sanity check: there has to be some data
        if(this%entry_size() > 0) then
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
        if(this%entry_size() > 0) then
            ! write the global det data to the buffer
            if(tAutoAdaptiveShift) call writeFFunc(&
                gdata_buf(this%fvals_start:this%fvals_end,:), ndets)
            if(tScaleBlooms) call write_max_ratio(&
                gdata_buf(this%max_ratio_start:this%max_ratio_end,:), ndets)
        end if
    end subroutine write_gdata
end module gdata_io
