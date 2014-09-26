#include "macros.h"

    module global_det_data

        use FciMCData, only: MaxWalkersPart
        use LoggingData, only: tRDMonFly, tExplicitAllRDM
        use constants
        use util_mod
        implicit none

        ! This is the tag for allocation/deallocation
        private :: glob_tag
        integer :: glob_tag = 0

        ! This is the data used to find the elements inside the storage array
        private :: pos_hel, len_hel

        ! The diagonal matrix element is always stored. As it is a real value it
        ! always has a length of 1 (never cplx). Therefore, encode these values
        ! as parameters to assist optimisation.
        integer, parameter :: pos_hel = 1, len_hel = 1
        
        ! Average sign and first occupation of iteration
        private :: pos_av_sgn, len_av_sgn, pos_iter_occ, len_iter_occ
        integer :: pos_av_sgn, len_av_sgn
        integer :: pos_iter_occ, len_iter_occ

        ! And somewhere to store the actual data
        real(dp), pointer :: global_determinant_data(:,:) => null()

        interface set_av_sgn
            module procedure set_av_sgn_sgl
            module procedure set_av_sgn_all
        end interface

        interface get_av_sgn
            module procedure get_av_sgn_sgl
            module procedure get_av_sgn_all
        end interface

        interface set_iter_occ
            module procedure set_iter_occ_sgl
            module procedure set_iter_occ_all
        end interface

        interface get_iter_occ
            module procedure get_iter_occ_sgl
            module procedure get_iter_occ_all
        end interface

    contains

        subroutine init_global_det_data ()

            ! Initialise the global storage of determinant specific persistent
            ! data
            !
            ! --> This is the data that should not be transmitted with each
            !     particle
            ! --> It is not stored in the bit representation

            integer :: tot_len
            integer :: ierr
            character(*), parameter :: t_r = 'init_global_det_data'

        ! The position and size of diagonal matrix elements in the array.
        ! This is set as a module wide parameter, rather than as runtime, as
        ! it is constant, and this aids optimisation.
        ! pos_hel = 1
        ! len_hel = 1

        ! If we are using calculating RDMs stochastically, need to include the
        ! average sign and the iteration on which it became occupied.
        if (tRDMonFly .and. .not. tExplicitAllRDM) then
            len_av_sgn = lenof_sign
            len_iter_occ = lenof_sign
            write(6, '(" The average current signs before death will be stored&
                       & for use in the RDMs.")')
        else
            len_av_sgn = 0
            len_iter_occ = 0
        end if

        ! Get the starting positions
        pos_av_sgn = pos_hel + len_hel
        pos_iter_occ = pos_av_sgn + len_av_sgn

        tot_len = len_hel + len_av_sgn + len_iter_occ

        ! Allocate and log the required memory (globally)
        allocate(global_determinant_data(tot_len, MaxWalkersPart), stat=ierr)
        log_alloc(global_determinant_data, glob_tag, ierr)

        write(6,'(a,f14.6,a)') &
            ' Determinant related persistent storage requires: ', &
            8_dp * tot_len / 1048576_dp, ' Mb / processor'

        ! As an added safety feature
        tot_len = 0_dp
                         
    end subroutine


    subroutine clean_global_det_data()

        character(*), parameter :: this_routine = 'clean_global_det_data'

        ! Cleanup the global storage for determinant specific data

        if (associated(global_determinant_data)) then
            deallocate(global_determinant_data)
            log_dealloc(glob_tag)
            glob_tag = 0
            nullify(global_determinant_data)
        end if

    end subroutine


    ! These are accessor functions to currentH --> they are data access/setting
    !
    ! Access the global_determinant_data structure by site index (j).

    subroutine set_det_diagH(j, hel_r)

        ! Diagonal matrix elements are real --> store a real value

        integer, intent(in) :: j
        real(dp), intent(in) :: hel_r

        global_determinant_data(pos_hel, j) = hel_r

    end subroutine

    function det_diagH(j) result(hel_r)
        
        integer, intent(in) :: j
        real(dp) :: hel_r

        hel_r = global_determinant_data(pos_hel, j)

    end function


    subroutine set_av_sgn_sgl (j, part, av_sgn)

        integer, intent(in) :: j, part
        real(dp), intent(in) :: av_sgn

        global_determinant_data(pos_av_sgn + part - 1, j) = av_sgn

    end subroutine

    subroutine set_av_sgn_all (j, av_sgn)

        integer, intent(in) :: j
        real(dp), intent(in) :: av_sgn(len_av_sgn)

        global_determinant_data(pos_av_sgn:pos_av_sgn + len_av_sgn - 1, j) = &
            av_sgn

    end subroutine

    function get_av_sgn_sgl (j, part) result(av_sgn)

        integer, intent(in) :: j, part
        real(dp) :: av_sgn

        av_sgn = global_determinant_data(pos_av_sgn + part - 1, j)

    end function

    function get_av_sgn_all (j) result(av_sgn)

        integer, intent(in) :: j
        real(dp) :: av_sgn(len_av_sgn)

        av_sgn = global_determinant_data(pos_av_sgn:&
                                       pos_av_sgn + len_av_sgn - 1, j)

    end function


    subroutine set_iter_occ_sgl (j, part, iter_occ)

        ! It is unusual, but all the RDM code uses real(dp) values for the
        ! current iteration. As floats store integers exactly up to a
        ! sensible limit, this is just fine, and simplifies this code.
        !
        ! But it is weird.

        integer, intent(in) :: j, part
        real(dp), intent(in) :: iter_occ

        global_determinant_data(pos_iter_occ + part - 1, j) = iter_occ

    end subroutine

    subroutine set_iter_occ_all (j, iter_occ)

        ! It is unusual, but all the RDM code uses real(dp) values for the
        ! current iteration. As floats store integers exactly up to a
        ! sensible limit, this is just fine, and simplifies this code.
        !
        ! But it is weird.

        integer, intent(in) :: j
        real(dp), intent(in) :: iter_occ(len_iter_occ)

        global_determinant_data(pos_iter_occ: &
                                pos_iter_occ + len_iter_occ - 1, j) = iter_occ

    end subroutine

    function get_iter_occ_sgl (j, part) result(iter_occ)

        integer, intent(in) :: j, part
        real(dp) :: iter_occ

        iter_occ = global_determinant_data(pos_iter_occ + part - 1, j)

    end function

    function get_iter_occ_all (j) result(iter_occ)

        integer, intent(in) :: j
        real(dp) :: iter_occ(len_iter_occ)

        iter_occ = global_determinant_data(pos_iter_occ: &
                                   pos_iter_occ + len_iter_occ - 1, j)

    end function


end module
