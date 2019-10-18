module csf_data
    use SystemData, only: nel
    implicit none

!    integer, parameter :: csf_orbital_mask = Z'1fffffff'
    integer, parameter :: csf_orbital_mask = 536870911
    integer, parameter :: csf_test_bit = 31
    integer, parameter :: csf_yama_bit = 30
    integer, parameter :: csf_ms_bit = 29

contains

    logical pure function iscsf (nI)

        ! Test if the specified determinant is a csf.
        !
        ! In:  nI    - Determinant in integer form
        ! Ret: iscsf - .true. if the determinant is a csf. .false. otherwise
        integer, dimension(:), intent(in) :: nI

        iscsf = btest(nI(1),csf_test_bit)
    end function

    subroutine csf_sort_det_block (dets, ndets, nopen)

        ! Sort a set of determinants, sorted as if they were a csf (closed
        ! orbitals followed by open ones) into normal order. This assumes that
        ! all the determinants being sorted have the same structure (ie they
        ! can all be sorted using moves from the first.
        !
        ! In:  ndets         - The number of determinants
        !      nopen         - The number of open shell electrons
        ! I/O: dets(ndets,:) - The determinant array to sort
        integer, intent(in) :: ndets, nopen
        integer, intent(inout) :: dets (nel, ndets)
        integer :: i, open_pos, tmp_dets (nel, ndets), nclosed, npos

        tmp_dets = dets
        nclosed = nel - nopen
        npos = 1
        open_pos = nclosed + 1
        ! Closed e- listed in pairs --> can jump pairs
        do i=1,nclosed-1,2
            do while (open_pos .le. nel)
                if (tmp_dets(Open_pos, 1) >= tmp_dets(i, 1)) exit
                dets(npos,:) = tmp_dets(open_pos,:)
                open_pos = open_pos + 1
                npos = npos + 1
            enddo
            dets(npos:npos+1,:) = tmp_dets(i:i+1,:)
            npos = npos + 2
        enddo
        ! Any remaining open electrons will be untouched.
    end subroutine

end module
