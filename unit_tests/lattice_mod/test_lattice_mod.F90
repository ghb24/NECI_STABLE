#include "macros.h"

! test my lattice module here, if the initializations and getter and
! setter functions work as intended

program test_lattice_mod

    use constants, only: dp, pi
    use lattice_mod
    use fruit
    use sort_mod, only: sort
    use SystemData, only: t_trans_corr_2body, t_k_space_hubbard

    implicit none

    integer :: failed_count

    call init_fruit()
    ! run my tests:
    call lattice_mod_test_driver
    call fruit_summary
    call fruit_finalize

    call get_failed_count(failed_count)
    if (failed_count /= 0) stop -1

contains

    subroutine lattice_mod_test_driver

        call run_test_case(test_init_lattice_chain, "test_init_lattice_chain")
        call run_test_case(test_init_lattice_rect, "test_init_lattice_rect")
        call run_test_case(test_sort_unique, "test_sort_unique")
        call run_test_case(test_init_lattice_tilted, "test_init_lattice_tilted")
        call run_test_case(test_init_lattice_tilted_transcorr, "test_init_lattice_tilted_transcorr")
        call run_test_case(test_init_lattice_triangular, "test_init_lattice_triangular")
        call run_test_case(test_init_lattice_hexagonal, "test_init_lattice_hexagonal")
        call run_test_case(test_init_lattice_kagome, "test_init_lattice_kagome")
        call run_test_case(inside_bz_2d_test, "inside_bz_2d_test")
        call run_test_case(on_line_2d_test, "on_line_2d_test")
        call run_test_case(apply_pbc_test, "apply_pbc_test")
        call run_test_case(apply_pbc_tilted_test, "apply_pbc_tilted_test")

    end subroutine lattice_mod_test_driver

    subroutine apply_pbc_test

        print *, ""
        print *, "testing: apply_pbc"
        print *, "TODO!"

    end subroutine apply_pbc_test

    subroutine on_line_2d_test

        print *, ""
        print *, "testing: on_line_2d "

        call assert_true(on_line_2d([1,1],[0,0],[1,1]))
        call assert_true(on_line_2d([2,2],[1,1],[0,0]))
        call assert_true(.not. on_line_2d([1,0],[0,0],[1,1]))

    end subroutine on_line_2d_test

    subroutine inside_bz_2d_test

        print *, ""
        print *, "testing: inside_bz_2d "

        call assert_true(inside_bz_2d(0,0, [-1,1],[-1,-1],[1,-1],[1,1]))
        call assert_true(inside_bz_2d(1,0, [-1,1],[-1,-1],[1,-1],[1,1]))
        call assert_true(inside_bz_2d(0,1, [-1,1],[-1,-1],[1,-1],[1,1]))
        call assert_true(inside_bz_2d(1,1, [-1,1],[-1,-1],[1,-1],[1,1]))

        call assert_true(inside_bz_2d(0,-1, [-1,1],[-1,-1],[1,-1],[1,1]))
        call assert_true(inside_bz_2d(1,-1, [-1,1],[-1,-1],[1,-1],[1,1]))
        call assert_true(.not.inside_bz_2d(2,-1, [-1,1],[-1,-1],[1,-1],[1,1]))
        call assert_true(.not.inside_bz_2d(2,0, [-1,1],[-1,-1],[1,-1],[1,1]))
        call assert_true(.not.inside_bz_2d(2,1, [-1,1],[-1,-1],[1,-1],[1,1]))
        call assert_true(.not.inside_bz_2d(3,2, [-1,1],[-1,-1],[1,-1],[1,1]))

        call assert_true(.not.inside_bz_2d(0,0,[3,3],[3,2],[4,2],[4,4]))

        call assert_true(inside_bz_2d(0,0, [-5,0],[0,-3],[3,0],[-2,3]))
        call assert_true(inside_bz_2d(-3,-1, [-5,0],[0,-3],[3,0],[-2,3]))
        call assert_true(inside_bz_2d(-1,2, [-5,0],[0,-3],[3,0],[-2,3]))
        call assert_true(.not.inside_bz_2d(0,3, [-5,0],[0,-3],[3,0],[-2,3]))


    end subroutine inside_bz_2d_test

    subroutine apply_pbc_tilted_test

        print *, ""
        print *, "testing: apply_pbc_tilted "
        print *, "TODO!"
    end subroutine apply_pbc_tilted_test

    subroutine test_sort_unique

        print *, ""
        print *, "testing sort_unique function"
        call assert_equals([1,2,3,4], sort_unique([3,2,1,4]), 4)
        call assert_equals([1,3,4], sort_unique([3,3,1,4,4]), 3)
        call assert_equals([-2,-1,0,1,2], sort_unique([2,1,0,0,-1,-2]), 5)
        call assert_equals([1], sort_unique([1,1,1,1]), 1)

    end subroutine test_sort_unique

    subroutine test_init_lattice_tilted_transcorr
        class(lattice), pointer :: ptr
        integer :: i
        real(dp) :: x(24)

        print *, ""
        print *, "initialize a 2x2 tilted square lattice with PBC with Gutwiller Transcorr"
        t_trans_corr_2body = .true.
        t_k_space_hubbard = .true.
        ptr => lattice('tilted', 2, 2, 1, .true., .true., .true., 'k-space')

        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(8, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))
        call assert_equals([3,5,7,8], ptr%get_neighbors(1),4)
        call assert_equals([3,5,7,8], ptr%get_neighbors(2),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(3),4)
        call assert_equals([3,5,7,8], ptr%get_neighbors(4),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(5),4)
        call assert_equals([3,5,7,8], ptr%get_neighbors(6),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(7),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(8),4)

        call assert_equals([1,3,7,11], ptr%get_spinorb_neighbors(15),4)

        call assert_equals(4.0_dp, ptr%dispersion_rel([0,0,0]))
        call assert_equals(-4.0_dp, ptr%dispersion_rel([2,0,0]))
        call assert_equals(0.0_dp, ptr%dispersion_rel([1,1,0]),1.e-10)
        call assert_equals(0.0_dp, ptr%dispersion_rel([-1,0,0]),1.e-10)

        x(1:8) = [(-ptr%dispersion_rel_orb(i), i = 1, 8)]
        call sort(x(1:8))

        do i = 1, 8
            print *, "e(k): ", x(i)
        end do

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 3x3 tilted square lattice with PBC"

        ptr => lattice('tilted', 3, 3, 1, .true., .true., .true., 'k-space')
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(18, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))
        call assert_equals([3,10,14,18], ptr%get_neighbors(1),4)
        call assert_equals([3,6,14,17], ptr%get_neighbors(2),4)
        call assert_equals([1,2,4,7], ptr%get_neighbors(3),4)
        call assert_equals([3,8,10,15], ptr%get_neighbors(4),4)
        call assert_equals([6,10,17,18], ptr%get_neighbors(5),4)
        call assert_equals([2,5,7,11], ptr%get_neighbors(6),4)
        call assert_equals([3,6,8,12], ptr%get_neighbors(7),4)
        call assert_equals([4,7,9,13], ptr%get_neighbors(8),4)
        call assert_equals([1,5,9,16], ptr%get_neighbors(18),4)
        call assert_equals([7,11,13,16], ptr%get_neighbors(12),4)

        call assert_equals(4.0_dp, ptr%dispersion_rel([0,0,0]))
        call assert_equals(2.0_dp, ptr%dispersion_rel([-1,0,0]),1.e-10)
        call assert_equals(1.0_dp, ptr%dispersion_rel([-1,-1,0]),1.e-10)
        call assert_equals(-1.0_dp, ptr%dispersion_rel([2,1,0]),1.e-10)
        call assert_equals(-2.0_dp, ptr%dispersion_rel([0,2,0]),1.e-10)
        call assert_equals(-4.0_dp, ptr%dispersion_rel([3,0,0]))

        x(1:18) = [(-ptr%dispersion_rel_orb(i), i = 1,18)]
        call sort(x(1:18))
        do i = 1, 18
            print *, "e(k): ", x(i)
        end do


        t_trans_corr_2body = .false.
        t_k_space_hubbard = .false.

    end subroutine test_init_lattice_tilted_transcorr

    subroutine test_init_lattice_tilted
        class(lattice), pointer :: ptr
        integer :: i
        real(dp) :: x(24)

        print *, ""
        print *, "initialize a 2x2 tilted square lattice with PBC"

        ptr => lattice('tilted', 2, 2, 1, .true., .true., .true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(8, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))
        call assert_equals([3,5,7,8], ptr%get_neighbors(1),4)
        call assert_equals([3,5,7,8], ptr%get_neighbors(2),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(3),4)
        call assert_equals([3,5,7,8], ptr%get_neighbors(4),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(5),4)
        call assert_equals([3,5,7,8], ptr%get_neighbors(6),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(7),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(8),4)

        call assert_equals([1,3,7,11], ptr%get_spinorb_neighbors(15),4)

        call assert_equals(4.0_dp, ptr%dispersion_rel([0,0,0]))
        call assert_equals(-4.0_dp, ptr%dispersion_rel([2,0,0]))
        call assert_equals(0.0_dp, ptr%dispersion_rel([1,1,0]),1.e-10)
        call assert_equals(0.0_dp, ptr%dispersion_rel([-1,0,0]),1.e-10)


        x(1:8) = [(-ptr%dispersion_rel_orb(i), i = 1, 8)]
        call sort(x(1:8))

        do i = 1, 8
            print *, "e(k): ", x(i)
        end do

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 2x2 tilted square lattice with closed BC"

        ptr => lattice('tilted', 2, 2, 1, .false., .false., .false.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(8, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( .not.ptr%is_periodic())
        call assert_true( .not.ptr%is_periodic(1))
        call assert_true( .not.ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))
        call assert_equals([3], ptr%get_neighbors(1),1)
        call assert_equals([3,5], ptr%get_neighbors(2),2)
        call assert_equals([1,2,4,6], ptr%get_neighbors(3),4)
        call assert_equals([3,7], ptr%get_neighbors(4),2)
        call assert_equals([2,6], ptr%get_neighbors(5),2)
        call assert_equals([3,5,7,8], ptr%get_neighbors(6),4)
        call assert_equals([4,6], ptr%get_neighbors(7),2)
        call assert_equals([6], ptr%get_neighbors(8),1)

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 2x2 tilted square lattice with closed BC in (1,-1)"

        ptr => lattice('tilted', 2, 2, 1, .true., .false., .false.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(8, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( .not.ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( .not.ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))
        call assert_equals([3,5], ptr%get_neighbors(1),2)
        call assert_equals([3,5], ptr%get_neighbors(2),2)
        call assert_equals([1,2,4,6], ptr%get_neighbors(3),4)
        call assert_equals([3,5,7,8], ptr%get_neighbors(4),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(5),4)
        call assert_equals([3,5,7,8], ptr%get_neighbors(6),4)
        call assert_equals([4,6], ptr%get_neighbors(7),2)
        call assert_equals([4,6], ptr%get_neighbors(8),2)

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 2x2 tilted square lattice with closed BC in (1,1)"

        ptr => lattice('tilted', 2, 2, 1, .false., .true., .false.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(8, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( .not.ptr%is_periodic())
        call assert_true( .not.ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))
        call assert_equals([3,7], ptr%get_neighbors(1),2)
        call assert_equals([3,5,7,8], ptr%get_neighbors(2),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(3),2)
        call assert_equals([3,7], ptr%get_neighbors(4),2)
        call assert_equals([2,6], ptr%get_neighbors(5),2)
        call assert_equals([3,5,7,8], ptr%get_neighbors(6),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(7),4)
        call assert_equals([2,6], ptr%get_neighbors(8),2)

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 3x3 tilted square lattice with closed BC"

        ptr => lattice('tilted', 3, 3, 1, .false., .false., .false.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(18, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( .not.ptr%is_periodic())
        call assert_true( .not.ptr%is_periodic(1))
        call assert_true( .not.ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))
        call assert_equals([3], ptr%get_neighbors(1),1)
        call assert_equals([3,6], ptr%get_neighbors(2),2)
        call assert_equals([1,2,4,7], ptr%get_neighbors(3),4)
        call assert_equals([3,8], ptr%get_neighbors(4),2)
        call assert_equals([6,10], ptr%get_neighbors(5),2)
        call assert_equals([2,5,7,11], ptr%get_neighbors(6),4)
        call assert_equals([3,6,8,12], ptr%get_neighbors(7),4)
        call assert_equals([4,7,9,13], ptr%get_neighbors(8),4)
        call assert_equals([16], ptr%get_neighbors(18),1)

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 3x3 tilted square lattice with PBC"

        ptr => lattice('tilted', 3, 3, 1, .true., .true., .true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(18, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))
        call assert_equals([3,10,14,18], ptr%get_neighbors(1),4)
        call assert_equals([3,6,14,17], ptr%get_neighbors(2),4)
        call assert_equals([1,2,4,7], ptr%get_neighbors(3),4)
        call assert_equals([3,8,10,15], ptr%get_neighbors(4),4)
        call assert_equals([6,10,17,18], ptr%get_neighbors(5),4)
        call assert_equals([2,5,7,11], ptr%get_neighbors(6),4)
        call assert_equals([3,6,8,12], ptr%get_neighbors(7),4)
        call assert_equals([4,7,9,13], ptr%get_neighbors(8),4)
        call assert_equals([1,5,9,16], ptr%get_neighbors(18),4)
        call assert_equals([7,11,13,16], ptr%get_neighbors(12),4)

        call assert_equals(4.0_dp, ptr%dispersion_rel([0,0,0]))
        call assert_equals(2.0_dp, ptr%dispersion_rel([-1,0,0]),1.e-10)
        call assert_equals(1.0_dp, ptr%dispersion_rel([-1,-1,0]),1.e-10)
        call assert_equals(-1.0_dp, ptr%dispersion_rel([2,1,0]),1.e-10)
        call assert_equals(-2.0_dp, ptr%dispersion_rel([0,2,0]),1.e-10)
        call assert_equals(-4.0_dp, ptr%dispersion_rel([3,0,0]))

        x(1:18) = [(-ptr%dispersion_rel_orb(i), i = 1,18)]
        call sort(x(1:18))
        do i = 1, 18
            print *, "e(k): ", x(i)
        end do

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 3x3 tilted square lattice with closed BC in (x,-x)"

        ptr => lattice('tilted', 3, 3, 1, .true., .false., .false.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(18, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( .not.ptr%is_periodic())
        call assert_true(  ptr%is_periodic(1))
        call assert_true( .not.ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))
        call assert_equals([3,10], ptr%get_neighbors(1),2)
        call assert_equals([3,6], ptr%get_neighbors(2),2)
        call assert_equals([1,2,4,7], ptr%get_neighbors(3),4)
        call assert_equals([3,8,10,15], ptr%get_neighbors(4),4)
        call assert_equals([6,10], ptr%get_neighbors(5),2)
        call assert_equals([2,5,7,11], ptr%get_neighbors(6),4)
        call assert_equals([3,6,8,12], ptr%get_neighbors(7),4)
        call assert_equals([4,7,9,13], ptr%get_neighbors(8),4)
        call assert_equals([9,16], ptr%get_neighbors(18),2)
        call assert_equals([1,4,5,11], ptr%get_neighbors(10),4)

        call lattice_deconstructor(ptr)

!         print *, ""
!         print *, "test also rectangular tilted lattices now!"
!         print *, "iniitalize a 1x2 4 site tilted! in k-space!"
!         ptr => lattice('tilted', 1, 2, 1, .true.,.true.,.true., 'k-space')
!
!         call assert_equals(2, ptr%get_ndim())
!         ! this would be nice: how do i implement that??
!         call assert_equals(1, ptr%get_length(1))
!         call assert_equals(2, ptr%get_length(2))
!         call assert_equals(4, ptr%get_nsites())
!         call assert_equals(4, ptr%get_nconnect_max())
!         call assert_true( ptr%is_periodic())
!         call assert_true( ptr%is_periodic(1))
!         call assert_true( ptr%is_periodic(2))
!         call assert_equals(1, ptr%get_site_index(1))
!         call assert_equals(2, ptr%get_site_index(2))
!         call assert_equals(3, ptr%get_site_index(3))
!         call assert_equals(4, ptr%get_site_index(4))
!         call assert_equals([2,4], ptr%get_neighbors(1),2)
!         call assert_equals([1,3], ptr%get_neighbors(2),2)
!         call assert_equals([2,4], ptr%get_neighbors(3),2)
!         print *, "neigh:1 ", ptr%get_neighbors(1)
!         print *, "neigh:2 ", ptr%get_neighbors(2)
!         print *, "neigh:3 ", ptr%get_neighbors(3)
!         print *, "neigh:4 ", ptr%get_neighbors(4)
!         call assert_equals([1,3], ptr%get_neighbors(4),2)
        ! also test the dispersion relation:
        ! i should also write a routine for the dispersion relation, where
        ! i just input the orbital(spin or spatial?) and which internally
        ! takes the correct k-vector  todo!
        ! i should also internally use real-space and k-vectors to better
        ! deal with periodicity and stuff!
!         call assert_equals([2.0_dp,0.0_dp], ptr%get_lat_vec(1),2)
!         call assert_equals([2.0_dp,-2.0_dp], ptr%get_lat_vec(2),2)
!         ! and from this i can calculate the k-vectors by r_i*k_j = 2pi\delta_{ij}
!         ! todo!
!         call assert_equals([0.0_dp,0.0_dp], ptr%get_rec_vec(1),2)
!         call assert_equals([0.0_dp,0.0_dp], ptr%get_rec_vec(2),2)
!
!         call assert_equals(0.0_dp, ptr%dispersion_rel(1))
!
!         print *, "initialize a 2x3 12-site tilted lattice: "
!         ptr => lattice('tilted', 2,3,1,.true.,.true.,.true.,'k-space')
!         call assert_equals(2, ptr%get_ndim())
!         ! this would be nice: how do i implement that??
!         call assert_equals(2, ptr%get_length(1))
!         call assert_equals(3, ptr%get_length(2))
!         call assert_equals(12, ptr%get_nsites())
!         call assert_equals(4, ptr%get_nconnect_max())
!         call assert_true( ptr%is_periodic())
!         call assert_true( ptr%is_periodic(1))
!         call assert_true( ptr%is_periodic(2))
!         call assert_equals(1, ptr%get_site_index(1))
!         call assert_equals(2, ptr%get_site_index(2))
!         call assert_equals(3, ptr%get_site_index(3))
!         call assert_equals(4, ptr%get_site_index(4))
!         call assert_equals([3,5,11,12], ptr%get_neighbors(1),2)
!         call assert_equals([3,5,11,12], ptr%get_neighbors(2),2)
!         call assert_equals([1,2,4,6], ptr%get_neighbors(3),2)
!         call assert_equals([3,5,7,9], ptr%get_neighbors(4),4)
!         call assert_equals([1,2,4,6], ptr%get_neighbors(5),4)
!         call assert_equals([3,5,7,9], ptr%get_neighbors(6),4)
!         call assert_equals([4,6,8,10], ptr%get_neighbors(7),4)
!         call assert_equals([7,9,11,12], ptr%get_neighbors(8),4)
!         call assert_equals([1,2,8,10], ptr%get_neighbors(12),4)
!
!         print *, "initialize a 3x2 12-site tilted lattice: "
!         print *, "due to symmetr the 2x3 is treated the dealt internally as"
!         print *, "the already working 2x3!"
!         ptr => lattice('tilted', 3,2,1,.true.,.true.,.true.,'k-space')
!         call assert_equals(2, ptr%get_ndim())
!         ! this would be nice: how do i implement that??
!         call assert_equals(2, ptr%get_length(1))
!         call assert_equals(3, ptr%get_length(2))
!         call assert_equals(12, ptr%get_nsites())
!         call assert_equals(4, ptr%get_nconnect_max())
!         call assert_true( ptr%is_periodic())
!         call assert_true( ptr%is_periodic(1))
!         call assert_true( ptr%is_periodic(2))
!         call assert_equals(1, ptr%get_site_index(1))
!         call assert_equals(2, ptr%get_site_index(2))
!         call assert_equals(3, ptr%get_site_index(3))
!         call assert_equals(4, ptr%get_site_index(4))
!         call assert_equals([3,5,11,12], ptr%get_neighbors(1),4)
!         call assert_equals([3,5,11,12], ptr%get_neighbors(2),4)
!         call assert_equals([1,2,4,6], ptr%get_neighbors(3),4)
!         call assert_equals([3,5,7,9], ptr%get_neighbors(4),4)
!         call assert_equals([1,2,4,6], ptr%get_neighbors(5),4)
!         call assert_equals([3,5,7,9], ptr%get_neighbors(6),4)
!         call assert_equals([4,6,8,10], ptr%get_neighbors(7),4)
!         call assert_equals([7,9,11,12], ptr%get_neighbors(8),4)
!         call assert_equals([1,2,8,10], ptr%get_neighbors(12),4)
!
!         print *, "initialize a 3x4 24-site tilted lattice: "
!         ptr => lattice('tilted', 3,4,1,.true.,.true.,.true.,'k-space')
!         call assert_equals(2, ptr%get_ndim())
!         ! this would be nice: how do i implement that??
!         call assert_equals(3, ptr%get_length(1))
!         call assert_equals(4, ptr%get_length(2))
!         call assert_equals(24, ptr%get_nsites())
!         call assert_equals(4, ptr%get_nconnect_max())
!         call assert_true( ptr%is_periodic())
!         call assert_true( ptr%is_periodic(1))
!         call assert_true( ptr%is_periodic(2))
!         call assert_equals(1, ptr%get_site_index(1))
!         call assert_equals(2, ptr%get_site_index(2))
!         call assert_equals(3, ptr%get_site_index(3))
!         call assert_equals(4, ptr%get_site_index(4))
!         call assert_equals([3,10,20,24], ptr%get_neighbors(1),4)
!         call assert_equals([3,6,20,23], ptr%get_neighbors(2),4)
!         call assert_equals([1,2,4,7], ptr%get_neighbors(3),4)
!         call assert_equals([3,8,10,16], ptr%get_neighbors(4),4)
!         call assert_equals([6,10,23,24], ptr%get_neighbors(5),4)
!         call assert_equals([2,5,7,11], ptr%get_neighbors(6),4)
!         call assert_equals([3,6,8,12], ptr%get_neighbors(7),4)
!         call assert_equals([13,17,19,22], ptr%get_neighbors(18),4)
!         call assert_equals([1,5,15,22], ptr%get_neighbors(24),4)
!
!         do i = 1, 24
!             print *, "k, e(k): ", ptr%get_k_vec(i), ptr%dispersion_rel_orb(i)
!         end do
!
!         x = [(-ptr%dispersion_rel_orb(i), i = 1,24)]
!         call sort(x)
! !         x = x(24:1:-1)
!
!         do i = 1,24
!             print *, x(i)
!         end do
!
!         print *, "initialize a 2x4 16-site tilted lattice: "
!         ptr => lattice('tilted', 2,4,1,.true.,.true.,.true.,'k-space')
!         call assert_equals(2, ptr%get_ndim())
!         ! this would be nice: how do i implement that??
!         call assert_equals(2, ptr%get_length(1))
!         call assert_equals(4, ptr%get_length(2))
!         call assert_equals(16, ptr%get_nsites())
!         call assert_equals(4, ptr%get_nconnect_max())
!         call assert_true( ptr%is_periodic())
!         call assert_true( ptr%is_periodic(1))
!         call assert_true( ptr%is_periodic(2))
!         call assert_equals(1, ptr%get_site_index(1))
!         call assert_equals(2, ptr%get_site_index(2))
!         call assert_equals(3, ptr%get_site_index(3))
!         call assert_equals(4, ptr%get_site_index(4))
!         call assert_equals([3,5,15,16], ptr%get_neighbors(1),4)
!         call assert_equals([3,5,15,16], ptr%get_neighbors(2),4)
!         call assert_equals([1,2,4,6], ptr%get_neighbors(3),4)
!         call assert_equals([3,5,7,9], ptr%get_neighbors(4),4)
!         call assert_equals([1,2,4,6], ptr%get_neighbors(5),4)
!         call assert_equals([3,5,7,9], ptr%get_neighbors(6),4)
!         call assert_equals([1,2,13,16], ptr%get_neighbors(16),4)
!         call assert_equals(0.0_dp, ptr%dispersion_rel_orb(1))
!
!         call stop_all("here","now")
!         call lattice_deconstructor(ptr)
!         print *, "iniitalize a 2x1 4 site tilted! in k-space!"
!         ptr => lattice('tilted', 2, 1, 1, .true.,.true.,.true., 'k-space')
!
!         call assert_equals(2, ptr%get_ndim())
!         ! this would be nice: how do i implement that??
!         call assert_equals(2, ptr%get_length(1))
!         call assert_equals(1, ptr%get_length(2))
!         call assert_equals(4, ptr%get_nsites())
!         call assert_equals(2, ptr%get_nconnect_max())
!         call assert_true( ptr%is_periodic())
!         call assert_true( ptr%is_periodic(1))
!         call assert_true( ptr%is_periodic(2))
!         call assert_equals(1, ptr%get_site_index(1))
!         call assert_equals(2, ptr%get_site_index(2))
!         call assert_equals(3, ptr%get_site_index(3))
!         call assert_equals(4, ptr%get_site_index(4))
!         call assert_equals([3,4], ptr%get_neighbors(1),2)
!         call assert_equals([3,4], ptr%get_neighbors(2),2)
!         call assert_equals([1,2], ptr%get_neighbors(3),2)
!         call assert_equals([1,2], ptr%get_neighbors(4),2)
!         ! also test the dispersion relation:
        ! i should also write a routine for the dispersion relation, where
        ! i just input the orbital(spin or spatial?) and which internally
        ! takes the correct k-vector  todo!
!
!         ! also test for the lattice vector and reciprocal vectors
!         call assert_equals([2.0_dp,2.0_dp], ptr%get_lat_vec(1),2)
!         call assert_equals([2.0_dp,0.0_dp], ptr%get_lat_vec(2),2)
!
!         call assert_equals([0.0_dp,0.0_dp], ptr%get_rec_vec(1),2)
!         call assert_equals([0.0_dp,0.0_dp], ptr%get_rec_vec(2),2)


!         call lattice_deconstructor(ptr)






    end subroutine test_init_lattice_tilted

    subroutine test_init_lattice_triangular
        class(lattice), pointer :: ptr

        print *, ""
        print *, "initialize a 2x2 triangular lattice with periodic boundary conditions"
        ptr => lattice('triangle',2,2,1,.true.,.true.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(4, ptr%get_nsites())
        call assert_equals(6, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))

        call assert_equals([2,3,4], ptr%get_neighbors(1),3)
        call assert_equals([1,3,4], ptr%get_neighbors(2),3)
        call assert_equals([1,2,4], ptr%get_neighbors(3),3)
        call assert_equals([1,2,3], ptr%get_neighbors(4),3)

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 2x3 triangular lattice with PBC"
        ptr => lattice('triangle',2,3,1,.true.,.true.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(6, ptr%get_nsites())
        call assert_equals(6, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))

        call assert_equals([2,3,4,5,6], ptr%get_neighbors(1),5)
        call assert_equals([1,3,4,5,6], ptr%get_neighbors(2),5)
        call assert_equals([1,2,4,5,6], ptr%get_neighbors(3),5)
        call assert_equals([1,2,3,4,5], ptr%get_neighbors(6),5)

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 3x2 triangular lattice with PBC"
        ptr => lattice('triangle',3,2,1,.true.,.true.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(6, ptr%get_nsites())
        call assert_equals(6, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))

        call assert_equals([2,3,4,5,6], ptr%get_neighbors(1),5)
        call assert_equals([1,3,4,5,6], ptr%get_neighbors(2),5)
        call assert_equals([1,2,4,5,6], ptr%get_neighbors(3),5)
        call assert_equals([1,2,3,4,5], ptr%get_neighbors(6),5)

        call lattice_deconstructor(ptr)


        print *, ""
        print *, "initialize a 3x3 triangular lattice with PBC"
        ptr => lattice('triangle',3,3,1,.true.,.true.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(9, ptr%get_nsites())
        call assert_equals(6, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))

        call assert_equals([2,3,4,5,7,9], ptr%get_neighbors(1),6)
        call assert_equals([1,3,5,6,7,8], ptr%get_neighbors(2),6)
        call assert_equals([1,2,4,6,8,9], ptr%get_neighbors(3),6)
        call assert_equals([1,2,4,6,8,9], ptr%get_neighbors(5),6)
        call assert_equals([1,3,5,6,7,8], ptr%get_neighbors(9),6)

        call lattice_deconstructor(ptr)


    end subroutine test_init_lattice_triangular

    subroutine test_init_lattice_rect

        class(lattice), pointer :: ptr

        print *, ""
        print *, "initialize a 2x2 square lattice with periodic boundary conditions"

        ptr => lattice('square',2,2,1,.true.,.true.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(4, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))

        call assert_equals([2,3], ptr%get_neighbors(1),2)
        call assert_equals([1,4], ptr%get_neighbors(2),2)
        call assert_equals([1,4], ptr%get_neighbors(3),2)
        call assert_equals([2,3], ptr%get_neighbors(4),2)

        call assert_equals(2, ptr%get_num_neighbors(1))
        call assert_equals(2, ptr%get_num_neighbors(2))
        call assert_equals(2, ptr%get_num_neighbors(3))
        call assert_equals(2, ptr%get_num_neighbors(4))
        call assert_equals(4.0_dp, ptr%dispersion_rel([0,0,0]))
        call assert_equals(0.0_dp, ptr%dispersion_rel([0,1,0]))
        call assert_equals(0.0_dp, ptr%dispersion_rel([1,0,0]))
        call assert_equals(-4.0_dp, ptr%dispersion_rel([1,1,0]))

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 3x3 square lattice with periodic boundary conditions"

        ptr => lattice('square',3,3,1,.true.,.true.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(9, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))
        call assert_equals(9, ptr%get_site_index(9))

        call assert_equals([2,3,4,7], ptr%get_neighbors(1),4)
        call assert_equals([1,3,5,8], ptr%get_neighbors(2),4)
        call assert_equals([1,2,6,9], ptr%get_neighbors(3),4)
        call assert_equals([1,5,6,7], ptr%get_neighbors(4),4)

        call assert_equals(4, ptr%get_num_neighbors(1))
        call assert_equals(4, ptr%get_num_neighbors(2))
        call assert_equals(4, ptr%get_num_neighbors(3))
        call assert_equals(4, ptr%get_num_neighbors(4))

        call assert_equals(4.0_dp, ptr%dispersion_rel([0,0,0]))
        call assert_equals(2.0_dp*(1.0_dp + cos(2*pi/3)), ptr%dispersion_rel([1,0,0]))
        call assert_equals(2.0_dp*(1.0_dp + cos(2*pi/3)), ptr%dispersion_rel([0,1,0]))


        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 2x2 square lattice with closed boundary conditions in y"

        ptr => lattice('square',2,2,1,.true.,.false.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(4, ptr%get_nsites())
        call assert_equals(3, ptr%get_nconnect_max())
        call assert_true( .not.ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( .not.ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))

        call assert_equals([2,3], ptr%get_neighbors(1),2)
        call assert_equals([1,4], ptr%get_neighbors(2),2)
        call assert_equals([1,4], ptr%get_neighbors(3),2)
        call assert_equals([2,3], ptr%get_neighbors(4),2)

        call assert_equals(2, ptr%get_num_neighbors(1))
        call assert_equals(2, ptr%get_num_neighbors(2))
        call assert_equals(2, ptr%get_num_neighbors(3))
        call assert_equals(2, ptr%get_num_neighbors(4))

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 2x2 square lattice with closed boundary conditions in x"

        ptr => lattice('square',2,2,1,.false.,.true.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(4, ptr%get_nsites())
        call assert_equals(3, ptr%get_nconnect_max())
        call assert_true( .not.ptr%is_periodic())
        call assert_true( .not.ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))

        call assert_equals([2,3], ptr%get_neighbors(1),2)
        call assert_equals([1,4], ptr%get_neighbors(2),2)
        call assert_equals([1,4], ptr%get_neighbors(3),2)
        call assert_equals([2,3], ptr%get_neighbors(4),2)

        call assert_equals(2, ptr%get_num_neighbors(1))
        call assert_equals(2, ptr%get_num_neighbors(2))
        call assert_equals(2, ptr%get_num_neighbors(3))
        call assert_equals(2, ptr%get_num_neighbors(4))

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 2x2 square lattice with closed boundary conditions "

        ptr => lattice('square',2,2,1,.false.,.false.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(4, ptr%get_nsites())
        call assert_equals(2, ptr%get_nconnect_max())
        call assert_true( .not.ptr%is_periodic())
        call assert_true( .not.ptr%is_periodic(1))
        call assert_true( .not.ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))

        call assert_equals([2,3], ptr%get_neighbors(1),2)
        call assert_equals([1,4], ptr%get_neighbors(2),2)
        call assert_equals([1,4], ptr%get_neighbors(3),2)
        call assert_equals([2,3], ptr%get_neighbors(4),2)

        call assert_equals(2, ptr%get_num_neighbors(1))
        call assert_equals(2, ptr%get_num_neighbors(2))
        call assert_equals(2, ptr%get_num_neighbors(3))
        call assert_equals(2, ptr%get_num_neighbors(4))

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 3x3 square lattice with closed boundary conditions in y"

        ptr => lattice('square',3,3,1,.true.,.false.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(9, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( .not.ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( .not.ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))
        call assert_equals(9, ptr%get_site_index(9))

        call assert_equals([2,4,7], ptr%get_neighbors(1),3)
        call assert_equals([1,3,5,8], ptr%get_neighbors(2),4)
        call assert_equals([2,6,9], ptr%get_neighbors(3),3)
        call assert_equals([1,5,7], ptr%get_neighbors(4),3)
        call assert_equals([2,4,6,8], ptr%get_neighbors(5),4)

        call assert_equals(3, ptr%get_num_neighbors(1))
        call assert_equals(4, ptr%get_num_neighbors(2))
        call assert_equals(3, ptr%get_num_neighbors(3))
        call assert_equals(3, ptr%get_num_neighbors(4))
        call assert_equals(4, ptr%get_num_neighbors(5))

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 3x3 square lattice with closed boundary conditions in x"

        ptr => lattice('square',3,3,1,.false.,.true.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(9, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( .not.ptr%is_periodic())
        call assert_true( .not.ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))
        call assert_equals(9, ptr%get_site_index(9))

        call assert_equals([2,3,4], ptr%get_neighbors(1),3)
        call assert_equals([1,3,5], ptr%get_neighbors(2),3)
        call assert_equals([1,2,6], ptr%get_neighbors(3),3)
        call assert_equals([1,5,6,7], ptr%get_neighbors(4),4)
        call assert_equals([2,4,6,8], ptr%get_neighbors(5),4)
        call assert_equals([6,7,8], ptr%get_neighbors(9),3)

        call assert_equals(3, ptr%get_num_neighbors(1))
        call assert_equals(3, ptr%get_num_neighbors(2))
        call assert_equals(3, ptr%get_num_neighbors(3))
        call assert_equals(4, ptr%get_num_neighbors(4))
        call assert_equals(4, ptr%get_num_neighbors(5))
        call assert_equals(3, ptr%get_num_neighbors(9))

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 3x3 square lattice with closed boundary conditions"

        ptr => lattice('square',3,3,1,.false.,.false.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(9, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( .not.ptr%is_periodic())
        call assert_true( .not.ptr%is_periodic(1))
        call assert_true( .not.ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))
        call assert_equals(9, ptr%get_site_index(9))

        call assert_equals([2,4], ptr%get_neighbors(1),2)
        call assert_equals([1,3,5], ptr%get_neighbors(2),3)
        call assert_equals([2,6], ptr%get_neighbors(3),2)
        call assert_equals([1,5,7], ptr%get_neighbors(4),3)
        call assert_equals([2,4,6,8], ptr%get_neighbors(5),4)
        call assert_equals([6,8], ptr%get_neighbors(9),2)

        call assert_equals(2, ptr%get_num_neighbors(1))
        call assert_equals(3, ptr%get_num_neighbors(2))
        call assert_equals(2, ptr%get_num_neighbors(3))
        call assert_equals(3, ptr%get_num_neighbors(4))
        call assert_equals(4, ptr%get_num_neighbors(5))
        call assert_equals(2, ptr%get_num_neighbors(9))

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 3x2 square lattice with closed boundary conditions"

        ptr => lattice('rectangle',3,2,1,.false.,.false.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(6, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( .not.ptr%is_periodic())
        call assert_true( .not.ptr%is_periodic(1))
        call assert_true( .not.ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))
        call assert_equals(6, ptr%get_site_index(6))

        call assert_equals([2,4], ptr%get_neighbors(1),2)
        call assert_equals([1,3,5], ptr%get_neighbors(2),3)
        call assert_equals([2,6], ptr%get_neighbors(3),2)
        call assert_equals([1,5], ptr%get_neighbors(4),2)
        call assert_equals([2,4,6], ptr%get_neighbors(5),3)

        call assert_equals(2, ptr%get_num_neighbors(1))
        call assert_equals(3, ptr%get_num_neighbors(2))
        call assert_equals(2, ptr%get_num_neighbors(3))
        call assert_equals(2, ptr%get_num_neighbors(4))
        call assert_equals(3, ptr%get_num_neighbors(5))

        call lattice_deconstructor(ptr)


        print *, ""
        print *, "initialize a 2x3 square lattice with closed boundary conditions in x"

        ptr => lattice('rectangle',2,3,1,.false.,.true.,.true.)
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(6, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( .not.ptr%is_periodic())
        call assert_true( .not.ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))
        call assert_equals(6, ptr%get_site_index(6))

        call assert_equals([2,3], ptr%get_neighbors(1),2)
        call assert_equals([1,4], ptr%get_neighbors(2),2)
        call assert_equals([1,4,5], ptr%get_neighbors(3),3)
        call assert_equals([2,3,6], ptr%get_neighbors(4),3)
        call assert_equals([3,6], ptr%get_neighbors(5),2)

        call assert_equals(2, ptr%get_num_neighbors(1))
        call assert_equals(2, ptr%get_num_neighbors(2))
        call assert_equals(3, ptr%get_num_neighbors(3))
        call assert_equals(3, ptr%get_num_neighbors(4))
        call assert_equals(2, ptr%get_num_neighbors(5))

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 4x2 square lattice with pbc"
        ptr => lattice('rectangle', 4, 2, 1, .true., .true., .true.)
        call assert_equals(4, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(8, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))
        call assert_equals(6, ptr%get_site_index(6))

        call assert_equals([2,4,5], ptr%get_neighbors(1),3)
        call assert_equals([1,3,6], ptr%get_neighbors(2),3)
        call assert_equals([2,4,7], ptr%get_neighbors(3),3)
        call assert_equals([1,3,8], ptr%get_neighbors(4),3)
        call assert_equals([1,6,8], ptr%get_neighbors(5),3)

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 4x3 square lattice with pbc"
        ptr => lattice('rectangle', 4, 3, 1, .true., .true., .true.)
        call assert_equals(4, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(12, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))
        call assert_equals(6, ptr%get_site_index(6))

        call assert_equals([2,4,5,9], ptr%get_neighbors(1),4)
        call assert_equals([1,3,6,10], ptr%get_neighbors(2),4)
        call assert_equals([2,4,7,11], ptr%get_neighbors(3),4)
        call assert_equals([1,3,8,12], ptr%get_neighbors(4),4)
        call assert_equals([1,6,8,9], ptr%get_neighbors(5),4)

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 2x4 square lattice with pbc"
        ptr => lattice('rectangle', 2, 4, 1, .true., .true., .true.)
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(4, ptr%get_length(2))
        call assert_equals(8, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))
        call assert_equals(6, ptr%get_site_index(6))

        call assert_equals([2,3,7], ptr%get_neighbors(1),3)
        call assert_equals([1,4,8], ptr%get_neighbors(2),3)
        call assert_equals([1,4,5], ptr%get_neighbors(3),3)
        call assert_equals([2,3,6], ptr%get_neighbors(4),3)
        call assert_equals([3,6,7], ptr%get_neighbors(5),3)

        call lattice_deconstructor(ptr)

    end subroutine test_init_lattice_rect

    subroutine test_init_lattice_chain
        ! test the specific initializers
        ! implicitly test the setter and getter functions or?
        class(lattice), pointer :: ptr

        print *, "initialize a periodic one-site 'chain' lattice: "
        ptr => lattice('chain', 1, 1, 1, .true., .true., .true.)

        call assert_equals(ptr%get_ndim(), 1)
        call assert_equals(ptr%get_nsites(), 1)
        call assert_equals(ptr%get_length(), 1)
        call assert_true(ptr%is_periodic())
        ! hm.. for this case the n_connect_max is not 2. since it is a special
        ! case.. for the 1 site and the 2 site non-periodic it is only 1 or 0
        ! but this special cases shouldnt matter.. or?
        ! leave it for now, so i am reminded that i might have to adjust that!
        call assert_equals(0, ptr%get_nconnect_max())

        ! todo: it would be really fancy to overload the round brackets
        ! to give the get function of the index value of our lattice!:
        ! for a lattice of length 1 we need some special initialization..
        ! thats good to have these edge cases!
        ! i need a public getter for the site indices..
        ! do i want to put it into the site type or in the lattice type?
        ! ptr%get_index(1) or ptr%sites(1)%get_index
        ! in the first i could check if the index is too high and i would
        ! not have to make so much public..
        call assert_equals(ptr%get_site_index(1), 1)
        ! and i want to have a get_neighbors routine
        ! apparently there is no assert equal for vectors of ints??
        ! thats BS! there is, but one has to give the additional number
        ! of elements input!
        call assert_equals(ptr%get_neighbors(1), [-1], 1)

        call assert_equals(2.0_dp, ptr%dispersion_rel([0,0,0]))
        call assert_equals(2.0_dp, ptr%dispersion_rel([1,0,0]))

        call lattice_deconstructor(ptr)

        call assert_true(.not.associated(ptr))

        print *, ""
        print *, "initialize a non-periodic two-site 'chain' lattice: "
        ptr => lattice('chain', 2, 1, 1, .false., .false., .false.)
        call assert_equals(ptr%get_ndim(), 1)
        call assert_equals(ptr%get_nsites(), 2)
        call assert_equals(ptr%get_length(), 2)
        call assert_true(.not. ptr%is_periodic())
        call assert_equals(1, ptr%get_nconnect_max())

        call assert_equals(ptr%get_site_index(1), 1)
        call assert_equals(ptr%get_site_index(2), 2)
        call assert_equals(ptr%get_neighbors(1), [2], size(ptr%get_neighbors(1)))
        call assert_equals(ptr%get_neighbors(2), [1], size(ptr%get_neighbors(2)))

        call assert_equals(1, ptr%get_num_neighbors(1))
        call assert_equals(1, ptr%get_num_neighbors(2))

        call lattice_deconstructor(ptr)

    call assert_true(.not. associated(ptr))

        print *, ""
        print *, "initialize a periodic 100 site 'chain' lattice: "
        ptr => lattice('chain', 0, 100, 1, .true., .true., .true.)
        call assert_equals(ptr%get_ndim(), 1)
        call assert_equals(ptr%get_nsites(), 100)
        call assert_equals(ptr%get_length(), 100)
        call assert_true(ptr%is_periodic())
        call assert_equals(2, ptr%get_nconnect_max())

        call assert_equals(ptr%get_site_index(1), 1)
        call assert_equals(ptr%get_site_index(100), 100)
        call assert_equals(ptr%get_site_index(50), 50)
        call assert_equals(ptr%get_neighbors(1), [100, 2], size(ptr%get_neighbors(1)))
        call assert_equals(ptr%get_neighbors(2), [1,3], size(ptr%get_neighbors(2)))
        call assert_equals(ptr%get_neighbors(100), [99,1], size(ptr%get_neighbors(2)))

        call assert_equals(2, ptr%get_num_neighbors(1))
        call assert_equals(2, ptr%get_num_neighbors(2))
        call assert_equals(2, ptr%get_num_neighbors(100))
        call assert_equals(2.0_dp, ptr%dispersion_rel([0,0,0]))
        call assert_equals(2.0_dp*cos(2*pi/100), ptr%dispersion_rel([1,0,0]))
        ! i actually do not need to have the common lattice type or?
        ! when i want to test a chain i could just use a chain or?

        call lattice_deconstructor(ptr)

        call assert_true(.not.associated(ptr))

    end subroutine test_init_lattice_chain

    subroutine test_init_lattice_hexagonal

        class(lattice), pointer :: ptr

        print *, ""
        print *, "testing: init_lattice_hexagonal"
        ptr => lattice('hexagonal', 1, 1, 1, .true.,.true.,.true.)

        call assert_equals(2, ptr%get_ndim())
        call assert_equals(8, ptr%get_nsites())
        call assert_equals(1, ptr%get_length(1))
        call assert_equals(1, ptr%get_length(2))
        call assert_true(ptr%is_periodic())
        call assert_equals(3, ptr%get_nconnect_max())

        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(8, ptr%get_site_index(8))
        call assert_equals(3, ptr%get_num_neighbors(1))
        call assert_equals(3, ptr%get_num_neighbors(8))

        call assert_equals([2,4,5], ptr%get_neighbors(1),3)
        call assert_equals([1,3,6], ptr%get_neighbors(2),3)
        call assert_equals([2,4,7], ptr%get_neighbors(3),3)
        call assert_equals([1,3,8], ptr%get_neighbors(4),3)
        call assert_equals([1,6,8], ptr%get_neighbors(5),3)
        call assert_equals([2,5,7], ptr%get_neighbors(6),3)
        call assert_equals([3,6,8], ptr%get_neighbors(7),3)
        call assert_equals([4,5,7], ptr%get_neighbors(8),3)

        call lattice_deconstructor(ptr)
        call assert_true(.not. associated(ptr))

        ptr => lattice('hexagonal', 2, 1,1,.true.,.true.,.true.)

        call assert_equals(2, ptr%get_ndim())
        call assert_equals(16, ptr%get_nsites())
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(1, ptr%get_length(2))
        call assert_true(ptr%is_periodic())
        call assert_equals(3, ptr%get_nconnect_max())

        call assert_equals([2,8,9], ptr%get_neighbors(1),3)
        call assert_equals([1,3,10], ptr%get_neighbors(2),3)
        call assert_equals([2,4,11], ptr%get_neighbors(3),3)
        call assert_equals([3,5,12], ptr%get_neighbors(4),3)
        call assert_equals([4,6,13], ptr%get_neighbors(5),3)
        call assert_equals([5,7,14], ptr%get_neighbors(6),3)
        call assert_equals([6,8,15], ptr%get_neighbors(7),3)
        call assert_equals([1,7,16], ptr%get_neighbors(8),3)

        ! the 1x2 is the first "new" lattice type since the
        ! Xx1 are like Xx2 square lattices with half the hopping t
        ptr => lattice('hexagonal', 1, 2,1,.true.,.true.,.true.)

        call assert_equals(2, ptr%get_ndim())
        call assert_equals(16, ptr%get_nsites())
        call assert_equals(1, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_true(ptr%is_periodic())
        call assert_equals(3, ptr%get_nconnect_max())

        call assert_equals([2,4,5], ptr%get_neighbors(1),3)
        call assert_equals([1,3,14], ptr%get_neighbors(2),3)
        call assert_equals([2,4,7], ptr%get_neighbors(3),3)
        call assert_equals([1,3,16], ptr%get_neighbors(4),3)
        call assert_equals([1,6,8], ptr%get_neighbors(5),3)
        call assert_equals([5,7,10], ptr%get_neighbors(6),3)
        call assert_equals([3,6,8], ptr%get_neighbors(7),3)
        call assert_equals([5,7,12], ptr%get_neighbors(8),3)

    end subroutine test_init_lattice_hexagonal

    subroutine test_init_lattice_kagome

        class(lattice), pointer :: ptr
        print *, ""
        print *, "testing: init_lattice_kagome:"

        ptr => lattice('kagome', 1,1,1,.true.,.true.,.true.)

        call assert_equals(2, ptr%get_ndim())
        call assert_equals(6, ptr%get_nsites())
        call assert_equals(1, ptr%get_length(1))
        call assert_equals(1, ptr%get_length(2))
        call assert_true(ptr%is_periodic())
        call assert_equals(4, ptr%get_nconnect_max())

        call assert_equals([2,4,6], ptr%get_neighbors(1),3)
        call assert_equals([1,3,4,5], ptr%get_neighbors(2),4)
        call assert_equals([2,5,6], ptr%get_neighbors(3),3)
        call assert_equals([1,2,6], ptr%get_neighbors(4),3)
        call assert_equals([2,3,6], ptr%get_neighbors(5),3)
        call assert_equals([1,3,4,5], ptr%get_neighbors(6),4)

        ptr => lattice('kagome', 2,1,1,.true.,.true.,.true.)

        call assert_equals(2, ptr%get_ndim())
        call assert_equals(12, ptr%get_nsites())
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(1, ptr%get_length(2))
        call assert_true(ptr%is_periodic())
        call assert_equals(4, ptr%get_nconnect_max())

        call assert_equals([3,5,7,10], ptr%get_neighbors(6),4)

        ptr => lattice('kagome', 1,2,1,.true.,.true.,.true.)

        call assert_equals(2, ptr%get_ndim())
        call assert_equals(12, ptr%get_nsites())
        call assert_equals(1, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_true(ptr%is_periodic())
        call assert_equals(4, ptr%get_nconnect_max())

        call assert_equals([1,3,10,11], ptr%get_neighbors(12),4)

        ptr => lattice('kagome', 2,2,1,.true.,.true.,.true.)

        call assert_equals(2, ptr%get_ndim())
        call assert_equals(24, ptr%get_nsites())
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_true(ptr%is_periodic())
        call assert_equals(4, ptr%get_nconnect_max())

        call assert_equals([1,10,15,23], ptr%get_neighbors(24),4)

    end subroutine test_init_lattice_kagome

end program test_lattice_mod
