#include "macros.h" 

! test my lattice module here, if the initializations and getter and 
! setter functions work as intended 

program test_lattice_mod 

    use constants, only: dp, pi
    use lattice_mod 
    use fruit 
    use sort_mod, only: sort

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
        call run_test_case(test_init_lattice_star, "test_init_lattice_star")
        call run_test_case(test_init_lattice_aim_chain, "test_init_lattice_aim_chain")
        call run_test_case(test_init_lattice_aim_star, "test_init_lattice_aim_star")
        call run_test_case(test_init_cluster_lattice_aim, "test_init_cluster_lattice_aim")
        call run_test_case(test_init_lattice_rect, "test_init_lattice_rect")
        call run_test_case(test_sort_unique, "test_sort_unique")
        call run_test_case(test_init_lattice_tilted, "test_init_lattice_tilted")
        call run_test_case(test_init_lattice_cube, "test_init_lattice_cube")
        call run_test_case(test_init_lattice_triangular, "test_init_lattice_triangular")
        call run_test_case(test_init_lattice_hexagonal, "test_init_lattice_hexagonal")
        call run_test_case(test_init_lattice_kagome, "test_init_lattice_kagome")
        call run_test_case(inside_bz_2d_test, "inside_bz_2d_test")
        call run_test_case(init_lattice_ole_test, "init_lattice_ole_test")
        call run_test_case(on_line_2d_test, "on_line_2d_test")
        call run_test_case(apply_pbc_test, "apply_pbc_test")
        call run_test_case(apply_pbc_tilted_test, "apply_pbc_tilted_test")

    end subroutine lattice_mod_test_driver 

    subroutine apply_pbc_test

        print *, ""
        print *, "testing: apply_pbc"
        print *, "TODO!"
!         call stop_all("apply_pbc_test", "todo")

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

    subroutine test_init_cluster_lattice_aim
        use OneEInts, only: gettmatel, tmat2d
        use SystemData, only: nbasis

        class(aim), pointer :: ptr 

        integer :: i 

        nbasis = 4

        allocate(tmat2d(4,4))
        tmat2d = 0.0
        tmat2d(1,3) = 1.0
        tmat2d(2,4) = 1.0
        tmat2d(3,1) = 1.0
        tmat2d(4,2) = 1.0

        

        print *, "" 
        print *, "initialize a 'cluster-lattice-aim' with 1 impurity and 1 bath site"
        ptr => aim('cluster', 1, 1)

        call assert_equals(1, ptr%get_ndim() )
        call assert_equals(1, ptr%get_length() ) 
        call assert_equals(2, ptr%get_nsites() ) 
        call assert_equals(1, ptr%get_nconnect_max() ) 
        call assert_true( .not. ptr%is_periodic() ) 
        call assert_equals(1, ptr%get_site_index(1) ) 
        call assert_equals(2, ptr%get_site_index(2) )
        call assert_equals([2], ptr%get_neighbors(1), 1)
        call assert_equals([1], ptr%get_neighbors(2), 1) 

        call assert_equals([3], ptr%get_spinorb_neighbors(1), 1)
        call assert_equals([4], ptr%get_spinorb_neighbors(2), 1) 
        call assert_equals([1], ptr%get_spinorb_neighbors(3), 1)
        call assert_equals([2], ptr%get_spinorb_neighbors(4), 1) 

        call assert_true( ptr%is_impurity_site(1) )
        call assert_true( .not. ptr%is_impurity_site(2) ) 
        call assert_true( ptr%is_bath_site(2) )
        call assert_true( .not. ptr%is_bath_site(1) )

        call assert_equals([1], ptr%get_impurities(), 1)
        call assert_equals([2], ptr%get_bath(), 1)

        call assert_equals(1, ptr%get_num_neighbors(1))
        call assert_equals(1, ptr%get_num_neighbors(2))

        call aim_deconstructor(ptr)

        call assert_true(.not.associated(ptr)) 

        deallocate(tmat2d)
        nbasis = -1

        nbasis = 202
        allocate(tmat2d(nbasis,nbasis))
        tmat2d = 0.0

        tmat2d(1,:) = 1.0
        tmat2d(2,:) = 1.0
        tmat2d(:,1) = 1.0
        tmat2d(:,2) = 1.0 
        print *, "" 
        print *, "initialize 1 impurity, 100-bath site 'aim-star' geometry"
        ptr => aim('cluster', 1, 100)

        call assert_equals(1, ptr%get_ndim() )
        call assert_equals(1, ptr%get_length() ) 
        call assert_equals(101, ptr%get_nsites() ) 
        call assert_equals(100, ptr%get_nconnect_max() ) 
        call assert_true( .not. ptr%is_periodic() ) 
        call assert_equals(1, ptr%get_site_index(1) ) 
        call assert_equals(2, ptr%get_site_index(2) )
        call assert_equals(101, ptr%get_site_index(101) )
        call assert_equals( [ (i, i = 2, 101) ], ptr%get_neighbors(1), 100)
        call assert_equals([1], ptr%get_neighbors(2), 1) 
        call assert_equals([1], ptr%get_neighbors(3), 1) 
        call assert_equals([1], ptr%get_neighbors(101), 1) 

        call assert_equals( [ (i, i = 3, 201, 2) ], ptr%get_spinorb_neighbors(1), 100)
        call assert_equals( [ (i, i = 4, 202, 2) ], ptr%get_spinorb_neighbors(2), 100)
        call assert_equals([1], ptr%get_spinorb_neighbors(3), 1) 
        call assert_equals([2], ptr%get_spinorb_neighbors(4), 1) 
        call assert_equals([1], ptr%get_spinorb_neighbors(5), 1) 
        call assert_equals([2], ptr%get_spinorb_neighbors(6), 1) 
        call assert_equals([1], ptr%get_spinorb_neighbors(201), 1) 
        call assert_equals([2], ptr%get_spinorb_neighbors(202), 1) 

        call assert_equals(100, ptr%get_num_neighbors(1))
        call assert_equals(1, ptr%get_num_neighbors(2))
        call assert_equals(1, ptr%get_num_neighbors(3))
        call assert_equals(1, ptr%get_num_neighbors(101))

        call assert_true( ptr%is_impurity_site(1) )
        call assert_true( .not. ptr%is_impurity_site(2) ) 
        call assert_true( .not. ptr%is_impurity_site(101) ) 
        call assert_true( ptr%is_bath_site(2) )
        call assert_true( .not. ptr%is_bath_site(1) )

        call assert_true( ptr%is_bath_site(3) )
        call assert_true( ptr%is_bath_site(101) )

        call assert_equals([1], ptr%get_impurities(), 1)
        call assert_equals([(i,i=2,101)], ptr%get_bath(), 1)


        call aim_deconstructor(ptr)

        call assert_true(.not.associated(ptr)) 

        deallocate(tmat2d) 
        nbasis = 6

        allocate(tmat2d(nbasis,nbasis))

        tmat2d = 0.0

        tmat2d(1,:) = 1.0
        tmat2d(2,:) = 1.0
        tmat2d(3,:) = 1.0
        tmat2d(4,:) = 1.0 
        tmat2d(5,:) = 1.0
        tmat2d(6,:) = 1.0

        print *, "" 
        print *, "initialize a 2 impurity 1 bath site geometry"
        ptr => aim('cluster', 2, 1) 
        call assert_equals(2, ptr%get_ndim() )
        call assert_equals(1, ptr%get_length() )
        call assert_equals(3, ptr%get_nsites() )
        call assert_equals(2, ptr%get_nconnect_max() )
        call assert_true( .not. ptr%is_periodic() )
        call assert_equals(1, ptr%get_site_index(1) ) 
        call assert_equals(2, ptr%get_site_index(2) )
        call assert_equals([2,3], ptr%get_neighbors(1), 2)
        call assert_equals([1,3], ptr%get_neighbors(2), 2)
        call assert_equals([1,2], ptr%get_neighbors(3), 2)

        call assert_equals(2, ptr%get_num_neighbors(1))
        call assert_equals(2, ptr%get_num_neighbors(2))
        call assert_equals(2, ptr%get_num_neighbors(3))

        call assert_true( ptr%is_impurity_site(1) )
        call assert_true( ptr%is_impurity_site(2) )
        call assert_true( .not. ptr%is_impurity_site(3) )
        call assert_true( ptr%is_bath_site(3) )
        call assert_true( .not. ptr%is_bath_site(2) )
        call assert_true( .not. ptr%is_bath_site(1) )

        call assert_equals([1,2], ptr%get_impurities() , 2)
        call assert_equals([3], ptr%get_bath(), 1)

        call aim_deconstructor(ptr)

        call assert_true(.not.associated(ptr)) 
        
        deallocate(tmat2d)
        nbasis = 208
        allocate(tmat2d(nbasis,nbasis))
        tmat2d = 0.0 
        tmat2d(1:8,:) = 1.0
        tmat2d(:,1:8) = 1.0

        print *, "" 
        print *, "initialize a 4 impurity 100 bath site geometry" 
        ptr => aim('cluster', 4, 100)
        call assert_equals(2, ptr%get_ndim() )
        call assert_equals(1, ptr%get_length() )
        call assert_equals(104, ptr%get_nsites() )
        call assert_equals(103, ptr%get_nconnect_max() )
        call assert_true(.not. ptr%is_periodic() )
        call assert_equals(1, ptr%get_site_index(1) )
        call assert_equals(2, ptr%get_site_index(2) )
        call assert_equals(104, ptr%get_site_index(104) )
        call assert_equals([(i, i = 2, 104)], ptr%get_neighbors(1), 103)
        call assert_equals([1, (i, i = 3, 104)], ptr%get_neighbors(2), 103)
        call assert_equals([1,2,3,4], ptr%get_neighbors(5), 4)
        call assert_equals([1,2,3,4], ptr%get_neighbors(104), 4)

        call assert_equals(103, ptr%get_num_neighbors(1))
        call assert_equals(103, ptr%get_num_neighbors(3))
        call assert_equals(4, ptr%get_num_neighbors(5))
        call assert_equals(4, ptr%get_num_neighbors(104))

        call assert_true( ptr%is_impurity_site(1) )
        call assert_true( ptr%is_impurity_site(4) )
        call assert_true( .not. ptr%is_impurity_site(5) )
        call assert_true( .not. ptr%is_impurity_site(104) )
        call assert_true( ptr%is_bath_site(5) )
        call assert_true( ptr%is_bath_site(104) )
        call assert_true( .not. ptr%is_bath_site(1) )
        call assert_true( .not. ptr%is_bath_site(4) )
        call assert_equals([1,2,3,4], ptr%get_impurities() , 4)
        call assert_equals([(i, i = 5, 104)], ptr%get_bath(), 100)

        call aim_deconstructor(ptr)

        deallocate(tmat2d)
        nbasis = -1

    end subroutine test_init_cluster_lattice_aim

    subroutine test_init_lattice_cube
        class(lattice), pointer :: ptr 

        print *, ""
        print *, "initialize a 2x2x2 cubic lattice with pbc"
        ptr => lattice('cube', 2,2,2, .true.,.true.,.true.)
        call assert_equals(3, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(2, ptr%get_length(3))
        call assert_equals(8, ptr%get_nsites())
        call assert_equals(6, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))
        call assert_true( ptr%is_periodic(3))

        ! now check the connectivity 
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))
        
        call assert_equals([2,3,5], ptr%get_neighbors(1), 3)
        call assert_equals([1,4,6], ptr%get_neighbors(2), 3)
        call assert_equals([1,4,7], ptr%get_neighbors(3), 3)
        call assert_equals([2,3,8], ptr%get_neighbors(4), 3)
        call assert_equals([1,6,7], ptr%get_neighbors(5), 3)
        call assert_equals([2,5,8], ptr%get_neighbors(6), 3)
        call assert_equals([3,5,8], ptr%get_neighbors(7), 3)
        call assert_equals([4,6,7], ptr%get_neighbors(8), 3)

        call assert_equals([8,12,14], ptr%get_spinorb_neighbors(16), 3)

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 2x2x3 cubic lattice with pbc"
        ptr => lattice('cube', 2,2,3, .true.,.true.,.true.)
        call assert_equals(3, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(2, ptr%get_length(2))
        call assert_equals(3, ptr%get_length(3))
        call assert_equals(12, ptr%get_nsites())
        call assert_equals(6, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))
        call assert_true( ptr%is_periodic(3))

        ! now check the connectivity 
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))
        
        call assert_equals([2,3,5,9], ptr%get_neighbors(1), 4)
        call assert_equals([1,4,6,10], ptr%get_neighbors(2), 4)
        call assert_equals([1,4,7,11], ptr%get_neighbors(3), 4)
        call assert_equals([2,3,8,12], ptr%get_neighbors(4), 4)
        call assert_equals([1,6,7,9], ptr%get_neighbors(5), 4)
        call assert_equals([2,5,8,10], ptr%get_neighbors(6), 4)
        call assert_equals([3,5,8,11], ptr%get_neighbors(7), 4)
        call assert_equals([4,6,7,12], ptr%get_neighbors(8), 4)
        call assert_equals([1,5,10,11], ptr%get_neighbors(9), 4)
        call assert_equals([2,6,9,12], ptr%get_neighbors(10), 4)
        call assert_equals([3,7,9,12], ptr%get_neighbors(11), 4)
        call assert_equals([4,8,10,11], ptr%get_neighbors(12), 4)

        call assert_equals([7,15,19,21], ptr%get_spinorb_neighbors(23), 4)

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 3x3x3 cubic lattice with pbc"
        ptr => lattice('cube', 3,3,3, .true.,.true.,.true.)
        call assert_equals(3, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(3, ptr%get_length(3))
        call assert_equals(27,ptr%get_nsites())
        call assert_equals(6, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))
        call assert_true( ptr%is_periodic(3))

        ! now check the connectivity 
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))
        
        call assert_equals([2,3,4,7,10,19], ptr%get_neighbors(1), 6)
        call assert_equals([1,3,5,8,11,20], ptr%get_neighbors(2), 6)
        call assert_equals([1,2,6,9,12,21], ptr%get_neighbors(3), 6)
        call assert_equals([1,5,6,7,13,22], ptr%get_neighbors(4), 6)
        call assert_equals([2,4,6,8,14,23], ptr%get_neighbors(5), 6)
        call assert_equals([3,4,5,9,15,24], ptr%get_neighbors(6), 6)
        call assert_equals([1,4,8,9,16,25], ptr%get_neighbors(7), 6)
        call assert_equals([2,5,7,9,17,26], ptr%get_neighbors(8), 6)
        call assert_equals([3,6,7,8,18,27], ptr%get_neighbors(9), 6)
        call assert_equals([1,11,12,13,16,19], ptr%get_neighbors(10), 6)
        call assert_equals([2,10,12,14,17,20], ptr%get_neighbors(11), 6)
        call assert_equals([3,10,11,15,18,21], ptr%get_neighbors(12), 6)
        call assert_equals([5,14,20,22,24,26], ptr%get_neighbors(23), 6)

        call assert_equals([10,28,40,44,48,52], ptr%get_spinorb_neighbors(46), 6)

        call lattice_deconstructor(ptr)

    end subroutine test_init_lattice_cube

    subroutine init_lattice_ole_test
        class(lattice), pointer :: ptr 
        real(dp) :: x(24)
        integer :: i

        print *, "" 
        print *, "testing: initialize a 3x5 Ole lattice" 
        ptr => lattice('ole', 3, 5, 1, .true.,.true.,.true.,'k-space')
        call assert_equals(2, ptr%get_ndim())
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(5, ptr%get_length(2))
        call assert_equals(24, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity 
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))
        call assert_equals([3,9,18,22], ptr%get_neighbors(1),4)
        call assert_equals([3,6,22,24], ptr%get_neighbors(2),4)
        call assert_equals([1,2,4,7], ptr%get_neighbors(3),4)
        call assert_equals([3,8,9,14], ptr%get_neighbors(4),4)
        call assert_equals([6,9,18,24], ptr%get_neighbors(5),4)
        call assert_equals([1,5,17,23], ptr%get_neighbors(18),4)
        call assert_equals([2,5,21,23], ptr%get_neighbors(24),4)

        x = [(-ptr%dispersion_rel_orb(i), i = 1,24 )]

        call sort(x)
        do i = 1, 24 
            print *, "e(x): ", x(i)
        end do

        print *, "check neighbors: "
        do i = 1, 24 
            print *, "i, neighbors", i, "|", ptr%get_neighbors(i)
        end do

        print *, "testing a 2x4 Ole lattice: "
        ptr => lattice('ole',2,4,1,.true.,.true.,.true.,'k-space')

        do i = 1, 12 
            print *, "i | neighbors:", i, "|", ptr%get_neighbors(i)
        end do
        call assert_equals(2, ptr%get_ndim())
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(4, ptr%get_length(2))
        call assert_equals(12, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))

        ! now check the connectivity 
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(8, ptr%get_site_index(8))

        call assert_equals([3,5,10,12], ptr%get_neighbors(1),4)
        call assert_equals([3,5,10,12], ptr%get_neighbors(2),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(3),4)
        call assert_equals([3,5,7,8], ptr%get_neighbors(4),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(5),4)
        call assert_equals([1,2,9,11], ptr%get_neighbors(10),4)
        call assert_equals([1,2,9,11], ptr%get_neighbors(12),4)

        x(1:12) = [(-ptr%dispersion_rel_orb(i), i = 1, 12)]
        call sort(x(1:12))

        print *, "e(x): "
        do i = 1, 12 
            print *, x(i)
        end do

    end subroutine init_lattice_ole_test

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

        print *, "initialize a 2x3 12-site tilted lattice: "
        ptr => lattice('tilted', 2,3,1,.true.,.true.,.true.,'k-space')
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(12, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))
        call assert_equals([3,5,11,12], ptr%get_neighbors(1),2)
        call assert_equals([3,5,11,12], ptr%get_neighbors(2),2)
        call assert_equals([1,2,4,6], ptr%get_neighbors(3),2)
        call assert_equals([3,5,7,9], ptr%get_neighbors(4),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(5),4)
        call assert_equals([3,5,7,9], ptr%get_neighbors(6),4)
        call assert_equals([4,6,8,10], ptr%get_neighbors(7),4)
        call assert_equals([7,9,11,12], ptr%get_neighbors(8),4)
        call assert_equals([1,2,8,10], ptr%get_neighbors(12),4)

        print *, "initialize a 3x2 12-site tilted lattice: "
        print *, "due to symmetr the 2x3 is treated the dealt internally as" 
        print *, "the already working 2x3!"
        ptr => lattice('tilted', 3,2,1,.true.,.true.,.true.,'k-space')
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(12, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))
        call assert_equals([3,5,11,12], ptr%get_neighbors(1),4)
        call assert_equals([3,5,11,12], ptr%get_neighbors(2),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(3),4)
        call assert_equals([3,5,7,9], ptr%get_neighbors(4),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(5),4)
        call assert_equals([3,5,7,9], ptr%get_neighbors(6),4)
        call assert_equals([4,6,8,10], ptr%get_neighbors(7),4)
        call assert_equals([7,9,11,12], ptr%get_neighbors(8),4)
        call assert_equals([1,2,8,10], ptr%get_neighbors(12),4)

        print *, "initialize a 3x4 24-site tilted lattice: " 
        ptr => lattice('tilted', 3,4,1,.true.,.true.,.true.,'k-space')
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(4, ptr%get_length(2))
        call assert_equals(24, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))
        call assert_equals([3,10,20,24], ptr%get_neighbors(1),4)
        call assert_equals([3,6,20,23], ptr%get_neighbors(2),4)
        call assert_equals([1,2,4,7], ptr%get_neighbors(3),4)
        call assert_equals([3,8,10,16], ptr%get_neighbors(4),4)
        call assert_equals([6,10,23,24], ptr%get_neighbors(5),4)
        call assert_equals([2,5,7,11], ptr%get_neighbors(6),4)
        call assert_equals([3,6,8,12], ptr%get_neighbors(7),4)
        call assert_equals([13,17,19,22], ptr%get_neighbors(18),4)
        call assert_equals([1,5,15,22], ptr%get_neighbors(24),4)

        do i = 1, 24 
            print *, "k, e(k): ", ptr%get_k_vec(i), ptr%dispersion_rel_orb(i)
        end do

        x = [(-ptr%dispersion_rel_orb(i), i = 1,24)]
        call sort(x)
!         x = x(24:1:-1)

        do i = 1,24 
            print *, x(i)
        end do

        print *, "initialize a 2x4 16-site tilted lattice: " 
        ptr => lattice('tilted', 2,4,1,.true.,.true.,.true.,'k-space')
        call assert_equals(2, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(2, ptr%get_length(1))
        call assert_equals(4, ptr%get_length(2))
        call assert_equals(16, ptr%get_nsites())
        call assert_equals(4, ptr%get_nconnect_max())
        call assert_true( ptr%is_periodic())
        call assert_true( ptr%is_periodic(1))
        call assert_true( ptr%is_periodic(2))
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(3, ptr%get_site_index(3))
        call assert_equals(4, ptr%get_site_index(4))
        call assert_equals([3,5,15,16], ptr%get_neighbors(1),4)
        call assert_equals([3,5,15,16], ptr%get_neighbors(2),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(3),4)
        call assert_equals([3,5,7,9], ptr%get_neighbors(4),4)
        call assert_equals([1,2,4,6], ptr%get_neighbors(5),4)
        call assert_equals([3,5,7,9], ptr%get_neighbors(6),4)
        call assert_equals([1,2,13,16], ptr%get_neighbors(16),4)
        call assert_equals(0.0_dp, ptr%dispersion_rel_orb(1))

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


        call lattice_deconstructor(ptr) 






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

    subroutine test_init_lattice_star
        class(lattice), pointer :: ptr

        integer :: i

        print *, ""
        print *, "initialize a 1 site 'star' geometry lattice: "
        ptr => lattice('star', 1, 1, 1, .false., .false., .false.)

        call assert_equals(1, ptr%get_ndim())
        call assert_equals(1, ptr%get_length())
        call assert_equals(1, ptr%get_nsites())
        call assert_equals(0, ptr%get_nconnect_max())
        call assert_true(.not. ptr%is_periodic())

        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(ptr%get_neighbors(1), [-1], 1)

        call lattice_deconstructor(ptr) 

        call assert_true(.not.associated(ptr)) 

        print *, "" 
        print *, "initialize a 2-site 'star' geometry lattice: "
        ptr => lattice('star', 2, 1, 1, .false., .false., .false.)
        call assert_equals(1, ptr%get_ndim())
        call assert_equals(1, ptr%get_length())
        call assert_equals(2, ptr%get_nsites())
        call assert_equals(1, ptr%get_nconnect_max())
        call assert_true(.not. ptr%is_periodic())

        call assert_equals(ptr%get_site_index(1), 1) 
        call assert_equals(ptr%get_site_index(2), 2) 
        call assert_equals(ptr%get_neighbors(1), [2], size(ptr%get_neighbors(1)))
        call assert_equals(ptr%get_neighbors(2), [1], size(ptr%get_neighbors(2)))

        call assert_equals(1, ptr%get_num_neighbors(1))
        call assert_equals(1, ptr%get_num_neighbors(2))

        call lattice_deconstructor(ptr) 

        call assert_true(.not.associated(ptr)) 

        print *, "" 
        print *, "initialize a 100 site 'star' geometry lattice: " 
        ptr => lattice('star', 100, 1, 1, .false., .false., .false.)
        call assert_equals(1, ptr%get_ndim())
        call assert_equals(1, ptr%get_length())
        call assert_equals(100, ptr%get_nsites())
        call assert_equals(99, ptr%get_nconnect_max())
        call assert_true(.not. ptr%is_periodic())

        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_equals(100, ptr%get_site_index(100))

        call assert_equals(1, ptr%get_num_neighbors(2))
        call assert_equals(1, ptr%get_num_neighbors(100))
        call assert_equals(99, ptr%get_num_neighbors(1))

        call assert_equals([1], ptr%get_neighbors(2), 1)
        call assert_equals([1], ptr%get_neighbors(100), 1)
        call assert_equals( [(i, i = 2, 100)], ptr%get_neighbors(1),99)

        call lattice_deconstructor(ptr) 

        call assert_true(.not.associated(ptr)) 

    end subroutine test_init_lattice_star
 
    subroutine test_init_lattice_aim_star()
        class(aim), pointer :: ptr 
        
        integer :: i 

        print *, "" 
        print *, "initialize 1 site 1 bath 'aim-star' geometry" 
        ptr => aim('aim-star', 1, 1) 

        call assert_equals(1, ptr%get_ndim() )
        call assert_equals(1, ptr%get_length() ) 
        call assert_equals(2, ptr%get_nsites() ) 
        call assert_equals(1, ptr%get_nconnect_max() ) 
        call assert_true( .not. ptr%is_periodic() ) 
        call assert_equals(1, ptr%get_site_index(1) ) 
        call assert_equals(2, ptr%get_site_index(2) )
        call assert_equals([2], ptr%get_neighbors(1), 1)
        call assert_equals([1], ptr%get_neighbors(2), 1) 

        call assert_equals(1, ptr%get_num_neighbors(1))
        call assert_equals(1, ptr%get_num_neighbors(2))

        call assert_true( ptr%is_impurity_site(1) )
        call assert_true( .not. ptr%is_impurity_site(2) ) 
        call assert_true( ptr%is_bath_site(2) )
        call assert_true( .not. ptr%is_bath_site(1) )
        call assert_equals([1], ptr%get_impurities(), 1)
        call assert_equals([2], ptr%get_bath(), 1)

        call aim_deconstructor(ptr)

        call assert_true(.not.associated(ptr)) 

        print *, "" 
        print *, "initialize 2-bath site 'aim-star' geometry"
        ptr => aim('star', 1, 2)

        call assert_equals(1, ptr%get_ndim() )
        call assert_equals(1, ptr%get_length() ) 
        call assert_equals(3, ptr%get_nsites() ) 
        call assert_equals(2, ptr%get_nconnect_max() ) 
        call assert_true( .not. ptr%is_periodic() ) 
        call assert_equals(1, ptr%get_site_index(1) ) 
        call assert_equals(2, ptr%get_site_index(2) )
        call assert_equals([2,3], ptr%get_neighbors(1), 2)
        call assert_equals([1], ptr%get_neighbors(2), 1) 
        call assert_equals([1], ptr%get_neighbors(3), 1) 

        call assert_equals(2, ptr%get_num_neighbors(1))
        call assert_equals(1, ptr%get_num_neighbors(2))
        call assert_equals(1, ptr%get_num_neighbors(3))

        call assert_true( ptr%is_impurity_site(1) )
        call assert_true( .not. ptr%is_impurity_site(2) ) 
        call assert_true( .not. ptr%is_impurity_site(3) ) 
        call assert_true( ptr%is_bath_site(2) )
        call assert_true( .not. ptr%is_bath_site(1) )

        call assert_true( ptr%is_bath_site(3) )

        call assert_equals([1], ptr%get_impurities(), 1)
        call assert_equals([2,3], ptr%get_bath(), 2)

        call aim_deconstructor(ptr)

        call assert_true(.not.associated(ptr)) 

        print *, "" 
        print *, "initialize 100-bath site 'aim-star' geometry"
        ptr => aim('star', 1, 100)

        call assert_equals(1, ptr%get_ndim() )
        call assert_equals(1, ptr%get_length() ) 
        call assert_equals(101, ptr%get_nsites() ) 
        call assert_equals(100, ptr%get_nconnect_max() ) 
        call assert_true( .not. ptr%is_periodic() ) 
        call assert_equals(1, ptr%get_site_index(1) ) 
        call assert_equals(2, ptr%get_site_index(2) )
        call assert_equals(101, ptr%get_site_index(101) )
        call assert_equals( [ (i, i = 2, 101) ], ptr%get_neighbors(1), 100)
        call assert_equals([1], ptr%get_neighbors(2), 1) 
        call assert_equals([1], ptr%get_neighbors(3), 1) 
        call assert_equals([1], ptr%get_neighbors(101), 1) 

        call assert_equals(100, ptr%get_num_neighbors(1))
        call assert_equals(1, ptr%get_num_neighbors(2))
        call assert_equals(1, ptr%get_num_neighbors(100))

        call assert_true( ptr%is_impurity_site(1) )
        call assert_true( .not. ptr%is_impurity_site(2) ) 
        call assert_true( .not. ptr%is_impurity_site(101) ) 
        call assert_true( ptr%is_bath_site(2) )
        call assert_true( .not. ptr%is_bath_site(1) )

        call assert_true( ptr%is_bath_site(3) )
        call assert_true( ptr%is_bath_site(101) )

        call assert_equals([1], ptr%get_impurities(), 1)
        call assert_equals([(i, i = 2, 101)], ptr%get_bath(), 100)

        call aim_deconstructor(ptr)

        call assert_true(.not.associated(ptr)) 

    end subroutine test_init_lattice_aim_star

    subroutine test_init_lattice_aim_chain
        class(aim), pointer :: ptr 

        integer :: i
        ! the question is: how do i deal with U and the hopping? 
        ! do i encode it in the lattice information? or do i just provide the 
        ! indices here
        print *, ""
        print *, "initialize a 'aim-chain' lattice with 1 bath site: "
        ptr => aim('aim-chain', 1, 1)

        call assert_equals(1, ptr%get_n_imps())
        call assert_equals(1, ptr%get_n_bath())
        call assert_true(.not. ptr%is_periodic())
        call assert_equals(1, ptr%get_ndim())
        call assert_equals(2, ptr%get_nsites())
        call assert_equals(2, ptr%get_length())
        call assert_equals(1, ptr%get_nconnect_max())
        call assert_equals(1, ptr%get_site_index(1))
        call assert_equals(2, ptr%get_site_index(2))
        call assert_true(ptr%is_impurity_site(1))
        call assert_true(ptr%is_bath_site(2))
        call assert_true(.not. ptr%is_bath_site(1))
        call assert_true(.not. ptr%is_impurity_site(2))
        call assert_equals([1], ptr%get_impurities(), 1)
        call assert_equals([2], ptr%get_bath(), 1)
        call assert_equals([2], ptr%get_neighbors(1), 1) 
        call assert_equals([1], ptr%get_neighbors(2), 1)
        
        call assert_equals(1, ptr%get_num_neighbors(1))
        call assert_equals(1, ptr%get_num_neighbors(2))

        call aim_deconstructor(ptr) 

        call assert_true(.not.associated(ptr)) 

        print *, "" 
        print *, "initialize a 'aim-chain' lattice with 100 bath sites: "
        ptr => aim('aim-chain', 1, 100) 

        call assert_equals(1, ptr%get_n_imps() )
        call assert_equals(100, ptr%get_n_bath() )
        call assert_true( .not. ptr%is_periodic() )
        call assert_equals(1 , ptr%get_ndim() )
        call assert_equals(101, ptr%get_nsites() )
        call assert_equals(101, ptr%get_length() )
        call assert_equals(2, ptr%get_nconnect_max() ) 
        call assert_equals(1, ptr%get_site_index(1) )
        call assert_equals(2, ptr%get_site_index(2) )
        call assert_equals(101, ptr%get_site_index(101) )
        call assert_true( ptr%is_impurity_site(1) )
        call assert_true( .not. ptr%is_bath_site(1) )
        call assert_true( ptr%is_bath_site(2) ) 
        call assert_true( ptr%is_bath_site(101) )
        call assert_true( .not. ptr%is_impurity_site(2) ) 
        call assert_true( .not. ptr%is_impurity_site(101) )
        call assert_equals([1], ptr%get_impurities(), 1) 
        call assert_equals( [ (i, i = 2,101) ], ptr%get_bath(), 100)
        call assert_equals([2], ptr%get_neighbors(1), 1) 
        call assert_equals([1,3], ptr%get_neighbors(2), 2) 
        call assert_equals([100], ptr%get_neighbors(101),1)

        call assert_equals(1, ptr%get_num_neighbors(1))
        call assert_equals(2, ptr%get_num_neighbors(2))
        call assert_equals(2, ptr%get_num_neighbors(50))
        call assert_equals(1, ptr%get_num_neighbors(101))

        call aim_deconstructor(ptr) 

        call assert_true(.not.associated(ptr)) 

    end subroutine test_init_lattice_aim_chain

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
