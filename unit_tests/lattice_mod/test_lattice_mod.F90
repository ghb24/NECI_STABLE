#include "macros.h" 

! test my lattice module here, if the initializations and getter and 
! setter functions work as intended 

program test_lattice_mod 

    use lattice_mod 
    use fruit 

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

    end subroutine lattice_mod_test_driver 

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

        call lattice_deconstructor(ptr)

        print *, ""
        print *, "initialize a 3x3x3 cubic lattice with pbc"
        ptr => lattice('cube', 3,3,3, .true.,.true.,.true.)
        call assert_equals(3, ptr%get_ndim())
        ! this would be nice: how do i implement that??
        call assert_equals(3, ptr%get_length(1))
        call assert_equals(3, ptr%get_length(2))
        call assert_equals(3, ptr%get_length(3))
        call assert_equals(27, ptr%get_nsites())
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

        call lattice_deconstructor(ptr)

    end subroutine test_init_lattice_cube

    subroutine test_init_lattice_tilted
        class(lattice), pointer :: ptr

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
        ! i actually do not need to have the common lattice type or? 
        ! when i want to test a chain i could just use a chain or? 

        call lattice_deconstructor(ptr) 

        call assert_true(.not.associated(ptr)) 

    end subroutine test_init_lattice_chain

end program test_lattice_mod
