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

    end subroutine lattice_mod_test_driver 

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

        deallocate(tmat2d)
        nbasis = -1

    end subroutine test_init_cluster_lattice_aim

    subroutine test_init_lattice_star
        class(lattice), pointer :: ptr

        integer :: i

        print *, ""
        print *, "initialize a 1 site 'star' geometry lattice: "
        ptr => lattice('star', 1, 1, .false., .false.)

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
        ptr => lattice('star', 2, 1, .false., .false.)
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
        ptr => lattice('star', 100, 1, .false., .false.) 
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
        ptr => lattice('chain', 1, 1, .true., .true.) 

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
        ptr => lattice('chain', 2, 1, .false., .false.) 
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
        ptr => lattice('chain', 0, 100, .true., .true.)
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
