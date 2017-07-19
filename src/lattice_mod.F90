#include "macros.h"

module lattice_mod 
    ! this will be the module to mimick oop in fortran and in the 
    ! lattice excitation generation implementation in neci. 
    ! my plan is to create classes of lattice implementations which 
    ! all come from a same base class "lattice" 
    ! and maybe also do the same for sites.. where in the AIM then eg. 
    ! the bath and impurity sites extend the base site class 

    implicit none 
    private 
    public :: chain, lattice, constructor

    type :: site 
        ! the basic site type for my lattice
        ! i guess here i need to store my neighbors and provide functionality 
        ! how to get them 
        ! and i think i want to store this data contigous in memory to 
        ! speed up cache, access
        private 
        integer :: n_neighbors = -1 
        ! ah.. here it is important: neighbors are also sites.. is this 
        ! possible? and i have to be flexible to allow different types of 
        ! site neighbors 
        class(site), pointer :: neighbors(:) 

    contains 
        private 
        procedure, public :: get_neighbors

    end type site 

    ! maybe make this abstract.. what are the benefits? 
    type :: lattice 
        private 
        ! base class of lattice 
        ! i think i want to try to store all this contigous in memory to 
        ! speed up chache access 
        integer :: n_sites = -1
        integer :: n_connect_max = -1 
        ! and in the end a lattice is a collection of sites 
        ! and all the topology could be stored in the connection of the 
        ! sites 
        type(site), allocatable :: sites(:)

    contains 
        private

        procedure, public :: initialize 
        procedure :: assign 
        generic :: assignment(=) => assign

    end type lattice 

    ! and the plan is than to extend this with the specific lattices 

    type, extends(lattice) :: chain
        private 
        integer :: n_dim = -1

    contains 
        private

    end type chain

    interface lattice
        procedure constructor 
    end interface

    interface chain 
        procedure test 
    end interface

    

contains 

    subroutine assign(to, from) 
        class(lattice), intent(out) :: to 
        class(lattice), intent(in) :: from 

!         allocate(to)
        to%n_sites = from%n_sites

    end subroutine assign


    subroutine initialize(this, n_sites) 
        ! and write the first dummy initialize 
        class(lattice) :: this 
        integer, intent(in) :: n_sites 
        character(*), parameter :: this_routine = "initialize"

        this%n_sites = n_sites 

        select type (this) 

        ! well i cannot init type is (lattice) if i choose to make it 
        ! abstract. since it is not allowed to ever be intantiated..

        class is (chain)
            this%n_dim = 1 

        class default 
            call stop_all(this_routine, "unexpected lattice type!")

        end select

    end subroutine initialize

    function get_neighbors(this) result(neighbors) 
        ! this is a generic routine to get the neighbors of a site 
        class(site) :: this 

        class(site), pointer :: neighbors(:)

    end function get_neighbors

    ! write a general public constructor for lattices 
    type(lattice) function constructor(n_sites) 
        integer, intent(in) :: n_sites
        call constructor%initialize(n_sites)
    end function constructor

    function test(n_sites) result(this)
        ! yes thats the way how the contructor should work in general! 
        ! but i guess it would be better to write a seperate module for each 
        ! of the "classes" so it is also better to make them publicly 
        ! available as the interfaced constructors.. 
        integer :: n_sites
        class(lattice), pointer :: this
        
        call this%initialize(n_sites)

    end function test

end module lattice_mod
