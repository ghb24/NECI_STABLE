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
    public :: lattice, deconstructor

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
        procedure :: nullify_neighbors

    end type site 

    ! maybe make this abstract.. what are the benefits? 
    type :: lattice 
        private 
        ! base class of lattice 
        ! i think i want to try to store all this contigous in memory to 
        ! speed up chache access 
        integer :: n_sites = -1
        integer :: n_connect_max = -1 
        integer :: n_dim = -1
        logical :: t_periodic = .true. 
        ! and i think additionally i want to store which type of lattice 
        ! this is in a string or? so i do not always have to 
        ! use the select type functionality 
        ! i would need constant expression.. so stick with select type!

        ! and in the end a lattice is a collection of sites 
        ! and all the topology could be stored in the connection of the 
        ! sites 
        ! and i just realized that if i want to use class(lattice) 
        ! generally in the whole program, i have to provide all the 
        ! functionality already for the lattice.. atleast in a dummy 
        ! way.. hm.. maybe there is a better way to do it.. 
        class(site), allocatable :: sites(:)

    contains 
        private

        ! i think i need some general interface here at the top of the 
        ! type definition, so all of those function can get called 
        ! but i need specifice ones then for each sub-class 
        ! how do i do that? 
        procedure, public :: initialize 
        procedure, public :: get_nsites
        procedure, public :: get_ndim
        procedure, public :: get_nconnect
        procedure, public :: is_periodic
        procedure, public :: get_length => get_length_lattice

        ! i definetly also want to have a print function! 
        procedure, public :: print 

        ! maybe i want set routines too? 
        ! but i guess i want the private, because there is no need of 
        ! them being used outside of this module 
        ! but this would make it flexible to make these function public
        procedure :: set_nsites 
        procedure :: set_ndim
        procedure :: set_nconnect_max
        procedure :: set_periodic 
        procedure :: set_length => set_length_lattice
        procedure :: calc_nsites => calc_nsites_lattice
        procedure :: allocate_sites
        procedure :: initialize_sites

        procedure :: deallocate_sites

    end type lattice 

    ! and the plan is than to extend this with the specific lattices 

    type, extends(lattice) :: chain
        private 
    
        ! now.. the chain has a concept of length how will i generalize 
        ! this.. and can i use this 
        integer :: length = -1

    contains 
        private

        procedure, public :: get_length => get_length_chain

        procedure :: set_length => set_length_chain
        procedure :: calc_nsites => calc_nsites_chain

    end type chain

    interface lattice
        procedure constructor 
    end interface

    integer, parameter :: DIM_CHAIN = 1, N_CONNECT_MAX_CHAIN = 2
contains 

    subroutine initialize(this, length_x, length_y, t_periodic) 
        ! and write the first dummy initialize 
        class(lattice) :: this 
        integer, intent(in) :: length_x, length_y
        logical, intent(in) :: t_periodic
        character(*), parameter :: this_routine = "initialize"

        integer :: n_sites

        select type (this) 

        ! well i cannot init type is (lattice) if i choose to make it 
        ! abstract. since it is not allowed to ever be intantiated..

        class is (chain)
            ! set some defaults for the chain lattice type 
            call this%set_ndim(DIM_CHAIN)
            call this%set_nconnect_max(N_CONNECT_MAX_CHAIN) 

            ! the type specific routine deal with the check of the 
            ! lenghts! 

            ! i should not call this set_nsites since this really should 
            ! just imply that it set the variable 
            ! introduce a second routine, which first determines the 
            ! number of sites depending on the lattice type 

        class default 
            call stop_all(this_routine, "unexpected lattice type!")

        end select

        ! and for the rest i can call general routines:
        n_sites = this%calc_nsites(length_x, length_y)

        call this%set_nsites(n_sites) 
        call this%set_periodic(t_periodic)
        call this%set_length(length_x, length_y)

        call this%allocate_sites(n_sites) 

    end subroutine initialize

    function constructor(lat_typ, length_x, length_y, t_periodic) result(this)
        ! write a general public constructor for lattices 
        ! the number of inputs are still undecided.. do we always have 
        ! the same number or differing number of inputs? 
        ! i guess, since we will be using global variables which are either 
        ! read in or set as default to init it.. 
        character(*), intent(in) :: lat_typ
        integer, intent(in) :: length_x, length_y
        logical, intent(in) :: t_periodic
        class(lattice), pointer :: this 

        ! depending on the string input defining lattice type 
        ! initialize corresponding lattice 
        select case (lat_typ) 
        case ('chain') 

            allocate(chain::this) 

        case default 
            ! stop here because a incorrect lattice type was given 
            call stop_all('lattice%constructor', & 
                'incorrect lattice type provided in lattice constructor!')

        end select 

        ! the initializer deals with the different types then.. 
        call this%initialize(length_x, length_y, t_periodic)

    end function constructor

    subroutine allocate_sites(this, n_sites) 
        ! do we want to 
        class(lattice) :: this 
        integer, intent(in) :: n_sites 
        character(*), parameter :: this_routine = "allocate_sites" 

        if (allocated(this%sites)) then 
            call stop_all(this_routine, "sites are already allocated!")
        end if
            
        if (n_sites < 1) then 
            call stop_all(this_routine, "0 or negative number of sites!")

        else 
            ! for now it is fine to allocate the sites just as "normal" 
            ! sites and not as bath/impurity sites.. 
            ! BUT: keep this in mind! 

            allocate(this%sites(n_sites)) 

        end if

    end subroutine allocate_sites


    subroutine deallocate_sites(this) 
        class(lattice) :: this
        integer :: i
        ! do i need an extra deallocater for the sites? 
        if (allocated(this%sites)) then 
            ! i have to run over all the sites and deallocate/nullify the 
            ! neighbor pointers! 
            do i = 1, this%get_nsites() 
                this%sites(1)%nullify_neighbors() 
            end do

            deallocate(this%sites)
        end if

    end subroutine deallocate_sites

    subroutine deconstructor(this) 
        ! routine to nullify the pointer to a lattice class 
        class(lattice), pointer :: this 

        ! first be sure that no sites are allocated 
        call this%deallocate_sites() 

        nullify(this) 

    end subroutine deconstructor
    

    function get_neighbors(this) result(neighbors) 
        ! this is a generic routine to get the neighbors of a site 
        class(site) :: this 

        class(site), pointer :: neighbors(:)

    end function get_neighbors

    function calc_nsites_chain(this, length_x, length_y) result(n_sites) 
        ! i acually do not want to rely on the previous calculated 
        ! length object in the class, since it is easy to recalc from 
        ! here and so i remove some dependencies.. 
        ! nah.. i can reuse this here in the set_length routine! 
        ! since the length equal the number of sites in the chain! 
        class(chain) :: this 
        integer, intent(in) :: length_x, length_y 
        integer :: n_sites 
        character(*), parameter :: this_routine = "calc_nsites_chain" 


        if (max(length_x,length_y) < 1 .or. min(length_x, length_y) > 1 .or. & 
            min(length_x,length_y) < 0) then 
            n_sites = -1 
            print *, "length_x: ", length_x
            print *, "length_y: ", length_y
            call stop_all(this_routine, "something went wrong in lenght input!")

        else 
            n_sites = max(length_x, length_y)

        end if

    end function calc_nsites_chain

    function calc_nsites_lattice(this, length_x, length_y) result(n_sites) 
        class(lattice) :: this 
        integer, intent(in) :: length_x, length_y
        integer :: n_sites 
        character(*), parameter :: this_routine = "calc_nsites_lattice"

        call stop_all(this_routine, &
            'type(lattice) should never be actually instantiated!') 

    end function calc_nsites_lattice

    ! Setter and getter routines for the private data of the lattice types
    function get_length_lattice(this) result(length)
        class(lattice) :: this 
        integer :: length
        character(*), parameter :: this_routine = "get_length_lattice"

        call stop_all(this_routine, &
            'type(lattice) should never be actually instantiated!') 

    end function get_length_lattice

    subroutine set_length_lattice(this, length_x, length_y) 
        class(lattice) :: this 
        integer, intent(in) :: length_x, length_y
        character(*), parameter :: this_routine = "set_length_lattice"

        call stop_all(this_routine, &
            'type(lattice) should never be actually instantiated!') 

    end subroutine set_length_lattice

    subroutine set_length_chain(this, length_x, length_y) 
        class(chain) :: this 
        integer, intent(in) :: length_x, length_y 
        character(*), parameter :: this_routine = "set_length_chain"

        ! the input checkin is all done in the calc_nsites routine!
        this%length = this%calc_nsites(length_x, length_y)

    end subroutine set_length_chain


    subroutine set_nsites(this, n_sites) 
        class(lattice) :: this 
        integer, intent(in) :: n_sites 

        this%n_sites = n_sites

    end subroutine set_nsites

    subroutine set_ndim(this, n_dim) 
        class(lattice) :: this 
        integer, intent(in) :: n_dim 

        this%n_dim = n_dim 

    end subroutine set_ndim

    subroutine set_nconnect_max(this, n_connect_max) 
        class(lattice) :: this 
        integer, intent(in) :: n_connect_max

        this%n_connect_max = n_connect_max

    end subroutine set_nconnect_max

    subroutine set_periodic(this, t_periodic) 
        class(lattice) :: this 
        logical, intent(in) :: t_periodic

        this%t_periodic = t_periodic

    end subroutine set_periodic

!     subroutine set_length(this, length)
!         class(chain) :: this 
!         integer, intent(in) :: length
! 
!         this%length = length 
! 
!     end subroutine set_length

    integer function get_nconnect(this) 
        class(lattice) :: this 
        
        get_nconnect = this%n_connect_max

    end function get_nconnect

    function get_length_chain(this) result(length)
        class(chain) :: this 
        integer :: length

        length = this%length

    end function get_length_chain

    logical function is_periodic(this) 
        class(lattice) :: this 

        is_periodic = this%t_periodic

    end function is_periodic

    integer function get_ndim(this) 
        class(lattice) :: this 

        get_ndim = this%n_dim

    end function get_ndim

    integer function get_nsites(this) 
        class(lattice) :: this 

        get_nsites = this%n_sites

    end function get_nsites

    subroutine print(this) 
        class(lattice) :: this 

        ! depending on the type print specific lattice information
        select type (this)
        class is (lattice) 

            call stop_all("lattice%print()", &
                "'lattice type should never be directly instantiated!")

        class is (chain) 
        
            print *, "Lattice type is: 'chain' "
            print *, " number of dimensions: ", this%get_ndim()
            print *, " max-number of neighbors: ", this%get_nconnect() 
            print *, " number of sites of chain: ", this%get_nsites() 
            print *, " is the chain periodic: ", this%is_periodic()

        end select 

    end subroutine print

    ! general non-type bound routines
!     subroutine get_lattice_type(this, string) 
!         class(lattice) :: this 
!         character(*), intent(out) :: string 
! 
!         select type (this) 
!         class is (chain) 
! 
!             string = 'chain'
! 
!         end select 
! 
!     end subroutine get_lattice_type
! 
!     function is_valid_lattice(this) result(flag) 
!         ! maybe it would be nice to have a routine which checks the basic 
!         ! necessities to be a proper lattice 
!         class(lattice) :: this 
!         logical :: flag
! 
!         character(*) :: string
! 
!         call this%get_lattice_type(string)
! 
!         select case (string)
!         case ('chain','lattice') 
!             flag = .true. 
! 
!         case default 
!             flag = .false. 
! 
!         end select 
! 
!     end function is_valid_lattice

end module lattice_mod
