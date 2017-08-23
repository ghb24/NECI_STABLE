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
    public :: lattice, lattice_deconstructor

    type :: site 
        ! the basic site type for my lattice
        ! i guess here i need to store my neighbors and provide functionality 
        ! how to get them 
        ! and i think i want to store this data contigous in memory to 
        ! speed up cache, access
        private 
        ! to i want to have an index, which gives me the number?
        ! i might not even need that? 
        ! in the mean-time use an index.. 
        integer :: ind = -1 
        integer :: n_neighbors = -1 
        ! ah.. here it is important: neighbors are also sites.. is this 
        ! possible? and i have to be flexible to allow different types of 
        ! site neighbors 
        ! and i am not even sure if i want pointers here.. maybe a vector 
        ! of neighbor indices is enough and even better.. 
        ! i cant really make it this way without another type since 
        ! fortran does not like arrays of pointers or atleast does not 
        ! intepret is as that from the beginning.. 
        ! but i think this would be perfect or? 
        ! recursive types are only allowed in the fortran 2008 standard.. 
        ! so i have to go over an intermediate type
        ! (which has the advantage to have an array of pointers.. 
        ! ok i realize this whole shabang is too much.. 
        ! so just make a list of indices of the neighbors
        integer, allocatable :: neighbors(:) 
!         class(neighbor), allocatable :: neighbors(:) 
        ! or can i point to the other sites here? hm.. 
        ! and maybe i want to store on-site repulsion here too.. lets see 

    contains 
        private 

        procedure :: allocate_neighbors
        procedure :: deallocate_neighbors

        procedure :: get_neighbors => get_neighbors_site
        ! maybe i do not need a initializer? 
        ! or maybe yes i do.. i would like to have something like a 
        ! with the lattice.. but now i want something optional.. 
        procedure :: initialize => init_site

        procedure :: set_index 
        procedure :: get_index
        procedure :: set_num_neighbors 
        procedure :: get_num_neighbors 
        procedure :: set_neighbors

        ! i could also use finalization routines instead of manually 
        ! deallocating everything.. 
        ! i need atleast gcc4.9.. which i am to lazy to update now..
        ! but will in the future!
!         final :: finalize_site

    end type site 

    ! maybe make this abstract.. what are the benefits? 
    type, abstract :: lattice 
        private 
        ! base class of lattice 
        ! i think i want to try to store all this contigous in memory to 
        ! speed up chache access 
        integer :: n_sites = -1
        integer :: n_connect_max = -1 
        integer :: n_dim = -1
        ! actually i want to have more flexibility: maybe periodic in x 
        ! but not y.. 
        logical :: t_periodic_x = .true. 
        logical :: t_periodic_y = .true. 
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
        ! maybe i have to use pointer attribute below to make it possible 
        ! to call an constructor of class(type) .. 
        ! well this also does not work as i like to have it.. 
        ! since it is not interpreted as an array of pointer, but as a 
        ! pointer to an array of class(sites), so redo this in the end!
        class(site), allocatable :: sites(:)

        ! this is just a small test if we can bring classic procedure
        ! pointers into the game.. but will be removed soon
        procedure(test), pointer :: a

        ! i just want to test if i can easily setup a lattice name
        ! no.. not yet supported in gcc4.8.. but in newer versions probably
!         character(*), allocatable :: lattice_name
    contains 
        private

        ! i think i need some general interface here at the top of the 
        ! type definition, so all of those function can get called 
        ! but i need specifice ones then for each sub-class 
        ! how do i do that? 
        procedure, public :: initialize => init_lattice
        procedure, public :: get_nsites
        procedure, public :: get_ndim
        procedure, public :: get_nconnect_max
        procedure, public :: is_periodic_x, is_periodic_y
        procedure(is_periodic_t), public, deferred :: is_periodic
        procedure(get_length_t), public, deferred :: get_length
!         procedure, public :: get_length => get_length_lattice
        procedure, public :: get_site_index
        ! make the get neighbors function public on the lattice level 
        procedure, public :: get_neighbors => get_neighbors_lattice

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
        procedure(set_length_t), deferred :: set_length
!         procedure :: set_length => set_length_lattice
        procedure(calc_nsites_t), deferred :: calc_nsites
!         procedure :: calc_nsites => calc_nsites_lattice
        procedure :: allocate_sites
        procedure(initialize_sites_t), deferred :: initialize_sites

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
        procedure, public :: is_periodic => is_periodic_chain

        procedure :: set_length => set_length_chain
        procedure :: calc_nsites => calc_nsites_chain
        procedure :: initialize_sites => init_sites_chain

    end type chain

    ! also implement a 'star' geometry, especially to deal with the AIM models
    type, extends(lattice) :: star
        private

    contains 
        private 

        procedure, public :: get_length => get_length_star
        ! i can just reuse the chain calc function..
        procedure, public :: calc_nsites => calc_nsites_star
        procedure, public :: is_periodic => is_periodic_star

        procedure :: initialize_sites => init_sites_star
        procedure :: set_length => set_length_star

    end type star

    ! i also want to have a rectangle.. 
    ! can i store every possibility in the rectangle type? 
    ! like also the tilted square lattice .. 
    ! the tilted essentially is only different boundary conditions.. 
    ! the rectangle type would be the basic 2D lattice with 4 nearest 
    ! neighbors.. 
    ! a square would be a special case of a rectangle with lx = ly 
    ! and a tilted would be a special case of a rectangle 
    ! and a tilted square? is it a special case of tilted or square? 
    ! but this is not of concern now.. actually i should finally 
    ! implement the rest of this code to actually work in neci! 
!     type, extends(lattice) :: rectangle 
!         private 
!         integer :: length_x, length_y
! 
!     contains 
! 
!     end type rectangle

    ! create the abstract interfaces for the deferred function in the base 
    ! abstract type: lattice
    abstract interface

        function get_length_t(this) result(length)
            import :: lattice 
            class(lattice) :: this 
            integer :: length
        end function get_length_t

        subroutine set_length_t(this, length_x, length_y) 
            import :: lattice 
            class(lattice) :: this 
            integer, intent(in) :: length_x, length_y
        end subroutine set_length_t

        function calc_nsites_t(this, length_x, length_y) result(n_sites) 
            import :: lattice
            class(lattice) :: this 
            integer, intent(in) :: length_x, length_y
            integer :: n_sites
        end function calc_nsites_t

        logical function is_periodic_t(this) 
            import :: lattice
            class(lattice) :: this 
        end function is_periodic_t

        subroutine initialize_sites_t(this) 
            import :: lattice 
            class(lattice) :: this 
        end subroutine initialize_sites_t

    end interface

    interface lattice
        procedure lattice_constructor 
!         procedure aim_lattice_constructor
    end interface

    ! can i make an abstract interface for the dummy procedures above?
    ! no not really.. since i would need class(lattice) in the interface 
    ! definition, but need the interface above already..or?
    ! maybe i do not have to declare the this argument here.. 
    ! since it is a proc ptr.. 
    ! ok.. this works now.. but i dont know for what yet.. 

    abstract interface 
        real function test(this) 
            import :: lattice
            class(lattice) :: this 
        end function test
    end interface

    interface site
        procedure site_constructor
    end interface

    interface assignment (=) 
        module procedure site_assign 
        module procedure lattice_assign
    end interface 
    integer, parameter :: DIM_CHAIN = 1, & 
                          DIM_STAR = 1, &
                          N_CONNECT_MAX_CHAIN = 2, &
                          STAR_LENGTH = 1

contains 

    subroutine set_num_neighbors(this, n_neighbors) 
        class(site) :: this 
        integer, intent(in) :: n_neighbors

        this%n_neighbors = n_neighbors 

    end subroutine set_num_neighbors

    pure integer function get_num_neighbors(this) 
        class(site), intent(in) :: this 

        get_num_neighbors = this%n_neighbors

    end function get_num_neighbors

    integer function get_site_index(this, ind)
        ! for now.. since i have not checked how efficient this whole 
        ! data-structure is to it in the nicest, safest way.. until i profile!
        class(lattice) :: this 
        integer, intent(in) :: ind
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_site_index"
#endif

        ASSERT(ind <= this%get_nsites()) 
        ASSERT(allocated(this%sites))

        get_site_index = this%sites(ind)%get_index()

    end function get_site_index

    subroutine init_sites_lattice(this) 
        ! these are the important routines.. which set up the connectivity 
        ! of the lattice sites. or atleast call the specific init-routines! 
        class(lattice) :: this

    end subroutine init_sites_lattice

    subroutine site_assign(lhs, rhs) 
        class(site), intent(out) :: lhs 
        class(site), intent(in), pointer :: rhs 

        ! argh.. on every change in the site constructor i have to make 
        ! the change here too.. annoying.. 
        ! but make i OOP style! 
        call lhs%set_index( rhs%get_index() )
        call lhs%set_num_neighbors( rhs%get_num_neighbors() )
        call lhs%allocate_neighbors( rhs%get_num_neighbors() ) 
        call lhs%set_neighbors( rhs%get_neighbors() )

        ! can i work with an allocatable statement here? 

    end subroutine site_assign

    subroutine lattice_assign(lhs, rhs) 
        ! specific lattice assigner to not have to use pointers.. 
        class(lattice), intent(out) :: lhs 
        class(lattice), intent(in), pointer :: rhs 

        character(*), parameter :: this_routine = "lattice_assign"

        ! here i have to copy all the specific values! 
        ! this is annoying but make the code more readable, and i do not 
        ! have to use pointers so much.. 
        ! todo although.. since i cannot use class(lattice) in the main 
        ! program without allocatable or pointer attribute anyway.. 
        ! i think there is no need of an assignment overload!
        call stop_all(this_routine, "not yet implemented!")
        
    end subroutine lattice_assign

    subroutine init_sites_star(this)
        ! create the lattice structure of the star geometry.. 
        ! with one special pivot site with index 1, which is connected 
        ! to all the other sites and the other sites are just connected to 
        ! the pivot site!
        class(star) :: this
        character(*), parameter :: this_routine = "init_sites_star"

        integer :: i 

        if (this%get_nsites() <= 0) then
            call stop_all(this_routine, & 
                "something went wrong: negative or 0 number of sites!")
        else if (this%get_nsites() == 1) then 
            this%sites(1) = site(ind = 1, n_neighbors = 0, neighbors = [integer :: ])

        else
            ! first to the special pivot site in the middle of the star
            this%sites(1) = site(ind = 1, n_neighbors = this%get_nconnect_max(), & 
                neighbors = [ (i, i = 2, this%get_nsites()) ])

            ! and all the others are just connected to the pivot
            do i = 2, this%get_nsites()
                this%sites(i) = site(ind = i, n_neighbors = 1, neighbors = [1])
            end do
        end if

    end subroutine init_sites_star

    subroutine init_sites_chain(this) 
        class(chain) :: this 

        integer :: i
        ! ok what exactly do i do here?? 
        ! and i think i want to write a constructor for the sites.. to do 
        ! nice initialization for AIM sites etc. 

        ! do i insist that other stuff is already set in the lattice type? 
        ! or do i take it as input here? 
        ! i think i inisist that it is set already! 
        ! and i have to deal with the edge case of a one-sited lattice 


        ! but how to i deal with the geometry here..
        ! i know it is a chain and i know how many sites and i know if it 
        ! is periodic or not.. 
        ! so just create the "lattice" 
        ! 1 - 2 - 3 - 4 - 5 - ... 
        ! maybe us an associate structure here to not always type so much.

        ! deal with the first and last sites specifically! 
        ! this initializes the site class with the index dummy variable:
        ! ok.. fortran does not really like arrays of pointers.. 
        ! i could make a workaround with yet another type or do i more 
        ! explicit here ..

        if (this%get_nsites() == 1) then 
            this%sites(1) = site(ind = 1, n_neighbors = 0, neighbors = [integer :: ])
            return 
        end if

        if (this%is_periodic()) then 

            ! use more concise site contructors!
            this%sites(1) = site(1, 2, [this%get_nsites(), 2])
            this%sites(this%get_nsites()) = site(this%get_nsites(), 2, & 
                [this%get_nsites() - 1, 1])

        else 
            ! open boundary conditions: 
            ! first site:
            this%sites(1) = site(1, 1, [2])

            ! last site: 
            this%sites(this%get_nsites()) = site(this%get_nsites(), 1, & 
                [this%get_nsites() - 1]) 

        end if 

        ! and do the rest inbetween which is always the same 
        do i = 2, this%get_nsites() - 1

            ! if periodic and first or last: already dealt with above
            this%sites(i) = site(i, N_CONNECT_MAX_CHAIN, [i - 1, i + 1])

        end do

    end subroutine init_sites_chain

    subroutine init_lattice(this, length_x, length_y, t_periodic_x, t_periodic_y) 
        ! and write the first dummy initialize 
        class(lattice) :: this 
        integer, intent(in) :: length_x, length_y
        logical, intent(in) :: t_periodic_x, t_periodic_y
        character(*), parameter :: this_routine = "init_lattice"

        integer :: n_sites


        n_sites = this%calc_nsites(length_x, length_y)

        ! and for the rest i can call general routines:
        call this%set_nsites(n_sites) 
        call this%set_periodic(t_periodic_x, t_periodic_y)

        select type (this) 

        ! well i cannot init type is (lattice) if i choose to make it 
        ! abstract. since it is not allowed to ever be intantiated..

        class is (chain)
            ! set some defaults for the chain lattice type 
            call this%set_ndim(DIM_CHAIN)

            call this%set_length(length_x, length_y)
            ! if incorrect lenght input it is caught in the calc_nsites above..
            if (this%get_length() == 1) then 
                call this%set_nconnect_max(0)
            else if (this%get_length() == 2 .and. (.not. this%is_periodic())) then 
                call this%set_nconnect_max(1)
            else 
                call this%set_nconnect_max(N_CONNECT_MAX_CHAIN) 
            end if

            ! the type specific routine deal with the check of the 
            ! lenghts! 

            ! i should not call this set_nsites since this really should 
            ! just imply that it set the variable 
            ! introduce a second routine, which first determines the 
            ! number of sites depending on the lattice type 

        class is (star) 
            call this%set_ndim(DIM_STAR)
            call this%set_nconnect_max(n_sites - 1)

            ! for the 'star' geometry the special point in the middle 
            ! is connected to all the others.. so i need to calc n_sites here.
            ! also check here if something went wrong in the input: 
            if (t_periodic_x .or. t_periodic_y) then 
                call stop_all(this_routine, &
                "incorrect initialization info: requested periodic 'star' geometry!")
            end if


        class default 
            call stop_all(this_routine, "unexpected lattice type!")

        end select
        ! do i want to allocate sites here or in the initializer? 
        ! well the specific site initializer will be different for all the 
        ! types of lattices.. so best would be to do everything which is 
        ! common to all routine here! 
        call this%allocate_sites(n_sites) 

        call this%initialize_sites() 

    end subroutine init_lattice

!     function aim_lattice_constructor() 
! 
!     end function aim_lattice_constructor

    function lattice_constructor(lat_typ, length_x, length_y, t_periodic_x , &
            t_periodic_y) result(this)
        ! write a general public lattice_constructor for lattices 
        ! the number of inputs are still undecided.. do we always have 
        ! the same number or differing number of inputs? 
        ! i guess, since we will be using global variables which are either 
        ! read in or set as default to init it.. 
        character(*), intent(in) :: lat_typ
        integer, intent(in) :: length_x, length_y
        logical, intent(in) :: t_periodic_x, t_periodic_y
        class(lattice), pointer :: this 

        ! depending on the string input defining lattice type 
        ! initialize corresponding lattice 
        select case (lat_typ) 
        case ('chain') 

            allocate(chain :: this) 

        case ('star') 

            allocate(star :: this) 

        case default 
            ! stop here because a incorrect lattice type was given 
            call stop_all('lattice_constructor', & 
                'incorrect lattice type provided in lattice_constructor!')

        end select 

        ! the initializer deals with the different types then.. 
        call this%initialize(length_x, length_y, t_periodic_x, t_periodic_y)

    end function lattice_constructor

    function site_constructor(ind, n_neighbors, neighbors, site_type) & 
            result(this) 
        ! similar to the lattice constructor i want to have a site 
        ! constructor, which depending on some input constructs the 
        ! specific sites on a lattice 
        ! for now the index is the only necessary input.. 
        ! include more in this site constructor here
        integer, intent(in) :: ind 
        integer, intent(in) :: n_neighbors 
        integer, intent(in) :: neighbors(n_neighbors)
        character(*), intent(in), optional :: site_type 
        ! i think i have to use pointers again.. 
        ! but maybe this is really bad to deal with in the rest of the code.. 
        class(site), pointer :: this 

        if (present(site_type)) then 
            ! not yet implementetd or to do.. so wait.. 

        else 
            ! this is the default case 
            allocate(site :: this) 

        end if

        call this%initialize(ind, n_neighbors, neighbors)

    end function site_constructor

    subroutine init_site(this, ind, n_neighbors, neighbors) 
        ! for now only use the index in the initalization 
        class(site) :: this 
        integer, intent(in) :: ind, n_neighbors
        integer, intent(in) :: neighbors(n_neighbors)

        ! for the beginning i do not need more than to set the index.. 
        ! independent of the type 
        call this%set_index(ind)
        call this%set_num_neighbors(n_neighbors)
        call this%allocate_neighbors(n_neighbors)
        call this%set_neighbors(neighbors)

    end subroutine init_site

    subroutine set_neighbors(this, neighbors) 
        class(site) :: this 
        integer, intent(in) :: neighbors(this%n_neighbors)

        this%neighbors = neighbors

    end subroutine set_neighbors

    subroutine allocate_neighbors(this, n_neighbors) 
        class(site) :: this 
        integer, intent(in) :: n_neighbors

        ! the procedure bound routine already checks if neighbors is 
        ! allocated.
        call this%deallocate_neighbors()

        allocate(this%neighbors(n_neighbors)) 

    end subroutine allocate_neighbors

    subroutine set_index(this, ind) 
        class(site) :: this 
        integer, intent(in) :: ind 

        this%ind = ind

    end subroutine set_index

    integer function get_index(this) 
        class(site) :: this 

        get_index = this%ind 

    end function get_index

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
                call this%sites(i)%deallocate_neighbors() 
            end do

            deallocate(this%sites)
        end if

    end subroutine deallocate_sites

    subroutine deallocate_neighbors(this) 
        class(site) :: this 

        if (allocated(this%neighbors)) deallocate(this%neighbors) 

    end subroutine deallocate_neighbors

    subroutine finalize_site(this) 
        type(site) :: this 
        if (allocated(this%neighbors)) deallocate(this%neighbors)
    end subroutine finalize_site

    subroutine lattice_deconstructor(this) 
        ! routine to nullify the pointer to a lattice class 
        class(lattice), pointer :: this 

        ! first be sure that no sites are allocated 
        call this%deallocate_sites() 

        nullify(this) 

    end subroutine lattice_deconstructor
    
    function get_neighbors_site(this) result(neighbors) 
        ! this is a generic routine to get the neighbors of a site 
        class(site) :: this 
        ! i need assumed array shape and size here or? check how i do that!
        ! can i use the stored number of neighbors in the type?
        ! i can use n_neighbors directly but not the function which gets me 
        ! n_neighbors.. strange and unfortunate.. 
        integer :: neighbors(this%n_neighbors)

        neighbors = this%neighbors

    end function get_neighbors_site

    function get_neighbors_lattice(this, ind) result(neighbors)
        class(lattice) :: this
        integer, intent(in) :: ind
        ! i can't really use the stored information here, since it is 
        ! input dependent and this can be out of bound.. exactly what i 
        ! want to avoid.. so allocate it here!
!         integer :: neighbors(this%sites(ind)%n_neighbors)
        integer, allocatable :: neighbors(:)
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_neighbors_lattice"
#endif 

        ! make all assert on a seperate line, so we exactly know what is 
        ! going wrong.. 
        ASSERT(ind <= this%get_nsites())
        ASSERT(ind > 0)
        ASSERT(allocated(this%sites))
        ASSERT(allocated(this%sites(ind)%neighbors))

        ! why doesn't this work?
!         neighbors = this%sites(ind)%get_neighbors() 

        ! apparently i have to access the data directly.. 
        if (this%get_nsites() == 1) then 
            ! output -1 to indicate that there are no neighbors, since the 
            ! lattice is too small and catch errors early.. or i could 
            ! users not allow to initialize such a lattice.. woudl also 
            ! make sense!
            allocate(neighbors(1)) 
            neighbors = -1
        else 
            allocate(neighbors(this%sites(ind)%get_num_neighbors())) 
            neighbors = this%sites(ind)%neighbors
        end if

    end function get_neighbors_lattice
! 
    function calc_nsites_star(this, length_x, length_y) result(n_sites) 
        ! the maximum of the input is used as the n_sites parameter! 
        ! this is the same function as the one for "chain" below.. 
        ! but i cannot use it somehow.. 
        class(star) :: this 
        integer, intent(in) :: length_x, length_y
        integer :: n_sites 
        character(*), parameter :: this_routine = "calc_nsites_star"

        if (max(length_x,length_y) < 1 .or. min(length_x, length_y) > 1 .or. & 
            min(length_x,length_y) < 0) then 
            n_sites = -1 
            call stop_all(this_routine, "something went wrong in lenght input!")

        else 
            n_sites = max(length_x, length_y)

        end if


    end function calc_nsites_star

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

    subroutine set_length_star(this, length_x, length_y)
        class(star) :: this 
        integer, intent(in) :: length_x, length_y
        character(*), parameter :: this_routine = "set_length_star"

        ! actually the length of a start is not really defined.. 
        ! maybe i should rethink if i make this part of the 
        ! original lattice class then.. 
        call stop_all(this_routine, "length not defined for 'star' geometry!")
    end subroutine set_length_star

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

    subroutine set_periodic(this, t_periodic_x, t_periodic_y) 
        class(lattice) :: this 
        logical, intent(in) :: t_periodic_x, t_periodic_y

        ! how do we decide which periodic flag to take? 
        ! well we just set it like that! the user has to be specific! 
        ! well.. but we do not want to ask if this x or y is periodic in 
        ! the 1dim case or? 
        ! periodic is the default.. and we want to turn off periodic by 
        ! inputting open-bc.. and in the case of the cain + open-bc this 
        ! means immediately that it is NOT periodic.. so only if both inputs 
        ! are true (which is the default case) the chain is set to be 
        ! periodic! (or in the is_periodic procedure, we check if 
        ! both are true!
        this%t_periodic_x = t_periodic_x
        this%t_periodic_y = t_periodic_y

    end subroutine set_periodic

    integer function get_nconnect_max(this) 
        class(lattice) :: this 
        
        get_nconnect_max = this%n_connect_max

    end function get_nconnect_max

    ! the star does not really have a concept of length.. 
    ! so always ouput STAR_LENGTH
    integer function get_length_star(this)
        class(star) :: this 
        get_length_star = STAR_LENGTH
    end function get_length_star

    integer function get_length_chain(this) 
        class(chain) :: this 

        get_length_chain = this%length

    end function get_length_chain

    logical function is_periodic_star(this)
        ! this is always false.. the star geometry can't be periodic
        class(star) :: this 
        is_periodic_star = .false.
    end function is_periodic_star

    logical function is_periodic_chain(this)
        class(chain) :: this 

        ! the chain is only treated as periodic if both the flags are set 
        ! to be periodic!
        is_periodic_chain = (this%is_periodic_x() .and. this%is_periodic_y())

    end function is_periodic_chain

    logical function is_periodic_lattice(this)
        class(lattice) :: this 
        character(*), parameter :: this_routine = "is_periodic_lattice"
        call stop_all(this_routine, &
            "lattice should not be directly instantiated!")


    end function is_periodic_lattice

    logical function is_periodic_x(this) 
        class(lattice) :: this 

        is_periodic_x = this%t_periodic_x

    end function is_periodic_x

    logical function is_periodic_y(this) 
        class(lattice) :: this 

        is_periodic_y = this%t_periodic_y

    end function is_periodic_y


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
            print *, " max-number of neighbors: ", this%get_nconnect_max() 
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
