#include "macros.h"

module lattice_mod 
    ! this will be the module to mimick oop in fortran and in the 
    ! lattice excitation generation implementation in neci. 
    ! my plan is to create classes of lattice implementations which 
    ! all come from a same base class "lattice" 
    ! and maybe also do the same for sites.. where in the AIM then eg. 
    ! the bath and impurity sites extend the base site class 

    use OneEInts, only: tmat2d
    use UmatCache, only: gtid
    use constants, only: dp, pi
    use SystemData, only: twisted_bc

    implicit none 
    private 
    public :: lattice, lattice_deconstructor, aim, aim_deconstructor, sort_unique, & 
              lat, determine_optimal_time_step, get_helement_lattice

    integer, parameter :: NAME_LEN = 13

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

        ! i could just hack into that i have a flag for impurities and 
        ! bath sites.. which i do not need for "normal" calculations but 
        ! for the AIM.. 
        logical :: t_impurity = .false. 
        logical :: t_bath = .false. 

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
        procedure :: get_num_neighbors => get_num_neighbors_site
        procedure :: set_neighbors

        procedure :: set_impurity 
        procedure :: is_impurity 
        procedure :: set_bath 
        procedure :: is_bath

        ! i could also use finalization routines instead of manually 
        ! deallocating everything.. 
        ! i need atleast gcc4.9.. which i am to lazy to update now..
        ! but will in the future!
!         final :: finalize_site

    end type site 

    ! can i do something like: 
    type, extends(site) :: impurity 
        private

    contains
        private

    end type impurity 

    type, extends(site) :: bath 
        private

    contains 
        private

    end type bath

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
        logical :: t_periodic(3) = .true. 

        ! i want to do a lattice type member also, do easier check, which 
        ! lattice we are looking at.. but i need to make this nice
        ! having a defered lenght character component would be the best 
        ! option: 
        !         character(:), allocatable, public :: type
        ! but this needs gfortran version >= 4.9, which i guess is not 
        ! available everywhere yet.. so do smth else for now: 
        character(NAME_LEN) :: name = ''
        ! and "just" use trim and align in the comparisons for lattice names! 

        ! and also add a flag, if we want a momentum space lattice 
        ! representation. 
        logical :: t_momentum_space = .false.

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
        type(site), allocatable :: sites(:)

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
        procedure :: initialize => init_lattice
        procedure, public :: get_nsites
        procedure, public :: get_ndim
        procedure, public :: get_nconnect_max
        procedure, public :: is_periodic_x
        procedure, public :: is_periodic_y
        procedure(is_periodic_t), public, deferred :: is_periodic
        procedure(get_length_t), public, deferred :: get_length
!         procedure, public :: get_length => get_length_lattice
        procedure, public :: get_site_index
        ! make the get neighbors function public on the lattice level 
        procedure, public :: get_neighbors => get_neighbors_lattice
        procedure, public :: get_num_neighbors => get_num_neighbors_lattice
        procedure, public :: get_spinorb_neighbors => get_spinorb_neighbors_lat

        procedure, public :: is_k_space 
        ! i definetly also want to have a print function! 
        procedure, public :: print 

        procedure :: set_name

        procedure, public :: get_name

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

        ! for the k-space implementations also implement a lattice 
        ! dependent dispersion relation function 
        procedure, public :: dispersion_rel => dispersion_rel_not_implemented

    end type lattice 
    ! and the plan is than to extend this with the specific lattices 

    ! i think it is better to extend lattice directly for the aim
    type, abstract, extends(lattice) :: aim 
        private 

        integer :: n_imps = -1 
        integer :: n_bath = -1 

        ! i think i want to store the bath and impurity indices.. 
        ! thats more efficient!
        integer, allocatable :: impurity_sites(:) 
        integer, allocatable :: bath_sites(:)

    contains 
        private 

!         procedure :: initialize => init_aim

        procedure :: set_n_imps
        procedure :: set_n_bath
        procedure :: calc_nsites => calc_nsites_aim

        procedure, public :: is_periodic => is_periodic_aim

        procedure, public :: get_n_imps
        procedure, public :: get_n_bath
        procedure, public :: is_impurity_site
        procedure, public :: is_bath_site

        procedure, public :: get_impurities
        procedure, public :: get_bath

    end type aim

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

        procedure, public :: dispersion_rel => dispersion_rel_chain

    end type chain

    
    type, extends(lattice) :: cube
        private 

        integer :: length(3) = -1 

    contains 
        private 

        procedure, public :: get_length => get_length_cube
        procedure, public :: is_periodic => is_periodic_cube
        procedure, public :: dispersion_rel => dispersion_rel_cube    

        procedure :: set_length => set_length_cube
        procedure :: calc_nsites => calc_nsites_cube
        procedure :: initialize_sites => init_sites_cube

    end type cube

    type, extends(lattice) :: rectangle 
        private 

        ! how to encode the length? 
        integer :: length(2) = -1
    contains 
        private 

        procedure, public :: get_length => get_length_rect
        procedure, public :: is_periodic => is_periodic_rect
        procedure, public :: dispersion_rel => dispersion_rel_rect

!         procedure, public :: get_length_x, get_length_y

        procedure :: set_length => set_length_rect
        procedure :: calc_nsites => calc_nsites_rect
        procedure :: initialize_sites => init_sites_rect


    end type rectangle

    type, extends(rectangle) :: kagome
        private
    contains
        private

!         procedure, public :: dispersion_rel => dispersion_rel_not_implemented

        procedure :: calc_nsites => calc_nsites_kagome
        procedure :: initialize_sites => init_sites_kagome


    end type kagome

    type, extends(rectangle) :: hexagonal 
        ! i found a unit cell for the hexagonal lattice. but this is 
        ! unfortunately 8 sites big alread.. anyway try it 
        private 

    contains 
        private 

!         procedure, public :: dispersion_rel => dispersion_rel_not_implemented

        procedure :: calc_nsites => calc_nsites_hexagonal 
        procedure :: initialize_sites => init_sites_hexagonal

    end type hexagonal

    type, extends(rectangle) :: triangular
        ! use the length quantitiy and periodicity of the rectangle
        private 

    contains 
        private 

!         procedure, public :: dispersion_rel => dispersion_rel_not_implemented

        ! number of sites is also the same! atleast in this definition of 
        ! the triangular lattice 
        ! so only init_sites must be made new 
        procedure :: initialize_sites => init_sites_triangular

    end type triangular

    type, extends(rectangle) :: tilted
        ! can i just extend rectangle and change the boundary conditions? 
        private

    contains
        private 

        procedure, public :: dispersion_rel => dispersion_rel_tilted
! 
        procedure :: calc_nsites => calc_nsites_tilted
        procedure :: initialize_sites => init_sites_tilted

    end type tilted

    ! can i just extend the chain class to make an impurity chain? 
    ! argh this is annoying without multiple inheritance.. 
    type, extends(aim) :: aim_chain 
        private 
    
        integer :: length = -1 

    contains 
        private 

        procedure, public :: get_length => get_length_aim_chain

!         procedure, public :: dispersion_rel => dispersion_rel_not_implemented

        procedure :: set_length => set_length_aim_chain

        procedure :: initialize_sites => init_sites_aim_chain

    end type aim_chain

    ! also implement a 'star' geometry, especially to deal with the AIM models
    type, extends(lattice) :: star
        private

    contains 
        private 

        procedure :: calc_nsites => calc_nsites_star
        procedure :: set_length => set_length_star
        procedure :: initialize_sites => init_sites_star

        procedure, public :: get_length => get_length_star
        procedure, public :: is_periodic => is_periodic_star

!         procedure, public :: dispersion_rel => dispersion_rel_not_implemented

    end type star

    type, extends(aim) :: aim_star 
        private 

    contains 
        private 

        procedure :: set_length => set_length_aim_star
        procedure :: calc_nsites => calc_nsites_aim_star
        procedure :: initialize_sites => init_sites_aim_star

        procedure, public :: is_periodic => is_periodic_aim_star
        procedure, public :: get_length => get_length_aim_star

!         procedure, public :: dispersion_rel => dispersion_rel_not_implemented

    end type aim_star

    ! write a general cluster aim class:
    ! this is a bit more tricky now, since the setup and the connectivity 
    ! in the impurity sites is input dependent
    ! so from the TMAT in the input i have to define the neighbors of each 
    ! impurity site within themselves 
    ! 
    type, extends(aim_star) :: cluster_aim 
        private 


    contains 
        private 

        procedure :: initialize_sites => init_sites_cluster_aim
!         procedure :: initialize_sites => init_sites_cluster_aim_test

!         procedure, public :: dispersion_rel => dispersion_rel_not_implemented

    end type cluster_aim 

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

        function get_length_t(this, dimen) result(length)
            import :: lattice 
            class(lattice) :: this 
            integer, intent(in), optional :: dimen 
            ! fix our self to 3 dimensions.. thats not good i know, but 
            ! atleast it gives us more flexibility
            integer :: length
        end function get_length_t

        subroutine set_length_t(this, length_x, length_y, length_z)
            import :: lattice 
            class(lattice) :: this 
            integer, intent(in) :: length_x, length_y
            integer, intent(in), optional :: length_z
        end subroutine set_length_t

        function calc_nsites_t(this, length_x, length_y, length_z) result(n_sites) 
            import :: lattice
            class(lattice) :: this 
            integer, intent(in) :: length_x, length_y
            integer, intent(in), optional :: length_z
            integer :: n_sites
        end function calc_nsites_t

        logical function is_periodic_t(this, dimen) 
            import :: lattice
            class(lattice) :: this 
            integer, intent(in), optional :: dimen
        end function is_periodic_t

        subroutine initialize_sites_t(this) 
            import :: lattice 
            class(lattice) :: this 
        end subroutine initialize_sites_t

        function dispersion_rel_t(this, k_vec) result(disp) 
            use constants, only: dp
            import :: lattice 
            class(lattice) :: this 
            integer, intent(in) :: k_vec(3) 
            real(dp) :: disp 
        end function dispersion_rel_t
    end interface

    interface lattice
        procedure lattice_constructor 
    end interface

    interface aim
        procedure aim_lattice_constructor
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
                          STAR_LENGTH = 1, &
                          DIM_RECT = 2, &
                          DIM_CUBE = 3

    ! define the global lattice class here or? this makes more sense
    class(lattice), pointer :: lat

    abstract interface 
        function get_helement_lattice_ex_mat_t(nI, ic, ex, tpar) result(hel)
            use SystemData, only: nel
            use constants, only: dp
            integer, intent(in) :: nI(nel), ic, ex(2,ic)
            logical, intent(in) :: tpar
            HElement_t(dp) :: hel 

        end function get_helement_lattice_ex_mat_t

        function get_helement_lattice_general_t(nI, nJ, ic_ret) result(hel) 
            use SystemData, only: nel 
            use constants, only: dp 
            integer, intent(in) :: nI(nel), nJ(nel) 
            integer, intent(inout), optional :: ic_ret 
            HElement_t(dp) :: hel 

        end function get_helement_lattice_general_t
    end interface 

    procedure(get_helement_lattice_ex_mat_t), pointer, public:: get_helement_lattice_ex_mat
    procedure(get_helement_lattice_general_t), pointer, public :: get_helement_lattice_general

    interface get_helement_lattice
        procedure get_helement_lattice_ex_mat
        procedure get_helement_lattice_general
    end interface get_helement_lattice

contains 

    subroutine init_sites_cluster_aim(this)
        use SystemData, only: nbasis
        use OneEInts, only: tmat2d, gettmatel
        use constants, only: EPS
        class(cluster_aim) :: this 
        character(*), parameter :: this_routine = "init_sites_cluster_aim"


        integer :: ind, orb, i, cnt, connections(nBasis/2), j, maximum

        maximum = -1
        ! now finally implement the actual cluster initializer, which 
        ! reads in a given TMAT and UMAT file, or atleast the information in it

        ! kai wrote a parser, but he has not pushed this information 
        ! and i guess this should be independent of other functionality 
        ! in neci.. so write it here 
        ! set it up so that 
        ASSERT(associated(tmat2d))

        associate(n_bath => this%get_n_bath(), n_imps => this%get_n_imps())
            ASSERT(n_imps + n_bath == nbasis/2)

            if (n_imps > 1) call this%set_ndim(2)

            ! no i have to do this better, loop over all the orbitals 
            ! first the impurities
            do i = 1, n_imps
                ! find all the connections in the tmat 
                connections = -1
                cnt = 0
                do j = 1, nBasis/2
                    ! no diagonal terms
                    if (i == j) cycle
                    if (abs(gettmatel(2*i, 2*j)) > EPS) then 
                        ! cound the number of connections
                        cnt = cnt + 1 
                        ! and maybe i should already store the connected 
                        ! states
                        connections(cnt) = j
                    end if
                    if (cnt > maximum) maximum = cnt
                end do

                this%sites(i) = site(i, cnt, connections(1:cnt-1), 'impurity')
            end do

            do i = n_imps + 1, nBasis/2
                connections = -1 
                cnt = 0
                do j = 1, nbasis/2
                    if (i == j) cycle 
                    if (abs(gettmatel(2*i,2*j)) > EPS) then 
                        cnt = cnt + 1
                        connections(cnt) = j 
                    end if
                    if (cnt > maximum) maximum = cnt
                end do
                this%sites(i) = site(i, cnt, &
                    connections(1:cnt), 'bath')
            end do

            call this%set_nconnect_max(maximum)
        end associate

    end subroutine init_sites_cluster_aim

    subroutine init_sites_cluster_aim_test(this)
        class(cluster_aim) :: this 
        
        integer :: i, imp_neighbors(this%n_imps - 1), j, k, l
        ! this is the tricky part.. i have to use the read-in UMAT and TMAT 
        ! files and setup the connectivity
        ! for now write one dummy function which assumes that all bath sites 
        ! are connected to each impurity and all the impurity sites are 
        ! connected within each-other
        
        ! remember if we have more than one impurity i still have to set 
        ! the dimension and n_connect_max here! 
        
        ! the bath sites are easier here! they are all connected to all 
        ! impurity sites! 
        associate(n_bath => this%get_n_bath(), n_imps => this%get_n_imps())

            if (n_imps > 1) call this%set_ndim(2)

            call this%set_nconnect_max(n_bath + n_imps - 1)

            do i = 1, n_bath
                this%sites(i + n_imps) = site(i + n_imps, n_imps, &
                    [ (i, i = 1, n_imps) ], 'bath')
            end do

            do i = 1, n_imps 
                k = 1
                do j = 1, n_imps
                    if (i /= j) then 
                        imp_neighbors(k) = j
                        k = k + 1
                    end if 
                end do

                this%sites(i) = site(i, n_imps - 1 + n_bath, & 
                    [imp_neighbors, (l + n_imps, l = 1, n_bath)], 'impurity')

            end do
        end associate

        ! for this test function assume every impurity site connected 
        ! within each other 


    end subroutine init_sites_cluster_aim_test

    integer function get_length_aim_star(this, dimen)
        class(aim_star) :: this 
        integer, intent(in), optional :: dimen

        get_length_aim_star = STAR_LENGTH 

    end function get_length_aim_star

    subroutine set_impurity(this, flag) 
        class(site) :: this
        logical, intent(in) :: flag 

        this%t_impurity = flag 

    end subroutine set_impurity

    subroutine set_bath(this, flag)
        class(site) :: this 
        logical, intent(in) :: flag 

        this%t_bath = flag 

    end subroutine set_bath

    logical function is_impurity(this) 
        class(site) :: this 

        is_impurity = this%t_impurity

    end function is_impurity

    logical function is_bath(this) 
        class(site) :: this 

        is_bath = this%t_bath

    end function is_bath

    subroutine init_sites_aim_chain(this)
        class(aim_chain) :: this
        character(*), parameter :: this_routine = "init_sites_aim_chain"
        ! now this is the important routine.. 

        integer :: i

        if (this%get_nsites() < 2) then 
            call stop_all(this_routine, & 
                "less than 2 sites!")
        end if
        ! for the chain we should assert that we only have one impurity! 
        if (this%get_n_imps() > 1) then 
            call stop_all(this_routine, "more than one impurity!")
        end if

        ! the first site is the impurity! 
        this%sites(1) = site(1, 1, [2], 'impurity')

        ! the bath sites are connected within each other, but not periodic!
        do i = 1, this%get_n_bath() - 1
            this%sites(i+1) = site(i + 1, 2, [i, i + 2], 'bath')
        end do

        ! and the last bath site only has one neighbor
        this%sites(this%get_nsites()) = site(this%get_nsites(), 1, &
            [this%get_nsites() - 1], 'bath')

    end subroutine init_sites_aim_chain

    function calc_nsites_aim(this, length_x, length_y, length_z) result(n_sites)
        class(aim) :: this
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        integer :: n_sites
        character(*), parameter :: this_routine = "calc_nsites_aim"

        ! for AIM systems assume first length input is number of impurity 
        ! sites and bath sites are number of path sites per impurity!!
        if (length_x <= 0) then 
            call stop_all(this_routine, &
                "incorrect aim_sites input <= 0!")
        end if
        if (length_y <= 0) then 
            call stop_all(this_routine, & 
                "incorrect bath_sites input <= 0!")
        end if

        n_sites = length_x * length_y + length_x

    end function calc_nsites_aim

    logical function is_bath_site(this, ind) 
        class(aim) :: this 
        integer, intent(in) :: ind 
        character(*), parameter :: this_routine = "is_bath_site" 
        class(site), allocatable :: temp

        ASSERT(ind > 0)
        ASSERT(ind <= this%get_nsites()) 

        is_bath_site = this%sites(ind)%is_bath()

        ! maybe i will use inherited sites(impurity and stuff)

!         allocate(temp, source=this%sites(ind))
! !         temp => this%sites(ind)
! !         associate(sites => this%sites(ind))
!             select type(temp)
! 
!             class is (impurity)
! 
!                 is_bath_site = .false. 
! 
!             class is (bath) 
! 
!                 is_bath_site = .false. 
! 
!             class default 
! 
!                 call stop_all(this_routine, &
!                     "something went wrong.. neither impurity nor bath..")
! 
!             end select
!         end associate

    end function is_bath_site

    logical function is_impurity_site(this, ind)
        class(aim) :: this 
        integer, intent(in) :: ind 
        character(*), parameter :: this_routine = "is_impurity_site"

        ASSERT(ind > 0)
        ASSERT(ind <= this%get_nsites()) 

        is_impurity_site = this%sites(ind)%is_impurity()

        ! just reuse is_bath_site function
!         is_impurity_site = .not. this%is_bath_site(ind) 

    end function is_impurity_site

    function get_bath(this) result(bath_sites) 
        class(aim) :: this 
        integer :: bath_sites(this%n_bath)
        character(*), parameter :: this_routine = "get_bath"

        integer :: i, j
        ! do i store the bath seperately or do i just loop here??

!         bath_sites = -1 
!         j = 1
!         do i = 1, this%get_nsites() 
!             if (this%is_bath_site(i)) then 
!                 bath_sites(j) = this%get_site_index(i)
!                 j = j + 1
!             end if
!         end do
! 
!         ASSERT(.not. any(bath_sites == -1))

        bath_sites = this%bath_sites

    end function get_bath

    function get_impurities(this) result(imp_sites) 
        class(aim) :: this 
        integer :: imp_sites(this%n_imps)
        character(*), parameter :: this_routine = "get_impurities"

        integer :: i, j
! 
!         imp_sites = -1 
! 
!         j = 1
!         do i = 1, this%get_nsites()
!             if (this%is_impurity_site(i)) then 
!                 imp_sites(j) = this%get_site_index(i)
!                 j = j + 1
!             end if
!         end do
! 
!         ASSERT(.not. any(imp_sites == -1))

        ! i guees it is better to store the impurity and the bath 
        ! indices
        imp_sites = this%impurity_sites

    end function get_impurities

    subroutine set_n_imps(this, n_imps)
        class(aim) :: this 
        integer, intent(in) :: n_imps
        character(*), parameter :: this_routine = "set_n_imps"

        ASSERT(n_imps > 0)

        this%n_imps = n_imps 

    end subroutine set_n_imps

    integer function get_n_imps(this) 
        class(aim) :: this 

        get_n_imps = this%n_imps 

    end function get_n_imps

    subroutine set_n_bath(this, n_bath) 
        class(aim) :: this 
        integer, intent(in) :: n_bath 
        character(*), parameter :: this_routine = "set_n_bath"

        ASSERT(n_bath > 0) 

        this%n_bath = n_bath

    end subroutine set_n_bath

    integer function get_n_bath(this) 
        class(aim) :: this 

        get_n_bath = this%n_bath

    end function get_n_bath

    logical function is_k_space(this) 
        class(lattice) :: this 

        is_k_space = this%t_momentum_space

    end function is_k_space 

    ! for the beginning set the aim periodicity to false all the time! 
    logical function is_periodic_aim(this, dimen)
        class(aim) :: this 
        integer, intent(in), optional :: dimen

        ! this function should never get called with dimension input or?
        is_periodic_aim = .false.

    end function is_periodic_aim

    subroutine set_num_neighbors(this, n_neighbors) 
        class(site) :: this 
        integer, intent(in) :: n_neighbors

        this%n_neighbors = n_neighbors 

    end subroutine set_num_neighbors

    pure integer function get_num_neighbors_site(this) 
        class(site), intent(in) :: this 

        get_num_neighbors_site = this%n_neighbors

    end function get_num_neighbors_site

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
        type(site), intent(out) :: lhs 
        type(site), intent(in) :: rhs 

        ! argh.. on every change in the site constructor i have to make 
        ! the change here too.. annoying.. 
        ! but make i OOP style! 
        call lhs%set_index( rhs%get_index() )
        call lhs%set_num_neighbors( rhs%get_num_neighbors() )
        call lhs%allocate_neighbors( rhs%get_num_neighbors() ) 
        call lhs%set_neighbors( rhs%get_neighbors() )
        call lhs%set_impurity( rhs%is_impurity() )
        call lhs%set_bath( rhs%is_bath() )

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

    subroutine init_sites_aim_star(this)
        class(aim_star) :: this 
        character(*), parameter :: this_routine = "init_sites_aim_star" 

        integer :: i 

        if (this%get_nsites() < 2) then 
            call stop_all(this_routine, & 
                "something went wrong: less than 2 sites!") 
        end if

        if (this%get_n_imps() > 1) then 
            call stop_all(this_routine, "more than one impurity not yet implemented!")
        end if

        ! do the first impurity 
        ! the impurity is connected to all the bath sites! 
        this%sites(1) = site(1, this%get_n_bath(), [(i, i = 2, this%get_nsites())], &
            'impurity')

        ! and all the bath sites are just connected to the impurity 
        do i = 2, this%get_nsites() 
            this%sites(i) = site(i, 1, [1], 'bath')
        end do

    end subroutine init_sites_aim_star

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

    subroutine init_sites_cube(this)
        class(cube) :: this
        character(*), parameter :: this_routine = "init_sites_cube" 
        integer :: temp_array(this%length(1),this%length(2),this%length(3))
        integer :: i, x, y, z, temp_neigh(6)
        integer :: up(this%length(1),this%length(2),this%length(3))
        integer :: down(this%length(1),this%length(2),this%length(3))
        integer :: left(this%length(1),this%length(2),this%length(3))
        integer :: right(this%length(1),this%length(2),this%length(3))
        integer :: in(this%length(1),this%length(2),this%length(3))
        integer :: out(this%length(1),this%length(2),this%length(3))
        integer, allocatable :: neigh(:) 

        ! minumum cube size:
        ASSERT(this%get_nsites() >= 8)

        ! encode the lattice with the fortran intrinsic ordering of matrices
        temp_array = reshape([(i, i = 1, this%get_nsites())], & 
            this%length)

        up = cshift(temp_array, -1, 1)
        down = cshift(temp_array, 1, 1)
        right = cshift(temp_array, 1, 2)
        left = cshift(temp_array, -1, 2)
        in = cshift(temp_array, -1, 3)
        out = cshift(temp_array, 1, 3)

        if (this%is_periodic()) then 
            do i = 1, this%get_nsites()
                ! find conversion to 3D matrix indices..
                ! i am not yet sure about that: but seems to work
                x = mod(i-1, this%length(1)) + 1
                z = (i-1)/(this%length(1)*this%length(2)) + 1
                y = mod((i-1)/this%length(1), this%length(2)) + 1

                temp_neigh = [up(x,y,z),down(x,y,z),left(x,y,z),right(x,y,z), &
                                in(x,y,z), out(x,y,z)]

                neigh = sort_unique(temp_neigh)

                this%sites(i) = site(i, size(neigh), neigh)

                deallocate(neigh) 

            end do
        else
            call stop_all(this_routine, &
                "closed boundary conditions not yet implemented for cubic lattice!")
        end if

    end subroutine init_sites_cube

    subroutine init_sites_kagome(this)
        ! the unit cell of the kagome i use is: 
        ! ___  
        ! \./  which is encoded as 1 - 4
        ! /_\                      2/- 5
        !    \                     3/  6
        ! 
        ! see below how this is implemented specifically! 
        class(kagome) :: this 
        character(*), parameter :: this_routine = "init_sites_kagome"

        ! and i think i will for the 6-site, 1x2 and 2x1 12 site and 2x2 24-site 
        ! hard encode the neighbors and stuff because this seems to be a pain 
        ! in the ass! 

        if (this%is_periodic()) then 
            if (this%length(1) == 1 .and. this%length(2) == 1) then 
                ! the smallest cluster 
                this%sites(1) = site(1, 3, [2,4,6])
                this%sites(2) = site(2, 4, [1,3,4,5])
                this%sites(3) = site(3, 3, [2,5,6])
                this%sites(4) = site(4, 3, [1,2,6])
                this%sites(5) = site(5, 3, [2,3,6])
                this%sites(6) = site(6, 4, [1,3,4,5])

            else if (this%length(1) == 1 .and. this%length(2) == 2) then 
                ! the 1x2 cluster with 12-sites: 
                this%sites(1) = site(1, 4, [2,4,10,12])
                this%sites(2) = site(2, 4, [1,3,4,5])
                this%sites(3) = site(3, 4, [2,5,11,12])
                this%sites(4) = site(4, 4, [1,2,6,7])
                this%sites(5) = site(5, 4, [2,3,6,9])
                this%sites(6) = site(6, 4, [4,5,7,9])
                this%sites(7) = site(7, 4, [4,6,8,10])
                this%sites(8) = site(8, 4, [7,9,10,11])
                this%sites(9) = site(9, 4, [5,6,8,11])
                this%sites(10) = site(10, 4, [1,7,8,12])
                this%sites(11) = site(11, 4, [3,8,9,12])
                this%sites(12) = site(12, 4, [1,3,10,11])

            else if (this%length(1) == 2 .and. this%length(2) == 1) then 
                ! the 2x1 12-site cluster: 
                this%sites(1) = site(1, 3, [2,4,12])
                this%sites(2) = site(2, 4, [1,3,4,5]) 
                this%sites(3) = site(3, 3, [2,5,6])
                this%sites(4) = site(4, 3, [1,2,12])
                this%sites(5) = site(5, 3, [2,3,6])
                this%sites(6) = site(6, 4, [3,5,7,10])
                this%sites(7) = site(7, 3, [6,8,10])
                this%sites(8) = site(8, 4, [7,9,10,11])
                this%sites(9) = site(9, 3, [8,11,12])
                this%sites(10) = site(10, 3, [6,7,8])
                this%sites(11) = site(11, 3, [8,9,12])
                this%sites(12) = site(12, 4, [1,4,9,11])

            else if (this%length(1) == 2 .and. this%length(2) == 2) then 
                ! the 2x2 24-site cluster.. puh.. thats a lot to do.. 
                this%sites(1) = site(1, 4, [2,4,10,24])
                this%sites(2) = site(2, 4, [1,3,4,5])
                this%sites(3) = site(3, 4, [2,5,11,12])
                this%sites(4) = site(4, 4, [1,2,7,18])
                this%sites(5) = site(5, 4, [2,3,6,9])
                this%sites(6) = site(6, 4, [5,9,16,19])
                this%sites(7) = site(7, 4, [4,8,10,18])
                this%sites(8) = site(8, 4, [7,9,10,11])
                this%sites(9) = site(9, 4, [5,6,8,11])
                this%sites(10) = site(10, 4, [1,7,8,24])
                this%sites(11) = site(11, 4, [3,8,9,12])
                this%sites(12) = site(12, 4, [3,11,13,22])
                this%sites(13) = site(13, 4, [12,14,16,22])
                this%sites(14) = site(14, 4, [13,15,16,17])
                this%sites(15) = site(15, 4, [14,17,23,24])
                this%sites(16) = site(16, 4, [6,13,14,19])
                this%sites(17) = site(17, 4, [13,15,18,21])
                this%sites(18) = site(18, 4, [4,7,17,21])
                this%sites(19) = site(19, 4, [6,16,20,21,22])
                this%sites(20) = site(20, 4, [19,21,22,23])
                this%sites(21) = site(21, 4, [17,18,20,23])
                this%sites(22) = site(22, 4, [12,13,19,20])
                this%sites(23) = site(23, 4, [15,20,21,24])
                this%sites(24) = site(24, 4, [1,10,15,23])

            else 
                call stop_all(this_routine, &
                    "only 1x1,1x2,2x1 and 2x2 clusters implemented for kagome yet!")
            end if

        else 
            call stop_all(this_routine, & 
                "closed boundary conditions not yet implemented for kagome lattice!")

        end if

    end subroutine init_sites_kagome

    subroutine init_sites_hexagonal(this) 
        ! the way to create the neighboring indices is based on the convention 
        ! that the lattice is interpreted in such a way: the unit cell is:
        !  __
        ! /  \   with encoding: 1 5 
        ! \__/                  2 6 
        ! /  \                  3 7
        !                       4 8

        ! which leads to an lattice:
        ! __/  \__/
        !   \__/  \_
        ! _ /  \__/
        !   \__/  \_
        ! __/  \__/
        !   \__/  \_

        ! so the up and down neighbors are easy to do again. 
        ! just the left and right are different: now each site only either 
        ! has left or right, not both and it alternates
        ! but this also depends on the colum of encoding one is in: 
        !   1 - 5    9 - 13 
        ! - 2   6 - 10   14 - 
        !   3 - 7   11 - 15
        ! - 4   8 - 12   16 - 
        class(hexagonal) :: this 
        integer :: temp_array(4*this%length(1),2*this%length(2))
        integer :: up(4*this%length(1),2*this%length(2))
        integer :: down(4*this%length(1),2*this%length(2)) 
        integer :: right(4*this%length(1),2*this%length(2)) 
        integer :: left(4*this%length(1),2*this%length(2))

        integer :: i , temp_neigh(3), x, y, special
        integer, allocatable :: neigh(:) 
        character(*), parameter :: this_routine = "init_sites_hexagonal"

        temp_array = reshape([(i, i = 1, this%get_nsites())], & 
            [4*this%length(1),2*this%length(2)]) 

        up = cshift(temp_array, -1, 1)
        down = cshift(temp_array, 1, 1) 
        left = cshift(temp_array, -1, 2)
        right = cshift(temp_array, 1, 2) 

        if (this%is_periodic()) then 
            do i = 1, this%get_nsites()
                ! columns and rows:
                x = mod(i-1, 4*this%length(1)) + 1 
                y = (i-1)/(4*this%length(1)) + 1 

                if (is_odd(y)) then 
                    if (is_odd(i)) then 
                        ! every odd number in a odd column has a right neighbor
                        special = right(x,y)
                    else 
                        ! otherwise left 
                        special = left(x,y) 
                    end if
                else 
                    ! for even columns it is the other way around 
                    if (is_odd(i)) then 
                        special = left(x,y) 
                    else 
                        special = right(x,y)
                    end if
                end if
                temp_neigh = [up(x,y), down(x,y), special] 

                neigh = sort_unique(temp_neigh) 

                this%sites(i) = site(i, size(neigh), neigh) 
            end do
        else 
            call stop_all(this_routine, &
                "closed boundary conditions not yet implemented for hexagonal lattice!")

        end if

    end subroutine init_sites_hexagonal

    subroutine init_sites_triangular(this)
        class(triangular) :: this 
        integer :: temp_array(this%length(1),this%length(2))
        integer :: down(this%length(1),this%length(2))
        integer :: up(this%length(1),this%length(2))
        integer :: left(this%length(1),this%length(2))
        integer :: right(this%length(1),this%length(2))
        integer :: lu(this%length(1),this%length(2))
        integer :: rd(this%length(1),this%length(2))
        integer :: i, temp_neigh(6), x, y
        integer, allocatable :: neigh(:)
        character(*), parameter :: this_routine = "init_sites_triangular"
        
        ASSERT(this%get_nsites() >= 4) 

        temp_array = reshape([(i, i = 1, this%get_nsites())], this%length)

        up = cshift(temp_array, -1, 1)
        down = cshift(temp_array, 1, 1)
        right = cshift(temp_array, 1, 2)
        left = cshift(temp_array, -1, 2)
        lu = cshift(up, -1, 2)
        rd = cshift(down, 1, 2) 

        if (this%is_periodic()) then 
            do i = 1, this%get_nsites()

                x = mod(i-1,this%length(1)) + 1
                y = (i-1)/this%length(1) + 1
            
                temp_neigh = [up(x,y),down(x,y),left(x,y),right(x,y),lu(x,y),rd(x,y)]

                neigh = sort_unique(temp_neigh)

                this%sites(i) = site(i, size(neigh), neigh)

                deallocate(neigh)

            end do
        else
            call stop_all(this_routine, &
                "closed boundary conditions not yet implemented for triangular lattice!")
        end if

    end subroutine init_sites_triangular

    subroutine init_sites_rect(this)
        class(rectangle) :: this 
        character(*), parameter :: this_routine = "init_sites_rect"

        integer :: i, temp_neigh(4), x, y, min_val, max_val, j
        integer :: temp_array(this%length(1),this%length(2))
        integer :: down(this%length(1),this%length(2))
        integer :: up(this%length(1),this%length(2))
        integer :: left(this%length(1),this%length(2))
        integer :: right(this%length(1),this%length(2))
        integer :: unique(4)
        integer, allocatable :: neigh(:)
        integer :: sort_array_3(3), sort_array_2(2), sort_array(4)

        ! this is the important routine.. 
        ! store lattice like that: 
        ! 1 4 7
        ! 2 5 8
        ! 3 6 9

        ASSERT(this%get_nsites() >= 4)

        ! use cshift intrinsic of fortran.. 
        ! how do i efficiently set that up? 
        temp_array =  reshape([(i, i = 1, this%get_nsites())], &
            this%length)

        up = cshift(temp_array, -1, 1)
        down = cshift(temp_array, 1, 1)
        right = cshift(temp_array, 1, 2)
        left = cshift(temp_array, -1, 2)

        if (this%is_periodic()) then 

            do i = 1, this%get_nsites() 
                ! create the neighbor list 
                x = mod(i-1,this%length(1)) + 1
                y = (i-1)/this%length(1) + 1

!                 print *, "i, (x,y): ", i, x, y
                temp_neigh = [up(x,y),down(x,y),left(x,y),right(x,y)]

                neigh = sort_unique(temp_neigh)

                this%sites(i) = site(i, size(neigh), neigh)

                deallocate(neigh)

            end do

        else if (this%is_periodic(1)) then 
            ! only periodic in x-direction
            do i = 1, this%get_nsites()

                x = mod(i-1,this%length(1)) + 1
                y = (i-1)/this%length(1) + 1

                ! now definetly always take the left and right neighbors 
                ! but up and down if we are not jumping boundaries 
                if (x  == 1) then 
                    ! dont take upper neighbor -> just repeat an neighbor, 
                    ! which will get removed from unique
                    temp_neigh  = [down(x,y), down(x,y), left(x,y), right(x,y)]
                else if (x == this%length(1)) then 
                    temp_neigh  = [up(x,y), up(x,y), left(x,y), right(x,y)]

                else
                    ! take all
                    temp_neigh = [up(x,y),down(x,y),left(x,y),right(x,y)]
                end if

                neigh = sort_unique(temp_neigh)

                this%sites(i) = site(i, size(neigh), neigh)

                deallocate(neigh)
            end do

        else if (this%is_periodic(2)) then 
            ! only periodic in the y-direction
            do i = 1, this%get_nsites()

                x = mod(i-1,this%length(1)) + 1
                y = (i-1)/this%length(1) + 1

                ! now definetly always take the left and right neighbors 
                ! but up and down if we are not jumping boundaries 
                if (y  == 1) then 
                    ! dont take upper neighbor -> just repeat an neighbor, 
                    ! which will get removed from unique
                    temp_neigh  = [up(x,y), down(x,y), right(x,y), right(x,y)]
                else if (y == this%length(2)) then 
                    temp_neigh  = [up(x,y), down(x,y), left(x,y), left(x,y)]

                else
                    ! take all
                    temp_neigh = [up(x,y),down(x,y),left(x,y),right(x,y)]
                end if

                neigh = sort_unique(temp_neigh)

                this%sites(i) = site(i, size(neigh), neigh)

                deallocate(neigh)
            end do

        else 
            ! non-periodic
            do i = 1, this%get_nsites()

                x = mod(i-1,this%length(1)) + 1
                y = (i-1)/this%length(1) + 1

                ! now definetly always take the left and right neighbors 
                ! but up and down if we are not jumping boundaries 
                if (x == 1 .and. y  == 1) then 
                    ! dont take upper neighbor -> just repeat an neighbor, 
                    ! which will get removed from unique
                    temp_neigh  = [down(x,y), down(x,y), right(x,y), right(x,y)]
                else if (y == this%length(2) .and. x == this%length(1)) then 
                    temp_neigh  = [up(x,y), up(x,y), left(x,y), left(x,y)]

                else if (x == this%length(1) .and. y == 1) then 
                    temp_neigh = [up(x,y),up(x,y),right(x,y),right(x,y)]
                else if (y == this%length(2) .and. x == 1) then 
                    temp_neigh = [down(x,y),down(x,y),left(x,y),left(x,y)]

                else if (x == 1) then 
                    temp_neigh = [left(x,y), right(x,y), down(x,y),down(x,y)]

                else if (x == this%length(1)) then 
                    temp_neigh = [left(x,y),right(x,y),up(x,y),up(x,y)]

                else if (y == 1) then 
                    temp_neigh = [right(x,y),right(x,y),up(x,y),down(x,y)]

                else if (y == this%length(2)) then 
                    temp_neigh = [left(x,y),left(x,y),up(x,y),down(x,y)]

                else
                    ! take all
                    temp_neigh = [up(x,y),down(x,y),left(x,y),right(x,y)]
                end if

                neigh = sort_unique(temp_neigh)

                this%sites(i) = site(i, size(neigh), neigh)

                deallocate(neigh)
            end do

        end if

    end subroutine init_sites_rect

    subroutine init_sites_tilted(this)
        class(tilted) :: this 
        character(*), parameter :: this_routine = "init_sites_tilted"

        integer :: temp_array(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: up(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: down(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: left(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: right(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: right_ul(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: right_ur(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: right_dl(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: right_dr(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: right_rr(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: right_ll(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: up_ul(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: up_ur(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: up_dl(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: up_dr(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: up_rr(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: up_ll(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: down_ul(-this%length(1):this%length(1),&
              -this%length(1):this%length(1)+1)
        integer :: down_ur(-this%length(1):this%length(1),&
              -this%length(1):this%length(1)+1)
        integer :: down_dl(-this%length(1):this%length(1),&
              -this%length(1):this%length(1)+1)
        integer :: down_dr(-this%length(1):this%length(1),&
             -this%length(1):this%length(1)+1)
        integer :: down_rr(-this%length(1):this%length(1),&
              -this%length(1):this%length(1)+1)
        integer :: down_ll(-this%length(1):this%length(1),&
              -this%length(1):this%length(1)+1)
        integer :: left_ul(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: left_ur(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: left_dl(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: left_dr(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: left_rr(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: left_ll(-this%length(1):this%length(1),&
                              -this%length(1):this%length(1)+1)
        integer :: i, j, k, l, pbc, temp_neigh(4)
        integer :: right_nn, left_nn, up_nn, down_nn
        integer, allocatable :: neigh(:)
        ! convention of lattice storage: 
        ! 
        !   2 5
        ! 1 3 6 8
        !   4 7

        ASSERT(this%get_nsites() >= 8) 

        ! set up the lattice indices, via the use of "k-vectors"
        temp_array(:,:) = 0

        k = 0
        l = 1
        do i = -this%length(1)+1, 0

            do j = -k, k

                temp_array(j,i) = l 

                l = l + 1
            end do
            k = k + 1
        end do

        k = k - 1
        do i = 1, this%length(1)

            do j = -k, k

                temp_array(j,i) = l

                l = l + 1

            end do

            k = k - 1
        end do

        up = cshift(temp_array, -1, 1)
        down = cshift(temp_array, 1, 1)
        right = cshift(temp_array, 1, 2)
        left = cshift(temp_array, -1, 2)

        ! apply the periodic boundary conditions to the neighbors
        pbc = this%length(1)
        up_ur = cshift(cshift(up, -pbc, 1), pbc, 2)
        up_dr= cshift(cshift(up, pbc, 1), pbc, 2)
        up_ul = cshift(cshift(up, -pbc, 1), -pbc, 2)
        up_dl = cshift(cshift(up, pbc, 1), -pbc, 2)
        up_rr = cshift(up, 2*pbc, 2)
        up_ll = cshift(up, -2*pbc, 2)


        down_ur = cshift(cshift(down, -pbc, 1), pbc, 2)
        down_dr= cshift(cshift(down, pbc, 1), pbc, 2)
        down_ul = cshift(cshift(down, -pbc, 1), -pbc, 2)
        down_dl = cshift(cshift(down, pbc, 1), -pbc, 2)
        down_rr = cshift(down, 2*pbc, 2)
        down_ll = cshift(down, -2*pbc, 2)

        right_ur = cshift(cshift(right, -pbc, 1), pbc, 2)
        right_dr= cshift(cshift(right, pbc, 1), pbc, 2)
        right_ul = cshift(cshift(right, -pbc, 1), -pbc, 2)
        right_dl = cshift(cshift(right, pbc, 1), -pbc, 2)
        right_rr = cshift(right, 2*pbc, 2)
        right_ll = cshift(right, -2*pbc, 2)

        left_ur = cshift(cshift(left, -pbc, 1), pbc, 2)
        left_dr= cshift(cshift(left, pbc, 1), pbc, 2)
        left_ul = cshift(cshift(left, -pbc, 1), -pbc, 2)
        left_dl = cshift(cshift(left, pbc, 1), -pbc, 2)
        left_rr = cshift(left, 2*pbc, 2)
        left_ll = cshift(left, -2*pbc, 2)

        k = 0
        l = 1
        ! now get the neighbors
        if (this%is_periodic()) then
            ! fully periodic case 
            do i = -this%length(1)+1, 0
                do j = -k, k 
                    ! make the neigbors list 
                    up_nn = maxval([up(j,i), up_ur(j,i), up_dr(j,i), up_ul(j,i), &
                                    up_dl(j,i)])

                    if (up_nn == 0) then
                        up_nn = maxval([up_rr(j,i), up_ll(j,i)])
                        if (up_nn == 0) then
                            print *, "smth wrong!"
                        end if
                    end if
                        
                    down_nn = maxval([down(j,i), down_ur(j,i), down_dr(j,i), down_ul(j,i), &
                        down_dl(j,i)])
                    
                    if (down_nn == 0) then 
                        down_nn = maxval([down_rr(j,i), down_ll(j,i)])

                        if (down_nn == 0) then 
                            print *, "smth wrong!"
                        end if
                    end if

                    right_nn = maxval([right(j,i), right_ur(j,i), right_dr(j,i), right_ul(j,i), &
                        right_dl(j,i)])
                    
                    if (right_nn == 0) then 
!                         right_nn = maxval([right_rr(j,i), right_ll(j,i)])
                        right_nn = right_ll(j,i)
                        if (right_nn == 0) then 
                            print *, "smth wrong!"
                        end if
                    end if

                    left_nn = maxval([left(j,i), left_ur(j,i), left_dr(j,i), left_ul(j,i), &
                        left_dl(j,i)])
                    
                    if (left_nn == 0) then 
!                         left_nn = maxval([left_rr(j,i), left_ll(j,i)])
                        left_nn = left_rr(j,i)
                        if (left_nn == 0) then 
                            print *, "smth wrong!"
                        end if
                    end if

                
                    neigh = sort_unique([up_nn, down_nn, left_nn, right_nn])

                    this%sites(l) = site(l, size(neigh), neigh)

                    l = l + 1

                    deallocate(neigh)

                end do
                k = k + 1
            end do

            k = k - 1
            do i = 1, this%length(1)
                do j = -k, k 
                    ! make the neigbors list 
                    up_nn = maxval([up(j,i), up_ur(j,i), up_dr(j,i), up_ul(j,i), &
                                    up_dl(j,i)])

                    if (up_nn == 0) then
                        up_nn = maxval([up_rr(j,i), up_ll(j,i)])
                        if (up_nn == 0) then
                            print *, "smth wrong!"
                        end if
                    end if
                        
                    down_nn = maxval([down(j,i), down_ur(j,i), down_dr(j,i), down_ul(j,i), &
                        down_dl(j,i)])
                    
                    if (down_nn == 0) then 
                        down_nn = maxval([down_rr(j,i), down_ll(j,i)])

                        if (down_nn == 0) then 
                            print *, "smth wrong!"
                        end if
                    end if

                    right_nn = maxval([right(j,i), right_ur(j,i), right_dr(j,i), right_ul(j,i), &
                        right_dl(j,i)])
                    
                    if (right_nn == 0) then 
!                         right_nn = maxval([right_rr(j,i), right_ll(j,i)])
                        right_nn = right_ll(j,i)
                        if (right_nn == 0) then 
                            print *, "smth wrong!"
                        end if
                    end if

                    left_nn = maxval([left(j,i), left_ur(j,i), left_dr(j,i), left_ul(j,i), &
                        left_dl(j,i)])
                    
                    if (left_nn == 0) then 
!                         left_nn = maxval([left_rr(j,i), left_ll(j,i)])
                        left_nn = left_rr(j,i)
                        if (left_nn == 0) then 
                            print *, "smth wrong!"
                        end if
                    end if

                    neigh = sort_unique([up_nn, down_nn, left_nn, right_nn])

                    this%sites(l) = site(l, size(neigh), neigh)

                    l = l + 1

                    deallocate(neigh)

                end do
                k = k -1
            end do
        else if (this%is_periodic(1)) then 
            ! only apply (x,x) periodicity 
            do i = -this%length(1) + 1, 0
                do j = -k, k

                    up_nn = maxval([up(j,i),up_ur(j,i),up_dl(j,i)])
                    down_nn = maxval([down(j,i),down_ur(j,i),down_dl(j,i)])
                    left_nn = maxval([left(j,i),left_ur(j,i),left_dl(j,i)])
                    right_nn = maxval([right(j,i),right_ur(j,i),right_dl(j,i)])

                    temp_neigh = [up_nn, down_nn, left_nn, right_nn]

                    neigh = sort_unique(pack(temp_neigh, temp_neigh > 0))

                    this%sites(l) = site(l, size(neigh), neigh)

                    l = l + 1 

                    deallocate(neigh)
                end do
                k = k + 1
            end do
            k = k - 1 
            
            do i = 1, this%length(1)
                do j = -k, k

                    up_nn = maxval([up(j,i),up_ur(j,i),up_dl(j,i)])
                    down_nn = maxval([down(j,i),down_ur(j,i),down_dl(j,i)])
                    left_nn = maxval([left(j,i),left_ur(j,i),left_dl(j,i)])
                    right_nn = maxval([right(j,i),right_ur(j,i),right_dl(j,i)])

                    temp_neigh = [up_nn, down_nn, left_nn, right_nn]

                    neigh = sort_unique(pack(temp_neigh, temp_neigh > 0))

                    this%sites(l) = site(l, size(neigh), neigh)

                    l = l + 1 

                    deallocate(neigh)
                end do
                k = k - 1
            end do

        else if (this%is_periodic(2)) then 
            ! only apply (x,-x) periodicity
            do i = -this%length(1) + 1, 0
                do j = -k, k

                    up_nn = maxval([up(j,i),up_ul(j,i),up_dr(j,i)])
                    down_nn = maxval([down(j,i),down_ul(j,i),down_dr(j,i)])
                    left_nn = maxval([left(j,i),left_ul(j,i),left_dr(j,i)])
                    right_nn = maxval([right(j,i),right_ul(j,i),right_dr(j,i)])

                    temp_neigh = [up_nn, down_nn, left_nn, right_nn]

                    neigh = sort_unique(pack(temp_neigh, temp_neigh > 0))

                    this%sites(l) = site(l, size(neigh), neigh)

                    l = l + 1 

                    deallocate(neigh)
                end do
                k = k + 1
            end do
            k = k - 1 
            
            do i = 1, this%length(1)
                do j = -k, k

                    up_nn = maxval([up(j,i),up_ul(j,i),up_dr(j,i)])
                    down_nn = maxval([down(j,i),down_ul(j,i),down_dr(j,i)])
                    left_nn = maxval([left(j,i),left_ul(j,i),left_dr(j,i)])
                    right_nn = maxval([right(j,i),right_ul(j,i),right_dr(j,i)])

                    temp_neigh = [up_nn, down_nn, left_nn, right_nn]

                    neigh = sort_unique(pack(temp_neigh, temp_neigh > 0))

                    this%sites(l) = site(l, size(neigh), neigh)

                    l = l + 1 

                    deallocate(neigh)
                end do
                k = k - 1
            end do

        else
            ! non-periodic case
            do i = -this%length(1) + 1, 0 
                do j = -k,k 
                    ! only a neighbor if the index is non-zero! 
                    temp_neigh = [up(j,i), down(j,i), left(j,i), right(j,i)] 

                    neigh = sort_unique(pack(temp_neigh, temp_neigh > 0)) 

                    this%sites(l) = site(l, size(neigh), neigh)

                    l = l + 1

                    deallocate(neigh)

                end do
                k = k + 1
            end do

            k = k - 1
            do i = 1, this%length(1)
                do j = -k, k
                    ! only a neighbor if the index is non-zero! 
                    temp_neigh = [up(j,i), down(j,i), left(j,i), right(j,i)] 

                    neigh = sort_unique(pack(temp_neigh, temp_neigh > 0)) 

                    this%sites(l) = site(l, size(neigh), neigh)

                    l = l + 1

                    deallocate(neigh)

                end do
                k = k - 1
            end do
        end if


    end subroutine init_sites_tilted

    function dispersion_rel_chain(this, k_vec) result(disp) 
        class(chain) :: this 
        integer, intent(in) :: k_vec(3) 
        real(dp) :: disp 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "dispersion_rel_chain" 
#endif


        ! for now only do it for periodic boundary conditions.. 
        ASSERT(this%is_periodic())

        ! and for now only for nearest neighbor interaction! 
        ! although this is just the nearest neighbor band.. 
        ! for next nearest and additional function should be implemented!

        ! i need to bring in the length of the chain and stuff.. 
        ! and i should consider twisted boundary conditions and nearest 
        ! neigbhors here too..? i think so.. 
        disp = 2.0_dp * cos(2*pi*(k_vec(1) + twisted_bc(1))/this%length)

    end function dispersion_rel_chain

    function dispersion_rel_rect(this, k_vec) result(disp)
        class(rectangle) :: this 
        integer, intent(in) :: k_vec(3)
        real(dp) :: disp 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "dispersion_rel_rect"
#endif

        ASSERT(this%is_periodic())

        disp = 2.0_dp * (cos(2*pi*(k_vec(1)+twisted_bc(1))/this%length(1)) & 
                        +cos(2*pi*(k_vec(2)+twisted_bc(2))/this%length(2)))

    end function dispersion_rel_rect

    function dispersion_rel_cube(this, k_vec) result(disp) 
        class(cube) :: this 
        integer, intent(in) :: k_vec(3)
        real(dp) :: disp 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "dispersion_rel_cube" 
#endif 

        ASSERT(this%is_periodic())

        disp = 2.0_dp * (cos(2*pi*(k_vec(1)+twisted_bc(1))/this%length(1)) &
                        +cos(2*pi*(k_vec(2)+twisted_bc(2))/this%length(2)) & 
                        +cos(2*pi*(k_vec(3)+twisted_bc(3))/this%length(3)))

    end function dispersion_rel_cube

    function dispersion_rel_tilted(this, k_vec) result(disp) 
        class(tilted) :: this 
        integer, intent(in) :: k_vec(3)
        real(dp) :: disp
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "dispersion_rel_tilted"
#endif

        ASSERT(this%is_periodic())

        disp = 2.0_dp * (cos(2*pi*((k_vec(1)+twisted_bc(1))*this%length(1) &
                                  +(k_vec(2)+twisted_bc(2))*this%length(2))/this%n_sites) & 
                        +cos(2*pi*((k_vec(1)+twisted_bc(1))*this%length(2) & 
                                  -(k_vec(2)+twisted_bc(2))*this%length(1))/this%n_sites))

    end function dispersion_rel_tilted

    function dispersion_rel_not_implemented(this, k_vec) result(disp)
        class(lattice) :: this 
        integer, intent(in) :: k_vec(3) 
        real(dp) :: disp 
        character(*), parameter :: this_routine = "dispersion_rel"

        call stop_all(this_routine, &
            "dispersion relation not yet implemented for this lattice type!")

    end function dispersion_rel_not_implemented

    function sort_unique(list) result(output)
        integer, intent(in) :: list(:)
        integer, allocatable :: output(:)

        integer :: i, min_val,  max_val, unique(size(list))
!         integer, allocatable :: unique(:)

!         allocate(unique(size(list)))

        unique = 0
        i = 0 
        min_val = minval(list) - 1
        max_val = maxval(list) 

        do while (min_val < max_val) 
            i = i + 1 
            min_val = minval(list, mask=list>min_val) 
            unique(i) = min_val 
        end do
        allocate(output(i), source=unique(1:i))

    end function sort_unique

    subroutine init_aim(this, length_x, length_y) 
        class(aim) :: this 
        integer, intent(in) :: length_x, length_y
        character(*), parameter :: this_routine = "init_aim"

    end subroutine init_aim

    subroutine init_lattice(this, length_x, length_y, length_z, & 
            t_periodic_x, t_periodic_y, t_periodic_z)
        ! and write the first dummy initialize 
        class(lattice) :: this 
        integer, intent(in) :: length_x, length_y, length_z
        logical, intent(in) :: t_periodic_x, t_periodic_y, t_periodic_z
        character(*), parameter :: this_routine = "init_lattice"

        integer :: n_sites, i

        n_sites = this%calc_nsites(length_x, length_y, length_z)

        ! and for the rest i can call general routines:
        call this%set_nsites(n_sites) 
        call this%set_periodic(t_periodic_x, t_periodic_y, t_periodic_z)

        select type (this) 

        ! well i cannot init type is (lattice) if i choose to make it 
        ! abstract. since it is not allowed to ever be intantiated..

        class is (chain)
            ! set some defaults for the chain lattice type 
            call this%set_ndim(DIM_CHAIN)

            call this%set_length(length_x, length_y)
            ! if incorrect length input it is caught in the calc_nsites above..
            if (this%get_length() == 1) then 
                call this%set_nconnect_max(0)
            else if (this%get_length() == 2 .and. (.not. this%is_periodic())) then 
                call this%set_nconnect_max(1)
            else 
                call this%set_nconnect_max(N_CONNECT_MAX_CHAIN) 
            end if

            ! the type specific routine deal with the check of the 
            ! length! 

            ! i should not call this set_nsites since this really should 
            ! just imply that it set the variable 
            ! introduce a second routine, which first determines the 
            ! number of sites depending on the lattice type 

        class is (rectangle) 

            call this%set_ndim(DIM_RECT)
            call this%set_length(length_x, length_y)

            if (this%get_length(1) == 2 .and. this%get_length(2) == 2) then 
                if (.not. this%is_periodic(1) .and. this%is_periodic(2)) then 
                    call this%set_nconnect_max(3)
                else if (.not. this%is_periodic(2) .and. this%is_periodic(1)) then 
                    call this%set_nconnect_max(3)
                else if (this%is_periodic()) then 
                    call this%set_nconnect_max(4)
                else if (.not. this%is_periodic()) then 
                    call this%set_nconnect_max(2)
                end if 
            else 
                call this%set_nconnect_max(4)
            end if 

        class is (tilted)
            
            call this%set_ndim(DIM_RECT) 
            call this%set_length(length_x, length_y) 
            call this%set_nconnect_max(4) 

        class is (cube) 
            call this%set_ndim(DIM_CUBE)
            call this%set_length(length_x, length_y, length_z)
            call this%set_nconnect_max(6)

        class is (triangular) 
            call this%set_ndim(DIM_RECT)
            call this%set_length(length_x, length_y)
            ! for a filling with triangles the maximum connection is 6! 

            call this%set_nconnect_max(6)

        class is (hexagonal) 

            call this%set_ndim(DIM_RECT) 
            call this%set_length(length_x, length_y, length_z) 
            call this%set_nconnect_max(3) 

        class is (kagome) 
            call this%set_ndim(DIM_RECT) 
            call this%set_length(length_x, length_y, length_z) 
            call this%set_nconnect_max(4)

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

        class is (aim_chain) 

            ! do stuff
            call this%set_ndim(DIM_CHAIN)
            call this%set_length(length_x, length_y)
            ! the neighbors is a bit complicated in this case.. 
            ! although it is a chain.. so it should not have more than 
            ! one impurity! 
            if (length_x > 1) then 
                call stop_all(this_routine, &
                    "more than 1 impurity taken in tha aim_chain setup!")
            end if
            if (length_y == 1) then 
                call this%set_nconnect_max(1) 
            else 
                call this%set_nconnect_max(N_CONNECT_MAX_CHAIN)
            end if
            
            ! i can use class specific routine in this block
            call this%set_n_imps(length_x)
            call this%set_n_bath(length_y)
            
            allocate(this%impurity_sites(length_x))
            allocate(this%bath_sites(length_y))

            this%impurity_sites = [(i, i = 1, length_x)]
            this%bath_sites = [(length_x + i, i = 1, length_y)]

        class is (aim_star) 
            ! this is the star with only 1 impurity for now.. 
            ! i still have to think how to efficiently setup up a 
            ! cluster impurity.. 
            ! i guess i have to decide on a lattice and a ab-initio 
            ! cluster impurity! thats good yeah 

            call this%set_ndim(DIM_STAR)
            ! number of bath sites is the maximal connectivity
            call this%set_nconnect_max(length_y) 

            ! also the one-site impurity can't be periodic! 
            if (t_periodic_x .or. t_periodic_y) then 
                call stop_all(this_routine, & 
                    "incorrect initialization info: requested periodic 'star' geometry!")
            end if

            if (length_x > 1) then 
                call stop_all(this_routine, &
                    "aim_star only implemented for one impurity!")
            end if

            call this%set_n_imps(length_x)
            call this%set_n_bath(length_y) 
 
            allocate(this%impurity_sites(length_x))
            allocate(this%bath_sites(length_y))

            this%impurity_sites = [(i, i = 1, length_x)]
            this%bath_sites = [(length_x + i, i = 1, length_y)]

        class is (cluster_aim) 

            ! now we have to be more specific.. 
            ! do i have to use the inputs already here? 
            if (length_x <= 0) then 
                call stop_all(this_routine, & 
                    "zero or negative impurity sites requested!")
            end if
            if (length_y <= 0) then 
                call stop_all(this_routine, & 
                    "zero or negative bath orbitals requested!")
            end if

            call this%set_n_imps(length_x)
            call this%set_n_bath(length_y) 

            if (length_x == 1) then 
                call this%set_ndim(DIM_STAR) 
                call this%set_nconnect_max(length_y)
            else 
                ! otherwise i have to do this later.. 
            end if
 
            allocate(this%impurity_sites(length_x))
            allocate(this%bath_sites(length_y))

            this%impurity_sites = [(i, i = 1, length_x)]
            this%bath_sites = [(length_x + i, i = 1, length_y)]

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

    function aim_lattice_constructor(lat_type, length_x, length_y)  result(this)
        character(*), intent(in) :: lat_type 
        integer, intent(in) :: length_x, length_y 
        class(aim), pointer :: this 
        character(*), parameter :: this_routine = "aim_lattice_constructor"

        select case(lat_type)
        case('chain','aim-chain', 'chain-aim')

            allocate(aim_chain :: this) 

        case ('star', 'aim-star', 'star-aim')

            allocate(aim_star :: this)

        case('cluster', 'aim-cluster', 'cluster-aim')

            allocate(cluster_aim :: this)

        case default 
            ! stop here because a incorrect lattice type was given 
            call stop_all(this_routine, & 
                'incorrect lattice type provided in lattice_constructor!')

        end select 

        ! the initializer deals with the different types then.. 
!         call this%initialize(length_x, length_y)
        call this%initialize(length_x, length_y, 1, .false., .false., .false.)

    end function aim_lattice_constructor

    function lattice_constructor(lattice_type, length_x, length_y, length_z, t_periodic_x , &
            t_periodic_y, t_periodic_z, space) result(this)
        ! write a general public lattice_constructor for lattices 
        ! the number of inputs are still undecided.. do we always have 
        ! the same number or differing number of inputs? 
        ! i guess, since we will be using global variables which are either 
        ! read in or set as default to init it.. 
        character(*), intent(in) :: lattice_type
        integer, intent(in) :: length_x, length_y, length_z
        logical, intent(in) :: t_periodic_x, t_periodic_y, t_periodic_z
        character(*), intent(in), optional :: space
        class(lattice), pointer :: this 
        character(*), parameter :: this_routine = "lattice_constructor"


        select case (lattice_type) 
        case ('chain') 

            allocate(chain :: this) 

        case ('star') 

            allocate(star :: this) 

!         case ('aim-chain') 
! 
!             allocate(aim_chain :: this)

        case ('square')
            ! i guess i want to make a seperate case for the tilted 
            ! square, although just the boundary conditions change, but also 
            ! the length input changes
            if (length_x /= length_y) then 
                call stop_all(this_routine, &
                    "incorrect length input for square lattice!")

            end if

            allocate(rectangle :: this)

        case('rectangle')

            allocate(rectangle :: this)

        case('tilted','tilted-square','square-tilted')

            if (length_x /= length_y) then 
                call stop_all(this_routine, &
                    "incorrect length_x /= length_y input for tilted lattice!")
            end if

            allocate(tilted :: this)

        case ('cube','cubic') 
            ! for the sake of no better name also use "cube" even if the 
            ! sides are not the same length

            if (any([length_x, length_y, length_z] < 2)) then 
                call stop_all(this_routine, &
                    "too short cube side lengths < 2!")
            end if

            allocate(cube :: this)

        case ('triangular','triangle')

            if (any([length_x, length_y] < 2)) then 
                call stop_all(this_routine, &
                    "too short lengths for triangular lattice! < 2!")
            end if

            allocate( triangular :: this )

        case ('hexagonal', 'hex', 'hexagon', 'honeycomb') 

            allocate( hexagonal :: this )

        case ('kagome') 
            allocate( kagome :: this) 

        case default 
            ! stop here because a incorrect lattice type was given 
            call stop_all(this_routine, & 
                'incorrect lattice type provided in lattice_constructor!')

        end select 

        call this%set_name(lattice_type)

        ! depending on the string input defining lattice type 
        ! initialize corresponding lattice 
        if (present(space)) then 
            select case(space) 
            case ('k-space') 
                this%t_momentum_space = .true. 

            case ('real-space')
                this%t_momentum_space = .false.

            case default 
                call stop_all(this_routine, "not recognized space!")

            end select 
        end if

        ! the initializer deals with the different types then.. 
        call this%initialize(length_x, length_y, length_z, &
            t_periodic_x, t_periodic_y, t_periodic_z)

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
        type(site) :: this 
        character(*), parameter :: this_routine = "site_constructor"

!         allocate(site :: this) 
        if (present(site_type)) then 
            ! not yet implementetd or to do.. so wait.. 
            select case (site_type)

            case ('impurity', 'imp')
                
                call this%set_impurity(.true.)

            case ('bath')

                call this%set_bath(.true.)

            case default 
                call stop_all(this_routine, &
                    "incorrect site type provided") 
                
            end select 

        else 
            ! this is the default case 

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

    subroutine set_name(this, lat_type)
        class(lattice) :: this
        character(*), intent(in) :: lat_type 

        this%name = lat_type

    end subroutine set_name

    function get_name(this) result(lattice_name) 
        class(lattice) :: this
        character(NAME_LEN) :: lattice_name

        lattice_name = this%name

    end function get_name

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

        select type (this)
        class is (aim) 
            deallocate(this%impurity_sites)
            deallocate(this%bath_sites)
        end select 

    end subroutine lattice_deconstructor
    
    subroutine aim_deconstructor(this) 
        class(aim), pointer :: this 

        call this%deallocate_sites()

        nullify(this)

    end subroutine aim_deconstructor

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

    function get_num_neighbors_lattice(this,ind) result(n_neighbors)
        class(lattice) :: this
        integer, intent(in) :: ind 
        integer :: n_neighbors
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_num_neighbors_lattice"
#endif
       
        ! make all assert on a seperate line, so we exactly know what is 
        ! going wrong.. 
        ASSERT(ind <= this%get_nsites())
        ASSERT(ind > 0)
        ASSERT(allocated(this%sites))
        ASSERT(allocated(this%sites(ind)%neighbors))

        n_neighbors = this%sites(ind)%get_num_neighbors()

    end function get_num_neighbors_lattice

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

    function get_spinorb_neighbors_lat(this, spinorb) result(neighbors)
        class(lattice) :: this 
        integer, intent(in) :: spinorb
        integer, allocatable :: neighbors(:) 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_spinorb_neighbors_lat"
#endif

        ASSERT(spinorb <= 2*this%get_nsites())
        ASSERT(spinorb > 0) 
        ASSERT(allocated(this%sites))
        ASSERT(allocated(this%sites(gtid(spinorb))%neighbors))

        neighbors = this%get_neighbors(gtid(spinorb))

        if (is_beta(spinorb)) then
            neighbors = 2 * neighbors - 1
        else 
            neighbors = 2 * neighbors
        end if

    end function get_spinorb_neighbors_lat

    function calc_nsites_aim_star(this, length_x, length_y, length_z) result(n_sites) 
        ! for the star geometry with maybe 
        class(aim_star) :: this 
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        integer :: n_sites
        character(*), parameter :: this_routine = "calc_nsites_aim_star"

        if (length_x < 1) then 
            call stop_all(this_routine, "n_imps < 1!") 
        end if
        if (length_y < 1) then 
            call stop_all(this_routine, "n_bath < 1!")
        end if

        n_sites = length_x + length_y

    end function calc_nsites_aim_star 

    function calc_nsites_star(this, length_x, length_y, length_z) result(n_sites) 
        ! the maximum of the input is used as the n_sites parameter! 
        ! this is the same function as the one for "chain" below.. 
        ! but i cannot use it somehow.. 
        class(star) :: this 
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        integer :: n_sites 
        character(*), parameter :: this_routine = "calc_nsites_star"

        if (max(length_x,length_y) < 1 .or. min(length_x, length_y) > 1 .or. & 
            min(length_x,length_y) < 0) then 
            n_sites = -1 
            call stop_all(this_routine, "something went wrong in length input!")

        else 
            n_sites = max(length_x, length_y)

        end if


    end function calc_nsites_star

    logical function is_periodic_aim_star(this, dimen)
        class(aim_star) :: this 
        integer, intent(in), optional :: dimen

        is_periodic_aim_star = .false.

    end function is_periodic_aim_star

    function calc_nsites_chain(this, length_x, length_y, length_z) result(n_sites) 
        ! i acually do not want to rely on the previous calculated 
        ! length object in the class, since it is easy to recalc from 
        ! here and so i remove some dependencies.. 
        ! nah.. i can reuse this here in the set_length routine! 
        ! since the length equal the number of sites in the chain! 
        class(chain) :: this 
        integer, intent(in) :: length_x, length_y 
        integer, intent(in), optional :: length_z
        integer :: n_sites 
        character(*), parameter :: this_routine = "calc_nsites_chain" 

        if (max(length_x,length_y) < 1 .or. min(length_x, length_y) > 1 .or. & 
            min(length_x,length_y) < 0) then 
            n_sites = -1 
            call stop_all(this_routine, "something went wrong in length input!")

        else 
            n_sites = max(length_x, length_y)

        end if

    end function calc_nsites_chain

    function calc_nsites_cube(this, length_x, length_y, length_z) result(n_sites)
        class(cube) :: this
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        integer :: n_sites
        character(*), parameter :: this_routine = "calc_nsites_cube"
        
        if (max(length_x, length_y, length_z) < 2) then 
            call stop_all(this_routine, "too small cube lengths specified! (< 2)")
        end if

        n_sites = length_x * length_y * length_z

    end function calc_nsites_cube

    function calc_nsites_hexagonal(this, length_x, length_y, length_z) result(n_sites) 
        class(hexagonal) :: this 
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        integer :: n_sites 
        character(*), parameter :: this_routine = "calc_nsites_hexagonal" 

        ! the length_x of the hexagonal is defined as the number of unit cells.. 
        ! and there are 8 sites in my hexagonal unit cell.. 
        ASSERT(length_x > 0) 
        ASSERT(length_y > 0) 

        n_sites = length_x * length_y * 8 

    end function calc_nsites_hexagonal

    function calc_nsites_kagome(this, length_x, length_y, length_z) result(n_sites) 
        class(kagome) :: this 
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        integer :: n_sites 
        character(*), parameter :: this_routine = "calc_nsites_kagome" 

        ! the length_x and length_y of the kagome are defined as the number of unit cells.. 
        ! and there are 8 sites in my kagome unit cell.. 
        ASSERT(length_x > 0) 
        ASSERT(length_y > 0) 

        n_sites = length_x * length_y * 6 

    end function calc_nsites_kagome

    function calc_nsites_rect(this, length_x, length_y, length_z) result(n_sites) 
        class(rectangle) :: this 
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        integer :: n_sites
        character(*), parameter :: this_routine = "calc_nsites_rect"

        if (length_x < 2 .or. length_y < 2) then 
            print *, "length_x: ", length_x
            print *, "length_y: ", length_y
            call stop_all(this_routine, "length input wrong for type rectangle!")

        else 
            n_sites = length_x * length_y

        end if

    end function calc_nsites_rect

    function calc_nsites_tilted(this, length_x, length_y, length_z) result(n_sites)
        class(tilted) :: this 
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        integer :: n_sites 
        character(*), parameter :: this_routine = "calc_nsites_tilted" 

        if (length_x < 2 .or. length_y < 2) then 
            print *, "length_x: ", length_x
            print *, "length_y: ", length_y
            call stop_all(this_routine, "length input wrong for type tilted!")

        else 
            n_sites = 2 * length_x * length_y
        end if

    end function calc_nsites_tilted

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

    subroutine set_length_aim_star(this, length_x, length_y, length_z)
        class(aim_star) :: this 
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        character(*), parameter :: this_routine = "set_length_aim_star"

        ! actually the length of a start is not really defined.. 
        ! maybe i should rethink if i make this part of the 
        ! original lattice class then.. 
        call stop_all(this_routine, "length not defined for 'star' geometry!")
    end subroutine set_length_aim_star

    subroutine set_length_star(this, length_x, length_y, length_z)
        class(star) :: this 
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        character(*), parameter :: this_routine = "set_length_star"

        ! actually the length of a start is not really defined.. 
        ! maybe i should rethink if i make this part of the 
        ! original lattice class then.. 
        call stop_all(this_routine, "length not defined for 'star' geometry!")
    end subroutine set_length_star

    subroutine set_length_chain(this, length_x, length_y, length_z) 
        class(chain) :: this 
        integer, intent(in) :: length_x, length_y 
        integer, intent(in), optional :: length_z
        character(*), parameter :: this_routine = "set_length_chain"

        ! the input checkin is all done in the calc_nsites routine!
        this%length = this%calc_nsites(length_x, length_y)

    end subroutine set_length_chain

    subroutine set_length_cube(this, length_x, length_y, length_z)
        class(cube) :: this 
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "set_length_cube"
#endif

        ASSERT(present(length_z))
        ASSERT(length_x > 1)
        ASSERT(length_y > 1)
        ASSERT(length_z > 1) 

        this%length(1) = length_x
        this%length(2) = length_y
        this%length(3) = length_z

    end subroutine set_length_cube

    subroutine set_length_rect(this, length_x, length_y, length_z) 
        class(rectangle) :: this 
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z

        this%length(1) = length_x
        this%length(2) = length_y

    end subroutine set_length_rect

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

    subroutine set_periodic(this, t_periodic_x, t_periodic_y, t_periodic_z)
        class(lattice) :: this 
        logical, intent(in) :: t_periodic_x, t_periodic_y
        logical, intent(in), optional :: t_periodic_z

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

        this%t_periodic(1) = t_periodic_x
        this%t_periodic(2) = t_periodic_y
        if (present(t_periodic_z)) then 
            this%t_periodic(3) = t_periodic_z
        end if

    end subroutine set_periodic

    integer function get_nconnect_max(this) 
        class(lattice) :: this 
        
        get_nconnect_max = this%n_connect_max

    end function get_nconnect_max

    ! the star does not really have a concept of length.. 
    ! so always ouput STAR_LENGTH
    integer function get_length_star(this, dimen)
        class(star) :: this 
        integer, intent(in), optional:: dimen 

        get_length_star = STAR_LENGTH

    end function get_length_star

    integer function get_length_chain(this, dimen)
        class(chain) :: this 
        integer, intent(in), optional :: dimen

        if (present(dimen)) then 
            if (dimen > 1) then 
                get_length_chain = 1
            else 
                get_length_chain = this%length
            end if
        else
            get_length_chain = this%length
        end if

    end function get_length_chain

    integer function get_length_cube(this, dimen) 
        class(cube) :: this 
        integer, intent(in), optional :: dimen
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_length_cube"
#endif

        ASSERT(present(dimen))
        ASSERT(dimen > 0)
        ASSERT(dimen <= 3)

        get_length_cube = this%length(dimen) 

    end function get_length_cube

    integer function get_length_rect(this, dimen)
        class(rectangle) :: this 
        integer, intent(in), optional :: dimen 
        character(*), parameter :: this_routine = "get_length_rect"

        ASSERT(present(dimen)) 
        ASSERT(dimen > 0) 

        if (dimen > 2) then 
            get_length_rect = 1 
        else 
            get_length_rect = this%length(dimen)
        end if

    end function get_length_rect
        
    integer function get_length_aim_chain(this, dimen)
        class(aim_chain) :: this 
        integer, intent(in), optional :: dimen

        if (present(dimen) .and. dimen > 1) then 
            get_length_aim_chain = 1 
        else 
            get_length_aim_chain = this%length
        end if

    end function get_length_aim_chain

    subroutine set_length_aim_chain(this, length_x, length_y, length_z)
        class(aim_chain) :: this 
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z

        ! as a definition make the length, even for multiple impurity chains 
        ! as bath_sites + 1
        this%length = length_y + 1

    end subroutine set_length_aim_chain

    logical function is_periodic_star(this, dimen)
        ! this is always false.. the star geometry can't be periodic
        class(star) :: this 
        integer, intent(in), optional :: dimen

        is_periodic_star = .false.

    end function is_periodic_star

    logical function is_periodic_chain(this, dimen)
        class(chain) :: this 
        integer, intent(in), optional :: dimen

        ! we do not want to deal with two dimensional flags for chains or? 
        is_periodic_chain = this%t_periodic(1)
        ! the chain is only treated as periodic if both the flags are set 
        ! to be periodic!
!         is_periodic_chain = (this%is_periodic_x() .and. this%is_periodic_y())

    end function is_periodic_chain

    logical function is_periodic_cube(this, dimen) 
        class(cube) :: this 
        integer, intent(in), optional :: dimen 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "is_periodic_cube"
#endif 

        if (present(dimen)) then 
            ASSERT(dimen > 0) 
            ASSERT(dimen <= 3) 

            is_periodic_cube = this%t_periodic(dimen)

        else 
            is_periodic_cube = all(this%t_periodic)

        end if

    end function is_periodic_cube

    logical function is_periodic_rect(this, dimen) 
        class(rectangle) :: this
        integer, intent(in), optional :: dimen 
        character(*), parameter :: this_routine = "is_periodic_rect"

        ! depending if we want to have a certain periodic flag or 
        ! a check for full periodicity 
        if (present(dimen)) then
            ! we only consider the first 2 dimensions here
            ASSERT(dimen > 0)
            ASSERT(dimen <= 3)

            is_periodic_rect = this%t_periodic(dimen)

        else 
            is_periodic_rect = all(this%t_periodic)

        end if

    end function is_periodic_rect
    
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

    function determine_optimal_time_step(time_step_death) result(time_step)
        use SystemData, only: nel, bhub, uhub, t_new_real_space_hubbard, & 
                              t_tJ_model, t_heisenberg_model, exchange_j, & 
                              nOccAlpha, nOccBeta, t_k_space_hubbard, omega, & 
                              nbasis

        real(dp), optional, intent(out) :: time_step_death
        real(dp) :: time_step
        ! move this time-step determination to this routine for the real
        ! space hubbard to have it fully conained
        character(*), parameter :: this_routine = "determine_optimal_time_step"

        real(dp) :: p_elec, p_hole, mat_ele, max_diag

        ! determine the optimal hubbard time-step for an optimized 
        ! hubbard excitation generation 
        ! the first electron is chosen at random 
        if (t_k_space_hubbard) then 
            p_elec = 1.0_dp / real(nOccAlpha * nOccBeta, dp) 
        else 
            p_elec = 1.0_dp / real(nel, dp)
        end if

        ! and for a picked electron the possible neighbors are looked for 
        ! open holes so the lowest probability is determined by the 
        ! maximum numbers of connections 
        ! do this in general in the same way for all types of 
        ! lattice models. thats not exact, but a good enough estimate
        if (t_k_space_hubbard) then 
            p_hole = 1.0_dp / real(nbasis - nel, dp) 
        else
            p_hole = 1.0_dp / real(lat%get_nconnect_max(), dp) 
        end if

        if (t_new_real_space_hubbard) then 

            mat_ele = real(abs(bhub), dp)

        else if (t_tJ_model) then 

            ! for the t-J take the maximum of hopping or exchange
            mat_ele = real(max(abs(bhub),abs(exchange_j)),dp)

        else if (t_heisenberg_model) then 

            mat_ele = real(abs(exchange_j), dp)

        else if (t_k_space_hubbard) then 
            mat_ele = abs(real(uhub,dp)/real(omega,dp))

        end if

        ! so the time-step is 
        time_step = p_elec * p_hole / mat_ele

        if (present(time_step_death)) then 
            ! the maximum death contribution is max(H_ii) - shift 
            ! and the maximum diagonal matrix element we know it is 
            ! the maximum possible number of doubly occupied sites times U 
            ! but this might be a too hard limit, since especially for high U 
            ! such an amount of double occupations is probably never reached 
            ! and the shift must also be included in the calculation.. 
            if (t_new_real_space_hubbard) then 
                max_diag = real(abs(uhub) * min(nOccAlpha, nOccBeta), dp)
            else if (t_tJ_model) then 
                print *, "todo!" 
            else if (t_heisenberg_model) then 
                ! we have to find the maximum number of non-frustrated 
                ! nearest neighbor spins! 
                print *, "todo: "
            else if (t_k_space_hubbard) then 
                print *, "todo" 

            end if

            time_step_death = 1.0_dp / max_diag
        end if

    end function determine_optimal_time_step


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
