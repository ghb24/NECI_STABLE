#include "macros.h"

module lattice_mod
    ! this will be the module to mimick oop in fortran and in the
    ! lattice excitation generation implementation in neci.
    ! my plan is to create classes of lattice implementations which
    ! all come from a same base class "lattice"
    ! and maybe also do the same for sites.. where in the AIM then eg.
    ! the bath and impurity sites extend the base site class

    ! for now disable the OneEInts usage, due to circular dependencies over
    ! the sym_mod!
    use constants, only: dp, pi, EPS
    use SystemData, only: twisted_bc, nbasis, basisfn, t_trans_corr_2body, &
                          symmetry, brr, t_input_order, orbital_order, &
                          t_k_space_hubbard, t_trans_corr_hop, &
                          t_new_real_space_hubbard
    use input_parser_mod, only: ManagingFileReader_t, TokenIterator_t
    use fortran_strings, only: to_upper, to_lower, to_int, to_realdp
    use util_mod, only: stop_all

    implicit none
    private
    public :: lattice, lattice_deconstructor, aim, aim_deconstructor, sort_unique, &
              lat, determine_optimal_time_step, get_helement_lattice, inside_bz_2d, &
              on_line_2d, epsilon_kvec, init_dispersion_rel_cache, setup_lattice_symmetry

    integer, parameter :: NAME_LEN = 13
    integer, parameter :: sdim = 3

    HElement_t(dp), allocatable, public :: dispersion_rel_cached(:)

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

        ! i think i also want to store the k-vectors of the sites here..
        ! since i will need them if I totally want to switch to my new
        ! implementation and also if i want to deal with the new type of
        ! periodic boundary conditions.
        integer :: k_vec(3) = 0
        ! i also need the real-space coordinates for the hopping
        ! transcorrelation!
        integer :: r_vec(3) = 0

        ! also use one integer to differentiate between the k-vectors!
        ! this makes it easier to access arrays..
        integer :: k_sym = -1

        ! do i also need the inverse k-vec in here??
        integer :: k_inv(3) = 0
        integer :: sym_inv = -1

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
        procedure :: set_k_vec
        procedure :: set_r_vec
        ! i could also use finalization routines instead of manually
        ! deallocating everything..
        ! i need atleast gcc4.9.. which i am to lazy to update now..
        ! but will in the future!

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
        ! Lookup table for momentum -> site index conversion
        integer, allocatable :: lu_table(:, :, :)
        ! Lookup table for first BZ (contains all sums of up to three momenta)
        logical, allocatable :: bz_table(:, :, :)
        ! size of the lookup tables
        integer :: kmin(sdim) = 0
        integer :: kmax(sdim) = 0
        ! also store an indexing for the real-space vectors
        integer :: r_min(sdim) = 0
        integer :: r_max(sdim) = 0

        ! actually i want to have more flexibility: maybe periodic in x
        ! but not y..
        logical :: t_periodic_x = .true.
        logical :: t_periodic_y = .true.
        logical :: t_periodic(3) = .true.

        logical :: t_bipartite_order = .false.

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

        integer :: lat_vec(3, 3) = 0

        ! also store k_vec ..
        integer :: k_vec(3, 3) = 0

        ! initialize the possible applicable basis vectors in the
        ! mapping to the BZ
        integer, allocatable :: basis_vecs(:, :)

        ! i also need a matrix mapping from the k-vectors to the
        ! k-symbols to quickly access them!
        integer, allocatable, public :: k_to_sym(:, :, :)

        ! and vice versa a mapping from the symbol to the k-vector
        ! or i could just use the orbital index? does this work with
        ! neci though?
        ! just use a matrix here and take the rows
        integer, allocatable, public :: sym_to_k(:, :)

        ! and also store a multiplication table in the lattice class..
        ! to make it consistend and store everything necessary in here..
        ! this just make use of the symbols!
        integer, allocatable, public :: mult_table(:, :)

        ! and also use an inverse table, which also just uses the
        ! symbols!
        integer, allocatable, public :: inv_table(:)

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
        procedure, public :: get_site_index
        ! make the get neighbors function public on the lattice level
        procedure, public :: get_neighbors => get_neighbors_lattice
        procedure, public :: get_num_neighbors => get_num_neighbors_lattice
        procedure, public :: get_spinorb_neighbors => get_spinorb_neighbors_lat

        procedure, public :: is_k_space
        ! i definetly also want to have a print function!
        procedure, public :: print_lat
        procedure, public :: add_k_vec
        procedure :: add_k_vec_symbol
        procedure, public :: inv_k_vec
        procedure :: inv_k_vec_symbol
        procedure, public :: get_sym
        procedure, public :: subtract_k_vec
        procedure, public :: get_sym_from_k
        procedure, public :: set_sym

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
        procedure(calc_nsites_t), deferred :: calc_nsites
        procedure :: allocate_sites
        procedure(initialize_sites_t), deferred :: initialize_sites

        procedure :: deallocate_sites

        ! for the k-space implementations also implement a lattice
        ! dependent dispersion relation function
        procedure, public :: dispersion_rel => dispersion_rel_not_implemented
        procedure, public :: dispersion_rel_orb
        procedure, public :: dispersion_rel_spin_orb

        procedure, public :: dot_prod => dot_prod_not_implemented
        procedure, public :: get_k_vec
        procedure, public :: get_r_vec
        procedure, public :: round_sym

        procedure, public :: map_k_vec
        procedure :: inside_bz
        procedure :: inside_bz_explicit
        procedure :: apply_basis_vector

        procedure, public :: get_orb_from_k_vec
        ! and procedures to initialize the site index lookup table and the
        ! matrix element lookup table
        procedure :: initialize_lu_table
        procedure :: fill_bz_table
        procedure :: fill_lu_table
        procedure :: get_lu_table_size
        procedure :: deallocate_caches
        ! actually i should make i deferred: todo
        procedure :: init_basis_vecs
        procedure, public :: init_hop_cache_bounds

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

        procedure, public :: dispersion_rel => dispersion_rel_chain_k
        procedure :: init_basis_vecs => init_basis_vecs_chain

        procedure, public :: dot_prod => dot_prod_chain

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

        procedure :: init_basis_vecs_rect_base

        procedure, public :: get_length => get_length_rect
        procedure, public :: is_periodic => is_periodic_rect
        procedure, public :: dispersion_rel => dispersion_rel_rect

        procedure :: set_length => set_length_rect
        procedure :: calc_nsites => calc_nsites_rect
        procedure :: initialize_sites => init_sites_rect
        procedure :: init_basis_vecs => init_basis_vecs_rect
        procedure, public :: dot_prod => dot_prod_rect

    end type rectangle

    type, extends(rectangle) :: kagome
        private
    contains
        private

        procedure :: calc_nsites => calc_nsites_kagome
        procedure :: initialize_sites => init_sites_kagome

    end type kagome

    type, extends(rectangle) :: hexagonal
        ! i found a unit cell for the hexagonal lattice. but this is
        ! unfortunately 8 sites big alread.. anyway try it
        private

    contains
        private

        procedure :: calc_nsites => calc_nsites_hexagonal
        procedure :: initialize_sites => init_sites_hexagonal
    end type hexagonal

    type, extends(rectangle) :: triangular
        ! use the length quantitiy and periodicity of the rectangle
        private

    contains
        private

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
        procedure :: calc_nsites => calc_nsites_tilted
        procedure :: initialize_sites => init_sites_tilted
        procedure :: init_basis_vecs => init_basis_vecs_tilted
        procedure, public :: dot_prod => dot_prod_tilted
    end type tilted

    type, extends(rectangle) :: sujun
        private

    contains
        private
        procedure :: calc_nsites => calc_nsites_sujun
        procedure :: initialize_sites => init_sites_sujun

    end type sujun

    type, extends(rectangle) :: ext_input
        private

    contains
        private

        procedure :: calc_nsites => read_lattice_n_sites
        procedure :: initialize_sites => read_sites

    end type ext_input

    type, extends(rectangle) :: ole
        private

    contains
        private
        procedure, public :: dispersion_rel => dispersion_rel_ole
        procedure :: calc_nsites => calc_nsites_ole
        procedure :: initialize_sites => init_sites_ole
        procedure :: find_periodic_neighbors => find_periodic_neighbors_ole

        procedure :: inside_bz => inside_bz_ole

    end type ole

    ! can i just extend the chain class to make an impurity chain?
    ! argh this is annoying without multiple inheritance..
    type, extends(aim) :: aim_chain
        private

        integer :: length = -1

    contains
        private

        procedure, public :: get_length => get_length_aim_chain

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

    end type aim_star

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

    ! create the abstract interfaces for the deferred function in the base
    ! abstract type: lattice
    abstract interface

        pure function get_length_t(this, dimen) result(length)
            import :: lattice
            class(lattice), intent(in) :: this
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

    interface assignment(=)
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
            integer, intent(in) :: nI(nel), ic, ex(2, ic)
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

    interface epsilon_kvec
        module procedure epsilon_kvec_vector
        module procedure epsilon_kvec_symmetry
        module procedure epsilon_kvec_orbital
    end interface epsilon_kvec

contains

    subroutine setup_lattice_symmetry
        ! since i need it also in the real-space lattice for the
        ! hopping transcorrelation move the symmetry setup for the
        ! k-spae hubbard model into the lattice_mod
#ifdef DEBUG_
        character(*), parameter :: this_routine = "setup_lattice_symmetry"
#endif
        integer :: i, kmin(3), kmax(3), j, k_i(3), k, l, ind
        ASSERT(associated(lat))

        if (allocated(lat%k_to_sym)) deallocate(lat%k_to_sym)
        if (allocated(lat%sym_to_k)) deallocate(lat%sym_to_k)
        if (allocated(lat%mult_table)) deallocate(lat%mult_table)
        if (allocated(lat%inv_table)) deallocate(lat%inv_table)

        allocate(lat%sym_to_k(lat%get_nsites(), 3))
        allocate(lat%mult_table(lat%get_nsites(), lat%get_nsites()))
        allocate(lat%inv_table(lat%get_nsites()))

        ! i have to setup the symlabels first ofc..
        do i = 1, lat%get_nsites()
            ! and also just encode the symmetry labels as integers, instead of
            ! 2^(k-1), to be able to treat more than 64 orbitals (in the old
            ! implementation, an integer overflow happened in this case!)
            ind = get_spatial(brr(2 * i))

            call lat%set_sym(ind, i)

        end do

        kmin = 0
        kmax = 0
        do i = 1, lat%get_nsites()
            k_i = lat%get_k_vec(i)
            do j = 1, lat%get_ndim()
                if (k_i(j) < kmin(j)) kmin(j) = k_i(j)
                if (k_i(j) > kmax(j)) kmax(j) = k_i(j)
            end do
        end do

        allocate(lat%k_to_sym(kmin(1):kmax(1), kmin(2):kmax(2), kmin(3):kmax(3)))

        lat%k_to_sym = 0

        ! now find the inverses:
        do i = 1, lat%get_nsites()

            ! find the orbital of -k
            j = lat%get_orb_from_k_vec(-lat%get_k_vec(i))

            lat%inv_table(lat%get_sym(i)) = lat%get_sym(j)

            lat%sym_to_k(lat%get_sym(i), :) = lat%get_k_vec(i)

            k_i = lat%get_k_vec(i)

            lat%k_to_sym(k_i(1), k_i(2), k_i(3)) = lat%get_sym(i)

            ! and create the symmetry product of (i) with every other symmetry
            do k = 1, lat%get_nsites()
                ! i just have to add the momenta and map it to the first BZ
                l = lat%get_orb_from_k_vec(lat%get_k_vec(i) + lat%get_k_vec(k))

                lat%mult_table(lat%get_sym(i), lat%get_sym(k)) = lat%get_sym(l)

            end do
        end do

    end subroutine setup_lattice_symmetry

    HElement_t(dp) function epsilon_kvec_vector(k_vec)
        ! and actually this function has to be defined differently for
        ! different type of lattices! TODO!
        ! actually i could get rid of this function and directly call
        ! the dispersion relation of the lattice..
        integer, intent(in) :: k_vec(3)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "epsilon_kvec_vector"
#endif

        ASSERT(associated(lat))

        ! i could save the basic lattice vectors for the lattice or even
        ! store the dispersion relation for each lattice type and call it
        ! with a given k-vector?
        ! change that to only access the cached result of the dispersion
        ! relation

        ! it is necessary to call this function only with k_vectors
        ! within the first BZ!!
        epsilon_kvec_vector = dispersion_rel_cached(lat%get_sym_from_k(k_vec))

    end function epsilon_kvec_vector

    HElement_t(dp) function epsilon_kvec_symmetry(sym)
        ! access the stored dispersion relation values through the symmetry
        ! symbol associated with a k-vector
        type(symmetry), intent(in) :: sym

        epsilon_kvec_symmetry = dispersion_rel_cached(sym%s)

    end function epsilon_kvec_symmetry

    HElement_t(dp) function epsilon_kvec_orbital(orb)
        ! access the stored dispersion relation values through the spatial
        ! orbital (orb)
        integer, intent(in) :: orb

        epsilon_kvec_orbital = dispersion_rel_cached(lat%get_sym(orb))

    end function epsilon_kvec_orbital

    subroutine init_dispersion_rel_cache()
        ! to avoid excessive calls to the cos() function cache the
        ! dispersion relation of the lattice and make them accessible
        ! through the symmetry label associated with the k-vectors!
        character(*), parameter :: this_routine = "init_dispersion_rel_cache"
        integer :: i, sym_min, sym_max, sym

        ASSERT(associated(lat))

        if (allocated(dispersion_rel_cached)) deallocate(dispersion_rel_cached)

        sym_min = 0
        sym_max = 0
        do i = 1, lat%get_nsites()
            sym = lat%get_sym(i)

            if (sym < sym_min) sym_min = sym
            if (sym > sym_max) sym_max = sym
        end do

        allocate(dispersion_rel_cached(sym_min:sym_max), source = h_cast(0.0_dp))

        do i = 1, lat%get_nsites()
            dispersion_rel_cached(lat%get_sym(i)) = &
                lat%dispersion_rel_orb(i)
        end do

    end subroutine init_dispersion_rel_cache

    subroutine set_sym(this, orb, sym)
        class(lattice) :: this
        integer, intent(in) :: orb, sym

        this%sites(orb)%k_sym = sym

    end subroutine set_sym

    pure function add_k_vec(this, k_1, k_2) result(k_out)
        class(lattice), intent(in) :: this
        integer, intent(in) :: k_1(3), k_2(3)
        integer :: k_out(3)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "add_k_vec"
#endif

        ASSERT(allocated(this%mult_table))
        ASSERT(allocated(this%k_to_sym))

        k_out = this%sym_to_k(this%mult_table( &
                              this%k_to_sym(k_1(1), k_1(2), k_1(3)), &
                              this%k_to_sym(k_2(1), k_2(2), k_2(3))), :)

    end function add_k_vec

    function add_k_vec_symbol(this, sym_1, sym_2) result(sym_out)
        class(lattice) :: this
        integer, intent(in) :: sym_1, sym_2
        integer :: sym_out
#ifdef DEBUG_
        character(*), parameter :: this_routine = "add_k_vec_symbol"
#endif

        ASSERT(allocated(this%mult_table))

        sym_out = this%mult_table(sym_1, sym_2)

    end function add_k_vec_symbol

    function inv_k_vec(this, k) result(k_inv)
        class(lattice) :: this
        integer, intent(in) :: k(3)
        integer :: k_inv(3)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "inv_k_vec"
#endif

        ASSERT(allocated(this%sym_to_k))
        ASSERT(allocated(this%inv_table))
        ASSERT(allocated(this%k_to_sym))

        k_inv = this%sym_to_k(this%inv_table(this%k_to_sym(k(1), k(2), k(3))), :)

    end function inv_k_vec

    pure function get_sym(this, orb) result(sym)
        ! gives the symmetry label associated with the k-vector of
        ! spatial orbital (orb)
        class(lattice), intent(in) ::  this
        integer, intent(in) :: orb
        integer :: sym
        unused_var(this)

        sym = this%sites(orb)%k_sym

    end function get_sym

    pure function get_sym_from_k(this, k) result(sym)
        ! the routine to get the symmetry label associated with the k-vector k
        class(lattice), intent(in) :: this
        integer, intent(in) :: k(3)
        integer :: sym

        sym = this%k_to_sym(k(1), k(2), k(3))

    end function get_sym_from_k

    ! do i also need a inv_k_vec_symbol function?
    function inv_k_vec_symbol(this, sym) result(inv_sym)
        class(lattice) :: this
        integer, intent(in) :: sym
        integer :: inv_sym
#ifdef DEBUG_
        character(*), parameter :: this_routine = "inv_k_vec_symbol"
#endif

        ASSERT(allocated(this%inv_table))

        inv_sym = this%inv_table(sym)

    end function inv_k_vec_symbol

    function subtract_k_vec(this, k_1, k_2) result(k_out)
        class(lattice) :: this
        integer, intent(in) :: k_1(3), k_2(3)
        integer :: k_out(3)
        unused_var(this)

        k_out = this%add_k_vec(k_1, this%inv_k_vec(k_2))

    end function subtract_k_vec

    function get_orb_from_k_vec(this, k_in, spin) result(orb)
        class(lattice) :: this
        integer, intent(in) :: k_in(3)
        integer, intent(in), optional :: spin
        integer :: orb
#ifdef DEBUG_
        character(*), parameter :: this_routine = "get_orb_from_k_vec"
#endif
        integer :: k_vec(3), i

        ! checking if it is in the first is not necessary anymore
        ! as the lookup table captures more than just the BZ

        ! the naive way would be to loop over all sites and check if the
        ! k-vector fits..
        ! but that would be too effortive, so we use the lookup table

        i = this%lu_table(k_in(1), k_in(2), k_in(3))

        ! and, if required, include the spin in the index
        if (present(spin)) then
            ASSERT(spin == 1 .or. spin == 2)
            if (spin == 1) then
                orb = 2 * i - 1
            else if (spin == 2) then
                orb = 2 * i
            end if
        else
            orb = i
        end if

    end function get_orb_from_k_vec

    elemental function round_sym(this, sym_in) result(sym_out)
        ! routine to map k-vectors outside first BZ back inside
        class(lattice), intent(in) :: this
        type(basisfn), intent(in) :: sym_in
        type(basisfn) :: sym_out

        integer :: k_vec(3)

        ! write a lattice specific routine, which checks if the k-vector is
        ! inside the first BZ of this lattice (which has to be defined by
        ! the input and implemented by me!)
        if (this%inside_bz(sym_in%k)) then
            ! then i have to do nothing
            sym_out = sym_in
        else
            ! otherwise map the k-vector back..
            k_vec = this%map_k_vec(sym_in%k)
            sym_out = sym_in
            sym_out%k = k_vec
        end if

    end function round_sym

    pure function map_k_vec(this, k_in) result(k_out)
        class(lattice), intent(in) :: this
        integer, intent(in) :: k_in(3)
        integer :: k_out(3)

        integer :: i

        k_out = k_in

        if (this%inside_bz(k_in)) then
            k_out = k_in

        else
            ! here i have to do something..
            ! should i store this matrix to setup the lattice within the
            ! lattice class? so i can reuse it here..
            ! or i apply the primitive vectors to the k_vec and check if
            ! a resulting vector lies within the first BZ..
            i = 1
            k_out = k_in
            do while (.not. this%inside_bz(k_out))
                ! apply all possible basis vectors of the lattice
                k_out = this%apply_basis_vector(k_in, i)
                i = i + 1
            end do
        end if

    end function map_k_vec

    logical pure function inside_bz(this, k_vec)
        class(lattice), intent(in) :: this
        integer, intent(in) :: k_vec(3)
        character(*), parameter :: this_routine = "inside_bz"

        ! this function should also be deferred!

        ! i think with Kais new BZ implementation we can write this function
        ! generally.

        ! I think this should be the approach for most lattices
        ! do a check if we have the bz-ishness of this vector stored
        if (all(k_vec <= this%kmax) .and. all(k_vec >= this%kmin)) then
            inside_bz = this%bz_table(k_vec(1), k_vec(2), k_vec(3))
        else
            ! if not, do the explicit check
            inside_bz = this%inside_bz_explicit(k_vec)
        end if

    end function inside_bz

    logical pure function inside_bz_explicit(this, k_vec)
        class(lattice), intent(in) :: this
        integer, intent(in) :: k_vec(sdim)

        integer :: i

        ! oles lattice is defined by four corner points and two lines
        inside_bz_explicit = .false.

        do i = 1, this%get_nsites()
            if (all(k_vec == this%get_k_vec(i))) inside_bz_explicit = .true.
        end do

    end function inside_bz_explicit

    logical pure function inside_bz_ole(this, k_vec)
        class(ole), intent(in) :: this
        integer, intent(in) :: k_vec(sdim)

        ! I think this should be the approach for most lattices
        ! do a check if we have the bz-ishness of this vector stored
        if (all(k_vec <= this%kmax) .and. all(k_vec >= this%kmin)) then
            inside_bz_ole = this%bz_table(k_vec(1), k_vec(2), k_vec(3))
        else
            ! if not, do the explicit check
            inside_bz_ole = this%inside_bz_explicit(k_vec)
        end if
    end function inside_bz_ole

    logical function inside_bz_chain(this, k_vec)
        class(chain) :: this
        integer, intent(in) :: k_vec(3)

        ! the chain goes as -Lx/2+1, +2, .. Lx/2
        inside_bz_chain = .false.

        if (k_vec(1) >= (-(this%length + 1) / 2 + 1) .and. &
            k_vec(1) <= this%length / 2) inside_bz_chain = .true.

    end function inside_bz_chain

    pure function apply_basis_vector_general(this, k_in) result(k_out)
        ! todo: i think i can write a general routine to apply the basis
        ! vectors.. so i do not have to write a specific one for all the
        ! different sorts of lattices..
        class(lattice), intent(in) :: this
        integer, intent(in) :: k_in(sdim)
        integer, allocatable :: k_out(:, :)

        character(*), parameter :: this_routine = "apply_basis_vector_general"

        unused_var(this)
        unused_var(k_in)

        call stop_all(this_routine, "not yet implemented!")

        k_out = 0

    end function apply_basis_vector_general

    pure function apply_basis_vector_cube(this, k_in, ind) result(k_out)
        class(cube), intent(in) :: this
        integer, intent(in) :: k_in(3)
        integer, intent(in), optional :: ind
        integer :: k_out(3)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "apply_basis_vector_cube"
#endif

        call stop_all("apply_basis_vector_cube", "not yet implemented!")
#ifdef WARNING_WORKAROUND_
        k_out = 0
#endif
        unused_var(this)
        unused_var(k_in)
        if (present(ind)) then
            unused_var(ind)
        end if

    end function apply_basis_vector_cube

    pure function apply_basis_vector_rect(this, k_in, ind) result(k_out)
        class(rectangle), intent(in) :: this
        integer, intent(in) :: k_in(3)
        integer, intent(in), optional :: ind
        integer :: k_out(3)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "apply_basis_vector_rect"
#endif

        integer :: basis_vec(8, 3), r1(3), r2(3)

        ASSERT(ind >= 0)
        ASSERT(ind <= 8)

        ! with negative signs we have in total 8 possibilities of
        ! vectors to apply:

        r1 = this%k_vec(:, 1)
        r2 = this%k_vec(:, 2)

        basis_vec(1, :) = r1
        basis_vec(2, :) = -r1
        basis_vec(3, :) = r2
        basis_vec(4, :) = -r2
        basis_vec(5, :) = r1 + r2
        basis_vec(6, :) = -(r1 + r2)
        basis_vec(7, :) = r1 - r2
        basis_vec(8, :) = -(r1 - r2)

        k_out = k_in + basis_vec(ind, :)

    end function apply_basis_vector_rect

    function apply_basis_vector_chain(this, k_in, ind) result(k_out)
        ! i think i have to make ind non-optional, as it is actually
        ! always needed.. what else should i do here?
        class(chain) :: this
        integer, intent(in) :: k_in(3)
        integer, intent(in), optional :: ind
        integer :: k_out(3)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "apply_basis_vector_chain"
#endif
        integer :: basis_vec(4, 3)

#ifdef DEBUG_
        if (t_trans_corr_2body) then
            ASSERT(ind > 0 .and. ind <= 4)
        else
            ASSERT(ind == 1 .or. ind == 2)
        end if
#endif

        k_out = k_in + this%basis_vecs(ind, :)

    end function apply_basis_vector_chain

    subroutine init_basis_vecs(this)
        class(lattice) :: this
        character(*), parameter :: this_routine = "init_basis_vecs"

        unused_var(this)
        ! K.G. 25.11.2019: Why is this a runtime check? Should be done compile-time
        call stop_all(this_routine, "this routine should always be deferred!")

    end subroutine init_basis_vecs

    subroutine init_basis_vecs_chain(this)
        class(chain) :: this

        integer :: i, j, k

        if (allocated(this%basis_vecs)) deallocate(this%basis_vecs)

        if (t_trans_corr_2body) then
            j = 2
            allocate(this%basis_vecs(5, 3))
        else
            j = 1
            allocate(this%basis_vecs(3, 3))
        end if

        this%basis_vecs = 0

        k = 0
        do i = -j, j
            k = k + 1
            this%basis_vecs(k, 1) = i * this%get_length(1)
        end do

    end subroutine init_basis_vecs_chain

    subroutine init_basis_vecs_rect(this)
        class(rectangle) :: this

        integer :: l

        if (t_trans_corr_2body) then
            l = 2
        else
            l = 1
        end if

        call this%init_basis_vecs_rect_base(l)
    end subroutine init_basis_vecs_rect

    subroutine init_basis_vecs_tilted(this)
        class(tilted) :: this

        ! Tilted lattices require more basis vectors stored (up to triple application of basis vector)
        call this%init_basis_vecs_rect_base(4)
    end subroutine init_basis_vecs_tilted

    !> Base function for setting up a the basis vector array for rectangular lattices (extracted from the previous init_basis_vecs_rect)
    !> @param[in] l  Maximal number of unit vectors to be combined into a basis vector
    subroutine init_basis_vecs_rect_base(this,l)
        class(rectangle), intent(inout) :: this
        integer, intent(in) :: l

        integer :: i,j,k

        if (allocated(this%basis_vecs)) deallocate(this%basis_vecs)
        allocate(this%basis_vecs((2*l+1)**2,3))
        this%basis_vecs = 0
        k = 0
        do i = -l, l
            do j = -l, l
                k = k + 1
                this%basis_vecs(k, :) = i * this%k_vec(:, 1) + j * this%k_vec(:, 2)
            end do
        end do

    end subroutine init_basis_vecs_rect_base

    pure function apply_basis_vector(this, k_in, ind) result(k_out)
        ! i have to specifically write this for every lattice type..
        class(lattice), intent(in) :: this
        integer, intent(in) :: k_in(3)
        integer, intent(in), optional :: ind
        integer :: k_out(3)
        character(*), parameter :: this_routine = "apply_basis_vector"

        ! i should only make this an abstract interface, since this function
        ! must be deffered!
        ASSERT(ind >= 0)
        ASSERT(ind <= size(this%basis_vecs, 1))

        k_out = k_in + this%basis_vecs(ind, :)

    end function apply_basis_vector

    function apply_basis_vector_ole(this, k_in, ind) result(k_out)
        class(ole) :: this
        integer, intent(in) :: k_in(3)
        integer, intent(in), optional :: ind
        integer :: k_out(3)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "apply_basis_vector_ole"
#endif

        integer :: basis_vec(8, 3), r1(3), r2(3)

        ASSERT(ind >= 0)
        ASSERT(ind <= 8)

        ! with negative signs we have in total 8 possibilities of
        ! vectors to apply:

        r1 = this%k_vec(:, 1)
        r2 = this%k_vec(:, 2)

        basis_vec(1, :) = r1
        basis_vec(2, :) = -r1
        basis_vec(3, :) = r2
        basis_vec(4, :) = -r2
        basis_vec(5, :) = r1 + r2
        basis_vec(6, :) = -(r1 + r2)
        basis_vec(7, :) = r1 - r2
        basis_vec(8, :) = -(r1 - r2)

        k_out = k_in + basis_vec(ind, :)

    end function apply_basis_vector_ole

    integer pure function get_length_aim_star(this, dimen)
        class(aim_star), intent(in) :: this
        integer, intent(in), optional :: dimen
        unused_var(this)
        if (present(dimen)) then
            unused_var(dimen)
        end if

        unused_var(this)
        unused_var(dimen)

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
        this%sites(1) = site(1, 1, [2], site_type='impurity')

        ! the bath sites are connected within each other, but not periodic!
        do i = 1, this%get_n_bath() - 1
            this%sites(i + 1) = site(i + 1, 2, [i, i + 2], site_type='bath')
        end do

        ! and the last bath site only has one neighbor
        this%sites(this%get_nsites()) = site(this%get_nsites(), 1, &
                                             [this%get_nsites() - 1], site_type='bath')

    end subroutine init_sites_aim_chain

    function calc_nsites_aim(this, length_x, length_y, length_z) result(n_sites)
        class(aim) :: this
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        integer :: n_sites
        character(*), parameter :: this_routine = "calc_nsites_aim"
        unused_var(this)
        if (present(length_z)) then
            unused_var(length_z)
        end if

        unused_var(this)
        unused_var(length_z)

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

    end function is_bath_site

    logical function is_impurity_site(this, ind)
        class(aim) :: this
        integer, intent(in) :: ind
        character(*), parameter :: this_routine = "is_impurity_site"

        ASSERT(ind > 0)
        ASSERT(ind <= this%get_nsites())

        is_impurity_site = this%sites(ind)%is_impurity()

    end function is_impurity_site

    function get_bath(this) result(bath_sites)
        class(aim) :: this
        integer :: bath_sites(this%n_bath)
        character(*), parameter :: this_routine = "get_bath"

        integer :: i, j

        bath_sites = this%bath_sites

    end function get_bath

    function get_impurities(this) result(imp_sites)
        class(aim) :: this
        integer :: imp_sites(this%n_imps)
        character(*), parameter :: this_routine = "get_impurities"

        integer :: i, j

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
        unused_var(this)
        if (present(dimen)) then
            unused_var(dimen)
        end if

        unused_var(this)
        unused_var(dimen)
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
#ifdef DEBUG_
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

        unused_var(this)

    end subroutine init_sites_lattice

    subroutine site_assign(lhs, rhs)
        type(site), intent(out) :: lhs
        type(site), intent(in) :: rhs

        ! argh.. on every change in the site constructor i have to make
        ! the change here too.. annoying..
        ! but make i OOP style!
        call lhs%set_index(rhs%get_index())
        call lhs%set_num_neighbors(rhs%get_num_neighbors())
        call lhs%allocate_neighbors(rhs%get_num_neighbors())
        call lhs%set_neighbors(rhs%get_neighbors())
        call lhs%set_impurity(rhs%is_impurity())
        call lhs%set_bath(rhs%is_bath())
        call lhs%set_k_vec(rhs%k_vec)
        call lhs%set_r_vec(rhs%r_vec)

        ! can i work with an allocatable statement here?

    end subroutine site_assign

    subroutine lattice_assign(lhs, rhs)
        ! specific lattice assigner to not have to use pointers..
        class(lattice), intent(out) :: lhs
        class(lattice), intent(in), pointer :: rhs

        character(*), parameter :: this_routine = "lattice_assign"

        unused_var(rhs)

        ! here i have to copy all the specific values!
        ! this is annoying but make the code more readable, and i do not
        ! have to use pointers so much..
        ! todo although.. since i cannot use class(lattice) in the main
        ! program without allocatable or pointer attribute anyway..
        ! i think there is no need of an assignment overload!
        call stop_all(this_routine, "Lattice assignment operator deleted")
        unused_var(rhs)

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
        this%sites(1) = site(1, this%get_n_bath(), [(i, i=2, this%get_nsites())], &
                             site_type='impurity')

        ! and all the bath sites are just connected to the impurity
        do i = 2, this%get_nsites()
            this%sites(i) = site(i, 1, [1], site_type='bath')
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
            this%sites(1) = site(ind=1, n_neighbors=0, neighbors=[integer ::])

        else
            ! first to the special pivot site in the middle of the star
            this%sites(1) = site(ind=1, n_neighbors=this%get_nconnect_max(), &
                                 neighbors=[(i, i=2, this%get_nsites())])

            ! and all the others are just connected to the pivot
            do i = 2, this%get_nsites()
                this%sites(i) = site(ind=i, n_neighbors=1, neighbors=[1])
            end do
        end if

    end subroutine init_sites_star

    subroutine init_sites_chain(this)
        class(chain) :: this

        integer :: i, vec(3), j, n
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

        if (this%t_bipartite_order) then
            if (t_input_order) then
                n = this%get_nsites()
                if (this%is_periodic()) then
                    vec = [-(this%length + 1) / 2 + orbital_order(1), 0, 0]
                    this%sites(orbital_order(1)) = site(orbital_order(1), 2, &
                        [orbital_order(n), orbital_order(2)], vec, vec)

                    vec = [this%length / 2, 0, 0]
                    this%sites(orbital_order(n)) = site(orbital_order(n), 2, &
                        [orbital_order(n-1), orbital_order(1)], vec, vec)


                else
                    vec = [-(this%length + 1) / 2 + orbital_order(1), 0, 0]
                    this%sites(orbital_order(1)) = site(orbital_order(1), 1, &
                        [orbital_order(2)], vec, vec)

                    vec = [this%length / 2, 0, 0]
                    this%sites(orbital_order(n)) = site(orbital_order(n), 1, &
                        [orbital_order(n-1)], vec, vec)

                end if

                ! and do the rest inbetween which is always the same
                do i = 2, this%get_nsites() - 1

                    vec = [-(this%length + 1) / 2 + orbital_order(i), 0, 0]
                    ! if periodic and first or last: already dealt with above
                    this%sites(orbital_order(i)) = site(orbital_order(i), &
                        N_CONNECT_MAX_CHAIN, [orbital_order(i-1),orbital_order(i + 1)], vec, vec)

                end do


            else
                if (this%is_periodic()) then

                    ! use more concise site contructors!
                    ! also encode k- and real-space vectors.. i have to get this right!
                    vec = [-(this%length + 1) / 2 + 1, 0, 0]

                    this%sites(1) = site(1, 2, &
                        [this%get_nsites(), this%get_nsites()/2 + 1], vec, vec)

                    vec = [this%length / 2, 0, 0]
                    this%sites(this%get_nsites()) = site(this%get_nsites(), 2, &
                                             [this%get_nsites()/2, 1], vec, vec)

                else
                    ! open boundary conditions:
                    ! first site:
                    this%sites(1) = site(1, 1, [this%get_nsites()/2 + 1], &
                                         [-(this%length + 1) / 2 + 1, 0, 0])

                    ! last site:
                    this%sites(this%get_nsites()) = site(this%get_nsites(), 1, &
                                         [this%get_nsites()/2], [this%length / 2, 0, 0])

                end if

                ! and do the rest inbetween which is always the same
                do i = 2, this%get_nsites()/2

                    vec = [-(this%length + 1) / 2 + 2 * i - 1, 0, 0]
                    ! if periodic and first or last: already dealt with above
                    this%sites(i) = site(i, N_CONNECT_MAX_CHAIN, &
                        [this%get_nsites()/2 + i - 1, this%get_nsites()/2 + i], vec, vec)

                end do

                j = 1
                do i = this%get_nsites()/2 + 1, this%get_nsites() - 1
                    vec = [-(this%length + 1) / 2 + 2*j , 0, 0]
                    ! if periodic and first or last: already dealt with above
                    this%sites(i) = site(i, N_CONNECT_MAX_CHAIN, &
                        [j, j + 1], vec, vec)

                    j = j + 1
                end do
            end if
        else
            if (this%get_nsites() == 1) then
                this%sites(1) = site(ind=1, n_neighbors=0, neighbors=[integer ::], &
                                     k_vec=[0, 0, 0], r_vec=[0, 0, 0])
                return
            end if

            if (this%is_periodic()) then
                ! use more concise site contructors!
                ! also encode k- and real-space vectors.. i have to get this right!
                vec = [-(this%length + 1) / 2 + 1, 0, 0]

                this%sites(1) = site(1, 2, [this%get_nsites(), 2], vec, vec)

                vec = [this%length / 2, 0, 0]
                this%sites(this%get_nsites()) = site(this%get_nsites(), 2, &
                                         [this%get_nsites() - 1, 1], vec, vec)

            else
                ! open boundary conditions:
                ! first site:
                this%sites(1) = site(1, 1, [2], &
                                     [-(this%length + 1) / 2 + 1, 0, 0])

                ! last site:
                this%sites(this%get_nsites()) = site(this%get_nsites(), 1, &
                             [this%get_nsites() - 1], [this%length / 2, 0, 0])

            end if

            ! and do the rest inbetween which is always the same
            do i = 2, this%get_nsites() - 1

                vec = [-(this%length + 1) / 2 + i, 0, 0]
                ! if periodic and first or last: already dealt with above
                this%sites(i) = site(i, N_CONNECT_MAX_CHAIN, [i - 1, i + 1], vec, vec)

            end do
        end if

    end subroutine init_sites_chain

    subroutine init_sites_cube(this)
        class(cube) :: this
        character(*), parameter :: this_routine = "init_sites_cube"
        integer :: temp_array(this%length(1), this%length(2), this%length(3))
        integer :: i, x, y, z, temp_neigh(6)
        integer :: up(this%length(1), this%length(2), this%length(3))
        integer :: down(this%length(1), this%length(2), this%length(3))
        integer :: left(this%length(1), this%length(2), this%length(3))
        integer :: right(this%length(1), this%length(2), this%length(3))
        integer :: in(this%length(1), this%length(2), this%length(3))
        integer :: out(this%length(1), this%length(2), this%length(3))
        integer, allocatable :: neigh(:)

        ! minumum cube size:
        ASSERT(this%get_nsites() >= 8)

        ! encode the lattice with the fortran intrinsic ordering of matrices
        temp_array = reshape([(i, i=1, this%get_nsites())], &
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
                x = mod(i - 1, this%length(1)) + 1
                z = (i - 1) / (this%length(1) * this%length(2)) + 1
                y = mod((i - 1) / this%length(1), this%length(2)) + 1

                temp_neigh = [up(x, y, z), down(x, y, z), left(x, y, z), right(x, y, z), &
                              in(x, y, z), out(x, y, z)]

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
        integer, allocatable :: order(:)
        integer :: i
        character(*), parameter :: this_routine = "init_sites_kagome"

        allocate(order(this%get_nsites()), source = 0)

        if (t_input_order .and. this%t_bipartite_order) then
            order = orbital_order
        else
            order = [(i, i = 1, this%get_nsites())]
        end if


        ! and i think i will for the 6-site, 1x2 and 2x1 12 site and 2x2 24-site
        ! hard encode the neighbors and stuff because this seems to be a pain
        ! in the ass!

        if (this%is_periodic()) then
            if (this%length(1) == 1 .and. this%length(2) == 1) then
                ! the smallest cluster
                this%sites(order(1)) = site(order(1), 3, order([2, 4, 6]))
                this%sites(order(2)) = site(order(2), 4, order([1, 3, 4, 5]))
                this%sites(order(3)) = site(order(3), 3, order([2, 5, 6]))
                this%sites(order(4)) = site(order(4), 3, order([1, 2, 6]))
                this%sites(order(5)) = site(order(5), 3, order([2, 3, 6]))
                this%sites(order(6)) = site(order(6), 4, order([1, 3, 4, 5]))

            else if (this%length(1) == 1 .and. this%length(2) == 2) then
                ! the 1x2 cluster with 12-sites:
                this%sites(1) = site(1, 4, [2, 4, 10, 12])
                this%sites(2) = site(2, 4, [1, 3, 4, 5])
                this%sites(3) = site(3, 4, [2, 5, 11, 12])
                this%sites(4) = site(4, 4, [1, 2, 6, 7])
                this%sites(5) = site(5, 4, [2, 3, 6, 9])
                this%sites(6) = site(6, 4, [4, 5, 7, 9])
                this%sites(7) = site(7, 4, [4, 6, 8, 10])
                this%sites(8) = site(8, 4, [7, 9, 10, 11])
                this%sites(9) = site(9, 4, [5, 6, 8, 11])
                this%sites(10) = site(10, 4, [1, 7, 8, 12])
                this%sites(11) = site(11, 4, [3, 8, 9, 12])
                this%sites(12) = site(12, 4, [1, 3, 10, 11])

            else if (this%length(1) == 2 .and. this%length(2) == 1) then
                ! the 2x1 12-site cluster:
                this%sites(1) = site(1, 3, [2, 4, 12])
                this%sites(2) = site(2, 4, [1, 3, 4, 5])
                this%sites(3) = site(3, 3, [2, 5, 6])
                this%sites(4) = site(4, 3, [1, 2, 12])
                this%sites(5) = site(5, 3, [2, 3, 6])
                this%sites(6) = site(6, 4, [3, 5, 7, 10])
                this%sites(7) = site(7, 3, [6, 8, 10])
                this%sites(8) = site(8, 4, [7, 9, 10, 11])
                this%sites(9) = site(9, 3, [8, 11, 12])
                this%sites(10) = site(10, 3, [6, 7, 8])
                this%sites(11) = site(11, 3, [8, 9, 12])
                this%sites(12) = site(12, 4, [1, 4, 9, 11])

            else if (this%length(1) == 2 .and. this%length(2) == 2) then
                ! the 2x2 24-site cluster.. puh.. thats a lot to do..
                this%sites(1) = site(1, 4,   [2, 4, 10, 24])
                this%sites(2) = site(2, 4,   [1, 3,  4,  5])
                this%sites(3) = site(3, 4,   [2, 5, 11, 12])
                this%sites(4) = site(4, 4,   [1, 2,  7, 18])
                this%sites(5) = site(5, 4,   [2, 3,  6,  9])
                this%sites(6) = site(6, 4,   [5, 9, 16, 19])
                this%sites(7) = site(7, 4,   [4, 8, 10, 18])
                this%sites(8) = site(8, 4,   [7, 9, 10, 11])
                this%sites(9) = site(9, 4,   [5, 6,  8, 11])
                this%sites(10) = site(10, 4, [1,  7, 8, 24])
                this%sites(11) = site(11, 4, [3,  8, 9, 12])
                this%sites(12) = site(12, 4, [3, 11,13, 22])
                this%sites(13) = site(13, 4, [12,14,16, 22])
                this%sites(14) = site(14, 4, [13,15,16, 17])
                this%sites(15) = site(15, 4, [14,17,23, 24])
                this%sites(16) = site(16, 4, [6, 13,14, 19])
                this%sites(17) = site(17, 4, [14,15,18, 21])
                this%sites(18) = site(18, 4, [4, 7, 17, 21])
                this%sites(19) = site(19, 4, [6, 16,20, 22])
                this%sites(20) = site(20, 4, [19,21,22, 23])
                this%sites(21) = site(21, 4, [17,18,20, 23])
                this%sites(22) = site(22, 4, [12,13,19, 20])
                this%sites(23) = site(23, 4, [15,20,21, 24])
                this%sites(24) = site(24, 4, [1, 10,15, 23])

            else if (this%length(1) == 3 .and. this%length(2) == 2) then
                ! the 3x2x6 36-site cluster
                this%sites( 1) = site( 1, 4, [ 2, 4,16,36])
                this%sites( 2) = site( 2, 4, [ 1, 3, 4, 5])
                this%sites( 3) = site( 3, 4, [ 2, 5,17,18])
                this%sites( 4) = site( 4, 4, [ 1, 2, 7,24])
                this%sites( 5) = site( 5, 4, [ 2, 3, 6, 9])
                this%sites( 6) = site( 6, 4, [ 5, 9,22,25])
                this%sites( 7) = site( 7, 4, [ 4, 8,10,24])
                this%sites( 8) = site( 8, 4, [ 7, 9,10,11])
                this%sites( 9) = site( 9, 4, [ 5, 6, 8,11])
                this%sites(10) = site(10, 4, [ 7, 8,13,30])
                this%sites(11) = site(11, 4, [ 8, 9,12,15])
                this%sites(12) = site(12, 4, [11,15,28,31])
                this%sites(13) = site(13, 4, [10,14,16,36])
                this%sites(14) = site(14, 4, [13,15,16,17])
                this%sites(15) = site(15, 4, [11,12,14,17])
                this%sites(16) = site(16, 4, [ 1,13,14,36])
                this%sites(17) = site(17, 4, [ 3,14,15,18])
                this%sites(18) = site(18, 4, [ 3,17,19,34])
                this%sites(19) = site(19, 4, [18,20,22,34])
                this%sites(20) = site(20, 4, [19,21,22,23])
                this%sites(21) = site(21, 4, [20,23,35,36])
                this%sites(22) = site(22, 4, [ 6,19,20,25])
                this%sites(23) = site(23, 4, [20,21,24,27])
                this%sites(24) = site(24, 4, [ 4, 7,23,27])
                this%sites(25) = site(25, 4, [ 6,22,26,28])
                this%sites(26) = site(26, 4, [25,27,28,29])
                this%sites(27) = site(27, 4, [23,24,26,29])
                this%sites(28) = site(28, 4, [12,25,26,31])
                this%sites(29) = site(29, 4, [26,27,30,33])
                this%sites(30) = site(30, 4, [10,13,29,33])
                this%sites(31) = site(31, 4, [12,28,32,34])
                this%sites(32) = site(32, 4, [31,33,34,35])
                this%sites(33) = site(33, 4, [29,30,32,35])
                this%sites(34) = site(34, 4, [18,19,31,32])
                this%sites(35) = site(35, 4, [21,32,33,36])
                this%sites(36) = site(36, 4, [ 1,16,21,35])

            else if (this%length(1) == 4 .and. this%length(2) == 2) then
                ! the 4x2x6 48-site cluster
                this%sites( 1) = site( 1, 4, [ 2, 4,22,48])
                this%sites( 2) = site( 2, 4, [ 1, 3, 4, 5])
                this%sites( 3) = site( 3, 4, [ 2, 5,23,24])
                this%sites( 4) = site( 4, 4, [ 1, 2, 7,30])
                this%sites( 5) = site( 5, 4, [ 2, 3, 6, 9])
                this%sites( 6) = site( 6, 4, [ 5, 9,28,31])
                this%sites( 7) = site( 7, 4, [ 4, 8,10,30])
                this%sites( 8) = site( 8, 4, [ 7, 9,10,11])
                this%sites( 9) = site( 9, 4, [ 5, 6, 8,11])
                this%sites(10) = site(10, 4, [ 7, 8,13,36])
                this%sites(11) = site(11, 4, [ 8, 9,12,15])
                this%sites(12) = site(12, 4, [11,15,34,37])
                this%sites(13) = site(13, 4, [10,14,16,36])
                this%sites(14) = site(14, 4, [13,15,16,17])
                this%sites(15) = site(15, 4, [11,12,14,17])
                this%sites(16) = site(16, 4, [13,14,19,42])
                this%sites(17) = site(17, 4, [14,15,18,21])
                this%sites(18) = site(18, 4, [17,21,40,43])
                this%sites(19) = site(19, 4, [16,20,22,42])
                this%sites(20) = site(20, 4, [19,21,22,23])
                this%sites(21) = site(21, 4, [17,18,20,23])
                this%sites(22) = site(22, 4, [ 1,19,20,48])
                this%sites(23) = site(23, 4, [ 3,20,21,24])
                this%sites(24) = site(24, 4, [ 3,23,25,46])
                this%sites(25) = site(25, 4, [24,26,28,46])
                this%sites(26) = site(26, 4, [25,27,28,29])
                this%sites(27) = site(27, 4, [26,29,47,48])
                this%sites(28) = site(28, 4, [ 6,25,26,31])
                this%sites(29) = site(29, 4, [26,27,30,33])
                this%sites(30) = site(30, 4, [ 4, 7,29,33])
                this%sites(31) = site(31, 4, [ 6,28,32,34])
                this%sites(32) = site(32, 4, [31,33,34,35])
                this%sites(33) = site(33, 4, [29,30,32,35])
                this%sites(34) = site(34, 4, [12,31,32,37])
                this%sites(35) = site(35, 4, [32,33,36,39])
                this%sites(36) = site(36, 4, [10,13,35,39])
                this%sites(37) = site(37, 4, [12,34,38,40])
                this%sites(38) = site(38, 4, [37,39,40,41])
                this%sites(39) = site(39, 4, [35,36,38,41])
                this%sites(40) = site(40, 4, [18,37,38,43])
                this%sites(41) = site(41, 4, [38,39,42,45])
                this%sites(42) = site(42, 4, [16,19,41,45])
                this%sites(43) = site(43, 4, [18,40,44,46])
                this%sites(44) = site(44, 4, [43,45,46,47])
                this%sites(45) = site(45, 4, [41,42,44,47])
                this%sites(46) = site(46, 4, [24,25,43,44])
                this%sites(47) = site(47, 4, [27,44,45,48])
                this%sites(48) = site(48, 4, [ 1,22,27,47])

            else
                call stop_all(this_routine, &
                              "only 1x1,1x2,2x1, 2x2, 3x2 and 4x2 clusters implemented for kagome yet!")
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
        integer :: temp_array(4 * this%length(1), 2 * this%length(2))
        integer :: up(4 * this%length(1), 2 * this%length(2))
        integer :: down(4 * this%length(1), 2 * this%length(2))
        integer :: right(4 * this%length(1), 2 * this%length(2))
        integer :: left(4 * this%length(1), 2 * this%length(2))

        integer :: i, temp_neigh(3), x, y, special
        integer, allocatable :: neigh(:), order(:)
        character(*), parameter :: this_routine = "init_sites_hexagonal"

        allocate(order(this%get_nsites()), source = 0)

        if (t_input_order .and. this%t_bipartite_order) then
            order = orbital_order
        else
            order = [(i, i = 1, this%get_nsites())]
        end if

        temp_array = reshape([(order(i), i=1, this%get_nsites())], &
                             [4 * this%length(1), 2 * this%length(2)])

        up = cshift(temp_array, -1, 1)
        down = cshift(temp_array, 1, 1)
        left = cshift(temp_array, -1, 2)
        right = cshift(temp_array, 1, 2)

        if (this%is_periodic()) then
            do i = 1, this%get_nsites()
                ! columns and rows:
                x = mod(i - 1, 4 * this%length(1)) + 1
                y = (i - 1) / (4 * this%length(1)) + 1

                if (is_odd(y)) then
                    if (is_odd(i)) then
                        ! every odd number in a odd column has a right neighbor
                        special = right(x, y)
                    else
                        ! otherwise left
                        special = left(x, y)
                    end if
                else
                    ! for even columns it is the other way around
                    if (is_odd(i)) then
                        special = left(x, y)
                    else
                        special = right(x, y)
                    end if
                end if
                temp_neigh = [up(x, y), down(x, y), special]

                neigh = sort_unique(temp_neigh)

                this%sites(order(i)) = site(order(i), size(neigh), neigh)
            end do
        else
            call stop_all(this_routine, &
                          "closed boundary conditions not yet implemented for hexagonal lattice!")

        end if

    end subroutine init_sites_hexagonal

    subroutine init_sites_triangular(this)
        class(triangular) :: this
        integer :: temp_array(this%length(1), this%length(2))
        integer :: down(this%length(1), this%length(2))
        integer :: up(this%length(1), this%length(2))
        integer :: left(this%length(1), this%length(2))
        integer :: right(this%length(1), this%length(2))
        integer :: lu(this%length(1), this%length(2))
        integer :: rd(this%length(1), this%length(2))
        integer :: i, temp_neigh(6), x, y
        integer, allocatable :: neigh(:), order(:)
        character(*), parameter :: this_routine = "init_sites_triangular"

        ASSERT(this%get_nsites() >= 4)

        allocate(order(this%get_nsites()), source = 0)
        if (t_input_order .and. this%t_bipartite_order) then
            order = orbital_order
        else
            order = [(i, i = 1, this%get_nsites())]
        end if

        temp_array = reshape([(order(i), i=1, this%get_nsites())], this%length)

        up = cshift(temp_array, -1, 1)
        down = cshift(temp_array, 1, 1)
        right = cshift(temp_array, 1, 2)
        left = cshift(temp_array, -1, 2)
        lu = cshift(up, -1, 2)
        rd = cshift(down, 1, 2)

        if (this%is_periodic()) then
            do i = 1, this%get_nsites()

                x = mod(i - 1, this%length(1)) + 1
                y = (i - 1) / this%length(1) + 1

                temp_neigh = [up(x, y), down(x, y), left(x, y), right(x, y), lu(x, y), rd(x, y)]

                neigh = sort_unique(temp_neigh)

                this%sites(order(i)) = site(order(i), size(neigh), neigh)

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
        integer :: temp_array(this%length(1), this%length(2))
        integer :: down(this%length(1), this%length(2))
        integer :: up(this%length(1), this%length(2))
        integer :: left(this%length(1), this%length(2))
        integer :: right(this%length(1), this%length(2))
        integer :: unique(4)
        integer, allocatable :: neigh(:)
        integer :: sort_array_3(3), sort_array_2(2), sort_array(4)
        integer :: k_vec(3), r_vec(3)
        integer, allocatable :: order(:)

        ! this is the important routine..
        ! store lattice like that:
        ! 1 4 7
        ! 2 5 8
        ! 3 6 9

        ASSERT(this%get_nsites() >= 4)

        ! use cshift intrinsic of fortran..
        ! how do i efficiently set that up?
        if (this%t_bipartite_order) then
            if (this%length(1) /= this%length(2)) then
                if (this%length(2) /= 2) then
                    call stop_all(this_routine, &
                        "ladder bipartite ordering is implemented with Ly == 2)")
                end if

                allocate(order(this%get_nsites()), source = 0)
                if (t_input_order) then
                    order = orbital_order
                else
                    do i = 1, this%length(1)/2
                        order(2*i-1) = i
                        order(2*i) = this%length(1) + i
                    end do
                    do i = this%length(1)/2 + 1, this%length(1)
                        order(2*i-1) = this%length(1) + i
                        order(2*i) = i
                    end do
                end if
            else
                if (this%get_nsites() == 16) then
                    allocate(order(16), source = 0)
                    if (t_input_order) then
                        order = orbital_order
                    else
                        order = [ 1,  9,  2, 10, &
                                 11,  3, 12,  4, &
                                  5, 13,  6, 14, &
                                 15,  7, 16,  8]
                    end if
                else if (this%get_nsites() == 36) then

                    allocate(order(36), source = 0)
                    if (t_input_order) then
                        order = orbital_order
                    else
                        order = [ 1, 19,  2, 20,  3, 21, &
                                 22,  4, 23,  5, 24,  6, &
                                  7, 25,  8, 26,  9, 27, &
                                 28, 10, 29, 11, 30, 12, &
                                 13, 31, 14, 32, 15, 33, &
                                 34, 16, 35, 17, 36, 18]
                    end if
                else
                    if (t_input_order) then
                        order = orbital_order
                    else
                        call stop_all(this_routine, &
                            "bipartite order for square only implemented for 4x4! and 6x6 for now!")
                    end if
                end if
            end if
        else
            allocate(order(this%get_nsites()), source = [(i, i = 1, this%get_nsites())])
        end if

        temp_array = reshape( order, this%length)

        up = cshift(temp_array, -1, 1)
        down = cshift(temp_array, 1, 1)
        right = cshift(temp_array, 1, 2)
        left = cshift(temp_array, -1, 2)

        if (this%is_periodic()) then

            do i = 1, this%get_nsites()
                ! create the neighbor list
                x = mod(i - 1, this%length(1)) + 1
                y = (i - 1) / this%length(1) + 1

                temp_neigh = [up(x, y), down(x, y), left(x, y), right(x, y)]

                neigh = sort_unique(temp_neigh)

                k_vec = [x - (this%length(1) + 1) / 2, y - (this%length(2) + 1) / 2, 0]
                r_vec = [x - (this%length(1) + 1) / 2, y - (this%length(2) + 1) / 2, 0]
                this%sites(order(i)) = site(order(i), size(neigh), neigh, k_vec, r_vec)

                deallocate(neigh)

            end do

        else if (this%is_periodic(1)) then
            ! only periodic in x-direction
            do i = 1, this%get_nsites()

                x = mod(i - 1, this%length(1)) + 1
                y = (i - 1) / this%length(1) + 1

                ! now definetly always take the left and right neighbors
                ! but up and down if we are not jumping boundaries
                if (x == 1) then
                    ! dont take upper neighbor -> just repeat an neighbor,
                    ! which will get removed from unique
                    temp_neigh = [down(x, y), down(x, y), left(x, y), right(x, y)]
                else if (x == this%length(1)) then
                    temp_neigh = [up(x, y), up(x, y), left(x, y), right(x, y)]

                else
                    ! take all
                    temp_neigh = [up(x, y), down(x, y), left(x, y), right(x, y)]
                end if

                neigh = sort_unique(temp_neigh)

                k_vec = [x - (this%length(1) + 1) / 2, y - (this%length(2) + 1) / 2, 0]
                r_vec = [x - (this%length(1) + 1) / 2, y - (this%length(2) + 1) / 2, 0]

                this%sites(order(i)) = site(order(i), size(neigh), neigh, k_vec, r_vec)

                deallocate(neigh)
            end do

        else if (this%is_periodic(2)) then
            ! only periodic in the y-direction
            do i = 1, this%get_nsites()

                x = mod(i - 1, this%length(1)) + 1
                y = (i - 1) / this%length(1) + 1

                ! now definetly always take the left and right neighbors
                ! but up and down if we are not jumping boundaries
                if (y == 1) then
                    ! dont take upper neighbor -> just repeat an neighbor,
                    ! which will get removed from unique
                    temp_neigh = [up(x, y), down(x, y), right(x, y), right(x, y)]
                else if (y == this%length(2)) then
                    temp_neigh = [up(x, y), down(x, y), left(x, y), left(x, y)]

                else
                    ! take all
                    temp_neigh = [up(x, y), down(x, y), left(x, y), right(x, y)]
                end if

                neigh = sort_unique(temp_neigh)

                k_vec = [x - (this%length(1) + 1) / 2, y - (this%length(2) + 1) / 2, 0]
                r_vec = [x - (this%length(1) + 1) / 2, y - (this%length(2) + 1) / 2, 0]

                this%sites(order(i)) = site(order(i), size(neigh), neigh, k_vec, r_vec)

                deallocate(neigh)
            end do

        else
            ! non-periodic
            do i = 1, this%get_nsites()

                x = mod(i - 1, this%length(1)) + 1
                y = (i - 1) / this%length(1) + 1

                ! now definetly always take the left and right neighbors
                ! but up and down if we are not jumping boundaries
                if (x == 1 .and. y == 1) then
                    ! dont take upper neighbor -> just repeat an neighbor,
                    ! which will get removed from unique
                    temp_neigh = [down(x, y), down(x, y), right(x, y), right(x, y)]
                else if (y == this%length(2) .and. x == this%length(1)) then
                    temp_neigh = [up(x, y), up(x, y), left(x, y), left(x, y)]

                else if (x == this%length(1) .and. y == 1) then
                    temp_neigh = [up(x, y), up(x, y), right(x, y), right(x, y)]
                else if (y == this%length(2) .and. x == 1) then
                    temp_neigh = [down(x, y), down(x, y), left(x, y), left(x, y)]

                else if (x == 1) then
                    temp_neigh = [left(x, y), right(x, y), down(x, y), down(x, y)]

                else if (x == this%length(1)) then
                    temp_neigh = [left(x, y), right(x, y), up(x, y), up(x, y)]

                else if (y == 1) then
                    temp_neigh = [right(x, y), right(x, y), up(x, y), down(x, y)]

                else if (y == this%length(2)) then
                    temp_neigh = [left(x, y), left(x, y), up(x, y), down(x, y)]

                else
                    ! take all
                    temp_neigh = [up(x, y), down(x, y), left(x, y), right(x, y)]
                end if

                neigh = sort_unique(temp_neigh)

                k_vec = [x - (this%length(1) + 1) / 2, y - (this%length(2) + 1) / 2, 0]
                r_vec = [x - (this%length(1) + 1) / 2, y - (this%length(2) + 1) / 2, 0]

                this%sites(order(i)) = site(order(i), size(neigh), neigh, k_vec, r_vec)

                deallocate(neigh)
            end do

        end if

    end subroutine init_sites_rect

    subroutine read_sites(this)

        class(ext_input):: this

        character(*), parameter :: this_routine = "read_sites"

        integer:: iunit, ios, i, n_site, n_neighbors, neigh
        integer, allocatable :: neighs(:)

        logical :: exists, leof

        CHARACTER(len=3) :: fmat
        CHARACTER(LEN=100) w
        type(ManagingFileReader_t) :: file_reader
        type(TokenIterator_t) :: tokens

        file_reader = ManagingFileReader_t("lattice.file")

        readsites: do while (file_reader%nextline(tokens))
            w = to_upper(tokens%next())
            select case (w)
            case ('SITE')
                n_site = to_int(tokens%next())

                n_neighbors = to_int(tokens%next())
                if (allocated(neighs)) deallocate(neighs)
                allocate(neighs(n_neighbors), source=0)
                do i = 1, size(neighs)
                    neighs(i) = to_int(tokens%next())
                end do
                this%sites(n_site) = site(n_site, n_neighbors, neighs)
            end select
        end do readsites


    end subroutine read_sites

    subroutine init_sites_sujun(this)
        ! order of the lattice sites
        !   4 7 10
        !   3 6 9
        ! 1 2 5 8
        !
        class(sujun) :: this
        character(*), parameter :: this_routine = "init_sites_sujun"
        this%sites(1) = site(1, 4, [2,4,8,10])
        this%sites(2) = site(2, 4, [1,3,5,7])
        this%sites(3) = site(3, 4, [2,4,6,8])
        this%sites(4) = site(4, 4, [1,3,7,9])
        this%sites(5) = site(5, 4, [2,6,8,10])
        this%sites(6) = site(6, 4, [3,5,7,9])
        this%sites(7) = site(7, 4, [2,4,6,10])
        this%sites(8) = site(8, 4, [1,3,5,9])
        this%sites(9) = site(9, 4, [4,6,8,10])
        this%sites(10) = site(10,4, [1,5,7,9])

    end subroutine init_sites_sujun

    subroutine init_sites_ole(this)
        class(ole) :: this
        character(*), parameter :: this_routine = "init_sites_ole"

        ! i think i can index an array in fortran in reversed order
        integer :: ind_array(-this%length(1):(this%length(1) + 1), &
                             -this%length(2):this%length(1))

        integer :: i, j, k, mat_ind(this%n_sites, 2), up, down, left, right, &
                   k_vec(3), A(2), B(2), C(2), D(2), k_vec_prep(24, 3)
        integer, allocatable :: neigh(:)

        ! how do i set up Ole cluster..
        ! in real and k-space.. this will be a pain i guess..

        ! i could make the same neighboring matrices as for the tilted
        k = 1
        ind_array = 0
        ! define the vertices of the parallelogram
        A = [-this%length(2), 0]
        B = [0, -this%length(1)]
        C = [this%length(1), 0]
        D = [this%length(1) - this%length(2), this%length(1)]

        ! i do not need to loop over the edge values, which belong to a
        ! different unit cell!
        do j = -this%length(2) + 1, this%length(1) - 1
            do i = this%length(1) - 1, -(this%length(1) + 1), -1

                ! i have to take 2 special points into account ! which are
                ! by definition on the edge of the (3,3) boundary
                if (inside_bz_2d(j, i, A, B, C, D) .and. .not. on_line_2d([j, i], A, D) &
                    .and. .not. on_line_2d([j, i], C, D)) then

                    ind_array(-i, j) = k
                    mat_ind(k, :) = [-i, j]

                    k = k + 1

                end if
            end do
        end do

        k_vec_prep(1, :) = [1, -3, 0]
        k_vec_prep(2, :) = [-2, 1, 0]
        k_vec_prep(3, :) = [-2, 0, 0]

        k_vec_prep(4, :) = [-1, 2, 0]
        k_vec_prep(5, :) = [-1, 1, 0]
        k_vec_prep(6, :) = [-1, 0, 0]
        k_vec_prep(7, :) = [-1, -1, 0]
        k_vec_prep(8, :) = [-1, -2, 0]

        k_vec_prep(9, :) = [0, 3, 0]
        k_vec_prep(10, :) = [0, 2, 0]
        k_vec_prep(11, :) = [0, 1, 0]
        k_vec_prep(12, :) = [0, 0, 0]
        k_vec_prep(13, :) = [0, -1, 0]
        k_vec_prep(14, :) = [0, -2, 0]
        k_vec_prep(15, :) = [0, -3, 0]

        k_vec_prep(16, :) = [1, 2, 0]
        k_vec_prep(17, :) = [1, 1, 0]
        k_vec_prep(18, :) = [1, 0, 0]
        k_vec_prep(19, :) = [1, -1, 0]
        k_vec_prep(20, :) = [1, -2, 0]

        k_vec_prep(21, :) = [2, 0, 0]
        k_vec_prep(22, :) = [2, -1, 0]
        k_vec_prep(23, :) = [2, -2, 0]

        k_vec_prep(24, :) = [3, -1, 0]

        ! now i want to get the neigbhors
        do i = 1, this%get_nsites()
            ! how to efficiently do this, and in a general way?
            ! better than in the other cases
            ! i am not going over the boundaries, due to the way i set up
            ! the matrix above.. hopefully
            up = ind_array(mat_ind(i, 1) - 1, mat_ind(i, 2))
            if (up == 0) then
                up = this%find_periodic_neighbors([mat_ind(i, 1) - 1, mat_ind(i, 2)], &
                                                  ind_array)
            end if

            down = ind_array(mat_ind(i, 1) + 1, mat_ind(i, 2))
            if (down == 0) then
                down = this%find_periodic_neighbors([mat_ind(i, 1) + 1, mat_ind(i, 2)], &
                                                    ind_array)
            end if

            left = ind_array(mat_ind(i, 1), mat_ind(i, 2) - 1)
            if (left == 0) then
                left = this%find_periodic_neighbors([mat_ind(i, 1), mat_ind(i, 2) - 1], &
                                                    ind_array)
            end if

            right = ind_array(mat_ind(i, 1), mat_ind(i, 2) + 1)
            if (right == 0) then
                right = this%find_periodic_neighbors([mat_ind(i, 1), mat_ind(i, 2) + 1], &
                                                     ind_array)
            end if

            neigh = sort_unique([up, down, left, right])

            ! i have to get the matrix indiced again, with the correct
            ! sign..
            if (this%get_nsites() == 24) then
                k_vec = k_vec_prep(i, :)
            else
                k_vec = [mat_ind(i, 2), -mat_ind(i, 1), 0]
            end if

            this%sites(i) = site(i, size(neigh), neigh, k_vec)

        end do

    end subroutine init_sites_ole

    integer function find_periodic_neighbors_ole(this, ind, A)
        ! function to give me a periodic neighbor of a specific lattice
        class(ole) :: this
        integer, intent(in) :: ind(2), A(:, :)
        integer :: ur(-this%length(1):(this%length(1) + 1), &
                      -this%length(2):this%length(1))
        integer :: dr(-this%length(1):(this%length(1) + 1), &
                      -this%length(2):this%length(1))
        integer :: ul(-this%length(1):(this%length(1) + 1), &
                      -this%length(2):this%length(1))
        integer :: dl(-this%length(1):(this%length(1) + 1), &
                      -this%length(2):this%length(1))
        integer :: rr(-this%length(1):(this%length(1) + 1), &
                      -this%length(2):this%length(1))
        integer :: ll(-this%length(1):(this%length(1) + 1), &
                      -this%length(2):this%length(1))
        integer :: temp(-this%length(1):(this%length(1) + 1), &
                        -this%length(2):this%length(1))

        integer :: unique, shift(4, 2), i
        ! i am not sure, if i have to specify the indices and size of the
        ! matrix inputted..

        ! get the lattice vectors:
        associate(r1 => this%lat_vec(1:2, 1), r2 => this%lat_vec(1:2, 2), &
                   x => ind(1), y => ind(2))

            shift = transpose(reshape([r1, r2, r1 + r2, r1 - r2], [2, 4]))

            find_periodic_neighbors_ole = -1

            do i = 1, 4
                ! apply all the periodic vectors one after the other
                ! negative and positive..
                temp = apply_pbc(A, shift(i, :))

                if (temp(x, y) /= 0) then
                    find_periodic_neighbors_ole = temp(x, y)
                    return
                end if

                temp = apply_pbc(A, -shift(i, :))

                if (temp(x, y) /= 0) then
                    find_periodic_neighbors_ole = temp(x, y)
                    return
                end if
            end do

        end associate

        find_periodic_neighbors_ole = unique

    end function find_periodic_neighbors_ole

    function apply_pbc(array, shift) result(s_array)
        integer, intent(in) :: array(:, :), shift(2)
        integer :: s_array(size(array, 1), size(array, 2))

        ! i have to be sure about the sign conventions here..
        s_array = eoshift(array, shift(1), dim=2)
        s_array = eoshift(s_array, -shift(2), dim=1)

    end function apply_pbc

    logical function inside_bz_2d(x, y, A, B, C, D)
        ! function to check if a point (x,y) is inside our outside the
        ! rectangle indicated by the 3 points A,B,C,D
        ! the definition is to provide the points in this fashion:
        !  A -- D
        !  |    |
        !  B -- C
        ! in a circular fashion, otherwise it does not work, since it would
        ! be a different rectangle
        ! idea from:
        ! https://stackoverflow.com/questions/2752725/finding-whether-a-point-lies-inside-a-rectangle-or-not
        integer, intent(in) :: x, y, A(2), B(2), C(2), D(2)

        integer :: vertex(4, 2), edges(4, 2), R(4), S(4), T(4), U(4), i
        inside_bz_2d = .false.

        vertex = transpose(reshape([A, B, C, D], [2, 4]))

        edges(1, :) = A - B
        edges(2, :) = B - C
        edges(3, :) = C - D
        edges(4, :) = D - A

        R = edges(:, 2)
        S = -edges(:, 1)
        T = -(R * vertex(:, 1) + S * vertex(:, 2))

        U = R * x + S * y + T

        if (all(U >= 0)) inside_bz_2d = .true.

    end function inside_bz_2d

    logical function on_line_2d(P, A, B)
        integer, intent(in) :: P(2), A(2), B(2)
        ! function to check if a point is on a line(for integers now only!)

        integer :: AB(2), AP(2)

        AB = B - A
        AP = P - A

        on_line_2d = .false.

        if (AB(1) * AP(2) - AB(2) * AP(1) == 0) on_line_2d = .true.

    end function on_line_2d

    subroutine init_sites_tilted(this)
        class(tilted) :: this
        character(*), parameter :: this_routine = "init_sites_tilted"

        integer :: temp_array(-this%length(1):this%length(2), &
                              -this%length(1):this%length(2) + 1)
        integer :: up(-this%length(1):this%length(2), &
                      -this%length(1):this%length(2) + 1)
        integer :: down(-this%length(1):this%length(2), &
                        -this%length(1):this%length(2) + 1)
        integer :: left(-this%length(1):this%length(2), &
                        -this%length(1):this%length(2) + 1)
        integer :: right(-this%length(1):this%length(2), &
                         -this%length(1):this%length(2) + 1)
        integer :: right_ul(-this%length(1):this%length(2), &
                            -this%length(1):this%length(2) + 1)
        integer :: right_ur(-this%length(1):this%length(2), &
                            -this%length(1):this%length(2) + 1)
        integer :: right_dl(-this%length(1):this%length(2), &
                            -this%length(1):this%length(2) + 1)
        integer :: right_dr(-this%length(1):this%length(2), &
                            -this%length(1):this%length(2) + 1)
        integer :: right_rr(-this%length(1):this%length(2), &
                            -this%length(1):this%length(2) + 1)
        integer :: right_ll(-this%length(1):this%length(2), &
                            -this%length(1):this%length(2) + 1)
        integer :: up_ul(-this%length(1):this%length(2), &
                         -this%length(1):this%length(2) + 1)
        integer :: up_ur(-this%length(1):this%length(2), &
                         -this%length(1):this%length(2) + 1)
        integer :: up_dl(-this%length(1):this%length(2), &
                         -this%length(1):this%length(2) + 1)
        integer :: up_dr(-this%length(1):this%length(2), &
                         -this%length(1):this%length(2) + 1)
        integer :: up_rr(-this%length(1):this%length(2), &
                         -this%length(1):this%length(2) + 1)
        integer :: up_ll(-this%length(1):this%length(2), &
                         -this%length(1):this%length(2) + 1)
        integer :: down_ul(-this%length(1):this%length(2), &
                           -this%length(1):this%length(2) + 1)
        integer :: down_ur(-this%length(1):this%length(2), &
                           -this%length(1):this%length(2) + 1)
        integer :: down_dl(-this%length(1):this%length(2), &
                           -this%length(1):this%length(2) + 1)
        integer :: down_dr(-this%length(1):this%length(2), &
                           -this%length(1):this%length(2) + 1)
        integer :: down_rr(-this%length(1):this%length(2), &
                           -this%length(1):this%length(2) + 1)
        integer :: down_ll(-this%length(1):this%length(2), &
                           -this%length(1):this%length(2) + 1)
        integer :: left_ul(-this%length(1):this%length(2), &
                           -this%length(1):this%length(2) + 1)
        integer :: left_ur(-this%length(1):this%length(2), &
                           -this%length(1):this%length(2) + 1)
        integer :: left_dl(-this%length(1):this%length(2), &
                           -this%length(1):this%length(2) + 1)
        integer :: left_dr(-this%length(1):this%length(2), &
                           -this%length(1):this%length(2) + 1)
        integer :: left_rr(-this%length(1):this%length(2), &
                           -this%length(1):this%length(2) + 1)
        integer :: left_ll(-this%length(1):this%length(2), &
                           -this%length(1):this%length(2) + 1)
        integer :: i, j, k, l, pbc, temp_neigh(4), k_min, k_max, offset, k_vec(3), m
        integer :: right_nn, left_nn, up_nn, down_nn, pbc_1(2), pbc_2(2), r_vec(3)
        integer, allocatable :: neigh(:)
        integer, allocatable :: order(:)
        ! convention of lattice storage:
        !
        !   2 5
        ! 1 3 6 8
        !   4 7
        ! update: we also want to have non-square tilted clusters.. how would
        ! that work. eg a 10-site 2x3 cluster:
        !   2 5
        ! 1 3 6 9
        !   4 7 10
        !     8
        ! can i also do a 1x2 4-site tilted? like:
        ! 1 2
        !   3 4
        ! and also a 2x1:
        !
        ! 1 3

        ASSERT(this%get_nsites() >= 4)

        ! set up the lattice indices, via the use of "k-vectors"
        temp_array(:, :) = 0
        if (this%t_bipartite_order) then
            if ( .not. (this%get_nsites() == 18 .or. this%get_nsites() == 8)) then
                call stop_all(this_routine, &
                    "bipartite only for 8 or 18 tilted sites for now")
            end if
            allocate(order(this%get_nsites()), source = 0)

            if (t_input_order) then
                order = orbital_order
            else
                if (this%get_nsites() == 18) then
                    order = [ 1, 2, 10, 3, 4, 11, 5, 12, 6, 13, 7, 14, 8, 15, 16, 9, 17, 18]
                else if (this%get_nsites() == 8) then
                    order = [1,2,5,3,6,4,7,8]
                end if
            end if
        else
            allocate(order(this%get_nsites()), source = [(i, i = 1, this%get_nsites())])
        end if


        k = 0
        l = 1
        m = this%get_nsites() / 2 + 1
        do i = -this%length(1) + 1, 0
            do j = -k, k

                temp_array(j, i) = order(l)
                l = l + 1

            end do
            k = k + 1
        end do

        ! here i need to change the k-vectors, differently, if it is a
        ! rectangular tilted lattice..
        ! and for now, until i have implemented it better over
        ! lattice vectors assume lx - ly <= 1! only one difference
        k = k - 1
        ! or should i do an inbetween-step if lx /= ly? this is also
        ! possible

        offset = abs(this%length(1) - this%length(2))
        k_min = -this%length(1) + 1
        k_max = this%length(2) - offset
        do i = 1, offset

            do j = k_min, k_max
                temp_array(j, i) = l
                l = l + 1
            end do

            ! shift the y indication by 1 up or down
            k_min = k_min + 1
            k_max = k_max + 1
        end do

        if (this%length(1) < this%length(2)) then
            k_min = k_min
            k_max = k_max - 1
        else if (this%length(1) > this%length(2)) then
            k_min = k_min + 1
            k_max = k_max
        else
            ! otherwise k_min and k_max where never defined
            k_min = -k
            k_max = k
        end if

        do i = offset + 1, this%length(2)


            ! if (this%t_bipartite_order) then
            !     do j = k_min, k_max, 2
            !
            !         temp_array(j, i) = m
            !
            !         m = m + 1
            !
            !     end do
            !     do j = k_min + 1, k_max - 1, 2
            !         temp_array(j, i) = l
            !         l = l + 1
            !     end do
            ! else
                do j = k_min, k_max

                    temp_array(j, i) = order(l)

                    l = l + 1

                end do
            ! end if

            ! k_min is negative
            k_min = k_min + 1
            k_max = k_max - 1
        end do

        up = cshift(temp_array, -1, 1)
        down = cshift(temp_array, 1, 1)
        right = cshift(temp_array, 1, 2)
        left = cshift(temp_array, -1, 2)

        ! apply the periodic boundary conditions to the neighbors
        pbc = this%length(1)
        ! for rectangular tilted lattices this is different of course

        if (this%length(1) == 1) then
            pbc_1 = [2, 0]
            pbc_2 = [this%length(2), -this%length(2)]
        else if (this%length(2) == 1) then
            pbc_1 = [this%length(1), this%length(1)]
            pbc_2 = [2, 0]
        else
            pbc_1 = [this%length(1), this%length(1)]
            pbc_2 = [this%length(2), -this%length(2)]
        end if

        ! do something like and do this generally maybe..
        call apply_pbc_tilted(up, pbc_1, pbc_2, up_ur, up_dr, up_ul, up_dl, up_rr, up_ll)
        call apply_pbc_tilted(down, pbc_1, pbc_2, down_ur, down_dr, down_ul, &
            down_dl, down_rr, down_ll)
        call apply_pbc_tilted(right, pbc_1, pbc_2, right_ur, right_dr, right_ul, &
            right_dl, right_rr, right_ll)
        call apply_pbc_tilted(left, pbc_1, pbc_2, left_ur, left_dr, left_ul, &
            left_dl, left_rr, left_ll)

        k = 0
        l = 1
        ! now get the neighbors
        if (this%is_periodic()) then
            ! fully periodic case
            do i = -this%length(1) + 1, 0
                do j = -k, k
                    ! make the neigbors list
                    up_nn = maxval([up(j, i), up_ur(j, i), up_dr(j, i), up_ul(j, i), &
                                    up_dl(j, i)])

                    if (up_nn == 0) then
                        up_nn = maxval([up_rr(j, i), up_ll(j, i)])
                        if (up_nn == 0) then
                            print *, " up: smth wrong!"
                        end if
                    end if

                    down_nn = maxval([down(j, i), down_ur(j, i), down_dr(j, i), &
                        down_ul(j, i), down_dl(j, i)])

                    if (down_nn == 0) then
                        down_nn = maxval([down_rr(j, i), down_ll(j, i)])

                        if (down_nn == 0) then
                            print *, "down: smth wrong!"
                        end if
                    end if

                    right_nn = maxval([right(j, i), right_ur(j, i), right_dr(j, i), &
                        right_ul(j, i), right_dl(j, i)])

                    if (right_nn == 0) then
                        right_nn = right_ll(j, i)
                        if (right_nn == 0) then
                            print *, "right: smth wrong!"
                        end if
                    end if

                    left_nn = maxval([left(j, i), left_ur(j, i), left_dr(j, i), &
                        left_ul(j, i), left_dl(j, i)])

                    if (left_nn == 0) then
                        left_nn = left_rr(j, i)
                        if (left_nn == 0) then
                            print *, "left: smth wrong!"
                        end if
                    end if

                    neigh = sort_unique([up_nn, down_nn, left_nn, right_nn])
                    ! also start to store the k-vector here!
                    ! have to be sure that i make it correct
                    k_vec = [i, j, 0]
                    r_vec = [j, i, 0]
                    this%sites(order(l)) = site(order(l), size(neigh), neigh, k_vec, r_vec)

                    l = l + 1

                    deallocate(neigh)

                end do
                k = k + 1
            end do

            k = k - 1

            k_min = -this%length(1) + 1
            k_max = this%length(2) - offset
            do i = 1, offset
                do j = k_min, k_max
                    ! make the neigbors list
                    up_nn = maxval([up(j, i), up_ur(j, i), up_dr(j, i), up_ul(j, i), &
                                    up_dl(j, i)])

                    if (up_nn == 0) then
                        up_nn = maxval([up_rr(j, i), up_ll(j, i)])
                        if (up_nn == 0) then
                            print *, "smth wrong!"
                        end if
                    end if

                    down_nn = maxval([down(j, i), down_ur(j, i), down_dr(j, i), &
                        down_ul(j, i), down_dl(j, i)])

                    if (down_nn == 0) then
                        down_nn = maxval([down_rr(j, i), down_ll(j, i)])

                        if (down_nn == 0) then
                            print *, "smth wrong!"
                        end if
                    end if

                    right_nn = maxval([right(j, i), right_ur(j, i), right_dr(j, i),&
                        right_ul(j, i), right_dl(j, i)])

                    if (right_nn == 0) then
                        right_nn = right_ll(j, i)
                        if (right_nn == 0) then
                            print *, "smth wrong!"
                        end if
                    end if

                    left_nn = maxval([left(j, i), left_ur(j, i), left_dr(j, i), &
                        left_ul(j, i), left_dl(j, i)])

                    if (left_nn == 0) then
                        left_nn = left_rr(j, i)
                        if (left_nn == 0) then
                            print *, "smth wrong!"
                        end if
                    end if

                    neigh = sort_unique([up_nn, down_nn, left_nn, right_nn])

                    k_vec = [i, j, 0]
                    r_vec = [j, 1, 0]
                    this%sites(order(l)) = site(order(l), size(neigh), neigh, k_vec, r_vec)

                    l = l + 1

                    deallocate(neigh)

                end do
                k_min = k_min + 1
                k_max = k_max + 1
            end do

            if (this%length(1) < this%length(2)) then
                k_min = k_min
                k_max = k_max - 1
            else if (this%length(1) > this%length(2)) then
                k_min = k_min + 1
                k_max = k_max
            else
                k_min = -k
                k_max = k
            end if

            do i = offset + 1, this%length(2)
                do j = k_min, k_max
                    ! make the neigbors list
                    up_nn = maxval([up(j, i), up_ur(j, i), up_dr(j, i), up_ul(j, i), &
                                    up_dl(j, i)])

                    if (up_nn == 0) then
                        up_nn = maxval([up_rr(j, i), up_ll(j, i)])
                        if (up_nn == 0) then
                            print *, "smth wrong!"
                        end if
                    end if

                    down_nn = maxval([down(j, i), down_ur(j, i), down_dr(j, i), &
                        down_ul(j, i), down_dl(j, i)])

                    if (down_nn == 0) then
                        down_nn = maxval([down_rr(j, i), down_ll(j, i)])

                        if (down_nn == 0) then
                            print *, "smth wrong!"
                        end if
                    end if

                    right_nn = maxval([right(j, i), right_ur(j, i), right_dr(j, i), &
                        right_ul(j, i), right_dl(j, i)])

                    if (right_nn == 0) then
                        right_nn = right_ll(j, i)
                        if (right_nn == 0) then
                            print *, "smth wrong!"
                        end if
                    end if

                    left_nn = maxval([left(j, i), left_ur(j, i), left_dr(j, i), &
                        left_ul(j, i), left_dl(j, i)])

                    if (left_nn == 0) then
                        left_nn = left_rr(j, i)
                        if (left_nn == 0) then
                            print *, "smth wrong!"
                        end if
                    end if

                    neigh = sort_unique([up_nn, down_nn, left_nn, right_nn])

                    k_vec = [i, j, 0]
                    r_vec = [j, i, 0]
                    this%sites(order(l)) = site(order(l), size(neigh), neigh, k_vec, r_vec)

                    l = l + 1

                    deallocate(neigh)

                end do
                k_min = k_min + 1
                k_max = k_max - 1

            end do
        else if (this%is_periodic(1)) then
            ! only apply (x,x) periodicity
            do i = -this%length(1) + 1, 0
                do j = -k, k

                    up_nn = maxval([up(j, i), up_ur(j, i), up_dl(j, i)])
                    down_nn = maxval([down(j, i), down_ur(j, i), down_dl(j, i)])
                    left_nn = maxval([left(j, i), left_ur(j, i), left_dl(j, i)])
                    right_nn = maxval([right(j, i), right_ur(j, i), right_dl(j, i)])

                    temp_neigh = [up_nn, down_nn, left_nn, right_nn]

                    neigh = sort_unique(pack(temp_neigh, temp_neigh > 0))

                    this%sites(order(l)) = site(order(l), size(neigh), neigh)

                    l = l + 1

                    deallocate(neigh)
                end do
                k = k + 1
            end do
            k = k - 1

            do i = 1, this%length(1)
                do j = -k, k

                    up_nn = maxval([up(j, i), up_ur(j, i), up_dl(j, i)])
                    down_nn = maxval([down(j, i), down_ur(j, i), down_dl(j, i)])
                    left_nn = maxval([left(j, i), left_ur(j, i), left_dl(j, i)])
                    right_nn = maxval([right(j, i), right_ur(j, i), right_dl(j, i)])

                    temp_neigh = [up_nn, down_nn, left_nn, right_nn]

                    neigh = sort_unique(pack(temp_neigh, temp_neigh > 0))

                    this%sites(order(l)) = site(order(l), size(neigh), neigh)

                    l = l + 1

                    deallocate(neigh)
                end do
                k = k - 1
            end do

        else if (this%is_periodic(2)) then
            ! only apply (x,-x) periodicity
            do i = -this%length(1) + 1, 0
                do j = -k, k

                    up_nn = maxval([up(j, i), up_ul(j, i), up_dr(j, i)])
                    down_nn = maxval([down(j, i), down_ul(j, i), down_dr(j, i)])
                    left_nn = maxval([left(j, i), left_ul(j, i), left_dr(j, i)])
                    right_nn = maxval([right(j, i), right_ul(j, i), right_dr(j, i)])

                    temp_neigh = [up_nn, down_nn, left_nn, right_nn]

                    neigh = sort_unique(pack(temp_neigh, temp_neigh > 0))

                    this%sites(order(l)) = site(order(l), size(neigh), neigh)

                    l = l + 1

                    deallocate(neigh)
                end do
                k = k + 1
            end do
            k = k - 1

            do i = 1, this%length(1)
                do j = -k, k

                    up_nn = maxval([up(j, i), up_ul(j, i), up_dr(j, i)])
                    down_nn = maxval([down(j, i), down_ul(j, i), down_dr(j, i)])
                    left_nn = maxval([left(j, i), left_ul(j, i), left_dr(j, i)])
                    right_nn = maxval([right(j, i), right_ul(j, i), right_dr(j, i)])

                    temp_neigh = [up_nn, down_nn, left_nn, right_nn]

                    neigh = sort_unique(pack(temp_neigh, temp_neigh > 0))

                    this%sites(order(l)) = site(order(l), size(neigh), neigh)

                    l = l + 1

                    deallocate(neigh)
                end do
                k = k - 1
            end do

        else
            ! non-periodic case
            do i = -this%length(1) + 1, 0
                do j = -k, k
                    ! only a neighbor if the index is non-zero!
                    temp_neigh = [up(j, i), down(j, i), left(j, i), right(j, i)]

                    neigh = sort_unique(pack(temp_neigh, temp_neigh > 0))

                    this%sites(order(l)) = site(order(l), size(neigh), neigh)

                    l = l + 1

                    deallocate(neigh)

                end do
                k = k + 1
            end do

            k = k - 1
            do i = 1, this%length(1)
                do j = -k, k
                    ! only a neighbor if the index is non-zero!
                    temp_neigh = [up(j, i), down(j, i), left(j, i), right(j, i)]

                    neigh = sort_unique(pack(temp_neigh, temp_neigh > 0))

                    this%sites(order(l)) = site(order(l), size(neigh), neigh)

                    l = l + 1

                    deallocate(neigh)

                end do
                k = k - 1
            end do
        end if

    end subroutine init_sites_tilted

    subroutine apply_pbc_tilted(array, pbc_1, pbc_2, ur, dr, ul, dl, rr, ll)
        ! intermediate routine to apply the periodic boundary conditions for
        ! the rectangular tilted lattice
        integer, intent(in) :: array(:, :), pbc_1(2), pbc_2(2)
        integer, intent(out) :: ur(:, :), dr(:, :), ul(:, :), dl(:, :), rr(:, :), ll(:, :)

        integer :: plus(2)

        plus = pbc_1 + pbc_2
        ! i have to do this properly, so it still works with the old
        ! tilted implementation:
        ! cshfting along the 2nd dimension is the x-axis shift on my
        ! derivation
        ! and along the 1-direction i should shift with -pbc to move it
        ! actually up..
        ! i think something is wrong with this routine still!
        rr = cshift(cshift(array, -plus(2), 1), plus(1), 2)
        ll = cshift(cshift(array, plus(2), 1), -plus(1), 2)

        ur = cshift(cshift(array, -pbc_1(2), 1), pbc_1(1), 2)
        dr = cshift(cshift(array, -pbc_2(2), 1), pbc_2(1), 2)

        ul = cshift(cshift(array, pbc_2(2), 1), -pbc_2(1), 2)
        dl = cshift(cshift(array, pbc_1(2), 1), -pbc_1(1), 2)

    end subroutine apply_pbc_tilted

    pure function get_k_vec(this, orb) result(k_vec)
        class(lattice), intent(in) :: this
        integer, intent(in) :: orb
        integer :: k_vec(3)

        k_vec = this%sites(orb)%k_vec

    end function get_k_vec

    function get_r_vec(this, orb) result(r_vec)
        class(lattice) :: this
        integer, intent(in) :: orb
        integer :: r_vec(3)

        r_vec = this%sites(orb)%r_vec

    end function get_r_vec

    function dispersion_rel_chain_k(this, k_vec) result(disp)
        class(chain) :: this
        integer, intent(in) :: k_vec(3)
        real(dp) :: disp
#ifdef DEBUG_
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
        disp = 2.0_dp * cos(2.0_dp * pi * (k_vec(1) + twisted_bc(1)) / this%length)

    end function dispersion_rel_chain_k

    function dispersion_rel_chain_orb(this, orb) result(disp)
        class(chain) :: this
        integer, intent(in) :: orb
        real(dp) :: disp

        integer :: k_vec(3)

        k_vec = this%get_k_vec(orb)

        disp = this%dispersion_rel(k_vec)

    end function dispersion_rel_chain_orb

    function dispersion_rel_rect(this, k_vec) result(disp)
        class(rectangle) :: this
        integer, intent(in) :: k_vec(3)
        real(dp) :: disp
#ifdef DEBUG_
        character(*), parameter :: this_routine = "dispersion_rel_rect"
#endif

        ASSERT(this%is_periodic())

        disp = 2.0_dp * (cos(2 * pi * (k_vec(1) + twisted_bc(1)) / this%length(1)) &
                         + cos(2 * pi * (k_vec(2) + twisted_bc(2)) / this%length(2)))

    end function dispersion_rel_rect

    function dispersion_rel_cube(this, k_vec) result(disp)
        class(cube) :: this
        integer, intent(in) :: k_vec(3)
        real(dp) :: disp
#ifdef DEBUG_
        character(*), parameter :: this_routine = "dispersion_rel_cube"
#endif

        ASSERT(this%is_periodic())

        disp = 2.0_dp * (cos(2 * pi * (k_vec(1) + twisted_bc(1)) / this%length(1)) &
                         + cos(2 * pi * (k_vec(2) + twisted_bc(2)) / this%length(2)) &
                         + cos(2 * pi * (k_vec(3) + twisted_bc(3)) / this%length(3)))

    end function dispersion_rel_cube

    function dispersion_rel_tilted(this, k_vec) result(disp)
        class(tilted) :: this
        integer, intent(in) :: k_vec(3)
        real(dp) :: disp
#ifdef DEBUG_
        character(*), parameter :: this_routine = "dispersion_rel_tilted"
#endif

        ASSERT(this%is_periodic())

        ! todo: i have to check if this also still holds for the
        ! rectangular tilted lattice!

        ! after some more consideration i believe this is the correct:
        ! although now i am not sure about the twist anymore... check that!
        disp = 2.0_dp * (cos(pi * ((k_vec(1) + twisted_bc(1)) / this%length(1) &
                                   + (k_vec(2) + twisted_bc(2)) / this%length(2))) &
                         + cos(pi * ((k_vec(1) + twisted_bc(1)) / this%length(1) &
                                     - (k_vec(2) + twisted_bc(2)) / this%length(2))))

    end function dispersion_rel_tilted

    function dispersion_rel_ole(this, k_vec) result(disp)
        class(ole) :: this
        integer, intent(in) :: k_vec(3)
        real(dp) :: disp
#ifdef DEBUG_
        character(*), parameter :: this_routine = "dispersion_rel_ole"
#endif

        ASSERT(this%is_periodic())

        ! i finally figured out how to do the non-orthogonal
        ! boundary conditions:
        ! although i am not sure about the twisted BC in this case!
        disp = 2.0_dp * (cos(2 * pi / real(sum(this%length(1:2)), dp) * &
                             ((k_vec(1) + twisted_bc(1)) - (k_vec(2) + twisted_bc(2)))) &
                         + cos(2 * pi / real(sum(this%length(1:2)), dp) * &
                               (this%length(2) / real(this%length(1), dp) * (k_vec(1) + twisted_bc(1)) &
                                + (k_vec(2) + twisted_bc(2)))))

    end function dispersion_rel_ole

    function dispersion_rel_not_implemented(this, k_vec) result(disp)
        class(lattice) :: this
        integer, intent(in) :: k_vec(3)
        real(dp) :: disp
        character(*), parameter :: this_routine = "dispersion_rel"

        unused_var(this)
        unused_var(k_vec)

        call stop_all(this_routine, &
                      "dispersion relation not yet implemented for this lattice type!")
#ifdef WARNING_WORKAROUND_
        disp = 0.0_dp
        unused_var(this)
        unused_var(k_vec)
#endif

    end function dispersion_rel_not_implemented

    function dispersion_rel_orb(this, orb) result(disp)
        class(lattice) :: this
        integer, intent(in) :: orb
        real(dp) :: disp

        integer :: k_vec(3)

        disp = this%dispersion_rel(this%get_k_vec(orb))

    end function dispersion_rel_orb

    function dispersion_rel_spin_orb(this, orb) result(disp)
        class(lattice) :: this
        integer, intent(in) :: orb
        real(dp) :: disp

        integer :: k_vec(3)

        disp = this%dispersion_rel(this%get_k_vec(get_spatial(orb)))

    end function dispersion_rel_spin_orb

    function dot_prod_not_implemented(this, k_vec, r_vec) result(dot)
        ! for the "fourier transform" implement the correct
        ! dot-product with all the factors of pi and n_sites implemented
        ! for each lattice
        class(lattice) :: this
        integer, intent(in) :: k_vec(3), r_vec(3)
        real(dp) :: dot
        character(*), parameter :: this_routine = "dot_prod_not_implemented"

        unused_var(this)
        unused_var(k_vec)
        unused_var(r_vec)

        call stop_all(this_routine, "not yet implemented for this lattice type!")
#ifdef WARNING_WORKAROUND_
        dot = 0.0_dp
        unused_var(this)
        unused_var(k_vec)
        unused_var(r_vec)
#endif

        dot = 0.0_dp

    end function dot_prod_not_implemented

    function dot_prod_chain(this, k_vec, r_vec) result(dot)
        class(chain) :: this
        integer, intent(in) :: k_vec(3), r_vec(3)
        real(dp) :: dot

        dot = 2.0_dp * PI / real(this%get_nsites(), dp) * &
              (k_vec(1) + twisted_bc(1)) * r_vec(1)

    end function dot_prod_chain

    function dot_prod_rect(this, k_vec, r_vec) result(dot)
        class(rectangle) :: this
        integer, intent(in) :: k_vec(3), r_vec(3)
        real(dp) :: dot

        dot = 2.0_dp * PI * ((k_vec(1) + twisted_bc(1)) * r_vec(1) / this%length(1) &
                             + (k_vec(2) + twisted_bc(2)) * r_vec(2) / this%length(2))

    end function dot_prod_rect

    function dot_prod_tilted(this, k_vec, r_vec) result(dot)
        class(tilted) :: this
        integer, intent(in) :: k_vec(3), r_vec(3)
        real(dp) :: dot

        dot = PI * (((k_vec(1) + twisted_bc(1)) / this%length(1) + &
                     (k_vec(2) + twisted_bc(2)) / this%length(2)) * r_vec(1) + &
                    ((k_vec(1) + twisted_bc(1)) / this%length(1) - &
                     (k_vec(2) + twisted_bc(2)) / this%length(2)) * r_vec(2))

    end function dot_prod_tilted

    function sort_unique(list) result(output)
        integer, intent(in) :: list(:)
        integer, allocatable :: output(:)

        integer :: i, min_val, max_val, unique(size(list))

        unique = 0
        i = 0
        min_val = minval(list) - 1
        max_val = maxval(list)

        do while (min_val < max_val)
            i = i + 1
            min_val = minval(list, mask=list > min_val)
            unique(i) = min_val
        end do
        allocate(output(i), source=unique(1:i))

    end function sort_unique

    subroutine init_aim(this, length_x, length_y)
        class(aim) :: this
        integer, intent(in) :: length_x, length_y
        character(*), parameter :: this_routine = "init_aim"

        unused_var(this)
        unused_var(length_x)
        unused_var(length_y)

    end subroutine init_aim

    subroutine init_lattice(this, length_x, length_y, length_z, &
                            t_periodic_x, t_periodic_y, t_periodic_z, t_bipartite_order)
        ! and write the first dummy initialize
        class(lattice) :: this
        integer, intent(in) :: length_x, length_y, length_z
        logical, intent(in) :: t_periodic_x, t_periodic_y, t_periodic_z
        logical, intent(in), optional :: t_bipartite_order
        character(*), parameter :: this_routine = "init_lattice"

        integer :: n_sites, i
        logical :: t_bipartite_order_
        def_default(t_bipartite_order_, t_bipartite_order, .false.)

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

            this%lat_vec(1, 1) = length_x
            ! the type specific routine deal with the check of the
            ! length!

            ! i should not call this set_nsites since this really should
            ! just imply that it set the variable
            ! introduce a second routine, which first determines the
            ! number of sites depending on the lattice type

            ! how should i define the lattice k_vectors..
            this%k_vec(1, 1) = length_x

            this%t_bipartite_order = t_bipartite_order_


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

            this%lat_vec(1, 1) = this%length(1)
            this%lat_vec(2, 2) = this%length(2)

            ! i also need to assign the lattice k-vectors..
            ! and i need to do it correctly..
            this%k_vec(1, 1) = this%length(1)
            this%k_vec(2, 2) = this%length(2)

            this%t_bipartite_order = t_bipartite_order_


        class is (tilted)

            call this%set_ndim(DIM_RECT)
            ! for the tilted we deal internally always with x as the
            ! lower of the two inputs. due to symmetry this does not
            ! make a difference
            ! and do not allow a 1xY or Yx1 lattice, since this implementation
            ! annoys me too much!
            if (length_x == 1 .or. length_y == 1) then
                call stop_all(this_routine, "incorrect size for tilted lattice!")
            end if
            call this%set_length(min(length_x, length_y), max(length_x, length_y))
            call this%set_nconnect_max(4)

            this%lat_vec(1:2, 1) = [this%length(1), this%length(1)]
            this%lat_vec(1:2, 2) = [-this%length(2), this%length(2)]

            this%k_vec(1:2, 1) = [this%length(1), this%length(1)]
            this%k_vec(1:2, 2) = [-this%length(2), this%length(2)]

            this%t_bipartite_order = t_bipartite_order_

        class is (ole)
            call this%set_ndim(DIM_RECT)

            if (length_x < 2 .or. length_y < 2 .or. length_x == length_y) then
                call stop_all(this_routine, "incorrect size for Oles Cluster")
            end if
            call this%set_length(min(length_x, length_y), max(length_x, length_y))
            call this%set_nconnect_max(4)

            this%lat_vec(1:2, 1) = [this%length(1), this%length(1)]
            this%lat_vec(1:2, 2) = [-this%length(2), this%length(1)]

            this%k_vec(1:2, 1) = [this%length(1), this%length(1)]
            this%k_vec(1:2, 2) = [-this%length(1), this%length(2)]

            if (t_bipartite_order_) then
                call stop_all(this_routine, &
                    "bipartite order not yet implemented for Ole lattice")
            end if

        class is (sujun)
            call this%set_ndim(DIM_RECT)

            if (length_x /= 1 .or. length_y /= 3) then
                call stop_all(this_routine, "incorrect size for Sujun cluster")
            end if

            call this%set_length(1,3)
            call this%set_nconnect_max(4)

            this%lat_vec(1:2, 1) = [1,3]
            this%lat_vec(1:2, 2) = [-3,1]

            ! k-vec todo..

        class is (ext_input)

            call read_lattice_struct(this)

        class is (cube)
            call this%set_ndim(DIM_CUBE)
            call this%set_length(length_x, length_y, length_z)
            call this%set_nconnect_max(6)

            this%lat_vec(1, 1) = this%length(1)
            this%lat_vec(2, 2) = this%length(2)
            this%lat_vec(3, 3) = this%length(3)

            this%k_vec(1, 1) = this%length(1)
            this%k_vec(2, 2) = this%length(2)
            this%k_vec(3, 3) = this%length(3)

            if (t_bipartite_order_) then
                call stop_all(this_routine, &
                    "bipartite order not yet implemented for cubic lattice")
            end if

        class is (triangular)
            call this%set_ndim(DIM_RECT)
            call this%set_length(length_x, length_y)
            ! for a filling with triangles the maximum connection is 6!

            call this%set_nconnect_max(6)


            ! todo: set lattice vector! and figure that out correctly!
            ! and write a more general routine to set the lattice
            ! vectors for all types of lattices!
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

            this%impurity_sites = [(i, i=1, length_x)]
            this%bath_sites = [(length_x + i, i=1, length_y)]

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

            this%impurity_sites = [(i, i=1, length_x)]
            this%bath_sites = [(length_x + i, i=1, length_y)]

        class default
            call stop_all(this_routine, "unexpected lattice type!")

        end select
        ! do i want to allocate sites here or in the initializer?
        ! well the specific site initializer will be different for all the
        ! types of lattices.. so best would be to do everything which is
        ! common to all routine here!
        call this%allocate_sites(n_sites)

        call this%initialize_sites()

        ! and fill the lookup table for the site index determination from k vectors
        if (t_k_space_hubbard .or. (t_trans_corr_hop .and. t_new_real_space_hubbard)) then
            call this%initialize_lu_table()
        end if

    end subroutine init_lattice

    subroutine init_hop_cache_bounds(this, r_min, r_max)
        ! initialize the lower and upper bounds of the cached
        ! hopping-transcorr factor, indexed by the involved r-vector
        ! maybe move this function to lattice_mod?!
        class(lattice) :: this
        integer, intent(out), optional :: r_min(3), r_max(3)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "init_hop_cache_bounds"
#endif
        integer :: n_sites, i, j, ri(3), rj(3), r_diff(3)

        if (.not. (all(this%r_min == this%r_max) .and. all(this%r_min == 0))) then

            if (present(r_min) .and. present(r_max)) then
                r_min = this%r_min
                r_max = this%r_max
            end if

            return
        end if

        n_sites = this%get_nsites()

        do i = 1, n_sites
            ri = this%get_r_vec(i)
            do j = 1, n_sites
                rj = this%get_r_vec(j)

                r_diff = ri - rj

                where (r_diff < this%r_min) this%r_min = r_diff
                where (r_diff > this%r_max) this%r_max = r_diff

            end do
        end do

        if (present(r_min) .and. present(r_max)) then
            r_min = this%r_min
            r_max = this%r_max
        end if

    end subroutine init_hop_cache_bounds

    subroutine initialize_lu_table(this)
        implicit none
        class(lattice) :: this

        ! first, get the dimension of the lookup tables
        call this%get_lu_table_size()
        call this%init_basis_vecs()
        ! and allocate the lookup tables accordingly
        allocate(this%lu_table(this%kmin(1):this%kmax(1), this%kmin(2):this%kmax(2) &
                                , this%kmin(3):this%kmax(3)))
        allocate(this%bz_table(this%kmin(1):this%kmax(1), this%kmin(2):this%kmax(2), &
                                this%kmin(3):this%kmax(3)))
        ! write(stdout, *) "Lookup table size is ", 2 * (this%kmax(1) - this%kmin(1) + 1) &
            ! * (this%kmax(2) - this%kmin(2) + 1) * 8 / 1024, " kB"
        ! and fill thee lookup table with the bz vectors
        call this%fill_bz_table()
        ! now, fill the lookup table with the indices of the states
        call this%fill_lu_table()

    end subroutine initialize_lu_table

    subroutine get_lu_table_size(this)
        implicit none
        class(lattice) :: this
        integer :: i, j, ki(sdim), kj(sdim), ka(sdim), nsites, m, a, k_sum(sdim)
        integer :: k, b, kk(sdim), kb(sdim)

        this%kmin = 0
        this%kmax = 0
        ! determine the maximum/minimum indices that appear in the lu table
        ! the lu table shall contain all reachable momenta ki + kj - ka
        nsites = this%get_nsites()
        ! for 2-body transcorrelation we need 5 momenta in total
        if (t_trans_corr_2body .and. t_k_space_hubbard) then
            do i = 1, nsites
                ki = this%get_k_vec(i)
                do j = 1, nsites
                    kj = this%get_k_vec(j)
                    do k = 1, nsites
                        kk = this%get_k_vec(k)
                        do a = 1, nsites
                            ka = this%get_k_vec(a)
                            do b = 1, nsites
                                kb = this%get_k_vec(b)

                                k_sum = ki + kj + kk - ka - kb

                                do m = 1, sdim
                                    if (k_sum(m) < this%kmin(m)) this%kmin(m) = k_sum(m)
                                    if (k_sum(m) > this%kmax(m)) this%kmax(m) = k_sum(m)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        else
            do i = 1, nsites
                ki = this%get_k_vec(i)
                do j = 1, nsites
                    kj = this%get_k_vec(j)
                    do a = 1, nsites
                        ka = this%get_k_vec(a)
                        k_sum = ki + kj - ka
                        do m = 1, sdim
                            if (k_sum(m) < this%kmin(m)) this%kmin(m) = k_sum(m)
                            if (k_sum(m) > this%kmax(m)) this%kmax(m) = k_sum(m)
                        end do
                    end do
                end do
            end do
        end if

    end subroutine get_lu_table_size

    subroutine fill_lu_table(this)
        implicit none
        class(lattice) :: this
        integer :: ki(sdim), kj(sdim), ka(sdim), i, j, m, a, k_check(sdim), k_sum(sdim), nsites
        integer :: k, b, kk(sdim), kb(sdim)

        nsites = this%get_nsites()
        !U.Ebling:
        !The older loop took a very long to finish for any lattice that is not super tiny.
        !I tried a 21x5x1 rectangle and it did not finish after 2 days
        !It over-counts a lot.
        !Below is my optimized version, which loops directly over momenta instead of orbitals
        !There is no need to distinguish the 1-body and 2-body transcorrelation terms, because it
        !uses the result of the subroutine get_lu_table_size
        do i=this%kmin(1),this%kmax(1)
            do j=this%kmin(2),this%kmax(2)
                do k=this%kmin(3),this%kmax(3)
                    k_sum(1)=i
                    k_sum(2)=j
                    k_sum(3)=k
                    k_check=this%map_k_vec(k_sum)
                    do m=1,nsites
                        if(all(k_check == this%get_k_vec(m))) then
                            this%lu_table(k_sum(1),k_sum(2),k_sum(3)) = m
                            exit
                        end if
                    end do
                end do
            end do
        end do

    end subroutine fill_lu_table

    subroutine fill_bz_table(this)
        implicit none
        class(lattice) :: this
        integer :: kx, ky, kz, k(sdim)

        ! here, we store for a bunch of important k vectors, whether they are in
        ! the first BZ or not

        this%bz_table = .false.
        ! check all vectors in the [kmin,kmax] range
        do kx = this%kmin(1), this%kmax(1)
            do ky = this%kmin(2), this%kmax(2)
                do kz = this%kmin(3), this%kmax(3)
                    k = (/kx, ky, kz/)
                    ! write down if k is in the BZ
                    ! the check is done by looping over all sites
                    this%bz_table(k(1), k(2), k(3)) = this%inside_bz_explicit(k)
                end do
            end do
        end do

    end subroutine fill_bz_table

    subroutine deallocate_caches(this)
        implicit none
        class(lattice) :: this

        deallocate(this%lu_table)
        deallocate(this%bz_table)
    end subroutine deallocate_caches

    function aim_lattice_constructor(lat_type, length_x, length_y) result(this)
        character(*), intent(in) :: lat_type
        integer, intent(in) :: length_x, length_y
        class(aim), pointer :: this
        character(*), parameter :: this_routine = "aim_lattice_constructor"

        select case (lat_type)
        case ('chain', 'aim-chain', 'chain-aim')

            allocate(aim_chain :: this)

        case ('star', 'aim-star', 'star-aim')

            allocate(aim_star :: this)

        case default
            ! stop here because a incorrect lattice type was given
            call stop_all(this_routine, &
                          'incorrect lattice type provided in lattice_constructor!')

        end select

        ! the initializer deals with the different types then..
        call this%initialize(length_x, length_y, 1, .false., .false., .false.)

    end function aim_lattice_constructor

    function lattice_constructor(lattice_type, length_x, length_y, length_z, t_periodic_x, &
                                 t_periodic_y, t_periodic_z, space, t_bipartite_order) result(this)
        ! write a general public lattice_constructor for lattices
        ! the number of inputs are still undecided.. do we always have
        ! the same number or differing number of inputs?
        ! i guess, since we will be using global variables which are either
        ! read in or set as default to init it..
        character(*), intent(in) :: lattice_type
        integer, intent(in) :: length_x, length_y, length_z
        logical, intent(in) :: t_periodic_x, t_periodic_y, t_periodic_z
        character(*), intent(in), optional :: space
        logical, intent(in), optional :: t_bipartite_order
        class(lattice), pointer :: this
        character(*), parameter :: this_routine = "lattice_constructor"

        select case (lattice_type)
        case ('chain')

            allocate(chain :: this)

        case ('star')

            allocate(star :: this)

        case ('square')
            ! i guess i want to make a seperate case for the tilted
            ! square, although just the boundary conditions change, but also
            ! the length input changes
            if (length_x /= length_y) then
                call stop_all(this_routine, &
                              "incorrect length input for square lattice!")

            end if

            allocate(rectangle :: this)

        case ('rectangle', 'rect')

            allocate(rectangle :: this)

        case ('tilted', 'tilted-square', 'square-tilted')

            allocate(tilted :: this)

        case ('cube', 'cubic')
            ! for the sake of no better name also use "cube" even if the
            ! sides are not the same length

            if (any([length_x, length_y, length_z] < 2)) then
                call stop_all(this_routine, &
                              "too short cube side lengths < 2!")
            end if

            allocate(cube :: this)

        case ('triangular', 'triangle', 'tri')

            if (any([length_x, length_y] < 2)) then
                call stop_all(this_routine, &
                              "too short lengths for triangular lattice! < 2!")
            end if

            allocate(triangular :: this)

        case ('hexagonal', 'hex', 'hexagon', 'honeycomb')

            allocate(hexagonal :: this)

        case ('kagome')
            allocate(kagome :: this)

        case ('ole')
            allocate(ole :: this)

        case ('sujun')
            allocate(sujun :: this)

        case ('ext_input')
            allocate(ext_input :: this)

        case default
            ! stop here because a incorrect lattice type was given
            call stop_all(this_routine, &
                          'incorrect lattice type provided in lattice_constructor!')

        end select

        call this%set_name(lattice_type)

        ! depending on the string input defining lattice type
        ! initialize corresponding lattice
        if (present(space)) then
            select case (space)
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
                             t_periodic_x, t_periodic_y, t_periodic_z, t_bipartite_order)

    end function lattice_constructor

    function site_constructor(ind, n_neighbors, neighbors, k_vec, r_vec, site_type) &
        result(this)
        ! similar to the lattice constructor i want to have a site
        ! constructor, which depending on some input constructs the
        ! specific sites on a lattice
        ! for now the index is the only necessary input..
        ! include more in this site constructor here
        integer, intent(in) :: ind
        integer, intent(in) :: n_neighbors
        integer, intent(in) :: neighbors(n_neighbors)
        integer, intent(in), optional :: k_vec(3)
        integer, intent(in), optional :: r_vec(3)
        character(*), intent(in), optional :: site_type
        ! i think i have to use pointers again..
        ! but maybe this is really bad to deal with in the rest of the code..
        type(site) :: this
        character(*), parameter :: this_routine = "site_constructor"

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

        if (present(k_vec)) then
            if (present(r_vec)) then
                call this%initialize(ind, n_neighbors, neighbors, k_vec, r_vec)
            else
                call this%initialize(ind, n_neighbors, neighbors, k_vec)
            end if
        else
            ASSERT(.not. present(r_vec))
            call this%initialize(ind, n_neighbors, neighbors)
        end if

    end function site_constructor

    subroutine init_site(this, ind, n_neighbors, neighbors, k_vec, r_vec)
        ! for now only use the index in the initalization
        class(site) :: this
        integer, intent(in) :: ind, n_neighbors
        integer, intent(in) :: neighbors(n_neighbors)
        integer, intent(in), optional :: k_vec(3)
        integer, intent(in), optional :: r_vec(3)

        ! for the beginning i do not need more than to set the index..
        ! independent of the type
        call this%set_index(ind)
        call this%set_num_neighbors(n_neighbors)
        call this%allocate_neighbors(n_neighbors)
        call this%set_neighbors(neighbors)

        if (present(k_vec)) then
            call this%set_k_vec(k_vec)
        end if

        if (present(r_vec)) then
            call this%set_r_vec(r_vec)
        end if

    end subroutine init_site

    subroutine set_k_vec(this, k_vec)
        class(site) :: this
        integer, intent(in) :: k_vec(3)

        this%k_vec = k_vec

    end subroutine set_k_vec

    subroutine set_r_vec(this, r_vec)
        class(site) :: this
        integer, intent(in) :: r_vec(3)

        this%r_vec = r_vec

    end subroutine set_r_vec

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

        nullify (this)

        select type (this)
        class is (aim)
            deallocate(this%impurity_sites)
            deallocate(this%bath_sites)
        end select

    end subroutine lattice_deconstructor

    subroutine aim_deconstructor(this)
        class(aim), pointer :: this

        call this%deallocate_sites()

        nullify (this)

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

    function get_num_neighbors_lattice(this, ind) result(n_neighbors)
        class(lattice) :: this
        integer, intent(in) :: ind
        integer :: n_neighbors
#ifdef DEBUG_
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
        integer, allocatable :: neighbors(:)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "get_neighbors_lattice"
#endif

        ! make all assert on a seperate line, so we exactly know what is
        ! going wrong..
        ASSERT(ind <= this%get_nsites())
        ASSERT(ind > 0)
        ASSERT(allocated(this%sites))
        ASSERT(allocated(this%sites(ind)%neighbors))

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
#ifdef DEBUG_
        character(*), parameter :: this_routine = "get_spinorb_neighbors_lat"
#endif

        ASSERT(spinorb <= 2 * this%get_nsites())
        ASSERT(spinorb > 0)
        ASSERT(allocated(this%sites))
        ASSERT(allocated(this%sites(get_spatial(spinorb))%neighbors))

        neighbors = this%get_neighbors(get_spatial(spinorb))

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

        unused_var(this)
        if (present(length_z)) then
            unused_var(length_z)
        end if

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
        unused_var(this)
        if (present(length_z)) then
            unused_var(length_z)
        end if

        unused_var(this)
        unused_var(length_z)

        if (max(length_x, length_y) < 1 .or. min(length_x, length_y) > 1 .or. &
            min(length_x, length_y) < 0) then
            n_sites = -1
            call stop_all(this_routine, "something went wrong in length input!")

        else
            n_sites = max(length_x, length_y)

        end if

    end function calc_nsites_star

    logical function is_periodic_aim_star(this, dimen)
        class(aim_star) :: this
        integer, intent(in), optional :: dimen
        unused_var(this)
        unused_var(dimen)
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

        unused_var(this)
        if (present(length_z)) then
            unused_var(length_z)
        end if

        if (max(length_x, length_y) < 1 .or. min(length_x, length_y) > 1 .or. &
            min(length_x, length_y) < 0) then
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

        unused_var(this)
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
        unused_var(this)
        if (present(length_z)) then
            unused_var(length_z)
        end if

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
        unused_var(this)
        if (present(length_z)) then
            unused_var(length_z)
        end if

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
        unused_var(this)
        if (present(length_z)) then
            unused_var(length_z)
        end if

        unused_var(this)
        unused_var(length_z)

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
        unused_var(this)
        if (present(length_z)) then
            unused_var(length_z)
        end if
        n_sites = 2 * length_x * length_y

    end function calc_nsites_tilted

    function calc_nsites_sujun(this, length_x, length_y, length_z) result(n_sites)
        class(sujun) :: this
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        integer :: n_sites
        character(*), parameter :: this_routine = "calc_nsites_sujun"
        unused_var(this)
        unused_var(length_x)
        unused_var(length_y)
        if (present(length_z)) then
            unused_var(length_z)
        end if

        n_sites = 10

    end function calc_nsites_sujun

    function calc_nsites_ole(this, length_x, length_y, length_z) result(n_sites)
        class(ole) :: this
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        integer :: n_sites
        character(*), parameter :: this_routine = "calc_nsites_ole"
        unused_var(this)
        if (present(length_z)) then
            unused_var(length_z)
        end if

        ! oles cluster we want to look at are defined by the vectors
        ! (x,x), (-y,x) and i also think  y = x + 2 is a requisite but i am
        ! not sure.. it is a non-orthogonal unit cell!
        ! check for the validity of the input
        if (length_x < 2 .or. length_y < 2 .or. length_x == length_y .or. &
            length_x > length_y) then
            ! and as a convention y > x is enforced!
            print *, "Lx: ", length_x
            print *, "Ly: ", length_y
            call stop_all(this_routine, "incorrect input for Ole Clusters")
        end if

        ! or shorter:
        n_sites = length_x * (length_x + length_y)

    end function calc_nsites_ole

    ! Setter and getter routines for the private data of the lattice types

    subroutine set_length_lattice(this, length_x, length_y)
        class(lattice) :: this
        integer, intent(in) :: length_x, length_y
        character(*), parameter :: this_routine = "set_length_lattice"

        unused_var(length_x)
        unused_var(length_y)
        unused_var(this)

        call stop_all(this_routine, &
                      'type(lattice) should never be actually instantiated!')

    end subroutine set_length_lattice

    subroutine set_length_aim_star(this, length_x, length_y, length_z)
        class(aim_star) :: this
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        character(*), parameter :: this_routine = "set_length_aim_star"

        unused_var(length_x)
        unused_var(length_y)
        unused_var(length_z)
        unused_var(this)

        ! actually the length of a start is not really defined..
        ! maybe i should rethink if i make this part of the
        ! original lattice class then..
        call stop_all(this_routine, "length not defined for 'star' geometry!")
        unused_var(length_x)
        unused_var(length_y)
        if (present(length_z)) then
            unused_var(length_z)
        end if
        unused_var(this)
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

        unused_var(length_x)
        unused_var(length_y)
        if (present(length_z)) then
            unused_var(length_z)
        end if
        unused_var(this)
    end subroutine set_length_star

    subroutine set_length_chain(this, length_x, length_y, length_z)
        class(chain) :: this
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        character(*), parameter :: this_routine = "set_length_chain"
        if (present(length_z)) then
            unused_var(length_z)
        end if

        unused_var(length_z)

        ! the input checkin is all done in the calc_nsites routine!
        this%length = this%calc_nsites(length_x, length_y)

    end subroutine set_length_chain

    subroutine set_length_cube(this, length_x, length_y, length_z)
        class(cube) :: this
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
#ifdef DEBUG_
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
        if (present(length_z)) then
            unused_var(length_z)
        end if

        unused_var(length_z)

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
    integer pure function get_length_star(this, dimen)
        class(star), intent(in) :: this
        integer, intent(in), optional:: dimen
        unused_var(dimen)
        unused_var(this)

        get_length_star = STAR_LENGTH

    end function get_length_star

    integer pure function get_length_chain(this, dimen)
        class(chain), intent(in) :: this
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

    integer pure function get_length_cube(this, dimen)
        class(cube), intent(in) :: this
        integer, intent(in), optional :: dimen
#ifdef DEBUG_
        character(*), parameter :: this_routine = "get_length_cube"
#endif

        ASSERT(present(dimen))
        ASSERT(dimen > 0)
        ASSERT(dimen <= 3)

        get_length_cube = this%length(dimen)

    end function get_length_cube

    integer pure function get_length_rect(this, dimen)
        class(rectangle), intent(in) :: this
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

    integer pure function get_length_aim_chain(this, dimen)
        class(aim_chain), intent(in) :: this
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
        unused_var(length_x)
        if (present(length_z)) then
            unused_var(length_z)
        end if

        unused_var(length_x)
        unused_var(length_z)

        ! as a definition make the length, even for multiple impurity chains
        ! as bath_sites + 1
        this%length = length_y + 1

    end subroutine set_length_aim_chain

    logical function is_periodic_star(this, dimen)
        ! this is always false.. the star geometry can't be periodic
        class(star) :: this
        integer, intent(in), optional :: dimen
        unused_var(dimen)
        unused_var(this)

        unused_var(dimen)
        unused_var(this)

        is_periodic_star = .false.

    end function is_periodic_star

    logical function is_periodic_chain(this, dimen)
        class(chain) :: this
        integer, intent(in), optional :: dimen
        unused_var(dimen)
        ! we do not want to deal with two dimensional flags for chains or?
        is_periodic_chain = this%t_periodic(1)
        ! the chain is only treated as periodic if both the flags are set
        ! to be periodic!

    end function is_periodic_chain

    logical function is_periodic_cube(this, dimen)
        class(cube) :: this
        integer, intent(in), optional :: dimen
#ifdef DEBUG_
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

    logical function is_periodic_x(this)
        class(lattice) :: this

        is_periodic_x = this%t_periodic_x

    end function is_periodic_x

    logical function is_periodic_y(this)
        class(lattice) :: this

        is_periodic_y = this%t_periodic_y

    end function is_periodic_y

    integer elemental function get_ndim(this)
        class(lattice), intent(in) :: this
        get_ndim = this%n_dim
    end function get_ndim

    integer elemental function get_nsites(this)
        class(lattice), intent(in) :: this
        get_nsites = this%n_sites
    end function get_nsites

    subroutine print_lat(this)
        class(lattice) :: this

        ! depending on the type print specific lattice information
        select type (this)
        class is (lattice)

            call stop_all("lattice%print()", &
                          "lattice type should never be directly instantiated!")

        class is (chain)

            print *, "Lattice type is: 'chain' "
            print *, " number of dimensions: ", this%get_ndim()
            print *, " max-number of neighbors: ", this%get_nconnect_max()
            print *, " number of sites of chain: ", this%get_nsites()
            print *, " is the chain periodic: ", this%is_periodic()

        class is (rectangle)

            if (trim(this%get_name()) == 'tilted') then
                print *, "Lattice type is 'tilted' "
            else
                print *, "Lattice type is 'rectangle' "
            end if

            print *, " number of dimensions: ", this%get_ndim()
            print *, " max-number of neighbors: ", this%get_nconnect_max()
            print *, " number of sites in the lattice: ", this%get_nsites()
            if (this%is_periodic()) then
                print *, " is the lattice periodic: True"
            end if
            print *, "primitive lattice vectors: (", this%lat_vec(1:2, 1), "), (", this%lat_vec(1:2, 2), ")"

        class is (ole)
            print *, "Lattice type is 'Ole' ;) "
            print *, "number of dimensions: ", this%get_ndim()
            print *, "max-number of neigbors: ", this%get_nconnect_max()
            print *, "number of sites: ", this%get_nsites()
            print *, "primitive lattice vectors: (", this%lat_vec(1:2, 1), "), (", this%lat_vec(1:2, 2), ")"
            print *, "TODO: more and better output! "

        class is (ext_input)
            print *, "Lattice read from lattice.file file "
            print *, "Lattice type is :",this%get_name()
            print *, "number of sites in the lattice: ", this%get_nsites()
            print *, "max-number of neigbors: ", this%get_nconnect_max()
            if (this%is_periodic()) then
                print *, " is the lattice periodic: True"
            end if

        end select

    end subroutine print_lat

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
            p_hole = 2.0_dp / real(nbasis - nel, dp)
        else
            p_hole = 1.0_dp / real(lat%get_nconnect_max(), dp)
        end if

        if (t_new_real_space_hubbard) then

            mat_ele = real(abs(bhub), dp)

        else if (t_tJ_model) then

            ! for the t-J take the maximum of hopping or exchange
            mat_ele = real(max(abs(bhub), abs(exchange_j)), dp)

        else if (t_heisenberg_model) then

            mat_ele = real(abs(exchange_j), dp)

        else if (t_k_space_hubbard) then
            mat_ele = abs(real(uhub, dp) / real(omega, dp))

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

    function read_lattice_n_sites(this, length_x, length_y, length_z) result(n_sites)

        class(ext_input):: this

        integer :: iunit, ios
        integer, intent(in) :: length_x, length_y
        integer, intent(in), optional :: length_z
        integer :: n_sites

        logical :: exists, leof

        CHARACTER(len=3) :: fmat
        CHARACTER(LEN=100) w
        character(*), parameter :: this_routine = "calc_nsites_gen"

        type(ManagingFileReader_t) :: file_reader
        type(TokenIterator_t) :: tokens

        unused_var(this)
        unused_var(length_x)
        unused_var(length_y)

        if (present(length_z)) then
            unused_var(length_z)
        end if

        file_reader = ManagingFileReader_t("lattice.file")

        lat: do while (file_reader%nextline(tokens))
            w = to_upper(tokens%next())
            select case (w)
            case ('N_SITES')
                n_sites = to_int(tokens%next())
            end select
        end do lat

    end function read_lattice_n_sites


    subroutine read_lattice_struct(this)

        class(ext_input):: this

        integer :: iunit, lat_dim, x1,y1, x2,y2, z1,z2, length_x, length_y, n_sites, ios

        logical :: exists, leof

        CHARACTER(len=3) :: fmat
        CHARACTER(LEN=100) w
        CHARACTER(LEN=NAME_LEN) lat_typ
        type(ManagingFileReader_t) :: file_reader
        type(TokenIterator_t) :: tokens

        file_reader = ManagingFileReader_t("lattice.file")

        lat: do while (file_reader%nextline(tokens))
            w = to_upper(tokens%next())

            select case (w)
            case ('DIM')
                call this%set_ndim(to_int(tokens%next()))

            case ('LATTICE_TYPE')
                this%name = to_upper(tokens%next())

            case ('LATTICE_PARAM')
                x1 = to_int(tokens%next())
                y1 = to_int(tokens%next())
                x2 = to_int(tokens%next())
                y2 = to_int(tokens%next())

                this%lat_vec(1:2, 1) = [x1, y1]
                this%lat_vec(1:2, 2) = [x2, y2]

            case ('N_CONNECT_MAX')
                this%n_connect_max = to_int(tokens%next())
            end select
        end do lat

        call this%set_length(1, 3)

    end subroutine read_lattice_struct


end module lattice_mod
