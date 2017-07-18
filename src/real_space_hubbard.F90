#include "macros.h"

! create a new module which covers all the real-space hubbard stuff in a more 
! clean and more efficient way. 
! I also want to implement triangular and kagome lattice types, so it is 
! good to get this better sorted. 
! also i definetly have to improve on the hubbard excitation generator
! i also probably should decouple all the UEG k-space hubbard at a certain 
! point. because it makes things super messy and also more inefficient how it 
! is done right now..
! and for the beginning i should probably stick to the already working 
! lattices in NECI
! and i also want to make it compatible to restart old runs with this new 
! implementation

module real_space_hubbard

    use SystemData, only: t_new_real_space_hubbard, lattice_type, length_x, &
                          length_y
    implicit none 


    ! some general needed quantities 
    integer :: n_dim, n_sites, n_connect_max

    ! i also want to have a linked list of neighbors for each site
    ! need a new type for that 
    type lattice_site
        ! what do i want to store here? 
        ! this could be optimized really neatly with bit-storing and 
        ! packing.. but do not optimize prematurely! 
        ! essentially i need a vector of vectors (of integers) 
        ! i will create this type nbasis times in a vector of types! 
        ! but also store the number of neighbors for each of these
        integer :: n_next_neighbors = 0
        integer, allocatable :: next_neighbor_list(:) 

        ! and it would be easy to generalize this for next-nearest neigbors 
        ! and stuff 
        ! and maybe for the AIM i could store flags here if this is a 
        ! impurity site or not.. 
        ! maybe i want to rename this type "lattice" 
        ! but do not implement unnecessary stuff for now! 

    end type lattice_site

    type lattice 
        ! create a type which contains general information about the 
        ! lattice structure 
        logical :: t_init = .false. 
        integer :: num_sites = 0 

        type(lattice_site), allocatable :: site(:)

    end type lattice

    ! and i guess i want to have a global lattice type.. to make it 
    ! easier accessible.. 
    type(lattice) :: global_lattice

contains 

    ! some brainstorming: 
    ! i want to change the input so that one specifies the lattice type 
    ! in the System input by 
    ! system hubbard [lattice-type] 

    subroutine init_real_space_hubbard() 
        ! routine, which does all of the necessary initialisation
        character(*), parameter :: this_routine = "init_real_space_hubbard"

        ! especially do some stop_all here so no wrong input is used 

    end subroutine init_real_space_hubbard

    ! then i have to think of how to set up the lattice.. 
    subroutine init_lattice()
        ! routine which sets up the lattice, like TMAT of nearest neighbors
        ! and creating the indexing of the neighbors of each site 
        ! lets break the convention of using a million of global 
        ! variables in neci.. and try to start using more and more 
        ! explicit variables, or atleast enable optional input 
        ! variables to unit-test the function more easily
        ! although i just realized that this bloats this whole 
        ! function way too much, with 3 variables for each input:
        ! a optional input one, a used on in this routine and 
        ! the global one which is used if no input is provided.. 
        ! argh.. fortran gets annoying..
        character(*), parameter :: this_routine = "init_lattice"

        ! what are the possible ones:
        ! the ones already in NECI (do them first!)
        ! CHAIN: just needs number of sites or length and if open-bc 
        ! SQUARE: need L_x and L_y and the boundary condition
        ! TILTED: do it as the input is done in the old implo, but maybe 
        !           enable L_x /= L_y 
        ! 
        ! after i check those above also implement: 
        ! TRIANGULAR: stick to equal sided triangles here and figure out 
        !               the boundary conditions 
        ! KAGOME: probably also stick to one length parameter and also 
        !           think about the BC

        select case (lattice_type) 
        case ("chain", "CHAIN")
            ! determine some basic values
            n_dim = 1
            n_sites = length_x
            ! the maximum connection is necessary for the time-step
            n_connect_max = 2

        case ("square", "SQUARE")
            n_dim = 2
            n_sites = length_x * length_y
            n_connect_max = 4

        case ("tilted", "TILTED")
            n_dim = 2 
            n_sites = 2 * length_x * length_y
            n_connect_max = 4

        case ("triangular", "TRIANGULAR") 
            n_dim = 2
!             n_sites = todo
            n_connect_max = 3

        case ("kagome", "KAGOME")
            n_dim = 2
!             n_sites = todo
            n_connect_max = 4

        case default 
            print *, "Incorrect lattice type provided! Choose:"
            print *, "chain, square, tilted, triangular, kagome"
            call stop_all(this_routine, "Incorrect lattice type!")
        end select

        ! i also have to set a lot of the other necessary stuff.. 
        ! but maybe i can avoid that by still using a lot of the 
        ! already provided initialisation and only set up the additional 
        ! necessary stuff here.

        ! i have to setup up the set_nearest neighbor routine
        if (global_lattice%t_init) then 
            ! if it has already been initialized before deinit all the 
            ! components 
        end if

        call create_neighbor_list(n_sites, 
    end subroutine init_lattice
    
    subroutine init_tmat() 
        ! i should create a new, more flexible routine which sets up the 
        ! TMAT for the different lattice types. although i am not sure if 
        ! we need this anymore 
        ! this also depends on the boundary conditions
        character(*), parameter :: this_routine = "init_tmat"

    end subroutine init_tmat

    subroutine create_neighbor_list(n_sit, lat_type, len_x, len_y, t_periodic)
        ! i want to index the neigbors of each lattice site with a linked 
        ! list. this will be helpful especially for the kagome lattice 
        ! where sites can have a differing amount of neighbors 
        ! and with open boundary conditions this also happens..
        integer, intent(in) :: n_sit
        character(*), intent(in) :: lat_type 
        integer, intent(in) :: len_x, len_y 
        logical, intent(in) :: t_periodic

        character(*), parameter :: this_routine = "create_neighbor_list"
    
        ! i have to get the data_type running 
        if (allocated(global_lattice)) deallocate(global_lattice)

    end subroutine create_neighbor_list

    subroutine determine_optimal_time_step()
        ! move this time-step determination to this routine for the real
        ! space hubbard to have it fully conained
        character(*), parameter :: this_routine = "determine_optimal_time_step"

    end subroutine determine_optimal_time_step

    subroutine real_space_excit_gen()
        ! i also want a new optimized real-space excitation generator
        character(*), parameter :: this_routine = "real_space_excit_gen"

    end subroutine real_space_excit_gen

    ! also optimize the matrix element calculation
    subroutine calc_diag_mat_ele_rsh() 
        ! the diagonal matrix element is essentialy just the number of 
        ! doubly occupied sites times U
        character(*), parameter :: this_routine = "calc_diag_mat_ele_rsh"

    end subroutine calc_diag_mat_ele_rsh

    subroutine calc_off_diag_mat_ele_rsh()
        ! in case we need it, the off-diagonal, except parity is just 
        ! -t if the hop is possible
        character(*), parameter :: this_routine = "calc_off_diag_mat_ele_rsh"

    end subroutine calc_off_diag_mat_ele_rsh
    
    ! what else?
    subroutine create_neel_state() 
        ! probably a good idea to have a routine which creates a neel state
        ! (or one of them if it is ambigous)
        character(*), parameter :: this_routine = "create_neel_state"

    end subroutine create_neel_state

end module real_space_hubbard
    
    
