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
                          length_y, uhub, nbasis, bhub, t_open_bc_x, &
                          t_open_bc_y, G1, ecore, nel, nOccAlpha, nOccBeta
    use lattice_mod, only: lattice
    use constants, only: dp
    use procedure_pointers, only: get_umat_el, generate_excitation
    use OneEInts, only: tmat2d
    use fcimcdata, only: pSingles, pDoubles, tsearchtau, tsearchtauoption
    use CalcData, only: t_hist_tau_search, t_hist_tau_search_option, tau
    use procedure_pointers, only: generate_excitation
    use tau_search, only: max_death_cpt
    use umatcache, only: gtid
    use dsfmt_interface, only: genrand_real2_dsfmt

    implicit none 

! and this is the global lattice class
    class(lattice), pointer :: lat

contains 

    ! some brainstorming: 
    ! i want to change the input so that one specifies the lattice type 
    ! in the System input by 
    ! system hubbard [lattice-type] 

    subroutine init_real_space_hubbard() 
        ! routine, which does all of the necessary initialisation
        character(*), parameter :: this_routine = "init_real_space_hubbard"
        ! especially do some stop_all here so no wrong input is used 

        ! first assert all the right input! 
        call check_real_space_hubbard_input() 

        ! which stuff do i need to initialize here? 
        get_umat_el => get_umat_el_hub

        ! i have to check if the lattice should be constructed from an fcidump 
        ! or created internally.. 
        if (trim(adjustl(lattice_type)) == 'read') then 
            ! then i have to construct tmat first 
            call init_tmat() 
            ! and then construct the lattice 
            lat => lattice(lattice_type, length_x, length_y, .not. t_open_bc_x, &
                .not. t_open_bc_y)
        else 
            ! otherwise i have to do it the other way around 
            lat => lattice(lattice_type, length_x, length_y, .not. t_open_bc_x, &
                .not. t_open_bc_y)
            call init_tmat(lat)

        end if
        
        ! i guess i have to setup G1 also.. argh.. i hate this! 
        allocate(G1(nbasis)) 
        G1(1:nbasis-1:2)%ms = -1
        G1(2:nbasis:2)%ms = 1


        ! Ecore should default to 0, but be sure anyway! 
        ecore = 0.0_dp

        ! and i have to point to the new hubbard excitation generator
        pSingles = 1.0_dp 
        pDoubles = 0.0_dp

        ! and i have to calculate the optimal time-step for the hubbard models. 
        ! where i need the connectivity of the lattice i guess? 
        generate_excitation => gen_excit_rs_hubbard
        
        ! i have to calculate the optimal time-step
        tau = determine_optimal_time_step()

        ! and i have to turn off the time-step search for the hubbard 
        tsearchtau = .false.
        ! set tsearchtauoption to true to use the death-tau search option
        tsearchtauoption = .true.

        t_hist_tau_search = .false. 
        t_hist_tau_search_option = .false. 
        
        max_death_cpt = 0.0_dp

    end subroutine init_real_space_hubbard

    subroutine check_real_space_hubbard_input() 
        use SystemData, only: tExch, tCSF, tReltvy, tUEG, tUEG2, tHub, & 
                              tKPntSym, tLatticeGens, tUEGNewGenerator, &
                tGenHelWeighted, tGen_4ind_weighted, tGen_4ind_reverse, &
                tUEGNewGenerator, tGen_4ind_part_exact, tGen_4ind_lin_exact, &
                tGen_4ind_2, tGen_4ind_2_symmetric, tGen_4ind_unbound, tStoreSpinOrbs, &
                tReal
        use umatcache, only : tTransGTid
        use OneEInts, only: tcpmdsymtmat, tOneelecdiag

        character(*), parameter :: this_routine = "check_real_space_hubbard_input"
        ! do all the input checking here, so no wrong input is used!

        if (tExch)            call stop_all(this_routine, "tExch set to true!")
        if (tCSF)             call stop_all(this_routine, "tCSF set to true!")
        if (tReltvy)          call stop_all(this_routine, "tReltvy set to true!")

        ! what else.. 
        if (tUEG)             call stop_all(this_routine, "tUEG set to true!")
        if (tUEG2)            call stop_all(this_routine, "tUEG2 set to true!")
        if (tHub)             call stop_all(this_routine, "tHub set to true!")
        if (tReal)            call stop_all(this_routine, "tReal set to true!")
        if (tKPntSym)         call stop_all(this_routine, "tKPntSym set to true!")
        if (tLatticeGens)     call stop_all(this_routine, "tLatticeGens set to true!")
        if (tUEGNewGenerator) call stop_all(this_routine, "tUEGNewGenerator set to true!")
        if (tGenHelWeighted)  call stop_all(this_routine, "tGenHelWeighted") 
        if (tGen_4ind_weighted) call stop_all(this_routine, "tGen_4ind_weighted") 
        if (tGen_4ind_reverse) call stop_all(this_routine, "tGen_4ind_reverse") 
        if (tGen_4ind_part_exact) call stop_all(this_routine, "tGen_4ind_part_exact") 
        if (tGen_4ind_2)        call stop_all(this_routine, "tGen_4ind_2") 
        if (tGen_4ind_2_symmetric) call stop_all(this_routine, "tGen_4ind_2_symmetric") 
        if (tGen_4ind_unbound)      call stop_all(this_routine, "tGen_4ind_unbound")
        if (tStoreSpinOrbs)     call stop_all(this_routine, "tStoreSpinOrbs")
        if (tTransGTid)         call stop_all(this_routine, "tTransGTid")
        if (tcpmdsymtmat)        call stop_all(this_routine, "tcpmdsymmat")
        if (tOneelecdiag)       call stop_all(this_routine, "tOneelecdiag")
            
    end subroutine check_real_space_hubbard_input

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

        class(lattice), pointer :: lat

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

        ! and especially for the Anderson Impurity models i need 
        ! Anderson-chain and Anderson-star

    end subroutine init_lattice
    
    subroutine init_tmat(lat)
        class(lattice), optional :: lat

        ! i should create a new, more flexible routine which sets up the 
        ! TMAT for the different lattice types. although i am not sure if 
        ! we need this anymore 
        ! this also depends on the boundary conditions
        character(*), parameter :: this_routine = "init_tmat"

        integer :: i, ind
        ! depending on the input i either create tmat2d here or is have to 
        ! set it up, so it can be used to create the lattice.. 
        ! but for the beginning i think i want to set it up here from the 
        ! lattice structure! 
        ! if the lattice class is alread set up and initialized, this indicates
        ! that it was build-in created and tmat has to be calculated from it 
        ! now!
        if (present(lat)) then 
            ! what do i need to do? 
            ! loop over the indices in the lattice and get the neighbors
            ! and i have to store it in spin-indices remember that!
            if (associated(tmat2d)) deallocate(tmat2d)
            allocate(tmat2d(nbasis,nbasis))
            tmat2d = 0.0_dp

            do i = 1, lat%get_nsites() 
                ind = lat%get_site_index(i)
                associate(next => lat%get_neighbors(i))
                    ! beta orbitals:
                    tmat2d(2*ind - 1, 2*next - 1) = bhub 
                    ! alpha: 
                    tmat2d(2*ind, 2*next) = bhub
                    
                    ASSERT(all(next > 0))
                    ASSERT(all(next <= nbasis/2))
                end associate
                ASSERT(lat%get_nsites() == nbasis/2)
                ASSERT(ind > 0) 
                ASSERT(ind <= nbasis/2)
                
            end do

        else
            ! this indicates that tmat has to be created from an fcidump 
            ! and the lattice is set up afterwards!
        end if

    end subroutine init_tmat

    function determine_optimal_time_step(time_step_death) result(time_step)
        real(dp), optional, intent(out) :: time_step_death
        real(dp) :: time_step
        ! move this time-step determination to this routine for the real
        ! space hubbard to have it fully conained
        character(*), parameter :: this_routine = "determine_optimal_time_step"

        real(dp) :: p_elec, p_hole, mat_ele, max_diag

        ! determine the optimal hubbard time-step for an optimized 
        ! hubbard excitation generation 
        ! the first electron is chosen at random 
        p_elec = 1.0_dp / real(nel, dp)

        ! and for a picked electron the possible neighbors are looked for 
        ! open holes so the lowest probability is determined by the 
        ! maximum numbers of connections 
        p_hole = 1.0_dp / real(lat%get_nconnect_max(), dp) 

        ! the matrix element is always just |t| 
        mat_ele = real(abs(bhub), dp)

        ! so the time-step is 
        time_step = p_elec * p_hole / mat_ele

        if (present(time_step_death)) then 
            ! the maximum death contribution is max(H_ii) - shift 
            ! and the maximum diagonal matrix element we know it is 
            ! the maximum possible number of doubly occupied sites times U 
            ! but this might be a too hard limit, since especially for high U 
            ! such an amount of double occupations is probably never reached 
            ! and the shift must also be included in the calculation.. 
            max_diag = real(abs(uhub) * min(nOccAlpha, nOccBeta), dp)

            time_step_death = 1.0_dp / max_diag
        end if

    end function determine_optimal_time_step

    !
    ! Generic excitaiton generator
    subroutine gen_excit_rs_hubbard (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)

        use SystemData, only: nel
        use bit_rep_data, only: NIfTot
        use FciMCData, only: excit_gen_store_type
        use constants, only: n_int, dp
        implicit none

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: run

        character(*), parameter :: this_routine = "gen_excit_rs_hubbard"

        integer :: iunused, ind , elec, id

        ! i only have single excitations in the hubbard model 
        ! the first plan is to choose an electron at random 
        ind = 1 + int(genrand_real2_dsfmt() * nel) 

        ! and then from the neighbors of this electron we pick an empty 
        ! spinorbital randomly, since all have the same matrix element 
        elec = nI(ind) 

        ! get the spatial index 
        id = gtid(elec) 


    end subroutine gen_excit_rs_hubbard

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

    function get_umat_el_hub(i,j,k,l) result(hel)
        integer, intent(in) :: i, j, k, l
        HElement_t(dp) :: hel 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_umat_el_hub"
#endif
        if (i == j .and. i == k .and. i == l) then 
            hel = h_cast(uhub)
        else 
            hel = h_cast(0.0_dp)
        end if

        ASSERT(i > 0)
        ASSERT(i <= nbasis/2)
        ASSERT(j > 0) 
        ASSERT(j <= nbasis/2)
        ASSERT(k > 0) 
        ASSERT(k <= nbasis/2)
        ASSERT(l > 0) 
        ASSERT(l <= nbasis/2)

    end function get_umat_el_hub

end module real_space_hubbard
    
    
