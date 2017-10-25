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
                          length_y, length_z, uhub, nbasis, bhub, t_open_bc_x, &
                          t_open_bc_y, t_open_bc_z, G1, ecore, nel, nOccAlpha, nOccBeta
    use lattice_mod, only: lattice
    use constants, only: dp, EPS, n_int
    use procedure_pointers, only: get_umat_el, generate_excitation
    use OneEInts, only: tmat2d, GetTMatEl
    use fcimcdata, only: pSingles, pDoubles, tsearchtau, tsearchtauoption
    use CalcData, only: t_hist_tau_search, t_hist_tau_search_option, tau
    use umatcache, only: gtid
    use dsfmt_interface, only: genrand_real2_dsfmt
    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet
    use bit_rep_data, only: NIfTot

    implicit none 

! and this is the global lattice class
    class(lattice), pointer :: lat

    real(dp) :: lat_tau_factor = 0.5_dp
    
    ! honjuns idea with the transcorrelated Hamiltonian we have a modified 
    ! hopping term: 
    ! t_ij^s = t exp[K(n_j^s' - n_i^s')] 
    ! so we need this K as an input parameter 
    real(dp) :: trans_corr_param = 1.0_dp 
    ! and a flag to start it 
    logical :: t_trans_corr = .false. 
    ! as one can see this modification is dependent on the current 
    ! occupation of the involved hopping orbitals! so it is not just a 
    ! change in the Hamiltonian 

    interface get_helement_rs_hub
        module procedure get_helement_rs_hub_ex_mat
        module procedure get_helement_rs_hub_general
    end interface get_helement_rs_hub

contains 

    ! some brainstorming: 
    ! i want to change the input so that one specifies the lattice type 
    ! in the System input by 
    ! system hubbard [lattice-type] 

    subroutine init_real_space_hubbard() 
        use SystemData, only: tExch, thub, treal
        ! routine, which does all of the necessary initialisation
        character(*), parameter :: this_routine = "init_real_space_hubbard"
        ! especially do some stop_all here so no wrong input is used 
        real(dp) :: tau_opt

        print *, "using new real-space hubbard implementation: "

        ! i do not need exchange integrals in the real-space hubbard model
        tExch = .false.
        ! after the whole setup i can set thub to false or? 
        thub = .false.
        ! and treal i can also set to false or? 
        treal = .false.
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
            lat => lattice(lattice_type, length_x, length_y, length_z, .not. t_open_bc_x, &
                .not. t_open_bc_y, .not. t_open_bc_z)
        else 
            ! otherwise i have to do it the other way around 
            lat => lattice(lattice_type, length_x, length_y, length_z, .not. t_open_bc_x, &
                .not. t_open_bc_y, .not. t_open_bc_z)
     
            ! if nbaiss was not yet provided:
            if (nbasis <= 0) then 
                nbasis = 2 * lat%get_nsites() 
            end if

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
        ! and maybe i have to be a bit more safe here and not be too near to 
        ! the optimal time-step
        tau_opt = determine_optimal_time_step()
        if (tau < EPS) then 
            print *, "setting time-step to optimally determined time-step: ", tau_opt
            print *, "times: ", lat_tau_factor
            tau = lat_tau_factor * tau_opt

        else 
            print *, "optimal time-step would be: ", tau_opt
            print *, "but tau specified in input!"
        end if

        ! and i have to turn off the time-step search for the hubbard 
        tsearchtau = .false.
        ! set tsearchtauoption to true to use the death-tau search option
        tsearchtauoption = .true.

        t_hist_tau_search = .false. 
        t_hist_tau_search_option = .false. 
        
        ! do not set that here, due to circular dependencies
!         max_death_cpt = 0.0_dp

    end subroutine init_real_space_hubbard

    subroutine check_real_space_hubbard_input() 
        use SystemData, only: tCSF, tReltvy, tUEG, tUEG2, tHub, & 
                              tKPntSym, tLatticeGens, tUEGNewGenerator, &
                tGenHelWeighted, tGen_4ind_weighted, tGen_4ind_reverse, &
                tUEGNewGenerator, tGen_4ind_part_exact, tGen_4ind_lin_exact, &
                tGen_4ind_2, tGen_4ind_2_symmetric, tGen_4ind_unbound, tStoreSpinOrbs, &
                tReal, tHPHF
        use umatcache, only : tTransGTid
        use OneEInts, only: tcpmdsymtmat, tOneelecdiag

        character(*), parameter :: this_routine = "check_real_space_hubbard_input"
        ! do all the input checking here, so no wrong input is used!

        if (thphf) call stop_all(this_routine, "hphf not yet implemented!")

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
        if (abs(bhub) < EPS)         call stop_all(this_routine, "bhub == 0!")
            
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
        use constants, only: n_int, dp, bits_n_int
        use get_excit, only: make_single
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

        integer :: iunused, ind , elec, id, src, orb, n_avail, n_orbs, i
        integer, allocatable :: neighbors(:), orbs(:)

        iunused = exflag; 

        ASSERT(associated(lat))

        ic = 1
        ! i only have single excitations in the hubbard model 
        ! the first plan is to choose an electron at random 
        elec = 1 + int(genrand_real2_dsfmt() * nel) 

        ! and then from the neighbors of this electron we pick an empty 
        ! spinorbital randomly, since all have the same matrix element 
        src = nI(elec) 

        ! get the spatial index 
        id = gtid(src) 

        ! now get neighbors
        n_orbs = lat%get_num_neighbors(id)
        allocate(neighbors(n_orbs)) 
        allocate(orbs(n_orbs))
        neighbors = lat%get_neighbors(id) 

        if (is_beta(src)) then 
            neighbors = 2 * neighbors - 1
        else 
            neighbors = 2 * neighbors
        end if

        n_avail = 0 
        orbs = -1 

        ! check which neighbors are empty 
        do i = 1, n_orbs 
            if (IsNotOcc(ilutI,neighbors(i))) then 
                n_avail = n_avail + 1 
                orbs(n_avail) = neighbors(i)
            end if
        end do

        if (n_avail == 0) then 
            ! no possible hoppings! so abort 
            nJ(1) = 0 
            pgen = 0.0_dp 
            return 
        end if

        ! otherwise pick a random orbital 
        ind = 1 + int(genrand_real2_dsfmt() * n_avail) 
        orb = orbs(ind) 
        pgen = 1.0_dp / real(n_avail * nel, dp) 

        call make_single(nI, nJ, elec, orb, ex, tParity) 

        ilutJ = ilutI 
        clr_orb(ilutJ, src)
        set_orb(ilutJ, orb)

    end subroutine gen_excit_rs_hubbard

    function get_helement_rs_hub_ex_mat(nI, nJ, ic, ex, tpar) result(hel)
        integer, intent(in) :: nI(nel), nJ(nel)
        integer, intent(in) :: ic, ex(2,2)
        logical, intent(in) :: tpar
        HElement_t(dp) :: hel 

        if (ic == 0) then 
            ! diagonal matrix element -> sum over doubly occupied orbitals! 
            hel = get_diag_helemen_rs_hub(nI) 

        else if (ic == 1) then 
            ! one-body operator: 
            ! here we need to make the distinction, if we are doing a 
            ! transcorrelated hamiltonian or not 
            hel = get_offdiag_helement_rs_hub(nI, ex(:,1), tpar)

        else 
            ! zero matrix element! 
            hel = h_cast(0.0_dp)

        end if

    end function get_helement_rs_hub_ex_mat

    function get_helement_rs_hub_general(nI, nJ, ic_ret) result(hel)
        integer, intent(in) :: nI(nel), nJ(nel) 
        integer, intent(inout), optional :: ic_ret 
        HElement_t(dp) :: hel
        
        integer :: ic, ex(2,2)
        logical :: tpar
        integer(n_int) :: ilutI(0:NIfTot), ilutJ(0:NIfTot)

        if (present(ic_ret)) then 
            if (ic_ret == 0) then 
                hel = get_diag_helemen_rs_hub(nI)

            else if (ic_ret == 1) then 
                ex(1,1) = 1
                call GetExcitation(nI, nJ, nel, ex, tpar)
                hel = get_offdiag_helement_rs_hub(nI, ex(:,1), tpar) 

            else if (ic_ret == -1) then 
                ! this indicates that ic_ret wants to get returned instead of 
                ! beeing calculated 
                ! its the same as if no ic_ret is present 
                call EncodeBitDet(nI, ilutI)
                call EncodeBitDet(nJ, ilutJ)

                ic_ret = FindBitExcitLevel(ilutI, ilutJ)

                if (ic_ret == 0) then 
                    hel = get_diag_helemen_rs_hub(nI)
                else if (ic_ret == 1) then 
                    ex(1,1) = 1
                    call GetBitExcitation(ilutI, ilutJ, ex, tpar)

                    hel = get_offdiag_helement_rs_hub(nI, ex(:,1), tpar)

                else 
                    hel = h_cast(0.0_dp)

                end if
            else
                hel = h_cast(0.0_dp)
            end if
        else
            call EncodeBitDet(nI, ilutI)
            call EncodeBitDet(nJ, ilutJ)

            ic = FindBitExcitLevel(ilutI, ilutJ)

            if (ic == 0) then 
                hel = get_diag_helemen_rs_hub(nI)
            else if (ic == 1) then 
                ex(1,1) = 1
                call GetBitExcitation(ilutI, ilutJ, ex, tpar)
                hel = get_offdiag_helement_rs_hub(nI, ex(:,1), tpar) 

            else
                hel = h_cast(0.0_dp)
            end if
        end if
    
    end function get_helement_rs_hub_general

    ! also optimize the matrix element calculation
    function get_diag_helemen_rs_hub(nI) result(hel)
        use double_occ_mod, only: count_double_orbs
        integer, intent(in) :: nI(nel)
        HElement_t(dp) :: hel

        integer(n_int) :: ilut(0:NIfTot)
        ! the diagonal matrix element is essentialy just the number of 
        ! doubly occupied sites times U

        call EncodeBitDet(nI, ilut)

        hel = h_cast(uhub * count_double_orbs(ilut))

    end function get_diag_helemen_rs_hub

    function get_offdiag_helement_rs_hub(nI, ex, tpar) result(hel)
        integer, intent(in) :: nI(nel), ex(2)
        logical, intent(in) :: tpar
        HElement_t(dp) :: hel

        real(dp) :: ni_opp, nj_opp

        ! in case we need it, the off-diagonal, except parity is just 
        ! -t if the hop is possible
        hel = GetTMatEl(ex(1),ex(2))

        if (tpar) hel = -hel

        ! put the transcorrelated stuff here for now, alhtough it would be 
        ! better to do it as a procedure pointer.. 
        if (t_trans_corr) then 
            if (is_beta(ex(1))) then 
                if (any(get_alpha(ex(1)) == nI)) then 
                    nj_opp = 1.0_dp
                end if 
                if (any(get_alpha(ex(2)) == nI)) then 
                    ni_opp = 1.0_dp
                end if
            else 
                if (any(get_beta(ex(1)) == nI)) then 
                    nj_opp = 1.0_dp
                end if 
                if (any(get_beta(ex(2)) == nI)) then 
                    ni_opp = 1.0_dp 
                end if
            end if

            hel = hel * exp(trans_corr_param * (nj_opp - ni_opp))

        end if

    end function get_offdiag_helement_rs_hub

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
        
        
