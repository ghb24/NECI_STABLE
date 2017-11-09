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
                          t_open_bc_y, t_open_bc_z, G1, ecore, nel, nOccAlpha, nOccBeta, & 
                          t_trans_corr, trans_corr_param
    use lattice_mod, only: lattice, determine_optimal_time_step, lat, &
                    get_helement_lattice, get_helement_lattice_ex_mat, & 
                    get_helement_lattice_general
    use constants, only: dp, EPS, n_int, bits_n_int
    use procedure_pointers, only: get_umat_el, generate_excitation
    use OneEInts, only: tmat2d, GetTMatEl
    use fcimcdata, only: pSingles, pDoubles, tsearchtau, tsearchtauoption
    use CalcData, only: t_hist_tau_search, t_hist_tau_search_option, tau
    use umatcache, only: gtid
    use dsfmt_interface, only: genrand_real2_dsfmt
    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet
    use bit_rep_data, only: NIfTot
    use util_mod, only: binary_search_first_ge
!     use fcimc_helper, only: update_run_reference
!     use Determinants, only: write_det

    implicit none 

    real(dp) :: lat_tau_factor = 0.5_dp

    ! create a flag which indicate to start in a neel state 
    logical :: t_start_neel_state = .false. 

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
        use SystemData, only: tExch, thub, treal, tHPHF
        ! routine, which does all of the necessary initialisation
        character(*), parameter :: this_routine = "init_real_space_hubbard"
        ! especially do some stop_all here so no wrong input is used 
        real(dp) :: tau_opt
        integer :: neel_state_ni(nel)
        integer(n_int) :: ilut_neel(0:NIfTot)

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

        ! also use the new lattice matrix elements

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
        if (.not. tHPHF) then
            generate_excitation => gen_excit_rs_hubbard
        end if
        
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
        
        if (t_start_neel_state) then 
!             neel_state_ni = create_neel_state(ilut_neel)

            print *, "starting from the Neel state: " 
            if (nel > nbasis/2) then 
                call stop_all(this_routine, &
                    "more than half-filling! does neel state make sense?")
            end if
!             call write_det(6, neel_state_ni, .true.)
!             call changerefdet(neel_state_ni)
!             call update_run_reference(ilut_neel, 1)

        end if
        ! do not set that here, due to circular dependencies
!         max_death_cpt = 0.0_dp

    end subroutine init_real_space_hubbard

    subroutine init_get_helement_hubbard
        get_helement_lattice_ex_mat => get_helement_rs_hub_ex_mat
        get_helement_lattice_general => get_helement_lattice_general
    end subroutine init_get_helement_hubbard

    subroutine check_real_space_hubbard_input() 
        use SystemData, only: tCSF, tReltvy, tUEG, tUEG2, tHub, & 
                              tKPntSym, tLatticeGens, tUEGNewGenerator, &
                tGenHelWeighted, tGen_4ind_weighted, tGen_4ind_reverse, &
                tUEGNewGenerator, tGen_4ind_part_exact, tGen_4ind_lin_exact, &
                tGen_4ind_2, tGen_4ind_2_symmetric, tGen_4ind_unbound, tStoreSpinOrbs, &
                tReal
        use umatcache, only : tTransGTid
        use OneEInts, only: tcpmdsymtmat, tOneelecdiag

        character(*), parameter :: this_routine = "check_real_space_hubbard_input"
        ! do all the input checking here, so no wrong input is used!

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
        real(dp), allocatable :: cum_arr(:)
        real(dp) :: cum_sum, elem, r, p_elec, p_orb

        iunused = exflag; 

        ASSERT(associated(lat))

        ic = 1
        ! i only have single excitations in the hubbard model 
        ! the first plan is to choose an electron at random 
        elec = 1 + int(genrand_real2_dsfmt() * nel) 

        p_elec = 1.0_dp / real(nel, dp)
        ! and then from the neighbors of this electron we pick an empty 
        ! spinorbital randomly, since all have the same matrix element 
        src = nI(elec) 

        ! get the spatial index 
        id = gtid(src) 

        ! now get neighbors
!         n_orbs = lat%get_num_neighbors(id)
        neighbors = lat%get_spinorb_neighbors(src)

!         allocate(neighbors(n_orbs)) 
!         allocate(orbs(n_orbs))
!         neighbors = lat%get_neighbors(id) 

!         if (is_beta(src)) then 
!             neighbors = 2 * neighbors - 1
!         else 
!             neighbors = 2 * neighbors
!         end if


        if (t_trans_corr) then 
            call create_cum_list_rs_hubbard(ilutI, src, neighbors, cum_arr, cum_sum)

            if (cum_sum < EPS) then 
                nJ(1) = 0
                pgen = 0.0_dp 
                return
            end if
            r = genrand_real2_dsfmt() * cum_sum 
            ind = binary_search_first_ge(cum_arr, r)

            if (ind == 1) then 
                p_orb = cum_arr(1) / cum_sum
            else
                p_orb = (cum_arr(ind) - cum_arr(ind-1)) / cum_sum 
            end if

            orb = neighbors(ind)
        else

            call create_avail_neighbors_list(ilutI, neighbors, orbs, n_avail)

            if (n_avail == 0) then 
                nJ(1) = 0
                pgen = 0.0_dp 
                return
            end if

            ind = 1 + int(genrand_real2_dsfmt() * n_avail)
            p_orb = 1.0_dp / real(n_avail, dp) 

            orb = orbs(ind) 

        end if

        pgen = p_elec * p_orb

        call make_single(nI, nJ, elec, orb, ex, tParity) 

        ilutJ = ilutI 
        clr_orb(ilutJ, src)
        set_orb(ilutJ, orb)

    end subroutine gen_excit_rs_hubbard

    function calc_pgen_rs_hubbard(nI, ilutI, ex, ic) result(pgen) 
        ! i also need a pgen recalculator.. specifically for the HPHF 
        ! implementation and i need to take the transcorrelated keyword 
        ! into account here! 
        integer, intent(in) :: nI(nel), ex(2,2), ic 
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        real(dp) :: pgen 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "calc_pgen_rs_hubbard"
#endif
        integer :: src, tgt, n_orbs
        real(dp) :: p_elec, p_orb, cum_sum
        real(dp), allocatable :: cum_arr(:)
        integer, allocatable :: orbs(:)

        ! only single excitations in the real-space hubbard
        if (ic /= 1) then 
            pgen = 0.0_dp 
            return
        end if

        src = ex(1,1)
        tgt = ex(2,1)

        ! can i assert the same spin of the 2 involved orbitals? 
        ! just return 0 if both have different spin? 
        ASSERT(is_beta(src) .eqv. is_beta(tgt))
        ! and assert that we actually take a valid excitation:
        ASSERT(any(tgt == lat%get_spinorb_neighbors(src)))
        ASSERT(IsOcc(ilutI, src))
        ASSERT(IsNotOcc(ilutI, tgt))

        p_elec = 1.0_dp / real(nel, dp) 

        if (t_trans_corr) then 
            call create_cum_list_rs_hubbard(ilutI, src, lat%get_spinorb_neighbors(src), &
                cum_arr, cum_sum, tgt, p_orb) 
            if (cum_sum < EPS) then 
                pgen = 0.0_dp 
                return
            end if

        else 
            ! i should also write a routine which gives me the 
            ! neighboring orbitals and the number of possible hops 
            call create_avail_neighbors_list(ilutI, lat%get_spinorb_neighbors(src), & 
                orbs, n_orbs) 

                p_orb = 1.0_dp / real(n_orbs, dp) 

        end if

        pgen = p_elec * p_orb

    end function calc_pgen_rs_hubbard

    subroutine create_cum_list_rs_hubbard(ilutI, src, neighbors, cum_arr, cum_sum, &
            tgt, cpt)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: src, neighbors(:)
        real(dp), intent(out) :: cum_sum 
        real(dp), intent(out), allocatable :: cum_arr(:)
        integer, intent(in), optional :: tgt
        real(dp), intent(out), optional :: cpt
#ifdef __DEBUG
        character(*), parameter :: this_routine = "create_cum_list_rs_hubbard"
#endif

        real(dp) :: elem
        integer :: i 

        ASSERT(IsOcc(ilutI,src))

        allocate(cum_arr(size(neighbors)))
        cum_arr = 0.0_dp
        cum_sum = 0.0_dp
        if (present(tgt)) then 
            ASSERT(present(cpt))
            do i = 1, ubound(neighbors,1)
                ASSERT(is_beta(src) .eqv. is_beta(neighbors(i)))
                if (IsNotOcc(ilutI, neighbors(i))) then 
                    cum_sum = cum_sum + abs(trans_corr_fac(ilutI, src, neighbors(i)))
                end if
                if (neighbors(i) == tgt) then 
                    cpt = abs(trans_corr_fac(ilutI, src, neighbors(i)))
                end if
            end do
            if (cum_sum < EPS) then
                cpt = 0.0
            else 
                cpt = cpt / cum_sum
            end if
        else
            do i = 1, ubound(neighbors,1)
            ASSERT(is_beta(src) .eqv. is_beta(neighbors(i)))
            if (IsNotOcc(ilutI,neighbors(i))) then 
                cum_sum = cum_sum + abs(trans_corr_fac(ilutI, src, neighbors(i)))
            end if
            cum_arr(i) = cum_sum
            end do
        end if
        
    end subroutine create_cum_list_rs_hubbard

    subroutine create_avail_neighbors_list(ilutI, neighbors, orbs, n_orbs)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: neighbors(:) 
        integer, intent(out), allocatable :: orbs(:)
        integer, intent(out) :: n_orbs 

        integer :: i, temp_orbs(size(neighbors))

        n_orbs = 0
        temp_orbs = 0

        do i = 1, ubound(neighbors,1)
            if (IsNotOcc(ilutI,neighbors(i))) then 
                n_orbs = n_orbs + 1
                temp_orbs(n_orbs) = neighbors(i)
            end if
        end do

        allocate(orbs(n_orbs), source=temp_orbs(1:n_orbs))

    end subroutine create_avail_neighbors_list

    function trans_corr_fac(ilutI, src, tgt) result(weight)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: src, tgt
        real(dp) :: weight 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "trans_corr_fac"
#endif
        real(dp) :: ni_opp, nj_opp 

        ! if the spins are not the same, something went wrong..
        ASSERT(is_beta(src) .eqv. is_beta(tgt))

        ni_opp = 0.0_dp
        nj_opp = 0.0_dp

        if (is_beta(src)) then 
            ! check if alpha orbital (i) and (j) in ilutI is occupied
            if (IsOcc(ilutI,get_alpha(src))) then 
                nj_opp = 1.0_dp
            end if
            if (IsOcc(ilutI,get_alpha(tgt))) ni_opp = 1.0_dp
        else 
            if (IsOcc(ilutI,get_beta(src))) nj_opp = 1.0_dp
            if (IsOcc(ilutI,get_beta(tgt))) ni_opp = 1.0_dp
        end if

        weight = exp(trans_corr_param*(nj_opp - ni_opp))

    end function trans_corr_fac

    function get_helement_rs_hub_ex_mat(nI, ic, ex, tpar) result(hel)
        integer, intent(in) :: nI(nel)
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
    function create_neel_state(ilut_neel) result(neel_state)
        ! probably a good idea to have a routine which creates a neel state
        ! (or one of them if it is ambigous)
        integer(n_int), intent(out), optional :: ilut_neel(0:NIfTot)
        integer :: neel_state(nel)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "create_neel_state"
#endif
        integer :: i, j, k, l, spin, ind
        

        ! make this independent of the lattice mod, to call it as early 
        ! as possible in the initialization
!         ASSERT(associated(lat))

        ! also assert that we have at most half-filling.. otherwise it does 
        ! not make so much sense to do this 
!         ASSERT(nel <= nBasis/2)
!         ASSERT(nel <= lat%get_nsites())
        ! there is no neel state for all lattice, but atleast something 
        ! similar to it.. 
        ! atleast for the lattice where there is a neel state, i should 
        ! create it automaticall and for the other a similar one 
        ! the lattice should already be intialized 
        ! i was thinking of putting this functionality into the 
        ! lattice_mod, but i guess i should keep it seperate, since it 
        ! actually has to do more with the system or?
        neel_state = 0

        select case (lattice_type) 

        case ('chain')
            neel_state = create_neel_state_chain()

        case ('square','rectangle','triangle','triangular')
            ! check if length_x is mod 2 
            if (mod(length_x, 2) == 0) then 

                neel_state = create_neel_state_chain()

                ! and flip the spins in every other column 
                do i = 1, length_y, 2
                    ! loop over every second column 
                    if (i*length_x >= nel) exit
                    do j = 1, length_x
                        ! whats the total index? 
                        ind = i*length_x + j
                        if (ind > nel) exit

                        ! flip the spins.. this should be way easier.. 
                        if (is_beta(neel_state(ind))) then
                            neel_state(ind) = neel_state(ind) + 1
                        else
                            neel_state(ind) = neel_state(ind) - 1
                        end if
                    end do
                end do
            else 
                ! here it is easy it is just like the chain case 
                neel_state = create_neel_state_chain()

            end if

        case ('cube') 
            ! not yet implemented
            ASSERT(.false.)

        case ('tilted','tilted-square','square-tilted')
            ! do is similar to the actual construction of the lattice 
            ! but way nicer as this stuff below..
            k = 0
            l = 1 
            spin = 1
            do i = -length_x+1, 0
                do j = -k, k

                    neel_state(l) = spin 
                    l = l + 1
                    if (l > nel) exit 

                    if (is_beta(spin)) then 
                        spin = spin + 3 
                    else 
                        spin = spin + 1 
                    end if 
                end do
                k = k + 1 
                if (l > nel) exit 

                ! in the first half we need the same starting spin 
                spin = spin - 1

            end do
            k = k - 1

            ! finished already?
            if (l > nel) return

            spin = spin + 1

            do i = 1, length_x 
                do j = -k, k

                    neel_state(l) = spin 

                    l = l + 1 

                    if (l > nel) exit 

                    if (is_beta(spin)) then 
                        spin = spin + 3 
                    else 
                        spin = spin + 1
                    end if
                end do
                k = k - 1 
                if (l > nel) exit 

                spin = spin + 1

            end do

        end select 

        if (present(ilut_neel)) then 
            call EncodeBitDet(neel_state, ilut_neel)
        end if
    end function create_neel_state

    function create_neel_state_chain() result(neel_state)
        integer :: neel_state(nel)

        integer :: i

        neel_state = [(i, i = 1, 2*nel-1,2)]
        neel_state(2:nel:2) = neel_state(2:nel:2) + 1

    end function create_neel_state_chain

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
        
        
