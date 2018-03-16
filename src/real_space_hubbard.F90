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
                          t_trans_corr, trans_corr_param, t_trans_corr_2body, & 
                          trans_corr_param_2body, tHPHF, t_trans_corr_new, & 
                          t_trans_corr_hop, t_uniform_excits

    use lattice_mod, only: lattice, determine_optimal_time_step, lat, &
                    get_helement_lattice, get_helement_lattice_ex_mat, & 
                    get_helement_lattice_general, init_dispersion_rel_cache, &
                    epsilon_kvec, setup_lattice_symmetry

    use constants, only: dp, EPS, n_int, bits_n_int, pi

    use procedure_pointers, only: get_umat_el, generate_excitation

    use OneEInts, only: tmat2d, GetTMatEl

    use fcimcdata, only: pSingles, pDoubles, tsearchtau, tsearchtauoption, &
                        excit_gen_store_type

    use CalcData, only: t_hist_tau_search, t_hist_tau_search_option, tau, & 
                        t_fill_frequency_hists, p_singles_input

    use dsfmt_interface, only: genrand_real2_dsfmt

    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet, ilut_lt, ilut_gt

    use bit_rep_data, only: NIfTot, nifd

    use util_mod, only: binary_search_first_ge, choose, swap, get_free_unit, &
                        binary_search

    use bit_reps, only: decode_bit_det

    use sort_mod, only: sort

    use back_spawn, only: make_ilutJ, get_ispn
    
    use get_excit, only: make_double, make_single

    use double_occ_mod, only: count_double_orbs

    implicit none 

    real(dp) :: lat_tau_factor = 0.5_dp

    ! create a flag which indicate to start in a neel state 
    logical :: t_start_neel_state = .false. 

    real(dp), allocatable :: umat_rs_hub_trancorr_hop(:,:,:,:)

    complex(dp) :: imag_unit = complex(0.0_dp,1.0_dp)

    interface get_helement_rs_hub
        module procedure get_helement_rs_hub_ex_mat
        module procedure get_helement_rs_hub_general
    end interface get_helement_rs_hub

    interface swap_excitations 
        module procedure swap_excitations_higher
        module procedure swap_excitations_singles
    end interface swap_excitations

    interface sum_hop_transcorr_factor
        module procedure sum_hop_transcorr_factor_vec
        module procedure sum_hop_transcorr_factor_orb
    end interface sum_hop_transcorr_factor
            
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
        integer :: neel_state_ni(nel)
        integer(n_int) :: ilut_neel(0:NIfTot)

        print *, "using new real-space hubbard implementation: "


        ! i do not need exchange integrals in the real-space hubbard model
        if (.not. t_trans_corr_hop) then
            tExch = .false.
        end if
        ! after the whole setup i can set thub to false or? 
        thub = .false.
        ! and treal i can also set to false or? 
        treal = .false.
        ! first assert all the right input! 
        call check_real_space_hubbard_input() 

        ! which stuff do i need to initialize here? 
        if (t_trans_corr_hop) then 
            get_umat_el => get_umat_rs_hub_trans
        else
            get_umat_el => get_umat_el_hub
        end if

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

        if (t_trans_corr_hop) then 
            ! we have double excitations with the hopping correlation! 
            pDoubles = p_singles_input 
            pSingles = 1.0_dp - pDoubles
        else
            ! and i have to point to the new hubbard excitation generator
            pSingles = 1.0_dp 
            pDoubles = 0.0_dp
        end if

        ! and i have to calculate the optimal time-step for the hubbard models. 
        ! where i need the connectivity of the lattice i guess? 
        if (t_trans_corr_hop) then 
            ASSERT(.not. tHPHF)
            if (t_uniform_excits) then 
                generate_excitation => gen_excit_rs_hubbard_transcorr_uniform
            else
                generate_excitation => gen_excit_rs_hubbard_transcorr
            end if
        else
            if (.not. tHPHF) then
                generate_excitation => gen_excit_rs_hubbard
            end if
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

        ! re-enable tau-search if we have transcorrelation 
        if (.not. (t_trans_corr_2body .or. t_trans_corr .or. t_trans_corr_hop)) then 
            ! and i have to turn off the time-step search for the hubbard 
            tsearchtau = .false.
            ! set tsearchtauoption to true to use the death-tau search option
            tsearchtauoption = .true.

            t_hist_tau_search = .false. 
            t_hist_tau_search_option = .false. 

            t_fill_frequency_hists = .false.
        end if
        
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

        ! i need to setup the necessary stuff for the new hopping 
        ! transcorrelated real-space hubbard! 
!         call init_hopping_transcorr() 

        call init_get_helement_hubbard()

    end subroutine init_real_space_hubbard

    subroutine init_hopping_transcorr()

        ! we also need the dispersion relation from the k-space hubbard! 
        ! should i move this to the lattice class?
        ! this also means i have to additionally initialize parts of the 
        ! k-space lattice for this transcorrelation! 
        ! and i finally have to get the real-space and k-space vector 
        ! relations fully correct! 
        integer :: i 

        call setup_lattice_symmetry()

        call init_dispersion_rel_cache() 

!         print *, "e(k)"
!         do i = 1, lat%get_nsites()
!             print *, i, "|", lat%get_k_vec(i), "|", epsilon_kvec(i)
!         end do

        ! i also need a umat array now! 
        call init_umat_rs_hub_transcorr()

    end subroutine init_hopping_transcorr

    real(dp) function hop_transcorr_factor(J, r_vec)
        ! compute \sum_p exp(i*p*r) exp(J,r) for the 2-body term in the 
        ! hopping transcorrelated real-space hubbard 
        real(dp), intent(in) :: J 
        integer, intent(in) :: r_vec(3) 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "hop_transcorr_factor"
#endif 
        complex(dp) :: temp
        integer :: i, n_sites, k_vec(3)

        ASSERT(associated(lat))
    
        n_sites = lat%get_nsites()

        temp = 0.0_dp 

        ! i have to have the correct dot-product of the real- and k-space 
        ! vectors... 
        do i = 1, n_sites
            k_vec = lat%get_k_vec(i)
!             temp = temp + exp(complex(0.0,1.0) * 2.0*pi/real(n_sites,dp) * &
!                 dot_product(k_vec, r_vec)) * exp(-J * epsilon_kvec(i))
            temp = temp + exp(imag_unit * lat%dot_prod(k_vec, r_vec)) * & 
                exp(-J * epsilon_kvec(i))

!             print *, "k: ", k_vec, "|", exp(complex(0.0,1.0) * 2.0*pi/real(n_sites,dp)* &
!                 dot_product(k_vec, r_vec))

        end do

        hop_transcorr_factor = real(temp) / real(n_sites,dp)

!         print *, "hop_transcorr_factor(r): ",r_vec, temp/real(n_sites,dp)
        ! will this be ever comlex? check the size 
        ASSERT(abs(aimag(temp)) < EPS) 


    end function hop_transcorr_factor

    real(dp) function sum_hop_transcorr_factor_vec(r1,r2,r3,r4)
        ! function to perform the summation over the four hopping 
        ! transcorrelation factors in the 2-body term of the hopping 
        ! transcorrelated real-space hubbard model 
        integer, intent(in) :: r1(3), r2(3), r3(3), r4(3)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "sum_hop_transcorr_factor_vec"
#endif
        integer i, m_vec(3)

        ASSERT(associated(lat))

        sum_hop_transcorr_factor_vec = 0.0_dp

        do i = 1, lat%get_nsites() 

            ! i need a routine to give me the real-space coordinates/vector
            ! and i also need to store that now! 
            m_vec = lat%get_r_vec(i)
            sum_hop_transcorr_factor_vec = sum_hop_transcorr_factor_vec + & 
                hop_transcorr_factor(trans_corr_param, r1 - m_vec) * &
                hop_transcorr_factor(trans_corr_param, r2 - m_vec) * &
                hop_transcorr_factor(-trans_corr_param, m_vec - r3) * &
                hop_transcorr_factor(-trans_corr_param, m_vec - r4)

        end do

    end function sum_hop_transcorr_factor_vec

    real(dp) function sum_hop_transcorr_factor_orb(i,j,k,l)
        ! function to perform the summation over the four hopping 
        ! transcorrelation factors in the 2-body term of the hopping 
        ! transcorrelated real-space hubbard model 
        integer, intent(in) :: i, j, k, l
#ifdef __DEBUG
        character(*), parameter :: this_routine = "sum_hop_transcorr_factor_orb"
#endif
        integer :: r1(3), r2(3), r3(3), r4(3)
        integer m, m_vec(3)

        ASSERT(associated(lat))

        sum_hop_transcorr_factor_orb = 0.0_dp

        r1 = lat%get_r_vec(i)
        r2 = lat%get_r_vec(j)
        r3 = lat%get_r_vec(k)
        r4 = lat%get_r_vec(l)

        do m = 1, lat%get_nsites() 

            ! i need a routine to give me the real-space coordinates/vector
            ! and i also need to store that now! 
            m_vec = lat%get_r_vec(m)
            sum_hop_transcorr_factor_orb = sum_hop_transcorr_factor_orb + & 
                hop_transcorr_factor(trans_corr_param, r1 - m_vec) * &
                hop_transcorr_factor(trans_corr_param, r2 - m_vec) * &
                hop_transcorr_factor(-trans_corr_param, m_vec - r3) * &
                hop_transcorr_factor(-trans_corr_param, m_vec - r4)

        end do

    end function sum_hop_transcorr_factor_orb

    subroutine init_umat_rs_hub_transcorr() 
        ! for now do it really stupidly without any concern for the symmetry 
        ! of the integrals 
        integer :: n_sites, i, j, k, l
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "init_umat_rs_hub_transcorr"
#endif
        real(dp) :: elem 
        integer :: iunit

        ASSERT(associated(lat))

        if (allocated(umat_rs_hub_trancorr_hop)) then 
            ! already initialized
            return
        end if

        n_sites = lat%get_nsites()

        ! create an fcidump file: 
        iunit = get_free_unit()
        open(iunit, file = 'UMAT')

        ! with the correct header

        ! and also try to allocate a umat_cache 
        allocate(umat_rs_hub_trancorr_hop(n_sites,n_sites,n_sites,n_sites))
        umat_rs_hub_trancorr_hop = 0.0_dp

        do i = 1, n_sites
            do j = 1, n_sites
                do k = 1, n_sites
                    do l = 1, n_sites
                        elem = uhub * sum_hop_transcorr_factor(i,j,k,l)

                        ! write to the dumpfile
                        write(iunit,*)  i,j,k,l, elem

                        ! and also store in the umat 
                        umat_rs_hub_trancorr_hop(i,j,k,l) = elem 
                    end do
                end do
            end do
        end do
        close(iunit)

    end subroutine init_umat_rs_hub_transcorr

    subroutine init_get_helement_hubbard

        get_helement_lattice_ex_mat => get_helement_rs_hub_ex_mat
        get_helement_lattice_general => get_helement_rs_hub_general
        if (t_trans_corr_hop) then 
            call init_hopping_transcorr()
        end if

        call init_tmat(lat)

    end subroutine init_get_helement_hubbard

    subroutine check_real_space_hubbard_input() 
        use SystemData, only: tCSF, tReltvy, tUEG, tUEG2, tHub, & 
                              tKPntSym, tLatticeGens, tUEGNewGenerator, &
                tGenHelWeighted, tGen_4ind_weighted, tGen_4ind_reverse, &
                tUEGNewGenerator, tGen_4ind_part_exact, tGen_4ind_lin_exact, &
                tGen_4ind_2, tGen_4ind_2_symmetric, tGen_4ind_unbound, tStoreSpinOrbs, &
                tReal
!         use umatcache, only : tTransGTid
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
!         if (tTransGTid)         call stop_all(this_routine, "tTransGTid")
        if (tcpmdsymtmat)        call stop_all(this_routine, "tcpmdsymmat")
        if (tOneelecdiag)       call stop_all(this_routine, "tOneelecdiag")
            
    end subroutine check_real_space_hubbard_input

    ! i could also create all determinants, not only the open-shells.. or? 
    function create_all_dets(n_orbs, n_alpha, n_beta) result(all_dets)
        integer, intent(in) :: n_orbs, n_alpha, n_beta
        integer(n_int), allocatable :: all_dets(:) 
        character(*), parameter :: this_routine = "create_all_dets" 

        integer, allocatable :: n_dets(:)
        integer(n_int), allocatable :: alpha_basis(:), beta_basis(:)
        integer :: i, n_states, n_total, j

        n_dets = calc_n_double(n_orbs, n_alpha, n_beta) 

        ! cannot deal with more than 32 orbitals yet.. since i am lazy! 
        ASSERT(n_orbs <= 32) 

        ! do it in a really simple way for now.. 
        ! just run over all the alpha and beta spin-bases and combine them 
        ! to the full-basis 

        alpha_basis = create_one_spin_basis(n_orbs, n_alpha) 
        beta_basis = create_one_spin_basis(n_orbs, n_beta) 

        ! keep an count on all the states
        n_states = 0

        ! and allocate the array, which keeps all the states 
        n_total = sum(n_dets) 

        allocate(all_dets(n_total))
        all_dets = 0_n_int 

        ! check if everything went correctly: 
        ASSERT(size(alpha_basis) * size(beta_basis) == n_total)

        do i = 1, size(alpha_basis) 
            do j = 1, size(beta_basis)
                n_states = n_states + 1
                all_dets(n_states) = general_product(alpha_basis(i),beta_basis(j),n_orbs)
            end do
        end do
        
        call sort(all_dets)

        ASSERT(n_states == n_total) 

    end function create_all_dets 

    subroutine create_hilbert_space_realspace(n_orbs, n_alpha, n_beta, & 
            n_states, state_list_ni, state_list_ilut) 
        integer, intent(in) :: n_orbs, n_alpha, n_beta
        integer, intent(out) :: n_states 
        integer, intent(out), allocatable :: state_list_ni(:,:)
        integer(n_int), intent(out), allocatable :: state_list_ilut(:,:) 

        integer(n_int), allocatable :: all_dets(:) 
        integer(n_int) :: temp_ilut(0:niftot)
        integer :: nJ(nel), i

        all_dets = create_all_dets(n_orbs, n_alpha, n_beta) 

        n_states = size(all_dets)

        allocate(state_list_ni(nel,n_states))
        allocate(state_list_ilut(0:niftot,n_states))

        do i = 1, size(all_dets) 
            temp_ilut = all_dets(i) 
            call decode_bit_det(nJ, temp_ilut)

            state_list_ni(:,i) = nJ
            state_list_ilut(:,i) = temp_ilut

        end do

    end subroutine create_hilbert_space_realspace

    function create_all_open_shell_dets(n_orbs, n_alpha, n_beta) result(open_shells)
        integer, intent(in) :: n_orbs, n_alpha, n_beta
        ! Alis idea for the use of the ADI option for the real-space hubbard 
        ! is to let all fully open shell dets, and excitations thereof 
        ! (dets with exactly one double occupancy in the case of half-filling) 
        ! be initiators 
        character(*), parameter :: this_routine = "create_all_open_shell_dets" 

        integer, allocatable :: n_dets(:)
        integer(n_int), allocatable :: open_shells(:), max_basis(:)
        integer :: n_max, n_min

        n_dets = calc_n_double(n_orbs, n_alpha, n_beta)

        ! only implement that for less than 32 orbitals, since i do not want to 
        ! deal with multiple integers to store the basis and also the size of 
        ! the hilbert space would be too big anyway.. 
        ASSERT(n_orbs <= 32)

        ! how should i create all of those purely open-shell dets? 
        ! in the case of half-filling it is easiest.. 
        
        ! i could first distribute the alpha spins and then the beta spins 
        ! afterwards.. 
        ! i should work with the bit-representation, then i get to know the 
        ! fortran intrinsic bit operations again.. 

        n_max = max(n_alpha, n_beta)
        n_min = min(n_alpha, n_beta)
        
        ! first create the spin-basis with more electrons: 
        max_basis = create_one_spin_basis(n_orbs, n_max) 

        open_shells = combine_spin_basis(n_orbs, n_max, n_min, n_dets(1),  max_basis, .true.)

    end function create_all_open_shell_dets

    function combine_spin_basis(n_orbs, n_first, n_second, n_total, first_basis, & 
            t_sort, n_doubles) result(spin_basis) 
        ! function to combine a already given spin-basis with the second one
        ! with the additional option to specify the number of doubly occupied 
        ! determinants.. 
        ! n_orbs   ...  is the number of spatial orbitals
        ! n_first  ...  number of spins, with which first_basis was created 
        ! n_second ...  number of spins with which basis is combined now 
        ! n_total  ...  total size of basis for given number of doubly occ. sites!
        ! first_basis . already contructed one-spin basis
        ! t_sort   ...  should the basis be sorted. (after combination i am not  sure it will be sorted..
        ! n_doubles ..  (optional:) number of doubly occupied sites(has to fit with n_total!)
        integer, intent(in) :: n_orbs, n_first, n_second, n_total
        integer(n_int), intent(in) :: first_basis(:) 
        logical, intent(in) :: t_sort
        integer, intent(in), optional :: n_doubles 
        integer(n_int) :: spin_basis(n_total) 
        character(*), parameter :: this_routine = "combine_spin_basis" 
#ifdef __DEBUG
        integer, allocatable :: n_dets(:)
#endif 
        integer :: n_doub, i, j, n_remain, n
        integer(n_int), allocatable :: second_basis(:)
    
        ! the half-filled case is the most easy one.. should we treat it 
        ! special? maybe.. 

        ASSERT(n_first >= n_second) 
#ifdef __DEBUG 
        ! be sure that the provided n_total fits with the n_doubles if 
        ! provided 
        n_dets = calc_n_double(n_orbs, n_first, n_second)
        if (present(n_doubles)) then 
            ASSERT(n_dets(n_doubles+1) == n_total)
        else 
            ! n_doubles = 0 is the default
            ASSERT(n_dets(1) == n_total)
        end if
#endif

        ! default value of n_doubles = 0 
        if (present(n_doubles)) then 
            n_doub = n_doubles
        else 
            n_doub = 0 
        end if

        ! do the n_double = 0 and half-filled case first, thats the 
        ! easiest and only necessary for now.. 
        select case (n_doub) 
        case (0) 

            if (n_first + n_second == n_orbs) then 
                ! half-filled case: easiest 
                ! since ms is not important, just convert the already 
                ! calculated spin-basis to a spin-orbital basis: 
                ! 0011 -> 00 00 01 01 eg. 
                do i = 1, n_total 
                    ! do have it ordered -> first set the beta spins on the left
                    spin_basis(i) = set_alpha_beta_spins(first_basis(i), n_orbs,.false.) 
                    ! and combine it with the alpha spins..
                    spin_basis(i) = ieor(spin_basis(i),set_alpha_beta_spins(not(first_basis(i)), n_orbs,.true.))
                end do
            else 
                ! we have to distribute the n_second remaining spins 
                ! across the n_orbs - n_first orbitals, so there are: 
                n_remain = int(choose(n_orbs - n_first, n_second))
                ! states per first_basis states 
                ASSERT(n_remain * size(first_basis) == n_total) 

                second_basis = create_one_spin_basis(n_orbs - n_first, n_second)

                ! i want to do some sort of tensor product here:
                ! |0011> x |01> = |00 10 01 01>
                ! |0011> x |10> = |10 00 01 01>  
                n = 1 
                do i = 1, size(first_basis) 
                    do j = 1, size(second_basis)
                        spin_basis(n) = open_shell_product(first_basis(i),second_basis(j),n_orbs)
                        n = n + 1
                    end do
                end do

            end if

        case default 

            ! i should not need a select case here.. this gets to 
            ! complicated.. 
            ! it is similar to non-half filled actually.. 
            if (n_first + n_second == n_orbs) then 
                ! then we have the half-filled case.. 
            else

            end if

            ! not yet implemented 
            call stop_all(this_routine, "not yet implemented!")

        end select

        ! do i want to sort it?? 
        if (t_sort) then 
            call sort(spin_basis) 
        end if


    end function combine_spin_basis

    function general_product(alpha, beta, n_orbs) result(basis) 
        ! this is a "general" tensor product of two spin basisfunctions.
        ! the ones in the alpha integer must be set at position 2*i 
        ! and the beta at 2*j-1 
        integer(n_int), intent(in) :: alpha, beta 
        integer, intent(in) :: n_orbs 
        integer(n_int) :: basis 

        ! it should be as simple: 
        basis = set_alpha_beta_spins(alpha, n_orbs, .false.) 
        basis = xor(basis, set_alpha_beta_spins(beta, n_orbs, .true.))

    end function general_product

    function open_shell_product(alpha, beta, n_orbs) result(basis)
        ! compute the tensor product of alpha and beta spin-bases to give 
        ! fully open-shell determinants, the first input is to be defined as 
        ! the alpha spin part
        ! and the beta input then is only as big as n_orbs - n_alpha and only 
        ! tells us in which empty orbitals of the alpha orbitals beta spins 
        ! should be set 
        ! sorting is not so easily ensured here anymore.. 
        integer(n_int), intent(in) :: alpha, beta
        integer, intent(in) :: n_orbs
        integer(n_int) :: basis

        integer, allocatable :: nZeros(:), nOnes(:)
        integer(n_int) :: mask_zeros
        integer :: i
        ! the setting of the alpha spins is the same or? 
        basis = set_alpha_beta_spins(alpha, n_orbs, .false.)

        ! i need the zeros, but just in the n_orbs range
        mask_zeros = iand(not(alpha), maskr(n_orbs))

        allocate(nZeros(popcnt(mask_zeros)))

        call decode_bit_det(nZeros, [mask_zeros])

        ! now i need the ones from beta
        allocate(nOnes(popcnt(beta)))
        call decode_bit_det(nOnes, [beta])

        ! and i have to set all beta-spins indicated by nZeros(nOnes) 
        ! so we need beta spin! 
        nZeros = 2*nZeros - 1
        do i = 1, size(nOnes) 
            basis = ibset(basis, nZeros(nOnes(i))-1)
        end do

    end function open_shell_product

    function set_alpha_beta_spins(beta_mask, n_orbs, t_beta) result(beta_spins) 
        ! a function which converts a flag of spatial beta-spins to a 
        ! spin orbita basis, eg.: 0011 -> 00 00 01 01
        integer(n_int), intent(in) :: beta_mask
        integer, intent(in) :: n_orbs 
        logical, intent(in) :: t_beta
        integer(n_int) :: beta_spins

        integer, allocatable :: nOnes(:) 
        integer :: i 

        allocate(nOnes(popcnt(iand(beta_mask, maskr(n_orbs))))) 

        call decode_bit_det(nOnes, [beta_mask])

        if (t_beta) then 
            ! then we want to set beta spins: 
            nOnes = 2*nOnes - 1
        else 
            ! otherwise we want to set alpha spins
            nOnes = 2*nOnes
        end if
                
        beta_spins = 0_n_int
        do i = 1, size(nOnes) 
            beta_spins = ibset(beta_spins, nOnes(i)-1)
        end do

    end function set_alpha_beta_spins

    function create_one_spin_basis(n_orbs, n_spins) result(one_spin_basis) 
        integer, intent(in) :: n_orbs, n_spins 
        integer(n_int), allocatable :: one_spin_basis(:) 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "create_one_spin_basis" 
#endif

        integer :: n_max_states, i, right_zero
        integer(n_int) :: count_mask, n_set_zero

        n_max_states = int(choose(n_orbs, n_spins))

        ! 
        allocate(one_spin_basis(n_max_states)) 

        ! create the first basis state: 
        one_spin_basis(1) = maskr(n_spins) 

        ! implement my matlab routine to create the basis states: 
        ! but now i have to do it for spin-orbital encoding!
        ! but do this afterwards so the mapping is easier for off-half-filling!
        do i = 2, n_max_states

            ! copy last state: 
            one_spin_basis(i) = one_spin_basis(i-1)

            ! find the right-most zero with atleast one 1 right of it
            right_zero = right_most_zero(one_spin_basis(i), n_orbs)

            ! if the right-most zero is bigger than n_orbs already, we should 
            ! have been finished already.. 
            ASSERT(right_zero <= n_orbs)

            ! i need to count the number of 1 right of this zero 
            ! so i want a mask where every bit right of right_zero is set 
            ! to one 
            count_mask = maskr(right_zero-1) 
            n_set_zero = popcnt(iand(one_spin_basis(i),count_mask)) 

            ! now i want to set the right_most zero to one  
            one_spin_basis(i) = ibset(one_spin_basis(i), right_zero-1)

            ! and everything to 0 right of it for now 
            one_spin_basis(i) = iand(one_spin_basis(i), not(count_mask)) 

            ! and then we want to set n_set_zero bits to the very right of 
            ! our bit-string 
            one_spin_basis(i) = merge_bits(one_spin_basis(i), maskr(n_set_zero-1,n_int), not(maskr(n_set_zero-1,n_int)))

        end do

    end function create_one_spin_basis

    integer function right_most_zero(i, n_orbs) 
        ! gives the position of the right most zero with atleast one 1 
        ! right of it in an integer bit representation, up to a input 
        ! dependent orbital n_orbs 
        integer(n_int), intent(in) :: i
        integer, intent(in) :: n_orbs
#ifdef __DEBUG
        character(*), parameter :: this_routine = "right_most_zero" 
#endif
        integer, allocatable :: nZeros(:), nOnes(:) 
        integer(n_int) :: j 

        ! first truncate the integer to be save.. 
        ! set all the ones left of n_orbs to 0 
        j = iand(i, maskr(n_orbs))

        ! be sure to avoid the edge case: 
        if (j == 0_n_int) then
            right_most_zero = -1 
            return
        end if

        ! i could misuse the decode_bit_det functionality 
        allocate(nZeros(popcnt(not(j))))
        allocate(nOnes(popcnt(j)))

        ! but exclude the unnecessary 0 left of n_orbs..
        call decode_bit_det(nZeros, [not(j)])
        call decode_bit_det(nOnes, [j])

        ! in this way to encode it.. the right most 0 is then the minumum 
        ! in nZeros which is bigger then the minumum in nOnes
        right_most_zero = minval(nZeros, nZeros > minval(nOnes))

    end function right_most_zero

    function calc_n_double(n_orbs, n_alpha, n_beta) result(n_double)
        ! routine to calculate the number of determinants with different 
        ! number of doubly occupied sites..
        integer, intent(in) :: n_orbs, n_alpha, n_beta
        integer, allocatable :: n_double(:)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "calc_n_double"
#endif

        real(dp) :: n_first
        integer :: n_max, n_min, i
        ! the number of possible open shell determinants is:
        ! for the half-filled ms=0 case, it is how often we can distribute 
        ! n/2 electrons in n orbitals (since the rest then gets filled up 
        ! by the other spin electrons -> nchoosek(n, n/2) 
        ! todo: merge the cc_amplitudes branch or atleast the 
        ! functionality in there! 
        ! for ms/=0 but still at half-filling it i still the binomial 
        ! coefficient of one spin-type 

        ! off-half filling, we get more then one possible distribution of 
        ! the electrons of the other spin for each distribution of the 
        ! first spin type.. 
        ! so we get some product of binomial of the remaining possible 
        ! distributions! 

        ! todo: write those formulas up! 

        ! todo: a way to create it would be like in my matlab code to 
        ! create the heisenberg basis states.. there i even did it in 
        ! a ordered fashion.

        ! the number of ways to distribute alpha spins on the lattice 
!         n_alpha = binomial(nbasis/2, nOccAlpha)

!         num_open_shell_dets = n_alpha * binomial(nbasis/2 - nOccAlpha, nOccBeta) 

        ! maybe the number of determinants with one double occupancy is 
        ! also of interest.. (those are the single excitations from the 
        ! fully open shell dets(in the case off-half-filling it could a
        ! single excitation can also lead to another fully-open-shell det, but
        ! these are already covered.
        ! one beta electron MUST reside in a orbital with a alpha spin (thats 
        ! just nOccAlpha) and the rest (nOccBeta-1) is on an empty orbital
!         num_one_double_occ = n_alpha * nOccAlpha * & 
!             binomial(nbasis/2 - nOccAlpha, nOccBeta - 1)

        ! just for the fun of it: i could calculate the size of each of 
        ! those space(1 double occ, 2 double occ.. etc. 

        ! only do that for less equal half-filling! 
!         ASSERT(n_alpha + n_beta <= n_orbs)

        n_max = max(n_alpha, n_beta) 
        n_min = min(n_alpha, n_beta)

        n_first = choose(n_orbs, n_max)

        allocate(n_double(n_min+1)) 

        do i = 0, n_min 
            n_double(i+1) = int(n_first * choose(n_max, i) * & 
                choose(n_orbs - n_max, n_min - i))
        end do

        ! make sure that the sum of basis states is the whole hilber space
        ASSERT(sum(n_double) == choose(n_orbs, n_alpha)*choose(n_orbs,n_beta))

    end function calc_n_double

    function get_optimal_correlation_factor() result(corr_factor) 
        ! Hongjuns derivation was for the k-space hubbard and in the low 
        ! density and U limit though.. 
        real(dp) :: corr_factor
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_optimal_correlation_factor" 
#endif

        ASSERT(associated(lat)) 

        ! the sign is not quite sure here.. which i need to take to 
        ! calculate the hermitian matrix elements..
        corr_factor = -log(abs(real(uhub,dp)/real(4 * lat%get_ndim() * bhub, dp)) + 1.0_dp)


    end function get_optimal_correlation_factor

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

!             if (t_trans_corr_hop) then 
                ! for the new transcorrelation term, we have modified 
                ! one-body matrix elements.. but these are dependent on the 
                ! determinant!! so we cannot setup the matrix element 
                ! beforehand

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

    subroutine gen_excit_rs_hubbard_transcorr_uniform (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)
        ! also create an uniform excitation generator for the hopping 
        ! transcorrelated hubbard. mainly to test where the instabilities 
        ! in the weighted come from 


        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: run
        character(*), parameter :: this_routine = "gen_excit_rs_hubbard_transcorr_uniform"

        integer :: iunused, elecs(2), orbs(2), src(2), spin
        real(dp) :: p_elec, p_orb

        iunused = exflag; 
        ilutJ = 0_n_int
        ic = 0
        nJ(1) = 0

        ASSERT(associated(lat))

        if (genrand_real2_dsfmt() < pDoubles) then 

            ic = 2 

            call pick_spin_opp_elecs(nI, elecs, p_elec)

            src = nI(elecs) 

            ASSERT(.not. same_spin(src(1),src(2)))

            call pick_spin_opp_holes(ilutI, orbs, p_orb)

            ASSERT(.not. same_spin(orbs(1),orbs(2)))

            if (any(orbs == 0)) then 
                nJ(1) = 0
                pgen = 0.0_dp 
                return 
            end if

            ! do i need a factor of 2? since orbitals could be switched 
            ! the other way around?
            pgen = p_elec * p_orb * pDoubles

            call make_double(nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), ex, tParity)
            ilutJ = make_ilutJ(ilutI, ex, 2) 

        else 

            ic = 1

            elecs(1) = 1 + int(genrand_real2_dsfmt() * nel) 

            src(1) = nI(elecs(1))

            p_elec = 1.0_dp / real(nel,dp)

            spin = get_spin(src(1)) - 1
            ! and now pick a spin-parallel hole! 
            call pick_random_hole(ilutI, orbs(1), p_orb, spin)

            if (orbs(1) == 0) then 
                nJ(1) = 0
                pgen = 0.0_dp 
                return 
            end if

            pgen = p_elec * p_orb * (1.0_dp - pDoubles)

            call make_single(nI, nJ, elecs(1), orbs(1), ex, tParity)
            ilutJ = make_ilutJ(ilutI, ex, 1)

        end if

    end subroutine gen_excit_rs_hubbard_transcorr_uniform

    function calc_pgen_rs_hubbard_transcorr_uniform(nI, ilutI, ex, ic) result(pgen)
        integer, intent(in) :: nI(nel), ex(2,2), ic
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        real(dp) :: pgen
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "calc_pgen_rs_hubbard_transcorr_uniform"
#endif 

        if (ic == 1) then 

            ASSERT(same_spin(ex(1,1),ex(2,1)))

            if (is_beta(ex(1,1))) then 
                pgen = 1.0_dp / real(nel * (nBasis/2 - nOccBeta), dp)
            else 
                pgen = 1.0_dp / real(nel * (nBasis/2 - nOccAlpha), dp)
            end if

        else if (ic == 2) then 

            ASSERT(.not. same_spin(ex(1,1), ex(1,2)))
            ASSERT(.not. same_spin(ex(2,1), ex(2,2)))

            pgen = 1.0_dp / real(nOccAlpha * nOccBeta * &
                (nBasis/2 - nOccAlpha) * (nBasis/2 - nOccBeta), dp)

        else 

            pgen = 0.0_dp 

        end if

    end function calc_pgen_rs_hubbard_transcorr_uniform

    subroutine pick_spin_opp_holes(ilutI, orbs, p_orb) 
        ! routine to pick two spin-opposite unoccupied orbitals 
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: orbs(2) 
        real(dp), intent(out) :: p_orb 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "pick_spin_opp_holes"
#endif

        integer, parameter :: max_trials = 100
        integer :: cnt 

        cnt = 0
        orbs = 0
        p_orb = 0.0_dp 

        do while (cnt < max_trials)
            orbs(1) = 1 + int(genrand_real2_dsfmt() * nBasis) 

            if (IsOcc(ilutI,orbs(1))) cycle 

            do 
                orbs(2) = 1 + int(genrand_real2_dsfmt() * nBasis) 

                if (IsOcc(ilutI,orbs(2))) cycle 

                if (orbs(1) /= orbs(2)) exit 

            end do

            if (.not. same_spin(orbs(1),orbs(2))) exit 

        end do

        orbs = [minval(orbs),maxval(orbs)]

        p_orb = 1.0_dp / real((nBasis/2 - nOccAlpha)*(nBasis/2 - nOccBeta),dp)

    end subroutine pick_spin_opp_holes

    subroutine gen_excit_rs_hubbard_transcorr (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)
        ! new excitation generator for the real-space hubbard model with 
        ! the hopping transcorrelation, which leads to double excitations 
        ! and long-range single excitations in the real-space hubbard.. 
        ! this complicates things alot! 

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: run

        character(*), parameter :: this_routine = "gen_excit_rs_hubbard_transcorr"

        integer :: iunused, ind , elec, src, orb
        real(dp) :: cum_arr(nBasis/2)
        real(dp) :: cum_sum, p_elec, p_orb

        iunused = exflag; 
        ilutJ = 0_n_int
        ic = 0
        nJ(1) = 0

        ASSERT(associated(lat))

        if (genrand_real2_dsfmt() < pDoubles) then 
            ic = 2 

            call gen_double_excit_rs_hub_transcorr(nI, ilutI, nJ, ilutJ, ex, tParity, pgen)

            pgen = pDoubles * pgen 

            if (nJ(1) == 0) then 
                pgen = 0.0_dp 
                return 
            end if

        else 
            ic = 1
            
            ! still choose an electron at random
            elec = 1 + int(genrand_real2_dsfmt() * nel) 

            p_elec = 1.0_dp / real(nel, dp)
            ! and then from the neighbors of this electron we pick an empty 
            ! spinorbital randomly, since all have the same matrix element 
            src = nI(elec) 

            ! now we can have more than only nearest neighbor hopping! 
            ! so implement a new cum-list creator
            call create_cum_list_rs_hubbard_transcorr_single(nI, ilutI, src, cum_arr, cum_sum)

            ! the rest stays the same i guess..
            if (cum_sum < EPS) then 
                nJ(1) = 0
                pgen = 0.0_dp 
                return
            end if

            call pick_from_cum_list(cum_arr, cum_sum, ind, p_orb)

            ! all orbitals are possible i guess, so make cum_arr for all 
            ! orbitals as ind already. we "just" have to fix the spin 
            if (is_beta(src)) then 
                orb = 2*ind - 1
            else 
                orb = 2*ind 
            end if
            
            pgen = p_elec * p_orb * (1.0_dp - pDoubles) 

            call make_single(nI, nJ, elec, orb, ex, tParity) 

        end if 

        ilutJ = make_ilutJ(ilutI, ex, ic)

    end subroutine gen_excit_rs_hubbard_transcorr

    subroutine create_cum_list_rs_hubbard_transcorr_single(nI, ilutI, src, &
            cum_arr, cum_sum, tgt, p_orb) 
        ! with transcorrelation use a different cum-list creator, due to 
        ! longer range single excitations possible now. 
        integer, intent(in) :: nI(nel), src 
        integer(n_int), intent(in) :: ilutI(0:NIfTot) 
        real(dp), intent(out) :: cum_arr(nBasis/2), cum_sum
        integer, intent(in), optional :: tgt 
        real(dp), intent(out), optional :: p_orb
#ifdef __DEBUG
        character(*), parameter :: this_routine = "create_cum_list_rs_hubbard_transcorr_single"
#endif
        integer :: spin, ex(2,2), nJ(nel), i, orb
        integer, allocatable :: ex2(:,:)
        real(dp) :: elem
        real(dp) :: temp

        ASSERT(IsOcc(ilutI, src))
        
        cum_arr = 0.0_dp
        cum_sum = 0.0_dp 

        ! 0.. alpha
        ! 1.. beta
        spin = get_spin(src) - 1

        ex = 0
        ex(1,1) = src

        if (present(tgt)) then 
            ASSERT(present(p_orb))
            ASSERT(same_spin(src,tgt))

            p_orb = 0.0_dp

            do i = 1, nBasis/2
                elem = 0.0_dp 

                ! take the same spin
                orb = 2 * i - spin

                ASSERT(same_spin(src,orb))

                if (IsNotOcc(ilutI,orb)) then 
                    
                    ex(2,1) = orb 
                    call swap_excitations(nI, ex, nJ, ex2)
                    elem = abs(get_single_helem_rs_hub_transcorr(nJ, ex2(:,1), .false.))

!                     temp = abs(get_single_helem_rs_hub_transcorr(nI, ex(:,1), .false.))
!                     if (abs(temp - elem) > EPS) then 
!                         print *, "singles do differ!"
!                     end if

                end if
                cum_sum = cum_sum + elem 

                if (tgt == orb) then 
                    p_orb = elem 
                end if
            end do
            if (cum_sum < EPS) then 
                p_orb = 0.0_dp 
            else
                p_orb = p_orb / cum_sum 
            end if
        else 
            do i = 1, nBasis/2 
                elem = 0.0_dp
                orb = 2 * i - spin 

                ASSERT(same_spin(src,orb))

                if (IsNotOcc(ilutI,orb)) then 
                    ex(2,1) = orb 
                    call swap_excitations(nI, ex, nJ, ex2)
                    elem = abs(get_single_helem_rs_hub_transcorr(nJ, ex2(:,1), .false.))
                    ! todo! i am not sure about the order of these matrix 

!                     temp = abs(get_single_helem_rs_hub_transcorr(nI, ex(:,1), .false.))
!                     if (abs(temp - elem) > EPS) then 
!                         print *, "singles do differ!"
!                     end if
                end if

                cum_sum = cum_sum + elem 
                cum_arr(i) = cum_sum
            end do
        end if
        
    end subroutine create_cum_list_rs_hubbard_transcorr_single

    subroutine gen_double_excit_rs_hub_transcorr(nI, ilutI, nJ, ilutJ, ex, tPar, pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2,2) 
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: tPar
        real(dp), intent(out) :: pgen 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "gen_double_excit_rs_hub_transcorr"
#endif
        integer :: elecs(2), orbs(2), src(2), ind
        real(dp) :: p_elec, cum_arr(nBasis/2), cum_sum, p_orb_a, p_orb_b, p_orb_switch

        call pick_spin_opp_elecs(nI, elecs, p_elec) 

        src = nI(elecs) 

        ! pick the first hole at random 
        call pick_random_hole(ilutI, orbs(1), p_orb_a) 

        if (orbs(1) == 0) then 
            nJ(1) = 0 
            pgen = 0.0_dp 
            return 
        end if

        ! create the cum-list for b 
        call create_cum_list_rs_hubbard_transcorr_double(nI, ilutI, src, orbs(1), cum_arr, & 
            cum_sum)

        if (cum_sum < EPS) then 
            nJ(1) = 0
            pgen = 0.0_dp
            return 
        end if

        call pick_from_cum_list(cum_arr, cum_sum, ind, p_orb_b)

        ! pick the right spin 
        if (is_beta(orbs(1))) then 
            orbs(2) = 2*ind 
        else 
            orbs(2) = 2*ind -1
        end if

        ! now pick the other way around 
        call create_cum_list_rs_hubbard_transcorr_double(nI, ilutI, src, orbs(2), cum_arr, &
            cum_sum, orbs(1), p_orb_switch) 

        ! if cum_sum can be 0 here i made something wrong with the cum_sum 
        ! check above! 
        ASSERT(cum_sum > EPS)

        ! no factor of 2, since we add the a|b b|a already!
        pgen = 1.0_dp * p_elec * p_orb_a * (p_orb_b + p_orb_switch) 

        call make_double(nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), ex, tPar)
        ilutJ = make_ilutJ(ilutI, ex, 2) 

    end subroutine gen_double_excit_rs_hub_transcorr

    subroutine create_cum_list_rs_hubbard_transcorr_double(nI, ilutI, src, orb_a, &
            cum_arr, cum_sum, tgt, p_orb) 
        ! routine to create the cum-list to pick the second orbital given 
        ! the electrons (src) and the first orbital (orb_a) 
        ! if second orbital (tgt) is given it calculates the probability (p_orb) 
        ! to have picked this orbital. 
        integer, intent(in) :: nI(nel), src(2), orb_a 
        integer(n_int), intent(in) :: ilutI(0:NIfTot) 
        real(dp), intent(out) :: cum_arr(nBasis/2), cum_sum 
        integer, intent(in), optional :: tgt 
        real(dp), intent(out), optional :: p_orb
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "create_cum_list_rs_hubbard_transcorr_double"
#endif
        integer :: ex(2,2), spin, b, nJ(nel), orb_b
        integer, allocatable :: ex2(:,:)
        real(dp) :: elem
        real(dp) :: temp

        ASSERT(.not. same_spin(src(1),src(2)))
        ASSERT(IsNotOcc(ilutI,orb_a))
        ASSERT(IsOcc(ilutI,src(1)))
        ASSERT(IsOcc(ilutI,src(2)))


        cum_arr = 0.0_dp
        cum_sum = 0.0_dp 

        ! make a spin factor for the orbital conversion
        ! 1...alpha
        ! 2...beta
        spin = get_spin(orb_a)

        ex(1,:) = src 
        ex(2,1) = orb_a

        if (present(tgt)) then 
            ASSERT(present(p_orb)) 
            ASSERT(.not. same_spin(orb_a, tgt))

            p_orb = 0.0_dp 

            do b = 0, nBasis/2 - 1
                elem = 0.0_dp 

                ! add the spin to get correct anti-parallel spin-orbtial to (a)
                ! if (a) is alpha, spin = 1 -> so add 
                orb_b = 2 * b + spin 

                if (IsNotOcc(ilutI, orb_b)) then 
                    ! with an occupancy everything is fine.. since a == b 
                    ! is not possible due to opposite spin 

                    ex(2,2) = orb_b 
                    call swap_excitations(nI, ex, nJ, ex2)
                    elem = abs(get_double_helem_rs_hub_transcorr(nJ, ex2, .false.))

!                     temp = abs(get_double_helem_rs_hub_transcorr(nI, ex, .false.))
!                     if (abs(elem - temp) > EPS) then 
!                         print *, "doubles do differ!"
!                     end if

                end if
                cum_sum = cum_sum + elem 
                
                if (tgt == orb_b) then 
                    p_orb = elem 
                end if
            end do
            if (cum_sum < EPS) then 
                p_orb = 0.0_dp 
            else 
                p_orb = p_orb / cum_sum
            end if
        else
            do b = 0, nBasis/2 - 1

                elem = 0.0_dp 
                orb_b = 2 * b + spin 

                if (IsNotOcc(ilutI, orb_b)) then 

                    ex(2,2) = orb_b 
                    call swap_excitations(nI, ex, nJ, ex2) 
                    elem = abs(get_double_helem_rs_hub_transcorr(nJ, ex2, .false.))

!                     temp = abs(get_double_helem_rs_hub_transcorr(nI, ex, .false.))
!                     if (abs(elem - temp) > EPS) then 
!                         print *, "doubles do differ!"
!                     end if
                end if

                cum_sum = cum_sum + elem 
                cum_arr(b+1) = cum_sum 
            end do
        end if
                    
    end subroutine create_cum_list_rs_hubbard_transcorr_double

    function get_single_helem_rs_hub_transcorr(nI, ex, tPar) result(hel) 
        ! function to get the new 1-body matrix elements in the real-space 
        ! hubbard with hopping transcorelation with influence from the 
        ! newly introduced double excitations 
        ! this is the standalone function outside get_offdiag_helement_rs_hub
        integer, intent(in) :: nI(nel), ex(2)
        logical, intent(in) :: tPar 
        HElement_t(dp) :: hel 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_single_helem_rs_hub_transcorr"
#endif 
        ASSERT(same_spin(ex(1),ex(2)))

        hel = GetTMatEl(ex(1),ex(2)) + get_2_body_contrib_transcorr_hop(nI,ex)

        if (tpar) hel = -hel

    end function get_single_helem_rs_hub_transcorr

    subroutine pick_random_hole(ilutI, orb, p_orb, spin) 
        ! routine to pick an unoccupied spin-orbital (orb) with probability 
        ! p_orb
        ! with an additional spin input we can restrict ourself to a 
        ! specific parallel spin! 
        ! spin = 1 -> beta orbital to be picked 
        ! spin = 0 -> alpha orbital to be picked
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: orb 
        real(dp), intent(out) :: p_orb
        integer, intent(in), optional :: spin 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "pick_random_hole" 
#endif
        integer, parameter :: max_trials = 100
        integer :: cnt, i

        orb = 0
        cnt = 0
        p_orb = 0.0_dp

        if (present(spin)) then 
            ASSERT(spin == 0 .or. spin == 1)
          
            do while (cnt < max_trials)
                ! create a random spin-orbital of parallel spin
                i = 2*(1 + int(genrand_real2_dsfmt() * nBasis/2)) - spin

                if (IsNotOcc(ilutI,i)) then 
                    orb = i
                    if (spin == 0) then 
                        p_orb = 1.0_dp / real(nBasis/2 - nOccAlpha,dp)
                    else if (spin == 1) then 
                        p_orb = 1.0_dp / real(nBasis/2 - nOccBeta,dp)
                    end if
                    return
                end if
            end do

        else

            do while (cnt < max_trials)
                cnt = cnt + 1
                i = 1 + int(genrand_real2_dsfmt() * nBasis)

                if (IsNotOcc(ilutI, i)) then 
                    orb = i 
                    p_orb = 1.0_dp / real(nBasis - nel, dp)
                    return 
                end if
            end do
        end if

    end subroutine pick_random_hole

    function calc_pgen_rs_hubbard_transcorr(nI, ilutI, ex, ic) result(pgen)
        integer, intent(in) :: nI(nel), ex(2,2), ic 
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        real(dp) :: pgen 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "calc_pgen_rs_hubbard_transcorr"
#endif
        integer :: src(2), tgt(2)
        real(dp) :: p_elec, p_orb, cum_arr(nBasis/2), cum_sum, p_hole_1, & 
                    p_orb_a, p_orb_b
        
        src = ex(1,:)
        tgt = ex(2,:)

        if (ic == 1) then 

            ASSERT(same_spin(src(1),tgt(1)))

            p_elec = 1.0_dp / real(nel, dp) 

            call create_cum_list_rs_hubbard_transcorr_single(nI, ilutI, src(1), &
                cum_arr, cum_sum,  tgt(1), p_orb) 

            if (cum_sum < EPS) then 
                pgen = 0.0_dp 
                return 
            end if
            
            pgen = p_elec * p_orb * (1.0_dp - pDoubles) 

        else if (ic == 2) then 

            if (same_spin(src(1),src(2)) .or. same_spin(tgt(1),tgt(2))) then 
                pgen = 0.0_dp
                return
            end if

            p_elec = 1.0_dp / real(nOccAlpha * nOccBeta, dp)

            ! we need two holes.. 
            ! pick the first at random.. 
            p_hole_1 = 1.0_dp / real(nBasis - nel, dp) 

            call create_cum_list_rs_hubbard_transcorr_double(nI, ilutI, src, tgt(1), & 
                cum_arr, cum_sum, tgt(2), p_orb_a) 

            if (cum_sum < EPS) then 
                pgen = 0.0_dp
                return
            end if

            ! and the other way around 
            call create_cum_list_rs_hubbard_transcorr_double(nI, ilutI, src, tgt(2), & 
                cum_arr, cum_sum, tgt(1), p_orb_b) 

            if (cum_sum < EPS) then 
                pgen = 0.0_dp
                return
            end if

            pgen = 1.0_dp * p_elec * p_hole_1 * (p_orb_a + p_orb_b) * pDoubles

        else 
            ! every other ic is 0
            pgen = 0.0_dp

        end if

    end function calc_pgen_rs_hubbard_transcorr

    ! Generic excitaiton generator
    subroutine gen_excit_rs_hubbard (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)

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

        integer :: iunused, ind , elec, src, orb, n_avail, n_orbs, i
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

        ! now get neighbors
        neighbors = lat%get_spinorb_neighbors(src)

        call create_cum_list_rs_hubbard(ilutI, src, neighbors, cum_arr, cum_sum)

        if (cum_sum < EPS) then 
            nJ(1) = 0
            pgen = 0.0_dp 
            return
        end if

        call pick_from_cum_list(cum_arr, cum_sum, ind, p_orb)

        orb = neighbors(ind)

        pgen = p_elec * p_orb

        call make_single(nI, nJ, elec, orb, ex, tParity) 

        ilutJ = make_ilutJ(ilutI, ex, 1)

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

        call create_cum_list_rs_hubbard(ilutI, src, lat%get_spinorb_neighbors(src), &
            cum_arr, cum_sum, tgt, p_orb) 
        if (cum_sum < EPS) then 
            pgen = 0.0_dp 
            return
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
        integer :: i, nI(nel), ex(2), ex2(2), nJ(nel)

        ASSERT(IsOcc(ilutI,src))

        call decode_bit_det(nI, ilutI)

        ex(1) = src 

        allocate(cum_arr(size(neighbors)))
        cum_arr = 0.0_dp
        cum_sum = 0.0_dp
        if (present(tgt)) then 
            ASSERT(present(cpt))
            do i = 1, ubound(neighbors,1)
                elem = 0.0_dp
                ASSERT(is_beta(src) .eqv. is_beta(neighbors(i)))
                if (IsNotOcc(ilutI, neighbors(i))) then 
                    ! change the order of determinants to reflect 
                    ! non-hermiticity correctly
                    ! old implo:
!                     elem = abs(get_offdiag_helement_rs_hub(nI,[src,neighbors(i)],.false.))
                    ex(2) = neighbors(i) 
                    call swap_excitations(nI, ex, nJ, ex2) 
                    elem = abs(get_offdiag_helement_rs_hub(nJ,ex2,.false.))
 
                end if
                cum_sum = cum_sum + elem
                if (neighbors(i) == tgt) then 
                    cpt = elem
                end if
            end do
            if (cum_sum < EPS) then
                cpt = 0.0
            else 
                cpt = cpt / cum_sum
            end if
        else
            do i = 1, ubound(neighbors,1)
                elem = 0.0_dp
                ASSERT(is_beta(src) .eqv. is_beta(neighbors(i)))
                if (IsNotOcc(ilutI,neighbors(i))) then 
                    ex(2) = neighbors(i) 
                    call swap_excitations(nI, ex, nJ, ex2) 
                    elem = abs(get_offdiag_helement_rs_hub(nJ,ex2,.false.))
                end if
                cum_sum = cum_sum + elem
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
        integer, intent(in) :: ic, ex(2,ic)
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

        else if (ic == 2 .and. t_trans_corr_hop) then 
            hel = get_double_helem_rs_hub_transcorr(nI, ex, tpar)

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

            else if (ic_ret == 2 .and. t_trans_corr_hop) then 

                ex(1,1) = 2 
                call GetExcitation(nI, nI, nel, ex, tpar)
                hel = get_double_helem_rs_hub_transcorr(nI, ex, tpar)

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

                else if (ic_ret == 2 .and. t_trans_corr_hop) then 
                    ex(1,1) = 2
                    call GetBitExcitation(ilutI, ilutJ, ex, tpar)

                    hel = get_double_helem_rs_hub_transcorr(nI, ex, tpar)

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

            else if (ic == 2 .and. t_trans_corr_hop) then 
                ex(1,1) = 2
                call GetBitExcitation(ilutI, ilutJ, ex, tpar)
                hel = get_double_helem_rs_hub_transcorr(nI, ex, tpar)
            else
                hel = h_cast(0.0_dp)
            end if
        end if
    
    end function get_helement_rs_hub_general

    ! also optimize the matrix element calculation
    function get_diag_helemen_rs_hub(nI) result(hel)
        integer, intent(in) :: nI(nel)
        HElement_t(dp) :: hel

        integer(n_int) :: ilut(0:NIfTot)
        ! the diagonal matrix element is essentialy just the number of 
        ! doubly occupied sites times U

        call EncodeBitDet(nI, ilut)

        if (t_trans_corr_hop) then 
            hel = get_diag_helemen_rs_hub_transcorr_hop(nI)
        else
            hel = h_cast(uhub * count_double_orbs(ilut))
        end if

    end function get_diag_helemen_rs_hub

    function get_diag_helemen_rs_hub_transcorr_hop(nI) result(hel)
        ! with the hopping transcorrelation also the diagonal matrix 
        ! elements change! 
        integer, intent(in) :: nI(nel)
        HElement_t(dp) :: hel

        integer :: i, j, id(nel), idX, idN

        hel = 0.0_dp 

        id = get_spatial(nI)

        ! now also n_iu n_jd contribute to the diagonal elements! 
        do i = 1, nel 
            do j = 1, nel 
                if (.not. same_spin(nI(i), nI(j))) then 

!                     idX = max(id(i),id(j))
!                     idN = min(id(i),id(j))
!                     ! is this sufficient? do i have a factor of 1/2?
!                     hel = hel + 0.5_dp * get_umat_rs_hub_trans(idN,idX,idN,idX)

                    hel = hel + 0.5_dp * get_umat_rs_hub_trans(id(i),id(j),id(j),id(i))

                end if
            end do
        end do

    end function get_diag_helemen_rs_hub_transcorr_hop

    function get_double_helem_rs_hub_transcorr(nI, ex, tpar) result(hel)
        ! newly introduced 2-body matrix element in the hopping 
        ! transcorrelated real-space hubbard model 
        integer, intent(in) :: nI(nel), ex(2,2) 
        logical, intent(in) :: tpar 
        HElement_t(dp) :: hel
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_double_helem_rs_hub_transcorr"
#endif
        integer :: src(2), tgt(2), ij(2), ab(2)
    
        ASSERT(t_trans_corr_hop) 

        if (same_spin(ex(1,1),ex(1,2)) .or. &
            same_spin(ex(2,1),ex(2,2))) then 
            hel = 0.0_dp
            return
        end if

        src = get_src(ex)
        tgt = get_tgt(ex)

        ij = get_spatial(src)
        ab = get_spatial(tgt)

        ! this weird combined sign convention.. 
        ! NOTE: i have to be careful due to the non-hermitian hamiltonian.. 
        ! i guess it matters now where i put the holes and electrons!
        ! and maybe it even matters where i put the indices for electrons 
        ! and their corresponding fitting hole.. to be tested! 
        if (same_spin(src(1),tgt(1)) .and. same_spin(src(2),tgt(2))) then 

            hel = get_umat_rs_hub_trans(ij(1),ij(2),ab(1),ab(2))

        else if (same_spin(src(1),tgt(2)) .and. same_spin(src(2),tgt(1))) then 

            hel = -get_umat_rs_hub_trans(ij(1),ij(2),ab(1),ab(2))
            
        end if

        if (tpar) hel = -hel

    end function get_double_helem_rs_hub_transcorr

    function get_offdiag_helement_rs_hub(nI, ex, tpar) result(hel)
        integer, intent(in) :: nI(nel), ex(2)
        logical, intent(in) :: tpar
        HElement_t(dp) :: hel

        integer(n_int) :: ilut(0:NIfTot)
        real(dp) :: n_i, n_j

        ! in case we need it, the off-diagonal, except parity is just 
        ! -t if the hop is possible
        hel = GetTMatEl(ex(1),ex(2))

        if (t_trans_corr_hop) then 
            hel = hel + get_2_body_contrib_transcorr_hop(nI, ex)
        end if

        if (tpar) hel = -hel

        ! put the transcorrelated stuff here for now, alhtough it would be 
        ! better to do it as a procedure pointer.. 
        if (t_trans_corr) then 
            call EncodeBitDet(nI, ilut)
            if (t_trans_corr_new) then
                n_j = get_opp_spin(ilut, ex(1))
                n_i = get_opp_spin(ilut, ex(2))

                ! try to go one step further before the transformation to 
                ! the momentum space.. and check if the results are still 
                ! correct.. 

                hel = hel * (1.0_dp + (exp(trans_corr_param) - 1.0_dp) * n_j + & 
                                      (exp(-trans_corr_param)- 1.0_dp) * n_i - & 
                                       2.0_dp*(cosh(trans_corr_param) - 1.0_dp)*n_i*n_j)

           else

                hel = hel * exp(trans_corr_param * & 
                    (get_opp_spin(ilut, ex(1)) - get_opp_spin(ilut, ex(2))))
            end if

        end if

        ! it is not really a 2-body trans, but still use the flag and 
        ! parameter
        if (t_trans_corr_2body) then 
            call EncodeBitDet(nI, ilut) 
            hel = hel * exp(trans_corr_param_2body * & 
                (get_spin_opp_neighbors(ilut, ex(1)) - get_spin_opp_neighbors(ilut,ex(2))))
        end if
    end function get_offdiag_helement_rs_hub

    function get_2_body_contrib_transcorr_hop(nI, ex) result(hel) 
        ! new single excitation matrix element calculation 
        ! in the hopping transcorrelation this has influence from the 
        ! 2-body term now. the original 1-body term is already calulated 
        ! outside this function
        integer, intent(in) :: nI(nel), ex(2)
        HElement_t(dp) :: hel
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_2_body_contrib_transcorr_hop"
#endif
        integer :: i, idX(2), idN

        ASSERT(same_spin(ex(1),ex(2)))

        hel = 0.0_dp 

        idX = get_spatial(ex)

        ! i have to loop over the occupied sites from the opposite spin 
        ! and retrieve the correct integral
        if (is_beta(ex(1))) then 
            do i = 1, nel 
                if (is_alpha(nI(i))) then 
                    ! get the correct indices. 
                    ! NOTE: i also have to think about the hermiticity here! 
                    ! so the order in umat counts i guess!!! also important 
                    ! for the setup of umat! 
                    idN = get_spatial(nI(i)) 
                    hel = hel + get_umat_rs_hub_trans(idX(1),idN,idN,idX(2))
!                     hel = hel + get_umat_rs_hub_trans(idX(1),idN,idX(2),idN)
                end if
            end do
        else 
            do i = 1, nel 
                if (is_beta(nI(i))) then 
                    idN = get_spatial(nI(i)) 
                    hel = hel + get_umat_rs_hub_trans(idX(1),idN,idN,idX(2))
!                     hel = hel + get_umat_rs_hub_trans(idX(1),idN,idX(2),idN)
                end if
            end do
        end if

    end function get_2_body_contrib_transcorr_hop

    function get_opp_spin(ilut, spin_orb) result(opp_spin) 
        integer(n_int), intent(in) :: ilut(0:NIfTot) 
        integer, intent(in) :: spin_orb
        real(dp) :: opp_spin 

        opp_spin = 0.0_dp

        if (is_beta(spin_orb)) then 
            if (IsOcc(ilut, spin_orb + 1)) opp_spin = 1.0_dp 
        else 
            if (IsOcc(ilut, spin_orb -1 )) opp_spin = 1.0_dp
        end if

    end function get_opp_spin

    function get_spin_opp_neighbors(ilut, spin_orb) result(spin_opp_neighbors) 
        ! function to give the number of opposite spin electron neighbors 
        integer(n_int), intent(in) :: ilut(0:niftot) 
        integer, intent(in) :: spin_orb 
        real(dp) :: spin_opp_neighbors
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_spin_opp_neighbors" 
#endif 
        integer :: i, flip 
        integer, allocatable :: neighbors(:) 

        ASSERT(associated(lat))

        spin_opp_neighbors = 0.0_dp 


        ! get the spin-opposite neigbhors 
        if (is_beta(spin_orb)) then 
            neighbors = lat%get_spinorb_neighbors(spin_orb) + 1
        else 
            neighbors = lat%get_spinorb_neighbors(spin_orb) - 1
        end if

        do i = 1, size(neighbors) 
            if (IsOcc(ilut, neighbors(i))) spin_opp_neighbors = spin_opp_neighbors + 1.0_dp
        end do

    end function get_spin_opp_neighbors


    ! what else?
    function create_neel_state(ilut_neel) result(neel_state)
        ! probably a good idea to have a routine which creates a neel state
        ! (or one of them if it is ambigous)
        integer(n_int), intent(out), optional :: ilut_neel(0:NIfTot)
        integer :: neel_state(nel)
        character(*), parameter :: this_routine = "create_neel_state"

        integer :: i, j, k, l, spin, ind
        

        if (tHPHF) then 
            call stop_all(this_routine, "remember to change for HPHF!")
        end if
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

        case ('square','rectangle','triangle','triangular','hexagonal','kagome')
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

        case default 
            call stop_all(this_routine, "unknown lattice type!")
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

    function get_umat_rs_hub_trans(i,j,k,l) result(hel) 
        ! do i need an explicit get_umat_rs_hub_trans? or can i just reuse 
        ! the one, whhich access the "normal" fcidump.. figure out! 
        integer, intent(in) :: i,j,k,l
        HElement_t(dp) :: hel 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_umat_rs_hub_trans"
#endif 

        hel = umat_rs_hub_trancorr_hop(i,j,k,l)

    end function get_umat_rs_hub_trans

    subroutine swap_excitations_higher(nI, ex, nJ, ex2) 
        ! routine to quickly, without make_double and make_triple 
        ! create the excited determinant nJ and ex2 to go from nJ -> nI 
        ! for the test of the order of matrix element calculation in the 
        ! transcorrelated approach. due to non-hermiticity 
        ! <i|H|j> /= <j|H|i> 
        ! also make this work for triple excitations! 
        integer, intent(in) :: nI(nel), ex(:,:) 
        integer, intent(out) :: nJ(nel) 
        integer, intent(out), allocatable :: ex2(:,:)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "swap_excitations" 
#endif

        ASSERT(size(ex,1) == 2) 
        ASSERT(size(ex,2) == 2 .or. size(ex,2) == 3)

        allocate(ex2(size(ex,1),size(ex,2)), source = ex)

        nJ = nI
        where (nJ == ex(1,1)) nJ = ex(2,1)
        where (nJ == ex(1,2)) nJ = ex(2,2)
        if (size(ex,2) == 3) then 
            where (nJ == ex(1,3)) nJ = ex(2,3)
        end if

        call sort(nJ)
        ex2 = ex 
        call swap(ex2(1,1),ex2(2,1))
        call swap(ex2(1,2),ex2(2,2))
        if (size(ex,2) == 3) then 
            call swap(ex2(1,3),ex2(2,3))
        end if

    end subroutine swap_excitations_higher

    subroutine swap_excitations_singles(nI, ex, nJ, ex2)
        integer, intent(in) :: nI(nel), ex(2) 
        integer, intent(out) :: nJ(nel), ex2(2)
        
        nJ = nI 

        where (nJ == ex(1)) nJ = ex(2)

        call sort(nJ) 

        ex2 = ex

        call swap(ex2(1),ex2(2))

    end subroutine swap_excitations_singles

    subroutine pick_spin_opp_elecs(nI, elecs, p_elec) 
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: elecs(2)
        real(dp), intent(out) :: p_elec
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_spin_opp_elecs"
#endif
        ! think of a routine to get the possible spin-opposite electron 
        ! pairs. i think i could do that way more efficiently, but do it in 
        ! the simple loop way for now 
        do 
            elecs(1) = 1 + int(genrand_real2_dsfmt() * nel)
            do 
                elecs(2) = 1 + int(genrand_real2_dsfmt() * nel) 

                if (elecs(1) /= elecs(2)) exit 

            end do
            if (get_ispn(nI(elecs)) == 2) exit
        end do

        ! output in ordered form 
        elecs = [minval(elecs),maxval(elecs)]

        ! actually the probability is twice that or? 
        ! or doesnt that matter, since it is the same
        p_elec = 1.0_dp / real(nOccBeta * nOccAlpha, dp)

    end subroutine pick_spin_opp_elecs

    subroutine pick_from_cum_list(cum_arr, cum_sum, ind, pgen) 
        real(dp), intent(in) :: cum_arr(:), cum_sum
        integer, intent(out) :: ind
        real(dp), intent(out) :: pgen 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_from_cum_list" 
#endif
        real(dp) :: r

        if (cum_sum < EPS) then 
            ind = -1 
            pgen = 0.0_dp
            return 
        end if

        r = genrand_real2_dsfmt() * cum_sum

        ind = binary_search_first_ge(cum_arr, r) 

        if (ind == 1) then 
            pgen = cum_arr(1)/cum_sum 
        else 
            pgen = (cum_arr(ind) - cum_arr(ind - 1)) / cum_sum
        end if

    end subroutine pick_from_cum_list

    subroutine gen_all_excits_r_space_hubbard(nI, n_excits, det_list) 
        ! for the purpose of excitation generation and time-step and 
        ! pDoubles determination create a routine to create all possible 
        ! excitations for a given determinant nI 
        ! the system specific variables need to be determined already 
        ! before! 
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:) 
        character(*), parameter :: this_routine = "gen_all_excits_r_space_hubbard"

        integer :: save_excits, n_doubles
        integer(n_int), allocatable :: double_dets(:,:), temp_dets(:,:)

        if (tHPHF) then 
            call stop_all(this_routine, &
                "not adapted to HPHF, since we create all the excitations!")
        end if

        call gen_all_singles_rs_hub(nI, n_excits, det_list) 

        ! if we have hopping transcorr we also have to compute the 
        ! possible doubles 
        ! we know does are different than any single
        if (t_trans_corr_hop) then 
            save_excits = n_excits

            call gen_all_doubles_rs_hub_hop_transcorr(nI, n_doubles, double_dets) 

            n_excits = n_excits + n_doubles

            allocate(temp_dets(0:niftot, save_excits), source = det_list(:,1:save_excits))

            deallocate(det_list) 

            allocate(det_list(0:niftot,n_excits)) 

            det_list(:,1:save_excits) = temp_dets 

            det_list(:,save_excits+1:n_excits) = double_dets

        end if

        call sort(det_list, ilut_lt, ilut_gt)

    end subroutine gen_all_excits_r_space_hubbard

    subroutine gen_all_doubles_rs_hub_hop_transcorr(nI, n_excits, det_list) 
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:) 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "gen_all_doubles_rs_hub_hop_transcorr"
#endif
        integer(n_int) :: ilut(0:NIfTot), ilutJ(0:NIfTot)
        integer :: n_bound, i, src(2), j, neigh, ex(2,2), pos, a, b
        integer(n_int), allocatable :: temp_list(:,:)
        integer, allocatable :: neighbors(:)

        call EncodeBitDet(nI, ilut) 

        n_excits = 1
        
        n_bound = nel*(nel-1)*(nbasis-nel)*(nbasis-nel-1)

        allocate(temp_list(0:NIfTot,n_bound))
        temp_list = 0_n_int

        ! we just have opposite spin doubles! 
        do i = 1, nel - 1
            do j = i + 1, nel 

                src = nI([i,j])

                if (same_spin(src(1),src(2))) cycle 

                ex(1,:) = src 
                ! and now loop over all possible empty spin opposite 
                ! orbitals in a unique way.. 
                do a = 1, nbasis - 1
                    if (IsOcc(ilut,a)) cycle 
                    do b = a + 1, nbasis 
                        if (IsOcc(ilut,b) .or. same_spin(a,b)) cycle 

                        ex(2,:) = [a,b] 

                        if (abs(get_double_helem_rs_hub_transcorr(nI, ex, .false.)) > EPS) then

                            ilutJ = make_ilutJ(ilut, ex, 2)
                            ! actually a search is not really necessary.. since 
                            ! all the single excitations are unique.. but 
                            ! just to be sure
                            pos = binary_search(temp_list(:,1:n_excits), ilutJ, nifd+1)

                            if (pos < 0) then 

                                temp_list(:,n_excits) = ilutJ
                                call sort(temp_list(:,1:n_excits),ilut_lt,ilut_gt)
                                n_excits = n_excits + 1
                                ! damn.. i have to sort everytime i guess..
                            end if
                        end if
                    end do
                end do
            end do
        end do

        n_excits = n_excits - 1
        allocate(det_list(0:NIfTot,n_excits), source = temp_list(:,1:n_excits))

        call sort(det_list, ilut_lt, ilut_gt)

    end subroutine gen_all_doubles_rs_hub_hop_transcorr

    subroutine gen_all_singles_rs_hub(nI, n_excits, det_list)
        ! create all single excitations in the real-space hubbard 
        ! without hopping transcorrelation this is quite easy.. 
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:) 
        character(*), parameter :: this_routine = "gen_all_singles_rs_hub" 

        if (t_trans_corr_hop) then 
            call gen_all_singles_rs_hub_hop_transcorr(nI, n_excits, det_list)
        else 
            call gen_all_singles_rs_hub_default(nI, n_excits, det_list) 
        end if

    end subroutine gen_all_singles_rs_hub

    subroutine gen_all_singles_rs_hub_hop_transcorr(nI, n_excits, det_list)
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:) 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_all_singles_rs_hub_hop_transcorr"
#endif 
        integer(n_int) :: ilut(0:NIfTot), ilutJ(0:NIfTot)
        integer :: n_bound, i, src, j, neigh, ex(2), pos, a
        integer(n_int), allocatable :: temp_list(:,:)
        integer, allocatable :: neighbors(:)
        real(dp) :: elem

        call EncodeBitDet(nI, ilut) 

        n_excits = 1

        n_bound = nel * (nbasis - nel)
        allocate(temp_list(0:NIfTot,n_bound))
        temp_list = 0_n_int
        ! loop over electrons 
        do i = 1, nel 
            ! but now loop over all orbitals 
            src = nI(i)
            ex(1) = src 

            do a = 1, nbasis

                ! only same-spin excitations possible
                if (same_spin(src,a) .and. IsNotOcc(ilut,a)) then 
                    ex(2) = a
                    elem = abs(get_single_helem_rs_hub_transcorr(nI, ex, .false.))

                    if (elem > EPS) then 

                        ilutJ = make_ilutJ(ilut, ex, 1)

                        ! actually a search is not really necessary.. since 
                        ! all the single excitations are unique.. but 
                        ! just to be sure
                        pos = binary_search(temp_list(:,1:n_excits), ilutJ, nifd+1)

                        if (pos < 0) then 

                            temp_list(:,n_excits) = ilutJ
                            call sort(temp_list(:,1:n_excits),ilut_lt,ilut_gt)
                            n_excits = n_excits + 1
                            ! damn.. i have to sort everytime i guess..
                        end if
                    end if
                end if
            end do
        end do

        n_excits = n_excits - 1
        allocate(det_list(0:NIfTot,n_excits), source = temp_list(:,1:n_excits))

        call sort(det_list, ilut_lt, ilut_gt)

    end subroutine gen_all_singles_rs_hub_hop_transcorr

    subroutine gen_all_singles_rs_hub_default(nI, n_excits, det_list) 
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:) 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_all_singles_rs_hub_default" 
#endif
        integer(n_int) :: ilut(0:NIfTot), ilutJ(0:NIfTot)
        integer :: n_bound, i, src, j, neigh, ex(2), pos
        integer(n_int), allocatable :: temp_list(:,:)
        integer, allocatable :: neighbors(:)
        real(dp) :: elem

        ASSERT(associated(lat))

        call EncodeBitDet(nI, ilut) 

        n_excits = 1

        n_bound = nel * (nbasis - nel)
        allocate(temp_list(0:NIfTot,n_bound))
        temp_list = 0_n_int

        ! loop over electrons
        do i = 1, nel 
            src = nI(i)

            ex(1) = src
            neighbors = lat%get_spinorb_neighbors(src)

            ! and loop over the neighboring sites of this electron
            do j = 1, size(neighbors)

                neigh = neighbors(j)

                ex(2) = neigh

                ASSERT(same_spin(src,neigh))

                ! if it is not occupied it should be a possible excitation
                if (IsNotOcc(ilut,neighbors(j))) then 
                    ! but to be sure, check the matrix element: 
                    elem = abs(get_offdiag_helement_rs_hub(nI,ex,.false.))

                    if (elem > EPS) then 

                        ilutJ = make_ilutJ(ilut, ex, 1)

                        ! actually a search is not really necessary.. since 
                        ! all the single excitations are unique.. but 
                        ! just to be sure
                        pos = binary_search(temp_list(:,1:n_excits), ilutJ, nifd+1)

                        if (pos < 0) then 

                            temp_list(:,n_excits) = ilutJ
                            call sort(temp_list(:,1:n_excits),ilut_lt,ilut_gt)
                            n_excits = n_excits + 1
                            ! damn.. i have to sort everytime i guess..
                        end if
                    end if
                end if
            end do
        end do

        n_excits = n_excits - 1
        allocate(det_list(0:NIfTot,n_excits), source = temp_list(:,1:n_excits))

        call sort(det_list, ilut_lt, ilut_gt)

    end subroutine gen_all_singles_rs_hub_default

end module real_space_hubbard
        
        
