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
                        t_fill_frequency_hists, p_singles_input, matele_cutoff

    use dsfmt_interface, only: genrand_real2_dsfmt

    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet, ilut_lt, ilut_gt

    use bit_rep_data, only: NIfTot, nifd

    use util_mod, only: binary_search_first_ge, choose, swap, get_free_unit, &
                        binary_search

    use bit_reps, only: decode_bit_det

    use sort_mod, only: sort

    use get_excit, only: make_double, make_single

    use double_occ_mod, only: count_double_orbs

    use lattice_models_utils, only: swap_excitations, pick_spin_opp_elecs, & 
                                    pick_from_cum_list, pick_spin_opp_holes, & 
                                    pick_random_hole, get_opp_spin, &
                                    get_spin_opp_neighbors, create_neel_state, &
                                    make_ilutJ, get_ispn

    implicit none 

    real(dp) :: lat_tau_factor = 0.5_dp

    ! create a flag which indicate to start in a neel state 
    logical :: t_start_neel_state = .false. 

    real(dp), allocatable :: umat_rs_hub_trancorr_hop(:,:,:,:)

    complex(dp), parameter :: imag_unit = cmplx(0.0_dp,1.0_dp)

    interface get_helement_rs_hub
        module procedure get_helement_rs_hub_ex_mat
        module procedure get_helement_rs_hub_general
    end interface get_helement_rs_hub

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
        if (t_trans_corr_hop .and. .not. tHPHF) then 
!             ASSERT(.not. tHPHF)
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

        root_print "initializing UMAT:"
        do i = 1, n_sites
            do j = 1, n_sites
                do k = 1, n_sites
                    do l = 1, n_sites
                        elem = uhub * sum_hop_transcorr_factor(i,j,k,l)
                        ! write to the dumpfile
                        if (abs(elem) > matele_cutoff) then 
                            write(iunit,*)  i,j,k,l, elem
                            ! and also store in the umat 
                            umat_rs_hub_trancorr_hop(i,j,k,l) = elem 
                        end if
                    end do
                end do
            end do
        end do
        close(iunit)
        root_print "Done"

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

            ! to be compatible with my test-suite actually calculate the 
            ! matrix element here.. 
            if (abs(get_double_helem_rs_hub_transcorr(nI, ex, .false.)) < EPS) then 
                nJ(1) = 0
                pgen = 0.0_dp
                return 
            end if
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

            if (abs(get_single_helem_rs_hub_transcorr(nI, ex, .false.)) < EPS) then 
                nJ(1) = 0
                pgen = 0.0_dp 
                return
            end if
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
                    
                    ! i am still not sure about the ordering of these weights..
                    ex(2,1) = orb 
                    call swap_excitations(nI, ex, nJ, ex2)
                    elem = abs(get_single_helem_rs_hub_transcorr(nJ, ex2(:,1), .false.))

!                     elem = abs(get_single_helem_rs_hub_transcorr(nI, ex(:,1), .false.))
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

!                     elem = abs(get_single_helem_rs_hub_transcorr(nI, ex(:,1), .false.))
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

!                     elem = abs(get_double_helem_rs_hub_transcorr(nI, ex, .false.))
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

!                     elem = abs(get_double_helem_rs_hub_transcorr(nI, ex, .false.))
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

end module real_space_hubbard
        
        
