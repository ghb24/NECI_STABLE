#include "macros.h" 

! finally write a specific routine for the k-space hubbard, where every 
! necessary functionality for the k-space/momentum space hubbard is brought 
! together. since i want to implement a transcorrelation factor for the 
! k-space hubbard, which would necessitate a cumulative sum exciation creation 
! i decided it is better to start a new branch, instead of hacking into the 
! old branches of the code. 
! in the end i also want to combine this fully with the lattice_mod implementation
! so i can easily decide which lattice to choose and then decide if we want to 
! use a real or momentum space basis (and in the future maybe even wavelets) 
module k_space_hubbard 
    use SystemData, only: t_lattice_model, t_k_space_hubbard, t_trans_corr, & 
                    trans_corr_param, t_trans_corr_2body, trans_corr_param_2body, & 
                    nel, tHPHF, nOccBeta, nOccAlpha, nbasis, tLatticeGens, tHub, &
                    omega, bhub, nBasisMax, G1, BasisFN, NullBasisFn, TSPINPOLAR, & 
                    treal, ttilt, tExch
    use lattice_mod, only: get_helement_lattice_ex_mat, get_helement_lattice_general, &
                           determine_optimal_time_step, lattice, sort_unique, lat
    use procedure_pointers, only: get_umat_el, generate_excitation
    use gen_coul_ueg_mod, only: get_hub_umat_el
    use constants, only: n_int, dp, EPS, bits_n_int
    use bit_rep_data, only: NIfTot
    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet
    use real_space_hubbard, only: lat_tau_factor
    use fcimcdata, only: tsearchtau, tsearchtauoption, pDoubles, pParallel
    use CalcData, only: tau, t_hist_tau_search, t_hist_tau_search_option
    use dsfmt_interface, only: genrand_real2_dsfmt
    use util_mod, only: binary_search_first_ge
    use back_spawn, only: make_ilutJ, get_orb_from_kpoints, is_allowed_ueg_k_vector, &
                          get_ispn
    use FciMCData, only: excit_gen_store_type
    use get_excit, only: make_double
    use UmatCache, only: gtid
    use OneEInts, only: GetTMatEl, tmat2d
    use sltcnd_mod, only: sltcnd_0
    use sym_mod, only: RoundSym, AddElecSym, SetupSym, lChkSym, mompbcsym, & 
                       TotSymRep, GenMolpSymTable
    use SymExcitDataMod, only: KPointToBasisFn

    implicit none 

    integer, parameter :: ABORT_EXCITATION = 0
    integer, parameter :: N_DIM = 3

    real(dp) :: p_triples = 0.0_dp 
    real(dp) :: three_body_prefac = 0.0_dp
    ! i especially need an interface for the matrix element calculation to 
    ! implement the transcorrelated hamiltonian 
    interface get_helement_k_space_hub
        module procedure get_helement_k_space_hub_ex_mat
        module procedure get_helement_k_space_hub_general
    end interface get_helement_k_space_hub

contains 

    subroutine init_k_space_hubbard() 

        real(dp) :: tau_opt
        print *, " new k-space hubbard implementation init:" 

        ! i have to set the incorrect excitaiton generator flags to false 
        tLatticeGens = .false.
        ! maybe i also need to turn off the hubbard keyword.. at this 
        ! point
        thub = .false. 

        call check_k_space_hubbard_input()

        get_umat_el => get_hub_umat_el

        call init_get_helement_k_space_hub()

        if (.not. tHPHF) generate_excitation => gen_excit_k_space_hub
        tau_opt = determine_optimal_time_step() 

        if (tau < EPS) then 
            print *, "setting time-step to optimally determined time-step: ", tau_opt
            print *, "times: ", lat_tau_factor
            tau = lat_tau_factor * tau_opt

        else 
            print *, "optimal time-step would be: ", tau_opt
            print *, "but tau specified in input!"
        end if

        tsearchtau = .false. 
        tsearchtauoption = .true.

        ! can i set exchange to false, if we have no transcorrelation
        ! although.. not sure..
!         if (.not. t_trans_corr_2body) tExch = .false. 
        t_hist_tau_search = .false. 
        t_hist_tau_search_option = .false.

        if (t_trans_corr_2body) then 
            three_body_prefac = 2.0_dp * (cosh(trans_corr_param_2body) - 1.0_dp) / real(omega**2,dp)
            ! i also have to set some generation probability parameters.. 
            p_triples = 1.0_dp - pDoubles
            pParallel = 0.2_dp
        end if
    end subroutine init_k_space_hubbard

    subroutine check_k_space_hubbard_input()
        character(*), parameter :: this_routine = "check_k_space_hubbard_input"

        print *, "checking input for k-space hubbard:" 
        !todo: find the incompatible input and abort here!

        print *, "input is fine!"

    end subroutine check_k_space_hubbard_input

    subroutine gen_excit_k_space_hub (nI, ilutI, nJ, ilutJ, exFlag, ic, &
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
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_excit_k_space_hub"
#endif
        real(dp) :: p_elec, p_orb
        integer :: elecs(2), orbs(2), src(2)

        ! i first have to choose an electron pair (ij) at random 
        ! but with the condition that they have to have opposite spin! 
        call pick_spin_opp_elecs(nI, elecs, p_elec) 

        src = nI(elecs)

        call pick_ab_orbitals_hubbard(nI, ilutI, src, orbs, p_orb)

        if (orbs(1) == ABORT_EXCITATION) then 
            nJ(1) = ABORT_EXCITATION
            pgen = 0.0_dp
            return 
        end if

        ic = 2 

        ! and make the excitation 
        call make_double(nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), ex, tParity)

        ilutJ = make_ilutJ(ilutI, ex, 2) 

        ! i think in both the electrons and the orbitals i have twice the 
        ! probability to pick them
        pgen = p_elec * p_orb * 2.0_dp

    end subroutine gen_excit_k_space_hub

    subroutine gen_excit_k_space_hub_transcorr (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)


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
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_excit_k_space_hub_transcorr"
#endif
        integer :: temp_ex(2,3) 

        if (genrand_real2_dsfmt() < pDoubles) then 
            if (genrand_real2_dsfmt() < pParallel) then 
                ! do a parallel triple excitation, coming from the triples..
                call gen_parallel_double_hubbard(nI, ilutI, nJ, ilutJ, ex, tParity, pgen)
                ic = 2
                pgen = pgen * pDoubles * pParallel
            else 
                ! do a "normal" hubbard k-space excitation 
                call gen_excit_k_space_hub (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)

                pgen = pgen * pDoubles * (1.0_dp - pParallel)
            end if 
        else 
            ! otherwise to a triple.. 
            call gen_triple_hubbard(nI, ilutI, nJ, ilutJ, temp_ex, tParity, pgen) 
            ic = 3 
            pgen = pgen * (1.0_dp - pDoubles)

        end if

    end subroutine gen_excit_k_space_hub_transcorr

    subroutine gen_triple_hubbard(nI, ilutI, nJ, ilutJ, ex, tParity, pgen) 
        ! i think i should calculat the matrix element in here already! 
        ! in this case.. otherwise i have to carry the tParity and ex 
        ! all the way through the rest of the code and this makes problems 
        ! i guess, since triples are actually never considered .. todo
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ex(2,3) 
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity 
        real(dp), intent(out) :: pgen 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_triple_hubbard" 
#endif
        integer :: elecs(3), orbs(3), src(3), sum_ms
        real(dp) :: p_elec, p_orb(2)

        ! first we pick 3 electrons in this case ofc.
        ! with the restriction, that they must not be all the same spin! 
        call pick_three_opp_elecs(nI, elecs, p_elec, sum_ms)

        src = nI(elecs) 
        ! then i pick 1 orbital? maybe.. 
        ! if this orbital is of the minority spin type, the other 2 orbitals 
        ! MUST be of parallel spin.. 
        ! if the orbital is of the majority spin type the remaining 
        ! orbitals must be of opposite spin type.. 
        ! can i take that into account with a tailored get_orb_from_kpoints() 
        ! function for 3 electrons or should i decide here if we always 
        ! want to do a specific picking order (restricting pgens, but making 
        ! it easier to handle algorithmically) or if we want full flexibility 
        ! (increasing pgens, but kind of making it a bit more difficult..) 
        ! NO: we decide to always pick the minority spin first! 
        call pick_a_orbital_hubbard(ilutI, orbs(1), p_orb(1), sum_ms) 

        ! and pick the remaining two orbitals (essentially it is only 
        ! one, since the last is restricted due to momentum conservation!)
        call pick_bc_orbitals_hubbard(nI, ilutI, src, orbs(1), orbs(2:3), p_orb(2))
        
        if (orbs(2) == 0) then 
            nJ(1) = ABORT_EXCITATION
            pgen = 0.0_dp
            return 
        end if

        ! so now.. did Robert make a routine like: 
        call make_triple(nI, nJ, elecs, orbs, ex, tParity)

        ilutJ = make_ilutJ(ilutI, ex, 3) 

        pgen = p_elec * product(p_orb) * 2.0_dp
        
    end subroutine gen_triple_hubbard

    subroutine pick_bc_orbitals_hubbard(nI, ilutI, src, orb_a, orbs, p_orb)
        ! this is the main routine, which picks the final (b,c) orbital for 
        ! the 3-body excitation
        integer, intent(in) :: nI(nel), src(3), orb_a
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: orbs(2)
        real(dp), intent(out) :: p_orb
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_bc_orbitals_hubbard"
        real(dp) :: test
#endif
        real(dp) :: cum_arr(nbasis/2), cum_sum
        integer :: orb_list(nbasis/2, 2), ind

        ! do it similar to the other routines.. 
        call create_bc_list_hubbard(nI, ilutI, src, orb_a, orb_list, cum_arr, cum_sum)

        if (cum_sum < EPS) then 
            orbs(1) = ABORT_EXCITATION
            return 
        end if

        call pick_from_cum_list(cum_arr, cum_sum, ind, p_orb) 

        orbs = orb_list(ind,:)

#ifdef __DEBUG 
        ! the influence of orb_a is important in the pgen recalc!!
        call create_bc_list_hubbard(nI, ilutI, src, orb_a, orb_list, cum_arr, cum_sum, & 
            orbs(2), test)

        if (abs(test - p_orb) > 1.e-8) then 
            print *, "pgen assumption wrong:!" 
            print *, "p_orb: ", p_orb
            print *, "test: ", test 
            print *, "orbs: ", orbs
        end if
#endif

    end subroutine pick_bc_orbitals_hubbard

    subroutine create_bc_list_hubbard(nI, ilutI, src, orb_a, orb_list, cum_arr, & 
                cum_sum, tgt, cpt) 
        integer, intent(in) :: nI(nel), src(3), orb_a
        integer(n_int), intent(in) :: ilutI(0:niftot) 
        integer, intent(out) :: orb_list(nbasis/2, 2) 
        real(dp), intent(out) :: cum_arr(nbasis/2), cum_sum
        integer, intent(in), optional :: tgt 
        real(dp), intent(out), optional :: cpt 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "create_bc_list_hubbard" 
#endif
        integer :: b, c, ex(2,3), spin, orb_b
        real(dp) :: elem

        orb_list = -1 
        cum_arr = 0.0_dp 
        cum_sum = 0.0_dp 
        ex(1,:) = src
        ex(2,1) = orb_a 

        ! we want to do a restriction! to make it easier to recalculate the 
        ! pgens and stuff! 
        ! the restriction is, that the first picked orbital (a) is of the 
        ! minority spin of the picked electrons. 
        ! so if we have to alpha and 1 beta electron picked, orbital (a) is 
        ! a beta orbital and the last 2 orbitals will be alpha and v.v.

        ! decide that the  first orbital orb_a, which is already picked is the 
        ! minority spin electron, so now pick two electrons of the opposite 
        ! spin 
        if (is_beta(orb_a)) then 
            ! then we want alpha orbitals
            spin = 0 
            ! and also be sure that we did the right thing until now, 
            ! otherwise it breaks 
            ASSERT(sum(get_spin_pn(src)) == 1)
        else 
            ! otherwise we want beta now
            spin = 1 
            ASSERT(sum(get_spin_pn(src)) == -1)
        end if

        if (present(tgt)) then 
            ASSERT(present(cpt))

            cpt = 0.0_dp 

            ! does the spin of tgt fit? 
            ! it must be opposite of orb_a!
            if (same_spin(orb_a, tgt)) return

            !TODO: we only need to consider one spin-type!!
            do b = 1, nbasis/2
                elem = 0.0_dp

                ! convert to the appropriate spin-orbitals 
                orb_b = 2 * b - spin

                if (IsNotOcc(ilutI,orb_b)) then 
                    ! get the appropriate momentum conserverd orbital c
                    c = get_orb_from_kpoints_three(src, orb_a, orb_b) 
                    
                    if (c /= orb_b .and. IsNotOcc(ilutI,c)) then 

                        ex(2,2:3) = [orb_b, c]

                        elem = abs(get_3_body_helement_ks_hub(nI, ex, .false.))

                    end if
                end if
                cum_sum = cum_sum + elem 
                if (tgt == orb_b) then 
                    cpt = elem 
                end if
            end do

            if (cum_sum < EPS) then 
                cpt = 0.0_dp 
            else 
                cpt = cpt / elem 
            end if
        else 
            do b = 1, nbasis/2 
                orb_b = 2 * b - spin 

                elem = 0.0_dp 

                if (IsNotOcc(ilutI, orb_b)) then 
                    c = get_orb_from_kpoints_three(src, orb_a, orb_b) 

                    if (c /= orb_b .and. IsNotOcc(ilutI,c)) then 

                        ex(2,2:3) = [orb_b, c] 

                        elem = abs(get_3_body_helement_ks_hub(nI, ex, .false.))

                    end if
                end if 
                cum_sum = cum_sum + elem 
                cum_arr(b) = cum_sum 
                orb_list(b,:) = [orb_b,c] 
            end do
        end if
                
    end subroutine create_bc_list_hubbard

    subroutine pick_a_orbital_hubbard(ilutI, orb, p_orb, sum_ms)
        ! hm... i think the first orbital can be picked totally random out 
        ! of the empty orbitals or? since every spin and every momentum 
        ! is allowed.. and i do not want to overdo the weighting i guess, 
        ! especially not for the beginning 
        integer(n_int), intent(in) :: ilutI(0:niftot) 
        integer, intent(out) :: orb 
        real(dp), intent(out) :: p_orb
        integer, intent(in), optional :: sum_ms
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_a_orbital_hubbard"
#endif
        integer :: spin

        ! if sum_ms is present, we pick the first orbital from the minority 
        ! spins in the picked electrons -> so it is the opposite 
        if (present(sum_ms)) then 
            if (sum_ms == -1) then 
                ! there a are 2 beta and one alpha electron picked -> 
                ! so pick alpha here! 
                spin = 0 
                p_orb = 1.0_dp / real(nbasis/2 - nOccAlpha, dp) 
            else if (sum_ms == 1) then 
                spin = -1
                p_orb = 1.0_dp / real(nbasis/2 - nOccBeta, dp) 
            end if

            do 
                orb = 2*(1 + int(genrand_real2_dsfmt() * nbasis/2)) + spin

                if (IsNotOcc(ilutI, orb)) exit 

            end do
        else 

            do 
                orb = 1 + int(genrand_real2_dsfmt() * nbasis) 

                if (IsNotOcc(ilutI, orb)) exit 

            end do

            p_orb = 1.0_dp/real(nbasis - nel, dp) 
        end if

    end subroutine pick_a_orbital_hubbard

    function get_orb_from_kpoints_three(src, orb_a, orb_b) result(orb_c) 
        ! momentum conservation for three-body terms 
        integer, intent(in) :: src(3), orb_a, orb_b
        integer :: orb_c 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_orb_from_kpoints_three"
#endif
        integer :: sum_ms, kc(3), spin_c, spin_ab

        ! implement that generally for also all-spin parallel excitation, which 
        ! might be necessary in the future.. 
        sum_ms = sum(get_spin_pn(src))

        ASSERT(sum_ms == -3 .or. sum_ms == -1 .or. sum_ms == 1 .or. sum_ms == 3)

        ! momentum conservation: ka + kb + kc = ki + kj + kk
        kc = G1(src(1))%k + G1(src(2))%k + G1(src(3))%k - G1(orb_a)%k - G1(orb_b)%k 

        ! perdiodic BC: 
        if (tHub .or. t_k_space_hubbard) then 
            call mompbcsym(kc, nBasisMax) 
        end if

        ! and now we want the spin: 
        if (sum_ms == -3) then 
            ! assert we can still reach the desired spin:
            ASSERT(get_spin_pn(orb_a) + get_spin_pn(orb_b) == -2)

            spin_c = 1 
        else if (sum_ms == 3) then 
            ASSERT(get_spin_pn(orb_a) + get_spin_pn(orb_b) == 2) 

            spin_c = 2

        else 
            spin_ab = get_ispn([orb_a,orb_b])

            ! we know we are not totally parallel so if: 
            if (spin_ab == 1) then 
                ! if we have already two beta, we need alpha
                spin_c = 2
            else if (spin_ab == 3) then 
                ! if we have two alpha we need a beta 
                spin_c = 1 

            else 
                ! in this case we need a spin to reach the final ms
                if (sum_ms == -1) then 
                    spin_c = 1
                else 
                    spin_c = 2
                end if
            end if
        end if

        orb_c = KPointToBasisFn(kc(1),kc(2),kc(3),spin_c) 

    end function get_orb_from_kpoints_three

    subroutine pick_three_opp_elecs(nI, elecs, p_elec, opt_sum_ms) 
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: elecs(3)
        real(dp), intent(out) :: p_elec 
        integer, intent(out), optional :: opt_sum_ms
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_three_opp_elecs"
#endif
        integer :: sum_ms

        ! this can be done way more effective than this simple implementation:
        first: do 
            elecs(1) = 1 + int(genrand_real2_dsfmt() * nel) 

            do 
                elecs(2) = 1 + int(genrand_real2_dsfmt() * nel) 

                if (elecs(1) /= elecs(2)) exit 
            end do

            if (get_ispn(nI(elecs(1:2))) /= 2) then 
                ! if the already picked electrons have same spin
                ! take opposite spin
                do 
                    elecs(3) = 1 + int(genrand_real2_dsfmt() * nel) 

                    if (elecs(1) /= elecs(3) .and. elecs(2) /= elecs(3) .and. &
                        .not. same_spin(nI(elecs(1)),nI(elecs(3)))) then 
                        exit first
                    end if
                end do
            else 
                ! here we can pick any spin electron 

                do 
                    elecs(3) = 1 + int(genrand_real2_dsfmt() * nel) 

                    if (elecs(1) /= elecs(3) .and. elecs(2) /= elecs(3)) then 
                        exit first
                    end if 
                end do
            end if
        end do first
                        
        ! i have to be careful how i could have gotten the electrons.. 
        ! because the order does not matter and there is no restriction set 

        ! the first two electrons are picked totally randon, independent of 
        ! spin just with the restriction i /= j 
        ! so p(i) = 1/nel, p(j|i) = 1/(nel-1) 
        ! the third electron has a spin restriction.. 
        ! if spin(i) = spin(j) then spin(k) must be -spin(i) giving 
        ! p(k|ij) = 1/n_opp_spin 
        ! if spin(i) /= spin(j) then k is freely except /= i,j 
        ! and since the order is not important the total probability is 
        ! p(i) * p(i|j) * [p(k|ij) + p(i|jk) + p(j|ik)] 
        ! which translates to 
        ! 1/nel * 1/(nel-1) * [1/n_opp + 2/(nel-2)] 
        ! so depending on the total mz, we get 
        sum_ms = sum(get_spin_pn(nI(elecs)))
        ASSERT(sum_ms == 1 .or. sum_ms == -1)
        if (sum_ms == 1) then 
            ! then we have 2 alpha and 2 beta electrons 
            p_elec = 1.0_dp/real(nel*(nel-1),dp) * &
                (1.0_dp/real(nOccBeta,dp) + 2.0_dp/real(nel-2,dp))
        else if (sum_ms == -1) then
            ! then we have 2 beta electrons and 1 alpha 
            p_elec = 1.0_dp/real(nel*(nel-1),dp) * & 
                (1.0_dp/real(nOccAlpha,dp) + 2.0_dp/real(nel-2,dp))
        end if

        if (present(opt_sum_ms)) opt_sum_ms = sum_ms

    end subroutine pick_three_opp_elecs

    subroutine gen_parallel_double_hubbard(nI, ilutI, nJ, ilutJ, ex, tParity, pgen)
        integer, intent(in) :: nI(nel) 
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ex(2,2) 
        integer, intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
#ifdef __DEBUG
        character(*), parameter :: this_routine = "gen_parallel_double_hubbard"
#endif
        real(dp) :: p_elec, p_orb
        integer :: elecs(2), orbs(2), src(2), ispn

        ! in the transcorrelated case we have to decide 
        ! i first have to choose an electron pair (ij) at random 
        ! but with the condition that they have to have opposite spin! 
        ! this is the only difference: i pick two spin-parallel electrons.. 
        ! the rest stays the same.. i just have to adjust the 
        ! get_orb_from_kpoints routine 
        ! and the matrix element calculation
        call pick_spin_par_elecs(nI, elecs, p_elec) 

        src = nI(elecs)

        ! i realised i could reuse the already implemented orbital picker, 
        ! but the other one does not use the fact that we know that we 
        ! have parallel spin in this case! so implement a parallel 
        ! spin-orbital picker!! (also better for matrix elements.. so we 
        ! must not check if the spins fit, if we only take the correct ones!)
        call pick_ab_orbitals_par_hubbard(nI, ilutI, src, orbs, p_orb)

        if (orbs(1) == ABORT_EXCITATION) then 
            nJ(1) = ABORT_EXCITATION
            pgen = 0.0_dp
            return 
        end if

        ! and make the excitation 
        call make_double(nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), ex, tParity)

        ilutJ = make_ilutJ(ilutI, ex, 2) 

        ! i think in both the electrons and the orbitals i have twice the 
        ! probability to pick them
        pgen = p_elec * p_orb * 2.0_dp

    end subroutine gen_parallel_double_hubbard

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

        ! actually the probability is twice that or? 
        ! or doesnt that matter, since it is the same
        p_elec = 1.0_dp / real(nOccBeta * nOccAlpha, dp)

    end subroutine pick_spin_opp_elecs

    subroutine pick_spin_par_elecs(nI, elecs, p_elec, opt_ispn)
        integer, intent(in) :: nI(nel) 
        integer, intent(out) :: elecs(2)
        real(dp), intent(out) :: p_elec 
        integer, intent(out), optional :: opt_ispn
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "pick_spin_par_elecs"
#endif 
        integer :: ispn

        do 
            elecs(1) = 1 + int(genrand_real2_dsfmt() * nel)

            do 
                elecs(2) = 1 + int(genrand_real2_dsfmt() * nel)
                
                if (elecs(1) /= elecs(2)) exit 

            end do
            ispn = get_ispn(nI(elecs))
            if (ispn /= 2) exit 

        end do

        ASSERT(ispn == 1 .or. ispn == 3) 

        if (ispn == 1) then 
            p_elec = 2.0_dp / real(nOccBeta*(nOccBeta-1),dp)
        else if (ispn == 3) then 
            p_elec = 2.0_dp / real(nOccAlpha*(nOccAlpha-1),dp)
        end if

        if (present(opt_ispn)) opt_ispn = ispn

    end subroutine pick_spin_par_elecs

    subroutine pick_ab_orbitals_par_hubbard(nI, ilutI, src, orbs, p_orb) 
        integer, intent(in) :: nI(nel), src(2)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: orbs(2) 
        real(dp), intent(out) :: p_orb
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_ab_orbitals_par_hubbard"
        real(dp) :: test
        integer :: ex(2,2)
#endif
        real(dp) :: cum_arr(nbasis/2)
        real(dp) :: cum_sum
        integer :: orb_list(nbasis/2, 2)
        integer :: ind

        ! without transcorrelation factor this is uniform, but with a 
        ! transcorrelation factor the matrix element might change and so also 
        ! the pgen should change. 
        call create_ab_list_par_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum)

        if (cum_sum < EPS) then 
            orbs(1) = ABORT_EXCITATION
            return
        end if

        ! this stuff is also written so often i should finally make a routine 
        ! out of that 
        call pick_from_cum_list(cum_arr, cum_sum, ind, p_orb)

        orbs = orb_list(ind,:)

#ifdef __DEBUG 
        ! check that the other way of picking the orbital has the same 
        ! probability.. 
        call create_ab_list_par_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum, & 
            orbs(2), test)

        if (abs(test - p_orb) > 1.e-8) then 
            print *, "pgen assumption wrong:!" 
            print *, "p_orb: ", p_orb
            print *, "test: ", test 
            print *, "orbs: ", orbs
        end if

        !todo: also call the calc_pgen_k_space_hubbard here and check 
        ! pgens 
        ex(1,:) = src
        ex(2,:) = orbs 

#endif

        ! do i have to recalc. the pgen the other way around? yes! 
        ! effectively reuse the above functionality
        ! i am pretty sure i just have to find the position in the 
        ! list.. OR: since in the hubbard it is just twice the 
        ! probability or? i am pretty sure yes.. but for all of them.. 
        ! so in the end it shouldnt matter again..
!         p_orb = 2.0_dp * p_orb

    end subroutine pick_ab_orbitals_par_hubbard

    subroutine pick_ab_orbitals_hubbard(nI, ilutI, src, orbs, p_orb) 
        ! depending on the already picked electrons (ij) pick an orbital 
        ! (a) and the connected orbital (b)
        integer, intent(in) :: nI(nel), src(2)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: orbs(2) 
        real(dp), intent(out) :: p_orb
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_ab_orbitals_hubbard"
        real(dp) :: test
        integer :: ex(2,2)
#endif
        real(dp) :: cum_arr(nbasis)
        real(dp) :: cum_sum
        integer :: orb_list(nbasis, 2)
        integer :: ind

        ! without transcorrelation factor this is uniform, but with a 
        ! transcorrelation factor the matrix element might change and so also 
        ! the pgen should change. 
        call create_ab_list_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum)

        if (cum_sum < EPS) then 
            orbs(1) = ABORT_EXCITATION
            return
        end if

        ! this stuff is also written so often i should finally make a routine 
        ! out of that 
        call pick_from_cum_list(cum_arr, cum_sum, ind, p_orb)

        orbs = orb_list(ind,:)

#ifdef __DEBUG 
        ! check that the other way of picking the orbital has the same 
        ! probability.. 
        call create_ab_list_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum, & 
            orbs(2), test)

        if (abs(test - p_orb) > 1.e-8) then 
            print *, "pgen assumption wrong:!" 
            print *, "p_orb: ", p_orb
            print *, "test: ", test 
            print *, "orbs: ", orbs
        end if

        !todo: also call the calc_pgen_k_space_hubbard here and check 
        ! pgens 
        ex(1,:) = src
        ex(2,:) = orbs 

#endif

        ! do i have to recalc. the pgen the other way around? yes! 
        ! effectively reuse the above functionality
        ! i am pretty sure i just have to find the position in the 
        ! list.. OR: since in the hubbard it is just twice the 
        ! probability or? i am pretty sure yes.. but for all of them.. 
        ! so in the end it shouldnt matter again..
!         p_orb = 2.0_dp * p_orb

    end subroutine pick_ab_orbitals_hubbard

    subroutine create_ab_list_par_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum, & 
            tgt, cpt) 
        integer, intent(in) :: nI(nel), src(2) 
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: orb_list(nbasis/2, 2)
        real(dp), intent(out) :: cum_arr(nbasis/2)
        real(dp), intent(out) :: cum_sum 
        integer, intent(in), optional :: tgt 
        real(dp), intent(out), optional :: cpt 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "create_ab_list_par_hubbard"
#endif
        integer :: a, b, ex(2,2), spin, orb_a
        real(dp) :: elem
        ! do the cum_arr for the k-space hubbard 
        ! i think here i might really use flags.. and not just do the 
        ! influence over the matrix elements.. since without transcorrelation 
        ! i would waste alot of effort if i calculate the matrix elements 
        ! here all the time.. 
        orb_list = -1 
        cum_arr = 0.0_dp 
        cum_sum = 0.0_dp 

        ex(1,:) = src
        
        ! this routine only checks for parallel spins, depending on src
        ASSERT(same_spin(src(1),src(2)))

        ! make a spin factor for the orbital conversion
        ! 0...alpha
        ! 1...beta
        spin = get_spin(src(1)) - 1

        ! and only loop over the correct spin
        if (present(tgt)) then 
            ASSERT(present(cpt))

            cpt = 0.0_dp 

            ! if target does not have the same spin, do we abort or return?
            if (.not. same_spin(src(1),tgt)) return

            do a = 1, nbasis/2
                elem = 0.0_dp 

                orb_a = 2 * a - spin

                if (IsNotOcc(ilutI,orb_a)) then 
                    ! modify get_orb_from_kpoints to give spin-parallel
                    ! it already does i think! 
                    b = get_orb_from_kpoints(src(1),src(2),orb_a)

                    if (b /= orb_a .and. IsNotOcc(ilutI,b)) then 

                        ex(2,:) = [orb_a,b] 

                        ! modify the matrix element calculation or 
                        ! write a new routine for the transcorrelated.. 
                        elem = abs(get_offdiag_helement_k_sp_hub(nI, ex, .false.))

                    end if
                end if
                cum_sum = cum_sum + elem
                if (tgt == orb_a) then 
                    cpt = elem 
                end if
            end do
            if (cum_sum < EPS) then 
                cpt = 0.0_dp
            else 
                cpt = cpt / cum_sum 
            end if
        else 
            do a = 1, nbasis/2
                orb_a = 2 * a - spin 

                elem = 0.0_dp 

                if (IsNotOcc(ilutI,orb_a)) then 
                    b = get_orb_from_kpoints(src(1),src(2), orb_a)

                    if (b /= orb_a .and. IsNotOcc(ilutI, orb_a)) then 
                        ex(2,:) = [orb_a, b]

                        elem = abs(get_offdiag_helement_k_sp_hub(nI, ex, .false.))

                    end if
                end if
                cum_sum = cum_sum + elem 
                cum_arr(a) = cum_sum 
                orb_list(a,:) = [orb_a,b]
            end do
        end if 
        
    end subroutine create_ab_list_par_hubbard
      
    subroutine create_ab_list_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum, & 
            tgt, cpt) 
        integer, intent(in) :: nI(nel), src(2) 
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: orb_list(nbasis, 2)
        real(dp), intent(out) :: cum_arr(nbasis)
        real(dp), intent(out) :: cum_sum 
        integer, intent(in), optional :: tgt 
        real(dp), intent(out), optional :: cpt 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "create_ab_list_hubbard"
#endif
        integer :: a, b, ex(2,2)
        real(dp) :: elem
        ! do the cum_arr for the k-space hubbard 
        ! i think here i might really use flags.. and not just do the 
        ! influence over the matrix elements.. since without transcorrelation 
        ! i would waste alot of effort if i calculate the matrix elements 
        ! here all the time.. 
        orb_list = -1 
        cum_arr = 0.0_dp 
        cum_sum = 0.0_dp 

        ex(1,:) = src

        ! we are also using this routine for the parallel excitations due to 
        ! the transcorrelation factor.. this is nice, since it is easily 
        ! reusable, but, loses alot of efficiency, since the we are looping 
        ! over all spin-orbital, although we know we only want to loop over a 
        ! certain spin! 
        ! todo
        if (present(tgt)) then 
            ASSERT(present(cpt))

            cpt = 0.0_dp

            do a = 1, nbasis
                elem = 0.0_dp
                ! if a is empty
                if (IsNotOcc(ilutI, a)) then 
                    ! i have to rewrite get_orb, so it gives me the same 
                    ! spin if src has the same spin! todo
                    ! to take into account spin-parallel double 
                    ! excitations!
                    b = get_orb_from_kpoints(src(1), src(2), a)

                    ! and b is empty and not a
                    if (b /= a .and. IsNotOcc(ilutI, b)) then
                        ! is it sure that we have opposite spin?
                        ! i have to do better asserts!
                        if (.not. t_trans_corr_2body) then
                            ASSERT(.not. same_spin(a,b))
                        end if

                        ex(2,:) = [a,b]
                        ! in the matrix element routine the check for 
                        ! transcorrelation is done.. although i could do it 
                        ! more efficiently out here.. todo
                        elem = abs(get_offdiag_helement_k_sp_hub(nI, ex, .false.))
                    end if
                end if
                cum_sum = cum_sum + elem 

                if (tgt == a)  then 
                    cpt = elem 
                end if
            end do
            if (cum_sum < EPS) then 
                cpt = 0.0_dp 
            else 
                cpt = cpt / cum_sum 
            end if
        else 
            do a = 1, nbasis
                elem = 0.0_dp 

                if (IsNotOcc(ilutI, a)) then 
                    b = get_orb_from_kpoints(src(1), src(2), a)

                    if (b /= a .and. IsNotOcc(ilutI, b)) then 

                        ex(2,:) = [a,b]
                        elem = abs(get_offdiag_helement_k_sp_hub(nI, ex, .false.))

                    end if
                end if
                cum_sum = cum_sum + elem 
                cum_arr(a) = cum_sum 
                orb_list(a,:) = [a,b] 
            end do
        end if

    end subroutine create_ab_list_hubbard

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

    function calc_pgen_k_space_hubbard_transcorr(nI, ilutI, ex, ic) result(pgen)
        ! this function i have to rewrite for the transcorrelated to take 
        ! the same-spin doubles and triples into account! 
        ! NOTE: ex could be of form (2,3) in the case of triples!
        integer, intent(in) :: nI(nel), ex(:,:), ic 
        integer(n_int), intent(in) :: ilutI(0:niftot) 
        real(dp) :: pgen 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "calc_pgen_k_space_hubbard_transcorr" 
#endif 
    
        pgen = 0.0_dp

        if (ic == 2) then 
            if (same_spin(ex(1,1),ex(1,2))) then 
                ! parallel double excitation
                ! the spins are checked within the function:
                pgen = calc_pgen_k_space_hubbard_par(nI,ilutI,ex,ic)

                pgen = pgen * pDoubles * pParallel

            else 
                ! "normal" opposite spin excitation
                ! the spins are checked within the function: 
                pgen = calc_pgen_k_space_hubbard(nI, ilutI, ex, ic) 

                pgen = pgen * pDoubles * (1.0_dp - pParallel) 
            end if
        else if (ic == 3) then 
            pgen = calc_pgen_k_space_hubbard_triples(nI, ilutI, ex, ic)

            pgen = pgen * (1.0_dp - pDoubles) 

        end if

    end function calc_pgen_k_space_hubbard_transcorr

    function calc_pgen_k_space_hubbard_triples(nI, ilutI, ex, ic) result(pgen)
        integer, intent(in) :: nI(nel), ex(:,:), ic
        integer(n_int), intent(in) :: ilutI(0:niftot) 
        real(dp) :: pgen 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "calc_pgen_k_space_hubbard_triples"
        real(dp) :: test
#endif
        real(dp) :: p_elec, p_orb(2), cum_arr(nbasis/2), cum_sum
        integer :: orb_list(nbasis/2,2), sum_ms, orb_a, orbs(2)

        if (ic /= 3) then 
            pgen = 0.0_dp
            return
        end if

        sum_ms = sum(get_spin_pn(ex(1,:)))

        ! check spins
        if (.not. (sum_ms == 1 .or. sum_ms == -1) .or. sum_ms /= sum(get_spin_pn(ex(2,:)))) then
            pgen = 0.0_dp
            return
        end if

        ! get the probabilites for the electrons and orbital (a)
        if (sum_ms == 1) then 
            p_elec = 1.0_dp / real(nel*(nel-1),dp) * & 
                (1.0_dp/real(nOccBeta,dp) + 2.0_dp / real(nel-2,dp))

            p_orb(1) = 1.0_dp / real(nbasis/2 - nOccBeta, dp) 

        else 
            p_elec = 1.0_dp / real(nel*(nel-1),dp) * & 
                (1.0_dp/real(nOccAlpha,dp) + 2.0_dp / real(nel-2,dp))

            p_orb(1) = 1.0_dp / real(nbasis/2 - nOccAlpha, dp)

        end if

        ! for this i need the minority spin orbital (a)
        orb_a = find_minority_spin(ex(2,:))

        orbs = pack(ex(2,:), ex(2,:) /= orb_a) 

        call create_bc_list_hubbard(nI, ilutI, ex(1,:), orb_a, orb_list, cum_arr, & 
            cum_sum, orbs(1), p_orb(2))

        pgen = p_elec * product(p_orb) * 2.0_dp 

#ifdef __DEBUG 
        call create_bc_list_hubbard(nI, ilutI, ex(1,:), orb_a, orb_list, cum_arr, & 
            cum_sum, orbs(2), test)

        if (abs(test - p_orb(2)) > 1.e-8) then 
            print *, "pgen assumption wrong:!" 
            print *, "p_orb: ", p_orb(2)
            print *, "test: ", test 
            print *, "ex(2,:): ", ex(2,:)
        end if
#endif
        
    end function calc_pgen_k_space_hubbard_triples

    function calc_pgen_k_space_hubbard_par(nI, ilutI, ex, ic) result(pgen) 
        integer, intent(in) :: nI(nel), ex(:,:), ic
        integer(n_int), intent(in) :: ilutI(0:niftot) 
        real(dp) :: pgen 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "calc_pgen_k_space_hubbard_par"
        real(dp) :: test
#endif
        real(dp) :: p_elec, p_orb, cum_arr(nbasis/2), cum_sum
        integer :: orb_list(nbasis/2,2)

        ! check ic:
        if (ic /= 2) then 
            pgen = 0.0_dp
            return
        end if

        ! check spin:
        if (.not.(same_spin(ex(1,1),ex(1,2)) .and. same_spin(ex(2,1),ex(2,2)) .and. & 
            same_spin(ex(1,1),ex(2,1)))) then 
            pgen = 0.0_dp
            return
        end if

        if (get_ispn(ex(1,1:2)) == 1) then 
            p_elec = 1.0_dp / real(nbasis/2 - nOccBeta, dp)
        else 
            p_elec = 1.0_dp / real(nbasis/2 - nOccAlpha, dp)
        end if

        call create_ab_list_par_hubbard(nI, ilutI, ex(1,1:2), orb_list, cum_arr, & 
            cum_sum, ex(2,1), p_orb)

        pgen = p_elec * p_orb * 2.0_dp

#ifdef __DEBUG 
        call create_ab_list_par_hubbard(nI, ilutI, ex(1,1:2), orb_list, cum_arr, & 
            cum_sum, ex(2,1), test)

        if (abs(test - p_orb) > 1.e-8) then 
            print *, "pgen assumption wrong:!" 
            print *, "p_orb: ", p_orb
            print *, "test: ", test 
            print *, "ex(2,:): ", ex(2,:)
        end if

#endif

    end function calc_pgen_k_space_hubbard_par

    function calc_pgen_k_space_hubbard(nI, ilutI, ex, ic) result(pgen) 
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(in) :: nI(nel), ex(2,2), ic
        real(dp) :: pgen
#ifdef __DEBUG
        character(*), parameter :: this_routine = "calc_pgen_k_space_hubbard"
        real(dp) :: test
#endif
        real(dp) :: p_elec, p_orb, cum_arr(nbasis), cum_sum
        integer :: orb_list(nbasis,2), src(2)

        if (ic /= 2) then 
            pgen = 0.0_dp 
            return 
        end if

        if (same_spin(ex(1,1),ex(1,2)) .or. same_spin(ex(2,1),ex(2,2))) then 
            pgen = 0.0_dp 
            return
        end if

        p_elec = 1.0_dp / real(nOccBeta * nOccAlpha, dp) 

        src = get_src(ex)

        call create_ab_list_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum, & 
                ex(2,1), p_orb) 

#ifdef __DEBUG
        call create_ab_list_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum, & 
                ex(2,2), test) 

        if (abs(test - p_orb) > 1.e-8) then 
            print *, "pgen assumption wrong:!" 
            print *, "p_orb: ", p_orb
            print *, "test: ", test 
            print *, "ex(2,:): ", ex(2,:)
        end if

#endif

        ! i do not need to recalc, the p(b|ij) since it is the same.. 
        ! but i need a factor of 2 somewhere.. figure that out!
        pgen = p_elec * p_orb * 2.0_dp
 
    end function calc_pgen_k_space_hubbard 

    subroutine init_get_helement_k_space_hub
        get_helement_lattice_ex_mat => get_helement_k_space_hub_ex_mat
        get_helement_lattice_general => get_helement_k_space_hub_general
        ! maybe i have to initialize more here, especially if we are using the 
        ! HPHF keyword I guess.. 

    end subroutine init_get_helement_k_space_hub

    function get_helement_k_space_hub_ex_mat(nI, ic, ex, tpar) result(hel)
        integer, intent(in) :: nI(nel), ic, ex(2,2)
        logical, intent(in) :: tpar 
        HElement_t(dp) :: hel
#ifdef __DEBUG 
        character(*), parameter :: this_routine ="get_helement_k_space_hub_ex_mat"
#endif

        !todo: if 2-body-transcorrelation, we can have triple excitations now..
        ! fix that here.. (and also in a lot of other parts in the code..)

        if (ic == 0) then 
            ! the diagonal is just the sum of the occupied one-particle 
            ! basis states 
            hel = sltcnd_0(nI) 

        else if (ic == 2) then 
            hel = get_offdiag_helement_k_sp_hub(nI, ex, tpar) 

        else 
            hel = h_cast(0.0_dp) 

        end if

    end function get_helement_k_space_hub_ex_mat

    function get_helement_k_space_hub_general(nI, nJ, ic_ret) result(hel) 
        integer, intent(in) :: nI(nel), nJ(nel) 
        integer, intent(inout), optional :: ic_ret
        HElement_t(dp) :: hel 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_helement_k_space_hub_general"
#endif
        integer :: ic, ex(2,2) 
        logical :: tpar 
        integer(n_int) :: ilutI(0:NIfTot), ilutJ(0:niftot)

        !todo: if 2-body-transcorrelation, we can have triple excitations now..
        ! fix that here.. (and also in a lot of other parts in the code..)
        if (present(ic_ret)) then 
            if (ic_ret == 0) then 
                hel = sltcnd_0(nI) 

            else if (ic_ret == 2) then 
                ex(1,1) = 2
                call GetExcitation(nI, nJ, nel, ex, tpar) 
                hel = get_offdiag_helement_k_sp_hub(nI, ex, tpar) 

            else if (ic_ret == -1) then 
                call EncodeBitDet(nI, ilutI) 
                call EncodeBitDet(nJ, ilutJ) 

                ic_ret = FindBitExcitLevel(ilutI, ilutJ) 

                if (ic_ret == 0) then 
                    hel = sltcnd_0(nI) 

                else if (ic_ret == 2) then 
                    ex(1,1) = 2 
                    call GetBitExcitation(ilutI, ilutJ, ex, tpar) 

                    hel = get_offdiag_helement_k_sp_hub(nI, ex, tpar) 

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
                hel = sltcnd_0(nI) 
            else if (ic == 2) then
                ex(1,1) = 2 
                call GetBitExcitation(ilutI, ilutJ, ex, tpar) 

                hel = get_offdiag_helement_k_sp_hub(nI, ex, tpar) 

            else 
                hel = h_cast(0.0_dp) 
            end if 
        end if

    end function get_helement_k_space_hub_general

    function get_diag_helement_k_sp_hub(nI) result(hel) 
        integer, intent(in) :: nI(nel) 
        HElement_t(dp) :: hel 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_diag_helement_k_sp_hub" 
#endif
        integer :: i, j, id(nel), idX, idN, spin 
        HElement_t(dp) :: hel_sing, hel_doub, hel_par, hel_opp

        ! todo: in the case of 2-body-transcorrelation, there are more 
        ! contributions.. 
        hel = h_cast(0.0_dp)

        if (t_trans_corr_2body) then 
            hel_sing = sum(GetTMatEl(nI,nI))

            id = gtID(nI) 

            hel_doub = h_cast(0.0_dp) 
            hel_par = h_cast(0.0_dp)
            hel_opp = h_cast(0.0_dp) 

            ! i do not need to run over the electrons, since all of this can 
            ! be calculated directly
            do i = 1, nel-1
                do j = i + 1, nel 

                    idX = max(id(i), id(j))
                    idN = min(id(i), id(j))

                    ! normal direct 
                    ! us the spin_restriction here directly! 
                    if (.not. same_spin(nI(i),nI(j))) then
                        hel_doub = hel_doub + get_hub_umat_el(idN,idX,idN,idX)
                    end if
                    
                    ! THEN WE do not need the exchange to cancel the 
                    ! incorrectly counted double excitations!
                    ! there is a bug in the hubbard i guess: the parallel 
                    ! spin excitations are taken into account.. for the 
                    ! diagonal contribution.. todo! 
!                     if (same_spin(nI(i),nI(j))) then 
!                         hel_doub = hel_doub - get_hub_umat_el(idN, idX, idX, idN)
!                     end if
                    ! and exchange terms 
                    ! actually for the "normal" double excitation, there is 
                    ! no exchange! 

                    ! we have the contribution from the parallel doubles now: 
                    ! this is really slow for now, i think most of that 
                    ! can be moved outside of the loop! 
                   if (is_beta(nI(i))) then 
                        spin = 1
                    else 
                        spin = -1
                    end if 

                    if (same_spin(nI(i),nI(j))) then 
                        hel_par = hel_par + 2.0_dp * three_body_prefac *  & 
                            get_one_body_diag(nI,spin) * & 
                            (1.0_dp - epsilon_kvec(G1(nI(i))%k - G1(nI(j))%k))

                    else 
                        ! take into account the opposite spin 3-body term
                        ! here i just have to be sure what is p and q in the 
                        ! formulasr.. 
                        hel_opp = hel_opp + three_body_transcorr_fac(nI, & 
                            G1(nI(i))%k, G1(nI(j))%k, [0,0,0], -spin)

                    end if 
                end do
            end do

            hel = hel_sing + hel_doub + hel_par + hel_opp

        else 
            hel = sltcnd_0(nI)
        end if

    end function get_diag_helement_k_sp_hub

    function get_one_body_diag(nI, spin) result(hel)
        integer, intent(in) :: nI(nel)
        integer, intent(in), optional :: spin 
        HElement_t(dp) :: hel 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_one_body_diag"
#endif 
        integer :: i

        ! the spin input: -1 is beta, +1 is alpha, 0 is both!
        ! if spin is not present, default is both!
        hel = h_cast(0.0_dp)

        if (present(spin)) then 
            ! either -1 or 1 input, if spin is given!
            ASSERT(spin == -1 .or. spin == 1)
            if (spin == -1) then 
                do i = 1, nel
                    if (is_beta(nI(i))) then 
                        ! is TMAT(nI) actually cos(k) ? or do i have to 
                        ! take this into account more specifically?
                        hel = hel + GetTMatEl(nI(i),nI(i))
                    end if
                end do
            else if (spin == 1) then
                do i = 1, nel 
                    if (is_alpha(nI(i))) then 
                        hel = hel + GetTMatEl(nI(i),nI(i))
                    end if 
                end do
            end if
        else 
            do i = 1, nel 
                hel = hel + GetTMatEl(nI(i),nI(i))
            end do
        end if

        ! remove the -t.. and add it afterwards is necesarry.. 
        ! this is just a quick hack. do it properly elsewhere!
        hel = hel / bhub

    end function get_one_body_diag

    function get_offdiag_helement_k_sp_hub(nI, ex, tpar) result(hel) 
        ! this routine is called for the double excitations in the 
        ! k-space hubbard. in case of transcorrelation, this can also be 
        ! spin-parallel excitations now. the triple excitation have a 
        ! seperate routine!
        integer, intent(in) :: nI(nel), ex(2,2)
        logical, intent(in) :: tpar 
        HElement_t(dp) :: hel 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_offdiag_helement_k_sp_hub"
#endif
        integer :: src(2), tgt(2), ij(2), ab(2), k_vec(3), spin

        src = get_src(ex)
        tgt = get_tgt(ex)
        if (.not. t_trans_corr_2body) then 
            if (same_spin(src(1),src(2)) .or. same_spin(tgt(1),tgt(2))) then 
                hel = h_cast(0.0_dp)
                return 
            end if
        else 
            ! if src has same spin but tgt has opposite spin -> 0 mat ele
            if (same_spin(src(1),src(2)) .and. (.not. same_spin(tgt(1),tgt(2)) &
                .or. .not. same_spin(src(1),tgt(1)))) then
                
                hel = h_cast(0.0_dp)
                return
            end if
        end if

        ij = gtid(src)
        ab = gtid(tgt) 
        ! that about the spin?? must spin(a) be equal spin(i) and same for 
        ! b and j? does this have an effect on the sign of the matrix element? 

        ! in the case of 2-body transcorrelation, parallel spin double exciattions 
        ! are possible todo: check if we get the coulomb and exchange contributions
        ! correct..

        ! the U part is still just the the spin-opposite part
        hel = get_hub_umat_el(ab(1),ab(2),ij(1),ij(2))

        ! if hel == 0, due to momentum conservation violation we can already 
        ! exit here, since this means this excitation is just no possible! 
        if (abs(hel) < EPS) return

        if (t_trans_corr) then 
            ! do something 
            hel = hel * exp(trans_corr_param/2.0_dp * & 
                (GetTMatEl(ij(1),ij(1)) + GetTMatEl(ij(2),ij(2))  & 
                - GetTMatEl(ab(1),ab(1)) - GetTMatEl(ab(2),ab(2)))) 
        end if

        if (t_trans_corr_2body) then 
            ! i need the k-vector of the transferred momentum.. 
            ! i am not sure if the orbitals involved in ex() are every 
            ! re-shuffled.. if yes, it is not so easy in the spin-parallel 
            ! case to reobtain the transferred momentum. although it must be 
            ! possible. for now just assume (ex(2,2)) is the final orbital b 
            ! with momentum k_i + k_j - k_a and we need the 
            ! k_j - k_a momentum
            
            
!             k_vec = get_transferred_momentum(ex)
            spin = get_spin_pn(src(1))
            if (same_spin(src(1),src(2))) then
                ! we need the spin of the excitation here if it is parallel

                ! in the same-spin case, this is the only contribution to the 
                ! matrix element
                ! and maybe i have to take the sign additionally into 
                ! account here?? or is this taken care of with tpar??

                ! thanks to Manu i have figured it out. we have to take 
                ! the momentum between the to equally possible excitations: 
                ! c^+_b c^+_a c_q c_p with W(q-a) 
                ! and 
                !-c^+_b c^+_a c_p c_q with W(p-a) 
                ! with one of the orbital spins. 
                ! i think it doesnt matter, which one. 
                ! although for the sign it maybe does.. check thate
                hel = same_spin_transcorr_factor(nI, G1(ex(1,1))%k - G1(ex(2,1))%k, spin) &
                    - same_spin_transcorr_factor(nI, G1(ex(1,2))%k - G1(ex(2,1))%k, spin)

            else 
                ! else we need the opposite spin contribution
                
                ! the two-body contribution needs two k-vector inputs. 
                ! figure out what momentum is necessary there! 
                ! i need the transferred momentum 
                ! and the momentum of other involved electron 
                ! which by definition of k, is always ex(1,2) todo: 
                ! check if this works as intented
                hel = hel + two_body_transcorr_factor(G1(ex(1,2))%k, k_vec)
                
                ! and now the 3-body contribution: 
                ! which also needs the third involved mometum, which then 
                ! again is ex(1,1)
                ! todo.. figure out spins! 
                hel = hel + three_body_transcorr_fac(nI, G1(ex(1,2))%k, & 
                    G1(ex(1,1))%k, k_vec, spin)

            end if
        end if

        if (tpar) hel = -hel 

    end function get_offdiag_helement_k_sp_hub

    function get_transferred_momentum(ex) result(k_vec) 
        ! routine to reobtain transferred momentum from a given excitation 
        ! for spin-opposite double excitations i am pretty sure how, but 
        ! for triple excitations and spin-parallel doubles not so much.. todo
        integer, intent(in) :: ex(:,:) 
        integer :: k_vec(3) 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_transferred_momentum"
#endif
        
        ASSERT(size(ex,1) == 2)
        ASSERT(size(ex,2) == 2 .or. size(ex,2) == 3)

        if (size(ex,2) == 2) then 
            ! double excitation
            if (same_spin(ex(1,1),ex(1,2))) then 
                ! spin-parallel excitation
                ASSERT(same_spin(ex(2,1),ex(2,2)))
                ASSERT(same_spin(ex(1,1),ex(2,1)))

                ! for now just take the momentum of ex(1,2) - ex(2,1) 
                call mompbcsym(G1(ex(1,2))%k - G1(ex(2,1))%k, nBasisMax)

            else 
                ! "normal" hubbard spin-opposite excitation
                ASSERT(.not. same_spin(ex(2,1),ex(2,2)))
                ! here it is easier, we need the momentum difference of the 
                ! same spin-electrons 
                ! the sign of k should be irrelevant or? todo!
                if (same_spin(ex(1,1),ex(2,1))) then 
                    call mompbcsym(G1(ex(1,1))%k - G1(ex(2,1))%k, nBasisMax)
                else 
                    call mompbcsym(G1(ex(1,1))%k - G1(ex(2,2))%k,nBasisMax)
                end if
            end if
        else 
            ! triple excitations..
            !todo
        end if

    end function get_transferred_momentum

    ! finally write the functions to setup up the pesky G1 and nBasisMax 
    ! quantities to be consistent with the rest of the old code 
    subroutine setup_g1(in_lat) 
        use SystemData, only: G1
        class(lattice), intent(in), optional :: in_lat
        character(*), parameter :: this_routine = "setup_g1"

        type(BasisFN) :: temp_g
        integer :: i,j,k,l,ind
        logical :: kallowed

        ! i think everything is in the System_neci file
        if (present(in_lat)) then 
            ! only do it if G1 has not been setup yet!
            if (.not. associated(G1)) then
                ! i need number of spin-orbitals
                allocate(G1(in_lat%get_nsites()*2))
                G1 = NullBasisFn
                
                ! should i rely on the already setup nBasisMax?
                if (all(nBasisMax == 0)) then 
                    call setup_nbasismax(in_lat)
                end if
                ind = 0
                do i = nBasisMax(1,1), nBasisMax(1,2)
                    do j = nBasisMax(2,1), nBasisMax(2,2)
                        do k = nBasisMax(3,1), nBasisMax(3,2)
                            do l  = nBasisMax(4,1), nBasisMax(4,2), 2
                               
                                temp_g%k = [i,j,k]
                                temp_g%ms = l 
                                if ((treal .and. .not. ttilt) .or. kallowed(temp_g, nBasisMax)) then
                                    ind = ind + 1 
                                    G1(ind)%k = [i,j,k] 
                                    G1(ind)%ms = l
                                    G1(ind)%Sym = TotSymRep()
                                    if (.not. in_lat%is_k_space()) then 
                                        ! turn off- symmetry in the hubbard case
                                        G1(ind)%sym%s = 0
                                    end if
                                end if
                            end do
                        end do
                    end do
                end do
                if (in_lat%is_k_space()) then 
                    call GenHubMomIrrepsSymTable(G1, in_lat%get_nsites()*2, nbasismax)
                else 
                    ! also to the rest of the symmetry stuff here: 
                    ! in case of real-space turn off symmetry completely: 
                    call GenMolpSymTable(1, G1, in_lat%get_nsites()*2)
                    ! and i have to redo the symmetry setting to 0 
                    do i = 1, in_lat%get_nsites()*2
                        G1(i)%sym%s = 0
                    end do
                end if

            end if
        else 
            ! not yet implemented!
            call Stop_All(this_routine, "not yet implemented")
        end if

    end subroutine setup_g1

    subroutine setup_nbasismax(in_lat) 
        use SystemData, only: nBasisMax
        class(lattice), intent(in), optional :: in_lat
        character(*), parameter :: this_routine =  "setup_nbasismax"

        integer :: dummy_size
        ! thats a fucking pain in the ass.. i do not want to do that now!
        if (present(in_lat)) then 
            if (all(nBasisMax == 0)) then 
                ! only do smth if nbasismax was not changed yet

                ! whatever spin-polar means: 
                if (TSPINPOLAR) then 
                    nBasisMax(4,1) = 1 
                    nBasisMax(2,3) = 1 
                else 
                    nBasisMax(4,1) = -1
                    if (nBasisMax(2,3) == 0) nBasisMax(2,3) = 2 
                end if

                ! this is never explained: 
                nBasisMax(4,2) = 1

                ! i should give lattice also a member type and a k-space flag..
                if (trim(in_lat%get_name()) == 'tilted') then 
                    ! how do i get nmaxx and the rest effectively?? 
!                     call SETBASISLIM_HUBTILT(nBasisMax, nmaxx, nmaxy, nmaxz, & 
!                         in_lat%gen_nsites()*2, in_lat%is_periodic(), itiltx, itilty))
                    ! if it is tilted the nmax stuff is usualy 1 or?? 
                    call SETBASISLIM_HUBTILT(nBasisMax, 1,1,1, dummy_size, & 
                        in_lat%is_periodic(), in_lat%get_length(1), in_lat%get_length(2))
                else 
                    call SETBASISLIM_HUB(nBasisMax, in_lat%get_length(1), & 
                        in_lat%get_length(2), in_lat%get_length(3), dummy_size, & 
                        in_lat%is_periodic(), .not. in_lat%is_k_space())
                end if
                if (thub .and. treal) then
                    ! apparently this allows integrals between different 
                    ! spins: so in the transcorrelated hubbard this should 
                    ! be changed also maybe? 
                    nBasisMax(2,3) = 1
                end if
                ASSERT(dummy_size == in_lat%get_nsites()*2)
            end if
        else 
            call Stop_All(this_routine, "not yet implemented")
        end if

    end subroutine setup_nbasismax

    ! create the necessary routines for the triple excitation in the 
    ! 2-body transcorrelated k-space hamiltonian 
    HElement_t(dp) function same_spin_transcorr_factor(nI, k_vec, spin)
        ! this is the term coming appearing in the spin-parallel 
        ! excitations coming from the k = 0 triple excitation
        integer, intent(in) :: nI(nel), k_vec(N_DIM), spin

        same_spin_transcorr_factor = three_body_prefac * get_one_body_diag(nI,-spin) * &
                                     epsilon_kvec(k_vec)

    end function same_spin_transcorr_factor

    HElement_t(dp) function epsilon_kvec(k_vec)
        ! and actually this function has to be defined differently for 
        ! different type of lattices! TODO!
        ! actually i could get rid of this function and directly call 
        ! the dispersion relation of the lattice.. 
        integer, intent(in) :: k_vec(N_DIM)
#ifdef __DEBUG
        character(*), parameter :: this_routine = "epsilon_kvec"
#endif

        ASSERT(associated(lat)) 

        ! i could save the basic lattice vectors for the lattice or even 
        ! store the dispersion relation for each lattice type and call it 
        ! with a given k-vector? 
        epsilon_kvec = h_cast(lat%dispersion_rel(k_vec))
        

!         epsilon_kvec = h_cast(2*sum(cos(real(k_vec,dp))))

    end function epsilon_kvec

    HElement_t(dp) function two_body_transcorr_factor(p,k) 
        integer, intent(in) :: p(N_DIM), k(N_DIM) 
        
        ! take out the part with U/2 since this is already covered in the 
        ! "normal" matrix elements
        two_body_transcorr_factor = real(bhub,dp)/real(omega,dp)*(&
            (exp(trans_corr_param_2body) - 1.0_dp) * epsilon_kvec(p - k) + &
            (exp(-trans_corr_param_2body) - 1.0_dp) * epsilon_kvec(p))

        ! thats it i gues.. 
    end function two_body_transcorr_factor

    HElement_t(dp) function three_body_transcorr_fac(nI, p, q, k, spin) 
        integer, intent(in) :: nI(nel), p(N_DIM), q(N_DIM), k(N_DIM), spin 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "three_body_transcorr_fac"
#endif
        real(dp) :: n_opp

        ASSERT(spin == 1 .or. spin == -1)

        ! i have to deside what i want to input here.. as spin sigma or -sigma..
        ! and here i want the number of electrons with the opposite spin
        if (spin == -1) then 
            n_opp = real(nOccAlpha,dp)
        else if (spin == 1) then 
            n_opp = real(nOccBeta, dp) 
        end if

        three_body_transcorr_fac = three_body_prefac * (&
            n_opp * (epsilon_kvec(p) + epsilon_kvec(p - k)) - & 
            get_one_body_diag(nI, -spin) * (epsilon_kvec(p-q-k) + epsilon_kvec(p+q)))

    end function three_body_transcorr_fac

    function get_3_body_helement_ks_hub(nI, ex, tpar) result(hel)
        ! the 3-body matrix element.. here i have to be careful about 
        ! the sign and stuff.. and also if momentum conservation is 
        ! fullfilled .. 
        integer, intent(in) :: nI(nel), ex(2,3)
        logical, intent(in) :: tpar
        HElement_t(dp) :: hel
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_3_body_helement_ks_hub"
#endif
        integer :: ms_elec, ms_orbs, opp_elec, opp_orb, par_elecs(2), par_orbs(2)

        hel = h_cast(0.0_dp)

        ms_elec = sum(get_spin_pn(ex(1,:)))
        ms_orbs = sum(get_spin_pn(ex(2,:)))

        ! check spin:
        if (.not.(ms_elec == ms_orbs)) return 
        if (.not.(ms_elec == 1 .or. ms_elec == -1)) return 

        ! check momentum conservation: 
        if (.not. check_momentum_sym(ex(1,:),ex(2,:))) return

        ! i have to get the correct momenta for the epsilon contribution 
        ! see the sheets for the k-vec relations. 
        ! we need the momentum p + (s-b) or variations thereof.. 
        ! p, we know, since it is the momentum of the minority spin-electron
        ! and (s-b) is the difference of an majority spin electron and the 
        ! fitting hole, so the total momentum conservation a + b + c = s + p + q
        ! is fulfilled
        ! the same ofc is a + (c-q) 
        ! is the minority hole always in ex(2,1)? otherwise we have to find it
        ! it is not! 

        ! i think i have figured it out with the help of Manu 
        ! the k-vector of the minority spin is always involved 
        ! but of the electron.. or can we transform it? 
        ! any way we have to calculate 
        ! W(k_p + k_s - k_b) - W(k_p + k_q - k_b) 
        ! for this we have to figure out what the minority and majority 
        ! electrons are! 
        opp_elec = find_minority_spin(ex(1,:))

        ! although i really can't be sure about the minority whole always 
        ! being at the first position in ex(2,:).. 
        opp_orb = find_minority_spin(ex(2,:)) 

        par_elecs = pack(ex(1,:),ex(1,:) /= opp_elec)
        par_orbs = pack(ex(2,:),ex(2,:) /= opp_orb)

        ! i hope it is fine if i always take par_orbs(1).. this has to do 
        ! with the overal sign i guess.. so maybe i should check if 
        ! if ex() is correctly sorted.. todo
        hel = three_body_prefac * (&
              epsilon_kvec(G1(opp_elec)%k + G1(par_elecs(1))%k - G1(par_orbs(1))%k) & 
            - epsilon_kvec(G1(opp_elec)%k + G1(par_elecs(2))%k - G1(par_orbs(1))%k))

        if (tpar) hel = -hel

    end function get_3_body_helement_ks_hub

    function find_minority_spin(spins) result(opp) 
        ! for now, given 3 spin orbitals, where 2 share a spin, this function
        ! returns the opposite spin
        integer, intent(in) :: spins(:) 
        integer :: opp
#ifdef __DEBUG
        character(*), parameter :: this_routine = "find_minority_spin" 
#endif 
        integer :: i

        ! for now, do it only for 3 inputted spins:
        ASSERT(size(spins) == 3)
        ASSERT(sum(get_spin_pn(spins)) == -1 .or. sum(get_spin_pn(spins)) == 1)

        opp = -1

        if (sum(get_spin_pn(spins)) == -1) then 
            ! then we have 2 beta and 1 alpha 
            do i = 1, size(spins)
                if (is_alpha(spins(i))) opp = spins(i)
            end do
        else
            ! we have 2 alpha and 1 beta 
            do i = 1, size(spins)
                if (is_beta(spins(i))) opp = spins(i)
            end do
        end if
    end function find_minority_spin

    logical function check_momentum_sym(elecs, orbs) 
        ! routine to check the momentum conservation for double and triple 
        ! spawns 
        ! although this could in fact be used for a general check of 
        ! symmetry adaptability
        integer, intent(in) :: elecs(:), orbs(:) 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "check_momentum_sym"
#endif 
        type(BasisFN) :: ka, kb 
        integer :: i

        ! make this function flexible..
!         ASSERT(size(elecs) == size(orbs))

        call SetupSym(ka)
        call SetupSym(kb) 

        do i = 1, size(elecs)
            call AddElecSym(elecs(i), G1, nBasisMax, ka)
        end do
        do i = 1, size(orbs)
            call AddElecSym(orbs(i), G1, nBasisMax, kb)
        end do
        
        ! apply periodic BC:
        call RoundSym(ka, nBasisMax)
        call RoundSym(kb, nBasisMax)

        ! and check sym:
        check_momentum_sym = (lChkSym(ka, kb)) 

    end function check_momentum_sym

    subroutine make_triple(nI, nJ, elecs, orbs, ex, tPar)
        integer, intent(in) :: nI(nel), elecs(3), orbs(3)
        integer, intent(out) :: nJ(nel), ex(2,3)
        logical, intent(out) :: tPar
#ifdef __DEBUG
        character(*), parameter :: this_routine = "make_triple"
#endif
        integer :: sort_elecs(3), sort_orbs(3), src(3), pos_moved, k, i

        ! i should also check if this excitation is actually possible! 
        ! and which spin i move to where or?? 

        ! figure out how to do triples efficiently.. 
        ! can we do a single and then a double? 

        ! NO: talk to Manu and do this with integer representation! not
        ! with nI and nJ, since this can be done way more effective as 
        ! via the occupied orbitals.. 

        ! TODO: thats a super strange convention, .. talk with Ali and 
        ! Simon about that.. but for now accept it as it is.. 

        sort_elecs = sort_unique(elecs)
        sort_orbs = sort_unique(orbs)
        
        src = nI(sort_elecs)

        ex(1,:) = src
        ex(2,:) = sort_orbs

        nJ = nI 

        ASSERT(sum(get_spin_pn(src)) == sum(get_spin_pn(orbs)))

        ! i should do some stuff depending on the order of src and tgt 
        ! we have to check if electrons hop over other electrons, so we 
        ! might have to change the indexing to adapt to the change in nJ! 

        ! or check it individually: 
        if (src(2) < sort_orbs(1)) then 
            ! then i hops over j: 
            sort_elecs(2) = sort_elecs(2) - 1 
        end if
        if (src(3) < sort_orbs(1)) then 
            ! then i hops over k 
            ! (note: this also implies that j hops over k, but treat that 
            ! seperately below, to also cover the case, where this if here 
            ! is not fullfilled!) 
            sort_elecs(3) = sort_elecs(3) - 1 
        end if 
        if (src(3) < sort_orbs(2)) then 
            ! then j hops over k 
            sort_elecs(3) = sort_elecs(3) - 1
        end if

        pos_moved = 0 

        do k = 1, 3 
            if (src(k) < sort_orbs(k)) then 
                if (sort_elecs(k) == nel) then 
                    ! this can only happen for k == 3 
                    i = nel + 1
                    nJ(nel) = sort_orbs(k) 
                else 
                    do i = sort_elecs(k) + 1, nel 
                        if (sort_orbs(k) < nJ(i)) then 
                            nJ(i-1) = sort_orbs(k)
                            exit 
                        else 
                            nJ(i-1) = nJ(i)
                        end if
                    end do
                    if (i == nel + 1) then 
                        nJ(nel) = sort_orbs(k)
                    end if
                end if
            else 
                if (sort_elecs(k) == 1) then 
                    i = 0
                    nJ(1) = sort_orbs(k)
                else 
                    do i = sort_elecs(k)-1, 1, -1
                        if (sort_orbs(k) > nJ(i)) then 
                            nJ(i+1) = sort_orbs(k)
                            exit 
                        else 
                            nJ(i+1) = nJ(i)
                        end if
                    end do
                    if (i == 0) then 
                        nJ(1) = sort_orbs(k)
                    end if
                end if
            end if

            pos_moved = pos_moved + sort_elecs(k) - i + 1

        end do

        tPar = btest(pos_moved, 0)
        
    end subroutine make_triple

    subroutine setup_tmat_k_space(in_lat)
        ! routine which sets up the (diagonal) t-matrix in the k-space 
        ! the dimensionality and connectivity of nearest and next-nearest 
        ! neighbors influences that! 
        class(lattice), intent(in), optional :: in_lat 
        character(*), parameter :: this_routine = "setup_tmat_k_space" 

        if (present(in_lat)) then 
            if (all(nBasisMax == 0)) then 
                call setup_nbasismax(in_lat)
            end if

            if (.not. associated(G1)) then 
                call setup_g1(in_lat) 
            end if
            ! else assume it is already setup correctly 

            if (.not. associated(tmat2d)) then 
                ! call the already implemented hubbard tmat calculator.. 
                ! for now only.. in the future we should do this standalone
                call CALCTMATHUB(in_lat%get_nsites()*2, nBasisMax, bhub, & 
                    ttilt,G1,.not. in_lat%is_k_space(), in_lat%is_periodic())

            end if

        else 
            call Stop_All(this_routine, "not yet implemented!")
        end if

    end subroutine setup_tmat_k_space

end module k_space_hubbard
