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
                    nel, tHPHF, nOccBeta, nOccAlpha, nbasis, tLatticeGens
    use lattice_mod, only: get_helement_lattice_ex_mat, get_helement_lattice_general, &
                           determine_optimal_time_step
    use procedure_pointers, only: get_umat_el, generate_excitation
    use gen_coul_ueg_mod, only: get_hub_umat_el
    use constants, only: n_int, dp, EPS, bits_n_int
    use bit_rep_data, only: NIfTot
    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet
    use real_space_hubbard, only: lat_tau_factor
    use fcimcdata, only: tsearchtau, tsearchtauoption
    use CalcData, only: tau, t_hist_tau_search, t_hist_tau_search_option
    use dsfmt_interface, only: genrand_real2_dsfmt
    use util_mod, only: binary_search_first_ge
    use back_spawn, only: make_ilutJ, get_orb_from_kpoints, is_allowed_ueg_k_vector, &
                          get_ispn
    use FciMCData, only: excit_gen_store_type
    use get_excit, only: make_double
    use UmatCache, only: gtid
    use OneEInts, only: GetTMatEl
    use sltcnd_mod, only: sltcnd_0

    implicit none 

    integer, parameter :: ABORT_EXCITATION = 0

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

        t_hist_tau_search = .false. 
        t_hist_tau_search_option = .false.

    end subroutine init_k_space_hubbard

    subroutine check_k_space_hubbard_input()
        character(*), parameter :: this_routine = "check_k_space_hubbard_input"

        print *, "checking input for k-space hubbard:" 
        !todo: find the incompatible input and abort here!

        print *, "input is fine!"

    end subroutine check_k_space_hubbard_input

    subroutine gen_excit_k_space_hub (nI, ilutI, nJ, ilutJ, exFlag, ic, &
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
        character(*), parameter :: this_routine = "gen_excit_k_space_hub"
#endif
        real(dp) :: p_elec, p_orbs
        integer :: elecs(2), orbs(2), src(2)

        ! i first have to choose an electron pair (ij) at random 
        ! but with the condition that they have to have opposite spin! 
        call pick_spin_opp_elecs(nI, elecs, p_elec) 

        src = nI(elecs)

        call pick_ab_orbitals_hubbard(nI, ilutI, src, orbs, p_orbs)

        if (orbs(1) == ABORT_EXCITATION) then 
            nJ(1) = ABORT_EXCITATION
            pgen = 0.0_dp
            return 
        end if

        ! and make the excitation 
        call make_double(nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), ex, tParity)

        ilutJ = make_ilutJ(ilutI, ex, 2) 

        pgen = p_elec * p_orbs

    end subroutine gen_excit_k_space_hub

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

    subroutine pick_ab_orbitals_hubbard(nI, ilutI, src, orbs, p_orbs) 
        ! depending on the already picked electrons (ij) pick an orbital 
        ! (a) and the connected orbital (b)
        integer, intent(in) :: nI(nel), src(2)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: orbs(2) 
        real(dp), intent(out) :: p_orbs
#ifdef __DEBUG
        character(*), parameter :: this_routine = "pick_ab_orbitals_hubbard"
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
        call pick_from_cum_list(cum_arr, cum_sum, ind, p_orbs)

        orbs = orb_list(ind,:)

        ! do i have to recalc. the pgen the other way around? yes! 
        ! effectively reuse the above functionality
        ! i am pretty sure i just have to find the position in the 
        ! list.. OR: since in the hubbard it is just twice the 
        ! probability or? i am pretty sure yes.. but for all of them.. 
        ! so in the end it shouldnt matter again..
        p_orbs = 2.0_dp * p_orbs

    end subroutine pick_ab_orbitals_hubbard

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

        if (present(tgt)) then 
            ASSERT(present(cpt))

            cpt = 0.0_dp

            do a = 1, nbasis
                elem = 0.0_dp
                ! if a is empty
                if (IsNotOcc(ilutI, a)) then 
                    b = get_orb_from_kpoints(src(1), src(2), a)

                    ! and b is empty and not a
                    if (b /= a .and. IsNotOcc(ilutI, b)) then
                        ! is it sure that we have opposite spin?
                        ASSERT(.not. same_spin(a,b))

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

    function calc_pgen_k_space_hubbard(nI, ilutI, ex, ic) result(pgen) 
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(in) :: nI(nel), ex(2,2), ic
        real(dp) :: pgen
#ifdef __DEBUG
        character(*), parameter :: this_routine = "calc_pgen_k_space_hubbard"
#endif
        real(dp) :: p_elec, p_orb, cum_arr(nbasis), cum_sum
        integer :: orb_list(nbasis,2), src(2)

        if (ic /= 2) then 
            pgen = 0.0_dp 
            return 
        end if

        p_elec = 1.0_dp / real(nOccBeta * nOccAlpha, dp) 

        src = get_src(ex)

        call create_ab_list_hubbard(nI, ilutI, src, orb_list, cum_arr, cum_sum, & 
                ex(2,1), p_orb) 

        ! i do not need to recalc, the p(b|ij) since it is the same.. 
        ! but i need a factor of 2 somewhere.. figure that out!
        pgen = p_elec * p_orb

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

        if (present(ic_ret)) then 
            if (ic_ret == 0) then 
                hel = sltcnd_0(nI) 

            else if (ic == 2) then 
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
        integer :: i, j

        ! just sum up the orbital energies of the occupied orbitals.. or?
        do i = 1, nel 
            hel = hel + GetTMatEl(nI(i),nI(i))
        end do

        ! so why not just use the sltcnd_0?
        ! is there a 2-body influence on the matrix elements? 

    end function get_diag_helement_k_sp_hub

    function get_offdiag_helement_k_sp_hub(nI, ex, tpar) result(hel) 
        integer, intent(in) :: nI(nel), ex(2,2)
        logical, intent(in) :: tpar 
        HElement_t(dp) :: hel 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "get_offdiag_helement_k_sp_hub"
#endif
        integer :: src(2), tgt(2), ij(2), ab(2)

        src = get_src(ex)
        tgt = get_tgt(ex)
        if (same_spin(src(1),src(2)) .or. same_spin(tgt(1),tgt(2))) then 
            hel = h_cast(0.0_dp)
            return 
        end if

        ij = gtid(src)
        ab = gtid(tgt) 
        ! that about the spin?? must spin(a) be equal spin(i) and same for 
        ! b and j? does this have an effect on the sign of the matrix element? 
        !todo!

        hel = get_hub_umat_el(ab(1),ab(2),ij(1),ij(2))

        if (t_trans_corr) then 
            ! do something 
            hel = hel * exp(trans_corr_param/2.0_dp * & 
                (GetTMatEl(ij(1),ij(1)) + GetTMatEl(ij(2),ij(2))  & 
                - GetTMatEl(ab(1),ab(1)) - GetTMatEl(ab(2),ab(2)))) 
        end if

        if (t_trans_corr_2body) then 
            ! do something 2-body.. 
        end if

        if (tpar) hel = -hel 

    end function get_offdiag_helement_k_sp_hub

end module k_space_hubbard
