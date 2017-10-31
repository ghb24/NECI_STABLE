#include "macros.h" 

module tJ_model 

    use SystemData, only: bhub, nel, nbasis
    use constants, only: dp, n_int, EPS, bits_n_int
    use real_space_hubbard, only: lat
    use bit_rep_data, only: NIfTot
    use umatcache, only: gtid
    use util_mod, only: binary_search_first_ge
    use OneEInts, only: GetTMatEl
    implicit none 

    logical :: t_tJ_model = .false.
    real(dp) :: exchange_j = 1.0_dp

    real(dp), allocatable :: exchange_matrix(:,:) 

contains 

    subroutine init_tJ_model 

        print *, "initializing tJ-model with parameters: "
        print *, "t: ", bhub
        print *, "J: ", exchange_j 

    end subroutine init_tJ_model 

    subroutine gen_excit_tJ_model (nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                      ex, tParity, pGen, hel, store, run)
        
        use SystemData, only: nel
        use bit_rep_data, only: NIfTot
        use FciMCData, only: excit_gen_store_type
        use constants, only: n_int, dp, bits_n_int
        use get_excit, only: make_single, make_double
        use back_spawn, only: make_ilutJ
        use dsfmt_interface, only: genrand_real2_dsfmt

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

        character(*), parameter :: this_routine = "gen_excit_tJ_model"

        integer :: elec, src, id, ind, elec_2, tgt_1, tgt_2
        real(dp) :: p_elec, p_orb, cum_sum, r, cum_sum_opp, cpt_opp
        integer, allocatable :: neighbors(:), ic_list(:), tmp_ic_list(:)
        real(dp), allocatable :: cum_arr(:), cum_arr_opp(:)

        ! the idea for the tJ excitation generator on a lattice is to 
        ! still pick the electron at random and then check for its 
        ! neighbors, if the neighboring site is empty, do a hop 
        ! and if there is a electron of opposite spin do a spin-flip 
        ! if it is occupied by an electron of same spin no excitation with 
        ! this neighbor is possible 

        ! use the lattice type like in the real-space hubbard implementation
        ASSERT(associated(lat))

        elec = 1 + int(genrand_real2_dsfmt() * nel) 

        p_elec = 1.0_dp / real(nel, dp) 

        src = nI(elec) 
        id = gtid(src) 

        neighbors = lat%get_neighbors(id) 

        call create_cum_list_tJ_model(ilutI, src, neighbors, cum_arr, cum_sum, & 
            ic_list) 

        if (cum_sum < EPS) then 
            nJ(1) = 0
            pgen = 0.0_dp
            return 
        end if

        r = genrand_real2_dsfmt() * cum_sum 
        
        ind = binary_search_first_ge(cum_arr, r) 

        ! i just realised that for spin-flip excitation we have to take the 
        ! opposite order of orbital picking into account too.. 
        ! because the same excitation could have been picked with the spin 
        ! opposite electron in the neighborhood of the first electron too.. 

        ic = ic_list(ind) 

        ! for spin-flips i have to add a contribution down below!
        if (ind == 1) then 
            p_orb = cum_arr(1) / cum_sum 
        else 
            p_orb = (cum_arr(ind) - cum_arr(ind-1)) / cum_sum 
        end if 

        if (ic == 1) then 
            if (is_beta(src)) then 
                tgt_1 = 2 * neighbors(ind) - 1
            else 
                tgt_1 = 2 * neighbors(ind)
            end if

            call make_single(nI, nJ, elec, tgt_1, ex, tParity) 

            ilutJ =  make_ilutJ(ilutI, ex, 1)

        else if (ic == 2) then 
            ASSERT(get_beta(src) .neqv. get_beta(spin_orb))

            ! here i have to recalc the contribution if i would have picked 
            ! the electron in orbital spin_orb first 

            if (is_beta(src)) then 
                ! need to the the index of electron 2 in nI 
                ! the second electron must be alpha 
                elec_2 = find_elec_in_ni(nI, 2*neighbors(ind))
                ! we need the orbital alpha of src 
                ! and the beta of the second orbital
                tgt_1 = get_alpha(src) 
                tgt_2 = 2 * neighbors(ind) - 1
            else 
                ! v.v here
                elec_2 = find_elec_in_ni(nI, 2*neighbors(ind) - 1)

                tgt_1 = get_beta(src)
                tgt_2 = 2 * neighbors(ind)

            end if

            call create_cum_list_tJ_model(ilutI, nI(elec_2), &
                lat%get_neighbors(gtid(nI(elec_2))), cum_arr_opp, cum_sum_opp, & 
                tmp_ic_list, tgt_1, cpt_opp)

            p_orb = p_orb + cpt_opp

            call make_double(nI, nJ, src, elec_2, tgt_1, tgt_2, ex, tParity)

            ilutJ = make_ilutJ(ilutI, ex, 2)

        else 
            ! something went wrong.. 
            call stop_all(this_routine, &
                "something went wrong ic > 2!")
        end if

        pgen = p_elec * p_orb
            
    end subroutine gen_excit_tJ_model

    subroutine create_cum_list_tJ_model(ilutI, src, neighbors, cum_arr, cum_sum, &
            ic_list, tgt, cpt)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: src, neighbors(:)
        real(dp), intent(out), allocatable :: cum_arr(:)
        real(dp), intent(out) :: cum_sum
        integer, intent(out), allocatable :: ic_list(:) 
        integer, intent(in), optional :: tgt
        real(dp), intent(out), optional :: cpt
#ifdef __DEBUG
        character(*), parameter :: this_routine = "create_cum_list_tJ_model"
#endif
        integer :: i
        integer, allocatable :: single_excits(:)
        integer, allocatable :: spin_flips(:)
        real(dp) :: elem
        logical :: t_single, t_flip

        ASSERT(IsOcc(ilutI,src))

        allocate(cum_arr(size(neighbors)))
        allocate(ic_list(size(neighbors)))
        allocate(single_excits(size(neighbors)))
        allocate(spin_flips(size(neighbors)))
        cum_arr = 0
        cum_sum = 0.0_dp
        ic_list = 0

        if (is_beta(src)) then 
            single_excits = 2*neighbors - 1
            spin_flips = 2 * neighbors
        else 
            single_excits = 2 * neighbors
            spin_flips = 2 * neighbors - 1
        end if

        if (present(tgt)) then 
            t_single = .false. 
            t_flip = .false. 

            ! find the probability of choosing orbital target
            if (is_beta(src) .eqv. is_beta(tgt)) then 
                ! then it was definetly a single excitation 
                t_single = .true.
            else
                t_flip = .true.
            end if

            ASSERT(present(cpt))
            do i = 1, ubound(neighbors,1)
                if (IsNotOcc(ilutI, single_excits(i)) .and. &
                    IsNotOcc(ilutI, spin_flips(i))) then 
                    ! just to be sure use the tmat, so both orbitals are 
                    ! definetly connected
                    elem = abs(GetTMatEl(src, single_excits(i)))

                else if (IsOcc(ilutI, spin_flips(i)) .and. &
                         IsNotOcc(ilutI, single_excits(i))) then 
                     elem = abs(get_heisenberg_exchange(src, spin_flips(i)))

                 else 
                     elem = 0.0_dp
                 end if
                 cum_sum = cum_sum + elem 

                 if (t_single .and. tgt == single_excits(i)) then 
                     cpt = elem
                 else if (t_flip .and. tgt == spin_flips(i)) then 
                     cpt = elem
                 end if
             end do
             if (cum_sum < EPS) then 
                 cpt = 0.0_dp
             else
                 cpt = cpt / cum_sum
             end if
        else
            ! create the list depending on the possible excitations 
            do i = 1, ubound(neighbors,1)
                if (IsNotOcc(ilutI,single_excits(i)) .and. &
                    IsNotOcc(ilutI, spin_flips(i))) then 
                    ! then the orbital is empty an we can do a hopping 
                    ! reuse the hubbard matrix elements..
                    ! or just the -t element? 
                    ! since we only need absolute value of matrix element 
                    ! it would be better to just use -t .. anyway.. keep it 
                    ! general
                    elem = abs(GetTMatEl(src, single_excits(i)))
                    ic_list(i) = 1

                else if (IsOcc(IlutI,spin_flips(i)) .and. &
                         IsNotOcc(IlutI,single_excits(i))) then 
                     ! then we can do a spin flip 
                     elem = abs(get_heisenberg_exchange(src, spin_flips(i)))
                     ic_list(i) = 2

                 else 
                     ! if the spin-parallel is occupied, no exciation 
                     ! possible, and also prohibit double occupancies
                     elem = 0.0_dp
                 end if

                    cum_sum = cum_sum + elem 
                    cum_arr(i) = cum_sum
             end do
         end if

    end subroutine create_cum_list_tJ_model

    function calc_pgen_tJ_model(nI, ilutI, ex, ic) result(pgen)
        integer, intent(in) :: nI(nel), ex(2,2), ic
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        real(dp) :: pgen
#ifdef __DEBUG
        character(*), parameter :: this_routine = "calc_pgen_tJ_model"
#endif

        integer :: src(2), tgt(2), tgt_orb
        real(dp) :: p_elec, p_orb, cum_sum, cpt_1, cpt_2
        real(dp), allocatable :: cum_arr(:)
        integer, allocatable :: tmp_list(:) 

        ASSERT(ic >= 0 )

        if (ic == 0 .or. ic > 2) then 
            pgen = 0.0_dp
            return
        end if

        src = get_src(ex)
        tgt = get_tgt(ex) 

        ASSERT(associated(lat))

        p_elec = 1.0_dp / real(nel, dp)

        if (ic == 1) then 
            ! here it is easy.. 
            ASSERT(is_beta(src(1)) .eqv. is_beta(tgt(1)))
            ASSERT(any(tgt(1) == lat%get_spinorb_neighbors(src(1))))
            ASSERT(IsOcc(ilutI, src(1)))
            ASSERT(IsNotOcc(ilutI,tgt(1)))

            call create_cum_list_tJ_model(ilutI, src(1), lat%get_neighbors(gtid(src(1))), & 
                cum_arr, cum_sum, tmp_list, tgt(1), p_orb)

        else if (ic == 2) then 
            ! here we have to find the correct orbital.. 
            ! and i just realised that we maybe have to take into account 
            ! of having picked the orbitals in a different order.. 
            ASSERT(is_beta(src(1)) .neqv. is_beta(src(2)))
            ASSERT(is_beta(tgt(1)) .neqv. is_beta(tgt(2)))
            ASSERT(is_in_pair(src(1),tgt(1)) .or. is_in_pair(src(1),tgt(2)))
            ASSERT(is_in_pair(src(2),tgt(1)) .or. is_in_pair(src(2),tgt(2)))

            if (is_beta(src(1)) .eqv. is_beta(tgt(1))) then 
                ! then those to orbitls were chosen 
                ASSERT(is_beta(src(2)) .eqv. is_beta(tgt(2)))

                call create_cum_list_tJ_model(ilutI, src(1), lat%get_neighbors(gtid(src(1))), &
                    cum_arr, cum_sum, tmp_list, tgt(1), cpt_1)
                call create_cum_list_tJ_model(ilutI, src(2), lat%get_neighbors(gtid(src(2))), &
                    cum_arr, cum_sum, tmp_list, tgt(2), cpt_2)

            else if (is_beta(src(1)) .eqv. is_beta(tgt(2))) then 
                ASSERT(is_beta(src(2)) .eqv. is_beta(tgt(1)))

                call create_cum_list_tJ_model(ilutI, src(1), lat%get_neighbors(gtid(src(1))), &
                    cum_arr, cum_sum, tmp_list, tgt(2), cpt_1)
                call create_cum_list_tJ_model(ilutI, src(2), lat%get_neighbors(gtid(src(2))), & 
                    cum_arr, cum_sum, tmp_list, tgt(1), cpt_2)
#ifdef __DEBUG
            else 
                ! something went wrong 
                call stop_all(this_routine, "something went wrong!")
#endif
            end if 

            p_orb = cpt_1 + cpt_2 

        end if

        pgen = p_elec * p_orb

    end function calc_pgen_tJ_model

    function find_elec_in_ni(nI, orb) result(elec)
        ! routine to find the number of the elctron in spin-orbital orb
        integer, intent(in) :: nI(nel), orb
        integer :: elec 
#ifdef __DEBUG
        character(*), parameter :: this_routine = "find_elec_in_ni"
#endif

        ASSERT(orb > 0)
        ASSERT(orb <= nbasis) 

        ! can i just reuse 
        elec = binary_search_first_ge(nI, orb)

        ! it already make the fail case.. or?
        ! and then make a fail-case: 
        if (nI(elec) /= orb) then 
            ! it is actually not found? 
            elec = -1 
        end if

    end function find_elec_in_ni

    function get_heisenberg_exchange(src, tgt) result(hel) 
        ! this is the wrapper function to get the heisenberg exchange 
        ! contribution. which will substitute get_umat in the 
        ! necessary places.. 
        ! NOTE: only in the debug mode the neighboring condition is tested
        ! also that that both orbitals src and tgt are occupied by opposite 
        ! spins have to be checked before this function call! 
        integer, intent(in) :: src, tgt 
        HElement_t(dp) :: hel 
#ifdef __DEBUG 
        character(*), parameter :: this_routine = "get_heisenberg_exchange"

        ASSERT(src > 0) 
        ASSERT(src <= nbasis)
        ASSERT(tgt > 0)
        ASSERT(tgt <= nbasis)

        ASSERT(associated(lat)) 
        ASSERT(is_beta(src) .neqv. is_beta(tgt)) 
        if (is_beta(src)) then 
            ASSERT(any(tgt == lat%get_spinorb_neighbors(src)-1))
        else 
            ASSERT(any(tgt == lat%get_spinorb_neighbors(src)+1))
        end if
        ASSERT(allocated(exchange_matrix))
#endif

        hel = h_cast(exchange_matrix(src, tgt))

    end function get_heisenberg_exchange

end module tJ_model
