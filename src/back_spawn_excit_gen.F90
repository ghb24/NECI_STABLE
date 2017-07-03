#include "macros.h"

module back_spawn_excit_gen

    use constants, only: dp, n_int, EPS, bits_n_int
    use SystemData, only: nel, G1, nbasis
    use bit_rep_data, only: niftot
    use SymExcitDataMod, only: excit_gen_store_type, SpinOrbSymLabel
    use bit_reps, only: test_flag, get_initiator_flag
    use FciMCData, only: pSingles, projedet, pDoubles
    use dSFMT_interface, only: genrand_real2_dSFMT
    use excit_gens_int_weighted, only: gen_single_4ind_ex, select_orb_sing, &
                                pick_weighted_elecs, get_paired_cc_ind, select_orb, &
                                pgen_select_orb
    use excit_gen_5, only: gen_double_4ind_ex2, pick_a_orb, pgen_select_a_orb
    use CalcData, only: t_back_spawn_flex, t_back_spawn_occ_virt, t_back_spawn, &
                        occ_virt_level
    use GenRandSymExcitNUMod, only: ClassCountInd, RandExcitSymLabelProd
    use back_spawn, only: check_electron_location, pick_virtual_electrons_double, & 
                          pick_occupied_orbital_single, pick_virtual_electron_single, &
                          pick_occupied_orbital, pick_second_occupied_orbital
    use get_excit, only: make_single, make_double

    implicit none


contains

    ! also write a wrapper-like routine for an excitation generator if 
    ! back-spawn is activated.. to not mess up all the old functions too much. 
    subroutine gen_excit_back_spawn(nI, ilutI, nJ, ilutJ, exFlag, ic, &
            ExcitMat, tParity, pgen, HelGen, store) 
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ic, ExcitMat(2,2) 
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity 
        real(dp), intent(out) :: pgen 
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store 
        character(*), parameter :: this_routine = "gen_excit_back_spawn"

        logical :: t_init 
        ! check the non-initiator criteria beforehand 
        t_init = test_flag(ilutI, get_initiator_flag(1)) 

        if (t_init) then 
            ! do the "normal" excitation type if it is an initiator
            if (genrand_real2_dSFMT() < pSingles) then 

                ic = 1 
                call gen_single_4ind_ex (nI, ilutI, nJ, ilutJ, ExcitMat, &
                                     tParity, pGen)
                pgen = pgen * pSingles
                
            else 

                ic = 2
                call gen_double_4ind_ex2 (nI, ilutI, nJ, ilutJ, ExcitMat, tParity, &
                                      pGen)
                pgen = pgen * pDoubles

            end if
        else 
            ! otherwise use custom made ones
            if (genrand_real2_dSFMT() < pSingles) then 
                
                ic = 1
                call gen_single_back_spawn(nI, ilutI, nJ, ilutJ, ExcitMat, &
                    tParity, pgen)
                pgen = pgen * pSingles 

            else 
                
                ic = 2 
                call gen_double_back_spawn(nI, ilutI, nJ, ilutJ, ExcitMat, &
                    tParity, pgen) 
                pgen = pgen * pDoubles 

            end if
        end if

    end subroutine gen_excit_back_spawn

    subroutine gen_single_back_spawn(nI, ilutI, nJ, ilutJ, ex, tPar, pgen) 
        ! specialised single excitation routine for the back-spawn method
        ! for the moment i still have to decide, which back-spawn method is 
        ! in use.. 
        integer, intent(in) :: nI(nel) 
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tPar
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "gen_single_back_spawn"

        integer :: elec, src, cc_index, loc, tgt
        real(dp) :: pgen_elec

        ! depending on the method we pick electrons accordingly
        if (t_back_spawn_flex) then 
            elec = 1 + floor(genrand_real2_dSFMT() * nel)

            call check_electron_location([nI(elec),0], 1, loc)

            pgen_elec = 1.0_dp/real(nel, dp)
        else
            call pick_virtual_electron_single(nI, elec, pgen_elec)
        end if

        src = nI(elec)

        cc_index = ClassCountInd (get_spin(src), SpinOrbSymLabel(src), &
                                  G1(src)%Ml)


        ! i have to make the logic easier here at some point.. 
        ! we just need to do some extensive testing and decide on one 
        ! back-spawn method and abandon the rest.. 
        ! i hope this logic is correct: otherwise check on the bottom
        if ((t_back_spawn_occ_virt) .or. (t_back_spawn_flex .and. (&
            (loc == 2 .and. occ_virt_level /= -1) .or. occ_virt_level == 2))) then 

            call pick_occupied_orbital_single(nI, src, cc_index, pgen, tgt)

        else 

            tgt = select_orb_sing(nI, ilutI, src, cc_index, pgen) 

        end if
! 
!         if (t_back_spawn_occ_virt) then 
!             call pick_occupied_orbital_single(nI, src, cc_index, pgen, tgt)
! 
!         else if (t_back_spawn_flex) then 
! 
!             if (loc == 2 .and. occ_virt_level /= -1) then 
!                 call pick_occupied_orbital_single(nI, src, cc_index, pgen, tgt) 
! 
!             else 
!                 if (occ_virt_level == 2) then 
!                     call pick_occupied_orbital_single(nI, src, cc_index, pgen, tgt) 
!                 else
!                     tgt = select_orb_sing(nI, ilutI, src, cc_index, pgen) 
!                 end if 
!             end if
!         else
!             tgt = select_orb_sing(nI, ilutI, src, cc_index, pgen) 
!         end if

        if (tgt == 0) then 
            nJ(1) = 0
            pgen = 0.0_dp 
            return
        end if

        call make_single(nI, nJ, elec, tgt, ex, tPar)

        ilutJ = ilutI
        clr_orb(ilutJ, src)
        set_orb(ilutJ, tgt)

        pgen = pgen * pgen_elec

    end subroutine gen_single_back_spawn

    subroutine gen_double_back_spawn(nI, ilutI, nJ, ilutJ, ex, tPar, pgen) 
        ! the double excitation routine for the back-spawn method 
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: tpar
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "gen_double_back_spawn"

        integer :: elecs(2), sym_product, ispn, sum_ml, src(2), loc, cc_a, cc_b,&
                   orbs(2)
        real(dp) :: int_cpt(2), cum_sum(2), sum_pair(2), cpt_pair(2), &
                    cum_arr(nbasis)

        if (t_back_spawn_flex) then 
            ! Pick the electrons in a weighted fashion
            call pick_weighted_elecs(nI, elecs, src, sym_product, ispn, sum_ml, &
                                 pgen)

            call check_electron_location(src, 2, loc)

        else
            call pick_virtual_electrons_double(nI, elecs, src, sym_product, ispn,&
                                                sum_ml, pgen)

            ! TODO!! thats the problem here! there are circular dependencies 
            ! when I call that from here.. so i guess i need to setup the 
            ! back_spawn excitation routines in a seperate file 
            sym_product = RandExcitSymLabelProd(SpinOrbSymLabel(src(1)), &
                                             SpinOrbSymLabel(src(2)))

            ! here it could be that there are no valid electrons.. 
            if (elecs(1) == 0) then 
                nJ(1) = 0
                pgen = 0.0_dp
                return
            end if

            ! indicate no flex option is used
            loc = -1

        end if

        ! at some point i have to make this whole logic easier here.. 
        ! in the end i think i have to remove some possibilities and stick to 
        ! one back-spawn method. 

        if ((t_back_spawn) .or. (t_back_spawn_flex .and. (&
            (loc == 1 .and. occ_virt_level /= -1) .or. (loc == 2) .or. &
            (loc == 0 .and. occ_virt_level >= 1)))) then 

            call pick_occupied_orbital(nI, src, ispn, int_cpt(1), cum_sum(1), &
                                        orbs(1))

        else 

            orbs(1) = pick_a_orb(ilutI, src, iSpn, int_cpt(1), cum_sum(1), cum_arr)
            
        end if

! 
!         if (t_back_spawn_temp ) then
!  
!             ! if we have one of the electrons in the occupied manifold 
!             ! pick atleast one hole also from this manifold to not increase 
!             ! the excitation level!
!            
!             call pick_occupied_orbital(nI, src, ispn, int_cpt(1), cum_sum(1), &
!                                         orbs(1))
! 
!         else if (t_back_spawn_flex ) then
!             ! now we have to decide on the flex-spawn + occ-virt implo:
!             if (loc == 1) then 
!                 ! the new option to excite on level
!                 if (occ_virt_level == -1) then
!                     orbs(1) = pick_a_orb(ilutI, src, iSpn, int_cpt(1), cum_sum(1), cum_arr)
!                 else
!                     call pick_occupied_orbital(nI, src, ispn, int_cpt(1), cum_sum(1), &
!                                         orbs(1))
!                 end if
! 
! 
!             else if (loc == 2) then 
!                 ! then we always pick an occupied orbital
!                 call pick_occupied_orbital(nI, src, ispn, int_cpt(1), cum_sum(1), &
!                                         orbs(1))
!             else 
!                 ! depending on the occ_virt_level
!                 if (occ_virt_level < 1) then 
!                     ! just pick normal:
!                     orbs(1) = pick_a_orb(ilutI, src, iSpn, int_cpt(1), cum_sum(1), cum_arr)
!                 else 
!                     ! otherwise if it is 1 or 2, we want to pick the (a) also
!                     ! from the occupied manifold
! 
!                     call pick_occupied_orbital(nI, src, ispn, int_cpt(1), cum_sum(1), &
!                                         orbs(1))
!                 end if
!             end if
!         else 
!             orbs(1) = pick_a_orb(ilutI, src, iSpn, int_cpt(1), cum_sum(1), cum_arr)
!         end if


        if (orbs(1) /= 0) then
            cc_a = ClasSCountInd(orbs(1))
            cc_b = get_paired_cc_ind(cc_a, sym_product, sum_ml, iSpn)

            if (t_back_spawn_flex .and.((loc == 2 .and. occ_virt_level /= -1) .or. &
                (occ_virt_level == 2))) then 

                call pick_second_occupied_orbital(nI, src, cc_b, orbs(1), ispn,&
                    int_cpt(2), cum_sum(2), orbs(2))

            else 

                orbs(2) = select_orb (ilutI, src, cc_b, orbs(1), int_cpt(2), &
                                  cum_sum(2))
            end if

            ASSERT((.not. (is_beta(orbs(2)) .and. .not. is_beta(orbs(1)))))
! 
!             if (t_back_spawn_flex  .and. .not. t_temp_init) then 
!                 ! in this case i have to pick the second orbital also from the 
!                 ! occupied list, but now also considering symmetries
!                 if (loc == 2) then 
!                     ! now also mix the occ-virt with this back-spawning
!                     ! for loc 2 always do it
!                     ! except specified by occ_virt_level= -1
!                     if (occ_virt_level == -1) then
!                         orbs(2) = select_orb (ilutI, src, cc_b, orbs(1), int_cpt(2), &
!                                   cum_sum(2))
!                     else
!                         call pick_second_occupied_orbital(nI, src, cc_b, orbs(1), ispn,&
!                             int_cpt(2), cum_sum(2), orbs(2))
!                     end if
! 
!                 else if (loc == 1 .and. occ_virt_level == 2) then 
!                     call pick_second_occupied_orbital(nI, src, cc_b, orbs(1), ispn,&
!                         int_cpt(2), cum_sum(2), orbs(2))
! 
!                 else if (loc == 0 .and. occ_virt_level == 2) then 
!                     call pick_second_occupied_orbital(nI, src, cc_b, orbs(1), ispn,&
!                         int_cpt(2), cum_sum(2), orbs(2))
!                 else
!                     orbs(2) = select_orb (ilutI, src, cc_b, orbs(1), int_cpt(2), &
!                                   cum_sum(2))
!                 end if
! 
!             else
! 
!                 orbs(2) = select_orb (ilutI, src, cc_b, orbs(1), int_cpt(2), &
!                                   cum_sum(2))
!             end if

        end if

        if (any(orbs == 0)) then
            nJ(1) = 0
            pgen = 0.0_dp
            return
        end if

        ! can i exit right away if this happens??
        ! i am pretty sure this means 
        if (any(cum_sum < EPS)) then 
           cum_sum = 1.0_dp
           int_cpt = 0.0_dp
        end if
        
        ! only on parallel excitations.. and symmetric exciation generator is 
        ! turned off for now in the back-spawning
        if (is_beta(orbs(1)) .eqv. is_beta(orbs(2)))  then
            if (t_back_spawn_occ_virt .or. (t_back_spawn_flex .and. (& 
                (loc == 1 .and. (occ_virt_level == 0 .or. occ_virt_level == 1)) &
                .or. (loc == 2 .and. occ_virt_level == -1) .or. &
                (loc == 0 .and. occ_virt_level == 1)))) then 

                if (any(orbs(2) == projedet(:,1))) then
                   ! if (b) is also in the occupied manifold i could have 
                    ! picked the other way around.. 
                    ! with the same uniform probability: 
                    cpt_pair(1) = int_cpt(1)
                    sum_pair(1) = cum_sum(1) 
                    ! and then (a) would have been picked according to the 
                    ! "normal" procedure
                    call pgen_select_orb(ilutI, src, orbs(2), orbs(1), &
                                 cpt_pair(2), sum_pair(2))
                else
                    ! if (b) is not in the occupied this does not work or? 
                    ! since i am forcing (a) to be in the occupied.. 
                    ! so remove this pgen:
                    cpt_pair = 0.0_dp
                    sum_pair = 1.0_dp
                end if

            else if (t_back_spawn_flex .and. ((loc == 0 .and. occ_virt_level == 2) & 
                .or. (loc == 1 .and. occ_virt_level == 2) .or. & 
                (loc == 2 .and. occ_virt_level /= -1))) then 

                cpt_pair = int_cpt 
                sum_pair = cum_sum 
                

            else

                ! otherwise "normal"
                call pgen_select_a_orb(ilutI, src, orbs(2), iSpn, cpt_pair(1), &
                                       sum_pair(1), cum_arr, .false.)
                call pgen_select_orb(ilutI, src, orbs(2), orbs(1), &
                                     cpt_pair(2), sum_pair(2))

            end if 

        else 

            cpt_pair = 0.0_dp
            sum_Pair = 1.0_dp
                
        end if

        if (any(sum_pair < EPS)) then 
            cpt_pair = 0.0_dp
            sum_pair = 1.0_dp
        end if

        pgen = pgen * (product(int_cpt) / product(cum_sum) + &
                       product(cpt_pair) / product(sum_pair))

        ! And generate the actual excitation
        call make_double (nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), &
                          ex, tpar)
        ilutJ = ilutI
        clr_orb (ilutJ, src(1))
        clr_orb (ilutJ, src(2))
        set_orb (ilutJ, orbs(1))
        set_orb (ilutJ, orbs(2))


    end subroutine gen_double_back_spawn

end module back_spawn_excit_gen
