#include "macros.h"

module back_spawn_excit_gen

    implicit none

contains

    ! also write a wrapper-like routine for an excitation generator if 
    ! back-spawn is activated.. to not mess up all the old functions too much. 
    subroutine gen_excit_back_spawn(nI, ilutI, nJ, ilutJ, exFlag, ic, &
            ExcitMat, tParity, pgen, HelGen, store) 
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ic, ExcitMat(2,2) 
        logical, intent(out) :: tParity 
        real(dp), intent(out) :: pgen 
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store 
        character(*), parameter :: this_routine = "gen_excit_back_spawn"

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
        integer, intent(out) : nJ(nel), ex(2,2)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tPar
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "gen_single_back_spawn"

        ! depending on the method we pick electrons accordingly
        if (t_back_spawn_flex) then 
            elec = i + floor(genrand_real2_dSFMT() * nel)

            call check_electron_location([nI(elec),0], 1, loc)

            pgen_elec = 1.0_dp/real(nel, dp)
        else
            call pick_virtual_electrons_double(nI, elec, pgen_elec)
        end if

        src = nI(elec)

        cc_index = ClassCountInd (get_spin(src), SpinOrbSymLabel(src), &
                                  G1(src)%Ml)


        ! i have to make the logic easier here at some point.. 
        ! we just need to do some extensive testing and decide on one 
        ! back-spawn method and abandon the rest.. 
        if (t_back_spawn_occ_virt) then 
            call pick_occupied_orbital_single(nI, src, cc_index, pgen, tgt)

        else if (t_back_spawn_flex) then 

            if (loc == 2 .and. occ_virt_level /= -1) then 
                call pick_occupied_orbital_single(nI, src, cc_index, pgen, tgt) 

            else 
                if (occ_virt_level == 2) then 
                    call pick_occupied_orbital_single(nI, src, cc_index, pgen, tgt) 
                else
                    tgt = select_orb_sing(nI, ilutI, src, cc_index, pgen) 
                end if 
            end if
        else
            tgt = select_orb_sing(nI, ilutI, src, cc_index, pgen) 
        end if

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

        if (t_back_spawn) then 
            call pick_virtual_electrons_double(nI, elecs, src, sym_product, ispn,&
                                                sum_ml, pgen)

            ! TODO!! thats the problem here! there are circular dependencies 
            ! when I call that from here.. so i guess i need to setup the 
            ! back_spawn excitation routines in a seperate file 
            sym_product = RandExcitSymLabelProd(SpinOrbSymLabel(src(1)), &
                                             SpinOrbSymLabel(src(2)))

    end subroutine gen_double_back_spawn

end module back_spawn_excit_gen
