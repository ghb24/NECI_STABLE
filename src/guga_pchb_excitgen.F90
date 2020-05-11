#include "macros.h"
module guga_pchb_excitgen

    implicit none

    ! start with one sampler for now!
    type(aliasSamplerArray_t) :: pchb_samplers(1)
    integer, allocatable :: tgt_orbs(:,:)

contains

    subroutine gen_rand_excit_pchb_guga(nI, ilutI, nJ, ilutJ, exFlag, ic, ex, &
            tpar, pgen, helgen, store, part_type)
        integer, intent(in) :: nI(nel), exFlag

        integer(n_int), intent(in) :: ilutI(0:IlutBits%len_tot)
        integer, intent(out) :: nJ(nel), ic, ex(2,maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:IlutBits%len_tot)
        logical, intent(out) :: tpar
        real(dp), intent(out) :: pgen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = "gen_rand_excit_pchb_guga"

        integer(n_int) :: ilut(0:GugaBits%len_tot), excitation(0:GugaBits%len_tot)
        integer :: excit_typ(2)

        unused_var(exFlag)
        unused_var(part_type)
        unused_var(store)
        unused_var(tpar)

        ! initialize output variables
        nJ = 0
        ilutJ = 0_n_int
        ic = -1
        helgen = h_cast(0.0_dp)
        ex = -1

        call convert_ilut_toGUGA(ilutI, ilut)
        ASSERT(isProperCSF_ilut(ilut))

        if (tNewDet) then
            ! use new setup function for additional CSF informtation
            ! instead of calculating it all seperately..
            call init_csf_information(ilut(0:nifd))

            ! then set tNewDet to false and only set it after the walker loop
            ! in FciMCPar
            tNewDet = .false.

        end if


        if (genrand_real2_dSFMT() < pSingles) then
            ic = 1
            call generate_single_pchb_guga(ilut, nI, excitation, pgen)
            pgen = pgen * pSingles
        else
            ic = 2
            call generate_double_pchb_guga(ilut, nI, excitation, pgen, excit_typ)

            pgen = pgen * pDoubles

            if (IsNullDet(nJ)) then
                nInvalidExcits = nInvalidExcits + 1
            else
                nValidExcits = nValidExcits + 1
            end if
        end if

        ! check if excitation generation was successful
        if (near_zero(pgen) .or. any(nj == 0)) then
            ! indicate NullDet to skip spawn step
            nJ(1) = 0
            pgen = 0.0_dp
            HElGen = h_cast(0.0_dp)
            return
        end if

        if (t_matele_cutoff .and. abs(HElGen) < matele_cutoff) then
            HElgen = h_cast(0.0_dp)
            nJ(1) = 0
            pgen = 0.0_dp
            return
        end if

        ! also store information on type of excitation for the automated
        ! tau-search for the non-weighted guga excitation generator in
        ! the excitMat variable
        ex(1,1:2) = excit_typ

        ! profile tells me this costs alot of time.. so remove it
        ! and only do it in debug mode..
        ! i just have to be sure that no wrong csfs are created..

        ASSERT(isProperCSF_ilut(excitation, .true.))
        ! otherwise extract H element and convert to 0

        call convert_ilut_toNECI(excitation, ilutJ, HElgen)
        call decode_bit_det(nJ, ilutJ)

    end subroutine gen_rand_excit_pchb_guga

    subroutine generate_double_pchb_guga(nI, ilutI, nJ, ilutJ, ex, pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:IlutBits%len_tot)
        integer, intent(out) :: nJ(nel), ex
        integer(n_int), intent(out) :: ilutJ(0:IlutBits%len_tot)
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "generate_double_pchb_guga"

        ! init
        nJ = 0
        ilutJ = 0_n_int
        ex = -1
        pgen = 0.0_dp

    end subroutine generate_double_pchb_guga

    subroutine init_guga_pchb_excitgen
        character(*), parameter :: this_routine = "init_guga_pchb_excitgen"

        ! try for now to use an uniform single excitation generator too!
        generate_single_pchb_guga => gen_uniform_guga_single

    end subroutine init_guga_pchb_excitgen

    subroutine gen_uniform_guga_single(ilut, nI, excitation, pgen)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(out) :: excitation(0:GugaBits%len_tot)
        real(dp) :: pgen

        unused_var(ilut)

        ! init for safety and so we can return on abort
        excitation = 0_n_int
        pgen = 0.0_dp

        ! pick random electron
        ! have to modify pgen for doubly occupied orbs! since double the
        ! chance!
        elec = 1 + floor(genrand_real2_dSFMT() * nel)

        call pick_uniform_spatial_hole(nI, elec, excitInfo, pgen)

    end subroutine gen_uniform_guga_single

    subroutine pick_uniform_spatial_hole(nI, elec, excitInfo, pgen)
        integer, intent(in) :: nI, elec
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        debug_function_name("pick_uniform_spatial_hole")

        pgen = 0.0_dp

        ! electron spin-orbital
        so_elec = nI(elec)

        ! electron spatial index
        s_elec = gtID(so_elec)

        ! get the symmetry index for this electron
        cc_i = ClassCountInd(1, SpinOrbSymLabel(so_elec), G1(so_elec)%ml)

        ! and the number of orbitals
        nOrb = OrbClassCount(cc_i)

        if (nOrb == 0) return

        ! get the symmetry index for later use
        sym_index = SymLabelCounts2(1, cc_i)

        ! now keep drawing from the symmetry orbitals until we pick an
        ! 'empty' (not doubly occupied!) one
        attempts = 0
        do while (attempts < 250)

            orb = 1 + floor(genrand_real2_dSFMT() * nOrb)

            s_orb = sym_label_list_spat(sym_index + orb - 1)

            ! if the spatial orbital is not doubly occupied an GUGA
            ! excitation (without taking into accound any other restrictions)
            ! can be possible
            ! and it must not be the original electron spatial index!
            if (current_stepvector(s_orb) /= 3 .and. s_orb /= s_elec) exit

            attempts = attempts + 1
#ifdef DEBUG_
            if (attempts > 200) then
                print *, "closing to 250 limit in random single orbital picking"
            end if
#endif
        end do

        if (s_orb < s_elec) then
            excitInfo = assign_excitinfo_values_single(gen_type%R, s_orb, s_elec, &
                s_orb, s_elec)
        else
            excitInfo = assign_excitinfo_values_single(gen_typ%L, s_orb, s_elec, &
                s_elec, s_orb)
        end if

        ! do i need nOrb now or the actual number of unoccupied in the
        ! CSF? i think the second..
        pgen = 1.0_dp / real(nOrb * nel, dp)

        if (current_stepvector(s_elec) == 3) pgen = 2.0_dp * pgen

    end subroutine pick_uniform_spatial_hole


end module guga_pchb_excitgen

