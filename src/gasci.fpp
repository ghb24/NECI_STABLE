#include "macros.h"
#:include "macros.fpph"
#:set ExcitationTypes = ['SingleExc_t', 'DoubleExc_t']

module gasci
    use SystemData, only: tNConservingGAS, tSpinConservingGAS, nBasis, nel
    use constants
    use util_mod, only: get_free_unit, binary_search_first_ge, &
        operator(.div.)
    use sort_mod, only : sort
    use bit_rep_data, only: NIfTot, NIfD
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: pDoubles
    use get_excit, only: make_double, make_single
    use Determinants, only: get_helement
    use excit_gens_int_weighted, only: pick_biased_elecs, pgen_select_orb
    use excitation_types, only: SingleExc_t, DoubleExc_t, UNKNOWN
    use sltcnd_mod, only: sltcnd_excit
    implicit none

    private
    public :: isValidExcit, loadGAS, generate_nGAS_excitation, clearGAS, &
        tSpinConservingGAS

    interface get_cumulative_list
        module procedure s_get_cumulative_list, d_get_cumulative_list
    end interface

    interface pick_weighted_hole
    #:for excitation_t in ExcitationTypes
        module procedure pick_weighted_hole_${excitation_t}$
    #:endfor
    end interface

    integer :: nGAS

    integer, parameter :: idx_alpha = 1, idx_beta = 2
#ifdef __INT64
    ! oddBits is the 2-base number: 101010101010....
    ! This corresponds to beta orbitals if interpreted as bitmask
    integer(n_int), parameter :: oddBits = 6148914691236517205_n_int
#else
    integer(n_int), parameter :: oddBits = 1431655765_n_int
#endif

    !> Bitmasks containing the active spaces
    !> (stored in the same format as an ilut)
    !> also have spin-resolved bitmasks
    integer(n_int), allocatable :: &
        !> The index of GAS_orbs is [0:NIfD, 1:nGAS]
        GAS_orbs(:, :), &
        !> Spin resolved version of GAS_orbs.
        !> The index is [0:NIfD, 1:nGAS, {idx_alpha, idx_beta}]
        GAS_spin_orbs(:, :, :)
    integer, allocatable :: &
        !> Integer list of spin-orbitals for each GAS space,
        !> that contains in continous order the occupied orbitals.
        !> The index is [i-th occupied orbital in given GAS space,
        !>               1:nGAS, {idx_alpha, idx_beta}].
        !> Note that a meaningful index is only [1:GAS_size(iGAS), iGAS, :]
        GAS_spin_orb_list(:, :, :), &
        !> Number of orbitals in each GAS space.
        GAS_size(:), &
        !> Lookup table containing the GAS space for each orbital.
        GAS_table(:)

contains

    subroutine loadGAS()
        integer :: GAS_unit, iOrb, nOrbs, GAS(nBasis .div. 2), iGAS

        nOrbs = nBasis .div. 2
        GAS_unit = get_free_unit()
        open(GAS_unit, file="GASOrbs", status='old')
            read(GAS_unit,  * ) GAS(1:nOrbs)
        close(GAS_unit)
        nGAS = maxval(GAS)
        allocate(GAS_orbs(0:NIfD, nGAS))
        GAS_orbs(:, :) = 0_n_int
        do iOrb = 1, nOrbs
            ! now write the orbitals read in to the current GAS
            ! set both the alpha- and the beta-spin orbital
            call setorb(GAS_orbs(:, GAS(iOrb)), 2 * iOrb)
            call setorb(GAS_orbs(:, GAS(iOrb)), 2 * iOrb - 1)
        end do

        ! now set up the auxiliary gas lookup tables
        allocate(GAS_spin_orbs(0:NIfD, nGAS, 2))
        ! GAS_spin_orbs is the same as GAS_orbs, but spin resolved.
        ! The order is beta, alpha, beta, alpha ...
        GAS_spin_orbs(:, :, idx_alpha) = iand(GAS_orbs, ishft(oddBits, 1))
        GAS_spin_orbs(:, :, idx_beta) = iand(GAS_orbs, oddBits)

        allocate(GAS_table(nBasis))
        ! gasTable contains the active space index for each spin orbital
        ! it is the same as GAS for spin orbitals
        GAS_table(1::2) = GAS(:)
        GAS_table(2::2) = GAS(:)

        allocate(GAS_spin_orb_list(nBasis, nGAS, 2))
        allocate(GAS_size(nGAS), source=0)
        associate (M => GAS_spin_orb_list)
            do iOrb = 1, nOrbs
                GAS_size(GAS(iOrb)) = GAS_size(GAS(iOrb)) + 1
                M(GAS_size(GAS(iOrb)), GAS(iOrb), idx_alpha) = 2 * iOrb
                M(GAS_size(GAS(iOrb)), GAS(iOrb), idx_beta) = 2 * iOrb - 1
            end do
        end associate

        do iGAS = 1, nGAS
            write (iout,  * ) "Number of orbs in GAS", iGAS, "is", &
                & sum(popCnt(GAS_orbs(:, iGAS)))
        end do

    contains

        ! One cannot use the macro here because of Fortran syntax rules.
        subroutine setOrb(ilut, orb)
            integer(n_int), intent(inout) :: ilut(0:NIfD)
            integer, intent(in) :: orb

            integer :: pos
            ! Get the bit index of the integer containing this orbital
            pos = (orb - 1) .div. bits_n_int
            ilut(pos) = ibset(ilut(pos), mod(orb - 1, bits_n_int))
        end subroutine setOrb
    end subroutine loadGAS

!----------------------------------------------------------------------------!

    subroutine clearGAS()
        if (allocated(GAS_size)) deallocate (GAS_size)
        if (allocated(GAS_spin_orb_list)) deallocate (GAS_spin_orb_list)
        if (allocated(GAS_table)) deallocate (GAS_table)
        if (allocated(GAS_spin_orbs)) deallocate (GAS_spin_orbs)
        if (allocated(GAS_orbs)) deallocate (GAS_orbs)
    end subroutine clearGAS

!----------------------------------------------------------------------------!

    pure function isValidExcit(ilut_i, ilut_j) result(valid)
        ! check if the excitation from ilut_i to ilut_j is valid within the GAS
        integer(n_int), intent(in) :: ilut_i(0:NIfTot), ilut_j(0:NIfTot)
        integer(n_int) :: GAS_ilut_i(0:NIfD), GAS_ilut_j(0:NIfD)
        logical :: valid

        integer :: iGAS

        valid = .true.
        ! safety check: do the GAS_orbs exist
        if (.not. allocated(GAS_orbs)) return
        do iGAS = 1, nGAS
            GAS_ilut_i = GAS_Component(ilut_i, iGAS)
            GAS_ilut_j = GAS_Component(ilut_j, iGAS)
            ! check if ilut_i and ilut_j have the same number of electrons in GAS-i
            valid = sum(popCnt(GAS_ilut_i)) == sum(popCnt(GAS_ilut_j))
            if (tSpinConservingGAS) then
                ! check if the number of odd bits set
                ! (=number of beta electrons) is the
                ! same in both determinants.
                valid = sum(popCnt(iand(GAS_ilut_i, oddBits))) &
                        &== sum(popCnt(iand(GAS_ilut_j, oddBits)))
            end if
        end do

    contains

        pure function GAS_Component(ilut, i) result(GAS_Ilut)
            integer(n_int), intent(in) :: ilut(0:NIfTot)
            integer, intent(in) :: i

            integer(n_int) :: GAS_Ilut(0:NIfD)

            GAS_Ilut = iand(ilut(0:NIfD), GAS_orbs(0:NIfD, i))

        end function GAS_Component

    end function isValidExcit

!----------------------------------------------------------------------------!

    subroutine generate_nGAS_excitation(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                        ex, tParity, pGen, hel, store, part_type)
        ! particle number conserving GAS excitation generator:
        ! we create only excitations, that preserver the number of electrons within
        ! each active space
        use SystemData, only: nel
        use FciMCData, only: excit_gen_store_type
        use constants

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, 2)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        real(dp) :: r

        @:unused_var(exFlag, store, part_type)

        ! exFlag and part_type are not used but part of interface for
        ! procedure pointers.
        hel = 0.0_dp

        ! single or double excitation?
        r = genrand_real2_dSFMT()
        if (r < pDoubles) then
            call generate_nGAS_double(nI, ilutI, nJ, ilutJ, ex, tParity, pgen)
            pgen = pgen * pDoubles
            ic = 2
        else
            call generate_nGAS_single(nI, ilutI, nJ, ilutJ, ex, tParity, pgen)
            pgen = pgen * (1.0_dp - pDoubles)
            ic = 1
        end if

    end subroutine generate_nGAS_excitation

!----------------------------------------------------------------------------!

    subroutine generate_nGAS_single(nI, ilutI, nJ, ilutJ, ex, par, pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2, 2)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen

        integer :: elec, tgt, src, spin_idx
        real(dp) :: r

        ! we assume that each active space has a possible single excitation
        ! (this is very mild because active spaces without such are trivial)
        ! -> pick any random electron
        r = genrand_real2_dSFMT()
        elec = int(r * nel) + 1
        src = nI(elec)
        ! adjust pgen
        pgen = 1.0_dp / nel
        ! from the same active space, get a hole
        spin_idx = get_spin(src)
        tgt = pick_weighted_hole( &
            nI, SingleExc_t(src), spin_idx, GAS_table(src), pgen)

        if (tgt == 0) then
            nJ(1) = 0
            ilutJ = 0_n_int
            return
        end if

        call make_single(nI, nJ, elec, tgt, ex, par)
        ilutJ = ilutI
        clr_orb(ilutJ, src)
        set_orb(ilutJ, tgt)

    end subroutine generate_nGAS_single

!----------------------------------------------------------------------------!

    subroutine generate_nGAS_double(nI, ilutI, nJ, ilutJ, ex, par, pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2, 2)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen

        integer :: elecs(2), src(2), sym_product, ispn, sum_ml, tgt(2)
        integer :: srcGAS(2)
        integer :: spin_idx, nJBase(nel)
        real(dp) :: r, pgen_pick1, pgen_pick2
        logical :: tExchange
        ! assuming that there are possible excitations within each active space,
        ! pick two random electrons (we would not include a full / empty space, so
        ! the assumption is very mild)
        call pick_biased_elecs(nI, elecs, src, sym_product, ispn, sum_ml, pgen)

        ! active spaces we consider
        srcGAS = GAS_table(src)

        tExchange = (ispn == 2) &
            .and. (.not. tSpinConservingGAS .or. srcGAS(1) == srcGAS(2))

        ! pick an empty orb from each active space chosen
        ! number of empty orbs in this active space
        r = genrand_real2_dSFMT()
        ! get the spin of the target orbital
        if (tExchange) then
            if (r > 0.5_dp) then
                spin_idx = 1
                r = 2 * (r - 0.5_dp)
            else
                spin_idx = 2
                r = 2 * (0.5_dp - r)
            end if
            pgen = pgen * 0.5_dp
        else
            spin_idx = get_spin(src(1))
        end if
        tgt = 0
        ! the first hole is chosen randomly from the first active space
        tgt(1) = pick_hole_from_active_space(&
                ilutI, nI, srcGAS(1), spin_idx, r, pgen_pick1)
        if (tgt(1) == 0) then
            call zeroResult()
            return
        end if

        if (tExchange) then
            ! we picked the first spin randomly,
            ! now the second elec has the opposite spin
            spin_idx = 3 - spin_idx
        else
            spin_idx = get_spin(src(2))
        end if
        ! the second hole is chosen in a weighted fashion
        tgt(2) = pick_weighted_hole(&
            nI, DoubleExc_t(src(1), src(2), tgt(1)), spin_idx, srcGAS(2), pgen_pick1)
        if (any(tgt == 0) .or. tgt(1) == tgt(2)) then
            call zeroResult()
            return
        end if

        ! if both excits are in the same GAS, we could also
        ! have picked them the other way around
        if (srcGAS(1) == srcGAS(2)) then
            nJBase = nI
            nJBase(elecs(1)) = tgt(2)
            pgen_pick2 = get_pgen_pick_weighted_hole(&
                                nI, src(1), src(2), tgt(2), tgt(1)) &
                         * get_pgen_pick_hole_from_active_space(&
                                ilutI, srcGAS(2), get_spin(tgt(2)))
            pgen = pgen * (pgen_pick1 + pgen_pick2)
        else
            pgen = pgen * pgen_pick1
        end if

        call make_double(nI, nJ, elecs(1), elecs(2), tgt(1), tgt(2), ex, par)
        ilutJ = ilutI
        clr_orb(ilutJ, src(1))
        clr_orb(ilutJ, src(2))
        set_orb(ilutJ, tgt(1))
        set_orb(ilutJ, tgt(2))

    contains

        subroutine zeroResult()
            integer :: src_copy(2)

            pgen = pgen * pgen_pick1
            src_copy(:) = src(:)
            call sort(src_copy)
            ex(1, :) = src_copy
            ex(2, :) = tgt
            nJ(1) = 0
            ilutJ = 0_n_int
        end subroutine zeroResult

    end subroutine generate_nGAS_double

!----------------------------------------------------------------------------!

    function get_pgen_pick_weighted_hole(nI, src1, src2, tgt1, tgt2) result(pgenVal)
        integer, intent(in) :: nI(nel)
        integer, intent(in) :: src1, src2, tgt1, tgt2

        real(dp) :: pgenVal
        real(dp) :: cSum(GAS_size(GAS_table(tgt2)))
        integer :: gasList(GAS_size(GAS_table(tgt2))), nOrbs
        integer :: gasInd2

        nOrbs = GAS_size(GAS_table(tgt2))
        gasList = GAS_spin_orb_list(1:nOrbs, GAS_table(tgt2), get_spin(tgt2))

        cSum = get_cumulative_list(gasList, nI, DoubleExc_t(src1, tgt1, src2))
        ! we know gasList contains tgt2, so we can look up its index with binary search
        gasInd2 = binary_search_first_ge(gasList, tgt2)
        if (gasInd2 == 1) then
            pgenVal = cSum(1)
        else
            pgenVal = (cSum(gasInd2) - cSum(gasInd2 - 1))
        end if

    end function get_pgen_pick_weighted_hole

!----------------------------------------------------------------------------!

    function get_pgen_pick_hole_from_active_space(ilut, srcGASInd, spin_idx) result(pgenVal)
        integer(n_int), intent(in) :: ilut(0:NIfD)
        integer, intent(in) :: srcGASInd, spin_idx

        real(dp) :: pgenVal
        integer :: nEmpty
        nEmpty = GAS_size(srcGASInd) &
                 - sum(popCnt(iand(ilut(0:NIfD), &
                                   GAS_spin_orbs(0:NIfD, srcGASInd, spin_idx))))
        ! if the excitation is not possible, pgen is void
        if (nEmpty == 0) then
            pgenVal = 1.0_dp
        else
            ! adjust pgen
            pgenVal = 1.0_dp / real(nEmpty, dp)
        end if
    end function get_pgen_pick_hole_from_active_space

!----------------------------------------------------------------------------!

    function s_get_cumulative_list(&
            Orb_list, occ_Orb_list, exc) result(cSum)
        integer, intent(in) :: Orb_list(:), occ_Orb_list(nel)
        type(SingleExc_t), intent(in) :: exc
        real(dp) :: cSum(size(Orb_list))

        type(SingleExc_t) :: trial_exc
        real(dp) :: previous
        integer :: i, tgt
        character(*), parameter :: this_routine = "s_get_cumulative_list"

        trial_exc = exc
        ASSERT(trial_exc%val(2) == UNKNOWN)

        previous = 0.0_dp
        do i = 1, size(Orb_list)
            tgt = Orb_list(i)
            if (all(occ_Orb_list /= tgt)) then
                trial_exc%val(2) = tgt
                cSum(i) = &
                    previous + abs(sltcnd_excit(occ_Orb_list, trial_exc, .false.))
            else
                cSum(i) = previous
            end if
            previous = cSum(i)
        end do

        cSum(:) = cSum(:) / cSum(size(cSum))
    end function

    function d_get_cumulative_list(&
            Orb_list, occ_Orb_list, exc) result(cSum)
        implicit none
        integer, intent(in) :: Orb_list(:), occ_Orb_list(nel)
        type(DoubleExc_t), intent(in) :: exc
        real(dp) :: cSum(size(Orb_list))

        type(DoubleExc_t) :: trial_exc
        real(dp) :: previous
        integer :: i, tgt
        character(*), parameter :: this_routine = "d_get_cumulative_list"

        trial_exc = exc
        ASSERT(trial_exc%val(2, 2) == UNKNOWN)

        previous = 0.0_dp
        do i = 1, size(Orb_list)
            tgt = Orb_list(i)
            if (all(occ_Orb_list /= tgt) .and. trial_exc%val(1, 2) /= tgt) then
                trial_exc%val(2, 2) = tgt
                cSum(i) = &
                    previous + abs(sltcnd_excit(occ_Orb_list, trial_exc, .false.))
            else
                cSum(i) = previous
            end if
            previous = cSum(i)
        end do

        cSum(:) = cSum(:) / cSum(size(cSum))
    end function


!----------------------------------------------------------------------------!

#:for excitation_t in ExcitationTypes
    function pick_weighted_hole_${excitation_t}$(&
            & nI, exc, spin_idx, iGAS, pgen) result(tgt)
        ! pick a hole of nI with spin ms from the active space with index
        ! srcGASInd the random number is to be supplied as r
        ! nI is the source determinant, nJBase the one from which we obtain
        ! the ket of the matrix element by single excitation
        integer, intent(in) :: nI(nel), spin_idx, iGAS
        type(${excitation_t}$), intent(in) :: exc
        real(dp), intent(inout) :: pgen

        integer :: tgt

        integer :: nOrbs, GAS_list(GAS_size(iGAS))
        real(dp) :: r, cSum(GAS_size(iGAS))

        ! initialize auxiliary variables
        nOrbs = GAS_size(iGAS)
        GAS_list = GAS_spin_orb_list(1:nOrbs, iGAS, spin_idx)
        ! build the cumulative list of matrix elements <src|H|tgt>
        cSum = get_cumulative_list(GAS_list, nI, exc)

        ! now, pick with the weight from the cumulative list
        r = genrand_real2_dSFMT()

        ! there might not be such an excitation
        if (cSum(nOrbs) > 0) then
            ! find the index of the target orbital in the gasList
            tgt = binary_search_first_ge(cSum, r)

            ! adjust pgen with the probability for picking tgt from the cumulative list
            if (tgt == 1) then
                pgen = pgen * cSum(1)
            else
                pgen = pgen * (cSum(tgt) - cSum(tgt - 1))
            end if

            ! convert to global orbital index
            tgt = GAS_list(tgt)
        else
            tgt = 0
        end if

    end function pick_weighted_hole_${excitation_t}$
#:endfor

!----------------------------------------------------------------------------!

    function pick_hole_from_active_space(ilutI, nI, srcGASInd, ms, r, pgen) result(tgt)
        ! pick a hole of ilutI with spin ms from the active space with index srcGASInd
        ! the random number is to be supplied as r
        integer(n_int), intent(in) :: ilutI(0:NIfD)
        integer, intent(in) :: nI(nel)
        integer, intent(in) :: ms, srcGASInd
        real(dp), intent(in) :: r
        real(dp), intent(out) :: pgen
        integer :: tgt
        integer :: nEmpty, nOrb

        ! this sum only converts an array of size 1 to a scalar
        nEmpty = GAS_size(srcGASInd) - sum(popCnt(iand(ilutI(0:NIfD), GAS_spin_orbs(0:NIfD, srcGASInd, ms))))

        ! if there are no empyty orbitals in this gas (can happen when allowing for
        ! spin-exchange), the excitation is invalid
        if (nEmpty == 0) then
            tgt = 0
            pGen = 1.0_dp
            return
        end if
        ! adjust pgen
        pgen = 1.0_dp / real(nEmpty, dp)
        ! TODO: check if there are enough empty orbs
        ! Re: We can safely assume there is always an empty orb in each active space

        ! index of the target orb
        nOrb = int(r * nEmpty) + 1
        call skipOrb(nOrb, GAS_spin_orbs(0:NIfD, srcGASInd, ms))
        nOrb = GAS_spin_orb_list(nOrb, srcGASInd, ms)
        tgt = nOrb

    contains

        subroutine skipOrb(nOrb, gasIlut)
            ! convert an unoccupied active space orbital index to an
            ! active space orbital index (i.e. correct for occ. orbs)
            integer, intent(inout) :: nOrb
            integer(n_int), intent(in) :: gasIlut(0:NIfD)
            integer :: j, globalOrbIndex

            do j = 1, nel
                globalOrbIndex = GAS_spin_orb_list(nOrb, srcGASInd, ms)
                if (nI(j) > globalOrbIndex) return
                ! check if an occ. orb is in the target active space
                if (validTarget(nI(j), gasIlut)) &
                    ! if yes, we skip this orbital, increase nOrb by 1
                    nOrb = nOrb + 1
            end do
        end subroutine skipOrb

        function validTarget(orb, gasIlut) result(valid)
            ! check if an orbital is a valid target for an excitation within an active space
            integer, intent(in) :: orb
            integer(n_int), intent(in) :: gasIlut(0:NIfD)

            logical :: valid

            ! is the orbital in the acitve space?
            valid = btest(gasIlut((orb - 1) / bits_n_int), mod(orb - 1, bits_n_int))
        end function validTarget
    end function pick_hole_from_active_space

end module gasci
