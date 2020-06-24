#include "macros.h"
#:include "macros.fpph"

#:set ExcitationTypes = ['SingleExc_t', 'DoubleExc_t']

module disconnected_gasci
    use SystemData, only: tGAS, tGASSpinRecoupling, nBasis, nel
    use constants
    use util_mod, only: get_free_unit, binary_search_first_ge, operator(.div.), &
        near_zero
    use sort_mod, only : sort
    use sets_mod, only: operator(.in.)
    use bit_rep_data, only: NIfTot, NIfD
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: pDoubles, pSingles
    use get_excit, only: make_double, make_single
    use Determinants, only: get_helement
    use excit_gens_int_weighted, only: pick_biased_elecs, pgen_select_orb
    use excitation_types, only: Excitation_t, SingleExc_t, DoubleExc_t, &
        get_last_tgt, set_last_tgt, defined, dyn_defined, UNKNOWN
    use sltcnd_mod, only: sltcnd_excit, dyn_sltcnd_excit
    use orb_idx_mod, only: SpinOrbIdx_t, calc_spin_raw, SpinProj_t, operator(==), operator(/=)
    use gasci, only: GASSpec_t
    implicit none

    private
    public :: init_disconnected_GAS, generate_nGAS_excitation, &
        clearGAS, oddBits, calc_pgen, dyn_calc_pgen


    #:for function_name in ['get_cumulative_list', 'get_mat_element', 'pick_weighted_hole', 'calc_pgen']
    interface ${function_name}$
        #:for Excitation_t in ExcitationTypes
            module procedure ${function_name}$_${Excitation_t}$
        #:endfor
    end interface
    #:endfor

    integer :: nGAS
    integer, parameter :: idx_alpha = 1, idx_beta = 2
#ifdef INT64_
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

    subroutine init_GAS(iGAS_per_spatorb)
        integer, intent(in) :: iGAS_per_spatorb(:)

        integer :: nOrbs, iOrb, iGAS

        associate(GAS => iGAS_per_spatorb)

        nOrbs = size(GAS)
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

        end associate

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
    end subroutine

    subroutine init_from_file()
        integer :: GAS_unit, nOrbs, GAS(nBasis .div. 2)

        nOrbs = nBasis .div. 2
        GAS_unit = get_free_unit()
        open(GAS_unit, file="GASOrbs", status='old')
            read(GAS_unit,  * ) GAS(1:nOrbs)
        close(GAS_unit)

        call init_GAS(GAS)
    end subroutine init_from_file

    subroutine init_disconnected_GAS(GAS_spec)
        type(GASSpec_t), intent(in) :: GAS_spec

        call init_GAS(GAS_spec%GAS_table(::2) )
    end subroutine

!----------------------------------------------------------------------------!

    subroutine clearGAS()
        if (allocated(GAS_size)) deallocate (GAS_size)
        if (allocated(GAS_spin_orb_list)) deallocate (GAS_spin_orb_list)
        if (allocated(GAS_table)) deallocate (GAS_table)
        if (allocated(GAS_spin_orbs)) deallocate (GAS_spin_orbs)
        if (allocated(GAS_orbs)) deallocate (GAS_orbs)
    end subroutine clearGAS

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
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type

        real(dp) :: r

        @:unused_var(exFlag, part_type, store)
#ifdef WARNING_WORKAROUND_
        hel = 0.0_dp
#endif

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
        integer :: ms, nJBase(nel)
        real(dp) :: r, pgen_pick1, pgen_pick2
        logical :: tExchange
        character(*), parameter :: this_routine = 'generate_nGAS_double'
        ! assuming that there are possible excitations within each active space,
        ! pick two random electrons (we would not include a full / empty space, so
        ! the assumption is very mild)
        call pick_biased_elecs(nI, elecs, src, sym_product, ispn, sum_ml, pgen)

        ! active spaces we consider
        srcGAS = GAS_table(src)

        @:ASSERT(tGASSpinRecoupling)

        tExchange = (ispn == 2) &
            .and. (tGASSpinRecoupling .or. srcGAS(1) == srcGAS(2))

        ! pick an empty orb from each active space chosen
        ! number of empty orbs in this active space
        r = genrand_real2_dSFMT()
        ! get the spin of the target orbital
        if (tExchange) then
            if (r > 0.5_dp) then
                ms = 1
                r = 2 * (r - 0.5_dp)
            else
                ms = 2
                r = 2 * (0.5_dp - r)
            end if
            pgen = pgen * 0.5_dp
        else
            ms = get_spin(src(1))
        end if
        tgt = 0
        ! the first hole is chosen randomly from the first active space
        tgt(1) = pick_hole_from_active_space(ilutI, nI, srcGAS(1), ms, genrand_real2_dSFMT(), pgen_pick1)
        if (tgt(1) == 0) then
            call zeroResult()
            return
        end if

        if (tExchange) then
            ! we picked the first spin randomly, now the second elec has the opposite spin
            ms = 3 - ms
        else
            ms = get_spin(src(2))
        end if
        ! the second hole is chosen in a weighted fashion
        tgt(2) = pick_weighted_hole(nI, DoubleExc_t(src(1), tgt(1), src(2)), ms, srcGAS(2), pgen_pick1)
        if (any(tgt == 0) .or. tgt(1) == tgt(2)) then
            call zeroResult()
            return
        end if

        ! if both excits are in the same GAS, we could also
        ! have picked them the other way around
        if (srcGAS(1) == srcGAS(2)) then
            pgen_pick2 = get_pgen_pick_weighted_hole(&
                                nI, DoubleExc_t(src(1), tgt(2), src(2), tgt(1))) &
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

    function get_pgen_pick_weighted_hole(nI, exc) result(pgenVal)
        integer, intent(in) :: nI(nel)
        type(DoubleExc_t), intent(in) :: exc
        character(*), parameter :: this_routine = 'get_pgen_pick_weighted_hole'

        real(dp) :: pgenVal

        @:ASSERT(defined(exc))

        associate (src1 => exc%val(1, 1), tgt1 => exc%val(2, 1), &
                   src2 => exc%val(1, 2), tgt2 => exc%val(2, 2))
        block
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
        end block
        end associate
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

    #:for Excitation_t in ExcitationTypes
        function get_cumulative_list_${Excitation_t}$(GAS_list, nI, incomplete_exc) result(cSum)
            integer, intent(in) :: GAS_list(:), nI(nel)
            type(${Excitation_t}$), intent(in) :: incomplete_exc
            real(dp) :: cSum(size(GAS_list))
            character(*), parameter :: this_routine = 'get_cumulative_list_${Excitation_t}$'

            real(dp) :: previous
            type(${Excitation_t}$) :: exc
            integer :: i, nOrbs

            nOrbs = size(GAS_list)

            exc = incomplete_exc
            ASSERT(get_last_tgt(exc) == UNKNOWN)

            ! build the cumulative list of matrix elements <src|H|tgt>
            previous = 0.0_dp
            do i = 1, nOrbs
                call set_last_tgt(exc, GAS_list(i))
                cSum(i) = get_mat_element(nI, exc) + previous
                previous = cSum(i)
            end do

            ! Normalize
            if (near_zero(cSum(nOrbs))) then
                cSum(:) = 0.0_dp
            else
                cSum(:) = cSum(:) / cSum(nOrbs)
            end if
        end function get_cumulative_list_${Excitation_t}$
    #:endfor

    function get_mat_element_SingleExc_t(nI, exc) result(res)
        integer, intent(in) :: nI(nEl)
        type(SingleExc_t), intent(in) :: exc
        real(dp) :: res

        associate (tgt => exc%val(2))
            if (all(nI /= tgt)) then
                res = abs(sltcnd_excit(nI, exc, .false.))
            else
                res = 0.0_dp
            end if
        end associate
    end function get_mat_element_SingleExc_t

    function get_mat_element_DoubleExc_t(nI, exc) result(res)
        integer, intent(in) :: nI(nEl)
        type(DoubleExc_t), intent(in) :: exc
        real(dp) :: res

        associate (src1 => exc%val(1, 1), tgt1 => exc%val(2, 1), &
                   src2 => exc%val(1, 2), tgt2 => exc%val(2, 2))
            if (tgt1 /= tgt2 .and. all(tgt2 /= nI)) then
                res = abs(sltcnd_excit(nI, exc, .false.))
            else
                res = 0.0_dp
            end if
        end associate

    end function get_mat_element_DoubleExc_t

!----------------------------------------------------------------------------!
    #:for Excitation_t in ExcitationTypes
        function pick_weighted_hole_${Excitation_t}$(nI, exc, spin_idx, iGAS, pgen) result(tgt)
            ! pick a hole of nI with spin ms from the active space with index
            ! srcGASInd the random number is to be supplied as r
            ! nI is the source determinant, nJBase the one from which we obtain
            ! the ket of the matrix element by single excitation
            integer, intent(in) :: nI(nel)
            type(${Excitation_t}$), intent(in) :: exc
            integer, intent(in) :: spin_idx, iGAS
            real(dp), intent(inout) :: pgen
            character(*), parameter :: this_routine = 'pick_weighted_hole_${Excitation_t}$'

            integer :: tgt, nOrbs, GAS_list(GAS_size(iGAS))
            real(dp) :: cSum(GAS_size(iGAS))


            ! initialize auxiliary variables
            nOrbs = GAS_size(iGAS)
            GAS_list = GAS_spin_orb_list(1:nOrbs, iGAS, spin_idx)
            ! build the cumulative list of matrix elements <src|H|tgt>
            ASSERT(get_last_tgt(exc) == UNKNOWN)
            cSum = get_cumulative_list(GAS_list, nI, exc)

            ! there might not be such an excitation
            if (cSum(nOrbs) > 0) then
                ! find the index of the target orbital in the gasList
                tgt = binary_search_first_ge(cSum, genrand_real2_dSFMT())

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

        end function pick_weighted_hole_${Excitation_t}$
    #:endfor

!----------------------------------------------------------------------------!


    function pick_hole_from_active_space(ilutI, nI, iGAS, ms, r, pgen) result(tgt)
        ! pick a hole of ilutI with spin ms from the active space with index srcGASInd
        ! the random number is to be supplied as r
        integer(n_int), intent(in) :: ilutI(0:NIfD)
        integer, intent(in) :: nI(nel)
        integer, intent(in) :: ms, iGAS
        real(dp), intent(in) :: r
        real(dp), intent(out) :: pgen
        integer :: tgt
        integer :: nEmpty, nOrb

        ! this sum only converts an array of size 1 to a scalar
        nEmpty = GAS_size(iGAS) - sum(popCnt(iand(ilutI(0:NIfD), GAS_spin_orbs(0:NIfD, iGAS, ms))))

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

        tgt = GAS_spin_orb_list(map_to_unoccupied(nOrb, iGAS, ms), iGAS, ms)

    contains

        !> Return the i-th unoccupied orbital in the iGAS space with spin ms.
        function map_to_unoccupied(i, iGAS, ms) result(new_orb_idx)
            integer, intent(in) :: i, iGAS, ms
            integer :: new_orb_idx

            integer :: j

            new_orb_idx = i
            do j = 1, nEl
                if (nI(j) > GAS_spin_orb_list(new_orb_idx, iGAS, ms)) return
                if (GAS_table(nI(j)) == iGAS .and. get_spin(nI(j)) == ms) then
                    ! if yes, we skip this orbital, increase by 1
                    new_orb_idx = new_orb_idx + 1
                end if
            end do
        end function
    end function pick_hole_from_active_space

    function calc_pgen_SingleExc_t(det_I, ilutI, exc) result(pgen)
        use FciMCData, only: pSingles
        type(SpinOrbIdx_t), intent(in) :: det_I
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        type(SingleExc_t), intent(in) :: exc
        character(*), parameter :: this_routine = 'get_pgen_pick_weighted_hole'

        real(dp) :: pgen

        real(dp) :: pgen_pick_elec, pgen_pick_weighted_hole

        @:ASSERT(defined(exc))
        @:unused_var(ilutI)

        pgen_pick_elec = 1.0_dp / nel

        associate(src1 => exc%val(1), tgt1 => exc%val(2))
        block
            real(dp) :: cSum(GAS_size(GAS_table(tgt1)))
            integer :: gasList(GAS_size(GAS_table(tgt1))), nOrbs
            integer :: j

            nOrbs = GAS_size(GAS_table(tgt1))
            gasList = GAS_spin_orb_list(1:nOrbs, GAS_table(tgt1), get_spin(tgt1))
            ! We know, that gasList contains tgt1
            j = binary_search_first_ge(gasList, tgt1)

            cSum = get_cumulative_list(gasList, det_I%idx, SingleExc_t(src1))

            if (j == 1) then
                pgen_pick_weighted_hole = cSum(1)
            else
                pgen_pick_weighted_hole = (cSum(j) - cSum(j - 1))
            end if
        end block
        end associate

        pgen = pSingles * pgen_pick_elec * pgen_pick_weighted_hole
    end function

    function calc_pgen_DoubleExc_t(det_I, ilutI, exc) result(pgen)
        use FciMCData, only: pDoubles, pParallel
        use SystemData, only: par_elec_pairs, AB_elec_pairs
        type(SpinOrbIdx_t), intent(in) :: det_I
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        type(DoubleExc_t), intent(in) :: exc
        character(*), parameter :: this_routine = 'calc_pgen_DoubleExc_t'

        real(dp) :: pgen
        real(dp) :: pgen_pick_elec, &
            ! These are arrays, because the pgen might be different if we
            ! pick AB or BA.
            ! pgen_first_pick == [p(A), p(B)]
            ! pgen_second_pick == [p(B | A), p(A | B)]
            pgen_first_pick(2), pgen_second_pick(2)

        @:ASSERT(defined(exc))

        associate (src1 => exc%val(1, 1), tgt1 => exc%val(2, 1), &
                   src2 => exc%val(1, 2), tgt2 => exc%val(2, 2))
        block
            real(dp) :: cSum(GAS_size(GAS_table(tgt1)))
            integer :: gasList(GAS_size(GAS_table(tgt1))), nOrbs
            integer :: j, iGAS, nEmpty, iOrb
            logical :: parallel_spin, tExchange

            parallel_spin = calc_spin_raw(src1) == calc_spin_raw(src2)
            if (parallel_spin) then
                pgen_pick_elec = pParallel / real(par_elec_pairs, dp)
            else
                pgen_pick_elec = (1.0_dp - pParallel) / real(AB_elec_pairs, dp)
            end if

            tExchange = .not. parallel_spin &
                        .and. (tGASSpinRecoupling .or. GAS_table(src1) == GAS_table(src2))

            iGAS = GAS_table(tgt1)
            nEmpty = GAS_size(iGAS) &
                     - sum(popCnt(iand(ilutI(0:NIfD), GAS_spin_orbs(0:NIfD, iGAS, get_spin(tgt1)))))
            if (nEmpty == 0) then
                pgen_first_pick(1) = 1.0_dp
            else
                ! adjust pgen
                pgen_first_pick(1) = 1.0_dp / real(nEmpty, dp)
            end if

            pgen_second_pick(1) = get_pgen_pick_weighted_hole(&
                    det_I%idx, DoubleExc_t(src1, tgt1, src2, tgt2))
            if (GAS_table(src1) == GAS_table(src2)) then
                pgen_first_pick(2) = get_pgen_pick_hole_from_active_space(&
                                    ilutI, GAS_table(src2), get_spin(tgt2))
                pgen_second_pick(2) = get_pgen_pick_weighted_hole(&
                    det_I%idx, DoubleExc_t(src1, tgt2, src2, tgt1))
            else
                pgen_first_pick(2) = 0.0_dp
                pgen_second_pick(2) = 0.0_dp
            end if
            pgen_first_pick = merge(0.5_dp, 1.0_dp, tExchange) * pgen_first_pick
        end block
        end associate

        pgen = pDoubles * pgen_pick_elec * sum(pgen_first_pick * pgen_second_pick)
    end function

    function dyn_calc_pgen(det_I, ilutI, exc) result(pgen)
        type(SpinOrbIdx_t), intent(in) :: det_I
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        class(Excitation_t), intent(in) :: exc
        character(*), parameter :: this_routine = 'dyn_calc_pgen'

        real(dp) :: pgen

        @:ASSERT(dyn_defined(exc))

        select type(exc)
        type is(SingleExc_t)
            pgen = calc_pgen(det_I, ilutI, exc)
        type is(DoubleExc_t)
            pgen = calc_pgen(det_I, ilutI, exc)
        end select
    end function

end module disconnected_gasci
