#include "macros.h"
#:include "macros.fpph"

#:set ExcitationTypes = ['SingleExc_t', 'DoubleExc_t']

module gasci_disconnected
    use SystemData, only: tGAS, nBasis, nel
    use constants
    use SymExcitDataMod, only: ScratchSize
    use util_mod, only: get_free_unit, binary_search_first_ge, operator(.div.), &
                        near_zero
    use sort_mod, only: sort
    use sets_mod, only: operator(.in.)
    use bit_rep_data, only: NIfTot, NIfD
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData, only: pDoubles, pSingles
    use get_excit, only: make_double, make_single
    use Determinants, only: get_helement
    use excit_gens_int_weighted, only: pick_biased_elecs, pgen_select_orb, get_pgen_pick_biased_elecs
    use excitation_types, only: Excitation_t, SingleExc_t, DoubleExc_t, &
                                get_last_tgt, set_last_tgt, defined, dyn_defined, UNKNOWN
    use sltcnd_mod, only: sltcnd_excit, dyn_sltcnd_excit
    use orb_idx_mod, only: SpinOrbIdx_t, calc_spin_raw, SpinProj_t, operator(==), operator(/=)
    use gasci, only: GASSpec_t
    use gasci_util, only: gen_all_excits

    use excitation_generators, only: ExcitationGenerator_t
    implicit none

    private
    public :: GAS_disc_ExcGenerator_t


    #:for function_name in ['get_cumulative_list', 'get_mat_element']
    interface ${function_name}$
        #:for Excitation_t in ExcitationTypes
        module procedure ${function_name}$_${Excitation_t}$
        #:endfor
    end interface
    #:endfor

    integer, parameter :: idx_alpha = 1, idx_beta = 2
#ifdef INT64_
    ! oddBits is the 2-base number: 101010101010....
    ! This corresponds to beta orbitals if interpreted as bitmask
    integer(n_int), parameter :: oddBits = 6148914691236517205_n_int
#else
    integer(n_int), parameter :: oddBits = 1431655765_n_int
#endif

    type, extends(ExcitationGenerator_t) :: GAS_disc_ExcGenerator_t
        private
        class(GASSpec_t), allocatable :: GAS_spec
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
        private
        procedure, public :: finalize => GAS_disc_finalize
        procedure, public :: gen_exc => GAS_disc_gen_exc
        procedure, public :: get_pgen => GAS_disc_get_pgen
        procedure, public :: gen_all_excits => GAS_disc_gen_all_excits

        procedure, public :: generate_nGAS_single
        procedure, public :: generate_nGAS_double

        procedure :: pick_hole_from_active_space
        procedure :: get_pgen_pick_weighted_hole
        procedure :: get_pgen_pick_hole_from_active_space

        #:for Excitation_t in ExcitationTypes
        procedure :: pick_weighted_hole_${Excitation_t}$
        generic :: pick_weighted_hole => pick_weighted_hole_${Excitation_t}$

        procedure :: calc_pgen_${Excitation_t}$
        generic :: calc_pgen => calc_pgen_${Excitation_t}$
        #:endfor
    end type

    interface GAS_disc_ExcGenerator_t
        module procedure construct_GAS_disc_ExcGenerator_t
    end interface


contains

    pure function construct_GAS_disc_ExcGenerator_t(GAS_spec) result(res)
        class(GASSpec_t), intent(in) :: GAS_spec
        type(GAS_disc_ExcGenerator_t) :: res

        !> GAS space for the i-th **spatial** orbital
        integer, allocatable :: spat_GAS(:)
        integer :: i, nOrbs, iOrb, nGAS

        res%GAS_spec = GAS_spec

        spat_GAS = GAS_spec%get_iGAS([(i, i = 1, GAS_spec%n_spin_orbs(), 2)])

        nOrbs = size(spat_GAS)
        nGAS = maxval(spat_GAS)

        allocate(res%GAS_orbs(0:NIfD, nGAS), source=0_n_int)
        do iOrb = 1, nOrbs
            ! now write the orbitals read in to the current GAS
            ! set both the alpha- and the beta-spin orbital
            call setorb(res%GAS_orbs(:, spat_GAS(iOrb)), 2 * iOrb)
            call setorb(res%GAS_orbs(:, spat_GAS(iOrb)), 2 * iOrb - 1)
        end do

        ! now set up the auxiliary gas lookup tables
        allocate(res%GAS_spin_orbs(0:NIfD, nGAS, 2))
        ! GAS_spin_orbs is the same as GAS_orbs, but spin resolved.
        ! The order is beta, alpha, beta, alpha ...
        res%GAS_spin_orbs(:, :, idx_alpha) = iand(res%GAS_orbs, ishft(oddBits, 1))
        res%GAS_spin_orbs(:, :, idx_beta) = iand(res%GAS_orbs, oddBits)

        allocate(res%GAS_table(nBasis))
        ! gasTable contains the active space index for each spin orbital
        ! it is the same as GAS for spin orbitals
        res%GAS_table(1::2) = spat_GAS(:)
        res%GAS_table(2::2) = spat_GAS(:)

        allocate(res%GAS_spin_orb_list(nBasis, nGAS, 2))
        allocate(res%GAS_size(nGAS), source=0)
        do iOrb = 1, nOrbs
            res%GAS_size(spat_GAS(iOrb)) = res%GAS_size(spat_GAS(iOrb)) + 1
            res%GAS_spin_orb_list(res%GAS_size(spat_GAS(iOrb)), spat_GAS(iOrb), idx_alpha) = 2 * iOrb
            res%GAS_spin_orb_list(res%GAS_size(spat_GAS(iOrb)), spat_GAS(iOrb), idx_beta) = 2 * iOrb - 1
        end do

    contains

        ! One cannot use the macro here because of Fortran syntax rules.
        pure subroutine setOrb(ilut, orb)
            integer(n_int), intent(inout) :: ilut(0:NIfD)
            integer, intent(in) :: orb

            integer :: pos
            ! Get the bit index of the integer containing this orbital
            pos = (orb - 1) .div. bits_n_int
            ilut(pos) = ibset(ilut(pos), mod(orb - 1, bits_n_int))
        end subroutine setOrb

    end function

    subroutine GAS_disc_finalize(this)
        class(GAS_disc_ExcGenerator_t), intent(inout) :: this
        if (allocated(this%GAS_size)) deallocate(this%GAS_size)
        if (allocated(this%GAS_spin_orb_list)) deallocate(this%GAS_spin_orb_list)
        if (allocated(this%GAS_table)) deallocate(this%GAS_table)
        if (allocated(this%GAS_spin_orbs)) deallocate(this%GAS_spin_orbs)
        if (allocated(this%GAS_orbs)) deallocate(this%GAS_orbs)
    end subroutine



    subroutine GAS_disc_gen_exc(this, nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                ex, tParity, pGen, hel, store, part_type)
        ! particle number conserving GAS excitation generator:
        ! we create only excitations, that preserver the number of electrons within
        ! each active space
        use SystemData, only: nel
        use FciMCData, only: excit_gen_store_type
        use constants
        class(GAS_disc_ExcGenerator_t), intent(inout) :: this
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
            call this%generate_nGAS_double(nI, ilutI, nJ, ilutJ, ex, tParity, pgen)
            pgen = pgen * pDoubles
            ic = 2
        else
            call this%generate_nGAS_single(nI, ilutI, nJ, ilutJ, ex, tParity, pgen)
            pgen = pgen * (1.0_dp - pDoubles)
            ic = 1
        end if
    end subroutine

    function GAS_disc_get_pgen(this, nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2) &
                result(pgen)
        class(GAS_disc_ExcGenerator_t), intent(inout) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit), ic
        integer, intent(in) :: ClassCount2(ScratchSize), ClassCountUnocc2(ScratchSize)
        real(dp) :: pgen
        character(*), parameter :: this_routine = 'GAS_disc_get_pgen'

        @:unused_var(ClassCount2, ClassCountUnocc2)

        select case(ic)
        case(1)
            pgen = pSingles * this%calc_pgen(SpinOrbIdx_t(nI), ilutI, SingleExc_t(ex(:, 1)))
        case(2)
            pgen = pDoubles * this%calc_pgen(SpinOrbIdx_t(nI), ilutI, DoubleExc_t(ex(:, :2)))
        case default
            call stop_all(this_routine, 'ic has to be 1 <= ic <= 2')
        end select
    end function


    subroutine GAS_disc_gen_all_excits(this, nI, n_excits, det_list)
        class(GAS_disc_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:)
        call gen_all_excits(this%GAS_spec, nI, n_excits, det_list)
    end subroutine


    subroutine generate_nGAS_single(this, nI, ilutI, nJ, ilutJ, ex, par, pgen)
        class(GAS_disc_ExcGenerator_t), intent(in) :: this
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
        tgt = this%pick_weighted_hole(nI, SingleExc_t(src), spin_idx, this%GAS_table(src), pgen)

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


    subroutine generate_nGAS_double(this, nI, ilutI, nJ, ilutJ, ex, par, pgen)
        class(GAS_disc_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ex(2, 2)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: par
        real(dp), intent(out) :: pgen

        integer :: elecs(2), src(2), sym_product, ispn, sum_ml, tgt(2)
        integer :: srcGAS(2)
        integer :: ms
        real(dp) :: r, pgen_pick1, pgen_pick2
        logical :: tExchange
        character(*), parameter :: this_routine = 'generate_nGAS_double'
        ! assuming that there are possible excitations within each active space,
        ! pick two random electrons (we would not include a full / empty space, so
        ! the assumption is very mild)
        call pick_biased_elecs(nI, elecs, src, sym_product, ispn, sum_ml, pgen)

        ! active spaces we consider
        srcGAS = this%GAS_table(src)

        tExchange = (ispn == 2) &
                    .and. (this%GAS_spec%recoupling() .or. srcGAS(1) == srcGAS(2))

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
        tgt(1) = this%pick_hole_from_active_space(ilutI, nI, srcGAS(1), ms, genrand_real2_dSFMT(), pgen_pick1)
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
        tgt(2) = this%pick_weighted_hole(nI, DoubleExc_t(src(1), tgt(1), src(2)), ms, srcGAS(2), pgen_pick1)
        if (any(tgt == 0) .or. tgt(1) == tgt(2)) then
            call zeroResult()
            return
        end if

        ! if both excits are in the same GAS, we could also
        ! have picked them the other way around
        if (srcGAS(1) == srcGAS(2)) then
            pgen_pick2 = &
                this%get_pgen_pick_weighted_hole(nI, DoubleExc_t(src(1), tgt(2), src(2), tgt(1))) &
                * this%get_pgen_pick_hole_from_active_space(ilutI, srcGAS(2), get_spin(tgt(2)))
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


    function get_pgen_pick_weighted_hole(this, nI, exc) result(pgenVal)
        class(GAS_disc_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        type(DoubleExc_t), intent(in) :: exc
        character(*), parameter :: this_routine = 'get_pgen_pick_weighted_hole'

        real(dp) :: pgenVal

        integer :: src1, tgt1, src2, tgt2
        real(dp) :: cSum(this%GAS_size(this%GAS_table(exc%val(2, 2))))
        integer :: gasList(this%GAS_size(this%GAS_table(exc%val(2, 2)))), nOrbs
        integer :: gasInd2

        @:ASSERT(defined(exc))
        src1 = exc%val(1, 1); tgt1 = exc%val(2, 1)
        src2 = exc%val(1, 2); tgt2 = exc%val(2, 2)

        nOrbs = this%GAS_size(this%GAS_table(tgt2))
        gasList = this%GAS_spin_orb_list(1:nOrbs, this%GAS_table(tgt2), get_spin(tgt2))

        cSum = get_cumulative_list(gasList, nI, DoubleExc_t(src1, tgt1, src2))
        ! we know gasList contains tgt2, so we can look up its index with binary search
        gasInd2 = binary_search_first_ge(gasList, tgt2)
        if (gasInd2 == 1) then
            pgenVal = cSum(1)
        else
            pgenVal = (cSum(gasInd2) - cSum(gasInd2 - 1))
        end if
    end function get_pgen_pick_weighted_hole


    function get_pgen_pick_hole_from_active_space(this, ilut, srcGASInd, spin_idx) result(pgenVal)
        class(GAS_disc_ExcGenerator_t), intent(in) :: this
        integer(n_int), intent(in) :: ilut(0:NIfD)
        integer, intent(in) :: srcGASInd, spin_idx

        real(dp) :: pgenVal
        integer :: nEmpty
        nEmpty = this%GAS_size(srcGASInd) &
                 - sum(popCnt(iand(ilut(0:NIfD), &
                                   this%GAS_spin_orbs(0:NIfD, srcGASInd, spin_idx) &
                                  ) &
                             ) &
                      )
        ! if the excitation is not possible, pgen is void
        if (nEmpty == 0) then
            pgenVal = 1.0_dp
        else
            ! adjust pgen
            pgenVal = 1.0_dp / real(nEmpty, dp)
        end if
    end function get_pgen_pick_hole_from_active_space


    #:for Excitation_t in ExcitationTypes
    function get_cumulative_list_${Excitation_t}$ (GAS_list, nI, incomplete_exc) result(cSum)
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

        if (all(nI /= exc%val(2))) then
            res = abs(sltcnd_excit(nI, exc, .false.))
        else
            res = 0.0_dp
        end if
    end function get_mat_element_SingleExc_t

    function get_mat_element_DoubleExc_t(nI, exc) result(res)
        integer, intent(in) :: nI(nEl)
        type(DoubleExc_t), intent(in) :: exc
        real(dp) :: res

        if (exc%val(2, 1) /= exc%val(2, 2) &
            .and. all(exc%val(2, 1) /= nI) .and. all(exc%val(2, 2) /= nI)) then
            res = abs(sltcnd_excit(nI, exc, .false.))
        else
            res = 0.0_dp
        end if

    end function get_mat_element_DoubleExc_t

    #:for Excitation_t in ExcitationTypes
    function pick_weighted_hole_${Excitation_t}$ (this, nI, exc, spin_idx, iGAS, pgen) result(tgt)
        ! pick a hole of nI with spin ms from the active space with index
        ! srcGASInd the random number is to be supplied as r
        ! nI is the source determinant, nJBase the one from which we obtain
        ! the ket of the matrix element by single excitation
        class(GAS_disc_ExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        type(${Excitation_t}$), intent(in) :: exc
        integer, intent(in) :: spin_idx, iGAS
        real(dp), intent(inout) :: pgen
        character(*), parameter :: this_routine = 'pick_weighted_hole_${Excitation_t}$'

        integer :: tgt, nOrbs, GAS_list(this%GAS_size(iGAS))
        real(dp) :: cSum(this%GAS_size(iGAS))

        ! initialize auxiliary variables
        nOrbs = this%GAS_size(iGAS)
        GAS_list = this%GAS_spin_orb_list(1:nOrbs, iGAS, spin_idx)
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


    function pick_hole_from_active_space(this, ilutI, nI, iGAS, ms, r, pgen) result(tgt)
        ! pick a hole of ilutI with spin ms from the active space with index srcGASInd
        ! the random number is to be supplied as r
        class(GAS_disc_ExcGenerator_t), intent(in) :: this
        integer(n_int), intent(in) :: ilutI(0:NIfD)
        integer, intent(in) :: nI(nel)
        integer, intent(in) :: ms, iGAS
        real(dp), intent(in) :: r
        real(dp), intent(out) :: pgen
        integer :: tgt
        integer :: nEmpty, nOrb

        ! this sum only converts an array of size 1 to a scalar
        nEmpty = this%GAS_size(iGAS) - sum(popCnt(iand(ilutI(0:NIfD), this%GAS_spin_orbs(0:NIfD, iGAS, ms))))

        ! if there are no empyty orbitals in this gas (can happen when allowing for
        ! spin-exchange), the excitation is invalid
        if (nEmpty == 0) then
            tgt = 0
            pGen = 1.0_dp
            return
        end if
        ! adjust pgen
        pgen = 1.0_dp / real(nEmpty, dp)
        ! We can safely assume there is always an empty orb in each active space

        ! index of the target orb
        nOrb = int(r * nEmpty) + 1

        tgt = this%GAS_spin_orb_list(map_to_unoccupied(nOrb, iGAS, ms), iGAS, ms)

    contains

        !> Return the i-th unoccupied orbital in the iGAS space with spin ms.
        function map_to_unoccupied(i, iGAS, ms) result(new_orb_idx)
            integer, intent(in) :: i, iGAS, ms
            integer :: new_orb_idx

            integer :: j

            new_orb_idx = i
            do j = 1, nEl
                if (nI(j) > this%GAS_spin_orb_list(new_orb_idx, iGAS, ms)) return
                if (this%GAS_table(nI(j)) == iGAS .and. get_spin(nI(j)) == ms) then
                    ! if yes, we skip this orbital, increase by 1
                    new_orb_idx = new_orb_idx + 1
                end if
            end do
        end function
    end function pick_hole_from_active_space

    function calc_pgen_SingleExc_t(this, det_I, ilutI, exc) result(pgen)
        class(GAS_disc_ExcGenerator_t), intent(in) :: this
        type(SpinOrbIdx_t), intent(in) :: det_I
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        type(SingleExc_t), intent(in) :: exc
        character(*), parameter :: this_routine = 'get_pgen_pick_weighted_hole'

        real(dp) :: pgen

        real(dp) :: pgen_pick_elec, pgen_pick_weighted_hole

        @:ASSERT(defined(exc))
        @:unused_var(ilutI)

        pgen_pick_elec = 1.0_dp / nel

        block
            real(dp) :: cSum(this%GAS_size(this%GAS_table(exc%val(2))))
            integer :: gasList(this%GAS_size(this%GAS_table(exc%val(2)))), nOrbs
            integer :: j

            nOrbs = this%GAS_size(this%GAS_table(exc%val(2)))
            gasList = this%GAS_spin_orb_list(1:nOrbs, this%GAS_table(exc%val(2)), get_spin(exc%val(2)))
            ! We know, that gasList contains the target
            j = binary_search_first_ge(gasList, exc%val(2))

            cSum = get_cumulative_list(gasList, det_I%idx, SingleExc_t(exc%val(1)))

            if (j == 1) then
                pgen_pick_weighted_hole = cSum(1)
            else
                pgen_pick_weighted_hole = (cSum(j) - cSum(j - 1))
            end if
        end block

        pgen = pgen_pick_elec * pgen_pick_weighted_hole
    end function

    function calc_pgen_DoubleExc_t(this, det_I, ilutI, exc) result(pgen)
        use FciMCData, only: pParallel
        use SystemData, only: par_elec_pairs, AB_elec_pairs
        class(GAS_disc_ExcGenerator_t), intent(in) :: this
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
        integer :: iGAS, nEmpty
        logical :: parallel_spin, tExchange
        integer :: src1, tgt1, src2, tgt2

        src1 = exc%val(1, 1); tgt1 = exc%val(2, 1)
        src2 = exc%val(1, 2); tgt2 = exc%val(2, 2)

        @:ASSERT(defined(exc))

        parallel_spin = calc_spin_raw(src1) == calc_spin_raw(src2)

        pgen_pick_elec = get_pgen_pick_biased_elecs(parallel_spin, pParallel, par_elec_pairs, AB_elec_pairs)

        tExchange = .not. parallel_spin &
                    .and. (this%GAS_spec%recoupling() .or. this%GAS_table(src1) == this%GAS_table(src2))

        iGAS = this%GAS_table(tgt1)
        nEmpty = this%GAS_size(iGAS) &
                 - sum(popCnt(iand(ilutI(0:NIfD), this%GAS_spin_orbs(0:NIfD, iGAS, get_spin(tgt1)))))
        if (nEmpty == 0) then
            pgen_first_pick(1) = 1.0_dp
        else
            ! adjust pgen
            pgen_first_pick(1) = 1.0_dp / real(nEmpty, dp)
        end if

        pgen_second_pick(1) = this%get_pgen_pick_weighted_hole( &
                              det_I%idx, DoubleExc_t(src1, tgt1, src2, tgt2))
        if (this%GAS_table(src1) == this%GAS_table(src2)) then
            pgen_first_pick(2) = this%get_pgen_pick_hole_from_active_space( &
                                 ilutI, this%GAS_table(src2), get_spin(tgt2))
            pgen_second_pick(2) = this%get_pgen_pick_weighted_hole( &
                                  det_I%idx, DoubleExc_t(src1, tgt2, src2, tgt1))
        else
            pgen_first_pick(2) = 0.0_dp
            pgen_second_pick(2) = 0.0_dp
        end if
        pgen_first_pick = merge(0.5_dp, 1.0_dp, tExchange) * pgen_first_pick

        pgen = pgen_pick_elec * sum(pgen_first_pick * pgen_second_pick)
    end function

end module gasci_disconnected
