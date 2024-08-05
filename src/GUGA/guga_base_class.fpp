#include "macros.h"
#:include "../macros.fpph"
#:include "../algorithms.fpph"

module guga_base_class
    !! This contains the base classes for the **orbital picking**.
    use constants, only: n_int, dp, maxExcit, int64, stdout, int_rdm
    use bit_rep_data, only: IlutBits, GugaBits, nIftot
    use SystemData, only: nel
    use FciMCData, only: pSingles, pDoubles
    use guga_data, only: ExcitationInformation_t, gen_type, excit_type, DistinctDouble_t, &
        ijkl_Index_t, distinct_doubles
    use guga_bitRepOps, only: CSF_Info_t
    use guga_excitations, only: get_excit_level_from_excitInfo
    use sets_mod, only: is_sorted
    use sort_mod, only: sort
    use error_handling_neci, only: stop_all
    use util_mod, only: EnumBase_t, operator(.isclose.), binary_search_int, near_zero, operator(.div.)

    better_implicit_none
    private
    public :: GUGABase_t, AbInitGUGABase_t, SinglesGUGABase_t, DoublesGUGABase_t, &
        GUGA_ExcGenerator_t, GUGA_ExcGenerator_vals_t, GUGA_exc_generator_vals, &
        GUGA_selected_exc_generator, test_pgen, Entry_t, analyze_entries, AnalyzedEntries_t, output


    type, extends(EnumBase_t) :: GUGA_ExcGenerator_t
    end type

    type :: GUGA_ExcGenerator_vals_t
        type(GUGA_ExcGenerator_t) :: &
            FCI_PCHB = GUGA_ExcGenerator_t(1), &
            GAS_PCHB = GUGA_ExcGenerator_t(2)
    end type

    type(GUGA_ExcGenerator_t), allocatable :: GUGA_selected_exc_generator
    type(GUGA_ExcGenerator_vals_t), parameter :: GUGA_exc_generator_vals = GUGA_ExcGenerator_vals_t()

    public :: init_guga_class, finalize_guga_class
    public :: class_pick_singles, class_pick_doubles
    public :: class_calc_orb_pgen_contr, class_calc_orb_pgen_contr_str, class_calc_orb_pgen_contr_end
    public :: class_calc_pgen
    public :: GUGA_exc_generator

    type, abstract :: Finalizable_t
    contains
        procedure(BoundFinalize_t), deferred :: finalize
    end type

    type, extends(Finalizable_t), abstract :: GUGABase_t
    contains
        procedure(PickOrbitals_t), deferred :: pickOrbitals_single
        procedure(PickOrbitals_t), deferred :: pickOrbitals_double

        procedure(CalcOrbitalPgenContr_t), deferred :: calc_orbital_pgen_contr_start
        procedure(CalcOrbitalPgenContr_t), deferred :: calc_orbital_pgen_contr_end
        procedure(calc_orbital_pgen_contr_t), deferred :: calc_orbital_pgen_contr

        procedure(CalcPgen_t), deferred :: calc_pgen

    end type

    type, extends(Finalizable_t), abstract :: SinglesGUGABase_t
    contains
        procedure(SinglesPickOrbitals_t), deferred :: pickOrbitals
        procedure(SinglesCalcPgen_t), deferred :: calc_pgen
    end type

    type, extends(Finalizable_t), abstract :: DoublesGUGABase_t
    contains
        procedure(DoublesPickOrbitals_t), deferred :: pickOrbitals
        procedure(DoublesCalcPgen_t), deferred :: calc_pgen


        procedure(DoublesCalcOrbitalPgenContr_t), deferred :: calc_orbital_pgen_contr_start
        procedure(DoublesCalcOrbitalPgenContr_t), deferred :: calc_orbital_pgen_contr_end
        procedure(Doublescalc_orbital_pgen_contr_t), deferred :: calc_orbital_pgen_contr

        procedure(Doubles_gen_all_distinct_doubles_t), deferred :: gen_all_distinct
    end type

    type, abstract, extends(GUGABase_t) :: AbInitGUGABase_t
        class(SinglesGUGABase_t), allocatable :: singles_generator
        class(DoublesGUGABase_t), allocatable :: doubles_generator
    contains
        procedure :: pickOrbitals_single => AbInitGUGABase_t_pickOrbitals_single
        procedure :: pickOrbitals_double  => AbInitGUGABase_t_pickOrbitals_double

        procedure :: calc_orbital_pgen_contr_start => AbInitGUGABase_t_calc_orbital_pgen_contr_start
        procedure :: calc_orbital_pgen_contr_end => AbInitGUGABase_t_calc_orbital_pgen_contr_end
        procedure :: calc_orbital_pgen_contr => AbInitGUGABase_t_calc_orbital_pgen_contr

        procedure :: calc_pgen => AbInitGUGABase_t_calc_pgen
    end type


    interface
        ! This is implemented in a submodule
        module subroutine init_guga_class()
        end subroutine
    end interface

    class(GUGABase_t), allocatable :: GUGA_exc_generator

    abstract interface
        subroutine PickOrbitals_t(this, nI, ilut, csf_i, excitInfo, pgen)
            import :: GUGABase_t, dp, n_int, GugaBits, CSF_Info_t, ExcitationInformation_t, nel
            implicit none
            class(GUGABase_t), intent(in) :: this
            integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
            integer, intent(in) :: nI(nel)
            type(CSF_Info_t), intent(in) :: csf_i
            type(ExcitationInformation_t), intent(out) :: excitInfo
            real(dp), intent(out) :: pgen
        end subroutine PickOrbitals_t

        subroutine SinglesPickOrbitals_t(this, nI, ilut, csf_i, excitInfo, pgen)
            import :: SinglesGUGABase_t, dp, n_int, GugaBits, CSF_Info_t, ExcitationInformation_t, nel
            implicit none
            class(SinglesGUGABase_t), intent(in) :: this
            integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
            integer, intent(in) :: nI(nel)
            type(CSF_Info_t), intent(in) :: csf_i
            type(ExcitationInformation_t), intent(out) :: excitInfo
            real(dp), intent(out) :: pgen
        end subroutine

        subroutine DoublesPickOrbitals_t(this, nI, ilut, csf_i, excitInfo, pgen)
            import :: DoublesGUGABase_t, dp, n_int, GugaBits, CSF_Info_t, ExcitationInformation_t, nel
            implicit none
            class(DoublesGUGABase_t), intent(in) :: this
            integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
            integer, intent(in) :: nI(nel)
            type(CSF_Info_t), intent(in) :: csf_i
            type(ExcitationInformation_t), intent(out) :: excitInfo
            real(dp), intent(out) :: pgen
        end subroutine

        subroutine CalcOrbitalPgenContr_t(this, nI, ilut, csf_i, occ_orbs, orb_a, orb_pgen)
            import :: GUGABase_t, dp, CSF_Info_t, nEl, GugaBits, n_int
            implicit none
            class(GUGABase_t), intent(in) :: this
            integer, intent(in) :: nI(nEl)
            integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: occ_orbs(2), orb_a
            real(dp), intent(out) :: orb_pgen
        end subroutine CalcOrbitalPgenContr_t

        subroutine calc_orbital_pgen_contr_t(this, nI, ilut, csf_i, occ_orbs, above_cpt, below_cpt)
            import :: GUGABase_t, dp, CSF_Info_t, nEl, GugaBits, n_int
            implicit none
            class(GUGABase_t), intent(in) :: this
            integer, intent(in) :: nI(nEl)
            integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: occ_orbs(2)
            real(dp), intent(out) :: above_cpt, below_cpt
        end subroutine calc_orbital_pgen_contr_t


        subroutine DoublesCalcOrbitalPgenContr_t(this, nI, ilut, csf_i, occ_orbs, orb_a, orb_pgen)
            import :: DoublesGUGABase_t, dp, CSF_Info_t, nEl, GugaBits, n_int
            implicit none
            class(DoublesGUGABase_t), intent(in) :: this
            integer, intent(in) :: nI(nEl)
            integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: occ_orbs(2), orb_a
            real(dp), intent(out) :: orb_pgen
        end subroutine

        subroutine Doublescalc_orbital_pgen_contr_t(this, nI, ilut, csf_i, occ_orbs, above_cpt, below_cpt)
            import :: DoublesGUGABase_t, dp, CSF_Info_t, nEl, GugaBits, n_int
            implicit none
            class(DoublesGUGABase_t), intent(in) :: this
            integer, intent(in) :: nI(nEl)
            integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: occ_orbs(2)
            real(dp), intent(out) :: above_cpt, below_cpt
        end subroutine

        pure function Doubles_gen_all_distinct_doubles_t(this, csf_i, sort_fused_key) result(res)
            !! This function returns all possible excitations for the orbital picker
            !!
            !! This means it considers only occupation, not spin recoupling.
            import :: DoublesGUGABase_t, CSF_Info_t, DistinctDouble_t
            implicit none
            class(DoublesGUGABase_t), intent(in) :: this
            type(CSF_Info_t), intent(in) :: csf_i
            logical, optional, intent(in) :: sort_fused_key
                !! The fused index of the excitations should be sorted.
            type(DistinctDouble_t), allocatable :: res(:)
        end function


        function CalcPgen_t(this, nI, ilutI, csf_i, excitInfo) result(pgen)
            import :: GUGABase_t, dp, CSF_Info_t, n_int, GUGABits, ExcitationInformation_t, nEl
            implicit none
            class(GUGABase_t), intent(in) :: this
            integer, intent(in) :: nI(nel)
            integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
            type(CSF_Info_t), intent(in) :: csf_i
            type(ExcitationInformation_t), intent(in) :: excitInfo
            real(dp) :: pgen
        end function

        function SinglesCalcPgen_t(this, nI, ilutI, csf_i, excitInfo) result(pgen)
            import :: SinglesGUGABase_t, dp, CSF_Info_t, n_int, GUGABits, ExcitationInformation_t, nEl
            implicit none
            class(SinglesGUGABase_t), intent(in) :: this
            integer, intent(in) :: nI(nel)
            integer(n_int), intent(in) :: ilutI(0 : GugaBits%len_tot)
            type(CSF_Info_t), intent(in) :: csf_i
            type(ExcitationInformation_t), intent(in) :: excitInfo
            real(dp) :: pgen
        end function

        function DoublesCalcPgen_t(this, nI, ilutI, csf_i, excitInfo) result(pgen)
            import :: DoublesGUGABase_t, dp, CSF_Info_t, n_int, GUGABits, ExcitationInformation_t, nEl
            implicit none
            class(DoublesGUGABase_t), intent(in) :: this
            integer, intent(in) :: nI(nel)
            integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
            type(CSF_Info_t), intent(in) :: csf_i
            type(ExcitationInformation_t), intent(in) :: excitInfo
            real(dp) :: pgen
        end function

        subroutine BoundFinalize_t(this)
            import :: Finalizable_t
            implicit none
            class(Finalizable_t), intent(inout) :: this
        end subroutine
    end interface

    type :: Entry_t
        integer(int64) :: key
        type(DistinctDouble_t) :: exc_operator
        integer(int64) :: accum = 0_int64
        real(dp) :: p_gen = -1._dp
        type(ExcitationInformation_t) :: exc_info = ExcitationInformation_t()
    end type

    type :: AnalyzedEntries_t
        integer(int64) :: key
        type(DistinctDouble_t) :: exc
        real(dp) :: p_gen
        real(dp) :: quality
        integer(int64) :: accum
    end type

    abstract interface
        logical pure function compare_entry_t(p1, p2)
            import :: AnalyzedEntries_t
            type(AnalyzedEntries_t), intent(in) :: p1, p2
        end function
    end interface

contains

    subroutine finalize_guga_class()
        if (allocated(GUGA_exc_generator)) then
            call GUGA_exc_generator%finalize()
            deallocate(GUGA_exc_generator)
        end if
    end subroutine

    subroutine class_pick_singles(nI, ilut, csf_i, excitInfo, pgen)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call GUGA_exc_generator%pickOrbitals_single(nI, ilut, csf_i, excitInfo, pgen)
    end subroutine

    subroutine class_pick_doubles(nI, ilut, csf_i, excitInfo, pgen)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call GUGA_exc_generator%pickOrbitals_double(nI, ilut, csf_i, excitInfo, pgen)
    end subroutine

    subroutine class_calc_orb_pgen_contr(nI, ilut, csf_i, occ_orbs, cpt_a, cpt_b)
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: cpt_a, cpt_b
        call GUGA_exc_generator%calc_orbital_pgen_contr(nI, ilut, csf_i, occ_orbs, cpt_a, cpt_b)
    end subroutine class_calc_orb_pgen_contr


    ! i think it would be better if i 'just' reimplement:
    subroutine class_calc_orb_pgen_contr_str(nI, ilut, csf_i, occ_orbs, a, orb_pgen)
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), a
        real(dp), intent(out) :: orb_pgen
        call GUGA_exc_generator%calc_orbital_pgen_contr_start(nI, ilut, csf_i, occ_orbs, a, orb_pgen)
    end subroutine

    subroutine class_calc_orb_pgen_contr_end(nI, ilut, csf_i, occ_orbs, a, orb_pgen)
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), a
        real(dp), intent(out) :: orb_pgen
        call GUGA_exc_generator%calc_orbital_pgen_contr_end(nI, ilut, csf_i, occ_orbs, a, orb_pgen)
    end subroutine class_calc_orb_pgen_contr_end

    function class_calc_pgen(nI, ilutI, csf_i, excitInfo) result(pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen
        pgen = GUGA_exc_generator%calc_pgen(nI, ilutI, csf_i, excitInfo)
    end function

    subroutine AbInitGUGABase_t_pickOrbitals_single(this, nI, ilut, csf_i, excitInfo, pgen)
        class(AbInitGUGABase_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call this%singles_generator%pickOrbitals(nI, ilut, csf_i, excitInfo, pgen)
    end subroutine

    subroutine AbInitGUGABase_t_pickOrbitals_double(this, nI, ilut, csf_i, excitInfo, pgen)
        class(AbInitGUGABase_t), intent(in) :: this
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call this%doubles_generator%pickOrbitals(nI, ilut, csf_i, excitInfo, pgen)
    end subroutine

    subroutine AbInitGUGABase_t_calc_orbital_pgen_contr_start(this, nI, ilut, csf_i, occ_orbs, orb_a, orb_pgen)
        class(AbInitGUGABase_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        call this%doubles_generator%calc_orbital_pgen_contr_start(nI, ilut, csf_i, occ_orbs, orb_a, orb_pgen)
    end subroutine

    subroutine AbInitGUGABase_t_calc_orbital_pgen_contr_end(this, nI, ilut, csf_i, occ_orbs, orb_a, orb_pgen)
        class(AbInitGUGABase_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        call this%doubles_generator%calc_orbital_pgen_contr_end(nI, ilut, csf_i, occ_orbs, orb_a, orb_pgen)
    end subroutine

    subroutine AbInitGUGABase_t_calc_orbital_pgen_contr(this, nI, ilut, csf_i, occ_orbs, above_cpt, below_cpt)
        class(AbInitGUGABase_t), intent(in) :: this
        integer, intent(in) :: nI(nEl)
        integer(n_int), intent(in) :: ilut(0 : GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: above_cpt, below_cpt
        call this%doubles_generator%calc_orbital_pgen_contr(nI, ilut, csf_i, occ_orbs, above_cpt, below_cpt)
    end subroutine

    function AbInitGUGABase_t_calc_pgen(this, nI, ilutI, csf_i, excitInfo) result(pgen)
        class(AbInitGUGABase_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen
        routine_name("AbInitGUGABase_t_calc_pgen")

        integer :: ic

        ic = get_excit_level_from_excitInfo(excitInfo)

        select case(ic)
        case(1)
            pgen = pSingles * this%singles_generator%calc_pgen(nI, ilutI, csf_i, excitInfo)
        case(2)
            pgen = pDoubles * this%doubles_generator%calc_pgen(nI, ilutI, csf_i, excitInfo)
        case(3:)
            pgen = 0.0_dp
        case default
            call stop_all(this_routine, 'Should not be here')
        end select
    end function


    subroutine test_pgen(doubles_generator, csf_i, n_iter, entries)
        class(DoublesGUGABase_t), intent(in) :: doubles_generator
        type(CSF_Info_t), intent(in) :: csf_i
        integer(int64), intent(in) :: n_iter
        type(Entry_t), intent(out), allocatable :: entries(:)
        routine_name("test_pgen")

        integer :: nI(nEl)
        integer(int64) :: i, L
        integer(n_int) :: ilut(0 : niftot)

        nI = csf_i%to_nI()
        ilut = csf_i%to_ilut()

        block
            type(DistinctDouble_t), allocatable :: exc_operators(:)
            exc_operators = doubles_generator%gen_all_distinct(csf_i, sort_fused_key=.true.)
            entries = [(Entry_t(exc_operators(i)%fuse(), exc_operators(i)), i = 1, size(exc_operators))]
            @:ASSERT(is_sorted(entries%key))
        end block

        do i = 1, size(entries)
            entries(i)%p_gen = doubles_generator%calc_pgen(nI, ilut, csf_i, entries(i)%exc_operator%get_repr_excit_info(1))
        end do

        L = n_iter .div. 100_int64
        write(stdout, '(A)') 'Running pgen tester'
        do i = 1_int64, n_iter
        block
            type(ExcitationInformation_t) :: excit_info
            type(DistinctDouble_t) :: exc
            integer(int64) :: pos
            real(dp) :: pgen
            if (mod(i, L) == 0_int64) then
                write(stdout, '(A)', advance='no') '#'
                flush(stdout)
            end if

            call doubles_generator%pickOrbitals(nI, ilut, csf_i, excit_info, pgen)
            if (.not. excit_info%valid) cycle

            if (.not. (pgen .isclose. doubles_generator%calc_pgen(nI, ilut, csf_i, excit_info))) then
                call stop_all(this_routine, 'p_gen does not match')
            end if

            associate(idx => ijkl_Index_t([excit_info%i, excit_info%j, excit_info%k, excit_info%l]))
                exc = DistinctDouble_t(idx%normalize())
                if (.not. (pgen .isclose. doubles_generator%calc_pgen(nI, ilut, csf_i, exc%get_repr_excit_info(1)))) then
                    call stop_all(this_routine, 'p_gen does not match')
            end if
            end associate

            if (exc%typ /= distinct_doubles%invalid) then
                pos = binary_search_int(entries%key, exc%fuse())
                if (pos == -1) then
                    call stop_all(this_routine, 'unexpected excitation')
                end if
                entries(pos)%accum = entries(pos)%accum + 1
                entries(pos)%exc_info = excit_info
                entries(pos)%p_gen = pgen
            end if
        end block
        end do
        write(stdout, *)

        contains
            logical elemental function leq_key(x, y)
            type(Entry_t), intent(in) :: x, y
            leq_key = x%key <= y%key
            end function
    end subroutine


    function analyze_entries(entries, n_iter, compare_entry) result(res)
        type(Entry_t), intent(in) :: entries(:)
        integer(int64), intent(in) :: n_iter
        procedure(compare_entry_t), optional :: compare_entry
        type(AnalyzedEntries_t) :: res(size(entries))
#ifdef IFORT_
        routine_name("analyze_entries")
#endif

        integer :: i
        integer(int64) :: total_hits

        total_hits = sum(entries%accum)


        do i = 1, size(entries)
            if (near_zero(entries(i)%p_gen)) then
                res(i) = AnalyzedEntries_t(entries(i)%key, entries(i)%exc_operator, &
                                entries(i)%p_gen, 1._dp, &
                                entries(i)%accum)
            else
                res(i) = AnalyzedEntries_t(entries(i)%key, entries(i)%exc_operator, &
                                entries(i)%p_gen, (entries(i)%accum / entries(i)%p_gen) / n_iter, &
                                entries(i)%accum)
            end if
        end do


        if (present(compare_entry)) then
#ifdef IFORT_
            call stop_all(this_routine, &
                "Does not work with Intel. &
                &At least with ifort 21 the compiler has an internal error.")
#else
            @:sort(AnalyzedEntries_t, res, rank=1, along=1, comp=compare_entry)
#endif
        end if

    end function

    subroutine output(file_id, analyzed_entries)
        integer, intent(in) :: file_id
        type(AnalyzedEntries_t), intent(in) :: analyzed_entries(:)

        integer :: i

        write(file_id, *) "=================================="
        write(file_id, *) "Problematic contribution List: "
        write(file_id, *) "=================================="
        write(file_id, '("|", A10, "|", A24, "|", A15, "|", A15, "|", A15, "|")') "Key", "Excitation", "Quality", "pgen", "n_accum"
        ! write(file_id, '("--------------------------------------------------------------------------")')
        write(file_id, '(A)') repeat('-', 10 + 24 + 15 + 15 + 15 + 6)
        do i = 1, size(analyzed_entries)
            write(file_id, '("|", I10, "|"I5, I5, " -> ", I5, I5, "|", F15.5, "|", E15.5, "|", I15"|")') &
                analyzed_entries(i)%key, analyzed_entries(i)%exc%idx%convert(), analyzed_entries(i)%quality, &
                analyzed_entries(i)%p_gen, analyzed_entries(i)%accum
        end do
    end subroutine

end module
