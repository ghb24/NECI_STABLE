#include "macros.h"
module guga_base_class
    use constants, only: n_int, dp, maxExcit, int64, stdout, int_rdm
    use bit_rep_data, only: IlutBits, GugaBits
    use SystemData, only: nel
    use FciMCData, only: pSingles, pDoubles
    use guga_data, only: ExcitationInformation_t, gen_type, excit_type
    use guga_bitRepOps, only: CSF_Info_t
    use guga_excitations, only: get_excit_level_from_excitInfo
    use error_handling_neci, only: stop_all
    use util_mod, only: EnumBase_t

    better_implicit_none
    private
    public :: GUGABase_t, AbInitGUGABase_t, SinglesGUGABase_t, DoublesGUGABase_t, &
        GUGA_ExcGenerator_t, GUGA_ExcGenerator_vals_t, GUGA_exc_generator_vals, &
        GUGA_selected_exc_generator


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
        subroutine PickOrbitals_t(this, ilut, nI, csf_i, excitInfo, pgen)
            import :: GUGABase_t, dp, n_int, GugaBits, CSF_Info_t, ExcitationInformation_t, nel
            implicit none
            class(GUGABase_t), intent(in) :: this
            integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
            integer, intent(in) :: nI(nel)
            type(CSF_Info_t), intent(in) :: csf_i
            type(ExcitationInformation_t), intent(out) :: excitInfo
            real(dp), intent(out) :: pgen
        end subroutine PickOrbitals_t

        subroutine SinglesPickOrbitals_t(this, ilut, nI, csf_i, excitInfo, pgen)
            import :: SinglesGUGABase_t, dp, n_int, GugaBits, CSF_Info_t, ExcitationInformation_t, nel
            implicit none
            class(SinglesGUGABase_t), intent(in) :: this
            integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
            integer, intent(in) :: nI(nel)
            type(CSF_Info_t), intent(in) :: csf_i
            type(ExcitationInformation_t), intent(out) :: excitInfo
            real(dp), intent(out) :: pgen
        end subroutine

        subroutine DoublesPickOrbitals_t(this, ilut, nI, csf_i, excitInfo, pgen)
            import :: DoublesGUGABase_t, dp, n_int, GugaBits, CSF_Info_t, ExcitationInformation_t, nel
            implicit none
            class(DoublesGUGABase_t), intent(in) :: this
            integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
            integer, intent(in) :: nI(nel)
            type(CSF_Info_t), intent(in) :: csf_i
            type(ExcitationInformation_t), intent(out) :: excitInfo
            real(dp), intent(out) :: pgen
        end subroutine

        subroutine CalcOrbitalPgenContr_t(this, csf_i, occ_orbs, orb_a, orb_pgen)
            import :: GUGABase_t, dp, CSF_Info_t
            implicit none
            class(GUGABase_t), intent(in) :: this
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: occ_orbs(2), orb_a
            real(dp), intent(out) :: orb_pgen
        end subroutine CalcOrbitalPgenContr_t

        subroutine calc_orbital_pgen_contr_t(this, csf_i, occ_orbs, above_cpt, below_cpt)
            import :: GUGABase_t, dp, CSF_Info_t
            implicit none
            class(GUGABase_t), intent(in) :: this
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: occ_orbs(2)
            real(dp), intent(out) :: above_cpt, below_cpt
        end subroutine calc_orbital_pgen_contr_t


        subroutine DoublesCalcOrbitalPgenContr_t(this, csf_i, occ_orbs, orb_a, orb_pgen)
            import :: DoublesGUGABase_t, dp, CSF_Info_t
            implicit none
            class(DoublesGUGABase_t), intent(in) :: this
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: occ_orbs(2), orb_a
            real(dp), intent(out) :: orb_pgen
        end subroutine

        subroutine Doublescalc_orbital_pgen_contr_t(this, csf_i, occ_orbs, above_cpt, below_cpt)
            import :: DoublesGUGABase_t, dp, CSF_Info_t
            implicit none
            class(DoublesGUGABase_t), intent(in) :: this
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: occ_orbs(2)
            real(dp), intent(out) :: above_cpt, below_cpt
        end subroutine


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
            integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
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

contains

    subroutine finalize_guga_class()
        call GUGA_exc_generator%finalize()
        deallocate(GUGA_exc_generator)
    end subroutine

    subroutine class_pick_singles(ilut, nI, csf_i, excitInfo, pgen)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call GUGA_exc_generator%pickOrbitals_single(ilut, nI, csf_i, excitInfo, pgen)
    end subroutine

    subroutine class_pick_doubles(ilut, nI, csf_i, excitInfo, pgen)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call GUGA_exc_generator%pickOrbitals_double(ilut, nI, csf_i, excitInfo, pgen)
    end subroutine

    subroutine class_calc_orb_pgen_contr(csf_i, occ_orbs, cpt_a, cpt_b)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: cpt_a, cpt_b
        call GUGA_exc_generator%calc_orbital_pgen_contr(csf_i, occ_orbs, cpt_a, cpt_b)
    end subroutine class_calc_orb_pgen_contr


    ! i think it would be better if i 'just' reimplement:
    subroutine class_calc_orb_pgen_contr_str(csf_i, occ_orbs, a, orb_pgen)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), a
        real(dp), intent(out) :: orb_pgen
        call GUGA_exc_generator%calc_orbital_pgen_contr_start(csf_i, occ_orbs, a, orb_pgen)
    end subroutine

    subroutine class_calc_orb_pgen_contr_end(csf_i, occ_orbs, a, orb_pgen)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), a
        real(dp), intent(out) :: orb_pgen
        call GUGA_exc_generator%calc_orbital_pgen_contr_end(csf_i, occ_orbs, a, orb_pgen)
    end subroutine class_calc_orb_pgen_contr_end

    function class_calc_pgen(nI, ilutI, csf_i, excitInfo) result(pgen)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen
        pgen = GUGA_exc_generator%calc_pgen(nI, ilutI, csf_i, excitInfo)
    end function

    subroutine AbInitGUGABase_t_pickOrbitals_single(this, ilut, nI, csf_i, excitInfo, pgen)
        class(AbInitGUGABase_t), intent(in) :: this
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call this%singles_generator%pickOrbitals(ilut, nI, csf_i, excitInfo, pgen)
    end subroutine

    subroutine AbInitGUGABase_t_pickOrbitals_double(this, ilut, nI, csf_i, excitInfo, pgen)
        class(AbInitGUGABase_t), intent(in) :: this
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        call this%doubles_generator%pickOrbitals(ilut, nI, csf_i, excitInfo, pgen)
    end subroutine

    subroutine AbInitGUGABase_t_calc_orbital_pgen_contr_start(this, csf_i, occ_orbs, orb_a, orb_pgen)
        class(AbInitGUGABase_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        call this%doubles_generator%calc_orbital_pgen_contr_start(csf_i, occ_orbs, orb_a, orb_pgen)
    end subroutine

    subroutine AbInitGUGABase_t_calc_orbital_pgen_contr_end(this, csf_i, occ_orbs, orb_a, orb_pgen)
        class(AbInitGUGABase_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen
        call this%doubles_generator%calc_orbital_pgen_contr_end(csf_i, occ_orbs, orb_a, orb_pgen)
    end subroutine

    subroutine AbInitGUGABase_t_calc_orbital_pgen_contr(this, csf_i, occ_orbs, above_cpt, below_cpt)
        class(AbInitGUGABase_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: above_cpt, below_cpt
        call this%doubles_generator%calc_orbital_pgen_contr(csf_i, occ_orbs, above_cpt, below_cpt)
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
end module
