#:include "../macros.fpph"
#include "macros.h"
module guga_gas_pchb_singles_main
    use util_mod, only: EnumBase_t, stop_all, operator(.implies.), &
        operator(.div.), operator(.isclose.)
    use SystemData, only: nEl, nSpatOrbs
    use excitation_generators, only: SingleExcitationGenerator_t
    use constants, only: n_int, dp, maxExcit, bits_n_int, stdout
    use dSFMT_interface, only: genrand_real2_dSFMT
    use SymExcitDataMod, only: ScratchSize
    use gasci, only: GASSpec_t
    use get_excit, only: make_single
    use gasci_supergroup_index, only: SuperGroupIndexer_t, lookup_supergroup_indexer
    use excitation_types, only: Excite_1_t, spin_allowed
    use FciMCData, only: excit_gen_store_type, GAS_PCHB_init_time
    use fortran_strings, only: to_upper
    use bit_rep_data, only: NIfTot, nIfD, GugaBits
    use bit_reps, only: decode_bit_det

    use guga_base_class, only: SinglesGUGABase_t
    use guga_bitrepops, only: CSF_Info_t
    use guga_data, only: ExcitationInformation_t, gen_type
    use guga_excitations, only: assign_excitinfo_values_single, &
        checkcompatibility_single
    better_implicit_none

    private
    public :: allocate_and_init, GAS_PCHB_SinglesOptions_t, GAS_PCHB_SinglesOptions_vals_t

    type, extends(EnumBase_t) :: GAS_PCHB_SinglesAlgorithm_t
        character(9) :: str
    end type

    type :: GAS_PCHB_SinglesAlgorithm_vals_t
        type(GAS_PCHB_SinglesAlgorithm_t) :: &
            UNIF_UNIF = GAS_PCHB_SinglesAlgorithm_t(1, 'UNIF:UNIF')
    contains
        procedure, nopass :: from_str => alg_from_str
    end type

    type, extends(EnumBase_t) :: GAS_PCHB_SinglesWeighting_t
        character(9) :: str
    end type


    type :: GAS_PCHB_SinglesWeighting_vals_t
        type(GAS_PCHB_SinglesWeighting_t) :: &
            UNIFORM = GAS_PCHB_SinglesWeighting_t(1, 'UNIFORM')
    contains
        procedure, nopass :: from_str => weight_from_str
    end type

    type :: GAS_PCHB_SinglesOptions_t
        type(GAS_PCHB_SinglesAlgorithm_t) :: algorithm
        type(GAS_PCHB_SinglesWeighting_t) :: weighting
    contains
        procedure :: is_valid
    end type

    type :: GAS_PCHB_SinglesOptions_vals_t
        type(GAS_PCHB_SinglesAlgorithm_vals_t) :: &
            algorithm = GAS_PCHB_SinglesAlgorithm_vals_t()
        type(GAS_PCHB_SinglesWeighting_vals_t) :: &
            weighting = GAS_PCHB_SinglesWeighting_vals_t()
    end type

    type(GAS_PCHB_SinglesOptions_vals_t), parameter :: &
        option_vals = GAS_PCHB_SinglesOptions_vals_t()


    !> The precomputed GAS uniform excitation generator
    type, extends(SinglesGUGABase_t) :: GAS_UniformExcGenerator_t
        private
        class(GASSpec_t), allocatable :: GAS_spec
        ! This is only a pointer because components cannot be targets
        ! otherwise. :-(
        type(SuperGroupIndexer_t), pointer :: indexer => null()
        logical, public :: use_lookup = .false.
            !! Use a lookup for the supergroup index in global_det_data.
        logical, public :: create_lookup = .false.
            !! Create **and** manage! the supergroup index lookup in global_det_data.

        integer(n_int), allocatable :: allowed_holes(:, :, :)
            !! Bitmask for the allowed holes for a single excitation
        integer :: L_spat_bits = -1
            !! The number of  integer(n_int) numbers to
            !! store information about spatial orbitals as bitmask
    contains
        private
        procedure, public :: init
        procedure, public :: finalize => GAS_singles_uniform_finalize

        procedure :: get_valid_orbs

        procedure, public :: calc_pgen => calc_pgen_uniform_singles
        procedure, public :: pickOrbitals => pickOrbitals_single
    end type

contains

    elemental logical function is_valid(this)
        class(GAS_PCHB_SinglesOptions_t), intent(in) :: this

        if (this%algorithm == option_vals%algorithm%UNIF_UNIF &
                .neqv. this%weighting == option_vals%weighting%UNIFORM) then
            is_valid = .false.
        else
            is_valid = .true.
        end if
    end function


    subroutine allocate_and_init(GAS_spec, options, use_lookup, create_lookup, generator)
        class(GASSpec_t), intent(in) :: GAS_spec
        type(GAS_PCHB_SinglesOptions_t), intent(in) :: options
        logical, intent(in) :: use_lookup, create_lookup
        class(SinglesGUGABase_t), allocatable, intent(inout) :: generator
        routine_name("allocate_and_init")

        @:ASSERT(options%is_valid())

        if (allocated(generator)) then
            call generator%finalize()
            deallocate(generator)
        end if
        if (options%algorithm == option_vals%algorithm%UNIF_UNIF) then
            allocate(GAS_UniformExcGenerator_t :: generator)
        end if
        select type(generator)
        type is(GAS_UniformExcGenerator_t)
            call generator%init(GAS_spec, use_lookup, create_lookup)
        end select
    end subroutine



    subroutine init(this, GAS_spec, use_lookup, create_lookup)
        class(GAS_UniformExcGenerator_t), intent(inout) :: this
        class(GASSpec_t), intent(in) :: GAS_spec
        logical, intent(in) :: use_lookup, create_lookup
        integer, allocatable :: supergroups(:, :)
        character(*), parameter :: this_routine = 'GAS_singles_uniform_init'

        integer :: i_sg, src, tgt

        this%GAS_spec = GAS_spec
        allocate(this%indexer, source=SuperGroupIndexer_t(GAS_spec, nEl))
        this%create_lookup = create_lookup
        if (create_lookup) then
            if (associated(lookup_supergroup_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the supergroup lookup.')
            else
                write(stdout, *) 'GAS singles is creating and managing the supergroup lookup'
                lookup_supergroup_indexer => this%indexer
            end if
        end if
        this%use_lookup = use_lookup
        if (use_lookup) write(stdout, *) 'GAS singles is using the supergroup lookup'
        ! possible supergroups
        supergroups = this%indexer%get_supergroups()

        this%L_spat_bits = nSpatOrbs .div. bits_n_int

        allocate(this%allowed_holes(0:this%L_spat_bits, nSpatOrbs, size(supergroups, 2)), &
                 source=0_n_int)


        do i_sg = 1, size(supergroups, 2)
            do src = 1, nSpatOrbs
                do tgt = 1, nSpatOrbs
                block
                    type(Excite_1_t) :: exc
                    exc = Excite_1_t(2 * src, 2 * tgt)
                    if (src /= tgt &
                            .and. this%GAS_spec%is_allowed(exc, supergroups(:, i_sg)) &
                            .and. symmetry_allowed(exc)) then
                        call my_set_orb(this%allowed_holes(:, src, i_sg), tgt)
                    end if
                end block
                end do
            end do
        end do

    contains

        ! Ugly Fortran syntax rules forces us to not use the macro.
        ! Even associate would not help here
        ! https://stackoverflow.com/questions/65734764/non-one-indexed-array-from-associate
        pure subroutine my_set_orb(ilut, orb)
            integer(n_int), intent(inout) :: ilut(0:this%L_spat_bits)
            integer, intent(in) :: orb
            set_orb(ilut, orb)
        end subroutine

        logical pure function symmetry_allowed(exc)
            use SymExcitDataMod, only: SpinOrbSymLabel
            type(Excite_1_t), intent(in) :: exc
            symmetry_allowed = SpinOrbSymLabel(exc%val(1, 1)) == SpinOrbSymLabel(exc%val(2, 1))
        end function
    end subroutine

    subroutine GAS_singles_uniform_finalize(this)
        class(GAS_UniformExcGenerator_t), intent(inout) :: this
        deallocate(this%allowed_holes)
    end subroutine


    subroutine pickOrbitals_single(this, ilut, nI, csf_i, excitInfo, pgen)
        debug_function_name("pickOrbitals_single")
        class(GAS_UniformExcGenerator_t), intent(in) :: this
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: src_orb, tgt_orb, i
        integer, allocatable :: occ_spat_orbs(:)
            !! The occupied spatial orbitals
        integer, allocatable :: valid_orbs(:)
        unused_var(ilut); unused_var(nI)

        occ_spat_orbs = pack([(i, i=1, size(csf_i%Occ_int))], csf_i%Occ_int /= 0)
        src_orb = occ_spat_orbs(1 + int(genrand_real2_dSFMT() * size(occ_spat_orbs)))

        valid_orbs = this%get_valid_orbs(csf_i, src_orb)

        if (size(valid_orbs) == 0) then
            pgen = 0.0_dp
            tgt_orb = 0
            excitInfo%valid = .false.
        else
            pgen = 1._dp / real(size(occ_spat_orbs) * size(valid_orbs), dp)
            tgt_orb = valid_orbs(1 + floor(genrand_real2_dSFMT() * size(valid_orbs)))


            if (tgt_orb < src_orb) then
                excitInfo = assign_excitinfo_values_single( &
                            gen_type%R, tgt_orb, src_orb, tgt_orb, src_orb)
            else
                excitInfo = assign_excitinfo_values_single( &
                            gen_type%L, tgt_orb, src_orb, src_orb, tgt_orb)
            end if
            excitInfo%valid = checkCompatibility_single(csf_i, excitInfo)
            if (excitInfo%valid) then
                @:ASSERT(pgen .isclose. this%calc_pgen(nI, ilut, csf_i, excitInfo))
            else
                pgen = 0._dp
            end if
        end if
    end subroutine pickOrbitals_single


    function get_valid_orbs(this, csf_i, src_orb) result(valid_orbs)
        class(GAS_UniformExcGenerator_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: src_orb
        integer, allocatable :: valid_orbs(:)
        character(*), parameter :: this_routine = 'get_valid_orbs'

        integer(n_int) :: &
            ilut_doubly_occupied(0:this%L_spat_bits), ilut_valid(0:this%L_spat_bits)
        integer :: i_sg

        if (this%use_lookup) then
            ! TODO implement
            call stop_all(this_routine, 'has to be implemented')
        else
            i_sg = this%indexer%idx_csf(csf_i)
        end if

        block
            integer, allocatable :: doubly_occupied(:)
            integer :: i
            doubly_occupied = pack([(i, i=1, size(csf_i%stepvector))], csf_i%stepvector == 3)
            ilut_doubly_occupied = 0_n_int
            do i = 1, size(doubly_occupied)
                set_orb(ilut_doubly_occupied, doubly_occupied(i))
            end do
        endblock

        ilut_valid = iand(this%allowed_holes(:, src_orb, i_sg), not(ilut_doubly_occupied))
        allocate(valid_orbs(sum(popcnt(ilut_valid))))
        if (size(valid_orbs) /= 0) then
            call decode_bit_det(valid_orbs, ilut_valid)
        end if
    end function get_valid_orbs


    function calc_pgen_uniform_singles(this, nI, ilutI, csf_i, excitInfo) result(pgen)
        debug_function_name("calc_pgen_uniform_singles")
        class(GAS_UniformExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen

        integer :: nOcc, nUnocc
        integer, allocatable :: valid_orbs(:)

        @:unused_var(nI, ilutI)

        ASSERT(1 <= excitInfo%i .and. excitInfo%i <= nSpatOrbs)
        ASSERT(1 <= excitInfo%j .and. excitInfo%j <= nSpatOrbs)

        if (excitInfo%i == excitInfo%j &
               .or. csf_i%stepvector(excitInfo%i) == 3 &
               .or. csf_i%stepvector(excitInfo%j) == 0) then
            pgen = 0.0_dp
        else
            nOcc = count(csf_i%Occ_int /= 0)
            valid_orbs = this%get_valid_orbs(csf_i, excitInfo%j)
            nUnocc = size(valid_orbs)
            pgen = 1.0_dp / real(nOcc * nUnocc, dp)
        end if
    end function calc_pgen_uniform_singles

    elemental function alg_from_str(str) result(res)
        character(*), intent(in) :: str
        type(GAS_PCHB_SinglesAlgorithm_t) :: res
        routine_name("alg_from_str")

        select case(to_upper(str))
        case(option_vals%algorithm%UNIF_UNIF%str)
            res = option_vals%algorithm%UNIF_UNIF
        case default
            call stop_all(this_routine, 'Invalid input: '//str)
        end select
    end function

    elemental function weight_from_str(str) result(res)
        character(*), intent(in) :: str
        type(GAS_PCHB_SinglesWeighting_t) :: res
        routine_name("alg_from_str")

        select case(to_upper(str))
        case(option_vals%weighting%UNIFORM%str)
            res = option_vals%weighting%UNIFORM
        case default
            call stop_all(this_routine, 'Invalid input: '//str)
        end select
    end function

end module
