#:include "../macros.fpph"
#include "macros.h"
module guga_gas_pchb_singles_main
    use util_mod, only: EnumBase_t, stop_all, operator(.implies.), &
        operator(.div.), operator(.isclose.), near_zero
    use SystemData, only: nEl, nSpatOrbs, nBasis
    use MPI_wrapper, only: root
    use constants, only: n_int, dp, bits_n_int, stdout
    use dSFMT_interface, only: genrand_real2_dSFMT
    use property_vector_index, only: AlsoGUGA_PropertyIndexer_t, lookup_property_indexer
    use excitation_types, only: Excite_1_t
    use fortran_strings, only: to_upper
    use bit_rep_data, only: nIfD, GugaBits
    use bit_reps, only: decode_bit_det, UnoccupiedGetter_t
    use aliasSampling, only: AliasSampler_1D_t, AliasSampler_2D_t, do_direct_calculation, AliasSampler_3D_t
    use UMatCache, only: gtID, get_umat_el

    use guga_base_class, only: SinglesGUGABase_t
    use guga_bitrepops, only: CSF_Info_t
    use guga_data, only: ExcitationInformation_t, gen_type
    use guga_pchb_singles_weights_mod, only: get_weight_t, PropVec_PCHB_SinglesWeighting_t, &
        PropVec_PCHB_SinglesWeighting_vals_t, set_get_weight_pointer
    use guga_excitations, only: assign_excitinfo_values_single, &
        checkcompatibility_single
    better_implicit_none


    private
    public :: allocate_and_init, PropVec_PCHB_SinglesOptions_t, PropVec_PCHB_SinglesOptions_vals_t, &
        scale_down_exch_singles

    type, extends(EnumBase_t) :: PropVec_PCHB_SinglesAlgorithm_t
        character(9) :: str
    end type

    real(dp) :: scale_down_exch_singles = 1._dp

    type :: PropVec_PCHB_SinglesAlgorithm_vals_t
        type(PropVec_PCHB_SinglesAlgorithm_t) :: &
            UNIF_UNIF = PropVec_PCHB_SinglesAlgorithm_t(1, 'UNIF:UNIF'), &
            FULL_FULL = PropVec_PCHB_SinglesAlgorithm_t(2, 'FULL:FULL'), &
            UNIF_FULL = PropVec_PCHB_SinglesAlgorithm_t(3, 'UNIF:FULL'), &
            UNIF_XNEW = PropVec_PCHB_SinglesAlgorithm_t(4, 'UNIF:XNEW')
    contains
        procedure, nopass :: from_str => alg_from_str
    end type

    type :: PropVec_PCHB_SinglesOptions_t
        type(PropVec_PCHB_SinglesAlgorithm_t) :: algorithm
        type(PropVec_PCHB_SinglesWeighting_t) :: weighting
    contains
        procedure :: is_valid
    end type

    type :: PropVec_PCHB_SinglesOptions_vals_t
        type(PropVec_PCHB_SinglesAlgorithm_vals_t) :: &
            algorithm = PropVec_PCHB_SinglesAlgorithm_vals_t()
        type(PropVec_PCHB_SinglesWeighting_vals_t) :: &
            weighting = PropVec_PCHB_SinglesWeighting_vals_t()
    end type

    type(PropVec_PCHB_SinglesOptions_vals_t), parameter :: &
        option_vals = PropVec_PCHB_SinglesOptions_vals_t()


    !> The precomputed PropVec uniform excitation generator
    type, extends(SinglesGUGABase_t) :: PropVec_UniformExcGenerator_t
        private
        class(AlsoGUGA_PropertyIndexer_t), allocatable :: indexer
        logical, public :: use_lookup = .false.
            !! Use a lookup for the prop_vec index in global_det_data.
        logical, public :: create_lookup = .false.
            !! Create **and** manage! the prop_vec index lookup in global_det_data.

        integer(n_int), allocatable :: allowed_holes(:, :, :)
            !! Bitmask for the allowed holes for a single excitation
        integer :: L_spat_bits = -1
            !! The number of  integer(n_int) numbers to
            !! store information about spatial orbitals as bitmask
    contains
        private
        procedure, public :: init => init_PropVec_UniformExcGenerator_t
        procedure, public :: finalize => PropVec_singles_uniform_finalize

        procedure :: get_valid_orbs

        procedure, public :: calc_pgen => calc_pgen_uniform_singles
        procedure, public :: pickOrbitals => pickOrbitals_single
    end type

    !> The precomputed PropVec uniform excitation generator
    type, extends(SinglesGUGABase_t), abstract :: PropVec_Weighted_t
        private
        class(AlsoGUGA_PropertyIndexer_t), allocatable :: indexer
        logical, public :: use_lookup = .false.
            !! Use a lookup for the prop_vec index in global_det_data.
        logical, public :: create_lookup = .false.
            !! Create **and** manage! the prop_vec index lookup in global_det_data.


        type(AliasSampler_1D_t) :: I_sampler
            !! The shape is (n_prop_vec) -> 2 * nSpatOrbs
            !!
            !! @note
            !! To make doubly occupied orbitals more probable we have doubly the
            !! amount of spatial orbitals and treat it as if there were spin-orbitals.
            !! @endnote
            !!
        type(AliasSampler_2D_t) :: A_sampler
            !! The shape is (nSpatOrbs, n_prop_vec) -> 2 * nSpatOrbs
            !!
            !! @note
            !! To make empty orbitals more probable than singly occpied ones,
            !! we have doubly the
            !! amount of spatial orbitals and treat it as if there were spin-orbitals.
            !! @endnote
            !!
        procedure(get_weight_t), pointer, nopass :: get_weight => null()

        type(UnoccupiedGetter_t) :: unoccupied

    contains
        private
        procedure, public :: init => init_PropVec_Weighted_t
        procedure, public :: finalize => finalize_PropVec_Weighted_t
    end type

    type, extends(PropVec_Weighted_t) :: SinglesPropVec_FullFull_t
    contains
        procedure, public :: calc_pgen => calc_pgen_FullFull
        procedure, public :: pickOrbitals => pickOrbitals_FullFull
    end type

    type, extends(PropVec_Weighted_t) :: SinglesPropVec_UnifFull_t
    contains
        procedure, public :: calc_pgen => calc_pgen_UnifFull
        procedure, public :: pickOrbitals => pickOrbitals_UnifFull
    end type

    !> The precomputed PropVec uniform excitation generator
    type, extends(SinglesGUGABase_t) :: UnifXnew_t
        private
        class(AlsoGUGA_PropertyIndexer_t), allocatable :: indexer
        logical, public :: use_lookup = .false.
        logical, public :: create_lookup = .false.


        type(AliasSampler_3D_t) :: A_sampler
            !! The shape is (nSpatOrbs: i, stepvector: i, n_prop_vec: i_sg) -> 2 * nSpatOrbs
        procedure(get_weight_t), pointer, nopass :: get_weight => null()

        type(UnoccupiedGetter_t) :: unoccupied

    contains
        private
        procedure, public :: init => init_UnifXnew_t
        procedure, public :: finalize => finalize_UnifXnew_t

        procedure, public :: pickOrbitals => pickOrbitals_UnifXnew_t
        procedure, public :: calc_pgen => calc_pgen_UnifXnew_t
    end type

contains

    elemental logical function is_valid(this)
        class(PropVec_PCHB_SinglesOptions_t), intent(in) :: this
        is_valid = (this%algorithm == option_vals%algorithm%UNIF_UNIF &
                            .implies. this%weighting == option_vals%weighting%UNIFORM)
    end function


    subroutine allocate_and_init(indexer, options, use_lookup, create_lookup, generator)
        class(AlsoGUGA_PropertyIndexer_t), intent(in) :: indexer
        type(PropVec_PCHB_SinglesOptions_t), intent(in) :: options
        logical, intent(in) :: use_lookup, create_lookup
        class(SinglesGUGABase_t), allocatable, intent(inout) :: generator
        routine_name("allocate_and_init")

        @:ASSERT(options%is_valid())

        if (allocated(generator)) then
            call generator%finalize()
            deallocate(generator)
        end if

        if (options%algorithm == option_vals%algorithm%UNIF_UNIF) then
            allocate(PropVec_UniformExcGenerator_t :: generator)
            select type(generator)
            type is(PropVec_UniformExcGenerator_t)
                call generator%init(indexer, use_lookup, create_lookup)
            end select
        else if (options%algorithm == option_vals%algorithm%FULL_FULL) then
            allocate(SinglesPropVec_FullFull_t :: generator)
            select type(generator)
            type is(SinglesPropVec_FullFull_t)
                call generator%init(options%weighting, indexer, use_lookup, create_lookup)
            end select
        else if (options%algorithm == option_vals%algorithm%UNIF_FULL) then
            allocate(SinglesPropVec_UnifFull_t :: generator)
            select type(generator)
            type is(SinglesPropVec_UnifFull_t)
                call generator%init(options%weighting, indexer, use_lookup, create_lookup)
            end select
        else if (options%algorithm == option_vals%algorithm%UNIF_XNEW) then
            allocate(UnifXnew_t :: generator)
            select type(generator)
            type is(UnifXnew_t)
                call generator%init(options%weighting, indexer, use_lookup, create_lookup, scale_down_exch_singles)
            end select
        else
            call stop_all(this_routine, "unsupported singles option.")
        end if
    end subroutine



    subroutine init_PropVec_UniformExcGenerator_t(this, indexer, use_lookup, create_lookup)
        class(PropVec_UniformExcGenerator_t), intent(inout) :: this
        class(AlsoGUGA_PropertyIndexer_t), intent(in) :: indexer
        logical, intent(in) :: use_lookup, create_lookup
        integer, allocatable :: prop_vecs(:, :)
        routine_name('init_PropVec_UniformExcGenerator_t')


        integer :: i_sg, src, tgt

        allocate(this%indexer, source=indexer)
        this%create_lookup = create_lookup
        if (create_lookup) then
            if (allocated(lookup_property_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the prop_vec lookup.')
            else
                write(stdout, *) 'PropVec singles is creating and managing the prop_vec lookup'
                allocate(lookup_property_indexer, source=this%indexer)
            end if
        end if
        this%use_lookup = use_lookup
        if (use_lookup) write(stdout, *) 'PropVec singles is using the prop_vec lookup'
        ! possible prop_vecs
        prop_vecs = this%indexer%get_prop_vecs()

        this%L_spat_bits = nSpatOrbs .div. bits_n_int

        allocate(this%allowed_holes(0:this%L_spat_bits, nSpatOrbs, size(prop_vecs, 2)), &
                 source=0_n_int)


        do i_sg = 1, size(prop_vecs, 2)
            do src = 1, nSpatOrbs
                do tgt = 1, nSpatOrbs
                block
                    type(Excite_1_t) :: exc
                    exc = Excite_1_t(2 * src, 2 * tgt)
                    if (src /= tgt &
                            .and. this%indexer%is_allowed(exc, prop_vecs(:, i_sg)) &
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

    end subroutine

    logical pure function symmetry_allowed(exc)
        use SymExcitDataMod, only: SpinOrbSymLabel
        type(Excite_1_t), intent(in) :: exc
        symmetry_allowed = SpinOrbSymLabel(exc%val(1, 1)) == SpinOrbSymLabel(exc%val(2, 1))
    end function

    subroutine PropVec_singles_uniform_finalize(this)
        class(PropVec_UniformExcGenerator_t), intent(inout) :: this
        deallocate(this%allowed_holes, this%indexer)
        if (this%create_lookup) deallocate(lookup_property_indexer)
    end subroutine


    subroutine pickOrbitals_single(this, nI, ilut, csf_i, excitInfo, pgen)
        debug_function_name("pickOrbitals_single")
        class(PropVec_UniformExcGenerator_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
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
        class(PropVec_UniformExcGenerator_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: src_orb
        integer, allocatable :: valid_orbs(:)

        integer(n_int) :: &
            ilut_doubly_occupied(0:this%L_spat_bits), ilut_valid(0:this%L_spat_bits)
        integer :: i_sg

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_prop_vec_idx(csf_I)
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
        class(PropVec_UniformExcGenerator_t), intent(in) :: this
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
        type(PropVec_PCHB_SinglesAlgorithm_t) :: res
        routine_name("alg_from_str")

        select case(to_upper(str))
        case(option_vals%algorithm%UNIF_UNIF%str)
            res = option_vals%algorithm%UNIF_UNIF
        case(option_vals%algorithm%FULL_FULL%str)
            res = option_vals%algorithm%FULL_FULL
        case(option_vals%algorithm%UNIF_FULL%str)
            res = option_vals%algorithm%UNIF_FULL
        case(option_vals%algorithm%UNIF_XNEW%str)
            res = option_vals%algorithm%UNIF_XNEW
        case default
            call stop_all(this_routine, 'Invalid input: '//str)
        end select
    end function

    subroutine init_PropVec_Weighted_t(this, weighting, indexer, use_lookup, create_lookup)
        class(PropVec_Weighted_t), intent(inout) :: this
        type(PropVec_PCHB_SinglesWeighting_t), intent(in) :: weighting
        class(AlsoGUGA_PropertyIndexer_t), intent(in) :: indexer
        logical, intent(in) :: use_lookup, create_lookup
        integer, allocatable :: prop_vecs(:, :)
        routine_name("init_PropVec_Weighted_t")

        integer :: i_sg, src, tgt, nBi
        real(dp), allocatable :: I_weights(:), A_weights(:)

        allocate(this%indexer, source=indexer)
        nBi = this%indexer%n_spin_orbs()
        @:ASSERT(nBi == 2 * nSpatOrbs)
        this%create_lookup = create_lookup
        if (create_lookup) then
            if (allocated(lookup_property_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the prop_vec lookup.')
            else
                write(stdout, *) 'PropVec singles is creating and managing the prop_vec lookup'
                allocate(lookup_property_indexer, source=this%indexer)
            end if
        end if
        this%use_lookup = use_lookup
        if (use_lookup) write(stdout, *) 'PropVec singles is using the prop_vec lookup'
        ! possible prop_vecs
        prop_vecs = this%indexer%get_prop_vecs()

        this%unoccupied = UnoccupiedGetter_t(nBasis)


        call set_get_weight_pointer(weighting, this%get_weight)

        call this%I_sampler%shared_alloc(size(prop_vecs, 2), nBi, 'PCHB')
        call this%A_sampler%shared_alloc([nSpatOrbs, size(prop_vecs, 2)], nBi, 'PCHB')

        allocate(I_weights(nBi), A_weights(nBi), source=0._dp)
        do i_sg = 1, size(prop_vecs, 2)
            I_weights(:) = 0._dp
            do src = 1, nSpatOrbs
                A_weights(:) = 0._dp
                do tgt = 1, nSpatOrbs
                block
                    type(Excite_1_t) :: exc
                    exc = Excite_1_t(2 * src, 2 * tgt)
                    if (src /= tgt &
                            .and. this%indexer%is_allowed(exc, prop_vecs(:, i_sg)) &
                            .and. symmetry_allowed(exc)) then
                        A_weights(2 * tgt) = abs(this%get_weight(src, tgt))
                        A_weights(2 * tgt - 1) = A_weights(2 * tgt)
                    end if
                end block
                end do
                I_weights(2 * src) = sum(A_weights)
                I_weights(2 * src - 1) = I_weights(2 * src)
                call this%A_sampler%setup_entry(src, i_sg, root, A_weights)
            end do
            call this%I_sampler%setup_entry(i_sg, root, I_weights)
        end do

    end subroutine

    subroutine finalize_PropVec_Weighted_t(this)
        class(PropVec_Weighted_t), intent(inout) :: this
        if (allocated(this%indexer)) then
            ! Yes I assume that either all or none are allocated
            call this%I_sampler%finalize()
            call this%A_sampler%finalize()
            nullify(this%get_weight)
            deallocate(this%indexer)
            if (this%create_lookup) deallocate(lookup_property_indexer)
        end if
    end subroutine


    subroutine pickOrbitals_FullFull(this, nI, ilut, csf_i, excitInfo, pgen)
        class(SinglesPropVec_FullFull_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        routine_name("pickOrbitals_PropVec_Weighted_t")

        integer :: i_sg, src, tgt
        real(dp) :: p_src, p_tgt

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_prop_vec_idx(csf_I)
        else
            i_sg = this%indexer%idx_nI(nI)
        end if

        select_particle: block
            real(dp) :: renorm_src
            integer :: spin_src, elec
            renorm_src = sum(this%I_sampler%get_prob(i_sg, nI))
            call this%I_sampler%constrained_sample(i_sg, nI, ilut(0 : nIfD), renorm_src, elec, spin_src, p_src)
            if (spin_src == 0) then
                call make_invalid()
                return
            end if
            @:ASSERT(nI(elec) == spin_src)
            src = gtId(spin_src)
            p_src = csf_I%occ_int(src) * p_src
        end block select_particle

        select_hole: block
            real(dp) :: renorm_tgt
            integer :: spin_tgt, dummy, unoccupied(nBasis - nEl)
            integer(n_int) :: ilut_unoccupied(0 : nIfD)
            call this%unoccupied%get(ilut(0 : nIfD), ilut_unoccupied, unoccupied)
            renorm_tgt = 1._dp - sum(this%A_sampler%get_prob(src, i_sg, nI))
            if (do_direct_calculation(renorm_tgt)) then
                renorm_tgt = sum(this%A_sampler%get_prob(src, i_sg, unoccupied))
            end if

            call this%A_sampler%constrained_sample(&
                 src, i_sg, unoccupied, ilut_unoccupied, renorm_tgt, dummy, spin_tgt, p_tgt)
            if (spin_tgt == 0) then
                call make_invalid()
                return
            end if
            tgt = gtId(spin_tgt)
            p_tgt = (2 - csf_I%occ_int(tgt)) * p_tgt
        end block select_hole
        pgen = p_src * p_tgt

        if (tgt < src) then
            excitInfo = assign_excitinfo_values_single( &
                        gen_type%R, tgt, src, tgt, src)
        else
            excitInfo = assign_excitinfo_values_single( &
                        gen_type%L, tgt, src, src, tgt)
        end if
        excitInfo%valid = checkCompatibility_single(csf_i, excitInfo)
        if (excitInfo%valid) then
            @:ASSERT(pgen .isclose. this%calc_pgen(nI, ilut, csf_i, excitInfo))
        else
            pgen = 0._dp
        end if

        contains

            subroutine make_invalid()
                pgen = 0.0_dp
                tgt = 0
                excitInfo%valid = .false.
            end subroutine
    end subroutine


    function calc_pgen_FullFull(this, nI, ilutI, csf_i, excitInfo) result(pgen)
        class(SinglesPropVec_FullFull_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen
        routine_name("calc_pgen_FullFull")

        real(dp) :: p_src, p_tgt
        integer :: i_sg

        @:ASSERT(1 <= excitInfo%i .and. excitInfo%i <= nSpatOrbs)
        @:ASSERT(1 <= excitInfo%j .and. excitInfo%j <= nSpatOrbs)

        if (excitInfo%i == excitInfo%j &
               .or. csf_i%stepvector(excitInfo%i) == 3 &
               .or. csf_i%stepvector(excitInfo%j) == 0) then
            pgen = 0.0_dp
            return
        end if

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_prop_vec_idx(csf_I)
        else
            i_sg = this%indexer%idx_csf(csf_i)
        end if

        select_particle: block
            real(dp) :: renorm_src
            integer :: spin_src
            spin_src = 2 * excitInfo%j - 1 + ((csf_I%stepvector(excitInfo%j) + 1) .div. 3)
            renorm_src = sum(this%I_sampler%get_prob(i_sg, nI))
            p_src = csf_I%occ_int(excitInfo%j) &
                    * this%I_sampler%constrained_getProb(i_sg, nI, renorm_src, spin_src)
        end block select_particle

        select_hole: block
            real(dp) :: renorm_tgt
            integer :: spin_tgt, unoccupied(nBasis - nEl)
            integer(n_int) :: ilut_unoccupied(0 : nIfD)
            spin_tgt = 2 * excitInfo%i - ((csf_I%stepvector(excitInfo%i) + 1) .div. 3)

            call this%unoccupied%get(ilutI(0 : nIfD), ilut_unoccupied, unoccupied)
            renorm_tgt = 1._dp - sum(this%A_sampler%get_prob(excitInfo%j, i_sg, nI))
            if (do_direct_calculation(renorm_tgt)) then
                renorm_tgt = sum(this%A_sampler%get_prob(excitInfo%j, i_sg, unoccupied))
            end if


            p_tgt = (2 - csf_I%occ_int(excitInfo%i)) &
                    * this%A_sampler%constrained_getProb(excitInfo%j, i_sg, unoccupied, renorm_tgt, spin_tgt)
        end block select_hole
        pgen = p_src * p_tgt
    end function

    subroutine pickOrbitals_UnifFull(this, nI, ilut, csf_i, excitInfo, pgen)
        class(SinglesPropVec_UnifFull_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        routine_name("pickOrbitals_PropVec_Weighted_t")

        integer :: i_sg, src, tgt
        real(dp) :: p_src, p_tgt

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_prop_vec_idx(csf_I)
        else
            i_sg = this%indexer%idx_nI(nI)
        end if

        select_particle: block
            src = gtID(nI(int(genrand_real2_dSFMT() * nEl) + 1))
            p_src = csf_I%occ_real(src) / real(nEl, dp)
        end block select_particle

        select_hole: block
            real(dp) :: renorm_tgt
            integer :: spin_tgt, dummy, unoccupied(nBasis - nEl)
            integer(n_int) :: ilut_unoccupied(0 : nIfD)
            call this%unoccupied%get(ilut(0 : nIfD), ilut_unoccupied, unoccupied)
            renorm_tgt = 1._dp - sum(this%A_sampler%get_prob(src, i_sg, nI))
            if (do_direct_calculation(renorm_tgt)) then
                renorm_tgt = sum(this%A_sampler%get_prob(src, i_sg, unoccupied))
            end if

            call this%A_sampler%constrained_sample(&
                 src, i_sg, unoccupied, ilut_unoccupied, renorm_tgt, dummy, spin_tgt, p_tgt)
            if (spin_tgt == 0) then
                call make_invalid()
                return
            end if
            tgt = gtId(spin_tgt)
            p_tgt = (2 - csf_I%occ_int(tgt)) * p_tgt
        end block select_hole
        pgen = p_src * p_tgt

        if (tgt < src) then
            excitInfo = assign_excitinfo_values_single( &
                        gen_type%R, tgt, src, tgt, src)
        else
            excitInfo = assign_excitinfo_values_single( &
                        gen_type%L, tgt, src, src, tgt)
        end if
        excitInfo%valid = checkCompatibility_single(csf_i, excitInfo)
        if (excitInfo%valid) then
            @:ASSERT(pgen .isclose. this%calc_pgen(nI, ilut, csf_i, excitInfo))
        else
            pgen = 0._dp
        end if

        contains

            subroutine make_invalid()
                pgen = 0.0_dp
                tgt = 0
                excitInfo%valid = .false.
            end subroutine
    end subroutine


    function calc_pgen_UnifFull(this, nI, ilutI, csf_i, excitInfo) result(pgen)
        class(SinglesPropVec_UnifFull_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen
        routine_name("calc_pgen_UnifFull")

        real(dp) :: p_src, p_tgt
        integer :: i_sg

        @:ASSERT(1 <= excitInfo%i .and. excitInfo%i <= nSpatOrbs)
        @:ASSERT(1 <= excitInfo%j .and. excitInfo%j <= nSpatOrbs)

        if (excitInfo%i == excitInfo%j &
               .or. csf_i%stepvector(excitInfo%i) == 3 &
               .or. csf_i%stepvector(excitInfo%j) == 0) then
            pgen = 0.0_dp
            return
        end if

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_prop_vec_idx(csf_I)
        else
            i_sg = this%indexer%idx_csf(csf_i)
        end if

        select_particle: block
            p_src = csf_I%occ_real(excitInfo%j) / real(nEl, dp)
        end block select_particle

        select_hole: block
            real(dp) :: renorm_tgt
            integer :: spin_tgt, unoccupied(nBasis - nEl)
            integer(n_int) :: ilut_unoccupied(0 : nIfD)
            spin_tgt = 2 * excitInfo%i - ((csf_I%stepvector(excitInfo%i) + 1) .div. 3)

            call this%unoccupied%get(ilutI(0 : nIfD), ilut_unoccupied, unoccupied)
            renorm_tgt = 1._dp - sum(this%A_sampler%get_prob(excitInfo%j, i_sg, nI))
            if (do_direct_calculation(renorm_tgt)) then
                renorm_tgt = sum(this%A_sampler%get_prob(excitInfo%j, i_sg, unoccupied))
            end if


            p_tgt = (2 - csf_I%occ_int(excitInfo%i)) &
                    * this%A_sampler%constrained_getProb(excitInfo%j, i_sg, unoccupied, renorm_tgt, spin_tgt)
        end block select_hole
        pgen = p_src * p_tgt

    end function



    subroutine init_UnifXnew_t(this, weighting, indexer, use_lookup, create_lookup, scale_down)
        class(UnifXnew_t), intent(inout) :: this
        type(PropVec_PCHB_SinglesWeighting_t), intent(in) :: weighting
        class(AlsoGUGA_PropertyIndexer_t), intent(in) :: indexer
        logical, intent(in) :: use_lookup, create_lookup
        real(dp), intent(in) :: scale_down
            !! Scale down the "exchange" single excitations
        integer, allocatable :: prop_vecs(:, :)
        routine_name("init_UnifXnew_t")

        integer :: i_sg, src, tgt, nBi
        real(dp), allocatable :: spat_w_A(:), A_weights(:)

        allocate(this%indexer, source=indexer)
        nBi = this%indexer%n_spin_orbs()
        @:ASSERT(nBi == 2 * nSpatOrbs)
        this%create_lookup = create_lookup
        if (create_lookup) then
            if (allocated(lookup_property_indexer)) then
                call stop_all(this_routine, 'Someone else is already managing the prop_vec lookup.')
            else
                lookup_property_indexer = this%indexer
            end if
        end if
        this%use_lookup = use_lookup

        prop_vecs = this%indexer%get_prop_vecs()

        this%unoccupied = UnoccupiedGetter_t(nBasis)

        call set_get_weight_pointer(weighting, this%get_weight)

        call this%A_sampler%shared_alloc([nSpatOrbs, 3, size(prop_vecs, 2)], nBi, 'PCHB')

        allocate(spat_w_A(nSpatOrbs), A_weights(nBi), source=0._dp)
        do i_sg = 1, size(prop_vecs, 2)
            do src = 1, nSpatOrbs
                A_weights(:) = 0._dp
                spat_w_A(:) = 0._dp
                do tgt = 1, nSpatOrbs
                block
                    type(Excite_1_t) :: exc
                    exc = Excite_1_t(2 * src, 2 * tgt)
                    if (src /= tgt &
                            .and. this%indexer%is_allowed(exc, prop_vecs(:, i_sg)) &
                            .and. symmetry_allowed(exc)) then
                        spat_w_A(tgt) = abs(this%get_weight(src, tgt))
                    end if
                end block
                end do

                start_from_3: block
                    A_weights(1::2) = spat_w_A(:)
                    A_weights(2::2) = spat_w_A(:)
                    @:ASSERT(near_zero(A_weights(2 * src - 1)) .and. near_zero(A_weights(2 * src)))
                    call this%A_sampler%setup_entry(src, 3, i_sg, root, A_weights)
                end block start_from_3

                start_from_1: block
                    A_weights(1::2) = spat_w_A(:)
                    A_weights(2::2) = spat_w_A(:) * scale_down
                    @:ASSERT(near_zero(A_weights(2 * src - 1)) .and. near_zero(A_weights(2 * src)))
                    if (src > 1) A_weights(2 * (src - 1)) = 0._dp
                    if (src < nSpatOrbs) A_weights(2 * (src + 1)) = 0._dp
                    call this%A_sampler%setup_entry(src, 1, i_sg, root, A_weights)
                end block start_from_1

                start_from_2: block
                    A_weights(1::2) = spat_w_A(:) * scale_down
                    A_weights(2::2) = spat_w_A(:)
                    @:ASSERT(near_zero(A_weights(2 * src - 1)) .and. near_zero(A_weights(2 * src)))
                    if (src > 1) A_weights(2 * (src - 1) - 1) = 0._dp
                    if (src < nSpatOrbs) A_weights(2 * (src + 1) - 1) = 0._dp
                    call this%A_sampler%setup_entry(src, 2, i_sg, root, A_weights)
                end block start_from_2
            end do
        end do
    end subroutine

    subroutine finalize_UnifXnew_t(this)
        class(UnifXnew_t), intent(inout) :: this
        call this%A_sampler%finalize()
        nullify(this%get_weight)
        deallocate(this%indexer)
        if (this%create_lookup) deallocate(lookup_property_indexer)
    end subroutine

    subroutine pickOrbitals_UnifXnew_t(this, nI, ilut, csf_i, excitInfo, pgen)
        class(UnifXnew_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        routine_name("pickOrbitals_PropVec_Weighted_t")

        integer :: i_sg, src, tgt
        real(dp) :: p_src, p_tgt
        real(dp) :: renorm_tgt
        integer :: spin_tgt, dummy, unoccupied(nBasis - nEl)
        integer(n_int) :: ilut_unoccupied(0 : nIfD)

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_prop_vec_idx(csf_I)
        else
            i_sg = this%indexer%idx_nI(nI)
        end if

        select_particle: block
            src = gtID(nI(int(genrand_real2_dSFMT() * nEl) + 1))
            p_src = csf_I%occ_real(src) / real(nEl, dp)
        end block select_particle

        select_hole: block

            call this%unoccupied%get(ilut(0 : nIfD), ilut_unoccupied, unoccupied)

            renorm_tgt = 1._dp - sum(this%A_sampler%get_prob(src, csf_i%stepvector(src), i_sg, nI))
            if (do_direct_calculation(renorm_tgt)) then
                renorm_tgt = sum(this%A_sampler%get_prob(src, csf_i%stepvector(src), i_sg, unoccupied))
            end if

            call this%A_sampler%constrained_sample(&
                 src, csf_i%stepvector(src), i_sg, unoccupied, ilut_unoccupied, &
                 renorm_tgt, dummy, spin_tgt, p_tgt)
            if (spin_tgt == 0) then
                call make_invalid()
                return
            end if
            tgt = gtId(spin_tgt)
            if (csf_I%occ_int(tgt) == 0) then
                if (mod(spin_tgt, 2) == 1) then
                    p_tgt = p_tgt &
                            + this%A_sampler%constrained_getProb(&
                                    src, csf_i%stepvector(src), i_sg, &
                                    unoccupied, renorm_tgt, 2 * tgt)
                else
                    p_tgt = p_tgt &
                            + this%A_sampler%constrained_getProb(&
                                    src, csf_i%stepvector(src), i_sg, &
                                    unoccupied, renorm_tgt, 2 * tgt - 1)
                end if
            end if
        end block select_hole
        pgen = p_src * p_tgt

        if (tgt < src) then
            excitInfo = assign_excitinfo_values_single( &
                        gen_type%R, tgt, src, tgt, src)
        else
            excitInfo = assign_excitinfo_values_single( &
                        gen_type%L, tgt, src, src, tgt)
        end if
        excitInfo%valid = checkCompatibility_single(csf_i, excitInfo)
        if (excitInfo%valid) then
            @:ASSERT(pgen .isclose. this%calc_pgen(nI, ilut, csf_i, excitInfo))
        else
            pgen = 0._dp
        end if

        contains

            subroutine make_invalid()
                pgen = 0.0_dp
                tgt = 0
                excitInfo%valid = .false.
            end subroutine
    end subroutine


    function calc_pgen_UnifXnew_t(this, nI, ilutI, csf_i, excitInfo) result(pgen)
        class(UnifXnew_t), intent(in) :: this
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:GugaBits%len_tot)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen
        routine_name("calc_pgen_UnifXnew_t")

        real(dp) :: p_src, p_tgt
        integer :: i_sg

        associate(tgt => excitInfo%i, src => excitInfo%j)

        @:ASSERT(1 <= src .and. src <= nSpatOrbs)
        @:ASSERT(1 <= tgt .and. tgt <= nSpatOrbs)


        if (src == tgt &
               .or. csf_i%occ_int(tgt) == 2 &
               .or. csf_i%occ_int(src) == 0) then
            pgen = 0.0_dp
            return
        end if

        if (this%use_lookup) then
            i_sg = this%indexer%lookup_prop_vec_idx(csf_I)
        else
            i_sg = this%indexer%idx_csf(csf_i)
        end if

        select_particle: block
            p_src = csf_I%occ_real(src) / real(nEl, dp)
        end block select_particle

        select_hole: block
            real(dp) :: renorm_tgt
            integer :: unoccupied(nBasis - nEl)
            integer(n_int) :: ilut_unoccupied(0 : nIfD)

            call this%unoccupied%get(ilutI(0 : nIfD), ilut_unoccupied, unoccupied)

            renorm_tgt = 1._dp - sum(this%A_sampler%get_prob(src, csf_i%stepvector(src), i_sg, nI))
            if (do_direct_calculation(renorm_tgt)) then
                renorm_tgt = sum(this%A_sampler%get_prob(src, csf_i%stepvector(src), i_sg, unoccupied))
            end if

            if (csf_i%occ_int(tgt) == 0) then
                p_tgt = this%A_sampler%constrained_getProb(src, csf_i%stepvector(src), i_sg, unoccupied, renorm_tgt, 2 * tgt - 1) &
                        + this%A_sampler%constrained_getProb(src, csf_i%stepvector(src), i_sg, unoccupied, renorm_tgt, 2 * tgt)
            else if (csf_i%stepvector(tgt) == 1) then
                p_tgt = this%A_sampler%constrained_getProb(src, csf_i%stepvector(src), i_sg, unoccupied, renorm_tgt, 2 * tgt)
            else if (csf_i%stepvector(tgt) == 2) then
                p_tgt = this%A_sampler%constrained_getProb(src, csf_i%stepvector(src), i_sg, unoccupied, renorm_tgt, 2 * tgt - 1)
            else
                call stop_all(this_routine, "should not be here.")
            end if

        end block select_hole
        pgen = p_src * p_tgt
        end associate

    end function

end module
