#include "macros.h"
#:include "../macros.fpph"
#:include "../algorithms.fpph"

module pq_property_index

    use constants, only: int64, n_int
    use util_mod, only: choose_i64, cumsum, binary_search_first_ge, custom_findloc, stop_all
    use bit_rep_data, only: nIfD
    use hash, only: hash_table_lookup
    use guga_bitRepOps, only: CSF_Info_t
    use growing_buffers, only: buffer_int_2D_t, buffer_int64_1D_t
    use composition_utils, only: n_compositions, get_compositions, composition_idx, &
                                 composition_from_idx, next_composition
    use excitation_types, only: Excite_1_t, Excite_2_t
    use guga_bitRepOps, only: CSF_Info_t
    use guga_data, only: DistinctDouble_t
    use property_vector_index, only: PropertyIndexer_t, AlsoGUGA_PropertyIndexer_t
    use gasci_supergroup_index, only: GAS_SuperGroupIndexer_t
    use gasci, only: GASSpec_t, LocalGASSpec_t, CumulGASSpec_t, FlexibleGASSpec_t, CAS_spec
    better_implicit_none

    private
    public :: PQ_indexer_t

    type, extends(AlsoGUGA_PropertyIndexer_t) :: PQ_indexer_t
        class(GAS_SuperGroupIndexer_t), allocatable :: total_space, P_space
        logical, allocatable :: contained_in_P(:)

    contains
        procedure :: get_contained_in_P
        procedure :: is_allowed_single
        procedure :: is_allowed_double
        procedure :: is_allowed_distinct_double

        procedure :: excite_single
        procedure :: excite_double
        generic :: excite => excite_single, excite_double

        procedure :: gen_all_excits => PQ_gen_all_excits

        procedure :: to_prop_vec_csf => PQ_to_prop_vec_csf
        procedure :: to_prop_vec_ni => PQ_to_prop_vec_ni
    end type

    interface PQ_indexer_t
        module procedure construct_PQ_indexer_t
    end interface
contains

    pure function get_contained_in_P(this) result(contained_in_P)
        class(PQ_indexer_t), intent(in) :: this
        integer :: i
        logical :: contained_in_P(this%total_space%n_prop_vecs())

        !TODO Loop through entries of total space and check if they are contained in P
        !Currently set to true to test constructor (everything is P space -> normal GAS)
        do i = 1, size(contained_in_P)
            contained_in_P(i) = .true.
        end do
        !Just compare composition indices?
    end function

    pure function excite_single(this, exc, prop_vec) result(res)
        class(PQ_indexer_t), intent(in) :: this
        type(Excite_1_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)
        integer :: res(size(prop_vec))

        res = this%total_space%excite(exc, prop_vec)
    end function

    pure function excite_double(this, exc, prop_vec) result(res)
        class(PQ_indexer_t), intent(in) :: this
        type(Excite_2_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)
        integer :: res(size(prop_vec))

        res = this%total_space%excite(exc, prop_vec)
    end function

    logical pure function is_allowed_single(this, exc, prop_vec)
        class(PQ_indexer_t), intent(in) :: this
        type(Excite_1_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)
        logical :: valid_excitation

        !We check if an excitation is GAS (total space) valid, otherwise there is no entry in contained_in_P
        valid_excitation = this%total_space%is_allowed(exc, prop_vec)

        if (valid_excitation) then
            is_allowed_single = (&
                                 this%contained_in_P(this%idx_prop_vec(prop_vec)) .or. &
                                 this%contained_in_P(this%idx_prop_vec(this%excite(exc,prop_vec))) &
                                )
        else
            is_allowed_single = .false.
        end if
    end function

    logical pure function is_allowed_double(this, exc, prop_vec)
        class(PQ_indexer_t), intent(in) :: this
        type(Excite_2_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)
        logical :: valid_excitation

        !We check if an excitation is GAS (total space) valid, otherwise there is no entry in contained_in_P
        valid_excitation = this%total_space%is_allowed(exc, prop_vec)

        if (valid_excitation) then
            is_allowed_double = (&
                                 this%contained_in_P(this%idx_prop_vec(prop_vec)) .or. &
                                 this%contained_in_P(this%idx_prop_vec(this%excite(exc,prop_vec))) &
                                )
        else
            is_allowed_double = .false.
        end if
    end function

    logical pure function is_allowed_distinct_double(this, exc, prop_vec)
        class(PQ_indexer_t), intent(in) :: this
        type(DistinctDouble_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)
        logical :: valid_excitation
        integer :: i, j, a, b
        call exc%idx%extract_ijkl(a, i, b, j)

        !We check if an excitation is GAS (total space) valid, otherwise there is no entry in contained_in_P
        valid_excitation = this%total_space%is_allowed(exc, prop_vec)

        if (valid_excitation) then
            is_allowed_distinct_double = (&
                                 this%contained_in_P(this%idx_prop_vec(prop_vec)) .or. &
                                 this%contained_in_P(this%idx_prop_vec(this%excite(Excite_2_t(2*i, 2*a, 2*j, 2*b),prop_vec))) &
                                )
        else
            is_allowed_distinct_double = .false.
        end if
    end function

    subroutine PQ_gen_all_excits(this, nI, n_excits, det_list, ic)
        class(PQ_indexer_t), intent(in) :: this
        integer, intent(in) :: nI(:)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:, :)
        integer, intent(in), optional :: ic

        call this%total_space%gen_all_excits(nI, n_excits, det_list, ic)
    end subroutine

    pure function PQ_to_prop_vec_csf(this, csf_i) result(supergroup)
        class(PQ_indexer_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        !integer :: supergroup(total_space%prop_vec_dim())
        integer :: supergroup(this%total_space%prop_vec_dim())

        supergroup = this%total_space%to_prop_vec(csf_i)

    end function

    pure function PQ_to_prop_vec_ni(this, nI) result(supergroup)
        class(PQ_indexer_t), intent(in) :: this
        integer, intent(in) :: nI(:)
        integer :: supergroup(this%total_space%prop_vec_dim())

        supergroup = this%total_space%to_prop_vec(nI)

    end function

    pure function construct_PQ_indexer_t(P_spec, GAS_spec) result(idxer)
        class(GASSpec_t), intent(in) :: P_spec, GAS_spec
        type(PQ_indexer_t) :: idxer
        integer(int64), allocatable :: allowed_indices(:)
        character(*), parameter :: this_routine = 'construct_PQ_indexer_t'

        @:pure_ASSERT(GAS_spec%is_valid())
        @:pure_ASSERT(P_spec%is_valid())

        idxer%total_space = GAS_SuperGroupIndexer_t(GAS_spec)
        idxer%P_space = GAS_SuperGroupIndexer_t(P_spec)

        @:pure_ASSERT(allocated(idxer%total_space))
        @:pure_ASSERT(allocated(idxer%P_space))

        if (GAS_spec%is_connected()) then
            allowed_indices = idxer%total_space%allowed_composition_indices
        else
            call stop_all(this_routine, "GAS Space cannot be disconnected")
        end if

        idxer%prop_vec_sum_ = idxer%total_space%prop_vec_sum()
        idxer%n_spin_orbs_ = idxer%total_space%n_spin_orbs()
        idxer%prop_vec_dim_ = idxer%total_space%prop_vec_dim()
        idxer%contained_in_P = get_contained_in_P(idxer)
        allocate (idxer%allowed_composition_indices, source=allowed_indices)

        @:pure_ASSERT(idxer%prop_vec_sum_ == GAS_spec%nEl())
        @:pure_ASSERT(idxer%n_spin_orbs_ == GAS_spec%n_spin_orbs())
        @:pure_ASSERT(idxer%prop_vec_dim_ == GAS_spec%nGAS())
        @:pure_ASSERT(size(idxer%contained_in_P) == size(idxer%total_space%allowed_composition_indices))
    end function
end module pq_property_index
