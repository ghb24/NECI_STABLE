#include "macros.h"
#:include "../macros.fpph"
#:include "../algorithms.fpph"

module pq_property_index

    use constants, only: int64, n_int
    use SystemData, only: nEl
    use util_mod, only: choose_i64, cumsum, binary_search_first_ge, custom_findloc, stop_all
    use sets_mod, only: subset
    use bit_rep_data, only: nIfD
    use hash, only: hash_table_lookup
    use guga_bitRepOps, only: CSF_Info_t
    use growing_buffers, only: buffer_int_2D_t, buffer_int64_1D_t
    use composition_utils, only: n_compositions, get_compositions, composition_idx, &
                                 composition_from_idx, next_composition
    use excitation_types, only: Excite_1_t, Excite_2_t, Excite_3_t
    use guga_bitRepOps, only: CSF_Info_t
    use guga_data, only: DistinctDouble_t
    use property_vector_index, only: PropertyIndexer_t, AlsoGUGA_PropertyIndexer_t
    use gasci_supergroup_index, only: GAS_SuperGroupIndexer_t
    use gasci, only: GASSpec_t, LocalGASSpec_t, CumulGASSpec_t, FlexibleGASSpec_t, CAS_spec
    use orb_idx_mod, only: SpinProj_t, calc_spin_raw

    better_implicit_none

    private
    public :: PQ_indexer_t, P_indexer, total_indexer

    type, extends(AlsoGUGA_PropertyIndexer_t) :: PQ_indexer_t

        ! This type extends the PropertyIndexer to perfrom PQ downfolding
        ! In order to do this we require a total space, P space and
        ! a vector contained_in_P that keeps track if a property vector is 
        ! in the P space or not.

        class(GAS_SuperGroupIndexer_t), allocatable :: total_space, P_space
        logical, allocatable :: contained_in_P(:)

    contains

        procedure :: get_contained_in_P
        procedure :: is_allowed_single
        procedure :: is_allowed_double
        procedure :: is_allowed_triple
        procedure :: is_allowed_configurations
        procedure :: is_allowed_distinct_double

        procedure :: excite_single
        procedure :: excite_double
        procedure :: excite_triple
        generic :: excite => excite_single, excite_double, excite_triple

        procedure :: gen_all_excits => PQ_gen_all_excits

        procedure, public :: to_prop_vec_csf => PQ_to_prop_vec_csf_t
        procedure, public :: to_prop_vec_ni => PQ_to_prop_vec_ni_t

        procedure, public :: write_to => PQ_write_to
    end type

    interface PQ_indexer_t
        module procedure construct_PQ_indexer_t
    end interface

    class(GAS_SuperGroupIndexer_t), allocatable :: P_indexer
    class(PropertyIndexer_t), allocatable :: total_indexer

contains

    pure function get_contained_in_P(this) result(contained_in_P)

        ! This function asserts for all composition indices in the total space
        ! wether they are contained in the P space or not.

        use sets_mod, only: operator(.in.)
        class(PQ_indexer_t), intent(in) :: this
        integer :: i
        logical :: contained_in_P(size(this%total_space%allowed_composition_indices))

        do i = 1, size(contained_in_P)
            contained_in_P(i) = this%total_space%allowed_composition_indices(i) .in. this%P_space%allowed_composition_indices
        end do
        
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

    pure function excite_triple(this, exc, prop_vec) result(res)
        class(PQ_indexer_t), intent(in) :: this
        type(Excite_3_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)
        integer :: res(size(prop_vec))

        res = this%total_space%excite(exc, prop_vec)

    end function

    logical pure function is_allowed_single(this, exc, prop_vec)

        ! This functions checks if a single excitation is allowed
        ! under GAS and PQ constraints

        class(PQ_indexer_t), intent(in) :: this
        type(Excite_1_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)
        logical :: valid_excitation
        integer :: exc_prop_vec(size(prop_vec))

        valid_excitation = this%total_space%is_allowed(exc, prop_vec)

        if (valid_excitation) then

            exc_prop_vec = this%excite(exc, prop_vec)

            is_allowed_single = (&
                                this%contained_in_P(this%idx_prop_vec(prop_vec)) &
                                .or. this%contained_in_P(this%idx_prop_vec(exc_prop_vec)) &
                                )
        else

            is_allowed_single = .false.

        end if

    end function

    logical pure function is_allowed_double(this, exc, prop_vec)

        ! This functions checks if a double excitation is allowed
        ! under GAS and PQ constraints

        class(PQ_indexer_t), intent(in) :: this
        type(Excite_2_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)
        logical :: valid_excitation
        integer :: exc_prop_vec(size(prop_vec))

        valid_excitation = this%total_space%is_allowed(exc, prop_vec)

        if (valid_excitation) then

            exc_prop_vec = this%excite(exc, prop_vec)

            is_allowed_double = this%contained_in_P(&
                                                    this%idx_prop_vec(prop_vec)) &
                                                    .or. this%contained_in_P(this%idx_prop_vec(exc_prop_vec) &
                                                    )

        else

            is_allowed_double = .false.

        end if

    end function

    logical pure function is_allowed_triple(this, exc, prop_vec)

        ! This functions checks if a triple excitation is allowed
        ! under GAS and PQ constraints

        class(PQ_indexer_t), intent(in) :: this
        type(Excite_3_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)
        logical :: valid_excitation
        integer :: exc_prop_vec(size(prop_vec))

        valid_excitation = this%total_space%is_allowed(exc, prop_vec)

        if (valid_excitation) then

            exc_prop_vec = this%excite(exc, prop_vec)

            is_allowed_triple = this%contained_in_P( &
                                                    this%idx_prop_vec(prop_vec)) &
                                                    .or. this%contained_in_P(this%idx_prop_vec(exc_prop_vec) &
                                                    )
        else

            is_allowed_triple = .false.

        end if

    end function

    logical pure function is_allowed_distinct_double(this, exc, prop_vec)

        ! This functions checks if a  distinct double excitation is allowed
        ! under GAS and PQ constraints

        class(PQ_indexer_t), intent(in) :: this
        type(DistinctDouble_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)
        integer ::  i, j, a, b

        call exc%idx%extract_ijkl(a, i, b, j)

        is_allowed_distinct_double = this%is_allowed_double(Excite_2_t(2*i, 2*a, 2*j, 2*b), prop_vec)

    end function

    logical pure function is_allowed_configurations(this, nI, nJ)
        !! Check if an excitation nI -> nJ is allowed.
        !!
        !! Is called once at initialization, so it does not have to be super fast.
        class(PQ_indexer_t), intent(in) :: this
        integer, intent(in) :: nI(nEl), nJ(nEl)
        is_allowed_configurations = &
            this%total_space%GAS_spec%contains_conf(nI) &
            .and. this%total_space%GAS_spec%contains_conf(nJ)
    end function

    subroutine PQ_gen_all_excits(this, nI, n_excits, det_list, ic)
        class(PQ_indexer_t), intent(in) :: this
        integer, intent(in) :: nI(:)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:, :)
        integer, intent(in), optional :: ic

        call this%total_space%gen_all_excits(nI, n_excits, det_list, ic)

    end subroutine

    pure function PQ_to_prop_vec_csf_t(this, csf_i) result(prop_vec)
        class(PQ_indexer_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
#ifndef IFORT_
        ! Please intel compilers...
        integer :: prop_vec(this%total_space%prop_vec_dim())
#else
        integer, allocatable :: prop_vec(:)
#endif
        prop_vec = this%total_space%to_prop_vec(csf_i)

    end function

    pure function PQ_to_prop_vec_nI_t(this, nI) result(prop_vec)
        class(PQ_indexer_t), intent(in) :: this
        integer, intent(in) :: nI(:)
#ifndef IFORT_
        ! Please intel compilers...
        integer :: prop_vec(this%total_space%prop_vec_dim())
#else
        integer, allocatable :: prop_vec(:)
#endif

        prop_vec = this%total_space%to_prop_vec(nI)

    end function

    pure function construct_PQ_indexer_t(P_indexer, total_indexer) result(idxer)

        ! This function constructs the PQ indexer

        use gasci_supergroup_index, only: get_allowed_composition_indices
        class(GAS_SuperGroupIndexer_t), intent(in) :: P_indexer, total_indexer
        type(PQ_indexer_t) :: idxer
        character(*), parameter :: this_routine = 'construct_PQ_indexer_t'

        ! Pass the definition of the total and P space to the indexer

        idxer%total_space = total_indexer
        idxer%P_space = P_indexer

        @:pure_ASSERT(allocated(idxer%total_space))
        @:pure_ASSERT(allocated(idxer%P_space))

        ! Check if the total space is connected
        ! Otherwise PQ downfolding is not possible

        if (total_indexer%GAS_spec%is_connected()) then

            idxer%allowed_composition_indices = idxer%total_space%allowed_composition_indices

        else

            call stop_all(this_routine, "The total GAS Space cannot be disconnected")

        end if

        ! idxer%P_space%allowed_composition_indices = [-1] for a disconnected P_specification
        ! Overwrite it with the composition_idx of the single property_vector
        ! Otherwise P_space can't be subset of total_space

        if (.not. P_indexer%GAS_spec%is_connected()) then

            idxer%P_space%allowed_composition_indices = get_allowed_composition_indices(P_indexer%GAS_spec)

        end if

        ! Assert if P is subset of total

        @:pure_ASSERT(subset(idxer%P_space%allowed_composition_indices,&
                             idxer%total_space%allowed_composition_indices))

        ! Pass properties of total space to PQ_indexer

        idxer%prop_vec_sum_ = idxer%total_space%prop_vec_sum()
        idxer%n_spin_orbs_ = idxer%total_space%n_spin_orbs()
        idxer%prop_vec_dim_ = idxer%total_space%prop_vec_dim()
        idxer%contained_in_P = idxer%get_contained_in_P()

        @:pure_ASSERT(idxer%prop_vec_sum_ == total_indexer%GAS_spec%nEl())
        @:pure_ASSERT(idxer%n_spin_orbs_ == total_indexer%GAS_spec%n_spin_orbs())
        @:pure_ASSERT(idxer%prop_vec_dim_ == total_indexer%GAS_spec%nGAS())

        ! Assert if contained_in_P has correct dimensions and if number of indices in P is consistent

        @:pure_ASSERT(size(idxer%contained_in_P) == size(idxer%total_space%allowed_composition_indices))
        @:pure_ASSERT(count(idxer%contained_in_P) == size(idxer%P_space%allowed_composition_indices))

    end function

    subroutine PQ_write_to(this, iunit)
        class(PQ_indexer_t), intent(in) :: this
        integer, intent(in) :: iunit

        write(iunit, '(/A)') 'PQ downfolding will be performed for the following GAS spaces:'

        write(iunit, '(/A)') 'Specification of the total space:'
        call this%total_space%write_to(iunit)

        write(iunit, '(/A)') 'Specification of the P space:'
        call this%P_space%write_to(iunit)

        write(iunit, '(/A20, 1x, A20, 1x)') 'n_prop_vecs (total)', 'n_prop_vecs (P)'
        write(iunit, '(A)') '------------------------------------------'
        write(iunit, '(I20, 1x, I20, 1x)') size(this%total_space%allowed_composition_indices), &
                                           size(this%P_space%allowed_composition_indices)
        write(iunit, '(A)') '------------------------------------------'

    end subroutine

end module pq_property_index
