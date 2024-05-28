#include "macros.h"
#:include "../macros.fpph"
#:include "../algorithms.fpph"


module property_vector_index
    !! This module generalizes the "supergroup" concept of GAS to arbitrary property vectors.
    !! This allows more general constraints such as Löwdin's perturbation theory via PQ space folding
    !!  or constraints that use e.g. the spin projection.
    !!
    !! A property vector is a vector of integers where each element encodes a property of a sub-space.
    !! In traditional GAS the property is just the particle number, while in general these can be e.g.
    !! `[N_1, Sz_1, N_2, Sz_2]`
    !! where `N_i` is the particle number in the i-th sub-space and `Sz_i` is the spin projection
    !! in the i-th sub-space (encoded as integer).
    !! The only assumed constraint of the property vectors is that the sum over them is preserved.
    !!
    !! The algorithm follows the ideas of stochastic GAS (https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00936),
    !! i.e. the property vectors are indexed by calculating the corresponding composition index.
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
    use gasci_supergroup_index, only: GAS_SuperGroupIndexer_t
    better_implicit_none

    private
    public :: PropertyIndexer_t, AlsoGUGA_PropertyIndexer_t, PQ_indexer_t,&
         indexer, lookup_property_indexer, &
         total_indexer, lookup_property_total_indexer, P_indexer, lookup_property_P_indexer

    type, abstract :: PropertyIndexer_t
        integer(int64), allocatable :: allowed_composition_indices(:)
        integer :: n_spin_orbs_
            !! The number of spin orbitals of underlying determinants.
        integer :: prop_vec_sum_
            !! The sum of the property vectors.
            !! In the case of GAS this is the particle number.
        integer :: prop_vec_dim_
            !! The dimension of the property vector
            !!   (corresponds to number of GAS spaces for pure GAS)
    contains
#ifndef IFORT_
        ! Stupid Ifort has problems with overwriting private methods.
        ! https://gitlab.mpcdf.mpg.de/elpa/elpa/-/issues/54
        private
#endif
        procedure, public :: idx_prop_vec
        procedure, public :: idx_nI => get_idx_prop_vec_det
        generic, public :: lookup_prop_vec_idx => lookup_prop_vec_idx_nI
        procedure :: lookup_prop_vec_idx_nI
        procedure, public :: n_prop_vecs => get_n_prop_vecs
        procedure, public :: get_prop_vecs
        procedure(is_allowed_single_t), deferred :: is_allowed_single
        procedure(is_allowed_double_t), deferred :: is_allowed_double
        generic, public :: is_allowed => is_allowed_single, is_allowed_double

        procedure(idxer_gen_all_excits_t), deferred, public :: gen_all_excits

        procedure(to_prop_vec_nI_t), deferred :: to_prop_vec_nI
        generic, public :: to_prop_vec => to_prop_vec_nI

        procedure, public :: prop_vec_dim
        procedure, public :: n_spin_orbs
        procedure, public :: prop_vec_sum
    end type

    type, abstract, extends(PropertyIndexer_t) :: AlsoGUGA_PropertyIndexer_t
    contains
#ifndef IFORT_
        ! Stupid Ifort has problems with overwriting private methods.
        ! https://gitlab.mpcdf.mpg.de/elpa/elpa/-/issues/54
        private
#endif
        procedure, public :: idx_csf => get_prop_vec_idx_csf
        generic, public :: lookup_prop_vec_idx => lookup_prop_vec_idx_CSF
        procedure :: lookup_prop_vec_idx_CSF
        procedure(is_allowed_distinct_double_t), deferred :: is_allowed_distinct_double
        generic, public :: is_allowed => is_allowed_distinct_double

        procedure(to_prop_vec_csf_t), deferred :: to_prop_vec_csf
        generic, public :: to_prop_vec => to_prop_vec_csf
    end type

    abstract interface

        logical pure function is_allowed_single_t(this, exc, prop_vec)
            import :: PropertyIndexer_t, Excite_1_t
            !! Check if a single excitation is allowed.
            !!
            !! Is called once at initialization, so it does not have to be super fast.
            class(PropertyIndexer_t), intent(in) :: this
            type(Excite_1_t), intent(in) :: exc
            integer, intent(in) :: prop_vec(:)
        end function

        logical pure function is_allowed_double_t(this, exc, prop_vec)
            import :: PropertyIndexer_t, Excite_2_t
            !! Check if a double excitation is allowed.
            !!
            !! Is called once at initialization, so it does not have to be super fast.
            class(PropertyIndexer_t), intent(in) :: this
            type(Excite_2_t), intent(in) :: exc
            integer, intent(in) :: prop_vec(:)
        end function

        logical pure function is_allowed_distinct_double_t(this, exc, prop_vec)
            import :: AlsoGUGA_PropertyIndexer_t, DistinctDouble_t
            class(AlsoGUGA_PropertyIndexer_t), intent(in) :: this
            type(DistinctDouble_t), intent(in) :: exc
            integer, intent(in) :: prop_vec(:)
        end function

        subroutine idxer_gen_all_excits_t(this, nI, n_excits, det_list, ic)
            import :: PropertyIndexer_t, n_int
            !!  Get all excitated determinants from
            !!  det_I that are allowed under GAS constraints.
            class(PropertyIndexer_t), intent(in) :: this
            integer, intent(in) :: nI(:)
                !!  Starting determinant
            integer, intent(out) :: n_excits
                !!  Number of determinants
            integer(n_int), intent(out), allocatable :: det_list(:,:)
                !!  Allocatable array of determinants in ilut format
            integer, intent(in), optional :: ic
                !!  Optional input for excitation level (ic=1 => singles, ic=2 => doubles)
                !!      If ommited generate all.
        end subroutine

        pure function to_prop_vec_nI_t(this, nI) result(prop_vec)
            import :: PropertyIndexer_t
            class(PropertyIndexer_t), intent(in) :: this
            integer, intent(in) :: nI(:)
#ifndef IFORT_
            ! Stupid Ifort has problems with defining abstract interfaces
            ! using stack arrays like this.
            ! The abstract interfaces are necessary for defining the deferred method.
            integer :: prop_vec(this%prop_vec_dim())
#else
            integer, allocatable :: prop_vec(:)
#endif
        end function

        pure function to_prop_vec_csf_t(this, csf_i) result(prop_vec)
            import :: AlsoGUGA_PropertyIndexer_t, CSF_Info_t
            class(AlsoGUGA_PropertyIndexer_t), intent(in) :: this
            type(CSF_Info_t), intent(in) :: csf_i
#ifndef IFORT_
            ! Stupid Ifort has problems with defining abstract interfaces
            ! using stack arrays like this.
            ! The abstract interfaces are necessary for defining the deferred method.
            integer :: prop_vec(this%prop_vec_dim())
#else
            integer, allocatable :: prop_vec(:)
#endif
        end function

    end interface

    type, extends(AlsoGUGA_PropertyIndexer_t) :: PQ_indexer_t
        class(GAS_SuperGroupIndexer_t), allocatable :: total_space, P_space
        logical, allocatable :: contained_in_P(:)

    contains
        procedure :: get_contained_in_P
        procedure :: is_allowed_single => PQ_is_allowed_single
        procedure :: is_allowed_double => PQ_is_allowed_double
        procedure :: is_allowed_distinct_double => PQ_is_allowed_distinct_double
        generic :: is_allowed => is_allowed_single, is_allowed_double, is_allowed_distinct_double

        procedure :: excite_single => PQ_excite_single
        procedure :: excite_double => PQ_excite_double
        generic :: excite => excite_single, excite_double

        procedure :: gen_all_excits => PQ_gen_all_excits 

        procedure :: to_prop_vec_csf => PQ_to_prop_vec_csf
        procedure :: to_prop_vec_ni => PQ_to_prop_vec_ni
    end type

    !interface
    !    module procedure construct_PQ_indexer
    !end interface


    class(PropertyIndexer_t), allocatable :: indexer, lookup_property_indexer
    
    !If we want to do the PQ Truncation we need an indexer for the P space, one for the full basis
    !And the indexer (see above) that is called by PCHB who contains both of them
    class(PropertyIndexer_t), allocatable :: total_indexer, lookup_property_total_indexer
    class(PropertyIndexer_t), allocatable :: P_indexer, lookup_property_P_indexer
contains

    pure function lookup_prop_vec_idx_nI(this, idet, nI) result(idx)
        !!  Use a precomputed prop_vec index from global_det_data.
        !!
        !!  This function heavily relies on correctly initialized global data
        !!  outside the control of this class.
        !!  Carefully make sure, that global_det_data is correctly initialized.
        use global_det_data, only: global_lookup => get_prop_vec_idx
        class(PropertyIndexer_t), intent(in) :: this
        integer, intent(in) :: idet
            !!  The index of nI in the FciMCData::CurrentDets array.
        integer, intent(in) :: nI(:)
            !!  The determinant for which the prop_vec index should be calculated.
        integer :: idx
        debug_function_name('lookup_prop_vec_idx')

        if (this%n_prop_vecs() > 1) then
            idx = global_lookup(idet)
            ! Assert that looked up and computed index agree.
            @:pure_ASSERT(idx == this%idx_nI(nI))
        else
            idx = 1
        end if
    end function

    pure function idx_prop_vec(this, prop_vec) result(idx)
        !! Get the index of a property vector.
        !!
        !! Corresponds to the supergroup index for GAS.
        class(PropertyIndexer_t), intent(in) :: this
        integer, intent(in) :: prop_vec(:)
        integer :: idx
        debug_function_name('idx_prop_vec')

        ! These ASSERTs are crucial for the index lookup via compositions.
        ! Otherwise the calculation of the composition index does not make sense.
        @:pure_ASSERT(all(prop_vec >= 0))
        @:pure_ASSERT(sum(prop_vec) == this%prop_vec_sum())
        @:pure_ASSERT(size(prop_vec) == this%prop_vec_dim())

        if (this%n_prop_vecs() > 1) then
            idx = int(binary_search_first_ge(this%allowed_composition_indices, composition_idx(prop_vec)))
        else
            idx = 1
        end if
        @:pure_ASSERT(idx /= -1)
    end function

    pure function get_idx_prop_vec_det(this, nI) result(idx)
        !!  Calculate the property vector index for a determinant nI
        class(PropertyIndexer_t), intent(in) :: this
        integer, intent(in) :: nI(:)
            !! The determinant for which the prop_vec index should be calculated.
        integer :: idx
        debug_function_name('get_idx_prop_vec_det')

        if (this%n_prop_vecs() > 1) then
            idx = this%idx_prop_vec(this%to_prop_vec(nI))
        else
            idx = 1
        end if
        @:pure_ASSERT(idx /= -1)
    end function

    pure function get_n_prop_vecs(this) result(res)
        !! Get the number of possible prop_vecs.
        !!
        !! GAS allowed compositions are called prop_vecs.
        class(PropertyIndexer_t), intent(in) :: this
        integer :: res
        res = size(this%allowed_composition_indices)
    end function

    pure function get_prop_vecs(this) result(res)
        !! Get the ordered compositions of n into k summands
        !!  constrained by cumulative minima and maxima.
        !!
        !! GAS allowed compositions are called prop_vecs.
        class(PropertyIndexer_t), intent(in) :: this
        integer, allocatable :: res(:, :)
        integer(int64) :: i

        allocate(res(this%prop_vec_dim(), this%n_prop_vecs()))
        do i = 1_int64, this%n_prop_vecs()
            res(:, i) = composition_from_idx(&
                this%prop_vec_dim(), this%prop_vec_sum(), this%allowed_composition_indices(i))
        end do
    end function

    elemental integer function prop_vec_dim(this)
        class(PropertyIndexer_t), intent(in) :: this
        prop_vec_dim = this%prop_vec_dim_
    end function

    elemental integer function n_spin_orbs(this)
        class(PropertyIndexer_t), intent(in) :: this
        n_spin_orbs = this%n_spin_orbs_
    end function

    elemental integer function prop_vec_sum(this)
        class(PropertyIndexer_t), intent(in) :: this
        prop_vec_sum = this%prop_vec_sum_
    end function

    pure function get_prop_vec_idx_csf(this, csf_i) result(idx)
        !!  Calculate the prop_vec index for a CSF `csf_i`
        class(AlsoGUGA_PropertyIndexer_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
            !! The CSF for which the prop_vec index should be calculated.
        integer :: idx
        character(*), parameter :: this_routine = 'get_prop_vec_idx_csf'

        if (this%n_prop_vecs() > 1) then
            idx = this%idx_prop_vec(this%to_prop_vec(csf_i))
        else
            idx = 1
        end if
        @:pure_ASSERT(idx /= -1)
    end function

    pure function lookup_prop_vec_idx_CSF(this, csf_i) result(idx)
        !!  Use a precomputed prop_vec index from global_det_data.
        !!
        !!  This function heavily relies on correctly initialized global data
        !!  outside the control of this class.
        !!  Carefully make sure, that global_det_data is correctly initialized.
        use global_det_data, only: global_lookup => get_prop_vec_idx
        class(AlsoGUGA_PropertyIndexer_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
            !! The CSF
        integer :: idx
        debug_function_name('lookup_prop_vec_idx')

        if (this%n_prop_vecs() > 1) then
            idx = global_lookup(csf_i%idx_curr_dets)
            ! Assert that looked up and computed index agree.
            @:pure_ASSERT(idx == this%idx_csf(csf_i))
        else
            idx = 1
        end if
    end function

    logical pure function get_contained_in_PQ(this) result(contained_in_P)
        class(PQ_indexer_t), intent(in) :: this
        logical, allocatable :: contained_in_P(:)
        integer :: i
        
        allocate(contained_in_P(size(this%total_space%allowed_composition_indices)))

        !TODO Loop through entries of total space and check if they are contained in P
        !do i=1,size(contained_in_P)
        !    contained_in_P(i) = .true.
        !end do
    end function
    
    pure function PQ_excite_single(this, exc, prop_vec)
        class(PQ_indexer_t), intent(in) :: this
        type(Excite_1_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)
        integer :: res(size(prop_vec))

        res = this%total_space%excite_single(exc, prop_vec)
    end function

    pure function PQ_excite_double(this, exc, prop_vec)
        class(PQ_indexer_t), intent(in) :: this
        type(Excite_2_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)
        integer :: res(prop_vec())

        res = this%total_space%excite_double(exc, prop_vec)
    end function

    logical pure function PQ_is_allowed_single(this, exc, prop_vec)
        class(PQ_indexer_t), intent(in) :: this
        type(Excite_1_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)

        !Query wether we have a P->P, P->Q, Q->P or Q->Q transition
        (this%contained_in_P(this%prop_vec_idx(prop_vec)) .or. &
         this%contained_in_P(this%prop_vec_idx(this%excite(exc, prop_vec)))) &
        .and. this%total_space%is_allowed(exc, prop_vec)
    end function
    
    logical pure function PQ_is_allowed_double(this, exc, prop_vec)
        class(PQ_indexer_t), intent(in) :: this
        type(Excite_2_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)

        !Query wether we have a P->P, P->Q, Q->P or Q->Q transition
        (this%contained_in_P(this%prop_vec_idx(prop_vec)) .or. &
         this%contained_in_P(this%prop_vec_idx(this%excite(exc, prop_vec)))) &
        .and. this%total_space%is_allowed(exc, prop_vec)
    end function
    
    logical pure function PQ_is_allowed_distinct_double(this, exc, prop_vec)
        class(PQ_indexer_t), intent(in) :: this
        type(DistinctDouble_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)

        integer :: i, j, a, b
        call exc%idx%extract_ijkl(a, i, b, j)

        !Query wether we have a P->P, P->Q, Q->P or Q->Q transition
        (this%contained_in_P(this%prop_vec_idx(prop_vec)) .or. &
         this%contained_in_P(this%prop_vec_idx(this%excite(Excite_2_t(2 * i, 2 * a, 2 * j, 2 * b), prop_vec)))) &
        .and. this%total_space%is_allowed_distinct_double(Excite_2_t(2 * i, 2 * a, 2 * j, 2 * b)c prop_vec)
    end function

    subroutine PQ_gen_all_excits(this, nI, n_excits, det_list, ic)
        class(PQ_indexer_t), intent(in) :: this
        integer, intent(in) :: nI(:)
        integer, intent(out) :: n_excits
        integer(n_int), intent(out), allocatable :: det_list(:,:)
        integer, intent(in), optional :: ic

        call this%total_space%gen_all_excits(nI, n_excits, det_list, ic)
    end subroutine
    
    pure function PQ_to_prop_vec_csf(this, csf_i) result(supergroup)
        class(PQ_indexer_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
        !integer :: supergroup(total_space%prop_vec_dim())
        integer, allocatable(:) :: supergroup

        supergroup = this%total_space%to_prop_vec_csf(csf_i)
        
    end function

    pure function PQ_to_prop_vec_ni(this, nI) result(supergroup)
        class(PQ_indexer_t), intent(in) :: this
        integer, intent(in) :: nI(:)
        !integer :: supergroup(total_space%prop_vec_dim())
        integer, allocatable(:) :: supergroup

        supergroup = this%total_space%to_prop_vec_ni(nI)
        
    end function
!    pure function construct_PQ_indexer_t(P_indexer, total_indexer) result(idxer)
!        class(GAS_SuperGroupIndexer_t) :: P_indexer
!        class(GAS_SuperGroupIndexer_t) :: total_indexer
!        type(PQ_indexer_t) :: idxer
!        integer(int64), allocatable :: allowed_indices(:)
!        character(*), parameter :. this_routine = 'construct_PQ_indexer_t'
!        
!        if (.not. total_indexer%GAS_spec%is_connected()) then
!            stop_all(t_r, "GAS Space cannot be disconnected")
!        end if
!        
!        idxer%contained_in_P
!    end function
        
end module
