#include "macros.h"
#:include "../macros.fpph"
#:include "../algorithms.fpph"


module property_vector_index
    !! This module generalizes the "supergroup" concept of PropVec to arbitrary property vectors.
    !! This allows more general constraints such as Löwdin's perturbation theory via PQ space folding
    !!  or constraints that use e.g. the spin projection.
    !!
    !! A property vector is a vector of integers where each element encodes a property of a sub-space.
    !! In traditional PropVec the property is just the particle number, while in general these can be e.g.
    !! `[N_1, Sz_1, N_2, Sz_2]`
    !! where `N_i` is the particle number in the i-th sub-space and `Sz_i` is the spin projection
    !! in the i-th sub-space (encoded as integer).
    !! The only assumed constraint of the property vectors is that the sum over them is preserved.
    !!
    !! The algorithm follows the ideas of stochastic PropVec (https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00936),
    !! i.e. the property vectors are indexed by calculating the corresponding composition index.
    use constants, only: int64, n_int
    use util_mod, only: choose_i64, cumsum, binary_search_first_ge, custom_findloc, stop_all
    use bit_rep_data, only: nIfD
    use hash, only: hash_table_lookup
    use guga_bitRepOps, only: CSF_Info_t
    use growing_buffers, only: buffer_int_2D_t, buffer_int64_1D_t
    use composition_utils, only: n_compositions, get_compositions, composition_idx, &
        composition_from_idx, next_composition
    use excitation_types, only: Excite_1_t, Excite_2_t, Excitation_t
    use guga_bitRepOps, only: CSF_Info_t
    use guga_data, only: DistinctDouble_t
!    use gasci_supergroup_index, only: GAS_SuperGroupIndexer_t
    better_implicit_none

    private
    public :: PropertyIndexer_t, AlsoGUGA_PropertyIndexer_t,&
         indexer, lookup_property_indexer

    type, abstract :: PropertyIndexer_t
        integer(int64), allocatable :: allowed_composition_indices(:)
        integer :: n_spin_orbs_
            !! The number of spin orbitals of underlying determinants.
        integer :: prop_vec_sum_
            !! The sum of the property vectors.
            !! In the case of PropVec this is the particle number.
        integer :: prop_vec_dim_
            !! The dimension of the property vector
            !!   (corresponds to number of PropVec spaces for pure PropVec)
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

        procedure, public :: dyn_is_allowed

        procedure(idxer_gen_all_excits_t), deferred, public :: gen_all_excits

        procedure(to_prop_vec_nI_t), deferred :: to_prop_vec_nI
        generic, public :: to_prop_vec => to_prop_vec_nI

        procedure, public :: prop_vec_dim
        procedure, public :: n_spin_orbs
        procedure, public :: prop_vec_sum

        procedure(write_to_t), public, deferred :: write_to
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
            !!  det_I that are allowed under PropVec constraints.
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

        subroutine write_to_t(this, iunit)
            !! Write a string representation of this PropVec specification to iunit
            import :: PropertyIndexer_t
            implicit none
            class(PropertyIndexer_t), intent(in) :: this
            integer, intent(in) :: iunit
        end subroutine

    end interface

    class(PropertyIndexer_t), allocatable :: indexer, lookup_property_indexer

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
        !! Corresponds to the supergroup index for PropVec.
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
        !! PropVec allowed compositions are called prop_vecs.
        class(PropertyIndexer_t), intent(in) :: this
        integer :: res
        res = size(this%allowed_composition_indices)
    end function

    pure function get_prop_vecs(this) result(res)
        !! Get the ordered compositions of n into k summands
        !!  constrained by cumulative minima and maxima.
        !!
        !! PropVec allowed compositions are called prop_vecs.
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


    logical pure function dyn_is_allowed(this, exc, prop_vec)
        !! Check if an excitation is allowed.
        !!
        !! Is called once at initialization, so it does not have to be super fast.
        class(PropertyIndexer_t), intent(in) :: this
        class(Excitation_t), intent(in) :: exc
        integer, intent(in) :: prop_vec(:)
        routine_name("dyn_is_allowed")

        select type(exc)
        type is(Excite_1_t)
            dyn_is_allowed = this%is_allowed(exc, prop_vec)
        type is(Excite_2_t)
            dyn_is_allowed = this%is_allowed(exc, prop_vec)
        class default
            call stop_all(this_routine, "invalid excitation passed in.")
        end select
    end function

end module
