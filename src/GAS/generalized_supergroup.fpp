#include "macros.h"
#:include "../macros.fpph"
#:include "../algorithms.fpph"

module generalized_supergroup
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

    better_implicit_none

    private
    public :: SuperGroupIndexer_t, AlsoGUGA_SuperGroupIndexer_t, &
         indexer, lookup_supergroup_indexer

    type, abstract :: SuperGroupIndexer_t
        integer(int64), allocatable :: allowed_composition_indices(:)
        integer :: n_El_, n_spin_orbs_, n_sg_dim_
            !! The particle number, the number of spin orbitals,
            !!    and the dimension of the supergroups.
            !!   (corresponds to number of GAS spaces for pure GAS)
    contains
        private
        procedure, public :: idx_supergroup => get_supergroup_idx
        procedure, public :: idx_nI => get_supergroup_idx_det
        generic, public :: lookup_supergroup_idx => lookup_supergroup_idx_nI
        procedure :: lookup_supergroup_idx_nI
        procedure, public :: n_supergroups => get_n_supergroups
        procedure, public :: get_supergroups => indexer_get_supergroups
        procedure(is_allowed_single_t), deferred :: is_allowed_single
        procedure(is_allowed_double_t), deferred :: is_allowed_double
        procedure(is_allowed_distinct_double_t), deferred :: is_allowed_distinct_double
        generic, public :: is_allowed => is_allowed_single, is_allowed_double, is_allowed_distinct_double

        procedure(idxer_gen_all_excits_t), deferred, public :: gen_all_excits

        procedure(to_supergroup_nI_t), deferred :: to_supergroup_nI
        generic, public :: to_supergroup => to_supergroup_nI

        procedure, public :: dim_sg => get_dim_sg
        procedure, public :: n_spin_orbs => get_nOrbs
        procedure, public :: n_El => get_nEl
    end type

    type, abstract, extends(SuperGroupIndexer_t) :: AlsoGUGA_SuperGroupIndexer_t
    contains
        private
        procedure, public :: idx_csf => get_supergroup_idx_csf
        generic, public :: lookup_supergroup_idx => lookup_supergroup_idx_CSF
        procedure :: lookup_supergroup_idx_CSF

        procedure(to_supergroup_csf_t), deferred :: to_supergroup_csf
        generic, public :: to_supergroup => to_supergroup_csf
    end type

    abstract interface

        logical pure function is_allowed_single_t(this, exc, supergroup)
            import :: SuperGroupIndexer_t, Excite_1_t
            !! Check if a single excitation is allowed.
            !!
            !! Is called once at initialization, so it does not have to be super fast.
            class(SuperGroupIndexer_t), intent(in) :: this
            type(Excite_1_t), intent(in) :: exc
            integer, intent(in) :: supergroup(:)
        end function

        logical pure function is_allowed_double_t(this, exc, supergroup)
            import :: SuperGroupIndexer_t, Excite_2_t
            !! Check if a double excitation is allowed.
            !!
            !! Is called once at initialization, so it does not have to be super fast.
            class(SuperGroupIndexer_t), intent(in) :: this
            type(Excite_2_t), intent(in) :: exc
            integer, intent(in) :: supergroup(:)
        end function

        logical pure function is_allowed_distinct_double_t(this, exc, supergroup)
            import :: SuperGroupIndexer_t, DistinctDouble_t
            class(SuperGroupIndexer_t), intent(in) :: this
            type(DistinctDouble_t), intent(in) :: exc
            integer, intent(in) :: supergroup(:)
        end function

        subroutine idxer_gen_all_excits_t(this, nI, n_excits, det_list, ic)
            import :: SuperGroupIndexer_t, n_int
            !!  Get all excitated determinants from
            !!  det_I that are allowed under GAS constraints.
            class(SuperGroupIndexer_t), intent(in) :: this
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

        pure function to_supergroup_nI_t(this, nI) result(supergroup)
            import :: SuperGroupIndexer_t
            class(SuperGroupIndexer_t), intent(in) :: this
            integer, intent(in) :: nI(:)
            integer :: supergroup(this%dim_sg())
        end function

        pure function to_supergroup_csf_t(this, csf_i) result(supergroup)
            import :: AlsoGUGA_SuperGroupIndexer_t, CSF_Info_t
            class(AlsoGUGA_SuperGroupIndexer_t), intent(in) :: this
            type(CSF_Info_t), intent(in) :: csf_i
            integer :: supergroup(this%dim_sg())
        end function

    end interface

    class(SuperGroupIndexer_t), allocatable :: indexer
    class(SuperGroupIndexer_t), pointer :: lookup_supergroup_indexer => null()

contains

    pure function lookup_supergroup_idx_nI(this, idet, nI) result(idx)
        !!  Use a precomputed supergroup index from global_det_data.
        !!
        !!  This function heavily relies on correctly initialized global data
        !!  outside the control of this class.
        !!  Carefully make sure, that global_det_data is correctly initialized.
        use global_det_data, only: global_lookup => get_supergroup_idx
        class(SuperGroupIndexer_t), intent(in) :: this
        integer, intent(in) :: idet
            !!  The index of nI in the FciMCData::CurrentDets array.
        integer, intent(in) :: nI(:)
            !!  The determinant for which the supergroup index should be calculated.
        integer :: idx
        debug_function_name('lookup_supergroup_idx')

        if (this%n_supergroups() > 1) then
            idx = global_lookup(idet)
            ! Assert that looked up and computed index agree.
            @:pure_ASSERT(idx == this%idx_nI(nI))
        else
            idx = 1
        end if
    end function

    pure function get_supergroup_idx(this, supergroup) result(idx)
        class(SuperGroupIndexer_t), intent(in) :: this
        integer, intent(in) :: supergroup(:)
        integer :: idx
        character(*), parameter :: this_routine = 'get_supergroup_idx'

        if (this%n_supergroups() > 1) then
            idx = int(binary_search_first_ge(this%allowed_composition_indices, composition_idx(supergroup)))
        else
            idx = 1
        end if
        @:pure_ASSERT(idx /= -1)
    end function

    pure function get_supergroup_idx_det(this, nI) result(idx)
        !!  Calculate the supergroup index for a determinant nI
        class(SuperGroupIndexer_t), intent(in) :: this
        integer, intent(in) :: nI(:)
            !! The determinant for which the supergroup index should be calculated.
        integer :: idx
        character(*), parameter :: this_routine = 'get_supergroup_idx_det'

        if (this%n_supergroups() > 1) then
            idx = this%idx_supergroup(this%to_supergroup(nI))
        else
            idx = 1
        end if
        @:pure_ASSERT(idx /= -1)
    end function

    pure function get_n_supergroups(this) result(res)
        !! Get the number of possible supergroups.
        !!
        !! GAS allowed compositions are called supergroups.
        class(SuperGroupIndexer_t), intent(in) :: this
        integer :: res
        res = size(this%allowed_composition_indices)
    end function

    pure function indexer_get_supergroups(this) result(res)
        !! Get the ordered compositions of n into k summands
        !!  constrained by cumulative minima and maxima.
        !!
        !! GAS allowed compositions are called supergroups.
        class(SuperGroupIndexer_t), intent(in) :: this
        integer, allocatable :: res(:, :)
        integer(int64) :: i

        allocate(res(this%dim_sg(), this%n_supergroups()))
        do i = 1_int64, this%n_supergroups()
            res(:, i) = composition_from_idx(&
                this%dim_sg(), this%n_El(), this%allowed_composition_indices(i))
        end do
    end function

    elemental integer function get_dim_sg(this)
        class(SuperGroupIndexer_t), intent(in) :: this
        get_dim_sg = this%n_sg_dim_
    end function

    elemental integer function get_nOrbs(this)
        class(SuperGroupIndexer_t), intent(in) :: this
        get_nOrbs = this%n_spin_orbs_
    end function

    elemental integer function get_nEl(this)
        class(SuperGroupIndexer_t), intent(in) :: this
        get_nEl = this%n_El_
    end function

    pure function get_supergroup_idx_csf(this, csf_i) result(idx)
        !!  Calculate the supergroup index for a CSF `csf_i`
        class(AlsoGUGA_SuperGroupIndexer_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
            !! The CSF for which the supergroup index should be calculated.
        integer :: idx
        character(*), parameter :: this_routine = 'get_supergroup_idx_csf'

        if (this%n_supergroups() > 1) then
            idx = this%idx_supergroup(this%to_supergroup(csf_i))
        else
            idx = 1
        end if
        @:pure_ASSERT(idx /= -1)
    end function

    pure function lookup_supergroup_idx_CSF(this, csf_i) result(idx)
        !!  Use a precomputed supergroup index from global_det_data.
        !!
        !!  This function heavily relies on correctly initialized global data
        !!  outside the control of this class.
        !!  Carefully make sure, that global_det_data is correctly initialized.
        use global_det_data, only: global_lookup => get_supergroup_idx
        class(AlsoGUGA_SuperGroupIndexer_t), intent(in) :: this
        type(CSF_Info_t), intent(in) :: csf_i
            !! The CSF
        integer :: idx
        debug_function_name('lookup_supergroup_idx')

        if (this%n_supergroups() > 1) then
            idx = global_lookup(csf_i%idx_curr_dets)
            ! Assert that looked up and computed index agree.
            @:pure_ASSERT(idx == this%idx_csf(csf_i))
        else
            idx = 1
        end if
    end function

end module
