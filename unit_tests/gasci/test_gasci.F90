module test_gasci_mod
    use fruit
    use constants, only: dp, n_int
    use bit_rep_data, only: NIfTot
    use SystemData, only: nEl
    use sets_mod, only: disjoint
    use util_mod, only: operator(.div.), operator(.isclose.)
    use sort_mod, only: sort
    use procedure_pointers, only: generate_excitation
    use orb_idx_mod, only: SpinOrbIdx_t, SpatOrbIdx_t, SpinProj_t, &
        size, operator(==), alpha, beta, sum, calc_spin, calc_spin_raw, &
        operator(-), to_ilut, write_det, operator(/=)
    use excitation_types, only: Excitation_t, SingleExc_t, DoubleExc_t, excite, dyn_excite
    use util_mod, only: cumsum

    use gasci, only: LocalGASSpec_t

    use unit_test_helper_excitgen, only: test_excitation_generator, &
        init_excitgen_test, finalize_excitgen_test, generate_random_integrals, &
        FciDumpWriter_t
    use DetBitOps, only: ilut_lt, ilut_gt
    implicit none
    private
    public :: test_igas, test_contains_det, test_particles_per_GAS, &
        test_is_valid, test_is_connected, test_get_orb_idx



contains

    subroutine test_igas()
        type(LocalGASSpec_t) :: GAS_spec
        GAS_spec = LocalGASSpec_t(n_min=[2, 2],  n_max=[2, 2], spat_GAS_orbs=[1, 1, 2, 2])

        call assert_equals([1, 1, 1, 1, 2, 2, 2, 2], &
                           GAS_spec%get_iGAS([1, 2, 3, 4, 5, 6, 7, 8]), 8)

        call assert_equals([1, 1, 2, 2], &
                           GAS_spec%get_iGAS([1, 4, 6, 7]), 4)
    end subroutine


    subroutine test_particles_per_GAS()
        type(LocalGASSpec_t) :: GAS_spec
        GAS_spec = LocalGASSpec_t(n_min=[2, 2],  n_max=[2, 2], spat_GAS_orbs=[1, 1, 2, 2])

        associate(expected => [2, 2], &
                  calculated => GAS_spec%count_per_GAS([1, 2, 5, 6]))
            call assert_equals(expected, calculated, size(expected))
        end associate

        associate(expected => [1, 3], &
                  calculated => GAS_spec%count_per_GAS([1, 5, 6, 7]))
            call assert_equals(expected, calculated, size(expected))
        end associate

    end subroutine

    subroutine test_contains_det()
        type(LocalGASSpec_t) :: GAS_spec
        GAS_spec = LocalGASSpec_t(n_min=[2, 2],  n_max=[2, 2], spat_GAS_orbs=[1, 1, 2, 2])

        call assert_true( GAS_spec%contains_conf([1, 2, 5, 6]))
        call assert_true( GAS_spec%contains_conf([1, 3, 5, 6]))
        call assert_false(GAS_spec%contains_conf([1, 2, 3, 4]))
        call assert_false(GAS_spec%contains_conf([1, 2, 3, 5]))
        call assert_false(GAS_spec%contains_conf([5, 6, 7, 8]))
        call assert_false(GAS_spec%contains_conf([1, 6, 7, 8]))
    end subroutine

    subroutine test_is_valid()
        type(LocalGASSpec_t) :: GAS_spec
        GAS_spec = LocalGASSpec_t(n_min=[2, 2], n_max=[2, 2], spat_GAS_orbs=[1, 1, 2, 2])
        call assert_true(GAS_spec%is_valid())
        GAS_spec = LocalGASSpec_t(n_min=[2, 2], n_max=[2, 2], spat_GAS_orbs=[1, 2, 2, 2])
        call assert_true(GAS_spec%is_valid())
        GAS_spec = LocalGASSpec_t(n_min=[1, 1], n_max=[3, 3], spat_GAS_orbs=[1, 1, 2, 2])
        call assert_true(GAS_spec%is_valid())
        GAS_spec = LocalGASSpec_t(n_min=[1, 1], n_max=[3, 3], spat_GAS_orbs=[1, 1, 2, 2])
        call assert_true(GAS_spec%is_valid(n_basis=8))
    end subroutine


    subroutine test_is_connected()
        integer, parameter :: GAS_table(4) = [1, 1, 2, 2]
        type(LocalGASSpec_t) :: GAS_spec
        GAS_spec = LocalGASSpec_t(n_min=[1, 1], n_max=[3, 3], spat_GAS_orbs=GAS_table)
        call assert_true(GAS_spec%is_connected())
        GAS_spec = LocalGASSpec_t(n_min=[1, 2], n_max=[3, 4], spat_GAS_orbs=GAS_table)
        call assert_true(GAS_spec%is_connected())

        GAS_spec = LocalGASSpec_t(n_min=[2, 2], n_max=[2, 2], spat_GAS_orbs=GAS_table)
        call assert_false(GAS_spec%is_connected())
        GAS_spec = LocalGASSpec_t(n_min=[2, 3], n_max=[2, 3], spat_GAS_orbs=GAS_table)
        call assert_false(GAS_spec%is_connected())
    end subroutine

    subroutine test_get_orb_idx
        type(LocalGASSpec_t) :: GAS_spec
        GAS_spec = LocalGASSpec_t(n_min=[2, 2], n_max=[2, 2], spat_GAS_orbs=[1, 1, 2, 2])

        call assert_equals(GAS_spec%get_orb_idx([1, 2], 1), [1, 2], 2)
        call assert_equals(GAS_spec%get_orb_idx([1, 2], 2), [5, 6], 2)
    end subroutine

end module test_gasci_mod

program test_gasci_program

    use mpi
    use fruit
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_gasci_mod, only: test_igas, test_is_connected, test_is_valid, &
        test_contains_det, test_particles_per_GAS, test_get_orb_idx

    implicit none
    integer :: failed_count
    block

        call MPIInit(.false.)

        call init_fruit()

        call test_gasci_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_gasci_program', 'failed_tests')

        call MPIEnd(.false.)
    end block

contains

    subroutine test_gasci_driver()
        call run_test_case(test_is_valid, "test_is_valid")
        call run_test_case(test_is_connected, "test_is_connected")
        call run_test_case(test_igas, "test_igas")
        call run_test_case(test_get_orb_idx, "test_get_orb_idx")
        call run_test_case(test_contains_det, "test_contains_det")
        call run_test_case(test_particles_per_GAS, "test_particles_per_GAS")
    end subroutine
end program test_gasci_program
