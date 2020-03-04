module test_gasci_mod
    use fruit
    use orb_idx_mod, only: SpinOrbIdx_t, SpatOrbIdx_t
    use util_mod, only: cumsum
    use gasci, only: GASSpec_t, get_iGAS, &
        contains_det, get_nGAS, particles_per_GAS, operator(.contains.), &
        is_valid, is_connected, get_possible_spaces
    implicit none
    private
    public :: test_igas_from_spatorb, test_igas_from_spinorb, &
        test_contains_det_spinorb, test_contains_det_spatorb, &
        test_particles_per_GAS_spatorb, test_particles_per_GAS_spinorb, &
        test_is_valid, test_is_connected, &
        test_get_possible_spaces_spinorb, test_get_possible_spaces_spatorb



contains

    subroutine test_igas_from_spatorb()
        type(GASSpec_t) :: GAS_spec
        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[2, 4],  n_max=[2, 4])
        call assert_equals([1, 1, 2, 2], &
                           get_iGAS(GAS_spec, SpatOrbIdx_t([1, 2, 3, 4])), 4)

    end subroutine

    subroutine test_igas_from_spinorb()
        type(GASSpec_t) :: GAS_spec
        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])

        call assert_equals([1, 1, 1, 1, 2, 2, 2, 2], &
                          get_iGAS(GAS_spec, SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 7, 8])), 8)

        call assert_equals([1, 2, 2, 2], &
                          get_iGAS(GAS_spec, SpinOrbIdx_t([1, 5, 6, 7])), 4)

    end subroutine

    subroutine test_particles_per_GAS_spatorb()
        type(GASSpec_t) :: GAS_spec
        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])
        call assert_equals([2, 2], &
                           particles_per_GAS(GAS_spec, SpatOrbIdx_t([1, 2, 3, 4])), 2)
        call assert_equals([1, 3], &
                           particles_per_GAS(GAS_spec, SpatOrbIdx_t([1, 3, 3, 4])), 2)

    end subroutine


    subroutine test_particles_per_GAS_spinorb()
        type(GASSpec_t) :: GAS_spec
        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])
        call assert_equals([2, 2], &
                           particles_per_GAS(GAS_spec, SpinOrbIdx_t([1, 2, 5, 6])), 2)
        call assert_equals([1, 3], &
                           particles_per_GAS(GAS_spec, SpinOrbIdx_t([1, 5, 6, 7])), 2)

    end subroutine

    subroutine test_contains_det_spinorb()
        type(GASSpec_t) :: GAS_spec
        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])

        call assert_true(GAS_spec .contains. SpinOrbIdx_t([1, 2, 5, 6]))
        call assert_true(GAS_spec .contains. SpinOrbIdx_t([1, 3, 5, 6]))
        call assert_false(GAS_spec .contains. SpinOrbIdx_t([1, 2, 3, 4]))
        call assert_false(GAS_spec .contains. SpinOrbIdx_t([1, 2, 3, 5]))
        call assert_false(GAS_spec .contains. SpinOrbIdx_t([5, 6, 7, 8]))
        call assert_false(GAS_spec .contains. SpinOrbIdx_t([1, 6, 7, 8]))
    end subroutine

    subroutine test_contains_det_spatorb()
        type(GASSpec_t) :: GAS_spec
        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])

        call assert_true(GAS_spec .contains. SpatOrbIdx_t([1, 1, 3, 3]))
        call assert_true(GAS_spec .contains. SpatOrbIdx_t([1, 2, 3, 4]))
        call assert_false(GAS_spec .contains. SpatOrbIdx_t([1, 1, 2, 3]))
        call assert_false(GAS_spec .contains. SpatOrbIdx_t([1, 1, 2, 2]))
        call assert_false(GAS_spec .contains. SpatOrbIdx_t([1, 3, 3, 4]))
        call assert_false(GAS_spec .contains. SpatOrbIdx_t([3, 3, 4, 4]))
    end subroutine

    subroutine test_is_valid()
        call assert_true(is_valid(GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])))
        call assert_true(is_valid(GASSpec_t(n_orbs=[1, 4], n_min=[2, 4], n_max=[2, 4])))
        call assert_true(is_valid(GASSpec_t(n_orbs=[2, 4], n_min=[1, 4], n_max=[3, 4])))
        call assert_true(is_valid(GASSpec_t(n_orbs=[2, 4], n_min=[1, 4], n_max=[3, 4]), &
                                  n_particles=4))
        call assert_true(is_valid(GASSpec_t(n_orbs=[2, 4], n_min=[1, 4], n_max=[3, 4]), &
                                  n_particles=4, n_basis=8))

        call assert_false(is_valid(GASSpec_t(n_orbs=[1, 4], n_min=[3, 4], n_max=[3, 4])))
        call assert_false(is_valid(GASSpec_t(n_orbs=[2, 4], n_min=[3, 3], n_max=[3, 4])))
        call assert_false(is_valid(GASSpec_t(n_orbs=[1, 2], n_min=[3, 5], n_max=[3, 5])))
        call assert_false(is_valid(GASSpec_t(n_orbs=[1, 2], n_min=[3, 5], n_max=[3, 5])))
        call assert_false(is_valid(GASSpec_t(n_orbs=[2, 4], n_min=[1, 4], n_max=[3, 4]), &
                                   n_particles=5))
        call assert_false(is_valid(GASSpec_t(n_orbs=[2, 4], n_min=[1, 4], n_max=[3, 4]), &
                                   n_particles=4, n_basis=5))
        call assert_false(is_valid(GASSpec_t(n_orbs=[2, 4], n_min=[1, 4], n_max=[3, 4]), &
                                   n_particles=5, n_basis=5))
    end subroutine


    subroutine test_is_connected()
        call assert_true(is_connected(GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])))
        call assert_true(is_connected(GASSpec_t(n_orbs=[2, 4], n_min=[2, 5], n_max=[2, 5])))

        call assert_false(is_connected(GASSpec_t(n_orbs=[2, 4], n_min=[1, 4], n_max=[3, 4])))
        call assert_false(is_connected(GASSpec_t(n_orbs=[2, 4], n_min=[1, 5], n_max=[3, 5])))
    end subroutine


    subroutine test_get_possible_spaces_spinorb()
        type(GASSpec_t) :: GAS_spec
        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])

        call assert_equals( &
            [0, 0], &
            get_possible_spaces(GAS_spec, SpinOrbIdx_t([1, 2, 5, 6])), &
            2)

        call assert_equals( &
            [1, 1], &
            get_possible_spaces(GAS_spec, SpinOrbIdx_t([1, 2, 5, 6]), &
                                additional_holes=SpinOrbIdx_t([1])), &
            2)

        call assert_equals( &
            [2, 2], &
            get_possible_spaces(GAS_spec, SpinOrbIdx_t([1, 2, 5, 6]), &
                                additional_holes=SpinOrbIdx_t([5])), &
            2)

        call assert_equals( &
            [0, 0], &
            get_possible_spaces(GAS_spec, SpinOrbIdx_t([1, 2, 5, 6]), &
                                additional_holes=SpinOrbIdx_t([1, 5]), &
                                n_particles=1), &
            2)

        call assert_equals( &
            [1, 2], &
            get_possible_spaces(GAS_spec, SpinOrbIdx_t([1, 2, 5, 6]), &
                                additional_holes=SpinOrbIdx_t([1, 5]), &
                                n_particles=2), &
            2)


        call assert_equals( &
            [1, 1], &
            get_possible_spaces(GAS_spec, SpinOrbIdx_t([1, 2, 5, 6]), &
                                additional_holes=SpinOrbIdx_t([1, 2]), &
                                n_particles=2), &
            2)

        call assert_equals( &
            [2, 2], &
            get_possible_spaces(GAS_spec, SpinOrbIdx_t([1, 2, 5, 6]), &
                                additional_holes=SpinOrbIdx_t([5, 6]), &
                                n_particles=2), &
            2)
    end subroutine


    subroutine test_get_possible_spaces_spatorb()
        type(GASSpec_t) :: GAS_spec
        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])

        call assert_equals( &
            [0, 0], &
            get_possible_spaces(GAS_spec, SpatOrbIdx_t([1, 1, 3, 3])), &
            2)

        call assert_equals( &
            [1, 1], &
            get_possible_spaces(GAS_spec, SpatOrbIdx_t([1, 1, 3, 3]), &
                                additional_holes=SpatOrbIdx_t([1])), &
            2)

        call assert_equals( &
            [2, 2], &
            get_possible_spaces(GAS_spec, SpatOrbIdx_t([1, 1, 3, 3]), &
                                additional_holes=SpatOrbIdx_t([3])), &
            2)

        call assert_equals( &
            [0, 0], &
            get_possible_spaces(GAS_spec, SpatOrbIdx_t([1, 1, 3, 3]), &
                                additional_holes=SpatOrbIdx_t([1, 3]), &
                                n_particles=1), &
            2)

        call assert_equals( &
            [1, 2], &
            get_possible_spaces(GAS_spec, SpatOrbIdx_t([1, 1, 3, 3]), &
                                additional_holes=SpatOrbIdx_t([1, 3]), &
                                n_particles=2), &
            2)


        call assert_equals( &
            [1, 1], &
            get_possible_spaces(GAS_spec, SpatOrbIdx_t([1, 1, 3, 3]), &
                                additional_holes=SpatOrbIdx_t([1, 1]), &
                                n_particles=2), &
            2)

        call assert_equals( &
            [2, 2], &
            get_possible_spaces(GAS_spec, SpatOrbIdx_t([1, 1, 3, 3]), &
                                additional_holes=SpatOrbIdx_t([3, 3]), &
                                n_particles=2), &
            2)
    end subroutine

end module test_gasci_mod

program test_gasci_program

    use mpi
    use fruit
    use test_gasci_mod, only: test_igas_from_spatorb, test_igas_from_spinorb, &
        test_contains_det_spinorb, test_contains_det_spatorb, &
        test_particles_per_GAS_spatorb, test_particles_per_GAS_spinorb, &
        test_is_valid, test_is_connected, &
        test_get_possible_spaces_spinorb, test_get_possible_spaces_spatorb


    implicit none
    integer :: failed_count, err

    integer :: n

    call mpi_init(err)

    call init_fruit()

    call test_gasci_driver()

    call fruit_summary()
    call fruit_finalize()
    call get_failed_count(failed_count)

    if (failed_count /= 0) call stop_all('test_gasci_program', 'failed_tests')

    call mpi_finalize(err)

contains

    subroutine test_gasci_driver()
        call run_test_case(test_igas_from_spatorb, "test_igas_from_spatorb")
        call run_test_case(test_igas_from_spinorb, "test_igas_from_spinorb")
        call run_test_case(test_contains_det_spinorb, "test_contains_det_spinorb")
        call run_test_case(test_contains_det_spatorb, "test_contains_det_spatorb")
        call run_test_case(test_particles_per_GAS_spatorb, "test_particles_per_GAS_spatorb")
        call run_test_case(test_particles_per_GAS_spinorb, "test_particles_per_GAS_spinorb")
        call run_test_case(test_is_valid, "test_is_valid")
        call run_test_case(test_is_connected, "test_is_connected")
        call run_test_case(test_get_possible_spaces_spinorb, "test_get_possible_spaces_spinorb")
        call run_test_case(test_get_possible_spaces_spatorb, "test_get_possible_spaces_spatorb")
    end subroutine
end program test_gasci_program
