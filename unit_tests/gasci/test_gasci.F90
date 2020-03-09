module test_gasci_mod
    use fruit
    use orb_idx_mod, only: SpinOrbIdx_t, SpatOrbIdx_t, Spin_t, &
        size, operator(==), spin => spin_values
    use excitation_types, only: SingleExc_t, excite
    use util_mod, only: cumsum
    use gasci, only: GASSpec_t, get_iGAS, &
        contains_det, get_nGAS, particles_per_GAS, operator(.contains.), &
        is_valid, is_connected, get_possible_spaces, get_possible_holes, &
        split_per_GAS
    implicit none(type, external)
    private
    public :: test_igas_from_spatorb, test_igas_from_spinorb, &
        test_contains_det_spinorb, test_contains_det_spatorb, &
        test_particles_per_GAS_spatorb, test_particles_per_GAS_spinorb, &
        test_is_valid, test_is_connected, &
        test_get_possible_spaces_spinorb, test_get_possible_spaces_spatorb, &
        test_possible_holes, test_split_per_GAS



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

        associate(expected => [2, 2], &
                  calculated => &
                    particles_per_GAS(split_per_GAS(GAS_spec, SpatOrbIdx_t([1, 2, 3, 4]))))
            call assert_equals(expected, calculated, size(expected))
        end associate

        associate(expected => [1, 3], &
                  calculated => &
                    particles_per_GAS(split_per_GAS(GAS_spec, SpatOrbIdx_t([1, 3, 3, 4]))))
            call assert_equals(expected, calculated, size(expected))
        end associate
    end subroutine


    subroutine test_particles_per_GAS_spinorb()
        type(GASSpec_t) :: GAS_spec
        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])

        associate(expected => [2, 2], &
                  calculated => &
                    particles_per_GAS(split_per_GAS(GAS_spec, SpinOrbIdx_t([1, 2, 5, 6]))))
            call assert_equals(expected, calculated, size(expected))
        end associate

        associate(expected => [1, 3], &
                  calculated => &
                    particles_per_GAS(split_per_GAS(GAS_spec, SpinOrbIdx_t([1, 5, 6, 7]))))
            call assert_equals(expected, calculated, size(expected))
        end associate

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
        type(SpinOrbIdx_t), allocatable :: splitted_det_I(:)
        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])
        splitted_det_I = split_per_GAS(GAS_spec, SpinOrbIdx_t([1, 2, 5, 6]))


        call assert_equals(0, size(get_possible_spaces(GAS_spec, splitted_det_I)))

        call assert_equals( &
            [1, 1], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpinOrbIdx_t([1])), &
            2)

        call assert_equals( &
            [2, 2], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpinOrbIdx_t([5])), &
            2)

        call assert_equals( &
            0, &
            size(get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpinOrbIdx_t([1, 5]), &
                                n_particles=1)))

        call assert_equals( &
            [1, 2], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpinOrbIdx_t([1, 5]), &
                                n_particles=2), &
            2)


        call assert_equals( &
            [1, 1], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpinOrbIdx_t([1, 2]), &
                                n_particles=2), &
            2)

        call assert_equals( &
            [2, 2], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpinOrbIdx_t([5, 6]), &
                                n_particles=2), &
            2)
    end subroutine


    subroutine test_get_possible_spaces_spatorb()
        type(GASSpec_t) :: GAS_spec
        type(SpatOrbIdx_t), allocatable :: splitted_det_I(:)
        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])
        splitted_det_I = split_per_GAS(GAS_spec, SpatOrbIdx_t([1, 1, 3, 3]))

        call assert_equals(0, size(get_possible_spaces(GAS_spec, splitted_det_I)))

        call assert_equals( &
            [1, 1], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpatOrbIdx_t([1])), &
            2)

        call assert_equals( &
            [2, 2], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpatOrbIdx_t([3])), &
            2)

        call assert_equals(0, &
            size(get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpatOrbIdx_t([1, 3]), &
                                n_particles=1)))

        call assert_equals( &
            [1, 2], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpatOrbIdx_t([1, 3]), &
                                n_particles=2), &
            2)


        call assert_equals( &
            [1, 1], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpatOrbIdx_t([1, 1]), &
                                n_particles=2), &
            2)

        call assert_equals( &
            [2, 2], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpatOrbIdx_t([3, 3]), &
                                n_particles=2), &
            2)
    end subroutine

    subroutine test_possible_holes
        type(GASSpec_t) :: GAS_spec
        type(SpinOrbIdx_t) :: reference
        type(SpinOrbIdx_t) :: expected, calculated
        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])
        reference = SpinOrbIdx_t([1, 2, 5, 6])

        expected = SpinOrbIdx_t([integer::])
        calculated = get_possible_holes(GAS_spec, reference)
        call assert_true(size(expected) == size(calculated))
        call assert_true(all(expected == calculated))

        expected = SpinOrbIdx_t([3, 4])
        calculated = get_possible_holes( &
              GAS_spec, reference, add_holes=SpinOrbIdx_t([1]))
        call assert_true(all(expected == calculated))

        expected = SpinOrbIdx_t([integer::])
        calculated = get_possible_holes( &
              GAS_spec, reference, add_holes=SpinOrbIdx_t([1, 2]))
        call assert_true(all(expected == calculated))

        expected = SpinOrbIdx_t([7, 8])
        calculated = get_possible_holes( &
                        GAS_spec, reference, add_holes=SpinOrbIdx_t([5]))
        call assert_true(size(expected) == size(calculated))
        call assert_true(all(expected == calculated))

        expected = SpinOrbIdx_t([7])
        calculated = get_possible_holes( &
              GAS_spec, reference, add_holes=SpinOrbIdx_t([5]), &
              m_s = Spin_t(spin%beta))
        call assert_true(size(expected) == size(calculated))
        call assert_true(all(expected == calculated))

        expected = SpinOrbIdx_t([8])
        calculated = get_possible_holes( &
              GAS_spec, reference, add_holes=SpinOrbIdx_t([5]), &
              m_s = Spin_t(spin%alpha))
        call assert_true(size(expected) == size(calculated))
        call assert_true(all(expected == calculated))

        GAS_spec = GASSpec_t(n_orbs=[2, 4, 6], n_min=[1, 3, 6], n_max=[3, 5, 6])

        reference = SpinOrbIdx_t([1, 2, 5, 6, 9, 10])
        expected = SpinOrbIdx_t([3, 4, 7, 8, 11, 12])
        calculated = get_possible_holes( &
              GAS_spec, reference, add_holes=SpinOrbIdx_t([5]))
        call assert_true(size(expected) == size(calculated))
        call assert_true(all(expected == calculated))

        reference = SpinOrbIdx_t([1, 5, 6, 7, 9, 10])
        expected = SpinOrbIdx_t([2, 3, 4])
        calculated = get_possible_holes( &
              GAS_spec, reference, add_holes=SpinOrbIdx_t([1]))
        call assert_true(size(expected) == size(calculated))
        call assert_true(all(expected == calculated))

        reference = SpinOrbIdx_t([1, 5, 6, 7, 9, 10])
        expected = SpinOrbIdx_t([2, 3, 4, 8, 11, 12])
        calculated = get_possible_holes( &
              GAS_spec, reference, add_holes=SpinOrbIdx_t([1, 5]), &
              n_particles=2)
        call assert_true(size(expected) == size(calculated))
        call assert_true(all(expected == calculated))

        expected = SpinOrbIdx_t([2, 3, 4])
        calculated = get_possible_holes( &
              GAS_spec, excite(reference, SingleExc_t(5, 11)), &
              add_holes=SpinOrbIdx_t([1]))
        call assert_true(size(expected) == size(calculated))
        call assert_true(all(expected == calculated))

        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[0, 4], n_max=[0, 4])
        reference = SpinOrbIdx_t([5, 6, 7, 8])
        expected = SpinOrbIdx_t([integer::])
        calculated = get_possible_holes( &
              GAS_spec, reference, add_holes=SpinOrbIdx_t([5]))
        call assert_true(size(expected) == size(calculated))
        call assert_true(all(expected == calculated))
    end subroutine

    subroutine test_split_per_GAS
        type(GASSpec_t) :: GAS_spec
        integer :: i
        GAS_spec = GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])

        associate(expected => [SpinOrbIdx_t([1, 2]), SpinOrbIdx_t([5, 6])], &
                  calculated => split_per_GAS(GAS_spec, SpinOrbIdx_t([1, 2, 5, 6])))
            call assert_true(all(expected(1) == calculated(1)))
            call assert_true(all(expected(2) == calculated(2)))
        end associate

        associate(expected => [SpinOrbIdx_t([integer::]), SpinOrbIdx_t([5, 6])], &
                  calculated => split_per_GAS(GAS_spec, SpinOrbIdx_t([5, 6])))
            call assert_true(all(expected(1) == calculated(1)))
            call assert_true(all(expected(2) == calculated(2)))
        end associate

        associate(expected => [SpinOrbIdx_t([1, 2, 3]), SpinOrbIdx_t([integer::])], &
                  calculated => split_per_GAS(GAS_spec, SpinOrbIdx_t([1, 2, 3])))
            call assert_true(all(expected(1) == calculated(1)))
            call assert_true(all(expected(2) == calculated(2)))
        end associate

        associate(expected => [SpatOrbIdx_t([1, 1]), SpatOrbIdx_t([3, 3])], &
                  calculated => split_per_GAS(GAS_spec, SpatOrbIdx_t([1, 1, 3, 3])))
            call assert_true(all(expected(1) == calculated(1)))
            call assert_true(all(expected(2) == calculated(2)))
        end associate

        associate(expected => [SpatOrbIdx_t([integer::]), SpatOrbIdx_t([3, 3])], &
                  calculated => split_per_GAS(GAS_spec, SpatOrbIdx_t([3, 3])))
            call assert_true(all(expected(1) == calculated(1)))
            call assert_true(all(expected(2) == calculated(2)))
        end associate

        associate(expected => [SpatOrbIdx_t([1, 1, 2]), SpatOrbIdx_t([integer::])], &
                  calculated => split_per_GAS(GAS_spec, SpatOrbIdx_t([1, 1, 2])))
            call assert_true(all(expected(1) == calculated(1)))
            call assert_true(all(expected(2) == calculated(2)))
        end associate


        GAS_spec = GASSpec_t(n_orbs=[2, 4, 6], n_min=[2, 4, 6], n_max=[2, 4, 6])
        associate(expected => [SpatOrbIdx_t([1, 1, 2]), SpatOrbIdx_t([3]), SpatOrbIdx_t([5, 6])], &
                  calculated => split_per_GAS(GAS_spec, SpatOrbIdx_t([1, 1, 2, 3, 5, 6])))
            do i = 1, get_nGAS(GAS_spec)
                call assert_true(all(expected(i) == calculated(i)))
            end do
        end associate


    end subroutine

end module test_gasci_mod

program test_gasci_program

    use mpi
    use fruit
    use test_gasci_mod, only: test_igas_from_spatorb, test_igas_from_spinorb, &
        test_contains_det_spinorb, test_contains_det_spatorb, &
        test_particles_per_GAS_spatorb, test_particles_per_GAS_spinorb, &
        test_is_valid, test_is_connected, &
        test_get_possible_spaces_spinorb, test_get_possible_spaces_spatorb, &
        test_possible_holes, test_split_per_GAS


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
        call run_test_case(test_is_valid, "test_is_valid")
        call run_test_case(test_is_connected, "test_is_connected")
        call run_test_case(test_igas_from_spatorb, "test_igas_from_spatorb")
        call run_test_case(test_igas_from_spinorb, "test_igas_from_spinorb")
        call run_test_case(test_split_per_GAS, "test_split_per_GAS")
        call run_test_case(test_contains_det_spinorb, "test_contains_det_spinorb")
        call run_test_case(test_contains_det_spatorb, "test_contains_det_spatorb")
        call run_test_case(test_particles_per_GAS_spatorb, "test_particles_per_GAS_spatorb")
        call run_test_case(test_particles_per_GAS_spinorb, "test_particles_per_GAS_spinorb")
        call run_test_case(test_get_possible_spaces_spinorb, "test_get_possible_spaces_spinorb")
        call run_test_case(test_get_possible_spaces_spatorb, "test_get_possible_spaces_spatorb")
        call run_test_case(test_possible_holes, "test_possible_holes")
    end subroutine
end program test_gasci_program
