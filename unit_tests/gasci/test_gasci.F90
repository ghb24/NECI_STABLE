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
    use disconnected_gasci, only: init_disconnected_GAS, &
        gen_disconnected => generate_nGAS_excitation, clearGAS, &
        dyn_calc_pgen


    use gasci, only: GASSpec_t, get_iGAS, &
        contains_det, get_nGAS, particles_per_GAS, operator(.contains.), &
        is_valid, is_connected, get_possible_spaces, get_possible_holes, &
        split_per_GAS, generate_nGAS_excitation, &
        get_available_singles, get_available_doubles

    use sltcnd_mod, only: dyn_sltcnd_excit
    use unit_test_helper_excitgen, only: test_excitation_generator, &
        init_excitgen_test, finalize_excitgen_test, generate_random_integrals, &
        FciDumpWriter_t
    use unit_test_helpers, only: run_excit_gen_tester
    use DetBitOps, only: ilut_lt, ilut_gt
    implicit none
    private
    public :: test_igas_from_spatorb, test_igas_from_spinorb, &
        test_contains_det_spinorb, test_contains_det_spatorb, &
        test_particles_per_GAS_spatorb, test_particles_per_GAS_spinorb, &
        test_is_valid, test_is_connected, &
        test_get_possible_spaces_spinorb, test_get_possible_spaces_spatorb, &
        test_possible_holes, test_split_per_GAS, test_available, test_pgen



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
                  calculated => particles_per_GAS(split_per_GAS(&
                                        GAS_spec, SpinOrbIdx_t([1, 2, 5, 6]))))
            call assert_equals(expected, calculated, size(expected))
        end associate

        associate(expected => [1, 3], &
                  calculated => particles_per_GAS(split_per_GAS(&
                                        GAS_spec, SpinOrbIdx_t([1, 5, 6, 7]))))
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
        call assert_true(is_connected(GASSpec_t(n_orbs=[2, 4], n_min=[1, 4], n_max=[3, 4])))
        call assert_true(is_connected(GASSpec_t(n_orbs=[2, 4], n_min=[1, 5], n_max=[3, 5])))

        call assert_false(is_connected(GASSpec_t(n_orbs=[2, 4], n_min=[2, 4], n_max=[2, 4])))
        call assert_false(is_connected(GASSpec_t(n_orbs=[2, 4], n_min=[2, 5], n_max=[2, 5])))
    end subroutine

    subroutine test_available()
        type(SpinOrbIdx_t), allocatable :: expect_singles(:), expect_doubles(:)
        type(SpinOrbIdx_t), allocatable :: singles_exc_list(:), doubles_exc_list(:)
        type(SpinOrbIdx_t) :: det_I
        type(GASSpec_t) :: GAS_spec
        integer :: i, j

        block
            GAS_spec = GASSpec_t(&
                n_orbs=[2, 4], &
                n_min=[2, 4], &
                n_max=[2, 4])
            det_I = SpinOrbIdx_t([1, 2, 5, 6])
            call assert_true(is_valid(GAS_spec))
            call assert_true(GAS_spec .contains. det_I)

            expect_singles = [&
               SpinOrbIdx_t([1, 2, 5, 8]), SpinOrbIdx_t([1, 2, 6, 7]), &
               SpinOrbIdx_t([1, 4, 5, 6]), SpinOrbIdx_t([2, 3, 5, 6])]

            expect_doubles = [&
               SpinOrbIdx_t([1, 2, 7, 8]), SpinOrbIdx_t([1, 3, 6, 8]), &
               SpinOrbIdx_t([1, 4, 5, 8]), SpinOrbIdx_t([1, 4, 6, 7]), &
               SpinOrbIdx_t([2, 3, 5, 8]), SpinOrbIdx_t([2, 3, 6, 7]), &
               SpinOrbIdx_t([2, 4, 5, 7]), SpinOrbIdx_t([3, 4, 5, 6])]

            singles_exc_list = get_available_singles(GAS_spec, det_I)
            doubles_exc_list = get_available_doubles(GAS_spec, det_I)


            call assert_true(size(expect_singles) == size(singles_exc_list))
            do i = 1, size(expect_singles)
                call assert_true(all(expect_singles(i) == singles_exc_list(i)))
            end do

            call assert_true(size(expect_doubles) == size(doubles_exc_list))
            do i = 1, size(expect_doubles)
                call assert_true(all(expect_doubles(i) == doubles_exc_list(i)))
            end do
        end block

        block
            GAS_spec = GASSpec_t(&
                n_orbs=[2, 4], &
                n_min=[0, 4], &
                n_max=[4, 4])
            det_I = SpinOrbIdx_t([1, 2, 5, 6])
            call assert_true(is_valid(GAS_spec))
            call assert_true(GAS_spec .contains. det_I)

            expect_singles = [&
                   SpinOrbIdx_t([1, 2, 3, 6]), SpinOrbIdx_t([1, 2, 4, 5]), &
                   SpinOrbIdx_t([1, 2, 5, 8]), SpinOrbIdx_t([1, 2, 6, 7]), &
                   SpinOrbIdx_t([1, 4, 5, 6]), SpinOrbIdx_t([1, 5, 6, 8]), &
                   SpinOrbIdx_t([2, 3, 5, 6]), SpinOrbIdx_t([2, 5, 6, 7])]

            expect_doubles = [&
                   SpinOrbIdx_t([1, 2, 3, 4]), SpinOrbIdx_t([1, 2, 3, 8]), &
                   SpinOrbIdx_t([1, 2, 4, 7]), SpinOrbIdx_t([1, 2, 7, 8]), &
                   SpinOrbIdx_t([1, 3, 4, 6]), SpinOrbIdx_t([1, 3, 6, 8]), &
                   SpinOrbIdx_t([1, 4, 5, 8]), SpinOrbIdx_t([1, 4, 6, 7]), &
                   SpinOrbIdx_t([1, 6, 7, 8]), SpinOrbIdx_t([2, 3, 4, 5]), &
                   SpinOrbIdx_t([2, 3, 5, 8]), SpinOrbIdx_t([2, 3, 6, 7]), &
                   SpinOrbIdx_t([2, 4, 5, 7]), SpinOrbIdx_t([2, 5, 7, 8]), &
                   SpinOrbIdx_t([3, 4, 5, 6]), SpinOrbIdx_t([3, 5, 6, 8]), &
                   SpinOrbIdx_t([4, 5, 6, 7]), SpinOrbIdx_t([5, 6, 7, 8])]



            singles_exc_list = get_available_singles(GAS_spec, det_I)
            doubles_exc_list = get_available_doubles(GAS_spec, det_I)

            call assert_true(size(expect_singles) == size(singles_exc_list))
            do i = 1, size(expect_singles)
                call assert_true(all(expect_singles(i) == singles_exc_list(i)))
            end do


            call assert_true(size(expect_doubles) == size(doubles_exc_list))
            do i = 1, size(expect_doubles)
                call assert_true(all(expect_doubles(i) == doubles_exc_list(i)))
            end do
        end block

        block
            GAS_spec = GASSpec_t(&
                n_orbs=[6, 12], &
                n_min=[6, 12], &
                n_max=[6, 12])
            det_I = SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18])
            call assert_true(is_valid(GAS_spec))
            call assert_true(GAS_spec .contains. det_I)

            expect_singles = [&
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 14, 15, 16, 17, 18])]

            expect_doubles = [&
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 19, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 19, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 19, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 20, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 20, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 21, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 21, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 22, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 23, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 17, 20, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 17, 20, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 17, 22, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 18, 19, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 18, 19, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 18, 19, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 18, 20, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 18, 20, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 18, 21, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 18, 21, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 18, 22, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 18, 23, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 17, 19, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 17, 19, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 17, 19, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 17, 20, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 17, 20, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 17, 21, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 17, 21, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 17, 22, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 17, 23, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 18, 19, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 18, 19, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 16, 18, 21, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 17, 18, 19, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 17, 18, 19, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 17, 18, 19, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 17, 18, 20, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 17, 18, 20, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 17, 18, 21, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 17, 18, 21, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 17, 18, 22, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 17, 18, 23, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 17, 20, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 17, 20, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 17, 22, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 18, 19, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 18, 19, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 18, 19, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 18, 20, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 18, 20, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 18, 21, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 18, 21, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 18, 22, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 16, 18, 23, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 17, 18, 20, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 17, 18, 20, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 15, 17, 18, 22, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 16, 17, 18, 19, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 16, 17, 18, 19, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 16, 17, 18, 19, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 16, 17, 18, 20, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 16, 17, 18, 20, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 16, 17, 18, 21, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 16, 17, 18, 21, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 16, 17, 18, 22, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 16, 17, 18, 23, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 17, 19, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 17, 19, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 17, 19, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 17, 20, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 17, 20, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 17, 21, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 17, 21, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 17, 22, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 17, 23, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 18, 19, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 18, 19, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 16, 18, 21, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 17, 18, 19, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 17, 18, 19, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 17, 18, 19, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 17, 18, 20, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 17, 18, 20, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 17, 18, 21, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 17, 18, 21, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 17, 18, 22, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 15, 17, 18, 23, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 16, 17, 18, 19, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 16, 17, 18, 19, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 14, 16, 17, 18, 21, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 20, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 20, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 21, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 21, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 22, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 23, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 7, 13, 14, 15, 16, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 7, 13, 14, 15, 16, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 7, 13, 14, 15, 16, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 7, 13, 14, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 7, 13, 14, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 7, 13, 14, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 7, 14, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 7, 14, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 7, 14, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 8, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 9, 13, 14, 15, 16, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 9, 13, 14, 15, 16, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 9, 13, 14, 15, 16, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 9, 13, 14, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 9, 13, 14, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 9, 13, 14, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 9, 14, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 9, 14, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 9, 14, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 10, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 11, 13, 14, 15, 16, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 11, 13, 14, 15, 16, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 11, 13, 14, 15, 16, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 11, 13, 14, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 11, 13, 14, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 11, 13, 14, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 11, 14, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 11, 14, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 11, 14, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 5, 12, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 7, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 8, 13, 14, 15, 16, 17, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 8, 13, 14, 15, 16, 17, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 8, 13, 14, 15, 16, 17, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 8, 13, 14, 15, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 8, 13, 14, 15, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 8, 13, 14, 15, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 8, 13, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 8, 13, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 8, 13, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 9, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 10, 13, 14, 15, 16, 17, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 10, 13, 14, 15, 16, 17, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 10, 13, 14, 15, 16, 17, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 10, 13, 14, 15, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 10, 13, 14, 15, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 10, 13, 14, 15, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 10, 13, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 10, 13, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 10, 13, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 11, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 12, 13, 14, 15, 16, 17, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 12, 13, 14, 15, 16, 17, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 12, 13, 14, 15, 16, 17, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 12, 13, 14, 15, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 12, 13, 14, 15, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 12, 13, 14, 15, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 12, 13, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 12, 13, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 6, 12, 13, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 7, 8, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 7, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 7, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 8, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 8, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 9, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 9, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 10, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 4, 11, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 7, 13, 14, 15, 16, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 7, 13, 14, 15, 16, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 7, 13, 14, 15, 16, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 7, 13, 14, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 7, 13, 14, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 7, 13, 14, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 7, 14, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 7, 14, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 7, 14, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 8, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 9, 13, 14, 15, 16, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 9, 13, 14, 15, 16, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 9, 13, 14, 15, 16, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 9, 13, 14, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 9, 13, 14, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 9, 13, 14, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 9, 14, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 9, 14, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 9, 14, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 10, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 11, 13, 14, 15, 16, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 11, 13, 14, 15, 16, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 11, 13, 14, 15, 16, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 11, 13, 14, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 11, 13, 14, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 11, 13, 14, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 11, 14, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 11, 14, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 11, 14, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 6, 12, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 8, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 8, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 5, 10, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 6, 7, 8, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 6, 7, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 6, 7, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 6, 8, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 6, 8, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 6, 9, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 6, 9, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 6, 10, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 3, 6, 11, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 7, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 8, 13, 14, 15, 16, 17, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 8, 13, 14, 15, 16, 17, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 8, 13, 14, 15, 16, 17, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 8, 13, 14, 15, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 8, 13, 14, 15, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 8, 13, 14, 15, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 8, 13, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 8, 13, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 8, 13, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 9, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 10, 13, 14, 15, 16, 17, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 10, 13, 14, 15, 16, 17, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 10, 13, 14, 15, 16, 17, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 10, 13, 14, 15, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 10, 13, 14, 15, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 10, 13, 14, 15, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 10, 13, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 10, 13, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 10, 13, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 11, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 12, 13, 14, 15, 16, 17, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 12, 13, 14, 15, 16, 17, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 12, 13, 14, 15, 16, 17, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 12, 13, 14, 15, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 12, 13, 14, 15, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 12, 13, 14, 15, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 12, 13, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 12, 13, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 6, 12, 13, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 7, 8, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 7, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 7, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 8, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 8, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 9, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 9, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 10, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 5, 11, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 6, 7, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 6, 7, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 4, 6, 9, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 5, 6, 7, 8, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 5, 6, 7, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 5, 6, 7, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 5, 6, 8, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 5, 6, 8, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 5, 6, 9, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 5, 6, 9, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 5, 6, 10, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 2, 5, 6, 11, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 7, 13, 14, 15, 16, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 7, 13, 14, 15, 16, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 7, 13, 14, 15, 16, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 7, 13, 14, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 7, 13, 14, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 7, 13, 14, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 7, 14, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 7, 14, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 7, 14, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 8, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 9, 13, 14, 15, 16, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 9, 13, 14, 15, 16, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 9, 13, 14, 15, 16, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 9, 13, 14, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 9, 13, 14, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 9, 13, 14, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 9, 14, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 9, 14, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 9, 14, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 10, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 11, 13, 14, 15, 16, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 11, 13, 14, 15, 16, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 11, 13, 14, 15, 16, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 11, 13, 14, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 11, 13, 14, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 11, 13, 14, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 11, 14, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 11, 14, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 11, 14, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 6, 12, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 8, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 8, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 5, 10, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 6, 7, 8, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 6, 7, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 6, 7, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 6, 8, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 6, 8, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 6, 9, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 6, 9, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 6, 10, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 4, 6, 11, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 5, 6, 8, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 5, 6, 8, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 3, 5, 6, 10, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 4, 5, 6, 7, 8, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 4, 5, 6, 7, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 4, 5, 6, 7, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 4, 5, 6, 8, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 4, 5, 6, 8, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 4, 5, 6, 9, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 4, 5, 6, 9, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 4, 5, 6, 10, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([1, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 7, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 8, 13, 14, 15, 16, 17, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 8, 13, 14, 15, 16, 17, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 8, 13, 14, 15, 16, 17, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 8, 13, 14, 15, 17, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 8, 13, 14, 15, 17, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 8, 13, 14, 15, 17, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 8, 13, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 8, 13, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 8, 13, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 9, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 10, 13, 14, 15, 16, 17, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 10, 13, 14, 15, 16, 17, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 10, 13, 14, 15, 16, 17, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 10, 13, 14, 15, 17, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 10, 13, 14, 15, 17, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 10, 13, 14, 15, 17, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 10, 13, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 10, 13, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 10, 13, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 14, 15, 16, 17, 20]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 14, 15, 16, 17, 22]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 14, 15, 16, 17, 24]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 14, 15, 16, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 14, 15, 16, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 14, 15, 16, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 14, 15, 17, 18, 20]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 14, 15, 17, 18, 22]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 14, 15, 17, 18, 24]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 14, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 14, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 14, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 15, 16, 17, 18, 20]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 15, 16, 17, 18, 22]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 13, 15, 16, 17, 18, 24]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 14, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 14, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 11, 14, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 12, 13, 14, 15, 16, 17, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 12, 13, 14, 15, 16, 17, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 12, 13, 14, 15, 16, 17, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 12, 13, 14, 15, 17, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 12, 13, 14, 15, 17, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 12, 13, 14, 15, 17, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 12, 13, 15, 16, 17, 18, 19]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 12, 13, 15, 16, 17, 18, 21]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 6, 12, 13, 15, 16, 17, 18, 23]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 7, 8, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 7, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 7, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 8, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 8, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 9, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 9, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 10, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 6, 7, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 6, 7, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 4, 6, 9, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 5, 6, 7, 8, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 5, 6, 7, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 5, 6, 7, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 5, 6, 8, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 5, 6, 8, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 5, 6, 9, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 5, 6, 9, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 5, 6, 10, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 3, 5, 6, 11, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 4, 5, 6, 7, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 4, 5, 6, 7, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([2, 4, 5, 6, 9, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([3, 4, 5, 6, 7, 8, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([3, 4, 5, 6, 7, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([3, 4, 5, 6, 7, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([3, 4, 5, 6, 8, 9, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([3, 4, 5, 6, 8, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([3, 4, 5, 6, 9, 10, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([3, 4, 5, 6, 9, 12, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([3, 4, 5, 6, 10, 11, 13, 14, 15, 16, 17, 18]), &
                    SpinOrbIdx_t([3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 18])]

            singles_exc_list = get_available_singles(GAS_spec, det_I)
            doubles_exc_list = get_available_doubles(GAS_spec, det_I)

            call assert_true(size(expect_singles) == size(singles_exc_list))
            do i = 1, size(expect_singles)
                call assert_true(all(expect_singles(i) == singles_exc_list(i)))
            end do

            call assert_true(size(expect_doubles) == size(doubles_exc_list))
            do i = 1, size(expect_doubles)
                call assert_true(all(expect_doubles(i) == doubles_exc_list(i)))
            end do
        end block
    end subroutine


    subroutine test_pgen()
        use gasci, only: global_GAS_spec => GAS_specification
        use SystemData, only: tGASSpinRecoupling
        use FciMCData, only: pSingles, pDoubles, pParallel
        type(GASSpec_t) :: GAS_spec
        type(SpinOrbIdx_t) :: det_I

        logical :: successful
        integer, parameter :: n_iters=10**5

        pParallel = 0.05_dp
        pSingles = 0.3_dp
        pDoubles = 1.0_dp - pSingles

        call assert_true(tGASSpinRecoupling)
        det_I = SpinOrbIdx_t([1, 2, 3, 7, 8, 10])

        GAS_spec = GASSpec_t(&
            n_orbs=[3, 6], &
            n_min=[3, size(det_I)], &
            n_max=[3, size(det_I)])
        call assert_true(is_valid(GAS_spec))
        call assert_true(GAS_spec .contains. det_I)
        global_GAS_spec = GAS_spec

        call init_excitgen_test(size(det_I), FciDumpWriter_t(write_small_Li2_FCIDUMP, 'FCIDUMP'))
        call run_excit_gen_tester( &
            generate_nGAS_excitation, 'general implementation, disconnected, Li2', &
            opt_nI=det_I%idx, &
            opt_n_iters=n_iters, &
            gen_all_excits=gen_all_excits, &
            successful=successful)
        call assert_true(successful)
        call finalize_excitgen_test()

        call init_excitgen_test(size(det_I), FciDumpWriter_t(write_small_Li2_FCIDUMP, 'FCIDUMP'))
        call init_disconnected_GAS(GAS_spec)
        call run_excit_gen_tester( &
            gen_disconnected, 'only disconnected implementation, disconnected, Li2', &
            opt_nI=det_I%idx, opt_n_iters=n_iters, &
            gen_all_excits=gen_all_excits, &
            calc_pgen=dyn_calc_pgen, &
            successful=successful)
        call assert_true(successful)
        call clearGAS()
        call finalize_excitgen_test()
    contains
        subroutine gen_all_excits(nI, n_excits, det_list)
            integer, intent(in) :: nI(nel)
            integer, intent(out) :: n_excits
            integer(n_int), intent(out), allocatable :: det_list(:,:)

            type(SpinOrbIdx_t) :: det_I
            type(SpinOrbIdx_t), allocatable :: singles(:), doubles(:)
            integer :: i, j, k

            det_I = SpinOrbIdx_t(nI)

            singles = get_available_singles(GAS_spec, det_I)
            doubles = get_available_doubles(GAS_spec, det_I)

            n_excits = size(singles) + size(doubles)
            allocate(det_list(0:niftot, n_excits))
            j = 1
            do i = 1, size(singles)
                det_list(:, j) = to_ilut(singles(i))
                j = j + 1
            end do

            do i = 1, size(doubles)
                det_list(:, j) = to_ilut(doubles(i))
                j = j + 1
            end do

            call sort(det_list, ilut_lt, ilut_gt)
        end subroutine gen_all_excits

        subroutine random_fcidump(iunit)
            integer, intent(in) :: iunit
            call generate_random_integrals(&
                iunit, n_el=size(det_I), n_spat_orb=GAS_spec%n_orbs(get_nGAS(GAS_spec)), &
                sparse=1.0_dp, sparseT=1.0_dp, total_ms=sum(calc_spin(det_I)))
        end subroutine

        subroutine write_large_Li2_FCIDUMP(unit_id)
            integer, intent(in) :: unit_id
            write(unit_id, '(A)') '  &FCI NORB= 10,NELEC=  6,MS2=  0,'
            write(unit_id, '(A)') '  ORBSYM= 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,'
            write(unit_id, '(A)') '  ISYM=0'
            write(unit_id, '(A)') ' &END'
            write(unit_id, '(A)') '     1.6327627870        1    1    1    1'
            write(unit_id, '(A)') '   -0.13816660653        2    1    1    1'
            write(unit_id, '(A)') '    0.17265361682E-01    2    1    2    1'
            write(unit_id, '(A)') '    0.34397810102        2    2    1    1'
            write(unit_id, '(A)') '   -0.53019371999E-02    2    2    2    1'
            write(unit_id, '(A)') '    0.25068860113        2    2    2    2'
            write(unit_id, '(A)') '    0.42029080786E-02    3    1    3    1'
            write(unit_id, '(A)') '    0.63322371809E-02    3    2    3    1'
            write(unit_id, '(A)') '    0.50916769052E-01    3    2    3    2'
            write(unit_id, '(A)') '    0.29004307384        3    3    1    1'
            write(unit_id, '(A)') '   -0.35754113716E-02    3    3    2    1'
            write(unit_id, '(A)') '    0.22795394153        3    3    2    2'
            write(unit_id, '(A)') '    0.22612217356        3    3    3    3'
            write(unit_id, '(A)') '    0.42029080786E-02    4    1    4    1'
            write(unit_id, '(A)') '    0.63322371809E-02    4    2    4    1'
            write(unit_id, '(A)') '    0.50916769052E-01    4    2    4    2'
            write(unit_id, '(A)') '    0.12372667071E-01    4    3    4    3'
            write(unit_id, '(A)') '    0.29004307384        4    4    1    1'
            write(unit_id, '(A)') '   -0.35754113716E-02    4    4    2    1'
            write(unit_id, '(A)') '    0.22795394153        4    4    2    2'
            write(unit_id, '(A)') '    0.20137683942        4    4    3    3'
            write(unit_id, '(A)') '    0.22612217356        4    4    4    4'
            write(unit_id, '(A)') '    0.23510693561E-01    5    1    1    1'
            write(unit_id, '(A)') '   -0.40301071636E-02    5    1    2    1'
            write(unit_id, '(A)') '    0.19240472786E-02    5    1    2    2'
            write(unit_id, '(A)') '    0.24144100660E-02    5    1    3    3'
            write(unit_id, '(A)') '    0.24144100660E-02    5    1    4    4'
            write(unit_id, '(A)') '    0.57751228916E-02    5    1    5    1'
            write(unit_id, '(A)') '   -0.24847734942E-01    5    2    1    1'
            write(unit_id, '(A)') '   -0.82183473553E-03    5    2    2    1'
            write(unit_id, '(A)') '   -0.80688684698E-02    5    2    2    2'
            write(unit_id, '(A)') '   -0.75495433303E-03    5    2    3    3'
            write(unit_id, '(A)') '   -0.75495433303E-03    5    2    4    4'
            write(unit_id, '(A)') '    0.43308083154E-02    5    2    5    1'
            write(unit_id, '(A)') '    0.29337139736E-01    5    2    5    2'
            write(unit_id, '(A)') '   -0.26878675252E-03    5    3    3    1'
            write(unit_id, '(A)') '    0.19061411205E-02    5    3    3    2'
            write(unit_id, '(A)') '    0.11681694683E-01    5    3    5    3'
            write(unit_id, '(A)') '   -0.26878675252E-03    5    4    4    1'
            write(unit_id, '(A)') '    0.19061411205E-02    5    4    4    2'
            write(unit_id, '(A)') '    0.11681694683E-01    5    4    5    4'
            write(unit_id, '(A)') '    0.30377515269        5    5    1    1'
            write(unit_id, '(A)') '   -0.51275560816E-02    5    5    2    1'
            write(unit_id, '(A)') '    0.22083460664        5    5    2    2'
            write(unit_id, '(A)') '    0.20326768821        5    5    3    3'
            write(unit_id, '(A)') '    0.20326768821        5    5    4    4'
            write(unit_id, '(A)') '    0.28444005389E-02    5    5    5    1'
            write(unit_id, '(A)') '    0.10673539613E-01    5    5    5    2'
            write(unit_id, '(A)') '    0.23222955692        5    5    5    5'
            write(unit_id, '(A)') '    0.75043896954E-02    6    1    1    1'
            write(unit_id, '(A)') '   -0.10585776208E-02    6    1    2    1'
            write(unit_id, '(A)') '    0.19681432522E-03    6    1    2    2'
            write(unit_id, '(A)') '    0.12946913459E-03    6    1    3    3'
            write(unit_id, '(A)') '    0.12946913459E-03    6    1    4    4'
            write(unit_id, '(A)') '    0.63388586124E-03    6    1    5    1'
            write(unit_id, '(A)') '    0.39256197191E-03    6    1    5    2'
            write(unit_id, '(A)') '    0.21635429134E-03    6    1    5    5'
            write(unit_id, '(A)') '    0.19848698782E-03    6    1    6    1'
            write(unit_id, '(A)') '   -0.97992463811E-02    6    2    1    1'
            write(unit_id, '(A)') '    0.41363792273E-03    6    2    2    1'
            write(unit_id, '(A)') '   -0.41124456245E-02    6    2    2    2'
            write(unit_id, '(A)') '   -0.40070999477E-02    6    2    3    3'
            write(unit_id, '(A)') '   -0.40070999477E-02    6    2    4    4'
            write(unit_id, '(A)') '    0.14415017417E-03    6    2    5    1'
            write(unit_id, '(A)') '    0.95133002289E-03    6    2    5    2'
            write(unit_id, '(A)') '   -0.24134011285E-02    6    2    5    5'
            write(unit_id, '(A)') '    0.58481893008E-03    6    2    6    1'
            write(unit_id, '(A)') '    0.44589483716E-02    6    2    6    2'
            write(unit_id, '(A)') '   -0.32605732421E-03    6    3    3    1'
            write(unit_id, '(A)') '   -0.23156361648E-02    6    3    3    2'
            write(unit_id, '(A)') '    0.79370563794E-03    6    3    5    3'
            write(unit_id, '(A)') '    0.32058144665E-03    6    3    6    3'
            write(unit_id, '(A)') '   -0.32605732421E-03    6    4    4    1'
            write(unit_id, '(A)') '   -0.23156361648E-02    6    4    4    2'
            write(unit_id, '(A)') '    0.79370563794E-03    6    4    5    4'
            write(unit_id, '(A)') '    0.32058144665E-03    6    4    6    4'
            write(unit_id, '(A)') '    0.85194148321E-02    6    5    1    1'
            write(unit_id, '(A)') '   -0.64975599013E-03    6    5    2    1'
            write(unit_id, '(A)') '    0.25818484021E-02    6    5    2    2'
            write(unit_id, '(A)') '    0.32449377861E-02    6    5    3    3'
            write(unit_id, '(A)') '    0.32449377861E-02    6    5    4    4'
            write(unit_id, '(A)') '    0.13754341276E-03    6    5    5    1'
            write(unit_id, '(A)') '    0.22918823421E-02    6    5    5    2'
            write(unit_id, '(A)') '    0.42153705731E-02    6    5    5    5'
            write(unit_id, '(A)') '   -0.32271116318E-03    6    5    6    1'
            write(unit_id, '(A)') '   -0.36448704566E-02    6    5    6    2'
            write(unit_id, '(A)') '    0.42347146970E-02    6    5    6    5'
            write(unit_id, '(A)') '    0.20051849601        6    6    1    1'
            write(unit_id, '(A)') '    0.41735324022E-02    6    6    2    1'
            write(unit_id, '(A)') '    0.18878583186        6    6    2    2'
            write(unit_id, '(A)') '    0.15749127301        6    6    3    3'
            write(unit_id, '(A)') '    0.15749127301        6    6    4    4'
            write(unit_id, '(A)') '   -0.12885230094E-02    6    6    5    1'
            write(unit_id, '(A)') '   -0.35609078534E-01    6    6    5    2'
            write(unit_id, '(A)') '    0.16750323527        6    6    5    5'
            write(unit_id, '(A)') '    0.75043896954E-02    6    6    6    1'
            write(unit_id, '(A)') '    0.66793041000E-01    6    6    6    2'
            write(unit_id, '(A)') '   -0.69455147750E-01    6    6    6    5'
            write(unit_id, '(A)') '     1.6327627870        6    6    6    6'
            write(unit_id, '(A)') '    0.66793041000E-01    7    1    1    1'
            write(unit_id, '(A)') '   -0.71513328264E-02    7    1    2    1'
            write(unit_id, '(A)') '    0.65121751160E-02    7    1    2    2'
            write(unit_id, '(A)') '    0.53241021228E-02    7    1    3    3'
            write(unit_id, '(A)') '    0.53241021228E-02    7    1    4    4'
            write(unit_id, '(A)') '    0.32575861986E-02    7    1    5    1'
            write(unit_id, '(A)') '    0.18623498631E-02    7    1    5    2'
            write(unit_id, '(A)') '    0.62023397312E-02    7    1    5    5'
            write(unit_id, '(A)') '    0.58481893008E-03    7    1    6    1'
            write(unit_id, '(A)') '   -0.53514496182E-03    7    1    6    2'
            write(unit_id, '(A)') '    0.10712954171E-02    7    1    6    5'
            write(unit_id, '(A)') '   -0.97992463811E-02    7    1    6    6'
            write(unit_id, '(A)') '    0.44589483716E-02    7    1    7    1'
            write(unit_id, '(A)') '   -0.24959346753E-01    7    2    1    1'
            write(unit_id, '(A)') '    0.79717594725E-03    7    2    2    1'
            write(unit_id, '(A)') '   -0.10346561190E-01    7    2    2    2'
            write(unit_id, '(A)') '   -0.78447694255E-02    7    2    3    3'
            write(unit_id, '(A)') '   -0.78447694255E-02    7    2    4    4'
            write(unit_id, '(A)') '    0.38311928397E-03    7    2    5    1'
            write(unit_id, '(A)') '    0.54438691459E-02    7    2    5    2'
            write(unit_id, '(A)') '   -0.40602638310E-02    7    2    5    5'
            write(unit_id, '(A)') '    0.21941615579E-04    7    2    6    1'
            write(unit_id, '(A)') '    0.28662672467E-04    7    2    6    2'
            write(unit_id, '(A)') '    0.10642128564E-02    7    2    6    5'
            write(unit_id, '(A)') '   -0.24959346753E-01    7    2    6    6'
            write(unit_id, '(A)') '    0.28662672467E-04    7    2    7    1'
            write(unit_id, '(A)') '    0.32376065304E-02    7    2    7    2'
            write(unit_id, '(A)') '   -0.21173245038E-03    7    3    3    1'
            write(unit_id, '(A)') '    0.38589709967E-03    7    3    3    2'
            write(unit_id, '(A)') '    0.13620273913E-02    7    3    5    3'
            write(unit_id, '(A)') '    0.16799583491E-03    7    3    6    3'
            write(unit_id, '(A)') '    0.12651804505E-02    7    3    7    3'
            write(unit_id, '(A)') '   -0.21173245038E-03    7    4    4    1'
            write(unit_id, '(A)') '    0.38589709967E-03    7    4    4    2'
            write(unit_id, '(A)') '    0.13620273913E-02    7    4    5    4'
            write(unit_id, '(A)') '    0.16799583491E-03    7    4    6    4'
            write(unit_id, '(A)') '    0.12651804505E-02    7    4    7    4'
            write(unit_id, '(A)') '    0.25571416458E-01    7    5    1    1'
            write(unit_id, '(A)') '   -0.21686924900E-02    7    5    2    1'
            write(unit_id, '(A)') '    0.69314080490E-02    7    5    2    2'
            write(unit_id, '(A)') '    0.83196872746E-02    7    5    3    3'
            write(unit_id, '(A)') '    0.83196872746E-02    7    5    4    4'
            write(unit_id, '(A)') '    0.15138436977E-02    7    5    5    1'
            write(unit_id, '(A)') '    0.12133194409E-01    7    5    5    2'
            write(unit_id, '(A)') '    0.16215012380E-01    7    5    5    5'
            write(unit_id, '(A)') '    0.26485037615E-03    7    5    6    1'
            write(unit_id, '(A)') '   -0.94215484865E-03    7    5    6    2'
            write(unit_id, '(A)') '    0.35601423246E-02    7    5    6    5'
            write(unit_id, '(A)') '   -0.36032118080E-01    7    5    6    6'
            write(unit_id, '(A)') '    0.31295027157E-02    7    5    7    1'
            write(unit_id, '(A)') '    0.20173358648E-02    7    5    7    2'
            write(unit_id, '(A)') '    0.12047629067E-01    7    5    7    5'
            write(unit_id, '(A)') '    0.41735324022E-02    7    6    1    1'
            write(unit_id, '(A)') '   -0.36901184632E-03    7    6    2    1'
            write(unit_id, '(A)') '    0.91458714190E-03    7    6    2    2'
            write(unit_id, '(A)') '    0.15869063848E-02    7    6    3    3'
            write(unit_id, '(A)') '    0.15869063848E-02    7    6    4    4'
            write(unit_id, '(A)') '    0.15611875539E-03    7    6    5    1'
            write(unit_id, '(A)') '    0.99885138838E-03    7    6    5    2'
            write(unit_id, '(A)') '    0.17217164860E-02    7    6    5    5'
            write(unit_id, '(A)') '   -0.10585776208E-02    7    6    6    1'
            write(unit_id, '(A)') '   -0.71513328264E-02    7    6    6    2'
            write(unit_id, '(A)') '    0.54753226077E-02    7    6    6    5'
            write(unit_id, '(A)') '   -0.13816660653        7    6    6    6'
            write(unit_id, '(A)') '    0.41363792273E-03    7    6    7    1'
            write(unit_id, '(A)') '    0.79717594725E-03    7    6    7    2'
            write(unit_id, '(A)') '    0.51169265382E-03    7    6    7    5'
            write(unit_id, '(A)') '    0.17265361682E-01    7    6    7    6'
            write(unit_id, '(A)') '    0.18878583186        7    7    1    1'
            write(unit_id, '(A)') '    0.91458714190E-03    7    7    2    1'
            write(unit_id, '(A)') '    0.16743661593        7    7    2    2'
            write(unit_id, '(A)') '    0.14955029106        7    7    3    3'
            write(unit_id, '(A)') '    0.14955029106        7    7    4    4'
            write(unit_id, '(A)') '    0.10855735869E-03    7    7    5    1'
            write(unit_id, '(A)') '   -0.22151564675E-01    7    7    5    2'
            write(unit_id, '(A)') '    0.15193566117        7    7    5    5'
            write(unit_id, '(A)') '    0.19681432522E-03    7    7    6    1'
            write(unit_id, '(A)') '    0.65121751160E-02    7    7    6    2'
            write(unit_id, '(A)') '   -0.12653081295E-01    7    7    6    5'
            write(unit_id, '(A)') '    0.34397810102        7    7    6    6'
            write(unit_id, '(A)') '   -0.41124456245E-02    7    7    7    1'
            write(unit_id, '(A)') '   -0.10346561190E-01    7    7    7    2'
            write(unit_id, '(A)') '   -0.25032329445E-01    7    7    7    5'
            write(unit_id, '(A)') '   -0.53019371999E-02    7    7    7    6'
            write(unit_id, '(A)') '    0.25068860113        7    7    7    7'
            write(unit_id, '(A)') '   -0.82333190944E-03    8    1    3    1'
            write(unit_id, '(A)') '   -0.18412795676E-02    8    1    3    2'
            write(unit_id, '(A)') '    0.28466432004E-03    8    1    5    3'
            write(unit_id, '(A)') '    0.14623348073E-03    8    1    6    3'
            write(unit_id, '(A)') '   -0.88597661874E-04    8    1    7    3'
            write(unit_id, '(A)') '    0.32058144665E-03    8    1    8    1'
            write(unit_id, '(A)') '   -0.35618532176E-03    8    2    3    1'
            write(unit_id, '(A)') '   -0.75979938319E-03    8    2    3    2'
            write(unit_id, '(A)') '   -0.18467178293E-02    8    2    5    3'
            write(unit_id, '(A)') '   -0.88597661874E-04    8    2    6    3'
            write(unit_id, '(A)') '    0.74379496651E-04    8    2    7    3'
            write(unit_id, '(A)') '    0.16799583491E-03    8    2    8    1'
            write(unit_id, '(A)') '    0.12651804505E-02    8    2    8    2'
            write(unit_id, '(A)') '   -0.79810529730E-02    8    3    1    1'
            write(unit_id, '(A)') '    0.40240805622E-03    8    3    2    1'
            write(unit_id, '(A)') '   -0.17890067998E-02    8    3    2    2'
            write(unit_id, '(A)') '   -0.10096059220E-02    8    3    3    3'
            write(unit_id, '(A)') '   -0.13256212443E-02    8    3    4    4'
            write(unit_id, '(A)') '   -0.51972466625E-03    8    3    5    1'
            write(unit_id, '(A)') '   -0.37198110259E-02    8    3    5    2'
            write(unit_id, '(A)') '   -0.29729353954E-02    8    3    5    5'
            write(unit_id, '(A)') '   -0.14350926511E-03    8    3    6    1'
            write(unit_id, '(A)') '   -0.77434722556E-03    8    3    6    2'
            write(unit_id, '(A)') '    0.42692007067E-03    8    3    6    5'
            write(unit_id, '(A)') '   -0.79810529730E-02    8    3    6    6'
            write(unit_id, '(A)') '   -0.77434722556E-03    8    3    7    1'
            write(unit_id, '(A)') '   -0.32701342776E-03    8    3    7    2'
            write(unit_id, '(A)') '   -0.17943412170E-02    8    3    7    5'
            write(unit_id, '(A)') '    0.40240805622E-03    8    3    7    6'
            write(unit_id, '(A)') '   -0.17890067998E-02    8    3    7    7'
            write(unit_id, '(A)') '    0.26045402532E-02    8    3    8    3'
            write(unit_id, '(A)') '    0.15800766114E-03    8    4    4    3'
            write(unit_id, '(A)') '    0.49796249053E-03    8    4    8    4'
            write(unit_id, '(A)') '   -0.67569386042E-03    8    5    3    1'
            write(unit_id, '(A)') '   -0.86875942289E-02    8    5    3    2'
            write(unit_id, '(A)') '    0.63085162628E-03    8    5    5    3'
            write(unit_id, '(A)') '    0.66296134760E-03    8    5    6    3'
            write(unit_id, '(A)') '   -0.62910978514E-03    8    5    7    3'
            write(unit_id, '(A)') '    0.72390456503E-03    8    5    8    1'
            write(unit_id, '(A)') '   -0.44464563393E-03    8    5    8    2'
            write(unit_id, '(A)') '    0.42563060127E-02    8    5    8    5'
            write(unit_id, '(A)') '    0.52578962662E-04    8    6    3    1'
            write(unit_id, '(A)') '    0.61941348667E-03    8    6    3    2'
            write(unit_id, '(A)') '   -0.27560857734E-03    8    6    5    3'
            write(unit_id, '(A)') '   -0.82333190944E-03    8    6    6    3'
            write(unit_id, '(A)') '   -0.35618532176E-03    8    6    7    3'
            write(unit_id, '(A)') '   -0.32605732421E-03    8    6    8    1'
            write(unit_id, '(A)') '   -0.21173245038E-03    8    6    8    2'
            write(unit_id, '(A)') '   -0.14673774416E-02    8    6    8    5'
            write(unit_id, '(A)') '    0.42029080786E-02    8    6    8    6'
            write(unit_id, '(A)') '    0.61941348667E-03    8    7    3    1'
            write(unit_id, '(A)') '    0.12480416371E-01    8    7    3    2'
            write(unit_id, '(A)') '   -0.42872299709E-02    8    7    5    3'
            write(unit_id, '(A)') '   -0.18412795676E-02    8    7    6    3'
            write(unit_id, '(A)') '   -0.75979938319E-03    8    7    7    3'
            write(unit_id, '(A)') '   -0.23156361648E-02    8    7    8    1'
            write(unit_id, '(A)') '    0.38589709967E-03    8    7    8    2'
            write(unit_id, '(A)') '   -0.12071158096E-01    8    7    8    5'
            write(unit_id, '(A)') '    0.63322371809E-02    8    7    8    6'
            write(unit_id, '(A)') '    0.50916769052E-01    8    7    8    7'
            write(unit_id, '(A)') '    0.15749127301        8    8    1    1'
            write(unit_id, '(A)') '    0.15869063848E-02    8    8    2    1'
            write(unit_id, '(A)') '    0.14955029106        8    8    2    2'
            write(unit_id, '(A)') '    0.13770968110        8    8    3    3'
            write(unit_id, '(A)') '    0.13284148553        8    8    4    4'
            write(unit_id, '(A)') '   -0.25169554469E-03    8    8    5    1'
            write(unit_id, '(A)') '   -0.20525696004E-01    8    8    5    2'
            write(unit_id, '(A)') '    0.13266789951        8    8    5    5'
            write(unit_id, '(A)') '    0.12946913459E-03    8    8    6    1'
            write(unit_id, '(A)') '    0.53241021228E-02    8    8    6    2'
            write(unit_id, '(A)') '   -0.10674136475E-01    8    8    6    5'
            write(unit_id, '(A)') '    0.29004307384        8    8    6    6'
            write(unit_id, '(A)') '   -0.40070999477E-02    8    8    7    1'
            write(unit_id, '(A)') '   -0.78447694255E-02    8    8    7    2'
            write(unit_id, '(A)') '   -0.23827492755E-01    8    8    7    5'
            write(unit_id, '(A)') '   -0.35754113716E-02    8    8    7    6'
            write(unit_id, '(A)') '    0.22795394153        8    8    7    7'
            write(unit_id, '(A)') '   -0.10096059220E-02    8    8    8    3'
            write(unit_id, '(A)') '    0.22612217356        8    8    8    8'
            write(unit_id, '(A)') '   -0.82333190944E-03    9    1    4    1'
            write(unit_id, '(A)') '   -0.18412795676E-02    9    1    4    2'
            write(unit_id, '(A)') '    0.28466432004E-03    9    1    5    4'
            write(unit_id, '(A)') '    0.14623348073E-03    9    1    6    4'
            write(unit_id, '(A)') '   -0.88597661874E-04    9    1    7    4'
            write(unit_id, '(A)') '    0.32058144665E-03    9    1    9    1'
            write(unit_id, '(A)') '   -0.35618532176E-03    9    2    4    1'
            write(unit_id, '(A)') '   -0.75979938319E-03    9    2    4    2'
            write(unit_id, '(A)') '   -0.18467178293E-02    9    2    5    4'
            write(unit_id, '(A)') '   -0.88597661874E-04    9    2    6    4'
            write(unit_id, '(A)') '    0.74379496651E-04    9    2    7    4'
            write(unit_id, '(A)') '    0.16799583491E-03    9    2    9    1'
            write(unit_id, '(A)') '    0.12651804505E-02    9    2    9    2'
            write(unit_id, '(A)') '    0.15800766114E-03    9    3    4    3'
            write(unit_id, '(A)') '    0.49796249053E-03    9    3    8    4'
            write(unit_id, '(A)') '    0.49796249053E-03    9    3    9    3'
            write(unit_id, '(A)') '   -0.79810529730E-02    9    4    1    1'
            write(unit_id, '(A)') '    0.40240805622E-03    9    4    2    1'
            write(unit_id, '(A)') '   -0.17890067998E-02    9    4    2    2'
            write(unit_id, '(A)') '   -0.13256212443E-02    9    4    3    3'
            write(unit_id, '(A)') '   -0.10096059220E-02    9    4    4    4'
            write(unit_id, '(A)') '   -0.51972466625E-03    9    4    5    1'
            write(unit_id, '(A)') '   -0.37198110259E-02    9    4    5    2'
            write(unit_id, '(A)') '   -0.29729353954E-02    9    4    5    5'
            write(unit_id, '(A)') '   -0.14350926511E-03    9    4    6    1'
            write(unit_id, '(A)') '   -0.77434722556E-03    9    4    6    2'
            write(unit_id, '(A)') '    0.42692007067E-03    9    4    6    5'
            write(unit_id, '(A)') '   -0.79810529730E-02    9    4    6    6'
            write(unit_id, '(A)') '   -0.77434722556E-03    9    4    7    1'
            write(unit_id, '(A)') '   -0.32701342776E-03    9    4    7    2'
            write(unit_id, '(A)') '   -0.17943412170E-02    9    4    7    5'
            write(unit_id, '(A)') '    0.40240805622E-03    9    4    7    6'
            write(unit_id, '(A)') '   -0.17890067998E-02    9    4    7    7'
            write(unit_id, '(A)') '    0.16086152722E-02    9    4    8    3'
            write(unit_id, '(A)') '   -0.13256212443E-02    9    4    8    8'
            write(unit_id, '(A)') '    0.26045402532E-02    9    4    9    4'
            write(unit_id, '(A)') '   -0.67569386042E-03    9    5    4    1'
            write(unit_id, '(A)') '   -0.86875942289E-02    9    5    4    2'
            write(unit_id, '(A)') '    0.63085162628E-03    9    5    5    4'
            write(unit_id, '(A)') '    0.66296134760E-03    9    5    6    4'
            write(unit_id, '(A)') '   -0.62910978514E-03    9    5    7    4'
            write(unit_id, '(A)') '    0.72390456503E-03    9    5    9    1'
            write(unit_id, '(A)') '   -0.44464563393E-03    9    5    9    2'
            write(unit_id, '(A)') '    0.42563060127E-02    9    5    9    5'
            write(unit_id, '(A)') '    0.52578962662E-04    9    6    4    1'
            write(unit_id, '(A)') '    0.61941348667E-03    9    6    4    2'
            write(unit_id, '(A)') '   -0.27560857734E-03    9    6    5    4'
            write(unit_id, '(A)') '   -0.82333190944E-03    9    6    6    4'
            write(unit_id, '(A)') '   -0.35618532176E-03    9    6    7    4'
            write(unit_id, '(A)') '   -0.32605732421E-03    9    6    9    1'
            write(unit_id, '(A)') '   -0.21173245038E-03    9    6    9    2'
            write(unit_id, '(A)') '   -0.14673774416E-02    9    6    9    5'
            write(unit_id, '(A)') '    0.42029080786E-02    9    6    9    6'
            write(unit_id, '(A)') '    0.61941348667E-03    9    7    4    1'
            write(unit_id, '(A)') '    0.12480416371E-01    9    7    4    2'
            write(unit_id, '(A)') '   -0.42872299709E-02    9    7    5    4'
            write(unit_id, '(A)') '   -0.18412795676E-02    9    7    6    4'
            write(unit_id, '(A)') '   -0.75979938319E-03    9    7    7    4'
            write(unit_id, '(A)') '   -0.23156361648E-02    9    7    9    1'
            write(unit_id, '(A)') '    0.38589709967E-03    9    7    9    2'
            write(unit_id, '(A)') '   -0.12071158096E-01    9    7    9    5'
            write(unit_id, '(A)') '    0.63322371809E-02    9    7    9    6'
            write(unit_id, '(A)') '    0.50916769052E-01    9    7    9    7'
            write(unit_id, '(A)') '    0.24340977871E-02    9    8    4    3'
            write(unit_id, '(A)') '    0.15800766114E-03    9    8    8    4'
            write(unit_id, '(A)') '    0.15800766114E-03    9    8    9    3'
            write(unit_id, '(A)') '    0.12372667071E-01    9    8    9    8'
            write(unit_id, '(A)') '    0.15749127301        9    9    1    1'
            write(unit_id, '(A)') '    0.15869063848E-02    9    9    2    1'
            write(unit_id, '(A)') '    0.14955029106        9    9    2    2'
            write(unit_id, '(A)') '    0.13284148553        9    9    3    3'
            write(unit_id, '(A)') '    0.13770968110        9    9    4    4'
            write(unit_id, '(A)') '   -0.25169554469E-03    9    9    5    1'
            write(unit_id, '(A)') '   -0.20525696004E-01    9    9    5    2'
            write(unit_id, '(A)') '    0.13266789951        9    9    5    5'
            write(unit_id, '(A)') '    0.12946913459E-03    9    9    6    1'
            write(unit_id, '(A)') '    0.53241021228E-02    9    9    6    2'
            write(unit_id, '(A)') '   -0.10674136475E-01    9    9    6    5'
            write(unit_id, '(A)') '    0.29004307384        9    9    6    6'
            write(unit_id, '(A)') '   -0.40070999477E-02    9    9    7    1'
            write(unit_id, '(A)') '   -0.78447694255E-02    9    9    7    2'
            write(unit_id, '(A)') '   -0.23827492755E-01    9    9    7    5'
            write(unit_id, '(A)') '   -0.35754113716E-02    9    9    7    6'
            write(unit_id, '(A)') '    0.22795394153        9    9    7    7'
            write(unit_id, '(A)') '   -0.13256212443E-02    9    9    8    3'
            write(unit_id, '(A)') '    0.20137683942        9    9    8    8'
            write(unit_id, '(A)') '   -0.10096059220E-02    9    9    9    4'
            write(unit_id, '(A)') '    0.22612217356        9    9    9    9'
            write(unit_id, '(A)') '    0.69455147750E-01   10    1    1    1'
            write(unit_id, '(A)') '   -0.54753226077E-02   10    1    2    1'
            write(unit_id, '(A)') '    0.12653081295E-01   10    1    2    2'
            write(unit_id, '(A)') '    0.10674136475E-01   10    1    3    3'
            write(unit_id, '(A)') '    0.10674136475E-01   10    1    4    4'
            write(unit_id, '(A)') '    0.15281018528E-02   10    1    5    1'
            write(unit_id, '(A)') '    0.14173834149E-03   10    1    5    2'
            write(unit_id, '(A)') '    0.11468694313E-01   10    1    5    5'
            write(unit_id, '(A)') '    0.32271116318E-03   10    1    6    1'
            write(unit_id, '(A)') '   -0.10712954171E-02   10    1    6    2'
            write(unit_id, '(A)') '    0.14450359076E-02   10    1    6    5'
            write(unit_id, '(A)') '   -0.85194148321E-02   10    1    6    6'
            write(unit_id, '(A)') '    0.36448704566E-02   10    1    7    1'
            write(unit_id, '(A)') '   -0.10642128564E-02   10    1    7    2'
            write(unit_id, '(A)') '    0.34494177544E-02   10    1    7    5'
            write(unit_id, '(A)') '    0.64975599013E-03   10    1    7    6'
            write(unit_id, '(A)') '   -0.25818484021E-02   10    1    7    7'
            write(unit_id, '(A)') '   -0.42692007067E-03   10    1    8    3'
            write(unit_id, '(A)') '   -0.32449377861E-02   10    1    8    8'
            write(unit_id, '(A)') '   -0.42692007067E-03   10    1    9    4'
            write(unit_id, '(A)') '   -0.32449377861E-02   10    1    9    9'
            write(unit_id, '(A)') '    0.42347146970E-02   10    1   10    1'
            write(unit_id, '(A)') '    0.36032118080E-01   10    2    1    1'
            write(unit_id, '(A)') '   -0.51169265382E-03   10    2    2    1'
            write(unit_id, '(A)') '    0.25032329445E-01   10    2    2    2'
            write(unit_id, '(A)') '    0.23827492755E-01   10    2    3    3'
            write(unit_id, '(A)') '    0.23827492755E-01   10    2    4    4'
            write(unit_id, '(A)') '   -0.23063334570E-03   10    2    5    1'
            write(unit_id, '(A)') '   -0.51591025887E-03   10    2    5    2'
            write(unit_id, '(A)') '    0.23554346201E-01   10    2    5    5'
            write(unit_id, '(A)') '   -0.26485037615E-03   10    2    6    1'
            write(unit_id, '(A)') '   -0.31295027157E-02   10    2    6    2'
            write(unit_id, '(A)') '    0.34494177544E-02   10    2    6    5'
            write(unit_id, '(A)') '   -0.25571416458E-01   10    2    6    6'
            write(unit_id, '(A)') '    0.94215484865E-03   10    2    7    1'
            write(unit_id, '(A)') '   -0.20173358648E-02   10    2    7    2'
            write(unit_id, '(A)') '    0.49557634117E-02   10    2    7    5'
            write(unit_id, '(A)') '    0.21686924900E-02   10    2    7    6'
            write(unit_id, '(A)') '   -0.69314080490E-02   10    2    7    7'
            write(unit_id, '(A)') '    0.17943412170E-02   10    2    8    3'
            write(unit_id, '(A)') '   -0.83196872746E-02   10    2    8    8'
            write(unit_id, '(A)') '    0.17943412170E-02   10    2    9    4'
            write(unit_id, '(A)') '   -0.83196872746E-02   10    2    9    9'
            write(unit_id, '(A)') '    0.35601423246E-02   10    2   10    1'
            write(unit_id, '(A)') '    0.12047629067E-01   10    2   10    2'
            write(unit_id, '(A)') '    0.14673774416E-02   10    3    3    1'
            write(unit_id, '(A)') '    0.12071158096E-01   10    3    3    2'
            write(unit_id, '(A)') '   -0.93801661712E-03   10    3    5    3'
            write(unit_id, '(A)') '   -0.72390456503E-03   10    3    6    3'
            write(unit_id, '(A)') '    0.44464563393E-03   10    3    7    3'
            write(unit_id, '(A)') '   -0.66296134760E-03   10    3    8    1'
            write(unit_id, '(A)') '    0.62910978514E-03   10    3    8    2'
            write(unit_id, '(A)') '   -0.37472223974E-02   10    3    8    5'
            write(unit_id, '(A)') '    0.67569386042E-03   10    3    8    6'
            write(unit_id, '(A)') '    0.86875942289E-02   10    3    8    7'
            write(unit_id, '(A)') '    0.42563060127E-02   10    3   10    3'
            write(unit_id, '(A)') '    0.14673774416E-02   10    4    4    1'
            write(unit_id, '(A)') '    0.12071158096E-01   10    4    4    2'
            write(unit_id, '(A)') '   -0.93801661712E-03   10    4    5    4'
            write(unit_id, '(A)') '   -0.72390456503E-03   10    4    6    4'
            write(unit_id, '(A)') '    0.44464563393E-03   10    4    7    4'
            write(unit_id, '(A)') '   -0.66296134760E-03   10    4    9    1'
            write(unit_id, '(A)') '    0.62910978514E-03   10    4    9    2'
            write(unit_id, '(A)') '   -0.37472223974E-02   10    4    9    5'
            write(unit_id, '(A)') '    0.67569386042E-03   10    4    9    6'
            write(unit_id, '(A)') '    0.86875942289E-02   10    4    9    7'
            write(unit_id, '(A)') '    0.42563060127E-02   10    4   10    4'
            write(unit_id, '(A)') '   -0.26848199120E-01   10    5    1    1'
            write(unit_id, '(A)') '   -0.70639097558E-03   10    5    2    1'
            write(unit_id, '(A)') '   -0.19913605550E-01   10    5    2    2'
            write(unit_id, '(A)') '   -0.15192773400E-01   10    5    3    3'
            write(unit_id, '(A)') '   -0.15192773400E-01   10    5    4    4'
            write(unit_id, '(A)') '    0.19822026018E-02   10    5    5    1'
            write(unit_id, '(A)') '    0.12872022466E-01   10    5    5    2'
            write(unit_id, '(A)') '   -0.12649454250E-01   10    5    5    5'
            write(unit_id, '(A)') '    0.47098413747E-03   10    5    6    1'
            write(unit_id, '(A)') '    0.15663120429E-02   10    5    6    2'
            write(unit_id, '(A)') '    0.69561376502E-03   10    5    6    5'
            write(unit_id, '(A)') '   -0.26848199120E-01   10    5    6    6'
            write(unit_id, '(A)') '    0.15663120429E-02   10    5    7    1'
            write(unit_id, '(A)') '    0.41587274600E-02   10    5    7    2'
            write(unit_id, '(A)') '    0.60423078748E-02   10    5    7    5'
            write(unit_id, '(A)') '   -0.70639097559E-03   10    5    7    6'
            write(unit_id, '(A)') '   -0.19913605550E-01   10    5    7    7'
            write(unit_id, '(A)') '   -0.28200553504E-02   10    5    8    3'
            write(unit_id, '(A)') '   -0.15192773400E-01   10    5    8    8'
            write(unit_id, '(A)') '   -0.28200553504E-02   10    5    9    4'
            write(unit_id, '(A)') '   -0.15192773400E-01   10    5    9    9'
            write(unit_id, '(A)') '   -0.69561376502E-03   10    5   10    1'
            write(unit_id, '(A)') '   -0.60423078748E-02   10    5   10    2'
            write(unit_id, '(A)') '    0.12067812377E-01   10    5   10    5'
            write(unit_id, '(A)') '    0.12885230094E-02   10    6    1    1'
            write(unit_id, '(A)') '   -0.15611875539E-03   10    6    2    1'
            write(unit_id, '(A)') '   -0.10855735869E-03   10    6    2    2'
            write(unit_id, '(A)') '    0.25169554469E-03   10    6    3    3'
            write(unit_id, '(A)') '    0.25169554469E-03   10    6    4    4'
            write(unit_id, '(A)') '    0.18322489151E-03   10    6    5    1'
            write(unit_id, '(A)') '    0.49183579199E-03   10    6    5    2'
            write(unit_id, '(A)') '    0.10274932039E-02   10    6    5    5'
            write(unit_id, '(A)') '   -0.63388586124E-03   10    6    6    1'
            write(unit_id, '(A)') '   -0.32575861986E-02   10    6    6    2'
            write(unit_id, '(A)') '    0.15281018528E-02   10    6    6    5'
            write(unit_id, '(A)') '   -0.23510693561E-01   10    6    6    6'
            write(unit_id, '(A)') '   -0.14415017417E-03   10    6    7    1'
            write(unit_id, '(A)') '   -0.38311928397E-03   10    6    7    2'
            write(unit_id, '(A)') '   -0.23063334570E-03   10    6    7    5'
            write(unit_id, '(A)') '    0.40301071636E-02   10    6    7    6'
            write(unit_id, '(A)') '   -0.19240472786E-02   10    6    7    7'
            write(unit_id, '(A)') '    0.51972466625E-03   10    6    8    3'
            write(unit_id, '(A)') '   -0.24144100660E-02   10    6    8    8'
            write(unit_id, '(A)') '    0.51972466625E-03   10    6    9    4'
            write(unit_id, '(A)') '   -0.24144100660E-02   10    6    9    9'
            write(unit_id, '(A)') '    0.13754341276E-03   10    6   10    1'
            write(unit_id, '(A)') '    0.15138436977E-02   10    6   10    2'
            write(unit_id, '(A)') '   -0.19822026018E-02   10    6   10    5'
            write(unit_id, '(A)') '    0.57751228916E-02   10    6   10    6'
            write(unit_id, '(A)') '    0.35609078534E-01   10    7    1    1'
            write(unit_id, '(A)') '   -0.99885138838E-03   10    7    2    1'
            write(unit_id, '(A)') '    0.22151564675E-01   10    7    2    2'
            write(unit_id, '(A)') '    0.20525696004E-01   10    7    3    3'
            write(unit_id, '(A)') '    0.20525696004E-01   10    7    4    4'
            write(unit_id, '(A)') '    0.49183579198E-03   10    7    5    1'
            write(unit_id, '(A)') '   -0.29654105531E-02   10    7    5    2'
            write(unit_id, '(A)') '    0.23840531387E-01   10    7    5    5'
            write(unit_id, '(A)') '   -0.39256197191E-03   10    7    6    1'
            write(unit_id, '(A)') '   -0.18623498631E-02   10    7    6    2'
            write(unit_id, '(A)') '    0.14173834149E-03   10    7    6    5'
            write(unit_id, '(A)') '    0.24847734942E-01   10    7    6    6'
            write(unit_id, '(A)') '   -0.95133002289E-03   10    7    7    1'
            write(unit_id, '(A)') '   -0.54438691459E-02   10    7    7    2'
            write(unit_id, '(A)') '   -0.51591025887E-03   10    7    7    5'
            write(unit_id, '(A)') '    0.82183473553E-03   10    7    7    6'
            write(unit_id, '(A)') '    0.80688684698E-02   10    7    7    7'
            write(unit_id, '(A)') '    0.37198110259E-02   10    7    8    3'
            write(unit_id, '(A)') '    0.75495433303E-03   10    7    8    8'
            write(unit_id, '(A)') '    0.37198110259E-02   10    7    9    4'
            write(unit_id, '(A)') '    0.75495433303E-03   10    7    9    9'
            write(unit_id, '(A)') '    0.22918823421E-02   10    7   10    1'
            write(unit_id, '(A)') '    0.12133194409E-01   10    7   10    2'
            write(unit_id, '(A)') '   -0.12872022466E-01   10    7   10    5'
            write(unit_id, '(A)') '    0.43308083154E-02   10    7   10    6'
            write(unit_id, '(A)') '    0.29337139736E-01   10    7   10    7'
            write(unit_id, '(A)') '    0.27560857734E-03   10    8    3    1'
            write(unit_id, '(A)') '    0.42872299709E-02   10    8    3    2'
            write(unit_id, '(A)') '   -0.24245682958E-02   10    8    5    3'
            write(unit_id, '(A)') '   -0.28466432004E-03   10    8    6    3'
            write(unit_id, '(A)') '    0.18467178293E-02   10    8    7    3'
            write(unit_id, '(A)') '   -0.79370563794E-03   10    8    8    1'
            write(unit_id, '(A)') '   -0.13620273913E-02   10    8    8    2'
            write(unit_id, '(A)') '   -0.93801661712E-03   10    8    8    5'
            write(unit_id, '(A)') '    0.26878675252E-03   10    8    8    6'
            write(unit_id, '(A)') '   -0.19061411205E-02   10    8    8    7'
            write(unit_id, '(A)') '    0.63085162628E-03   10    8   10    3'
            write(unit_id, '(A)') '    0.11681694683E-01   10    8   10    8'
            write(unit_id, '(A)') '    0.27560857734E-03   10    9    4    1'
            write(unit_id, '(A)') '    0.42872299709E-02   10    9    4    2'
            write(unit_id, '(A)') '   -0.24245682958E-02   10    9    5    4'
            write(unit_id, '(A)') '   -0.28466432004E-03   10    9    6    4'
            write(unit_id, '(A)') '    0.18467178293E-02   10    9    7    4'
            write(unit_id, '(A)') '   -0.79370563794E-03   10    9    9    1'
            write(unit_id, '(A)') '   -0.13620273913E-02   10    9    9    2'
            write(unit_id, '(A)') '   -0.93801661712E-03   10    9    9    5'
            write(unit_id, '(A)') '    0.26878675252E-03   10    9    9    6'
            write(unit_id, '(A)') '   -0.19061411205E-02   10    9    9    7'
            write(unit_id, '(A)') '    0.63085162628E-03   10    9   10    4'
            write(unit_id, '(A)') '    0.11681694683E-01   10    9   10    9'
            write(unit_id, '(A)') '    0.16750323527       10   10    1    1'
            write(unit_id, '(A)') '    0.17217164860E-02   10   10    2    1'
            write(unit_id, '(A)') '    0.15193566117       10   10    2    2'
            write(unit_id, '(A)') '    0.13266789951       10   10    3    3'
            write(unit_id, '(A)') '    0.13266789951       10   10    4    4'
            write(unit_id, '(A)') '   -0.10274932039E-02   10   10    5    1'
            write(unit_id, '(A)') '   -0.23840531387E-01   10   10    5    2'
            write(unit_id, '(A)') '    0.13499295253       10   10    5    5'
            write(unit_id, '(A)') '    0.21635429134E-03   10   10    6    1'
            write(unit_id, '(A)') '    0.62023397312E-02   10   10    6    2'
            write(unit_id, '(A)') '   -0.11468694313E-01   10   10    6    5'
            write(unit_id, '(A)') '    0.30377515269       10   10    6    6'
            write(unit_id, '(A)') '   -0.24134011285E-02   10   10    7    1'
            write(unit_id, '(A)') '   -0.40602638310E-02   10   10    7    2'
            write(unit_id, '(A)') '   -0.23554346201E-01   10   10    7    5'
            write(unit_id, '(A)') '   -0.51275560816E-02   10   10    7    6'
            write(unit_id, '(A)') '    0.22083460664       10   10    7    7'
            write(unit_id, '(A)') '   -0.29729353954E-02   10   10    8    3'
            write(unit_id, '(A)') '    0.20326768821       10   10    8    8'
            write(unit_id, '(A)') '   -0.29729353954E-02   10   10    9    4'
            write(unit_id, '(A)') '    0.20326768821       10   10    9    9'
            write(unit_id, '(A)') '   -0.42153705731E-02   10   10   10    1'
            write(unit_id, '(A)') '   -0.16215012380E-01   10   10   10    2'
            write(unit_id, '(A)') '   -0.12649454250E-01   10   10   10    5'
            write(unit_id, '(A)') '   -0.28444005389E-02   10   10   10    6'
            write(unit_id, '(A)') '   -0.10673539613E-01   10   10   10    7'
            write(unit_id, '(A)') '    0.23222955692       10   10   10   10'
            write(unit_id, '(A)') '    -5.0104776095        1    1    0    0'
            write(unit_id, '(A)') '    0.16825602123        2    1    0    0'
            write(unit_id, '(A)') '    -1.3444694987        2    2    0    0'
            write(unit_id, '(A)') '    -1.1641813794        3    3    0    0'
            write(unit_id, '(A)') '    -1.1641813794        4    4    0    0'
            write(unit_id, '(A)') '   -0.55111212653E-01    5    1    0    0'
            write(unit_id, '(A)') '    0.10744608355        5    2    0    0'
            write(unit_id, '(A)') '    -1.2064263890        5    5    0    0'
            write(unit_id, '(A)') '   -0.21959082572E-01    6    1    0    0'
            write(unit_id, '(A)') '   -0.17000642444        6    2    0    0'
            write(unit_id, '(A)') '    0.24383954093        6    5    0    0'
            write(unit_id, '(A)') '    -5.0104776095        6    6    0    0'
            write(unit_id, '(A)') '   -0.17000642444        7    1    0    0'
            write(unit_id, '(A)') '    0.83304412852E-02    7    2    0    0'
            write(unit_id, '(A)') '    0.11133814911        7    5    0    0'
            write(unit_id, '(A)') '    0.16825602123        7    6    0    0'
            write(unit_id, '(A)') '    -1.3444694987        7    7    0    0'
            write(unit_id, '(A)') '   -0.50817865497E-02    8    3    0    0'
            write(unit_id, '(A)') '    -1.1641813794        8    8    0    0'
            write(unit_id, '(A)') '   -0.50817865497E-02    9    4    0    0'
            write(unit_id, '(A)') '    -1.1641813794        9    9    0    0'
            write(unit_id, '(A)') '   -0.24383954093       10    1    0    0'
            write(unit_id, '(A)') '   -0.11133814911       10    2    0    0'
            write(unit_id, '(A)') '    0.16764058901       10    5    0    0'
            write(unit_id, '(A)') '    0.55111212653E-01   10    6    0    0'
            write(unit_id, '(A)') '   -0.10744608355       10    7    0    0'
            write(unit_id, '(A)') '    -1.2064263890       10   10    0    0'
            write(unit_id, '(A)') '    -2.4825000000        1    0    0    0'
            write(unit_id, '(A)') '    -2.4750000000        2    0    0    0'
            write(unit_id, '(A)') '   -0.30682000000        3    0    0    0'
            write(unit_id, '(A)') '   -0.84108000000E-01    4    0    0    0'
            write(unit_id, '(A)') '     3.1292000000        5    0    0    0'
            write(unit_id, '(A)') '     3.1292000000        6    0    0    0'
            write(unit_id, '(A)') '     3.2137000000        7    0    0    0'
            write(unit_id, '(A)') '     3.2397000000        8    0    0    0'
            write(unit_id, '(A)') '     3.2397000000        9    0    0    0'
            write(unit_id, '(A)') '     3.4306000000       10    0    0    0'
            write(unit_id, '(A)') '     1.7817547824        0    0    0    0'
        end subroutine write_large_Li2_FCIDUMP


        subroutine write_small_Li2_FCIDUMP(unit_id)
            integer, intent(in) :: unit_id
            write(unit_id, '(A)') '  &FCI NORB=  6,NELEC=  6,MS2=  0,'
            write(unit_id, '(A)') '  ORBSYM= 1, 1, 1, 1, 1, 1,'
            write(unit_id, '(A)') '  ISYM=0'
            write(unit_id, '(A)') ' &END'
            write(unit_id, '(A)') '     1.6327627870        1    1    1    1'
            write(unit_id, '(A)') '   -0.13816660653        2    1    1    1'
            write(unit_id, '(A)') '    0.17265361682E-01    2    1    2    1'
            write(unit_id, '(A)') '    0.34397810102        2    2    1    1'
            write(unit_id, '(A)') '   -0.53019371999E-02    2    2    2    1'
            write(unit_id, '(A)') '    0.25068860113        2    2    2    2'
            write(unit_id, '(A)') '    0.42029080786E-02    3    1    3    1'
            write(unit_id, '(A)') '    0.63322371809E-02    3    2    3    1'
            write(unit_id, '(A)') '    0.50916769052E-01    3    2    3    2'
            write(unit_id, '(A)') '    0.29004307384        3    3    1    1'
            write(unit_id, '(A)') '   -0.35754113716E-02    3    3    2    1'
            write(unit_id, '(A)') '    0.22795394153        3    3    2    2'
            write(unit_id, '(A)') '    0.22612217356        3    3    3    3'
            write(unit_id, '(A)') '    0.75043896954E-02    4    1    1    1'
            write(unit_id, '(A)') '   -0.10585776208E-02    4    1    2    1'
            write(unit_id, '(A)') '    0.19681432522E-03    4    1    2    2'
            write(unit_id, '(A)') '    0.12946913459E-03    4    1    3    3'
            write(unit_id, '(A)') '    0.19848698782E-03    4    1    4    1'
            write(unit_id, '(A)') '   -0.97992463811E-02    4    2    1    1'
            write(unit_id, '(A)') '    0.41363792273E-03    4    2    2    1'
            write(unit_id, '(A)') '   -0.41124456245E-02    4    2    2    2'
            write(unit_id, '(A)') '   -0.40070999477E-02    4    2    3    3'
            write(unit_id, '(A)')  '0.58481893008E-03    4    2    4    1'
            write(unit_id, '(A)') '    0.44589483716E-02    4    2    4    2'
            write(unit_id, '(A)') '   -0.32605732421E-03    4    3    3    1'
            write(unit_id, '(A)') '   -0.23156361648E-02    4    3    3    2'
            write(unit_id, '(A)') '    0.32058144665E-03    4    3    4    3'
            write(unit_id, '(A)') '    0.20051849601        4    4    1    1'
            write(unit_id, '(A)') '    0.41735324022E-02    4    4    2    1'
            write(unit_id, '(A)') '    0.18878583186        4    4    2    2'
            write(unit_id, '(A)') '    0.15749127301        4    4    3    3'
            write(unit_id, '(A)') '    0.75043896954E-02    4    4    4    1'
            write(unit_id, '(A)') '    0.66793041000E-01    4    4    4    2'
            write(unit_id, '(A)') '     1.6327627870        4    4    4    4'
            write(unit_id, '(A)') '    0.66793041000E-01    5    1    1    1'
            write(unit_id, '(A)') '   -0.71513328264E-02    5    1    2    1'
            write(unit_id, '(A)') '    0.65121751160E-02    5    1    2    2'
            write(unit_id, '(A)') '    0.53241021228E-02    5    1    3    3'
            write(unit_id, '(A)') '    0.58481893008E-03    5    1    4    1'
            write(unit_id, '(A)') '   -0.53514496182E-03    5    1    4    2'
            write(unit_id, '(A)') '   -0.97992463811E-02    5    1    4    4'
            write(unit_id, '(A)') '    0.44589483716E-02    5    1    5    1'
            write(unit_id, '(A)') '   -0.24959346753E-01    5    2    1    1'
            write(unit_id, '(A)') '    0.79717594725E-03    5    2    2    1'
            write(unit_id, '(A)') '   -0.10346561190E-01    5    2    2    2'
            write(unit_id, '(A)') '   -0.78447694255E-02    5    2    3    3'
            write(unit_id, '(A)') '    0.21941615579E-04    5    2    4    1'
            write(unit_id, '(A)') '    0.28662672467E-04    5    2    4    2'
            write(unit_id, '(A)') '   -0.24959346753E-01    5    2    4    4'
            write(unit_id, '(A)') '    0.28662672467E-04    5    2    5    1'
            write(unit_id, '(A)') '    0.32376065304E-02    5    2    5    2'
            write(unit_id, '(A)') '   -0.21173245038E-03    5    3    3    1'
            write(unit_id, '(A)') '    0.38589709967E-03    5    3    3    2'
            write(unit_id, '(A)') '    0.16799583491E-03    5    3    4    3'
            write(unit_id, '(A)') '    0.12651804505E-02    5    3    5    3'
            write(unit_id, '(A)') '    0.41735324022E-02    5    4    1    1'
            write(unit_id, '(A)') '   -0.36901184632E-03    5    4    2    1'
            write(unit_id, '(A)') '    0.91458714190E-03    5    4    2    2'
            write(unit_id, '(A)') '    0.15869063848E-02    5    4    3    3'
            write(unit_id, '(A)') '   -0.10585776208E-02    5    4    4    1'
            write(unit_id, '(A)') '   -0.71513328264E-02    5    4    4    2'
            write(unit_id, '(A)') '   -0.13816660653        5    4    4    4'
            write(unit_id, '(A)') '    0.41363792273E-03    5    4    5    1'
            write(unit_id, '(A)') '    0.79717594725E-03    5    4    5    2'
            write(unit_id, '(A)') '    0.17265361682E-01    5    4    5    4'
            write(unit_id, '(A)') '    0.18878583186        5    5    1    1'
            write(unit_id, '(A)') '    0.91458714190E-03    5    5    2    1'
            write(unit_id, '(A)') '    0.16743661593        5    5    2    2'
            write(unit_id, '(A)') '    0.14955029106        5    5    3    3'
            write(unit_id, '(A)') '    0.19681432522E-03    5    5    4    1'
            write(unit_id, '(A)') '    0.65121751160E-02    5    5    4    2'
            write(unit_id, '(A)') '    0.34397810102        5    5    4    4'
            write(unit_id, '(A)') '   -0.41124456245E-02    5    5    5    1'
            write(unit_id, '(A)') '   -0.10346561190E-01    5    5    5    2'
            write(unit_id, '(A)') '   -0.53019371999E-02    5    5    5    4'
            write(unit_id, '(A)') '    0.25068860113        5    5    5    5'
            write(unit_id, '(A)') '   -0.82333190944E-03    6    1    3    1'
            write(unit_id, '(A)') '   -0.18412795676E-02    6    1    3    2'
            write(unit_id, '(A)') '    0.14623348073E-03    6    1    4    3'
            write(unit_id, '(A)') '   -0.88597661874E-04    6    1    5    3'
            write(unit_id, '(A)') '    0.32058144665E-03    6    1    6    1'
            write(unit_id, '(A)') '   -0.35618532176E-03    6    2    3    1'
            write(unit_id, '(A)') '   -0.75979938319E-03    6    2    3    2'
            write(unit_id, '(A)') '   -0.88597661874E-04    6    2    4    3'
            write(unit_id, '(A)') '    0.74379496651E-04    6    2    5    3'
            write(unit_id, '(A)') '    0.16799583491E-03    6    2    6    1'
            write(unit_id, '(A)') '    0.12651804505E-02    6    2    6    2'
            write(unit_id, '(A)') '   -0.79810529730E-02    6    3    1    1'
            write(unit_id, '(A)') '    0.40240805622E-03    6    3    2    1'
            write(unit_id, '(A)') '   -0.17890067998E-02    6    3    2    2'
            write(unit_id, '(A)') '   -0.10096059220E-02    6    3    3    3'
            write(unit_id, '(A)') '   -0.14350926511E-03    6    3    4    1'
            write(unit_id, '(A)') '   -0.77434722556E-03    6    3    4    2'
            write(unit_id, '(A)') '   -0.79810529730E-02    6    3    4    4'
            write(unit_id, '(A)') '   -0.77434722556E-03    6    3    5    1'
            write(unit_id, '(A)') '   -0.32701342776E-03    6    3    5    2'
            write(unit_id, '(A)') '    0.40240805622E-03    6    3    5    4'
            write(unit_id, '(A)') '   -0.17890067998E-02    6    3    5    5'
            write(unit_id, '(A)') '    0.26045402532E-02    6    3    6    3'
            write(unit_id, '(A)') '    0.52578962662E-04    6    4    3    1'
            write(unit_id, '(A)') '    0.61941348667E-03    6    4    3    2'
            write(unit_id, '(A)') '   -0.82333190944E-03    6    4    4    3'
            write(unit_id, '(A)') '   -0.35618532176E-03    6    4    5    3'
            write(unit_id, '(A)') '   -0.32605732421E-03    6    4    6    1'
            write(unit_id, '(A)') '   -0.21173245038E-03    6    4    6    2'
            write(unit_id, '(A)') '    0.42029080786E-02    6    4    6    4'
            write(unit_id, '(A)') '    0.61941348667E-03    6    5    3    1'
            write(unit_id, '(A)') '    0.12480416371E-01    6    5    3    2'
            write(unit_id, '(A)') '   -0.18412795676E-02    6    5    4    3'
            write(unit_id, '(A)') '   -0.75979938319E-03    6    5    5    3'
            write(unit_id, '(A)') '   -0.23156361648E-02    6    5    6    1'
            write(unit_id, '(A)') '    0.38589709967E-03    6    5    6    2'
            write(unit_id, '(A)') '    0.63322371809E-02    6    5    6    4'
            write(unit_id, '(A)') '    0.50916769052E-01    6    5    6    5'
            write(unit_id, '(A)') '    0.15749127301        6    6    1    1'
            write(unit_id, '(A)') '    0.15869063848E-02    6    6    2    1'
            write(unit_id, '(A)') '    0.14955029106        6    6    2    2'
            write(unit_id, '(A)') '    0.13770968110        6    6    3    3'
            write(unit_id, '(A)') '    0.12946913459E-03    6    6    4    1'
            write(unit_id, '(A)') '    0.53241021228E-02    6    6    4    2'
            write(unit_id, '(A)') '    0.29004307384        6    6    4    4'
            write(unit_id, '(A)') '   -0.40070999477E-02    6    6    5    1'
            write(unit_id, '(A)') '   -0.78447694255E-02    6    6    5    2'
            write(unit_id, '(A)') '   -0.35754113716E-02    6    6    5    4'
            write(unit_id, '(A)') '    0.22795394153        6    6    5    5'
            write(unit_id, '(A)') '   -0.10096059220E-02    6    6    6    3'
            write(unit_id, '(A)') '    0.22612217356        6    6    6    6'
            write(unit_id, '(A)') '    -5.0104776095        1    1    0    0'
            write(unit_id, '(A)') '    0.16825602123        2    1    0    0'
            write(unit_id, '(A)') '    -1.3444694987        2    2    0    0'
            write(unit_id, '(A)') '    -1.1641813794        3    3    0    0'
            write(unit_id, '(A)') '   -0.21959082572E-01    4    1    0    0'
            write(unit_id, '(A)') '   -0.17000642444        4    2    0    0'
            write(unit_id, '(A)') '    -5.0104776095        4    4    0    0'
            write(unit_id, '(A)') '   -0.17000642444        5    1    0    0'
            write(unit_id, '(A)') '    0.83304412852E-02    5    2    0    0'
            write(unit_id, '(A)') '    0.16825602123        5    4    0    0'
            write(unit_id, '(A)') '    -1.3444694987        5    5    0    0'
            write(unit_id, '(A)') '   -0.50817865497E-02    6    3    0    0'
            write(unit_id, '(A)') '    -1.1641813794        6    6    0    0'
            write(unit_id, '(A)') '    -2.4825000000        1    0    0    0'
            write(unit_id, '(A)') '    -2.4750000000        2    0    0    0'
            write(unit_id, '(A)') '   -0.30682000000        3    0    0    0'
            write(unit_id, '(A)') '   -0.84108000000E-01    4    0    0    0'
            write(unit_id, '(A)') '     3.1292000000        5    0    0    0'
            write(unit_id, '(A)') '     3.1292000000        6    0    0    0'
            write(unit_id, '(A)') '     1.7817547824        0    0    0    0'
        end subroutine write_small_Li2_FCIDUMP
    end subroutine test_pgen


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
                                n_total=1)))

        call assert_equals( &
            [1, 2], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpinOrbIdx_t([1, 5]), &
                                n_total=2), &
            2)


        call assert_equals( &
            [1, 1], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpinOrbIdx_t([1, 2]), &
                                n_total=2), &
            2)

        call assert_equals( &
            [2, 2], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpinOrbIdx_t([5, 6]), &
                                n_total=2), &
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
                                n_total=1)))

        call assert_equals( &
            [1, 2], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpatOrbIdx_t([1, 3]), &
                                n_total=2), &
            2)


        call assert_equals( &
            [1, 1], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpatOrbIdx_t([1, 1]), &
                                n_total=2), &
            2)

        call assert_equals( &
            [2, 2], &
            get_possible_spaces(GAS_spec, splitted_det_I, &
                                add_holes=SpatOrbIdx_t([3, 3]), &
                                n_total=2), &
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
              excess=alpha)
        call assert_true(size(expected) == size(calculated))
        call assert_true(all(expected == calculated))

        expected = SpinOrbIdx_t([8])
        calculated = get_possible_holes( &
              GAS_spec, reference, add_holes=SpinOrbIdx_t([5]), &
              excess=beta)
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
              n_total=2)
        call assert_true(size(expected) == size(calculated))
        call assert_true(all(expected == calculated))

        expected = SpinOrbIdx_t([2, 3, 4])
        calculated = get_possible_holes( &
              GAS_spec, excite(reference, SingleExc_t(5, 11)), &
              add_holes=SpinOrbIdx_t([1]))
        call assert_true(size(expected) == size(calculated))
        call assert_true(all(expected == calculated))

        expected = SpinOrbIdx_t([2, 3, 4])
        calculated = get_possible_holes( &
              GAS_spec, reference, &
              add_particles=SpinOrbIdx_t([11]), &
              add_holes=SpinOrbIdx_t([1, 5]))
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


        associate(calculated => split_per_GAS(GAS_spec, SpinOrbIdx_t([1, 2, 5, 6, 9, 10])), &
                  expected => [SpinOrbIdx_t([1, 2]), SpinOrbIdx_t([5, 6]), SpinOrbIdx_t([9, 10])])
            do i = 1, get_nGAS(GAS_spec)
                call assert_true(all(expected(i) == calculated(i)))
            end do
        end associate

    end subroutine

end module test_gasci_mod

program test_gasci_program

    use mpi
    use fruit
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_gasci_mod, only: test_igas_from_spatorb, test_igas_from_spinorb, &
        test_contains_det_spinorb, test_contains_det_spatorb, &
        test_particles_per_GAS_spatorb, test_particles_per_GAS_spinorb, &
        test_is_valid, test_is_connected, &
        test_get_possible_spaces_spinorb, test_get_possible_spaces_spatorb, &
        test_possible_holes, test_split_per_GAS, test_available, test_pgen


    implicit none
    integer :: failed_count, err

    integer :: n
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
        call run_test_case(test_available, "test_available")
        call run_test_case(test_pgen, "test_pgen")
    end subroutine
end program test_gasci_program
