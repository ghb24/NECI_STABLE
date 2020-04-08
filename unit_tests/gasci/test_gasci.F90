module test_gasci_mod
    use fruit
    use constants, only: dp, n_int
    use bit_rep_data, only: NIfTot
    use SystemData, only: nEl
    use sets_mod, only: disjoint
    use util_mod, only: operator(.div.)
    use sort_mod, only: sort
    use procedure_pointers, only: generate_excitation
    use orb_idx_mod, only: SpinOrbIdx_t, SpatOrbIdx_t, SpinProj_t, &
        size, operator(==), alpha, beta, sum, calc_spin, calc_spin_raw, &
        operator(-), to_ilut
    use excitation_types, only: SingleExc_t, DoubleExc_t, excite
    use util_mod, only: cumsum
    use gasci, only: GASSpec_t, get_iGAS, &
        contains_det, get_nGAS, particles_per_GAS, operator(.contains.), &
        is_valid, is_connected, get_possible_spaces, get_possible_holes, &
        split_per_GAS, generate_nGAS_excitation, &
        get_available_singles, get_available_doubles
    use unit_test_helper_excitgen, only: test_excitation_generator, &
        init_excitgen_test, finalize_excitgen_test
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

    subroutine test_available()
        type(SpinOrbIdx_t), allocatable :: expect_singles(:), expect_doubles(:)
        type(SpinOrbIdx_t) :: det_I
        type(GASSpec_t) :: GAS_spec
        integer :: i

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

        block
            type(SpinOrbIdx_t), allocatable :: singles_exc_list(:), doubles_exc_list(:)

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


        block
            type(SpinOrbIdx_t), allocatable :: singles_exc_list(:), doubles_exc_list(:)

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
        use FciMCData, only: pSingles, pDoubles, pParallel
        type(GASSpec_t) :: GAS_spec
        type(SpinOrbIdx_t) :: det_I
        integer, parameter :: n_spat_orbs = 12, n_iters=10**6


        det_I = SpinOrbIdx_t([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18])
        call init_excitgen_test(size(det_I), n_spat_orbs, &
                                1.0_dp, 1.0_dp, sum(calc_spin(det_I)))
        pParallel = 0.5_dp
        pSingles = 0.1_dp
        pDoubles = 1.0_dp - pSingles


        ! two benzenes stacked: disconnected spaces
        GAS_spec = GASSpec_t(&
            n_orbs=[6, 12], &
            n_min=[6, size(det_I)], &
            n_max=[6, size(det_I)])
        call assert_true(is_valid(GAS_spec))
        call assert_true(GAS_spec .contains. det_I)
        global_GAS_spec = GAS_spec
        call run_excit_gen_tester( &
            generate_nGAS_excitation, 'generate_nGAS_excitation', &
            opt_nI=det_I%idx, opt_n_iters=n_iters, &
            gen_all_excits=gen_all_excits)

!         two benzenes stacked: 1exc in both directions
        GAS_spec = GASSpec_t(&
            n_orbs=[6, 12], &
            n_min=[5, size(det_I)], &
            n_max=[7, size(det_I)])
        call assert_true(is_valid(GAS_spec))
        call assert_true(GAS_spec .contains. det_I)
        global_GAS_spec = GAS_spec
        call run_excit_gen_tester( &
            generate_nGAS_excitation, 'generate_nGAS_excitation', &
            opt_nI=det_I%idx, opt_n_iters=n_iters, &
            gen_all_excits=gen_all_excits)

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

        call mpi_init(err)

        call init_fruit()

        call test_gasci_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_gasci_program', 'failed_tests')

        call mpi_finalize(err)
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
