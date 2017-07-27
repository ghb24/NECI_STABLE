! finally use some unit-testing on my newly implemented features 

program test_back_spawn

    use back_spawn
    use fruit 

    implicit none 

    integer :: failed_count

    call init_fruit() 
    call back_spawn_test_driver() 
    call fruit_summary() 
    call fruit_finalize() 

    call get_failed_count(failed_count) 
    if (failed_count /= 0) stop -1 

contains 

    subroutine back_spawn_test_driver() 
        ! main test-case driver for the back-spawn procedures 

        call run_test_case(test_is_in_ref, "test_is_in_ref")
        call run_test_case(test_is_in_virt_mask, "test_is_in_virt_mask")
        call run_test_case(test_init_back_spawn, "test_init_back_spawn")
        call run_test_case(test_setup_virtual_mask, "test_setup_virtual_mask")
        call run_test_case(test_check_electron_location, "test_check_electron_location")
        call run_test_case(test_pick_virtual_electrons_double, "test_pick_virtual_electrons_double")
        call run_test_case(test_pick_occupied_orbital_single, "test_pick_occupied_orbital_single")
        call run_test_case(test_pick_occupied_orbital_hubbard, "test_pick_occupied_orbital_hubbard")
        call run_test_case(test_pick_occupied_orbital, "test_pick_occupied_orbital") 
        call run_test_case(test_pick_second_occupied_orbital, "test_pick_second_occupied_orbital")
        call run_test_case(test_pick_virtual_electrons_double_hubbard, &
            "test_pick_virtual_electrons_double_hubbard")
        call run_test_case(test_pick_virtual_electron_single, "test_pick_virtual_electron_single")
        call run_test_case(test_get_ispn, "test_get_ispn")

    end subroutine back_spawn_test_driver 

    subroutine test_is_in_ref
        use bit_rep_data, only: niftot 
        use fcimcdata, only: ilutref
        use constants, only: n_int
        use DetBitOps, only: EncodeBitDet
        use SystemData, only: nel

        integer :: run

        niftot = 1
        allocate(ilutref(0:niftot,2)); ilutref = 0_n_int

        print *, "" 
        print *, "testing: is_in_ref()" 
        print *, "with necessary global data: " 
        print *, "ilutref: ", ilutref
        print *, "niftot: ", niftot
        print *, "n_int: ", n_int
        print *, "EncodeBitDet() "

        ! also do some tests with multiple replicas: 
#ifdef __CMPLX
        run = 3
#else
        run = 2
#endif

        call assert_true(.not. is_in_ref(1))
        call assert_true(.not. is_in_ref(1,1))
        call assert_true(.not. is_in_ref(1,run))

        nel = 1
        call EncodeBitDet([1], ilutref(:,1))
        call EncodeBitDet([2], ilutref(:,2))

        call assert_true(is_in_ref(1))
        call assert_true(.not. is_in_ref(2))

        call assert_true(is_in_ref(2,run))
        call assert_true(.not. is_in_ref(2,1))

        ilutref = 0_n_int 

        nel = 2
        call EncodeBitDet([1,2],ilutref)

        call assert_true(is_in_ref(1))
        call assert_true(is_in_ref(2))
        call assert_true(.not. is_in_ref(3))

        deallocate(ilutref)
        nel = -1

    end subroutine test_is_in_ref

    subroutine test_is_in_virt_mask
        use bit_rep_data, only: niftot
        use constants, only: n_int
        use DetBitOps, only: EncodeBitDet
        use SystemData, only: nel

        integer :: run
        niftot = 1
        allocate(mask_virt_ilut(0:niftot,2)); mask_virt_ilut = 0_n_int

        ! also do some tests with multiple replicas: 
#ifdef __CMPLX
        run = 3
#else
        run = 2
#endif

        print *, "" 
        print *, "testing: is_in_virt_mask() "
        print *, "with necessary global data: " 
        print *, "mask_virt_ilut", mask_virt_ilut
        print *, "niftot: ", niftot
        print *, "n_int: ", n_int
        print *, "EncodeBitDet() "

        call assert_true(.not. is_in_virt_mask(1))
        call assert_true(.not. is_in_virt_mask(2,1))
        call assert_true(.not. is_in_virt_mask(3,run))

        nel = 1 
        call EncodeBitDet([1],mask_virt_ilut(:,1))
        call EncodeBitDet([2],mask_virt_ilut(:,2))

        call assert_true(is_in_virt_mask(1))
        call assert_true(is_in_virt_mask(1,1))
        call assert_true(is_in_virt_mask(2,2))

        call assert_true(.not. is_in_virt_mask(2,1))

        nel = -1
        niftot = -1

        deallocate(mask_virt_ilut)

    end subroutine test_is_in_virt_mask

    subroutine test_init_back_spawn 
        ! check if everything is set in the initializer.. 
        ! and only use the global variables inside this routine, so i know 
        ! which variables and functions are necessary on top of each routine!
        use SystemData, only: nBasis, nel
        use constants, only: inum_runs
        use fcimcdata, only: max_calc_ex_level, projedet, ilutref
        use bit_rep_data, only: niftot
        use detbitops, only: EncodeBitDet

        integer :: i, j

        nBasis = 12
        nel = 6 
        niftot = 1
        ! depending on the executable i might have to set inum_runs..
        ! but it seems we do not have the preprocessor flags here.. 
        ! thats a shame.. try to get them here..
        ! ok i cant get them running yet.. i could work with if statements.. 
        ! but thats just an intermediate solution!

#ifdef __PROG_NUMRUNS
        inum_runs = 2
#endif
        allocate(projedet(nel,inum_runs)) 
        allocate(ilutref(0:niftot,inum_runs))
        do j = 1, inum_runs
            projedet(:,j) = [(i, i = 1, nel)]
            call EncodeBitDet(projedet, ilutref)
        end do


        tTruncInitiator = .true. 
        tGen_4ind_2 = .true. 
        tGen_4ind_2_symmetric = .false. 

        print *, "" 
        print *, "testing init_back_spawn on System with required global data: " 
        print *, "nBasis: ", nBasis 
        print *, "nel: ", nel 
        print *, "inum_runs: ", inum_runs
        print *, "projedet: ", projedet
        print *, "tTruncInitiator: ", tTruncInitiator
        print *, "tGen_4ind_2: ", tGen_4ind_2
        print *, "tGen_4ind_2_symmetric: ", tGen_4ind_2_symmetric
        print *, "ilutref: ", ilutref
        print *, "niftot: ", niftot
        print *, "EncodeBitDet()"

        call init_back_spawn()

        call assert_true(allocated(mask_virt_ni))
        call assert_equals(size(mask_virt_ni, 1), nBasis - nel)
        call assert_equals([nBasis - nel, inum_runs], shape(mask_virt_ni), 2)
        do j = 1, inum_runs
            call assert_equals(mask_virt_ni(:,j), [7, 8, 9, 10, 11, 12], 2)
        end do
        call assert_equals(max_calc_ex_level, nel)

        ! maybe i have to deallocate all the above initialized global data..
        ! i guess so.. yes i do! 
        nBasis = -1 
        nel = -1
        niftot = -1

        deallocate(projedet)
        deallocate(mask_virt_ni)
        deallocate(ilutref)

#ifdef __PROG_NUMRUNS
        inum_runs = -1
#endif
    end subroutine test_init_back_spawn

    subroutine test_setup_virtual_mask() 
        use fcimcdata, only: projedet, ilutref
        use SystemData, only: nel, nBasis
        use constants, only: inum_runs, dp
        use bit_rep_data, only: niftot
        use detbitops, only: EncodeBitDet

        HElement_t(dp) :: test

        integer :: i, j
        nel = 6
        nBasis = 12 
        niftot = 2

        ! no.. if statements do not work.. since it is on compile-time 
        ! checked if inum_runs is a parameter..
        ! ok works now.. i just had to do it in a fresh build_dir.. duh..
#ifdef __PROG_NUMRUNS
        inum_runs = 2
#endif

        allocate(projedet(nel, inum_runs)) 
        allocate(ilutref(0:niftot,inum_runs))
        do j = 1, inum_runs
            projedet(:,j) = [(i, i = 1, nel)]
            call EncodeBitDet(projedet, ilutref)
        end do

        print *, "" 
        print *, "testing: setup_virtual_mask() on a system with required global data: " 
        print *, "nel: ", nel
        print *, "nBasis: ", nBasis
        print *, "projedet: ", projedet
        print *, "inum_runs: ", inum_runs
        print *, "ilutref: ", ilutref
        print *, "niftot: ", niftot
        print *, "EncodeBitDet()"

        ! for some reason the allocatio/deallocation is not in the setup 
        ! routine.. or atleast was.. change that!
        call setup_virtual_mask()

        call assert_true(allocated(mask_virt_ni))
        do j = 1, inum_runs
            call assert_equals(mask_virt_ni(:,j), [(i,i=7,12)], 2)
        end do

        ! nullify used global data.. so no leak to other test cases happen
        nBasis = -1 
        nel = -1 
        niftot = -1
        deallocate(projedet)
        deallocate(mask_virt_ni)
        deallocate(ilutref)
        deallocate(mask_virt_ilut)

#ifdef __PROG_NUMRUNS
        inum_runs = -1
#endif
    end subroutine test_setup_virtual_mask

    subroutine test_check_electron_location() 
        ! i need a projedet in this case.. 
        use fcimcdata, only: projedet, ilutref
        use bit_rep_data, only: niftot
        use SystemData, only: nel
        use DetBitOps, only: EncodeBitDet

        integer :: i, run2

        niftot = 1
        nel = 4

        allocate(projedet(4,2))
        projedet(:,1) = [(i,i=1,4)]
        projedet(:,2) = [(i,i=5,8)]
        allocate(ilutref(0:niftot,2))

        call EncodeBitDet(projedet(:,1), ilutref(:,1))
        call EncodeBitDet(projedet(:,2), ilutref(:,2))
        
        ! this test does not yet work with the complex code..
        ! maybe i misuse the run variable.. this is good to know, 
        ! since this may cause other parts of the code to break.. 
        print *, "" 
        print *, "testing: check_electron_location() on a system with required global data:"
        print *, "projedet(:,1): ", projedet(:,1)
        print *, "projedet(:,2): ", projedet(:,2)
        print *, "ilutref: ", ilutref
        print *, "niftot: ", niftot
        print *, "nel: ", nel 
        print *, "EncodeBitDet() "


        ! i have to take into account if it is complex.. since 
        ! in the excitation generator i loop over the lenof_sign which 
        ! usually is 2 * inum_runs .. and 2 is mapped to 1 in the 
        ! check_electron_location function.. 
        call assert_equals(check_electron_location([1,2],1,1), 2)
        call assert_equals(check_electron_location([1,2],2,1), 2)
        call assert_equals(check_electron_location([1,5],2,1), 1)

#ifdef __CMPLX 
        run2 = 3
#else 
        run2 = 2
#endif
        call assert_equals(check_electron_location([1,2],1,run2), 0)
        call assert_equals(check_electron_location([1,2],2,run2), 0)
        call assert_equals(check_electron_location([1,5],2,run2), 1)

        deallocate(projedet)
        deallocate(ilutref)
        nel = -1 
        niftot = -1 

    end subroutine test_check_electron_location

    subroutine test_pick_virtual_electrons_double() 
        ! this is a stochastic routine.. but by forcing the choices we can 
        ! test it anyway.. and also maybe only test the pgens not the 
        ! actual picked orbitals. but alot of global information needed.. 
        use constants, only: dp
        use SystemData, only: nel, G1
        use dSFMT_interface, only: dSFMT_init
        use bit_rep_data, only: niftot
        use DetBitOps, only: EncodeBitDet


        real(dp) :: pgen
        integer :: run, elecs(2), src(2), ispn, sum_ml
        integer, allocatable :: nI(:)

        ! i can make really stupid masks here too.. i just want to test 
        ! the pick_virtual_electrons_double functionality
        allocate(mask_virt_ni(1,1))
        mask_virt_ni(1,1) = 1
        nel = 1
        allocate(nI(nel))
        nI(1) = 1
        niftot = 1
        allocate(mask_virt_ilut(0:niftot,1))
        call EncodeBitDet(mask_virt_ni, mask_virt_ilut)

        ! for this code i might not need G1 and dSFMT_init..

        print *, "" 
        print *, "testing: pick_virtual_electrons_double() with necessary global data: "
        print *, "nel: ", nel
        print *, "mask_virt_ni: ", mask_virt_ni

        call pick_virtual_electrons_double(nI, 1, elecs, src, ispn, sum_ml, pgen)

        call assert_equals(elecs, [0,0], 2) 
        call assert_equals(src, [0,0], 2) 
        call assert_equals(ispn, -1)
        call assert_equals(sum_ml, -1)
        call assert_equals(pgen, 0.0_dp)

        ! now a case where it does not exit prematurely
        nel = 2 
        deallocate(nI)
        allocate(ni(nel))
        ! although this case is stupid.. and actually should not be called 
        ! with.. 
        nI = 1

        call dSFMT_init(123)
        allocate(G1(1)) 
        G1(1)%ml = 0

        print *, "" 
        print *, "now with additional global data: "
        print *, "dSFMT_init() "
        print *, "G1: ", G1 

        call pick_virtual_electrons_double(nI, 1, elecs, src, ispn, sum_ml, pgen)

        call assert_equals(elecs, [2,1], 2) 
        call assert_equals(src, [1,1], 2)
        call assert_equals(ispn, 1) 
        call assert_equals(sum_ml, 0)
        call assert_equals(pgen, 1.0_dp)

        ! todo more test cases in the future 
        deallocate(mask_virt_ni)
        deallocate(nI)
        nullify(G1)
        deallocate(mask_virt_ilut)
        niftot = 1
        nel = -1

    end subroutine test_pick_virtual_electrons_double

    ! now the tricky things start.. with a lot ot dependencies on global 
    ! variables.. 
    subroutine test_pick_occupied_orbital_single 
        use SymExcitDataMod, only: OrbClassCount, SymLabelCounts2, SymLabelList2
        use SymExcitDataMod, only: SpinOrbSymLabel
        use fcimcdata, only: projedet, ilutref
        use constants, only: dp, n_int
        use dSFMT_interface, only: dSFMT_init
        use bit_rep_data, only: niftot
        use SystemData, only: nel
        use DetBitOps, only: EncodeBitDet

        integer :: nI(1) = 0, src = 1, cc_index = 1, run = 1
        real(dp) :: pgen 
        integer :: orb
        integer(n_int), allocatable :: ilut(:) 

        niftot = 1
        ! allocate and initialite the necessary data
        allocate(OrbClassCount(1));     OrbClassCount(1) = 1
        allocate(SymLabelCounts2(1,1)); SymLabelCounts2(1,1) = 1
        allocate(SymLabelList2(1));     SymLabelList2(1) = 1
        allocate(SpinOrbSymLabel(1));   SpinOrbSymLabel(1) = 1

        allocate(projedet(1,1));        projedet(1,1) = 2
        allocate(ilutref(0:niftot,1));  ilutref = 0_n_int

        call EncodeBitDet(projedet, ilutref)

        nel = 1

        print *, "" 
        print *, "testing: test_pick_occupied_orbital_single"
        print *, "with necessary global data: " 
        print *, "OrbClassCount: ", OrbClassCount
        print *, "SymLabelCounts2: ", SymLabelCounts2
        print *, "SymLabelList2: ", SymLabelList2
        print *, "SpinOrbSymLabel: ", SpinOrbSymLabel
        print *, "projedet: ", projedet
        print *, "nel: ", nel
        ! first make a test where it definetly fails 
        allocate(ilut(0:niftot));   ilut = 0_n_int

        call pick_occupied_orbital_single(nI, ilut, src, cc_index, run, pgen, orb)

        call assert_equals(pgen, 0.0_dp)
        call assert_equals(orb, 0)

        ! and then do a test where we need the random numbers.. 
        print *, ""
        print *, "now with additional global data: "
        print *, "dSFMT_init() "
        call dSFMT_init(123)

        projedet(1,1) = 1
        call EncodeBitDet(projedet, ilutref)

        ! but to a test where there is a fixed result
        call pick_occupied_orbital_single(nI, ilut, src, cc_index, run, pgen, orb)

        call assert_equals(pgen, 1.0_dp) 
        call assert_equals(orb, 1)

        deallocate(projedet)
        deallocate(SpinOrbSymLabel)
        deallocate(OrbClassCount)
        deallocate(SymLabelCounts2)
        deallocate(SymLabelList2)
        deallocate(ilutref)
        nel = -1 
        niftot = -1


    end subroutine test_pick_occupied_orbital_single

    subroutine test_pick_occupied_orbital_hubbard
        use SystemData, only: nel 
        use fcimcdata, only: projedet, ilutref
        use constants, only: dp, n_int
        use dSFMT_interface, only: dSFMT_init
        use bit_rep_data, only: niftot
        use DetBitOps, only: EncodeBitDet

        integer :: nI(1) = 1, run = 1
        integer(n_int), allocatable :: ilutI(:)
        real(dp) :: pgen 
        integer :: orb

        nel = 1
        niftot = 1
        allocate(projedet(1,1)); projedet(1,1) = 1
        allocate(ilutref(0:niftot,1))

        call EncodeBitDet(projedet, ilutref)

        allocate(ilutI(0:niftot)); ilutI = 0_n_int
        call EncodeBitDet(nI, ilutI)

        print *, "" 
        print *, "testing: pick_occupied_orbital_hubbard() " 
        print *, "with necessary global data: "
        print *, "nel: ", nel 
        print *, "projedet: ", projedet
        call pick_occupied_orbital_hubbard(nI, ilutI, run, pgen, orb) 

        ! those 2 fail not:
        call assert_equals(0.0_dp, pgen)
        call assert_equals(0, orb)

        print *, "now with additional global data: "
        print *, "dSFMT_init()"
        call dSFMT_init(123)

        nI = 2
        call EncodeBitDet(nI, ilutI)

        call pick_occupied_orbital_hubbard(nI, ilutI, run, pgen, orb) 
        call assert_equals(pgen, 1.0_dp) 
        call assert_equals(orb, 1)

        deallocate(projedet)
        deallocate(ilutref)
        niftot = -1
    end subroutine test_pick_occupied_orbital_hubbard

    subroutine test_pick_occupied_orbital
        use SystemData, only: nel 
        use fcimcdata, only: projedet, ilutref
        use constants, only: dp, n_int
        use dSFMT_interface, only: dSFMT_init
        use bit_rep_data, only: niftot
        use DetBitOps, only: EncodeBitDet

        integer :: src(2) = 0, ispn = 0, run = 1
        integer, allocatable :: nI(:)
        real(dp) :: pgen, cum_sum
        integer :: orb
        integer(n_int), allocatable :: ilutI(:)

        nel = 1
        niftot = 1
        allocate(projedet(1,1)); projedet(1,1) = 1
        allocate(ilutref(0:niftot,1)) 
        call EncodeBitDet(projedet, ilutref)

        print *, "" 
        print *, "testing: pick_occupied_orbital() "
        print *, "with necessary global data: " 
        print *, "nel: ", nel 
        print *, "projedet: ", projedet

        ! do one test with no valid orbs first..
        ! although there are various combinations with ispn and such, which 
        ! can lead to non-possible orbs.. which i should probably all test..
        allocate(nI(nel)); nI = 1
        allocate(ilutI(0:niftot))
        call EncodeBitDet(nI, ilutI)

        call pick_occupied_orbital(nI, ilutI, src, ispn, run, pgen, cum_sum, orb)

        call assert_equals(orb, 0)
        call assert_equals(pgen, 0.0_dp)
        call assert_equals(cum_sum, 1.0_dp)

        print *, "now with additional global data: "
        print *, "dSFMT_init()"
        call dSFMT_init(123)
        ! i have to consider ispn here too and that the first picked 
        ! orbital always should be a beta orbital!
        ! apparently even numbers are beta orbitals! since when? 
        ! check that again.. because i am pretty sure odd are beta orbitals..
        ! nah i am just entering the first if-statement of the routine with 
        ! these settings.. 
        projedet(1,1) = 2
        call EncodeBitDet(projedet, ilutref)

        ispn = 1
        call pick_occupied_orbital(nI, ilutI, src, ispn, run, pgen, cum_sum, orb)
        call assert_equals(2, orb)
        call assert_equals(1.0_dp, pgen)
        call assert_equals(1.0_dp, cum_sum)

        ! do also a test with ispn = 2 
        ispn = 2 
        call pick_occupied_orbital(nI, ilutI, src, ispn, run, pgen, cum_sum, orb)
        call assert_equals(0, orb)
        call assert_equals(0.0_dp, pgen)
        call assert_equals(1.0_dp, cum_sum)

        ! and a succesful one.. 
        projedet(1,1) = 1
        call EncodeBitDet(projedet, ilutref)

        nI = 2
        call EncodeBitDet(nI, ilutI)
        call pick_occupied_orbital(nI, ilutI, src, ispn, run, pgen, cum_sum, orb)
        call assert_equals(1, orb)
        call assert_equals(1.0_dp, pgen)
        call assert_equals(1.0_dp, cum_sum)

        ! todo: i should do more tests here
        nel = -1
        niftot = -1
        deallocate(projedet)
        deallocate(ilutref)

    end subroutine test_pick_occupied_orbital

    subroutine test_pick_second_occupied_orbital
        use SymExcitDataMod, only: SymLabelCounts2, OrbClassCount, SymLabelList2
        use fcimcdata, only: projedet, ilutref
        use constants, only: dp, n_int
        use dSFMT_interface, only: dSFMT_init
        use SystemData, only: nel 
        use bit_rep_data, only: niftot
        use DetBitOps, only: EncodeBitDet

        integer, allocatable :: nI(:) 
        integer :: src(2) = 0, cc_b = 1, orb_a = 1, ispn = 0, run = 1
        real(dp) :: cpt, cum_sum 
        integer :: orb
        integer(n_int), allocatable :: ilutI(:)

        nel = 1
        niftot = 1

        ! allocate and initialite the necessary data
        allocate(OrbClassCount(1));     OrbClassCount(1) = 1
        allocate(SymLabelCounts2(1,1)); SymLabelCounts2(1,1) = 1
        allocate(SymLabelList2(1));     SymLabelList2(1) = 1

        allocate(projedet(1,1));        projedet(1,1) = 1

        allocate(ilutref(0:niftot,1))
        call EncodeBitDet(projedet, ilutref)

        print *, "" 
        print *, "testing: pick_second_occupied_orbital() " 
        print *, "with necessary global data: "
        print *, "nel: ", nel 
        print *, "SymLabelCounts2: ", SymLabelCounts2
        print *, "OrbClassCount: ", OrbClassCount
        print *, "SymLabelList2: ", SymLabelList2
        print *, "projedet: ", projedet
        allocate(nI(nel)); nI = 1
        allocate(ilutI(0:niftot)) 
        call EncodeBitDet(nI, ilutI)

        call pick_second_occupied_orbital(nI, ilutI, src, cc_b, orb_a, ispn, run, &
            cpt, cum_sum, orb)

        call assert_equals(0, orb) 
        call assert_equals(0.0_dp, cpt) 
        call assert_equals(1.0_dp, cum_sum)

        ! now do a test with one fixed output: 
        nI = 2
        orb_a = 3

        call EncodeBitDet(nI, ilutI)
        print *, "now with additional global data: "
        print *, "dSFMT_init()" 
        call dSFMT_init(123)

        call pick_second_occupied_orbital(nI, ilutI, src, cc_b, orb_a, ispn, run, &
            cpt, cum_sum, orb)

        call assert_equals(1, orb) 
        call assert_equals(1.0_dp, cpt)
        call assert_equals(1.0_dp, cum_sum)


        ! deallocate:
        deallocate(projedet)
        deallocate(OrbClassCount)
        deallocate(SymLabelCounts2)
        deallocate(SymLabelList2)
        deallocate(ilutref)
        nel = -1 
        niftot = -1

    end subroutine test_pick_second_occupied_orbital
    
    subroutine test_pick_virtual_electrons_double_hubbard
        use SystemData, only: nel
        use constants, only: dp
        use dSFMT_interface, only: dSFMT_init
        use bit_rep_data, only: niftot
        use DetBitOps, only: EncodeBitDet

        integer, allocatable :: nI(:)
        integer :: run = 1, elecs(2), src(2), ispn
        real(dp) :: pgen 

        nel = 1
        niftot = 1
        allocate(mask_virt_ni(1,1)); mask_virt_ni(1,1) = 1
        allocate(mask_virt_ilut(0:niftot,1))
        call EncodeBitDet(mask_virt_ni, mask_virt_ilut)

        print *, ""
        print *, "testing: pick_virtual_electrons_double_hubbard() "
        print *, "with necessary global data: "
        print *, "nel: ", nel 
        print *, "mask_virt_ni: ", mask_virt_ni

        allocate(nI(nel)); nI = 1

        call pick_virtual_electrons_double_hubbard(nI, run, elecs, src, ispn, pgen)
        call assert_equals([0,0], elecs, 2)
        call assert_equals([0,0], src, 2) 
        call assert_equals(0, ispn)
        call assert_equals(0.0_dp, pgen)

        print *, "now with additional global data: "
        print *, "dSFMT_init()" 
        call dSFMT_init(123)

        deallocate(nI)
        nel = 2
        allocate(nI(nel)); nI = [1,2] 

        deallocate(mask_virt_ni)
        allocate(mask_virt_ni(2,1)); mask_virt_ni(:,1) = [1,2]
        call EncodeBitDet(mask_virt_ni, mask_virt_ilut)

        call pick_virtual_electrons_double_hubbard(nI, run, elecs, src, ispn, pgen)

        call assert_equals(2, ispn)
        call assert_equals(1.0_dp, pgen) 
        ! the electrons can be in any order now.. how to test that?
        call assert_true( all(elecs == [1,2]) .or. all(elecs == [2,1]))
        call assert_true( all(src == [1,2]) .or. all(src == [2,1]))

        ! todo: more tests.. 
        nel = -1 
        niftot = -1
        deallocate(mask_virt_ni)
        deallocate(mask_virt_ilut)

    end subroutine test_pick_virtual_electrons_double_hubbard

    subroutine test_pick_virtual_electron_single
        use SystemData, only: nel 
        use constants, only: dp, n_int
        use dSFMT_interface, only: dSFMT_init
        use bit_rep_data, only: niftot
        use DetBitOps, only: EncodeBitDet

        integer, allocatable :: nI(:) 
        integer :: run = 1, elec
        real(dp) :: pgen

        nel = 1
        niftot = 1
        allocate(mask_virt_ni(1,1)); mask_virt_ni(1,1) = 2
        allocate(mask_virt_ilut(0:niftot,1))
        call EncodeBitDet(mask_virt_ni, mask_virt_ilut)

        print *, "" 
        print *, "testing: pick_virtual_electron_single() " 
        print *, "with necessary global data: "
        print *, "nel: ", nel 
        print *, "mask_virt_ni: ", mask_virt_ni

        allocate(nI(nel)); nI = 1

        call pick_virtual_electron_single(nI, run, elec, pgen) 
        
        call assert_equals(0, elec)
        call assert_equals(0.0_dp, pgen) 

        print *, "now with additional global data: "
        print *, "dSFMT_init()" 
        call dSFMT_init(123)

        nI = 2 
        call pick_virtual_electron_single(nI, run, elec, pgen) 

        call assert_equals(1, elec) 
        call assert_equals(1.0_dp, pgen)

        nel = -1 
        niftot = -1
        deallocate(mask_virt_ni)
        deallocate(mask_virt_ilut)


    end subroutine test_pick_virtual_electron_single

    subroutine test_get_ispn

        print *, "" 
        print *, "testing: get_ispn()" 

        call assert_equals(1, get_ispn([1,1]))
        call assert_equals(2, get_ispn([1,2]))
        call assert_equals(3, get_ispn([2,2]))

    end subroutine test_get_ispn
end program test_back_spawn
