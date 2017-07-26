program test_back_spawn_excit_gen

    use back_spawn_excit_gen
    use fruit 

    implicit none 

    integer :: failed_count

    call init_fruit() 
    call back_spawn_excit_gen_test_driver() 
    call fruit_summary() 
    call fruit_finalize() 

    call get_failed_count(failed_count) 
    if (failed_count /= 0) stop -1 

contains 

    subroutine back_spawn_excit_gen_test_driver

        call run_test_case(calc_pgen_back_spawn_ueg_test, "calc_pgen_back_spawn_ueg_test")

    end subroutine back_spawn_excit_gen_test_driver

    subroutine calc_pgen_back_spawn_ueg_test
        ! actually this routine does not take too much effort to test 
        ! maybe i can even test the ueg and hubbard excitation generator..
        ! although there is annoying symmetry stuff to be set up.. 
        ! indirectly i can also test calcpgenlattice with this unit test.. 
        ! since it is part of calc_pgen_back_spawn_ueg
        use SystemData, only: nel, tUEG, nOccBeta, nOccAlpha, tNoFailAb, nBasis
        use FciMCData, only: projedet
        use dSFMT_interface, only: dSFMT_init
        use constants, only: dp, n_int
        use bit_reps, only: set_flag, get_initiator_flag, test_flag
        use bit_rep_data, only: niftot, noffflag, tUseflags

        integer, allocatable :: nI(:) 
        integer(n_int), allocatable :: ilut(:)
        integer :: ex(2,2), ic = 1, run = 1

        nel = 2 
        nOccBeta = 1
        nOccAlpha = 1
        tUEG = .true. 
        tNoFailAb = .false. 
        niftot = 2
        noffflag = 1
        tUseflags = .true. 
        nBasis = 4

        allocate(projedet(1,1)); projedet(1,1) = 1
 
        print *, ""
        print *, "testing: calc_pgen_back_spawn_ueg "
        print *, "with necessary global data: " 
        print *, "nel: ", nel
        print *, "nOccBeta: ", nOccBeta
        print *, "nOccAlpha: ", nOccAlpha
        print *, "tUEG: ", tUEG
        print *, "tNoFailAb: ", tNoFailAb
        print *, "projedet: ", projedet
        print *, "dSFMT_init() "
        print *, "niftot: ", niftot
        print *, "n_int: ", n_int
        print *, "set_flag() "
        print *, "get_initiator_flag()"
        print *, "noffflag: ", noffflag
        print *, "tUseflags: ", tUseflags
        print *, "nBasis: ", nBasis

        call dSFMT_init(123)

        allocate(nI(nel));          nI = [1, 2]
        allocate(ilut(0:niftot));   ilut = 0_n_int

        ex(1,:) = [1,2]
        ex(2,:) = [3,4]

        ! test with initiators first.. 
        call set_flag(ilut, get_initiator_flag(run), .true.)
        ! first get 0 pgen:
        call assert_equals(0.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))

        ic = 2 
        call assert_equals(1.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))
        ! apparently thats it.. 


    end subroutine calc_pgen_back_spawn_ueg_test

end program test_back_spawn_excit_gen

