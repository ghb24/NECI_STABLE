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
        call run_test_case(calc_pgen_back_spawn_hubbard_test, "calc_pgen_back_spawn_hubbard_test")

    end subroutine back_spawn_excit_gen_test_driver

    subroutine calc_pgen_back_spawn_ueg_test
        ! actually this routine does not take too much effort to test 
        ! maybe i can even test the ueg and hubbard excitation generator..
        ! although there is annoying symmetry stuff to be set up.. 
        ! indirectly i can also test calcpgenlattice with this unit test.. 
        ! since it is part of calc_pgen_back_spawn_ueg
        use SystemData, only: nel, tUEG, nOccBeta, nOccAlpha, tNoFailAb, nBasis
        use FciMCData, only: projedet, ilutref
        use dSFMT_interface, only: dSFMT_init
        use constants, only: dp, n_int
        use bit_reps, only: set_flag, get_initiator_flag, test_flag
        use bit_rep_data, only: niftot, noffflag, tUseflags
        use detbitops, only: encodebitdet
        use CalcData, only: occ_virt_level

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
        occ_virt_level = 0

        allocate(projedet(nel,1)); projedet(:,1) = [1,2]
        allocate(ilutref(0:niftot,1))
        call encodebitdet(projedet, ilutref)
 
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

        ! now do a bit more advanced tests.. 
        nel = 4 
        nBasis = 8 
        nOccBeta = 2 
        nOccAlpha = 2

        ! i am not quite sure why this can be called since nI is too short 
        ! now actually..  but not used.. so maybe thats the reason
        call assert_equals(1.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))

        ! change to alpha excitation: 
        ex(1,1) = 2
        call assert_equals(1.0_dp/6.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))
        ! no beta:
        ex(1,:) = 1
        call assert_equals(1.0_dp/6.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))

        ! and now test with back-spawn:
        call set_flag(ilut, get_initiator_flag(run), .false.)

        ! first check with ic = 1
        ic = 1 
        call assert_equals(0.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))

        ! and then with 
        ! electrons outside the occupied manifold 
        ic = 2
        deallocate(projedet); allocate(projedet(nel,1)); projedet(:,1) = [3,4,5,6]
        deallocate(ilutref);  allocate(ilutref(0:niftot,1)); 
        call encodebitdet(projedet, ilutref)

        ! so this is now same as the above beta pgen
        call assert_equals(1.0_dp/6.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))
        ex(1,1) = 2 
        call assert_equals(1.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))
        ex(1,2) = 2
        call assert_equals(1.0_dp/6.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))

        ! and now actually test with recalculated pgen.. 
        ! and now think hard.. is the formula actually 
        ! accurate why is there a factor of 4 for the UEG?
        ! i think it is because of the possibility to pick both electrons 
        ! in either order.. and the second 2 is because p(b|ij) = p(a|ij) 
        ! but this is actually maybe not even true in the UEG and hubbard 
        ! case.. it depends 

        ! i guess i got it right now.. so get it done.. 
        ! i want both electrons in the occupied manifold 
        ! i know these setting make no sense but atleast this way i can 
        ! check it
        ex(2,:) = 5
        ! i need ni and ilut 
        deallocate(nI); allocate(nI(nel)); nI = [1,2,3,4]
        call EncodeBitDet(nI, ilut)

        ex(1,:) = 3
        call assert_equals(4.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))

        ex(1,:) = 4 
        call assert_equals(4.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))

        ex(1,:) = [3,4]
        call assert_equals(2.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))


        ! and now with b outside occupied range
        ex(2,:) = 7

        ex(1,:) = 3
        call assert_equals(2.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))

        ex(1,:) = 4 
        call assert_equals(2.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))

        ex(1,:) = [3,4]
        call assert_equals(1.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))

        ! this mimicks a beta-beta excitation with orb b also in the 
        ! occupied manifold.. 
        ! oh damn.. i just realized.. i do not need the restriction to pick 
        ! a beta orbital first in the UEG case! damn.. 
        ! so i need a extra ueg orbital picker.. 
        ! ok.. so are we finished with the prep? 
        deallocate(projedet)
        deallocate(ilutref) 
        nel = -1 
        nOccBeta = -1
        nOccAlpha = -1
        niftot = -1 
        noffflag = -1
        nBasis = -1
        occ_virt_level = 0
        tUseflags = .false.
        tUEG = .false.

    end subroutine calc_pgen_back_spawn_ueg_test

    subroutine calc_pgen_back_spawn_hubbard_test
        use SystemData, only: nel, tHub, nOccBeta, nOccAlpha, nBasis
        use FciMCData, only: projedet, ilutref
        use dSFMT_interface, only: dSFMT_init
        use constants, only: dp, n_int
        use bit_reps, only: set_flag, get_initiator_flag, test_flag
        use bit_rep_data, only: niftot, noffflag, tUseflags
        use detbitops, only: encodebitdet
        use CalcData, only: occ_virt_level, t_back_spawn

        integer, allocatable :: nI(:) 
        integer(n_int), allocatable :: ilut(:)
        integer :: ex(2,2), ic = 1, run = 1

        print *, "" 
        print *, "testing: calc_pgen_back_spawn_hubbard" 

        nel = 2 
        nOccBeta = 1
        nOccAlpha = 1
        tHub = .true. 
        niftot = 2
        noffflag = 1
        tUseflags = .true. 
        nBasis = 4
        occ_virt_level = 0
        t_back_spawn = .false.

        allocate(projedet(nel,1)); projedet(:,1) = [1,2]
        allocate(ilutref(0:niftot,1))
        call encodebitdet(projedet, ilutref)
 
        print *, ""
        print *, "testing: calc_pgen_back_spawn_hubbard "
        print *, "with necessary global data: " 
        print *, "nel: ", nel
        print *, "nOccBeta: ", nOccBeta
        print *, "nOccAlpha: ", nOccAlpha
        print *, "tHub: ", tHub
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
        call assert_equals(0.0_dp, calc_pgen_back_spawn_hubbard(nI, ilut, ex, ic, run))

        ic = 2 
        call assert_equals(1.0_dp, calc_pgen_back_spawn_hubbard(nI, ilut, ex, ic, run))

        ! now do a bit more advanced tests.. 
        nel = 4 
        nBasis = 8 
        nOccBeta = 2 
        nOccAlpha = 2

        ! i am not quite sure why this can be called since nI is too short 
        ! now actually..  but not used.. so maybe thats the reason
        ! a conversion happens! not good actually.. 
        call assert_equals(1.0_dp/8.0_dp, calc_pgen_back_spawn_hubbard(nI, ilut, ex, ic, run))

        ! only ispn == 2 excitation in the hubbard.. 

        ! and now test with back-spawn:
        call set_flag(ilut, get_initiator_flag(run), .false.)

        ! first check with ic = 1
        ic = 1 
        call assert_equals(0.0_dp, calc_pgen_back_spawn_hubbard(nI, ilut, ex, ic, run))

        ! and then with 
        ! electrons outside the occupied manifold 
        ic = 2
        deallocate(projedet); allocate(projedet(nel,1)); projedet(:,1) = [3,4,5,6]
        deallocate(ilutref);  allocate(ilutref(0:niftot,1)); 
        call encodebitdet(projedet, ilutref)

        call assert_equals(1.0_dp/8.0_dp, calc_pgen_back_spawn_hubbard(nI, ilut, ex, ic, run))

        ! and now actually test with recalculated pgen.. 
        ! and now think hard.. is the formula actually 
        ! accurate why is there a factor of 4 for the UEG?
        ! i think it is because of the possibility to pick both electrons 
        ! in either order.. and the second 2 is because p(b|ij) = p(a|ij) 
        ! but this is actually maybe not even true in the UEG and hubbard 
        ! case.. it depends 

        ! i guess i got it right now.. so get it done.. 
        ! i want both electrons in the occupied manifold 
        ! i know these setting make no sense but atleast this way i can 
        ! check it
        ex(2,:) = 5
        ! i need ni and ilut 
        deallocate(nI); allocate(nI(nel)); nI = [1,2,3,4]
        call EncodeBitDet(nI, ilut)

        ! so test if we have both electrons in the occupied manifold
        ! and beta also in the occupied manifold.. 
        ! so i have 2 * 1/4 * 1/2
        ex(1,:) = [3,4]
        call assert_equals(1.0_dp/4.0_dp, calc_pgen_back_spawn_hubbard(nI, ilut, ex, ic, run))

        ! and now with b outside occupied range
        ex(2,:) = 7

        call assert_equals(1.0_dp/8.0_dp, calc_pgen_back_spawn_hubbard(nI, ilut, ex, ic, run))


        ! todo: test also the functionality with the old back-spawn and 
        ! also with different combinations of occ_virt_level.. 

        ! this mimicks a beta-beta excitation with orb b also in the 
        ! occupied manifold.. 
        ! oh damn.. i just realized.. i do not need the restriction to pick 
        ! a beta orbital first in the UEG case! damn.. 
        ! so i need a extra ueg orbital picker.. 
        ! ok.. so are we finished with the prep? 
        deallocate(projedet)
        deallocate(ilutref) 
        nel = -1 
        nOccBeta = -1
        nOccAlpha = -1
        niftot = -1 
        noffflag = -1
        nBasis = -1
        occ_virt_level = 0
        tHub = .false.

    end subroutine calc_pgen_back_spawn_hubbard_test

end program test_back_spawn_excit_gen

