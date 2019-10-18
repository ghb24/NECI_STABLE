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
        call run_test_case(calc_pgen_back_spawn_ueg_new_test, "calc_pgen_back_spawn_ueg_new_test")

    end subroutine back_spawn_excit_gen_test_driver

    subroutine calc_pgen_back_spawn_ueg_test
        ! actually this routine does not take too much effort to test
        ! maybe i can even test the ueg and hubbard excitation generator..
        ! although there is annoying symmetry stuff to be set up..
        ! indirectly i can also test calcpgenlattice with this unit test..
        ! since it is part of calc_pgen_back_spawn_ueg
        use SystemData, only: nel, tUEG, nOccBeta, nOccAlpha, tNoFailAb, nBasis, &
                              G1, nmaxx, nmaxy, nmaxz, tOrbECutoff
        use FciMCData, only: projedet, ilutref
        use dSFMT_interface, only: dSFMT_init
        use constants, only: dp, n_int
        use bit_reps, only: set_flag, get_initiator_flag, test_flag
        use bit_rep_data, only: niftot, noffflag
        use detbitops, only: encodebitdet
        use CalcData, only: occ_virt_level
        use symexcitdatamod, only: kpointtobasisfn

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
        nBasis = 4
        occ_virt_level = 0

        allocate(G1(nBasis))
        nmaxx = 2
        nmaxy = 2
        nmaxz = 2
        allocate(KPointToBasisFn(-nmaxx:nmaxx, -nmaxy:nmaxy, -nmaxz:nmaxz, 2))
        tOrbECutoff = .false.

        allocate(projedet(nel,1)); projedet(:,1) = [1,2]
        allocate(ilutref(0:niftot,1))
        call encodebitdet(projedet(:,1), ilutref(:,1))

        G1(1)%k = [1,0,0]
        G1(2)%k = [0,1,0]
        G1(3)%k = [1,1,0]
        G1(4)%k = [0,0,0]

        ! i also have to set the ms value in G1
        G1(1)%ms = -1
        G1(2)%ms = 1
        G1(3)%ms = -1
        G1(4)%ms = 1

        KPointToBasisFn(0,1,0,2) = 2 ! i should get kb = [0,1,0]
        KPointToBasisFn(1,0,0,1) = 3
        KPointToBasisFn(0,0,0,2) = 4
        KPointToBasisFn(1,1,0,1) = 1


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
        ! do that later..

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
        call encodebitdet(projedet(:,1), ilutref(:,1))

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

        ! do those tests later.. too much setup hussle..
!         ex(1,:) = 3
!         call assert_equals(4.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))
!
!         ex(1,:) = 4
!         call assert_equals(4.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))
!
!         ex(1,:) = [3,4]
!         call assert_equals(2.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))
!
!
!         ! and now with b outside occupied range
!         ex(2,:) = 7
!
!         ex(1,:) = 3
!         call assert_equals(2.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))
!
!         ex(1,:) = 4
!         call assert_equals(2.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))
!
!         ex(1,:) = [3,4]
!         call assert_equals(1.0_dp/12.0_dp, calc_pgen_back_spawn_ueg(nI, ilut, ex, ic, run))

        ! this mimicks a beta-beta excitation with orb b also in the
        ! occupied manifold..
        ! oh damn.. i just realized.. i do not need the restriction to pick
        ! a beta orbital first in the UEG case! damn..
        ! so i need a extra ueg orbital picker..
        ! ok.. so are we finished with the prep?
        deallocate(projedet)
        deallocate(ilutref)
        deallocate(G1)
        deallocate(kpointtobasisfn)
        nel = -1
        nOccBeta = -1
        nOccAlpha = -1
        niftot = -1
        noffflag = -1
        nBasis = -1
        occ_virt_level = 0
        tUEG = .false.
        nmaxx = -1
        nmaxy = -1
        nmaxz = -1

    end subroutine calc_pgen_back_spawn_ueg_test

    subroutine calc_pgen_back_spawn_hubbard_test
        use SystemData, only: nel, tHub, nOccBeta, nOccAlpha, nBasis
        use FciMCData, only: projedet, ilutref
        use dSFMT_interface, only: dSFMT_init
        use constants, only: dp, n_int
        use bit_reps, only: set_flag, get_initiator_flag, test_flag
        use bit_rep_data, only: niftot, noffflag
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
        nBasis = 4
        occ_virt_level = 0
        t_back_spawn = .false.

        allocate(projedet(nel,1)); projedet(:,1) = [1,2]
        allocate(ilutref(0:niftot,1))
        call encodebitdet(projedet(:,1), ilutref(:,1))

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
        call encodebitdet(projedet(:,1), ilutref(:,1))

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

    subroutine calc_pgen_back_spawn_ueg_new_test
        use SystemData, only: nel, nBasis, G1, nmaxx, nmaxy, nmaxz, ElecPairs, &
                              tOrbECutoff
        use bit_rep_data, only: niftot, noffflag
        use constants, only: dp, n_int
        use detbitops, only: encodebitdet
        use dSFMT_interface, only: dSFMT_init
        use bit_reps, only: set_flag, get_initiator_flag, test_flag
        use symexcitdatamod, only: kpointtobasisfn
        use procedure_pointers, only: get_umat_el
        use CalcData, only: t_back_spawn_flex, occ_virt_level
        use ueg_excit_gens, only: create_ab_list_ueg
        use FciMCData, only: ilutref, projedet

        integer, allocatable :: nI(:)
        integer(n_int), allocatable :: ilut(:)
        integer :: ex(2,2), ic, run
        real(dp), allocatable :: cum_arr(:)
        real(dp) :: cum_sum

        nel = 2
        nBasis = 4
        niftot = 1
        noffflag = 0
        allocate(cum_arr(nBasis))
        get_umat_el => get_umat_test
        ElecPairs = 1
        allocate(G1(nBasis))

        allocate(nI(nel)); nI = [1,2]
        allocate(ilut(0:niftot)); call EncodeBitDet(nI, ilut)

        allocate(ilutref(0:niftot,1))
        allocate(projedet(nel,1)); projedet(:,1) = [1,2]
        call encodebitdet(projedet(:,1), ilutref(:,1))
        nmaxx = 2
        nmaxy = 2
        nmaxz = 2
        tOrbECutoff = .false.
        niftot = 1
        allocate(KPointToBasisFn(-nmaxx:nmaxx, -nmaxy:nmaxy, -nmaxz:nmaxz, 2))

        t_back_spawn_flex = .true.
        occ_virt_level = 0

        G1(1)%k = [1,0,0]
        G1(2)%k = [0,1,0]
        G1(3)%k = [1,1,0]
        G1(4)%k = [0,0,0]

        ! i also have to set the ms value in G1
        G1(1)%ms = -1
        G1(2)%ms = 1
        G1(3)%ms = -1
        G1(4)%ms = 1

        KPointToBasisFn(0,1,0,2) = 2 ! i should get kb = [0,1,0]
        KPointToBasisFn(1,0,0,1) = 3
        KPointToBasisFn(0,0,0,2) = 4
        KPointToBasisFn(1,1,0,1) = 1

        ex(1,:) = [1,2]
        ex(2,:) = [3,4]
        ic = 1
        run = 1

        print *, ""
        print *, "testing: calc_pgen_back_spawn_ueg"
        print *, "with necessary global data: "
        print *, "nel: ", nel
        print *, "nBasis: ", nbasis
        print *, "niftot: ", niftot
        print *, "encodebitdet() "
        print *, "dSFMT_init "
        call dSFMT_init(123)

        ! here ic = 1
        call assert_equals(0.0_dp, calc_pgen_back_spawn_ueg_new(nI, ilut, ex, ic, run))

        call set_flag(ilut, get_initiator_flag(run), .true.)
        ic = 2


        ! here ic = 2 and inititator
        call assert_equals(1.0_dp, calc_pgen_back_spawn_ueg_new(nI, ilut, ex, ic, run))

        ex(2,:) = [1,2]

        ilut = 0_n_int
        call set_flag(ilut, get_initiator_flag(run), .true.)

        call create_ab_list_ueg(ilut, [1,2], cum_arr, cum_sum)

        ! here it should give the same result as the ueg test.. which is
        ! 0.5.. why doesnt it?
        call assert_equals(1.0_dp/3.0_dp, calc_pgen_back_spawn_ueg_new(nI, ilut, ex, ic, run))

        ! although this test is dangerous, since actually orb_b > orb_a is
        ! enforced in the excitation generator.. but anyway..
        ex(2,:) = [2,1]

        ! first try if calc_pgen_ueg gets called correctly for inits
        call assert_equals(1.0_dp/3.0_dp, calc_pgen_back_spawn_ueg_new(nI, ilut, ex, ic, run))


        ! and now we want to test the actual back-spawn part of it..
        call set_flag(ilut, get_initiator_flag(run), .false.)

        ! but i need even more setup for that!
        ! ilutref and stuff
        ! start with
        ex(2,:) = [3,4]

        call assert_equals(1.0_dp/2.0_dp, calc_pgen_back_spawn_ueg_new(nI, ilut, ex, ic, run))

        ex(1,:) = [3,4]
        ! what if both electrons are outside?
        ! it should be the same as above..
        call assert_equals(1.0_dp/3.0_dp, calc_pgen_back_spawn_ueg_new(nI, ilut, ex, ic, run))

        ! those 2 fail now with the new implementation:
        ex(1,:) = [1,3]
        call assert_equals(0.0_dp, calc_pgen_back_spawn_ueg_new(nI, ilut, ex, ic, run))

        projedet(:,1) = [1,3]
        call assert_equals(0.0_dp, calc_pgen_back_spawn_ueg_new(nI, ilut, ex, ic, run))

        ex(1,:) = [2,4]
        call assert_equals(0.0_dp, calc_pgen_back_spawn_ueg_new(nI, ilut, ex, ic, run))

        get_umat_el => null()
        nel = -1
        niftot = -1
        ElecPairs = -1
        nBasis = -1
        deallocate(G1)
        deallocate(kpointtobasisfn)
        deallocate(ilutref)
        nmaxx = -1
        nmaxy = -1
        nmaxz = -1


    end subroutine calc_pgen_back_spawn_ueg_new_test

    function get_umat_test(i,j,k,l) result(hel)
        use constants, only: dp
        implicit none
        integer, intent(in) :: i,j,k,l
        HElement_t(dp) :: hel

        hel = 1.0_dp
    end function get_umat_test
end program test_back_spawn_excit_gen

