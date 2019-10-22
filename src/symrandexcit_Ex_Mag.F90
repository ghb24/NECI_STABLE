#include "macros.h"

module symrandexcit_Ex_mag

    ! This is another version of the excitation generators. It creates
    ! random excitations with a calculable, but non-uniform, probability.
    !
    ! Motivation:
    !     i) Use fewer random numbers than symrandexcit2
    !    ii) Generate (a bit) more uniform generation probabilities.

    use SystemData, only: nel, tFixLz, G1, ElecPairs, tUEG, tHub, &
                          tLatticeGens, tNoBrillouin, tUseBrillouin, &
                          tNoSymGenRandExcits, nOccAlpha, nOccBeta
    use SymExcitDataMod, only: ScratchSize, SpinOrbSymLabel, SymInvLabel, &
                               SymLabelList2, SymLabelCounts2, pDoubNew, &
                               pSingNew, pSing_spindiff1_new, pDoub_spindiff1_new, pDoub_spindiff2_new, &
                               OrbClassCount, ScratchSize1, ScratchSize2, &
                               ScratchSize3
    use SymData, only: nSymLabels
    use dSFMT_interface, only: genrand_real2_dSFMT
    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd, ClassCountInd, &
                                    CreateSingleExcit, CreateExcitLattice, &
                                    init_excit_gen_store,clean_excit_gen_store
    use FciMCData, only: pDoubles, pSingles, iter, excit_gen_store_type, &
                         pDoub_spindiff1, pDoub_spindiff2, pSing_spindiff1
    use bit_reps, only: niftot, decode_bit_det_lists, getExcitationType
    use constants, only: dp, n_int, bits_n_int
    use sym_general_mod, only: SymAllowedExcit
    use timing_neci
    use Parallel_neci
    use util_mod, only: binary_search_first_ge
    use symrandexcit3, only: pick_elec_pair, count_orb_pairs, select_syms, select_orb_pair, &
                             create_excit_det2, construct_class_counts
    use symexcit3, only: GenSingleExcit
    implicit none

contains


    subroutine gen_rand_excit_Ex_mag (nI, ilutI, nJ, ilutJ, exFlag, IC, ExcitMat, &
                                tParity, pGen, HElGen, store, part_type)

        ! This generator _requires_ store to have already been filled. This
        ! involves calling decode_bit_det_lists.

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        real(dp) :: r
        character(*), parameter :: this_routine = 'gen_rand_excit_Ex_Mag'
        integer :: excitType, i

        logical tAllExcitFound, tij_lt_ab_only, tSpinRestrict
        integer doubleExcitsFound


        ! Just in case
        ilutJ(0) = -1
        HElGen = 0.0_dp


!        tAllExcitFound = .false.
!        tij_lt_ab_only = .true.
!        tSpinRestrict = .true.

!        doubleExcitsFound = -1

 !       ExcitMat(1,1)=0
 !       do while(.not. tAllExcitFound)
 !           call GenSingleExcit(nI,iLutI,nJ,1,ExcitMat,tParity,tAllExcitFound,tij_lt_ab_only)
 !           doubleExcitsFound = doubleExcitsFound + 1
 !       enddo

!        write(*,*) "single excits found", doubleExcitsFound
!        call stop_all(this_routine, "found excits")





        ! UEG and Hubbard interjection for now
        ! TODO: This should be made into its own fn-pointered case.
        if ((tUEG .and. tLatticeGens) .or. (tHub .and. tLatticeGens)) then
            call CreateExcitLattice (nI, iLutI, nJ, tParity, ExcitMat, pGen)
            IC = 2
            return
        endif

        !if (.not. store%tFilled)    ** n.b. not needed anymore **
        !    call construct_class_counts (~)

        ! If exFlag is 3, select singles or doubles randomly, according
        ! to the value in pDoubles. Otherwise exFlag = 1 gives a single,
        ! and exFlag = 2 gives a double.


ASSERT(exFlag<=3.and.exFlag>=1)
        IC = exFlag
        select case(IC)
        case(1)
            ! normalise single probabilities
            pSingNew = pSingles / (pSingles + pSing_spindiff1)
            pSing_spindiff1_new = pSing_spindiff1 / (pSingles + pSing_spindiff1)
            pDoubNew = 0
            pDoub_spindiff1_new = 0
            pDoub_spindiff2_new = 0
        case(2)
            pSingNew = 0
            pSing_spindiff1_new = 0
            ! normalise double probabilities
            pDoubNew = pDoubles / (pDoubles + pDoub_spindiff1 + pDoub_spindiff2)
            pDoub_spindiff1_new = pDoub_spindiff1 / (pDoubles + pDoub_spindiff1 + pDoub_spindiff2)
            pDoub_spindiff2_new = pDoub_spindiff2 / (pDoubles + pDoub_spindiff1 + pDoub_spindiff2)
        case(3)
            pSingNew = pSingles
            pSing_spindiff1_new = pSing_spindiff1
            pDoubNew = pDoubles
            pDoub_spindiff1_new = pDoub_spindiff1
            pDoub_spindiff2_new = pDoub_spindiff2
        end select

        call select_spin_diff(excitType, IC)

        ! Call the actual single/double excitation generators.

        if (excitType==1 .or. excitType==3) then
            pGen = gen_single (nI, nJ, ilutI, ExcitMat, tParity, &
                               store, IC, excitType)

        else
            pGen = gen_double (nI, nJ, iLutI, ExcitMat, tParity, store, excitType)
        endif


        if (nJ(1)/=0 .and. excitType.ne.getExcitationType(ExcitMat, IC)) then
            write(iout,*) "NI", ni
            write(iout,*) "NJ", nj
            write(iout,*) "--- excit type wanted", excitType
            write(iout,*) "--- generated excit type", getExcitationType(ExcitMat, IC)
            write(iout,*) "--- ic", ic
            write(iout,*) pSingles
            write(iout,*) pSing_spindiff1
            write(iout,*) pDoubles
            write(iout,*) pDoub_spindiff1
            write(iout,*) pDoub_spindiff2
            write(iout,*) "first excit", excitMat(:,1)
            write(iout,*) "second excit", excitMat(:,2)

            call stop_all(this_routine, "invalid excitation generated")
        endif

ASSERT(nJ(1)==0 .or. excitType == getExcitationType(ExcitMat, IC))


    end subroutine


    subroutine select_spin_diff(excitType, IC)
        integer, intent(inout) :: excitType, IC
        real(dp) :: r, ptot


        r = genrand_real2_dSFMT()
        select case(IC)
        case(1)
            ! single excitation has been selected
            if (r<pSingles) then
                excitType = 1
            else
                excitType = 3
            endif
        case(2)
            ! double excitation has been selected
            if (r<pDoubles) then
                excitType = 2
                return
            elseif (r<pDoubles+pDoub_spindiff1) then
                excitType = 4
                return
            else
                excitType = 5
            endif
        case(3)
            ptot = pSingles
            if(r<pTot) then
                excitType = 1
                IC = 1
                return
            endif
            ptot = ptot+pSing_spindiff1
            if(r<ptot) then
                excitType = 3
                IC = 1
                return
            endif
            ptot = ptot+pDoubles
            if(r<ptot) then
                excitType = 2
                IC = 2
                return
            endif
            ptot = ptot+pDoub_spindiff1
            if(r<ptot) then
                excitType = 4
                IC = 2
                return
            endif
            excitType = 5
            IC = 2

        end select
    end subroutine


    function gen_double (nI, nJ, iLutI, ExcitMat, tParity, store, excitType) result(pGen)

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: nJ(nel)
        integer(n_int), intent(in) :: iLutI(0:niftot)
        integer, intent(out) :: ExcitMat(2,2)
        logical, intent(out) :: tParity
        integer, intent(in) :: excitType
        type(excit_gen_store_type), intent(inout), target :: store
        real(dp) :: r
        real(dp) :: pGen

        real(dp) :: pElecs
        integer :: elecs(2), elec_spn(2), virt_spn(2), orbs(2), sym_inds(2)
        integer :: sym_prod, rint, tot_pairs, realised_elecpairs
        integer :: pair_list(0:nSymLabels-1)

        if (excitType == 5) then
            ! Pick a pair of electrons with the constraint that they have a commom spin
            call pick_likespin_elec_pair(nI, elecs, sym_prod, elec_spn, realised_elecpairs, store)
        else
            ! Pick an unbiased, distinct, electron pair.
            call pick_elec_pair (nI, elecs, sym_prod, elec_spn)
            realised_elecpairs = elecpairs
        endif


        ! Pick a pair of symmetries, such that
        virt_spn = elec_spn
        select case(excitType)
        ! f(n) = 3-n maps 1->2 and 2->1
        case(4) !S->S+1
            if (nint(genrand_real2_dSFMT())==1) then
                virt_spn(1) = 3-virt_spn(1)
            else
                virt_spn(2) = 3-virt_spn(2)
            endif
        case(5) !S->S+2
            virt_spn(1) = 3-virt_spn(1)
            virt_spn(2) = 3-virt_spn(2)
        end select

        tot_pairs = count_orb_pairs (sym_prod, virt_spn, store%ClassCountUnocc, pair_list)

        ! If there are no possible excitations for the electron pair picked,
        ! then we abort the excitation
        if (tot_pairs == 0) then
            nJ(1) = 0
            pGen = 0.0_dp
            return
        endif

        r=genrand_real2_dSFMT()
        rint = 1 + int(r * tot_pairs)
        ! Given a random number, the remainder of the generation is entirely
        ! deterministic



        ! Select a pair of symmetries to choose from
        call select_syms(rint, sym_inds, sym_prod, virt_spn, pair_list)


        ! Select a pair of orbitals from the symmetries above.
        call select_orb_pair (rint, sym_inds, orbs, store%ClassCountUnocc, &
                              store%virt_list)

        ! Generate the final determinant.
        call create_excit_det2 (nI, nJ, tParity, ExcitMat, elecs, orbs)

        if (excitType==4) then
            if (virt_spn(1)==virt_spn(2) .or. (virt_spn(1)/=virt_spn(2) .and. orbs(1)==orbs(2) )) then
                tot_pairs = 2*tot_pairs
            endif
        endif


        ! Return the final probability
        select case(excitType)
        case(2)
            pGen = pDoubNew / real(realised_elecPairs * tot_pairs, dp)
        case(4)
            pGen = pDoub_spindiff1_new / real(realised_elecPairs * tot_pairs, dp)
        case(5)
            pGen = pDoub_spindiff2_new / real(realised_elecPairs * tot_pairs, dp)
        end select

    end function

    subroutine pick_likespin_elec_pair (nI, elecs, sym_prod, spn, nPairs, store)
        integer :: nPairs_alpha, nPairs_beta, nel_beta
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: elecs(2), sym_prod, spn(2), nPairs
        integer :: ind, orbs(2)
        type(excit_gen_store_type), intent(inout), target :: store

        ! Elecpairs is normally given by the number of elements in a triangular indexing system
        ! Nel * (Nel-1)/2 pairs.
        ! Here we demand that each electron in the pair has the same spin
        ! e.g. for 6 electrons we would normally have the mapping:
        !
        !                   21                     1
        !               32  31                 3   2
        !           43  42  41  ==>        6   5   4
        !       54  53  52  51         10  9   8   7
        !   65  64  63  62  61     15  14  13  12  11
        !
        ! But for 3 alpha electrons and 3 beta electrons in the parallel spin picking scheme:
        !
        !                   21              1
        !               32  31          3   2
        !                       ==>
        !                   21              4
        !               32  31          6   5
        !

        nel_beta = nel - store%nel_alpha
        nPairs_alpha = (store%nel_alpha * (store%nel_alpha-1))/2
        nPairs_beta = (nel_beta*(nel_beta-1))/2
        nPairs = nPairs_alpha + nPairs_beta

        ! Generate a random integer 1 <= i <= nPairs
        ! then shift by the appropriate pair number and use the normal mapping from random number to indices

        ind = 1 + int(nPairs * genrand_real2_dSFMT())
        if (ind <= nPairs_alpha) then
            elecs(1) = ceiling((1 + sqrt(1 + 8*real(ind,dp))) / 2)
            elecs(2) = ind - ((elecs(1) - 1) * (elecs(1) - 2)) / 2
            orbs = store%nI_alpha(elecs)
            elecs = store%nI_alpha_inds(elecs)
            ! alpha = 1
        else
            ind = ind - nPairs_alpha
            elecs(1) = ceiling((1 + sqrt(1 + 8*real(ind,dp))) / 2)
            elecs(2) = (ind - ((elecs(1) - 1) * (elecs(1) - 2)) / 2)
            orbs = store%nI_beta(elecs)
            elecs = store%nI_beta_inds(elecs)
            ! beta = 2
        endif
        spn = get_spin(orbs)
        sym_prod = RandExcitSymLabelProd (SpinOrbSymLabel(orbs(1)), &
                                          SpinOrbSymLabel(orbs(2)))

    end subroutine


    ! note: should tidy this up so that ccocc, ccunocc, pair_list, occ_list and virt_list are
    ! referenced directly from the store variable
    function gen_single (nI, nJ, iLutI, ExcitMat,  tParity, &
                         store, IC, excitType) result(pGen)

        integer, intent(in) :: nI(nel), IC
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ExcitMat(2,2)
        logical, intent(out) :: tParity
        ! pair_list is in store%scratch3
        type(excit_gen_store_type), intent(inout), target :: store
        ! encode_child interface uses value of IC in this scope,
        ! so intent is changed here to inout in order to avoid having to change the interface.
        integer, intent(inout) :: excitType
        real(dp) :: pGen
        character(*), parameter :: this_routine = 'gen_single'

        integer :: npairs, rint, ind, src, tgt, i

        ! We still do not work with lz symmetry
        ASSERT(.not. tFixLz)
        ! Find the number of available pairs in each symmetry & overall

        ! If we are changing spin, we want the number of pairs is the product
        ! of the number of occupied orbitals of class (spin, symlabel) and the
        ! number of unoccupied orbitals of class (3-spin, symlabel)


        ! This outer loop has been tentatively removed, since we now need two different
        ! pairs lists to accomodate magnetic excitations.

!        if (pair_list(1) == -1) then
            ! this is the first run for the current det
            if (excitType==1) then
                store%scratch3(1) = store%ClassCountOcc(1) * store%ClassCountUnocc(1)
                do i = 2, ScratchSize
                    store%scratch3(i) = store%scratch3(i-1) + (store%ClassCountOcc(i) * store%ClassCountUnocc(i))
                enddo
            else ! magnetic excitation
                store%scratch3(1) = store%ClassCountOcc(1) * store%ClassCountUnocc(2)
                do i = 2, ScratchSize
                    ! if i is even, we want the unocc count of class i-1
                    ! if i is odd, we want that of class i+1
                    ! mod(i,2) is either (0, 1) (i even, i odd)
                    ! 2*mod(i,2)-1 is either -1 or 1
                    ! therefore i+2*mod(i,2)-1 = i-1 (i even) or i+1 (i odd) as required
                    store%scratch3(i) = store%scratch3(i-1) + (store%ClassCountOcc(i) * store%ClassCountUnocc( i+2*mod(i,2)-1 ))
                enddo
            endif
!        endif
        npairs = store%scratch3(ScratchSize)

        ! If there are no possible singles, then abandon.
        if (npairs == 0) then
            nJ(1) = 0
            pGen = 0.0_dp
            return
        endif

        ! Pick a pair
        rint = int(1.0_dp + (genrand_real2_dSFMT() * real(npairs,dp)),sizeof_int)

        ! Select which symmetry/spin category we want for the currently  occupied orbital
        !ind = binary_search_first_ge (pair_list, rint)
        do ind = 1, ScratchSize
            if (store%scratch3(ind) >= rint) exit
        enddo
        ASSERT(ind <= ScratchSize)

        ! We are selecting one from the occupied list, and one from the
        ! unoccupied list
        ! --> There must be no overlap, so use a rectangular selection.
        if (ind > 1) rint = rint - store%scratch3(ind - 1)
        ASSERT(1.le. rint .and. rint.le.store%scratch3(ind))
        ! we now have a random number between 1 and the number of pairs in the selected class
        ! virt_list has the indices
        !src = mod(rint - 1, CCOcc(ind)) + 1

        src = mod((rint-1),store%ClassCountOcc(ind)) + 1
        tgt = (rint-1)/store%ClassCountOcc(ind) + 1

        ! Find the index of the src orbital in the list and the target orbital
        ! (the tgt'th vacant orbital in the given symmetry)
        ExcitMat(1,1) = store%occ_list(src, ind)
        if (excitType==1) then
            ExcitMat(2,1) = store%virt_list(tgt, ind)
        elseif (excitType==3) then
            ExcitMat(2,1) = store%virt_list(tgt, ind+2*mod(ind,2)-1)
        endif

        ! Generate the new determinant
        nJ = nI
        call FindExcitDet (ExcitMat, nJ, 1, tParity)
#ifdef __DEBUG
        ! For debugging purposes only (O[N] operation).
        if (.not. SymAllowedExcit(nI, nJ, 1, ExcitMat) .or. ExcitMat(1,1)*ExcitMat(2,1)==0) then
            write(iout,*) "ccocc(1)", store%ClassCountOcc(1)
            write(iout,*) "ccocc(2)", store%ClassCountOcc(2)
            write(iout,*) "ccunocc(1)", store%ClassCountUnocc(1)
            write(iout,*) "ccunocc(2)", store%ClassCountUnocc(2)
            write(iout,*) "ind", ind
            write(iout,*) "pair_list", store%scratch3
            write(iout,*) "alpha count", store%nel_alpha
            write(iout,*) "src", src, "tgt", tgt
            write(iout,*) "excitmat", excitmat(1,:)
            call stop_all(this_routine, "invalid single excitation generated")
        endif
#endif
        ! Return the generation probability
        if (excitType==1) then
            pGen = pSingNew / real(npairs, dp)
        else
            pGen = pSing_spindiff1_new / real(npairs, dp)
        endif

    end function



    subroutine test_sym_excit_ExMag (nI, iterations)

    use SystemData, only: NEl, nBasis, G1, nBasisMax, LzTot, tUEG, &
                          tLatticeGens, tHub,tKPntSym, tFixLz
    use GenRandSymExcitNUMod, only: gen_rand_excit, ScratchSize
    Use SymData , only : nSymLabels
!    use soft_exit , only : ChangeVars
    use DetBitOps , only : EncodeBitDet, FindExcitBitDet
    use GenRandSymExcitNUMod, only: IsMomentumAllowed
    use constants, only: n_int
    use bit_reps, only: NIfTot,decode_bit_det_spinsep
    use sym_mod, only: mompbcsym, GetLz
    use SymExcit4, only: CountExcitations4
    use neci_intfce
    IMPLICIT NONE
    INTEGER :: i,Iterations,exFlag,nI(NEl),nJ(NEl),IC,ExcitMat(2,2),kx,ky,kz,ktrial(3)
    real(dp) :: pDoub,pGen,AverageContrib,AllAverageContrib
    INTEGER(KIND=n_int) :: iLutnJ(0:NIfTot),iLut(0:NIfTot)
    INTEGER :: iExcit
    LOGICAL :: tParity,IsMomAllowedDet,test

    ! Accumulator arrays. These need to be allocated on the heap, or we
    ! get a segfault by overflowing the stack using ifort
    real(dp), allocatable :: DoublesHist(:,:,:,:)
    real(dp), allocatable :: AllDoublesHist(:,:,:,:)
    real(dp), allocatable :: SinglesHist(:,:)
    real(dp), allocatable :: AllSinglesHist(:,:)
    integer, allocatable :: DoublesCount(:,:,:,:)
    integer, allocatable :: AllDoublesCount(:,:,:,:)
    integer, allocatable :: SinglesCount(:,:)
    integer, allocatable :: AllSinglesCount(:,:)

    INTEGER , ALLOCATABLE :: EXCITGEN(:)
    INTEGER :: ierr,Ind1,Ind2,Ind3,Ind4,iMaxExcit,nStore(6),nExcitMemLen(1),j,k,l,DetNum,DetNumS
    integer :: DetNumD0, DetNumD1, DetNumD2, DetNumS1, DetNumS0
    INTEGER :: Lz,excitcount,ForbiddenIter,error, iter_tmp,nSing,nSing_1,nDoub,nDoub_1,nDoub_2
    HElement_t(dp) :: HElGen
    type(excit_gen_store_type) :: store
    logical :: brillouin_tmp(2)
    type(timer), save :: test_timer
    character(*), parameter :: t_r = 'test_sym_excit_ExMag'

    WRITE(6,*) nI(:)
    WRITE(6,*) Iterations
    WRITE(6,*) "nSymLabels: ",nSymLabels
    CALL neci_flush(6)

    call CountExcitations4(nI, 1, 1, 0, 0, nSing)
    call CountExcitations4(nI, 1, 1, 1, 1, nSing_1)
    call CountExcitations4(nI, 1, 1, 0, 0, nDoub)
    call CountExcitations4(nI, 1, 1, 1, 1, nDoub_1)
    call CountExcitations4(nI, 1, 1, 2, 2, nDoub_2)

    write(6,*) "nSing: ",nSing
    write(6,*) "nSing_1: ",nSing_1
    write(6,*) "nDoub: ",nDoub
    write(6,*) "nDoub_1: ",nDoub_1
    write(6,*) "nDoub_2: ",nDoub_2

    excitcount = nSing + nSing_1 + nDoub + nDoub_1 + nDoub_2

    pSingles = real(nSing,dp)/real(excitcount,dp)
    pSing_spindiff1 = real(nSing_1,dp)/real(excitcount,dp)
    pDoubles = real(nDoub,dp)/real(excitcount,dp)
    pDoub_spindiff1 = real(nDoub_1,dp)/real(excitcount,dp)
    pDoub_spindiff2 = real(nDoub_2,dp)/real(excitcount,dp)

    WRITE(6,*) "Determinant has ",excitcount," total excitations from it."
    CALL neci_flush(6)

    ! Allocate the accumulators
    allocate (DoublesHist(nbasis, nbasis, nbasis, nbasis))
    allocate (AllDoublesHist(nbasis, nbasis, nbasis, nbasis))
    allocate (SinglesHist(nbasis, nbasis))
    allocate (AllSinglesHist(nbasis, nbasis))
    allocate (DoublesCount(nbasis, nbasis, nbasis, nbasis))
    allocate (AllDoublesCount(nbasis, nbasis, nbasis, nbasis))
    allocate (SinglesCount(nbasis, nbasis))
    allocate (AllSinglesCount(nbasis, nbasis))

    ! Initialise the excitation generator store
    call init_excit_gen_store (store)

    ! Zero the accumulators
    DoublesHist = 0
    SinglesHist = 0
    AllDoublesHist = 0
    AllSinglesHist = 0
    DoublesCount = 0
    SinglesCount = 0
    AllDoublesCount = 0
    AllSinglesCount = 0

    CALL EncodeBitDet(nI, iLut)
    CALL decode_bit_det_spinsep(nJ,iLut,store)

    do i = 1,NEl
        if(nI(i).ne.nJ(i)) stop 'Error with det'
    enddo

    AverageContrib=0.0_dp
    AllAverageContrib=0.0_dp
    ForbiddenIter=0
!    pDoub=1.0_dp
!    IF(iProcIndex.eq.0) OPEN(9,FILE="AvContrib",STATUS="UNKNOWN")

    test_timer%timer_name = 'test_symrandexcit3'
    call set_timer (test_timer)
    iter_tmp = iter
    do i=1,Iterations
        iter = i

        IF(mod(i,1).eq.0) THEN
        !IF(mod(i,40000).eq.0) THEN
            !WRITE(6,"(A,I10)") "Iteration: ",i
            CALL neci_flush(6)
        ENDIF

        call gen_rand_excit_Ex_Mag (nI, iLut, nJ, iLutnJ, 3, IC, ExcitMat, &
                             tParity, pGen, HElGen, store)
        IF(nJ(1).eq.0) THEN
!            ForbiddenIter=ForbiddenIter+1
            CYCLE
        ENDIF
        AverageContrib=AverageContrib+1.0_dp/pGen

!        CALL EncodeBitDet(nJ,iLutnJ)
!        IF(IC.eq.1) THEN
!            WRITE(6,*) ExcitMat(1,1),ExcitMat(2,1)
!        ELSE
!            WRITE(6,*) "Double Created"
!            WRITE(6,*) ExcitMat(1,1),ExcitMat(1,2),ExcitMat(2,1),ExcitMat(2,2)
!        ENDIF

        IF(IC.eq.1) THEN
            SinglesHist(ExcitMat(1,1),ExcitMat(2,1))=SinglesHist(ExcitMat(1,1),ExcitMat(2,1))+(1.0_dp/pGen)
            SinglesCount(ExcitMat(1,1), ExcitMat(2,1)) = &
                SinglesCount(ExcitMat(1,1), ExcitMat(2,1)) + 1
!            SinglesNum(ExcitMat(1,1),ExcitMat(2,1))=SinglesNum(ExcitMat(1,1),ExcitMat(2,1))+1
        ELSE
!Have to make sure that orbitals are in the same order...
            IF(ExcitMat(1,1).gt.ExcitMat(1,2)) THEN
                Ind1=ExcitMat(1,2)
                Ind2=ExcitMat(1,1)
            ELSE
                Ind1=ExcitMat(1,1)
                Ind2=ExcitMat(1,2)
            ENDIF
            IF(ExcitMat(2,1).gt.ExcitMat(2,2)) THEN
                Ind3=ExcitMat(2,2)
                Ind4=ExcitMat(2,1)
            ELSE
                Ind3=ExcitMat(2,1)
                Ind4=ExcitMat(2,2)
            ENDIF
            DoublesHist(Ind1,Ind2,Ind3,Ind4)=DoublesHist(Ind1,Ind2,Ind3,Ind4)+(1.0_dp/pGen)
            DoublesCount(ind1,ind2,ind3,ind4) = &
                DoublesCount(ind1,ind2,ind3,ind4) + 1
        ENDIF
!        IF(mod(i,iWriteEvery).eq.0) THEN
!            AllAverageContrib=0.0_dp
!#ifdef PARALLEL
!            CALL MPI_AllReduce(AverageContrib,AllAverageContrib,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
!#else
!            AllAverageContrib=AverageContrib
!#endif
!            IF(iProcIndex.eq.0) THEN
!                WRITE(9,*) i,AllAverageContrib/(REAL(i,8)*excitcount*nProcessors)
!            ENDIF
!!            CALL ChangeVars(tDummy,tSoftExitFound,tDummy2)
!!            IF(tSoftExitFound) EXIT
!        ENDIF

!Check excitation
        if (.not.SymAllowedExcit(nI, nJ, ic, excitmat)) then
            write(6,*) "nI: ",nI
            write(6,*) "nJ: ",nJ
            write(6,*) "IC: ",IC
            write(6,*) "electrons ",excitmat(1,:)
            write(6,*) "holes: ",excitmat(2,:)
            call stop_all(t_r, 'Invalid determinant')
        endif

    enddo
    iter = iter_tmp
    call halt_timer (test_timer)

!    IF(iProcIndex.eq.0) CLOSE(9)

#ifdef PARALLEL
    call MPIBarrier(error)
    call MPIAllReduce (DoublesHist, MPI_SUM, AllDoublesHist)
    call MPIAllReduce (SinglesHist, MPI_SUM, AllSinglesHist)
    call MPIAllReduce (DoublesCount, MPI_SUM, AllDoublesCount)
    call MPIAllReduce (SinglesCount, MPI_SUM, AllSinglesCount)
#else
    AllDoublesHist = DoublesHist
    AllSinglesHist = SinglesHist
    AllDoublesCount = DoublesCount
    AllSinglesCount = SinglesCount
#endif

!Now run through arrays normalising them so that numbers are more managable.
    IF(iProcIndex.eq.0) THEN
        OPEN(8,FILE="DoublesHist3",STATUS="UNKNOWN")
        DetNum=0
        DetNumD0=0
        DetNumD1=0
        DetNumD2=0
        do i=1,nBasis-1
            do j=i+1,nBasis
                do k=1,nBasis-1
                    do l=k+1,nBasis
                        IF(AllDoublesHist(i,j,k,l).gt.0.0_dp) THEN
    !                        DoublesHist(i,j,k,l)=DoublesHist(i,j,k,l)/real(Iterations,8)

                            DetNum=DetNum+1

                            ExcitMat(1,1)=i
                            ExcitMat(1,2)=j
                            ExcitMat(2,1)=k
                            ExcitMat(2,2)=l


                            select case(getExcitationType(excitMat, 2))
                            case(2)
                                DetNumD0 = DetNumD0+1
                            case(4)
                                DetNumD1 = DetNumD1+1
                            case(5)
                                DetNumD2 = DetNumD2+1
                            end select



                            CALL FindExcitBitDet(iLut,iLutnJ,2,ExcitMat)
                            write(8,"(i12,f20.12,2i5,'->',2i5,2i15,i12)") DetNum,&
                                AllDoublesHist(i,j,k,l) / (real(Iterations,dp)&
                                                        * nProcessors), &
                                i, j, k, l, iLutnJ(0),AllDoublesCount(i,j,k,l),getExcitationType(excitMat,2)
!                            WRITE(6,*) DetNum,DoublesHist(i,j,k,l),i,j,"->",k,l
                            IF(tHub.or.tUEG) THEN
                                write(8,*) "#",G1(i)%k(1),G1(i)%k(2)
                                write(8,*) "#",G1(j)%k(1),G1(j)%k(2)
                                write(8,*) "#",G1(k)%k(1),G1(k)%k(2)
                                write(8,*) "#",G1(l)%k(1),G1(l)%k(2)
                            ENDIF
                        ENDIF
                    enddo
                enddo
            enddo
        enddo
        CLOSE(8)
        WRITE(6,*) DetNum," Total double excitations found from nI"
        WRITE(6,*) DetNumD0," Double spindiff0 excitations found from nI"
        WRITE(6,*) DetNumD1," Double spindiff1 excitations found from nI"
        WRITE(6,*) DetNumD2," Double spindiff2 excitations found from nI"
        OPEN(9,FILE="SinglesHist3",STATUS="UNKNOWN")
        DetNumS=0
        DetNumS0=0
        DetNumS1=0
        do i=1,nBasis
            do j=1,nBasis
                IF(AllSinglesHist(i,j).gt.0.0_dp) THEN
                    DetNumS=DetNumS+1
                    select case(getExcitationType(excitMat, 1))
                    case(1)
                        DetNumS0 = DetNumS0+1
                    case(3)
                        DetNumS1 = DetNumS1+1
                    end select
                    ExcitMat(1,1)=i
                    ExcitMat(2,1)=j
                    CALL FindExcitBitDet(iLut,iLutnJ,1,ExcitMat)
                    write(9,*) DetNumS, AllSinglesHist(i,j) / &
                                        (real(Iterations,dp) * nProcessors), &
                               i, "->", j, ALlSinglesCount(i, j)
!                    WRITE(6,*) DetNumS,AllSinglesHist(i,j),i,"->",j
                ENDIF
            enddo
        enddo
        CLOSE(9)
        WRITE(6,*) DetNumS," Single excitations found from nI"
        WRITE(6,*) DetNumS0," Single spindiff0 excitations found from nI"
        WRITE(6,*) DetNumS1," Single spindiff1 excitations found from nI"

        WRITE(6,*) DetNumS+DetNum," Total excitations found from nI"
        IF((DetNum+DetNumS).ne.ExcitCount) THEN
            CALL construct_class_counts(nI, store%ClassCountOcc, &
                                            store%ClassCountUnocc, &
                                            store%scratch3)
            WRITE(6,*) "Total determinants = ", ExcitCount
            WRITE(6,*) "ClassCount2(:)= ", store%ClassCountOcc
            WRITE(6,*) "***"
            WRITE(6,*) "ClassCountUnocc2(:)= ", store%ClassCountUnocc
            CALL Stop_All("TestGenRandSymExcitNU","Not all excitations accounted for...")
        ENDIF
    ENDIF

    CALL MPIBarrier(error)

    ! Deallocate the accumulators
    deallocate (DoublesHist, AllDoublesHist, &
                SinglesHist, AllSinglesHist, &
                DoublesCount, AllDoublesCount, &
                SinglesCount, AllSinglesCount)
    call clean_excit_gen_store (store)

    END SUBROUTINE

end module



! N.B. This is outside the module *sigh*
!subroutine virt_uniform_sym_setup ()

!    use SymExcitDataMod, only: ScratchSize, ScratchSize3
!    implicit none

    ! We use the third scratch array to store data for single
    ! excitations

!    call SpinOrbSymSetup ()

!    ScratchSize3 = ScratchSize

!end subroutine

