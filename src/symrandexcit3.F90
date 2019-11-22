#include "macros.h"

module symrandexcit3

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
                               OrbClassCount, ScratchSize1, ScratchSize2, &
                               ScratchSize3
    use SymData, only: nSymLabels
    use dSFMT_interface, only: genrand_real2_dSFMT
    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd, ClassCountInd, &
                                    CreateSingleExcit, CreateExcitLattice, &
                                    init_excit_gen_store,clean_excit_gen_store
    use FciMCData, only: pDoubles, iter, excit_gen_store_type
    use bit_reps, only: niftot, decode_bit_det_lists
    use constants, only: dp, n_int, bits_n_int
    use sym_general_mod, only: SymAllowedExcit
    use timing_neci
    use Parallel_neci
    use util_mod, only: binary_search_first_ge
    implicit none

contains

    subroutine gen_rand_excit3 (nI, ilutI, nJ, ilutJ, exFlag, IC, ExcitMat, &
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
        character(*), parameter :: this_routine = 'gen_rand_excit3'

        unused_var(part_type)

        ! Just in case
        ilutJ(0) = -1
        HElGen = 0.0_dp

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
            pDoubNew = 0
        case(2)
            pDoubNew = 1
        case(3)
            r = genrand_real2_dSFMT()
            if (r < pDoubles) then
                IC = 2
            else
                IC = 1
            endif
            pDoubNew = pDoubles
        end select

        ! Call the actual single/double excitation generators.
        if (IC == 2) then
            pGen = gen_double (nI, nJ, ExcitMat, tParity, &
                               store%ClassCountUnocc, store%virt_list)
        else
            pGen = gen_single (nI, nJ, ExcitMat, tParity, &
                               store%ClassCountOcc, store%ClassCountUnocc, &
                               store%scratch3, store%occ_list, &
                               store%virt_list)
        endif

    end subroutine


    function gen_double (nI, nJ, ExcitMat, tParity, CCUnocc, &
                         virt_list) result(pGen)

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: nJ(nel)
        integer, intent(in) :: CCUnocc(ScratchSize)
        integer, intent(in) :: virt_list(:,:)
        integer, intent(out) :: ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp) :: pGen

        real(dp) :: pElecs
        integer :: elecs(2), spn(2), orbs(2), sym_inds(2)
        integer :: sym_prod, rint, tot_pairs
        integer :: pair_list(0:nSymLabels-1)

        ! Pick and unbiased, distinct, electron pair.
        call pick_elec_pair (nI, elecs, sym_prod, spn)

        ! Pick a pair of symmetries, such that
        tot_pairs = count_orb_pairs (sym_prod, spn, CCUnocc, pair_list)

        ! If there are no possible excitations for the electron pair picked,
        ! then we abort the excitation
        if (tot_pairs == 0) then
            nJ(1) = 0
            pGen = 0.0_dp
            return
        endif

        ! Given a random number, the remainder of the generation is entirely
        ! deterministic
        rint = 1 + int(genrand_real2_dSFMT() * tot_pairs)

        ! Select a pair of symmetries to choose from
        call select_syms(rint, sym_inds, sym_prod, spn, pair_list)

        ! Select a pair of orbitals from the symmetries above.
        call select_orb_pair (rint, sym_inds, orbs, CCUnocc, virt_list)

        ! Generate the final determinant.
        call create_excit_det2 (nI, nJ, tParity, ExcitMat, elecs, orbs)

        ! Return the final probability
        pGen = pDoubNew / real(ElecPairs * tot_pairs, dp)

    end function


    subroutine pick_elec_pair (nI, elecs, sym_prod, spn)

        ! Use a triangular indexing system.
        !  --> Only need one random number to pick two distinct electrons
        !      from N(N-1)/2 pairs.
        !
        ! i.e.              21                     1
        !               32  31   ==>           3   2
        !           43  42  41             6   5   4
        !       54  53  52  51         10  9   8   7
        !
        ! We can obtain the first index, A, by considering the largest
        ! integer, i, which can give an element on that row. For an integer
        ! 1 <= i <= npair:
        !
        !   --> A = ceil((1 + sqrt(1 + 8*i)) / 2)
        !   --> B = i - (A-1)(A-2)/2

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: elecs(2), sym_prod, spn(2)

        integer :: ind, orbs(2)

        ! Generate a random integer 1 <= i <= nel(nel-1)/2 (ElecPairs)
        ind = 1 + int(ElecPairs * genrand_real2_dSFMT())

        ! Generate the two indices, and obtain the associated orbitals
        elecs(1) = ceiling((1 + sqrt(1 + 8*real(ind,dp))) / 2)
        elecs(2) = ind - ((elecs(1) - 1) * (elecs(1) - 2)) / 2
        orbs = nI(elecs)

        ! Obtain the symmetry product label
        sym_prod = RandExcitSymLabelProd (SpinOrbSymLabel(orbs(1)), &
                                          SpinOrbSymLabel(orbs(2)))

        ! Obtain spins
        spn = get_spin(orbs)

    end subroutine

    function count_orb_pairs (sym_prod, spn, CCUnocc, num_pairs) &
                              result(tot_pairs)

        integer, intent(in) :: sym_prod
        integer, intent(inout) :: spn(2)
        integer, intent(in) :: CCUnocc(ScratchSize)
        integer, intent(inout) :: num_pairs(0:nSymLabels-1)
        integer :: tot_pairs
        character(*), parameter :: this_routine = 'count_orb_pairs'

        integer :: symA, symB, indA, indB, rint, tmp_tot

        ! TODO: Are we going to be able to store this one in ScratchSize
        !       arrays as well?
        tot_pairs = 0
        if (spn(1) == spn(2)) then
            indA = spn(1)
            do symA = 0, nSymLabels - 1
                symB = RandExcitSymLabelProd(SymInvLabel(symA), sym_prod)

                if (symA == symB) then
                    tot_pairs = tot_pairs + (CCUnocc(indA) * &
                                             max(CCUnocc(indA) - 1,0)) / 2
                elseif (symB > symA) then
                    indB = ClassCountInd(spn(1), symB, -1)
                    tot_pairs = tot_pairs + (CCUnocc(indA) * CCUNocc(indB))
                endif

                num_pairs(symA) = tot_pairs
                indA = indA + 2
            enddo

        else ! spn(1) /= spn(2)

            indA = spn(1)
            do symA = 0, nSymLabels - 1
                symB = RandExcitSymLabelProd(SymInvLabel(symA), sym_prod)

                if (symA == symB) then
                    tot_pairs = tot_pairs + (CCUnocc(indA) * &
                                             CCUNocc(ieor(indA-1,1)+1))
                else
                    ! Don't restrict to A < B. Use the B < A case to count
                    ! the equivalent with the spins swapped.
                    indB = ClassCountInd(spn(2), symB, -1)
                    tot_pairs = tot_pairs + (CCUnocc(indA) * CCUnocc(indB))
                endif

                num_pairs(symA) = tot_pairs
                indA = indA + 2
            enddo
        endif

    end function

    subroutine select_syms (rint, sym_inds, sym_prod, spn, num_pairs)

        integer, intent(inout) :: rint, spn(2)
        integer, intent(out) :: sym_inds(2)
        integer, intent(in) :: sym_prod
        integer, intent(in) :: num_pairs(0:nSymLabels-1)

        integer :: syms(2), npairs, inds(2), tmp, symA

        ! Select a symA/symB pair biased by the number of possible
        ! excitations which can be made into them.
        do symA = 0, nSymLabels - 1
            if (num_pairs(symA) >= rint) then
                syms(1) = symA
                syms(2) = RandExcitSymLabelProd(SymInvLabel(symA), sym_prod)
                exit
            endif
        enddo

        ! Modify rint such that it now specifies which of the orbital pairs
        ! within the selected symmetry categories is desired.
        if (symA /= 0) rint = rint - num_pairs(symA - 1)

        ! Return the symmetry indices, rather than the symmetry labels
        ! as that is what we will need for the selections.
        sym_inds = ClassCountInd(spn, syms, -1)

    end subroutine


    subroutine select_orb_pair (rint, sym_inds, orbs, CCUnocc, &
                                virt_list)

        integer, intent(in) :: rint, sym_inds(2)
        integer, intent(in) :: CCUnocc(ScratchSize)
        integer, intent(in) :: virt_list(:,:)
        integer, intent(out) :: orbs(2)
        character(*), parameter :: this_routine = 'select_orb_pair'

        integer :: i

        if (sym_inds(1) == sym_inds(2)) then
            ! We are picking two orbitals from the same category
            ! --> Use the triangular scheme previously for selecting
            !     electrons.

            ! Select the positions of the two orbitals in the vacant list.
            orbs(1) = ceiling((1 + sqrt(1 + 8*real(rint,dp))) / 2)
            orbs(2) = rint - ((orbs(1) - 1) * (orbs(1) - 2) / 2)
        else
            ! We are picking two orbitals from different categories
            ! --> use a 'rectangular', mapping scheme.
            orbs(1) = mod(rint - 1, CCUnocc(sym_inds(1))) + 1
            orbs(2) = floor((real(rint,dp) - 1) / CCUnocc(sym_inds(1))) + 1
        endif

        ! Extract the orbitals from the vacant list.
        orbs(1) = virt_list (orbs(1), sym_inds(1))
        orbs(2) = virt_list (orbs(2), sym_inds(2))

    end subroutine


    subroutine create_excit_det2 (nI, nJ, tParity, ExcitMat, elecs, orbs)

        integer, intent(in) :: nI(nel), elecs(2), orbs(2)
        integer, intent(out) :: nJ(nel), ExcitMat(2,2)
        logical, intent(out) :: tParity

        ExcitMat(1,:) = elecs
        ExcitMat(2,:) = orbs
        nJ = nI

        ! TODO FindExcitDet (excit.F) is not modularised/interfaced yet
        call FindExcitDet(ExcitMat, nJ, 2, tParity)

    end subroutine

    subroutine construct_class_counts (nI, CCOcc, CCUnocc, pair_list)

        ! Return two arrays of length ScratchSize, containing information on
        ! the number of orbitals - occupied and unoccupied - in each symmetry.
        !
        ! The arrays are indexed via the indices returned by ClassCountInd
        ! n.b. this is O[nel], so we should store this if we can.

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: CCOcc(ScratchSize), CCUnocc(ScratchSize)
        integer, intent(out) :: pair_list(ScratchSize)

        integer :: ind_alpha, ind_beta, i, ind, tot

        CCOcc = 0
        if (tNoSymGenRandExcits) then
            ind_alpha = ClassCountInd(1,0,0)
            ind_beta = ClassCountInd(2,0,0)
            CCOcc(ind_alpha) = nOccAlpha
            CCOcc(ind_beta) = nOccBeta
        else
            do i = 1, nel
                ind = ClasscountInd(nI(i))
                CCOcc(ind) = CCOcc(ind) + 1
            enddo
        endif
        CCUnocc = OrbClassCount - CCOcc

        ! Store a -1 to indicate to the singles routine that this
        ! structure hasn't been filled in yet.
        pair_list(1) = -1

    end subroutine

    function gen_single (nI, nJ, ExcitMat,  tParity, CCOcc, CCUnocc, &
                         pair_list, occ_list, virt_list) result(pGen)

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: nJ(nel), ExcitMat(2,2)
        logical, intent(out) :: tParity
        integer, intent(in) :: CCOcc(ScratchSize), CCUnocc(ScratchSize)
        integer, intent(out) :: pair_list(ScratchSize)
        integer, intent(in) :: occ_list(:,:), virt_list(:,:)
        real(dp) :: pGen
        character(*), parameter :: this_routine = 'gen_single'

        integer :: npairs, rint, ind, src, tgt, i

        ! We still do not work with lz symmetry
        ASSERT(.not. tFixLz)

        ! Find the number of available pairs in each symmetry & overall
        if (pair_list(1) == -1) then
            pair_list(1) = CCOcc(1) * CCUnocc(1)
            do i = 2, ScratchSize
                pair_list(i) = pair_list(i-1) + (CCOcc(i) * CCUnocc(i))
            enddo
        endif
        npairs = pair_list(ScratchSize)

        ! If there are no possible singles, then abandon.
        if (npairs == 0) then
            nJ(1) = 0
            pGen = 0.0_dp
            return
        endif

        ! Pick a pair
        rint = int(1.0_dp + (genrand_real2_dSFMT() * real(npairs,dp)),sizeof_int)

        ! Select which symmetry/spin category we want.
        !ind = binary_search_first_ge (pair_list, rint)
        do ind = 1, ScratchSize
            if (pair_list(ind) >= rint) exit
        enddo
        ASSERT(ind <= ScratchSize)

        ! We are selecting one from the occupied list, and one from the
        ! unoccupied list
        ! --> There must be no overlap, so use a rectangular selection.
        if (ind > 1) rint = rint - pair_list(ind - 1)
        src = mod(rint - 1, CCOcc(ind)) + 1
        tgt = floor((real(rint,dp) - 1) / CCOcc(ind)) + 1

        ! Find the index of the src orbital in the list and the target orbital
        ! (the tgt'th vacant orbital in the given symmetry)
        ExcitMat(1,1) = occ_list(src, ind)
        ExcitMat(2,1) = virt_list(tgt, ind)

        ! Generate the new determinant
        nJ = nI
        call FindExcitDet (ExcitMat, nJ, 1, tParity)

#ifdef DEBUG_
        ! For debugging purposes only (O[N] operation).
        if (.not. SymAllowedExcit(nI, nJ, 1, ExcitMat)) &
            call stop_all(this_routine, 'Invalid excitation generated')
#endif

        ! Return the generation probability
        pGen = (1 - pDoubNew) / real(npairs, dp)

    end function



    subroutine test_sym_excit3 (nI, iterations, pDoub, exFlag)
    use SystemData, only: NEl, nBasis, G1, nBasisMax, LzTot, tUEG, &
                          tLatticeGens, tHub,tKPntSym, tFixLz
    use GenRandSymExcitNUMod, only: gen_rand_excit, ScratchSize
    Use SymData , only : nSymLabels
!    use soft_exit , only : ChangeVars
    use DetBitOps , only : EncodeBitDet, FindExcitBitDet
    use GenRandSymExcitNUMod, only: IsMomentumAllowed
    use constants, only: n_int
    use bit_reps, only: NIfTot
    use sym_mod, only: mompbcsym, GetLz
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
    INTEGER :: Lz,excitcount,ForbiddenIter,error, iter_tmp
    HElement_t(dp) :: HElGen
    type(excit_gen_store_type) :: store
    logical :: brillouin_tmp(2)
    type(timer), save :: test_timer
    character(*), parameter :: t_r = 'test_sym_excit3'

    WRITE(6,*) nI(:)
    WRITE(6,*) Iterations,pDoub,exFlag
    WRITE(6,*) "nSymLabels: ",nSymLabels
    CALL neci_flush(6)

    ! The old excitation generator will not generate singles from the HF
    ! unless tNoBrillouin is set
    brillouin_tmp(1) = tNoBrillouin
    brillouin_tmp(2) = tUseBrillouin
    tNoBrillouin = .true.
    tUseBrillouin = .false.

!Find the number of symmetry allowed excitations there should be by looking at the full excitation generator.
!Setup excit generators for this determinant
    iMaxExcit=0
    nStore(1:6)=0
    CALL GenSymExcitIt2(nI,NEl,G1,nBasis,.TRUE.,nExcitMemLen,nJ,iMaxExcit,nStore,exFlag)
    ALLOCATE(EXCITGEN(nExcitMemLen(1)),stat=ierr)
    IF(ierr.ne.0) CALL Stop_All("SetupExcitGen","Problem allocating excitation generator")
    EXCITGEN(:)=0
    CALL GenSymExcitIt2(nI,NEl,G1,nBasis,.TRUE.,EXCITGEN,nJ,iMaxExcit,nStore,exFlag)
!    CALL GetSymExcitCount(EXCITGEN,DetConn)
    excitcount=0

lp2: do while(.true.)
        CALL GenSymExcitIt2(nI,nEl,G1,nBasis,.false.,EXCITGEN,nJ,iExcit,nStore,exFlag)
        IF(nJ(1).eq.0) exit lp2
        IF(tUEG.or.tHub) THEN
            IF (IsMomentumAllowed(nJ)) THEN
                excitcount=excitcount+1
                CALL EncodeBitDet(nJ,iLutnJ)
                IF(iProcIndex.eq.0) WRITE(25,*) excitcount,iExcit,iLutnJ(0)
            ENDIF
        ELSEIF(tFixLz) THEN

            CALL GetLz(nJ,NEl,Lz)
            IF(Lz.eq.LzTot) THEN
                excitcount=excitcount+1
                CALL EncodeBitDet(nJ,iLutnJ)
                IF(iProcIndex.eq.0) WRITE(25,*) excitcount,iExcit,iLutnJ(0)
            ENDIF
        ELSEIF(tKPntSym) THEN
            IF(IsMomAllowedDet(nJ)) THEN
                excitcount=excitcount+1
                CALL EncodeBitDet(nJ,iLutnJ)
                IF(iProcIndex.eq.0) WRITE(25,*) excitcount,iExcit,nJ(:)
            ENDIF
        ELSE
            excitcount=excitcount+1
            CALL EncodeBitDet(nJ,iLutnJ)
            IF(iProcIndex.eq.0) WRITE(25,*) excitcount,iExcit,iLutnJ(0)
        ENDIF
    enddo lp2
    tNoBrillouin = brillouin_tmp(1)
    tUseBrillouin = brillouin_tmp(2)

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

    ! Build the lists we need
    store%tFilled = .false.
    store%ClassCountOcc = 0
    store%ClassCountUnocc = 0
    store%scratch3 = 0
    call decode_bit_det_lists (nI, ilut, store)

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


        IF(mod(i,400000).eq.0) THEN
            WRITE(6,"(A,I10)") "Iteration: ",i
            CALL neci_flush(6)
        ENDIF

        call gen_rand_excit3 (nI, iLut, nJ, iLutnJ, exFlag, IC, ExcitMat, &
                             tParity, pGen, HElGen, store)
        IF(nJ(1).eq.0) THEN
!            ForbiddenIter=ForbiddenIter+1
            CYCLE
        ENDIF
        IF(tKPntSym) THEN
            test=IsMomAllowedDet(nJ)
        ENDIF
        ! This is implemented for the old excitation generators, that could only handle momentum conservation under
        ! zero momentum conditions
        IF(tUEG.and.(.not.tLatticeGens)) THEN
            kx=0
            ky=0
            kz=0
            do j=1,NEl
                kx=kx+G1(nJ(j))%k(1)
                ky=ky+G1(nJ(j))%k(2)
                kz=kz+G1(nJ(j))%k(3)
            enddo
            IF(.not.(kx.eq.0.and.ky.eq.0.and.kz.eq.0)) THEN
                CYCLE
            ENDIF
        ELSEIF(tHub.and.(.not.tLatticeGens)) THEN
            kx=0
            ky=0
            kz=0
            do j=1,NEl
                kx=kx+G1(nJ(j))%k(1)
                ky=ky+G1(nJ(j))%k(2)
                kz=kz+G1(nJ(j))%k(3)
            enddo
            ktrial=(/kx,ky,0/)
            CALL MomPbcSym(ktrial,nBasisMax)
            IF(.not.(ktrial(1).eq.0.and.ktrial(2).eq.0.and.kz.eq.0)) THEN
                CYCLE
            ENDIF
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
        if (SymAllowedExcit(nI, nJ, ic, excitmat)) &
            call stop_all(t_r, 'Invalid determinant')

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
                            CALL FindExcitBitDet(iLut,iLutnJ,2,ExcitMat)
                            write(8,"(i12,f20.12,2i5,'->',2i5,2i15)") DetNum,&
                                AllDoublesHist(i,j,k,l) / (real(Iterations,dp)&
                                                        * nProcessors), &
                                i, j, k, l, iLutnJ(0),AllDoublesCount(i,j,k,l)
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
        WRITE(6,*) DetNum," Double excitations found from nI"
        OPEN(9,FILE="SinglesHist3",STATUS="UNKNOWN")
        DetNumS=0
        do i=1,nBasis
            do j=1,nBasis
                IF(AllSinglesHist(i,j).gt.0.0_dp) THEN
                    DetNumS=DetNumS+1
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
subroutine virt_uniform_sym_setup ()

    use SymExcitDataMod, only: ScratchSize, ScratchSize3
    implicit none

    ! We use the third scratch array to store data for single
    ! excitations

    call SpinOrbSymSetup ()

    ScratchSize3 = ScratchSize

end subroutine

