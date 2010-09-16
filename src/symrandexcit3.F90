#include "macros.h"

module symrandexcit3

    ! This is another version of the excitation generators. It creates 
    ! random excitations with a calculable, but non-uniform, probability.
    !
    ! Motivation:
    !     i) Use fewer random numbers than symrandexcit2
    !    ii) Generate (a bit) more uniform generation probabilities.

    use SystemData, only: nel, tFixLz, G1, ElecPairs, tUEG, tHub, &
                          tLatticeGens, tNoBrillouin, tUseBrillouin
    use SymExcitDataMod, only: ScratchSize, SpinOrbSymLabel, SymInvLabel, &
                               SymLabelList2, SymLabelCounts2, pDoubNew
    use SymData, only: nSymLabels
    use dSFMT_interface, only: genrand_real2_dSFMT
    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd, ClassCountInd, &
                                    construct_class_counts,CreateSingleExcit,&
                                    CreateExcitLattice
    use FciMCData, only: pDoubles, iter
    use bit_reps, only: niftot
    use constants, only: dp, n_int, bits_n_int
    use Determinants, only: write_det
    implicit none

contains

    subroutine gen_rand_excit3 (nI, ilutI, nJ, ilutJ, exFlag, IC, ExcitMat, &
                                tParity, pGen, HElGen, tFilled, CCOcc, &
                                CCUnocc, pair_list)

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), IC, ExcitMat(2,2)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
        integer, intent(inout) :: CCOcc(ScratchSize)
        integer, intent(inout) :: CCUnocc(ScratchSize)
        integer, intent(inout) :: pair_list(ScratchSize)
        logical, intent(inout) :: tFilled
        HElement_t, intent(out) :: HElGen

        real(dp) :: r
        character(*), parameter :: this_routine = 'gen_rand_excit3'

        ! Just in case
        ilutJ(0) = -1

        ! UEG and Hubbard interjection for now
        ! TODO: This should be made into its own fn-pointered case.
        if ((tUEG .and. tLatticeGens) .or. (tHub .and. tLatticeGens)) then
            call CreateExcitLattice (nI, iLutI, nJ, tParity, ExcitMat, pGen)
            IC = 2
            return
        endif

        ! Count occupied/unoccupied orbitals in each symmetry class. This
        ! is an O[nel] operation. For efficiency, store these arrays
        ! between invocations of the excitation generator.
        if (.not. tFilled) then
            call construct_class_counts (nI, CCOcc, CCUnocc)
            tFilled = .true.
        endif

        ! If exFlag is 3, select singles or doubles randomly, according
        ! to the value in pDoubles. Otherwise exFlag = 1 gives a single,
        ! and exFlag = 2 gives a double.
        ASSERT(exFlag <= 3 .and. exFlag >= 1)
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
            pGen = gen_double (nI, nJ, iLutI, ExcitMat, tParity, CCUnocc, &
                               pair_list)
        else
            pGen = gen_single (nI, nJ, ilutI, ExcitMat, tParity, CCocc, &
                               CCUnocc, pair_list)
        endif

    end subroutine


    function gen_double (nI, nJ, iLutI, ExcitMat, tParity, CCUnocc, &
                         pair_list) result(pGen)

        integer, intent(in) :: nI(nel)
        integer, intent(out) :: nJ(nel)
        integer, intent(in) :: CCUnocc(ScratchSize)
        integer, intent(inout) :: pair_list(ScratchSize)
        integer(n_int), intent(in) :: iLutI(0:niftot)
        integer, intent(out) :: ExcitMat(2,2)
        logical, intent(out) :: tParity
        real(dp) :: pGen

        real(dp) :: pElecs
        integer :: elecs(2), spn(2), orbs(2), sym_inds(2)
        integer :: sumMl, sym_prod, rint, tot_pairs

        ! Pick and unbiased, distinct, electron pair.
        pElecs = pick_elec_pair (nI, elecs, sym_prod, spn, sumMl)

        ! Pick a pair of symmetries, such that 
        tot_pairs = count_orb_pairs (sym_prod, spn, sumMl, CCUnocc, pair_list)

        ! If there are no possible excitations for the electron pair picked, 
        ! then we abort the excitation
        if (tot_pairs == 0) then
            nJ(1) = 0
            return
        endif

        ! Given a random number, the remainder of the generation is entirely
        ! deterministic
        rint = 1 + int(genrand_real2_dSFMT() * tot_pairs)

        ! Select a pair of symmetries to choose from
        call select_syms(rint, sym_inds, sym_prod, spn, sumMl, CCUnocc, &
                         pair_list)

        ! Select a pair of orbitals from the symmetries above.
        call select_orb_pair (rint, sym_inds, ilutI, orbs, CCUnocc)

        ! Generate the final determinant.
        call create_excit_det2 (nI, nJ, tParity, ExcitMat, elecs, orbs)
       
        ! Return the final probability
        pGen = pDoubNew * pElecs / real(tot_pairs, dp)

    end function


    function pick_elec_pair (nI, elecs, sym_prod, spn, sumMl) result(pElec)

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
        integer, intent(out) :: elecs(2), sym_prod, spn(2), sumMl
        real(dp) :: pElec

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

        ! Store the Lz value to preserve
        if (tFixLz) sumMl = product(G1(orbs)%Ml)

        ! Return the generation probability
        pElec = 2 / real(nel*(nel-1), dp)

    end function
    
    function count_orb_pairs (sym_prod, spn, sumMl, CCUnocc, num_pairs) &
                              result(tot_pairs)

        integer, intent(in) :: sym_prod, spn(2), sumMl
        integer, intent(in) :: CCUnocc(ScratchSize)
        integer, intent(inout) :: num_pairs(0:nSymLabels-1)
        integer :: tot_pairs
        character(*), parameter :: this_routine = 'count_orb_pairs'

        integer :: symA, symB, indA, indB, rint, tmp_tot

        ! TODO: Are we going to be able to store this one in ScratchSize
        !       arrays as well?
        num_pairs  = 0
        do symA = 0, nSymLabels - 1

            ! Get the required B symmetry to complement the chosen A sym.
            symB = RandExcitSymLabelProd(SymInvLabel(symA), sym_prod)

            ! We must avoid double counting A/B pairs
            if (symB < symA) cycle

            ASSERT(.not. tFixLz)
            indA = ClassCountInd(spn(1), symA, -1)
            indB = ClassCountInd(spn(2), symB, -1)
            if (symA == symB) then
                if (CCUnocc(indA) > 0 .and. CCUnocc(indB) > 0) then
                    if (spn(1) /= spn(2)) then
                        num_pairs(symA) = CCUnocc(indA) * CCUnocc(indB)
                    else
                        num_pairs(symA) = (CCUnocc(indA) * &
                                           (CCUnocc(indB) - 1)) / 2
                    endif
                endif
            else
                if (CCUnocc(indA) > 0 .and. CCUnocc(indB) > 0) &
                    num_pairs(symA) = CCUnocc(indA) * CCUnocc(indB)
                if (spn(1) /= spn(2)) then
                    indA = ClassCountInd(spn(2), symA, -1)
                    indB = ClassCountInd(spn(1), symB, -1)
                    num_pairs(symA) = num_pairs(symA) + &
                                      (CCUnocc(indA) * CCUnocc(indB)) 
                endif
            endif
        enddo

        ! Count the total number of pairs.
        tot_pairs = sum(num_pairs)

    end function

    subroutine select_syms (rint, sym_inds, sym_prod, spn, sumMl, CCUnocc, &
                            num_pairs)

        integer, intent(inout) :: rint, spn(2)
        integer, intent(out) :: sym_inds(2)
        integer, intent(in) :: sumMl, sym_prod
        integer, intent(in) :: CCUnocc(ScratchSize), num_pairs(0:nSymLabels-1)

        integer :: tmp_tot, syms(2), npairs, inds(2), tmp, symA

        ! Select a symA/symB pair biased by the number of possible 
        ! excitations which can be made into them.
        tmp_tot = 0
        do symA = 0, nSymLabels - 1
            if (tmp_tot + num_pairs(symA) >= rint) then
                syms(1) = symA
                syms(2) = RandExcitsymLabelProd(SymInvLabel(symA), sym_prod)
                exit
            endif
            tmp_tot = tmp_tot + num_pairs(symA)
        enddo

        ! Modify rint such that it now specifies which of the orbital pairs
        ! within the selected symmetry categories is desired.
        rint = rint - tmp_tot

        ! If there are more than one symmetry index corresponding to these
        !
        ! n.b. need to look in SymLabelList2(SymLabelCounts2(1,Ind)+i)
        !      which is a list of length SymLabelCounts2(2,ind)
        if (spn(1) /= spn(2) .and. syms(1) /= syms(2)) then
            inds = ClassCountInd(spn, syms, -1)
            npairs = CCUnocc(inds(1)) * CCUnocc(inds(2))
            if (rint > npairs) then
                ! Swap the spins around --> the other possibility, and adjust
                ! the random number to select a pair in this set.
                tmp = spn(1)
                spn(1) = spn(2)
                spn(2) = tmp
                rint = rint - npairs
            endif
        endif

        ! Return the symmetry indices, rather than the symmetry labels
        ! as that is what we will need for the selections.        
        sym_inds = ClassCountInd(spn, syms, -1)

    end subroutine


    subroutine select_orb_pair (rint, sym_inds, ilutI, orbs, CCUnocc)

        integer, intent(in) :: rint, sym_inds(2)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(in) :: CCUnocc(ScratchSize)
        integer :: orbs(2)
        character(*), parameter :: this_routine = 'select_orb_pair'

        integer :: a, b, i, nvac, pos, orb, offset, norbs

        if (sym_inds(1) == sym_inds(2)) then
            ! We are picking two orbitals from the same category
            ! --> Use the triangular scheme previously for selecting
            !     electrons.

            ! Select the positions of the two orbitals in the vacant list.
            a = ceiling((1 + sqrt(1 + 8*real(rint,dp))) / 2)
            b = rint - ((a - 1) * (a - 2) / 2)
            orbs(1) = min(a, b)
            orbs(2) = max(a, b)

            nvac = 0
            pos = 1
            offset = SymLabelCounts2(1, sym_inds(1))
            norbs = SymLabelCounts2(2, sym_inds(1))
            do i = 0, norbs
                orb = SymLabelList2(offset + i)
                if (IsNotOcc(ilutI, orb)) then
                    nvac = nvac + 1
                    if (nvac == orbs(pos)) then
                        orbs(pos) = orb
                        pos = pos + 1
                        if (pos > 2) exit
                    endif
                endif
            enddo

            if (i == norbs) &
                call stop_all (this_routine, 'Unable to find enough vacant &
                                             &orbitals in SymLabelList2')
        else
            ! We are picking two orbitals from different categories
            ! --> use a 'rectangular', mapping scheme.
            orbs(1) = mod(rint - 1, CCUnocc(sym_inds(1))) + 1
            orbs(2) = floor((real(rint,dp) - 1) / CCUnocc(sym_inds(1))) + 1

            do pos = 1, 2
                nvac = 0
                offset = SymLabelCounts2(1, sym_inds(pos))
                norbs = SymLabelCounts2(2, sym_inds(pos))
                do i = 0, norbs-1
                    orb = SymLabelList2(offset + i)
                    if (IsNotOcc(ilutI, orb)) then
                        nvac = nvac + 1
                        if (nvac == orbs(pos)) then
                            orbs(pos) = orb
                            exit
                        endif
                    endif
                enddo

                if (i == norbs) &
                    call stop_all (this_routine, 'Unable to find enough &
                                          &vacant orbitals in SymLabelList2')
            enddo

        endif

    end subroutine


    subroutine create_excit_det2 (nI, nJ, tParity, ExcitMat, elecs, orbs)

        integer, intent(in) :: nI(nel), elecs(2), orbs(2)
        integer, intent(out) :: nJ(nel), ExcitMat(2,2)
        logical, intent(out) :: tParity

        ! TODO: Update the bit representations here as well?
        ExcitMat(1,:) = elecs
        ExcitMat(2,:) = orbs
        nJ = nI

        ! TODO FindExcitDet (excit.F) is not modularised/interfaced yet
        call FindExcitDet(ExcitMat, nJ, 2, tParity)

    end subroutine

    function gen_single (nI, nJ, iLutI, ExcitMat,  tParity, CCOcc, CCUnocc, &
                         pair_list) result(pGen)

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ExcitMat(2,2)
        logical, intent(out) :: tParity
        integer, intent(in) :: CCOcc(ScratchSize), CCUnocc(ScratchSize)
        integer, intent(out) :: pair_list(ScratchSize)
        real(dp) :: pGen
        character(*), parameter :: this_routine = 'gen_single'

        integer :: npairs, rint, ind, tot, src, tgt, i, cnt
        integer :: offset, norbs, orb

        ! We still do not work with lz symmetry
        ASSERT(.not. tFixLz)

        ! Find the number of available pairs
        pair_list = CCOcc * CCUnocc
        npairs = sum(pair_list)

        ! If there are no possible singles, then abandon.
        if (npairs == 0) then
            nJ(1) = 0
            return
        endif

        ! Pick a pair
        rint = 1 + (genrand_real2_dSFMT() * npairs)
        
        ! Select which symmetry/spin category we want.
        tot = 0
        do ind = 1, ScratchSize
            if (tot + pair_list(ind) >= rint) exit
            tot = tot + pair_list(ind)
        enddo
        ASSERT(ind <= ScratchSize) ! Ensure we haven't overflowed.

        ! We are selecting one from the occupied list, and one from the
        ! unoccupied list
        ! --> There must be no overlap, so use a rectangular selection.
        rint = rint - tot
        src = mod(rint - 1, CCOcc(ind)) + 1
        tgt = floor((real(rint,dp) - 1) / CCOcc(ind)) + 1

        ! Find the index of the src orbital in the list
        cnt = 0
        do i = 1, nel
            if (ClassCountInd(nI(i)) == ind) then
                cnt = cnt + 1
                if (cnt == src) exit
            endif
        enddo
        ASSERT(i <= nel) ! Ensure we haven't overflowed.
        src = i

        ! Find the target orbital, as the tgt'th vacant orbital under the
        ! given symmetry index.
        cnt = 0
        offset = SymLabelCounts2(1, ind)
        norbs = SymLabelCounts2(2, ind)
        do i = 0, norbs - 1
            orb = SymLabelList2(offset + i)
            if (IsNotOcc(ilutI, orb)) then
                cnt = cnt + 1
                if (cnt == tgt) then
                    tgt = orb
                    exit
                endif
            endif
        enddo
        ASSERT(i < norbs) ! Ensure we haven't overflowed.

        ! Generate the new determinant
        nJ = nI
        ExcitMat(1,1) = src
        ExcitMat(2,1) = tgt
        call FindExcitDet (ExcitMat, nJ, 1, tParity)

#ifdef __DEBUG
        ! For debugging purposes only (O[N] operation).
        call IsSymAllowedExcit (nI, nJ, 1, ExcitMat)
#endif

        ! Return the generation probability
        pGen = (1 - pDoubNew) / real(npairs, dp)

    end function



    subroutine test_sym_excit3 (nI, iterations, pDoub, exFlag)
    use SystemData, only: NEl, nBasis, G1, nBasisMax, LzTot, tUEG, &
                          tLatticeGens, tHub,tKPntSym, tFixLz
    use GenRandSymExcitNUMod, only: gen_rand_excit, construct_class_counts,ScratchSize
    Use SymData , only : nSymLabels
    use Parallel
!    use soft_exit , only : ChangeVars 
    use DetBitOps , only : EncodeBitDet, FindExcitBitDet
    use GenRandSymExcitNUMod, only: IsMomentumAllowed
    use constants, only: n_int
    use bit_reps, only: NIfTot
    use sym_mod, only: mompbcsym, GetLz
    IMPLICIT NONE
    INTEGER :: i,Iterations,exFlag,nI(NEl),nJ(NEl),IC,ExcitMat(2,2),kx,ky,kz,ktrial(3)
    REAL*8 :: pDoub,pGen,AverageContrib,AllAverageContrib
    INTEGER :: ClassCount2(ScratchSize),Scratch1(ScratchSize),Scratch2(ScratchSize),scratch3(scratchsize)
    INTEGER(KIND=n_int) :: iLutnJ(0:NIfTot),iLut(0:NIfTot)
    INTEGER :: ClassCountUnocc2(ScratchSize),iExcit
    LOGICAL :: tParity,tFilled,IsMomAllowedDet,test
    REAL*8 , ALLOCATABLE :: DoublesHist(:,:,:,:),SinglesHist(:,:),AllDoublesHist(:,:,:,:),AllSinglesHist(:,:)
    INTEGER , ALLOCATABLE :: EXCITGEN(:)
    INTEGER :: ierr,Ind1,Ind2,Ind3,Ind4,iMaxExcit,nStore(6),nExcitMemLen,j,k,l,DetNum,DetNumS,Lz,excitcount,ForbiddenIter,error, iter_tmp
    HElement_t :: HElGen
    logical :: brillouin_tmp(2)

    WRITE(6,*) nI(:)
    WRITE(6,*) Iterations,pDoub,exFlag
    WRITE(6,*) "nSymLabels: ",nSymLabels
    CALL FLUSH(6)

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
    ALLOCATE(EXCITGEN(nExcitMemLen),stat=ierr)
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
    CALL FLUSH(6)

!Allocate memory for histogramming determinants
    ALLOCATE(DoublesHist(nBasis,nBasis,nBasis,nBasis),stat=ierr)
    ALLOCATE(AllDoublesHist(nBasis,nBasis,nBasis,nBasis),stat=ierr)
    IF(ierr.ne.0) THEN
        CALL Stop_All("TestGenRandSymExcitNU","Not possible to allocate memory to do histogramming")
    ENDIF
    ALLOCATE(SinglesHist(nBasis,nBasis),stat=ierr)
    ALLOCATE(AllSinglesHist(nBasis,nBasis),stat=ierr)
    IF(ierr.ne.0) THEN
        CALL Stop_All("TestGenRandSymExcitNU","Not possible to allocate memory to do histogramming")
    ENDIF
    DoublesHist(:,:,:,:)=0.D0
    SinglesHist(:,:)=0.D0
    AllDoublesHist(:,:,:,:)=0.D0
    AllSinglesHist(:,:)=0.D0

    CALL EncodeBitDet(nI,iLut)

    tFilled=.false.
    Scratch1(:)=0
    Scratch2(:)=0

    AverageContrib=0.D0
    AllAverageContrib=0.D0
    ForbiddenIter=0
!    pDoub=1.D0
!    IF(iProcIndex.eq.0) OPEN(9,FILE="AvContrib",STATUS="UNKNOWN")

    iter_tmp = iter
    do i=1,Iterations
        iter = i

    
        IF(mod(i,40000).eq.0) THEN
            WRITE(6,"(A,I10)") "Iteration: ",i
            CALL FLUSH(6)
        ENDIF

        call gen_rand_excit3 (nI, iLut, nJ, iLutnJ, exFlag, IC, ExcitMat, &
                             tParity, pGen, HElGen, tFilled, Scratch1, Scratch2, &
                             Scratch3)
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
        AverageContrib=AverageContrib+1.D0/pGen

!        CALL EncodeBitDet(nJ,iLutnJ)
!        IF(IC.eq.1) THEN
!            WRITE(6,*) ExcitMat(1,1),ExcitMat(2,1)
!        ELSE
!            WRITE(6,*) "Double Created"
!            WRITE(6,*) ExcitMat(1,1),ExcitMat(1,2),ExcitMat(2,1),ExcitMat(2,2)
!        ENDIF

        IF(IC.eq.1) THEN
            SinglesHist(ExcitMat(1,1),ExcitMat(2,1))=SinglesHist(ExcitMat(1,1),ExcitMat(2,1))+(1.D0/pGen)
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
            DoublesHist(Ind1,Ind2,Ind3,Ind4)=DoublesHist(Ind1,Ind2,Ind3,Ind4)+(1.D0/pGen)
        ENDIF
!        IF(mod(i,iWriteEvery).eq.0) THEN
!            AllAverageContrib=0.D0
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
        CALL IsSymAllowedExcit(nI,nJ,IC,ExcitMat)

    enddo
    iter = iter_tmp

!    IF(iProcIndex.eq.0) CLOSE(9)

#ifdef PARALLEL
    CALL MPI_BARRIER(MPI_COMM_WORLD,error)
    CALL MPI_AllReduce(DoublesHist,AllDoublesHist,nBasis**4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
    CALL MPI_AllReduce(SinglesHist,AllSinglesHist,nBasis**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
#else
    AllDoublesHist=DoublesHist
    AllSinglesHist=SinglesHist
#endif

!Now run through arrays normalising them so that numbers are more managable.
    IF(iProcIndex.eq.0) THEN
        OPEN(8,FILE="DoublesHist3",STATUS="UNKNOWN")
        DetNum=0
        do i=1,nBasis-1
            do j=i+1,nBasis
                do k=1,nBasis-1
                    do l=k+1,nBasis
                        IF(AllDoublesHist(i,j,k,l).gt.0.D0) THEN
    !                        DoublesHist(i,j,k,l)=DoublesHist(i,j,k,l)/real(Iterations,8)
                            DetNum=DetNum+1
                            ExcitMat(1,1)=i
                            ExcitMat(1,2)=j
                            ExcitMat(2,1)=k
                            ExcitMat(2,2)=l
                            CALL FindExcitBitDet(iLut,iLutnJ,2,ExcitMat)
                            WRITE(8,"(I12,F20.12,2I5,A,2I5,I15)") DetNum,AllDoublesHist(i,j,k,l)/(real(Iterations,8)*nProcessors),i,j,"->",k,l,iLutnJ(0)
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
                IF(AllSinglesHist(i,j).gt.0.D0) THEN
                    DetNumS=DetNumS+1
                    ExcitMat(1,1)=i
                    ExcitMat(2,1)=j
                    CALL FindExcitBitDet(iLut,iLutnJ,1,ExcitMat)
                    WRITE(9,*) DetNumS,AllSinglesHist(i,j)/(real(Iterations,8)*nProcessors),i,"->",j
!                    WRITE(6,*) DetNumS,AllSinglesHist(i,j),i,"->",j
                ENDIF
            enddo
        enddo
        CLOSE(9)
        WRITE(6,*) DetNumS," Single excitations found from nI"
        IF((DetNum+DetNumS).ne.ExcitCount) THEN
            CALL construct_class_counts(nI,ClassCount2,ClassCountUnocc2)
            WRITE(6,*) "Total determinants = ", ExcitCount
            WRITE(6,*) "ClassCount2(:)= ",ClassCount2(:)
            WRITE(6,*) "***"
            WRITE(6,*) "ClassCountUnocc2(:)= ",ClassCountUnocc2(:)
            CALL Stop_All("TestGenRandSymExcitNU","Not all excitations accounted for...")
        ENDIF
    ENDIF
    CALL MPIBarrier(error)

    END SUBROUTINE

end module

