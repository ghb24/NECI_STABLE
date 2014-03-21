#include "macros.h"

MODULE ISKRandExcit 
!ISK (Inversion-Symmetry K-point) wavefunctions are a linear combinations of two determinants, 
!where all orbitals are exchanged with 
!the equivalent orbital where the k-points of all occupied orbitals have been switched.
!In simple notation we will consider an excitation where determinants i and j are in the original ISK,
!and determinants a and b in the excited ISK.
!As with HPHF wavefunctions, a large simplification occurs when it is realised that P(i->a) = P(j->b),
!therefore all excitation are connected to the determinantal excitations of just one of the constituent determinants.
!This means that only one constituent determinant will be considered in the space.

    use SystemData, only: nel, tCSF
    use GenRandSymExcitNUMod, only: gen_rand_excit, ScratchSize 
    use HPHFRandExcitMod, only: CalcNonUniPGen
    use DetBitOps, only: DetBitLT, DetBitEQ, FindExcitBitDet, &
                         FindBitExcitLevel, TestClosedShellDet
    use FciMCData, only: pDoubles, excit_gen_store_type
    use constants, only: dp,n_int,bits_n_int
    use HElem
    use bit_reps, only: NIfD, NIfDBO, NIfTot
    use SymExcitDataMod, only: SpinOrbSymLabel,SymTableLabels,SymInvLabel, &
                               KPntInvSymOrb
    IMPLICIT NONE

    contains

    subroutine gen_ISK_excit (nI, iLutnI, nJ, iLutnJ, exFlag, IC, ExcitMat, &
                               tParity, pGen, HEl, store)
        use FciMCData, only: tGenMatHEl
        integer, intent(in) :: nI(nel) 
        integer(kind=n_int), intent(in) :: iLutnI(0:niftot)
        integer, intent(in) :: exFlag
        integer, intent(out) :: nJ(nel)
        integer(kind=n_int), intent(out) :: iLutnJ(0:niftot)
        integer, intent(out) :: IC, ExcitMat(2,2)
        logical, intent(out) :: tParity ! Not used
        type(excit_gen_store_type), intent(inout), target :: store
        real(dp), intent(out) :: pGen
        HElement_t, intent(out) :: HEl
        character(*), parameter :: this_routine='gen_ISK_excit'
        logical :: tSignOrig,tSwapped,tSame_ISK,tCrossConnected,tSignCross
        integer(n_int) :: iLutnJSym(0:NIfTot)
        integer :: nJSym(NEl),CrossIC,CrossEx(2,2)
        real(dp) :: pGen2

        !First, generate a random excitation from the determinant which is given in the argument
        call gen_rand_excit (nI, iLutnI, nJ, iLutnJ, exFlag, IC, ExcitMat, &
                             tSignOrig, pGen, HEl, store)
        if(IsNullDet(nJ)) return    !No excitation created

!Create bit representation of excitation - iLutnJ
        call FindExcitBitDet(iLutnI,iLutnJ,IC,ExcitMat)

        if(is_self_inverse(nJ,iLutnJ)) then
!There is only one way which we could have generated the excitation nJ since it has no determinant partner.
!Also, we will always return the 'correct' determinant.
            if(tGenMatHEl) then
                !Generate matrix element of excitation.
                call stop_all(this_routine,"tGenMatHEl not yet implemented")
            endif
        else
!This excitation *does* have a inverse determinant. Could we have generated this instead?
!Find the inverse
            call returned_invsym(nJ,nJSym,iLutnJ,iLutnJSym,.true.,tSwapped)

!Calculate whether the 'cross-term' is non-zero. i.e. is the original determinant connected to the 
!inverse determinant of the excitation.
            if(tSwapped) then
                !We have swapped the excitation, so that iLutnJ now contains the cross determinant
                call ISK_cross_det_conn(iLutnI,iLutnJ,tSame_ISK,tCrossConnected,CrossIC,CrossEx,tSignCross)
            else
                call ISK_cross_det_conn(iLutnI,iLutnJSym,tSame_ISK,tCrossConnected,CrossIC,CrossEx,tSignCross)
            endif

            if(tSame_ISK) then
                !We have created the same ISK - return null ISK
                nJ(1)=0
                if(tGenMatHEl) HEl=0.0_dp
            elseif(tCrossConnected) then
                !The cross-term is connected. Calculate the probability that we created this det instead (in same ISK)
                CALL CalcNonUniPGen(nI, ilutnI, CrossEx, CrossIC, &
                                    store%ClassCountOcc, &
                                    store%ClassCountUnocc,pDoubles,pGen2)
                pGen=pGen+pGen2

                if(tGenMatHEl) then
                    !calculate matrix element to open shell excitation.
                endif
            else
                !The cross-term is not connected. Therefore, there was only one way we could have generated this ISK.
                !nI *CANNOT* be a self-inverse ISK here, otherwise we should have had a cross-term.
                !Check this, but remove check once we know it is working.
                if(is_self_inverse(nI,iLutnI)) then
                    write(6,*) "no cross term, therefore should not be self-inv, but is..."
                    write(6,*) "nI: ",nI(:)
                    write(6,*) "nJ: ",nJ(:)
                    write(6,*) "nJSym: ",nJSym(:)
                    write(6,*) tSwapped
                    call stop_all(this_routine,"Dodgy logic - this should not be self-inv")
                endif
                
                if(tGenMatHEl) then
                endif
            endif
        endif

        ! Eliminate compiler warnings
        tParity = tParity

    end subroutine gen_ISK_excit


    !determine whether two determinants are connected.
    !Return whether the ISKs are actually the same
    !If the cross term is found to exist, return the excitation matrix and parity of the excitation
    !TODO: Optimisation - do we need to calculate the excitation matrix, or can we just pass it in?
    subroutine ISK_cross_det_conn(iLutnI,iLutnJSym,tSame_ISK,tcross_conn,CrossIC,ExCross,tExSign)
        integer(n_int), intent(in) :: iLutnI(0:NIfTot),iLutnJSym(0:NIfTot)
        logical, intent(out) :: tcross_conn,tSame_ISK,tExSign
        integer, intent(out) :: ExCross(2,2),CrossIC
        integer :: SymB,ijSymProd

        tSame_ISK=.false.   !Is the ISK generated actually the one we started with?

        !Comment out these lines when we are sure above call is working!
        CrossIC = FindBitExcitLevel(iLutnI, iLutnJSym, 2)

        if(CrossIC.eq.0) then
            !We have generated the same ISK!!
            !The inverse of the excitation is the determinant that we started with!
            !This does not want to be allowed.
            tSame_ISK=.true.
!            write(6,*) "Generated same det"
            return
        elseif(CrossIC.gt.2) then
            tcross_conn=.false.
!            write(6,*) "Det is more than double excit",CrossIC
            return
        endif

!        write(6,*) "Det is allowed by excitlevel: ",CrossIC

        !cross-term is allowed by excitation level. Is it allowed by momentum conservation?
        !calculate excitation matrix
        ExCross(1,1)=CrossIC
        call GetBitExcitation(iLutnI,iLutnJSym,ExCross,tExSign)
        !i = ex(1,1)
        !j = ex(1,2)
        if(CrossIC.eq.1) then
            if(SpinOrbSymLabel(ExCross(1,1)).ne.SpinOrbSymLabel(ExCross(2,1))) then
                !symmetry forbidden
                tcross_conn=.false.
            else
                tcross_conn=.true.
            endif
        else
            !Excitlevel *must* equal 2
            ijSymProd=SymTableLabels(SpinOrbSymLabel(ExCross(1,1)),SpinOrbSymLabel(ExCross(1,2)))
            SymB=SymTableLabels(ijSymProd,SymInvLabel(SpinOrbSymLabel(ExCross(2,1))))
            if(SymB.ne.SpinOrbSymLabel(ExCross(2,2))) then
                !symmetry forbidden
                tcross_conn=.false.
            else
                tcross_conn=.true.
            endif
        endif

    end subroutine ISK_cross_det_conn

!Function to determine whether a determinant is its own self inverse when inverted.
!If it is (returns true), then there is no other determinant in the ISK function.
!Input: the determinant to test, in both natural ordered, and bit representations.
!TODO: This could be sped up by a factor of two, since we only actually need to 
!search through half of the electrons, since if it is symmetric, then the other 
!half should be inverses of ones that we've already tested!
    pure function is_self_inverse(nI,iLut) result(tSelfInv)
        integer, intent(in) :: nI(Nel)
        integer(n_int), intent(in) :: iLut(0:niftot)
        logical :: tSelfInv
        integer :: pos,i

        tSelfInv=.true.

        do i=1,NEl
            !run through the electrons, find inverse and test whether it is in the bit representation of the determinant.
            pos = KPntInvSymOrb(nI(i))-1

            if(.not. (btest(iLut(pos/bits_n_int),mod(pos,bits_n_int)))) then
                !the inverse orbital is not found in the determinant.
                !therefore the determinant can not be a self inverse.
                tSelfInv=.false.
                return
            endif
        enddo

    end function is_self_inverse

!This routine will take a determinant (nI and nISym) and return in the same place the correct
!determinant out of the two inverse determinants which make up the ISK. If this is actually its
!partner determinant, then tSwapped is returned as true, and nISym and iLutnISym are returned
!as the original determinant.
!The flag tCalcnISym let the routine know whether to calculate the symmetry
!partner of nI initially, or whether it is passed in correctly.
!Whether one determinant or its partner is returned as nI is purely based on which bit representation is 'largest'.
    subroutine returned_invsym(nI,nISym,iLutnI,iLutnISym,tCalcnISym,tSwapped)
        integer, intent(inout) :: nI(NEl),nISym(NEl)
        integer :: nTemp(NEl),i
        integer(n_int) , intent(inout) :: iLutnI(0:NIfTot),iLutnISym(0:NIfTot)
        integer(n_int) :: iLutTemp(0:NIfTot)
        logical, intent(in) :: tCalcnISym
        logical, intent(out) :: tSwapped

        if(tCalcnISym) then
            call find_invsym_det(nI,nISym,iLutnISym)
        endif

!        write(6,*) "Original det: ",nI(:)
!        write(6,*) "Inv det: ",nISym(:)

        ! iLutnI is 'less' than iLutSym, so iLutSym is the determinant with 
        ! the first open-shell = alpha. Swap them around.
        i=DetBitLT(iLutnI,iLutnISym,NIfD)
        if(i.eq.1) then
            iLutTemp(:)=iLutnI(:)
            iLutnI(:)=iLutnISym(:)
            iLutnISym(:)=iLutTemp(:)
            nTemp(:)=nI(:)
            nI(:)=nISym(:)
            nISym(:)=nTemp(:)
            tSwapped=.true.
        elseif(i.eq.0) then
            CALL Stop_All("Returned_InvSym","Shouldn't have self-inverse determinants in here")
        else
            tSwapped=.false.
        endif

    end subroutine returned_invsym

!When inverting a determinant, the ordering of the orbitals can change.
!By ordering, we can introduce a phase to the calculation, which we have to
!calculate. Here, we take a non-ordered determinant, and the same ordered
!determinant, and calculate the number of permutations required to line up the
!determinants.
!Other routines for finding permutations look at excitations between normal-ordered
!determinants which are excitations of each other, and are therefore no use here.
    subroutine find_invsym_permut(nINonOrder,nIOrder,tPermute)
        use util_mod, only : swap 
        integer, intent(in) :: nINonOrder(NEl),nIOrder(Nel)
        logical, intent(out) :: tPermute
        integer :: nITemp(NEl),permute,i,j,orb
        logical :: tFoundOrb

        !Initially, do this is a very inefficient way.
        nITemp=nINonOrder
        permute=0

        do i=1,NEl
            if(nITemp(i).eq.nIOrder(i)) then
                cycle
            endif

            !Need to find nIOrder(i) in nITemp and count permutations
            tFoundOrb=.false.
            do j=Nel,i+1,-1
                !Orbital must be higher
                !Search from the end of the list to find orbital we want.
                !Once found, cycle back to where we started, swapping orbitals as we go.
                !If optimisation needed, then we should be able to skip the rest of the loop
                !once the correct orbital has been found.
                if((.not.tFoundOrb).and.(nITemp(j).ne.nIOrder(i))) then
                    cycle
                else
                    !We have found the orbital / it is at position j.
                    tFoundOrb=.true.
                    !Move it to j-1
                    call swap(nITemp(j),nITemp(j-1))
                    permute=permute+1
                endif
            enddo

            !Now, the correct orbital should be at i.
            if(nIOrder(i).ne.nITemp(i)) then
                call stop_all("find_invsym_permut","Error in finding permutation between differently ordered determinants")
            endif
        enddo

        tpermute = mod(permute,2) == 1  !Odd or even number of permutations?

    end subroutine find_invsym_permut

!Find the inverse of a determinant *without* reordering the new determinant into natural ordering.
    subroutine find_invsym_det_noorder(nI,nISym,iLutnISym)
        integer , intent(in) :: nI(NEl)
        integer(n_int) , intent(out) :: iLutnISym(0:NIfTot)
        integer , intent(out) :: nISym(NEl)
        integer :: orb,pos,i
        iLutnISym(:)=0
        do i=1,NEl
            orb=KPntInvSymOrb(nI(i))
            nISym(i)=orb
            pos=(orb-1)/bits_n_int
            iLutnISym(pos)=ibset(iLutnISym(pos),mod(orb-1,bits_n_int))
        enddo
    end subroutine find_invsym_det_noorder

    subroutine find_invsym_det(nI,nISym,iLutnISym)
!        use systemdata, only: nBasis
        use sort_mod
        integer , intent(in) :: nI(NEl)
        integer(n_int) , intent(out) :: iLutnISym(0:NIfTot)
        integer , intent(out) :: nISym(NEl)
        integer :: orb,pos,i

        iLutnISym(:)=0

!        write(6,*) "KPntInvSymOrb: "
!        do i=1,nBasis
!            write(6,*) i,KPntInvSymOrb(i)
!        enddo
        do i=1,NEl
            orb=KPntInvSymOrb(nI(i))
            nISym(i)=orb
            pos=(orb-1)/bits_n_int
            iLutnISym(pos)=ibset(iLutnISym(pos),mod(orb-1,bits_n_int))
        enddo

        call sort(nISym)    !Need to sort the determinant, since even if it is a self-inverse, it will change natural ordering.

    end subroutine find_invsym_det

!This routine will take an ISK nI, and find Iterations number of excitations of it. It will then histogram these, 
!summing in 1/pGen for every occurance of
!the excitation. This means that all excitations should be 0 or 1 after enough iterations. It will then count the 
!excitations and compare the number to the
!number of excitations generated using the full enumeration excitation generation.
    SUBROUTINE TestGenRandISKExcit(nI,Iterations,pDoub)
        Use SystemData , only : NEl,nBasis,G1,nBasisMax,tUseBrillouin
        use DetBitOps, only: EncodeBitDet
        use bit_reps, only: decode_bit_det
        use GenRandSymExcitNuMod, only: scratchsize
        use FciMCData, only: tGenMatHEl
        use util_mod, only: get_free_unit
        use sort_mod
        use HPHFRandExcitMod, only: BinSearchListHPHF
        use Parallel_neci
        IMPLICIT NONE
        INTEGER :: nIX(NEl)
        INTEGER :: i,Iterations,nI(NEl),nJ(NEl),DetConn,nI2(NEl),nJ2(NEl),DetConn2
        INTEGER :: iUniqueHPHF,iUniqueBeta,PartInd,ierr,iExcit
        real(dp) :: pDoub, pGen
        logical :: Unique, Die, tSwapped
        INTEGER(KIND=n_int) :: iLutnI(0:NIfTot),iLutnJ(0:NIfTot),iLutnI2(0:NIfTot),iLutSym(0:NIfTot)
        INTEGER(KIND=n_int), ALLOCATABLE :: ConnsAlpha(:,:),ConnsBeta(:,:),UniqueHPHFList(:,:)
        INTEGER , ALLOCATABLE :: ExcitGen(:)
        real(dp) , ALLOCATABLE :: Weights(:),AllWeights(:)
        INTEGER :: iMaxExcit,nStore(6),nExcitMemLen(1),j,k,l, iunit
        integer :: icunused, exunused(2,2)
        logical :: tParityunused, tTmp
        type(excit_gen_store_type) :: store
        HElement_t :: HEl

        tUseBrillouin=.false.

        CALL EncodeBitDet(nI,iLutnI)
        !Returns both forms
!        call find_invsym_det(nI,nI2,iLutnI2)
!        CALL FindDetSpinSym(nI,nI2,NEl)
!        CALL EncodeBitDet(nI2,iLutnI2)
        IF(is_self_inverse(nI,iLutnI)) THEN
            call find_invsym_det(nI,nI2,iLutnI2)
            IF(.not.DetBitEQ(iLutnI,iLutnI2,NIfDBO)) THEN
                WRITE(6,*) "nI: ",nI(:)
                WRITE(6,*) ""
                WRITE(6,*) "nISym: ",nI2(:)
                WRITE(6,*) ""
                WRITE(6,*) "iLutnI: ",iLutnI(:)
                WRITE(6,*) "iLutnISym: ",iLutnI2(:)
                WRITE(6,*) "***"
                CALL Stop_All("TestGenRandISKExcit","Self-inverse determinant entered, but coupled dets different...")
            ENDIF
        ELSE
            call returned_invsym(nI,nI2,iLutnI,iLutnI2,.true.,tSwapped)
            if(tSwapped) then
                write(6,*) "Swapped nI with inverse determinant to maintain excitations from correct det."
            endif
        ENDIF
        WRITE(6,*) "nI: ",nI(:)
        WRITE(6,*) ""
        WRITE(6,*) "nISym: ",nI2(:)
        WRITE(6,*) ""
        WRITE(6,*) "iLutnI: ",iLutnI(:)
        WRITE(6,*) "iLutnISym: ",iLutnI2(:)
        WRITE(6,*) "***"
        WRITE(6,*) Iterations,pDoub
!        WRITE(6,*) "nSymLabels: ",nSymLabels
        CALL neci_flush(6)

!First, we need to enumerate all possible ISK wavefunctions from each inverse-pair of determinants.
!These need to be stored in an array
!Find the number of symmetry allowed excitations there should be by looking at the full excitation generator.
!Setup excit generators for this determinant
        write(6,*) "Generating ALL excitations from first determinant..."
        iMaxExcit=0
        nStore(1:6)=0
        CALL GenSymExcitIt2(nI,NEl,G1,nBasis,.TRUE.,nExcitMemLen,nJ,iMaxExcit,nStore,3)
        ALLOCATE(EXCITGEN(nExcitMemLen(1)),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All("SetupExcitGen","Problem allocating excitation generator")
        EXCITGEN(:)=0
        CALL GenSymExcitIt2(nI,NEl,G1,nBasis,.TRUE.,EXCITGEN,nJ,iMaxExcit,nStore,3)
        CALL GetSymExcitCount(EXCITGEN,DetConn)
        WRITE(6,*) "Stored determinant (dubbed alpha determinant) has ",DetConn," total excitations:"
        ALLOCATE(ConnsAlpha(0:NIfTot,DetConn))
        i=1
        lp2: do while(.true.)
            CALL GenSymExcitIt2(nI,NEl,G1,nBasis,.false.,EXCITGEN,nJ,iExcit,nStore,3)
            IF(IsNullDet(nJ)) exit lp2
            if((iExcit.eq.1).and.(nJ(NEl).eq.105)) then
                write(6,*) "Generated single: "
                write(6,*) nJ(:)
            endif
            CALL EncodeBitDet(nJ,iLutnJ)
            IF(.not.is_self_inverse(nJ,iLutnJ)) THEN
!                CALL ReturnAlphaOpenDet(nJ,nJ2,iLutnJ,iLutSym,.true.,.true.,tSwapped)
                call returned_invsym(nJ,nJ2,iLutnJ,iLutSym,.true.,tSwapped)
            ENDIF
!            WRITE(6,"(4I4,A,I4,A,I13)") nJ(:), " *** ", iExcit, " *** ", iLutnJ(:)
            ConnsAlpha(0:NIfD,i)=iLutnJ(:)
            i=i+1
        enddo lp2

!Now we also need to store the excitations from the other spin-coupled determinant.
        iMaxExcit=0
        nStore(1:6)=0
        DEALLOCATE(EXCITGEN)
        
        CALL GenSymExcitIt2(nI2,NEl,G1,nBasis,.TRUE.,nExcitMemLen,nJ,iMaxExcit,nStore,3)
        ALLOCATE(EXCITGEN(nExcitMemLen(1)),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All("SetupExcitGen","Problem allocating excitation generator")
        EXCITGEN(:)=0
        CALL GenSymExcitIt2(nI2,NEl,G1,nBasis,.TRUE.,EXCITGEN,nJ,iMaxExcit,nStore,3)
        CALL GetSymExcitCount(EXCITGEN,DetConn2)
        WRITE(6,*) "Inverse determinant has ",DetConn2," total excitations"
        ALLOCATE(ConnsBeta(0:NIfTot,DetConn2))
        i=1
        lp: do while(.true.)
            CALL GenSymExcitIt2(nI2,NEl,G1,nBasis,.false.,EXCITGEN,nJ,iExcit,nStore,3)
            IF(IsNullDet(nJ)) exit lp
            CALL EncodeBitDet(nJ,iLutnJ)
            IF(.not.is_self_inverse(nJ,iLutnJ)) THEN
!                CALL ReturnAlphaOpenDet(nJ,nJ2,iLutnJ,iLutSym,.true.,.true.,tSwapped)
                call returned_invsym(nJ,nJ2,iLutnJ,iLutSym,.true.,tSwapped)
            ENDIF
!            IF(.not.TestClosedShellDet(iLutnJ)) THEN
!                CALL ReturnAlphaOpenDet(nJ,nJ2,iLutnJ,iLutSym,.true.,.true.,tSwapped)
!            ENDIF
!            WRITE(6,"(4I4,A,I4,A,I13)") nJ(:), " *** ",iExcit," *** ",iLutnJ(:)
            ConnsBeta(:,i)=iLutnJ(:)
            i=i+1
        enddo lp
        DEALLOCATE(EXCITGEN)

!Now we need to find how many HPHF functions there are.
        iUniqueHPHF=0
        do j=1,DetConn
!Run though all HPHF in the first array
            Unique=.true.
            do k=j-1,1,-1
!Run backwards through the array to see if this HPHF has come before
                IF(DetBitEQ(ConnsAlpha(0:NIfTot,k),ConnsAlpha(0:NIfTot,j),NIfDBO)) THEN
!This HPHF has already been counted before...
                    Unique=.false.
                    EXIT
                ENDIF
            enddo
            IF(Unique) THEN
!Unique HPHF found, count it
                iUniqueHPHF=iUniqueHPHF+1
            ENDIF
        enddo

        iUniqueBeta=0

!Now look through all the excitations for the spin-coupled determinant from the original HPHF...
        do j=1,DetConn2
!Run though all excitations in the first array, *and* up to where we are in the second array
            Unique=.true.
            do k=1,DetConn
                IF(DetBitEQ(ConnsAlpha(:,k),ConnsBeta(:,j),NIfDBO)) THEN
                    Unique=.false.
                    EXIT
                ENDIF
            enddo
            IF(Unique) THEN
!Need to search backwards through the entries we've already looked at in this array...
                do k=j-1,1,-1
                    IF(DetBitEQ(ConnsBeta(0:NIfTot,k),ConnsBeta(0:NIfTot,j),NIfDBO)) THEN
                        Unique=.false.
                        EXIT
                    ENDIF
                enddo
            ENDIF
            IF(Unique) THEN
                iUniqueHPHF=iUniqueHPHF+1
                iUniqueBeta=iUniqueBeta+1
            ENDIF
        enddo

        WRITE(6,*) "There are ",iUniqueHPHF," unique ISK wavefunctions from the ISK given."
        
        WRITE(6,*) "There are ",iUniqueBeta," unique ISK wavefunctions from the inverted determinant, " &
        & //"which are not in the alpha version."
        IF(iUniqueBeta.ne.0) THEN
            WRITE(6,*) "ISK from beta, but not from alpha!"
            CALL neci_flush(6)
            STOP
        ENDIF

        ALLOCATE(UniqueHPHFList(0:NIfTot,iUniqueHPHF))
        UniqueHPHFList(:,:)=0
!Now fill the list of HPHF Excitations.
        iUniqueHPHF=0
        do j=1,DetConn
!Run though all HPHF in the first array
            Unique=.true.
            do k=j-1,1,-1
!Run backwards through the array to see if this HPHF has come before
                IF(DetBitEQ(ConnsAlpha(0:NIfTot,k),ConnsAlpha(0:NIfTot,j),NIfDBO)) THEN
!This HPHF has already been counted before...
                    Unique=.false.
                    EXIT
                ENDIF
            enddo
            IF(Unique) THEN
!Unique HPHF found, count it
                iUniqueHPHF=iUniqueHPHF+1
                UniqueHPHFList(:,iUniqueHPHF)=ConnsAlpha(:,j)
            ENDIF
        enddo

!Now look through all the excitations for the spin-coupled determinant from the original HPHF...
        do j=1,DetConn2
!Run though all excitations in the first array, *and* up to where we are in the second array
            Unique=.true.
            do k=1,DetConn
                IF(DetBitEQ(ConnsAlpha(:,k),ConnsBeta(:,j),NIfDBO)) THEN
                    Unique=.false.
                    EXIT
                ENDIF
            enddo
            IF(Unique) THEN
!Need to search backwards through the entries we've already looked at in this array...
                do k=j-1,1,-1
                    IF(DetBitEQ(ConnsBeta(:,k),ConnsBeta(:,j),NIfDBO)) THEN
                        Unique=.false.
                        EXIT
                    ENDIF
                enddo
            ENDIF
            IF(Unique) THEN
                iUniqueHPHF=iUniqueHPHF+1
                UniqueHPHFList(:,iUniqueHPHF)=ConnsBeta(0:NIfD,j)
            ENDIF
        enddo

!Now sort the list, so that it can be easily binary searched.
        ALLOCATE(ExcitGen(iUniqueHPHF))
        call sort (UniqueHPHFList(:,1:iUniqueHPHF), ExcitGen)
        DEALLOCATE(ExcitGen)

        WRITE(6,*) "Unique ISK wavefunctions are: ",iUniqueHPHF, " written to unit 78"
        do i=1,iUniqueHPHF
            WRITE(78,*) UniqueHPHFList(0:NIfTot,i)
        enddo

        ALLOCATE(Weights(iUniqueHPHF))
        ALLOCATE(AllWeights(iUniqueHPHF))
        AllWeights(:)=0.0_dp
        Weights(:)=0.0_dp
        store%tFilled = .false.

        write(6,*) "Generating ISK random excitations..."

        do i=1,Iterations

            IF((mod(i,10000).eq.0).or.(i.le.200)) WRITE(6,"(A,I10)") "Iteration: ",i

!            tTmp = tGenMatHEl
!            tGenMatHel = .false.
            call gen_ISK_excit (nI, iLutnI, nJ, iLutnJ, 3, icunused, &
                                 exunused, tparityunused, pGen, HEl, store)
!            tGenMatHel = tTmp
!            CALL GenRandSymExcitNU(nI,iLut,nJ,pDoub,IC,ExcitMat,TParity,exFlag,pGen)

!Search through the list of HPHF wavefunctions to find slot.
            CALL BinSearchListHPHF(iLutnJ,UniqueHPHFList(:,1:iUniqueHPHF),iUniqueHPHF,1,iUniqueHPHF,PartInd,Unique)

            IF(.not.Unique) THEN
                write(6,*) "iLutnJ: ",iLutnJ(:)
                write(6,*) "nJ: ",nJ(:)
                CALL Stop_All("TestGenRandISKExcit","Cannot find excitation in list of allowed excitations")
            ENDIF

            Weights(PartInd)=Weights(PartInd)+(1.0_dp/pGen)
             
!Check excitation
!            CALL IsSymAllowedExcit(nI,nJ,IC,ExcitMat)

        enddo

        call MPISumAll(Weights,AllWeights)

        if(iProcIndex.eq.Root) then
        
            iunit = get_free_unit()
            OPEN(iunit,FILE="PGenHist",STATUS="UNKNOWN")

    !normalise excitation probabilities
            Die=.false.
            do i=1,iUniqueHPHF
                AllWeights(i)=AllWeights(i)/(real(Iterations,dp)*real(nNodes,dp))
                IF(abs(AllWeights(i)-1.0_dp).gt.0.1_dp) THEN
                    WRITE(6,*) "Error here!",i
                    Die=.true.
                ENDIF
    !            WRITE(6,*) i,UniqueHPHFList(0:NIfTot,i),Weights(i)
                call decode_bit_det (nIX, UniqueHPHFList(0:NIfTot,i))
    !            WRITE(6,*) nIX(:)
                WRITE(iunit,"(I8,G25.10)",advance='no') i,AllWeights(i)
                do j=0,NIfTot-1
                    WRITE(iunit,"(I24)",advance='no') UniqueHPHFList(j,i)
                enddo
                WRITE(iunit,"(I24)") UniqueHPHFList(NIfTot,i)
                write(79,*) i,"***"
                do j=1,NEl-1
                    write(79,"(I5)",advance='no') nIX(j)
                enddo
                write(79,"(I5)") nIX(NEl)
            enddo

            CLOSE(iunit)
            IF(Die) THEN
                CALL Stop_All("IUB","TestFail")
            ENDIF
        endif
        call mpiBarrier(ierr)

    END SUBROUTINE TestGenRandISKExcit


END MODULE ISKRandExcit 
