#include "macros.h"
module MomInvRandExcit
    use constants, only: dp,n_int
    use SystemData, only: NEl, G1, nBasis, tAntisym_MI, modk_offdiag
    use DetBitOps, only: DetBitEQ, FindExcitBitDet, FindBitExcitLevel
    use sltcnd_mod, only: sltcnd, sltcnd_excit
    use bit_reps, only: NIfD, NIfTot, NIfDBO, decode_bit_det
    use MomInv, only: IsBitMomSelfInv, InvertMomDet, ReturnMomAllowedDet, InvertMomBitDet 
    use MomInv, only: IsMomSelfInv, CalcMomAllowedBitDet
    use SymExcitDataMod, only:  excit_gen_store_type
    use FciMCData, only: tGenMatHEl,pDoubles
    use GenRandSymExcitNUMod, only: gen_rand_excit, CalcNonUniPGen
    implicit none

    contains
    
    
    subroutine gen_MI_excit (nI, iLutnI, nJ, iLutnJ, exFlag, IC, ExcitMat, &
                               parity_unused, pGen, HEl, store)

        ! Generate an MI excitation using only one of the determinants in 
        ! the source MI function.
        !
        ! --> Relies on both determinants in the MI function being connected
        !     to all excited MI functions.
        ! --> nI will always need to be a unique choice of determinant within
        !     each MI function, and then we nevver need to (explicitly) 
        !     refer to its spin-coupled partner.
        ! --> If tGenMatEl is true, the Hamiltonian matrix element between the
        !     two determinants will be calculated, and returned in Hel

        integer, intent(in) :: nI(nel) 
        integer(kind=n_int), intent(in) :: iLutnI(0:niftot)
        integer, intent(in) :: exFlag
        integer, intent(out) :: nJ(nel)
        integer(kind=n_int), intent(out) :: iLutnJ(0:niftot)
        integer, intent(out) :: IC, ExcitMat(2,2)
        integer, intent(out) :: parity_unused ! Not used
        real(dp), intent(out) :: pGen
        HElement_t, intent(out) :: HEl
        type(excit_gen_store_type), intent(inout), target :: store

        integer(kind=n_int) :: iLutnJ2(0:niftot)
        integer :: nJ2(nel), ex2(2,2), excitLevel, parity, parity_orig
        real(dp) :: pGen2
        HElement_t :: MatEl, MatEl2
        logical :: TestClosedShellDet
        logical :: tSwapped
        character(len=*) , parameter :: t_r='gen_MI_excit'

        ! Avoid warnings
        parity_unused = 1

        call gen_rand_excit (nI, iLutnI, nJ, iLutnJ, exFlag, IC, ExcitMat, &
                             parity_orig, pGen, HEl, store)

!Create excitation of uniquely chosen determinant in this HPHF function.
        IF(IsNullDet(nJ)) RETURN
!Create bit representation of excitation - iLutnJ
        CALL FindExcitBitDet(iLutnI,iLutnJ,IC,ExcitMat)
            
        IF(IsMomSelfInv(nJ,iLutnJ)) THEN
!There is only one way which we could have generated the excitation nJ since it has no momentum-partner. 
!Also, we will always return the 'correct' version.
            IF(tGenMatHEl) THEN
!Generate matrix element -> MI function to closed shell det.
                IF(IsMomSelfInv(nI,iLutnI)) THEN
                    !Self-inv -> Self-inv
                    HEl = sltcnd_excit (nI, IC, ExcitMat, parity_orig)
                ELSE
                    !Mom coupled -> Self-inv
                    if(tAntisym_MI) then
                        HEl=0.0_dp
                    else
                        MatEl = sltcnd_excit (nI, IC, ExcitMat, parity_orig)
                        HEl=MatEl*SQRT(2.0_dp)
                    endif
                ENDIF
                if (ic /= 0 .and. modk_offdiag) hel = -abs(hel)
            ENDIF
        ELSE
!Momentum paired excitation - could we have generated the momentum inverted determinant instead?

!Find the momentum inverted version. 
            CALL CalcMomAllowedBitDet(nJ,nJ2,iLutnJ,iLutnJ2,.true.,.true.,tSwapped)

!Try and find if spin-coupled determinant from excitation is attached.
            IF(tSwapped) THEN
                ExcitLevel = FindBitExcitLevel(iLutnI, iLutnJ, 2)
            ELSE
                ExcitLevel = FindBitExcitLevel(iLutnI, iLutnJ2, 2)
            ENDIF

            IF((ExcitLevel.eq.2).or.(ExcitLevel.eq.1)) THEN     
                !Momentum inverted determinant is connected to root det.
                !This is the connection to the determinant we didn't generate...

                Ex2(1,1)=ExcitLevel

                IF(tSwapped) THEN
                    CALL GetBitExcitation(iLutnI,iLutnJ,Ex2,parity)
                ELSE
                    CALL GetBitExcitation(iLutnI,iLutnJ2,Ex2,parity)
                ENDIF
                !Calc probability of generating it
                CALL CalcNonUniPGen(nI, Ex2, ExcitLevel, store%ClassCountOcc,&
                                    store%ClassCountUnocc, pDoubles, pGen2)
!!We cannot guarentee that the pGens are going to be the same - in fact, generally, they wont be.
                pGen=pGen+pGen2

                IF(tGenMatHEl) THEN
!Generate matrix element to MI excitation
                    IF(IsMomSelfInv(nI,iLutnI)) THEN    
                        !Self-inv det -> Momentum paired func : Want to sum in SQRT(2)* Hij
                        if(tAntisym_MI) then
                            !Must be zero matrix element
                            HEl=0.0_dp
                        else
                            !Could actually use either of these - they should be the same?!
                            IF(tSwapped) THEN
                                MatEl = sltcnd_excit (nI, IC, Ex2, parity)
                            ELSE
                                MatEl = sltcnd_excit (nI, IC, ExcitMat, &
                                                      parity_orig)
                            ENDIF
                            HEl=MatEl*SQRT(2.0_dp)
                        endif
                    ELSE     !MI function -> MI function

!First find nI -> nJ. If nJ has swapped, then this will be different.
                        IF(tSwapped) THEN
                            MatEl = sltcnd_excit (nI, ExcitLevel, Ex2, &
                                                  parity)
                        ELSE
                            MatEl = sltcnd_excit (nI, IC, ExcitMat, &
                                                  parity_orig)
                        ENDIF

                        !now nI2 -> nJ - cross term present
                        IF((ExcitLevel.eq.2).or.(ExcitLevel.eq.1)) THEN
                            
!                            CALL CalcOpenOrbs(iLutnJ,OpenOrbsJ)
!                            CALL CalcOpenOrbs(iLutnI,OpenOrbsI)

                            IF(tSwapped) THEN
!                                IF((OpenOrbsJ+OpenOrbsI).eq.3) tSignOrig=.not.tSignOrig  !I.e. J odd and I even or vice versa, but since these can only be at max quads, then they can only have 1/2 open orbs
                                MatEl2 = sltcnd_excit (nI, IC, ExcitMat, &
                                                       parity_orig)
                            ELSE
!                                IF((OpenOrbsJ+OpenOrbsI).eq.3) tSign=.not.tSign     !I.e. J odd and I even or vice versa, but since these can only be at max quads, then they can only have 1/2 open orbs
                                MatEl2 = sltcnd_excit (nI,  ExcitLevel, &
                                                       Ex2, parity)
                            ENDIF

                            if(tAntisym_MI) then
                                MatEl=MatEl-MatEl2
                            else
                                MatEl=MatEl+MatEl2
                            endif

                        ENDIF
                        HEl=MatEl
                    
                    ENDIF   !Endif from MomInv/not det

                    if (ic /= 0 .and. modk_offdiag) hel = -abs(hel)
                ENDIF   !Endif want to generate matrix element

            ELSEIF(ExcitLevel.eq.0) THEN
!We have generated the same MI function. MatEl wants to be zero.
                nJ(1)=0
                IF(tGenMatHEl) THEN
                    HEl=0.0_dp
                ENDIF

            ELSE    !MI func to MI func, but with no cross-connection.
                
                IF(tGenMatHEl) THEN
!iLutnI MUST be an MI func here, since otherwise it would have been connected to iLutnJ2. Also, we know the cross connection (i.e. MatEl2 = 0)
                    IF(tSwapped) THEN
                        if(tAntisym_MI) parity_orig = -parity_orig
                        MatEl = sltcnd_excit(nI,  IC, ExcitMat, parity_orig)
                    ELSE
                        MatEl = sltcnd_excit (nI, IC, ExcitMat, parity_orig)
                    ENDIF

                    HEl=MatEl
                    if (ic /= 0 .and. modk_offdiag) hel = -abs(hel)
                        
                ENDIF


            ENDIF
            
        ENDIF

    end subroutine gen_MI_excit


!    !This routine checks the properties of the MI functions, so that we know we can use the more efficient
!    !HPHF-type excitation generation and integral evaluation routines.
!    !This relies on the fact that all connected MI functions can be generated from the determinant excitations
!    !of just one determinant in an MI function.
!    subroutine TestMomInvInts()
!        use util_mod, only: get_free_unit
!        use Determinants, only: get_helement
!        use symexcitdatamod, only: ScratchSize
!        use GenRandSymExcitNUMod, only: construct_class_counts,CalcNonUniPGen
!        implicit none
!        character(len=*), parameter :: t_r='TestMomInvInts'
!        integer :: PairedDetsUnit,nPairedDets,i,TempnI(NEl),TempnI2(NEl),IC,ICSym,ICConnect,SelfInvDetsUnit,nSelfDets
!        real(dp) :: State1,State2,pGen_ia,pGen_ib,pGen_ja,pGen_jb
!        integer(n_int) :: iLut(0:NIfTot),iLut2(0:NIfTot)
!        integer(n_int), allocatable :: PairedDets(:,:),SelfInvDets(:,:)
!        logical :: tSwapped,tSign
!        integer :: j,IC_ia,IC_ib,IC_ja,IC_jb,Ex_ia(2,2),Ex_ib(2,2),Ex_ja(2,2),Ex_jb(2,2)
!        integer(n_int) :: iLutnI(0:NIfTot),iLutnJ(0:NIfTot),iLutnI2(0:NIfTot),iLutnJ2(0:NIfTot)
!        HElement_t :: Hel_ia,Hel_ib,Hel_ja,Hel_jb
!        integer :: nI(NEl),nJ(NEl),nI2(NEl),nJ2(NEl)
!        integer :: CC_ia(ScratchSize),CC_ib(ScratchSize),CC_ja(ScratchSize),CC_jb(ScratchSize)  
!        integer :: CCU_ia(ScratchSize),CCU_ib(ScratchSize),CCU_ja(ScratchSize),CCU_jb(ScratchSize)  
!
!        write(6,*) "Checking properties of MI functions..."
!
!        !Open files with determinants written out from diagonalisation
!        PairedDetsUnit=get_free_unit()
!
!        open(PairedDetsUnit,file='LzPairedDets',status='old',action='read')
!        nPairedDets=0
!        iLut=0
!        iLut2=0
!        do while(.true.)
!            read(PairedDetsUnit,"(2I16,15I4,2G19.8)",end=99) iLut(0),iLut2(0),TempnI(1:NEl),TempnI2(1:NEl),IC,ICSym,ICConnect,State1,State2
!            nPairedDets=nPairedDets+1
!        enddo
!99      continue
!        rewind(PairedDetsUnit)
!        allocate(PairedDets(0:NIfTot,nPairedDets))
!        PairedDets=0
!        i=0
!        do while(.true.)
!            read(PairedDetsUnit,"(2I16,15I4,2G19.8)",end=199) iLut(0),iLut2(0),TempnI(1:NEl),TempnI2(1:NEl),IC,ICSym,ICConnect,State1,State2
!            i=i+1
!            call ReturnMomAllowedDet(TempnI,TempnI2,iLut,iLut2,tSwapped)
!!            write(6,*) tSwapped
!            PairedDets(:,i)=iLut(:) !Store the allowed one only
!        enddo
!199     continue
!        close(PairedDetsUnit)
!
!        !Now get self-inverse dets
!        SelfInvDetsUnit=get_free_unit()
!        open(SelfInvDetsUnit,file='SelfInvDet',status='old',action='read')
!        nSelfDets=0
!        do while(.true.)
!            read(SelfInvDetsUnit,"(I16,6I4,G19.8)",end=98) iLut(0),TempnI(1:NEl),State1
!            nSelfDets=nSelfDets+1
!        enddo
!98      continue
!        rewind(SelfInvDetsUnit)
!        allocate(SelfInvDets(0:NIfTot,nSelfDets))
!        SelfInvDets=0
!        i=0
!        do while(.true.)
!            i=i+1
!            read(SelfInvDetsUnit,"(I16,6I4,G19.8)",end=198) SelfInvDets(0,i),TempnI(1:NEl),State1
!        enddo
!198     continue
!        close(SelfInvDetsUnit)
!
!        !We now have two arrays with all "allowed" MomInvDets in them - SelfInvDets and PairedDets
!        !We first want to test the connections between *All* MI functions.
!        write(6,*) "Running through all pairs of MI functions..."
!        WRITE(6,*) "CHECKING CONNECTIVITY..."
!        write(6,*) "CHECKING MATRIX ELEMENTS..."
!        write(6,*) "CHECKING PGENS..."
!        do i=1,nPairedDets
!
!            do j=1,nPairedDets
!
!                !For this i,j pair, isolate all four hamiltonian matrix element between determinants.
!                iLutnI(:)=PairedDets(:,i)
!                iLutnJ(:)=PairedDets(:,j)
!
!                call InvertMomBitDet(iLutnI,iLutnI2)
!                call InvertMomBitDet(iLutnJ,iLutnJ2)
!
!                IC_ia = FindBitExcitLevel(iLutnI,iLutnJ,NEl)
!                IC_ib = FindBitExcitLevel(iLutnI,iLutnJ2,NEl)
!                IC_ja = FindBitExcitLevel(iLutnI2,iLutnJ,NEl)
!                IC_jb = FindBitExcitLevel(iLutnI2,iLutnJ2,NEl)
!
!                if(IC_ia.ne.IC_jb) call stop_all(t_r,"I,A excitation level not the same as J,B excitation")
!                if(IC_ib.ne.IC_ja) call stop_all(t_r,"I,B excitation level not the same as J,A excitation")
!
!                call decode_bit_det(nI,iLutnI)
!                call decode_bit_det(nI2,iLutnI2)
!                call decode_bit_det(nJ,iLutnJ)
!                call decode_bit_det(nJ2,iLutnJ2)
!
!                !Calculate matrix elements between them...
!                Hel_ia=get_helement(nI,nJ,IC_ia,iLutnI,iLutnJ)
!                Hel_ib=get_helement(nI,nJ2,IC_ib,iLutnI,iLutnJ2)
!                Hel_ja=get_helement(nI2,nJ,IC_ja,iLutnI2,iLutnJ)
!                Hel_jb=get_helement(nI2,nJ2,IC_jb,iLutnI2,iLutnJ2)
!
!                if(abs(HEl_ia-HEl_jb).gt.1.0e-7_dp_dp) then
!                    call stop_all(t_r,"Matrix element ia .ne. jb")
!                endif
!                if(abs(HEl_ib-HEl_ja).gt.1.0e-7_dp_dp) then
!                    call stop_all(t_r,"Matrix element ib .ne. ja")
!                endif
!
!                !Calculate PGens...
!                !First need Ex...
!                if(IC_ia.le.2) then
!                    Ex_ia=0
!                    Ex_ia(1,1)=IC_ia
!                    call GetBitExcitation(iLutnI,iLutnJ,Ex_ia,tSign)
!                    call construct_class_counts(nI,CC_ia,CCU_ia)
!                    call CalcNonUniPGen(nI,Ex_ia,IC_ia,CC_ia,CCU_ia,0.9_8,pGen_ia)
!                else
!                    pGen_ia=0.0_dp
!                endif
!                if(IC_ib.le.2) then
!                    Ex_ib=0
!                    Ex_ib(1,1)=IC_ib
!                    call GetBitExcitation(iLutnI,iLutnJ2,Ex_ib,tSign)
!                    call construct_class_counts(nI,CC_ib,CCU_ib)
!                    call CalcNonUniPGen(nI,Ex_ib,IC_ib,CC_ib,CCU_ib,0.9_8,pGen_ib)
!                else
!                    pGen_ib=0.0_dp
!                endif
!                if(IC_ja.le.2) then
!                    Ex_ja=0
!                    Ex_ja(1,1)=IC_ja
!                    call GetBitExcitation(iLutnI2,iLutnJ,Ex_ja,tSign)
!                    call construct_class_counts(nI2,CC_ja,CCU_ja)
!                    call CalcNonUniPGen(nI,Ex_ja,IC_ja,CC_ja,CCU_ja,0.9_8,pGen_ja)
!                else
!                    pGen_ja=0.0_dp
!                endif
!                if(IC_jb.le.2) then
!                    Ex_jb=0
!                    Ex_jb(1,1)=IC_jb
!                    call GetBitExcitation(iLutnI2,iLutnJ2,Ex_jb,tSign)
!                    call construct_class_counts(nI2,CC_jb,CCU_jb)
!                    call CalcNonUniPGen(nI2,Ex_jb,IC_jb,CC_jb,CCU_jb,0.9_8,pGen_jb)
!                else
!                    pGen_jb=0.0_dp
!                endif
!                
!                if(abs(PGen_ia-PGen_jb).gt.1.0e-7_dp_dp) then
!                    call stop_all(t_r,"PGen ia .ne. jb")
!                endif
!                if(abs(PGen_ib-PGen_ja).gt.1.0e-7_dp_dp) then
!                    call stop_all(t_r,"PGen ib .ne. ja")
!                endif
!
!            enddo
!        enddo
!
!        write(6,*) "Now testing matrix elements from closed shell MI functions..."
!
!        do i=1,nSelfDets
!            do j=1,nPairedDets
!                
!                !For this i,j pair, isolate all four hamiltonian matrix element between determinants.
!                iLutnI(:)=SelfInvDets(:,i)
!                iLutnJ(:)=PairedDets(:,j)
!
!                call InvertMomBitDet(iLutnJ,iLutnJ2)
!
!                IC_ia = FindBitExcitLevel(iLutnI,iLutnJ,NEl)
!                IC_ib = FindBitExcitLevel(iLutnI,iLutnJ2,NEl)
!
!                if(IC_ia.ne.IC_ib) call stop_all(t_r,"I,A excitation level not the same as I,B excitation")
!
!                call decode_bit_det(nI,iLutnI)
!                call decode_bit_det(nJ,iLutnJ)
!                call decode_bit_det(nJ2,iLutnJ2)
!
!                !Calculate matrix elements between them...
!                Hel_ia=get_helement(nI,nJ,IC_ia,iLutnI,iLutnJ)
!                Hel_ib=get_helement(nI,nJ2,IC_ib,iLutnI,iLutnJ2)
!
!                if(abs(HEl_ia-HEl_ib).gt.1.0e-7_dp_dp) then
!                    call stop_all(t_r,"Matrix element ia .ne. ib")
!                endif
!
!                !Calculate PGens...
!                !First need Ex...
!                if(IC_ia.le.2) then
!                    Ex_ia=0
!                    Ex_ia(1,1)=IC_ia
!                    call GetBitExcitation(iLutnI,iLutnJ,Ex_ia,tSign)
!                    call construct_class_counts(nI,CC_ia,CCU_ia)
!                    call CalcNonUniPGen(nI,Ex_ia,IC_ia,CC_ia,CCU_ia,0.9_8,pGen_ia)
!                else
!                    pGen_ia=0.0_dp
!                endif
!                if(IC_ib.le.2) then
!                    Ex_ib=0
!                    Ex_ib(1,1)=IC_ib
!                    call GetBitExcitation(iLutnI,iLutnJ2,Ex_ib,tSign)
!                    call construct_class_counts(nI,CC_ib,CCU_ib)
!                    call CalcNonUniPGen(nI,Ex_ib,IC_ib,CC_ib,CCU_ib,0.9_8,pGen_ib)
!                else
!                    pGen_ib=0.0_dp
!                endif
!                
!                if(abs(PGen_ia-PGen_ib).gt.1.0e-7_dp_dp) then
!                    call stop_all(t_r,"PGen ia .ne. ib")
!                endif
!
!            enddo
!        enddo
!
!        deallocate(PairedDets,SelfInvDets)
!
!    end subroutine TestMomInvInts

end module MomInvRandExcit
