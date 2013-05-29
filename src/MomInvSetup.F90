! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
#include "macros.h"

Module MomInv
    use constants, only: dp,n_int,end_n_int
    use bit_reps, only: NIfD, NIfTot, NIfDBO, decode_bit_det
    use DetBitOps, only: DetBitEQ,DetBitLT
    use SystemData, only: NEl,G1,nBasis,Brr,Arr,tFixLz,LzTot,G1,tAntisym_MI
    use Parallel_neci
    use sort_mod
    use SymExcitDataMod, only: MomInvSymOrb  

    implicit none

    contains

    subroutine SetupMomInv()
        use util_mod, only: get_free_unit
        implicit none
        character(len=*), parameter :: t_r='SetupMomInv'
        logical :: exists,tOrbTaken(nBasis)
        integer :: LzMap_unit,Orb,OrbPartner,i,j,ierr
    
        if(.not.tFixLz) then
            call stop_all(t_r,"Cannot use MomInv functions without Lz symmetry")
        endif
        if(LzTot.ne.0) then
            call stop_all(t_r,"Cannot use MomInv functions with non-zero total momentum")
        endif

        !This array will contain the inverse orbitals
        allocate(MomInvSymOrb(1:nBasis),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Allocation error")
        MomInvSymOrb=0

        if(iProcIndex.eq.Root) then
            inquire(file='LzMapping',exist=exists)
            if(exists) then
                write(6,*) "LzMapping file found for momentum inversion orbitals"
                LzMap_unit = get_free_unit()
                open(LzMap_unit,file='LzMapping',status='old',action='read')
                do i=1,nBasis
                    read(LzMap_unit) Orb,OrbPartner
                    if(MomInvSymOrb(Orb).ne.0) then
                        call stop_all(t_r,"Orbital already assigned")
                    endif
                    MomInvSymOrb(Orb)=OrbPartner
                    if((G1(Orb)%Ml).ne.-((G1(MomInvSymOrb(Orb))%Ml))) then
                        call stop_all(t_r,"Error in LzMapping - Lz")
                    endif
                    if((G1(Orb)%Ms).ne.(G1(MomInvSymOrb(Orb))%Ms)) then
                        call stop_all(t_r,"Error in LzMapping - Ms")
                    endif
                enddo

                close(LzMap_unit)

            else
                !Create mapping

                tOrbTaken(:)=.false.
                do i=1,nBasis
                    !Find partner for orbital i, just by finding the next symmetry allowed pair
                    !This will need to be refined for molecular systems

                    if(tOrbTaken(i)) cycle  !Already found pair

                    if(G1(i)%Ml.eq.0) then
                        tOrbTaken(i)=.true. !Self inverse orbital
                        MomInvSymOrb(i)=i
                    else

                        do j=i+2,nBasis,2
                            !Loop over remaining orbitals of correct spin
!                            write(6,*) i,j,tOrbTaken(j),G1(i)%Ml,-G1(j)%Ml
                            if((G1(i)%Ml.eq.-G1(j)%Ml).and.(.not.tOrbTaken(j))) then
                                !Found pair
                                MomInvSymOrb(i)=j
                                MomInvSymOrb(j)=i
                                tOrbTaken(i)=.true.
                                tOrbTaken(j)=.true.
                                exit
                            endif
                        enddo
                        
                        if(j.gt.nBasis) then
                            call stop_all(t_r,"Could not find correct pair for orbital")
                        endif

                    endif

                enddo

                do i=1,nBasis
                    if(.not.tOrbTaken(i).or.(MomInvSymOrb(i).eq.0)) then
                        call stop_all(t_r,"Not all orbitals paired")
                    endif
                enddo


            endif

            !Check Mapping
            do i=1,nBasis
                if(MomInvSymOrb(i).eq.0) then
                    call stop_all(t_r,"Orbital not assigned")
                endif
                if((G1(i)%Ms).ne.(G1(MomInvSymOrb(i))%Ms)) then
                    call stop_all(t_r,"Ms values not matches")
                endif
                if((G1(i)%Ml).ne.(-(G1(MomInvSymOrb(i))%Ml))) then
                    call stop_all(t_r,"Ml pairs not assigned correctly")
                endif

                !Pairs should be (approximately) degenerate
                if(abs(Arr(i,2)-Arr(MomInvSymOrb(i),2)).gt.2.0e-4_dp) then
                    call stop_all(t_r,"Ml pairs not degenerate")
                endif
            enddo
        endif

        call MPIBCast(MomInvSymOrb)

        if(tAntisym_MI) then
            write(6,*) "Hilbert space will be constructed from antisymmetric momentum-coupled determinants"
        else
            write(6,*) "Hilbert space will be constructed from symmetric momentum-coupled determinants"
        endif
        write(6,*) "Sucessfully paired orbitals by momentum."
        write(6,*) "Orbital     Lz_Paired Orbital"
        do i=1,nBasis
            write(6,"(2I10)") i,MomInvSymOrb(i)
        enddo
                
    end subroutine SetupMomInv

    !Create momentum paired nI
    subroutine InvertMomDet(nI,MomSymDet)
        implicit none
        integer , intent(in) :: nI(NEl)
        integer , intent(out) :: MomSymDet(NEl)
        integer :: i

        do i=1,nel
            MomSymDet(i)=MomInvSymOrb(nI(i))
        enddo
        call sort(MomSymDet)

    end subroutine InvertMomDet

    !Create momentum paired iLut
    subroutine InvertMomBitDet(iLut,MomSymiLut, nI)
        implicit none
        integer(n_int) , intent(in) :: iLut(0:NIfTot)
        integer(n_int) , intent(out) :: MomSymiLut(0:NIfTot)
        integer :: i,j,nITmp(nel)
        integer, intent(in), optional :: nI(NEl)
        
        MomSymiLut(:)=0
        if(present(nI)) then
            do i=1,nel
                j=MomInvSymOrb(nI(i))
                set_orb(MomSymiLut,j)
            enddo
        else
            call decode_bit_det(nITmp,iLut)
            !No need to sort, since we are encoding straight away
            do i=1,nel
                j=MomInvSymOrb(nITmp(i))
                set_orb(MomSymiLut,j)
            enddo
        endif


    end subroutine InvertMomBitDet

    !Return the inverse determinant in both representations
    subroutine InvertMomDetBoth(nI,nISym,iLut,iLutSym)
        implicit none
        integer, intent(in) :: nI(NEl)
        integer, intent(out) :: nISym(NEl)
        integer(n_int), intent(in) :: iLut(0:NIfTot)
        integer(n_int), intent(out) :: iLutSym(0:NIfTot)
        integer :: i
            
        iLutSym(:)=0
        do i=1,nel
            nISym(i)=MomInvSymOrb(nI(i))
            set_orb(iLutSym,nISym(i))
        enddo
        call sort(nISym)

    end subroutine InvertMomDetBoth

    !Return the allowed Momentum Inverse determinant, with optional calculation of either/both
    !of the symmetry determinants.
    subroutine CalcMomAllowedBitDet(nI,nISym,iLut,iLutSym,tCalcnISym,tCalciLutSym,tSwapped)
        implicit none
        logical, intent(in) :: tCalcnISym,tCalciLutSym
        logical, intent(out) :: tSwapped
        integer, intent(inout) :: nI(NEl)
        integer, intent(inout) :: nISym(NEl)
        integer(n_int), intent(inout) :: iLut(0:NIfTot)
        integer(n_int), intent(inout) :: iLutSym(0:NIfTot)
        integer :: i,nITemp(NEl)
        integer(n_int) :: iLutTemp(0:NIfTot)

        if(tCalcnISym.and.tCalciLutSym) then
            !Calculate both representations
            call InvertMomDetBoth(nI,nISym,iLut,iLutSym)
        elseif(tCalcnISym) then
            !Calculate nI representation
            call decode_bit_det(nISym,iLutSym)
        elseif(tCalciLutSym) then
            !Calculation iLut representation
            iLutSym(:)=0
            do i=1,nel
                set_orb(iLutSym,nISym(i))
            enddo
        endif
            
        !Use ilut representation to return correct determinant
        i=DetBitLT(iLut,iLutSym,NIfD)
        if(i.eq.1) then
            !swap
            iLutTemp(:)=iLut(:)
            iLut(:)=iLutSym(:)
            iLutSym(:)=iLutTemp(:)
            nITemp(:)=nI(:)
            nI(:)=nISym(:)
            nISym(:)=nITemp(:)
            tSwapped=.true.
        elseif(i.eq.0) then
            call stop_all("CalcMomAllowedBitDet","Shouldn't have self-inverse in here")
        else
            tswapped=.false.
        endif

    end subroutine CalcMomAllowedBitDet

    subroutine ReturnMomAllowedBitDet(iLut,iLutSym,tSwapped)
        implicit none
        integer(n_int), intent(inout) :: iLut(0:NIfTot),iLutSym(0:NIfTot)
        logical, intent(out) :: tSwapped
        integer(n_int) :: iLutTemp(0:NIfTot)
        integer :: i
        
        i=DetBitLT(iLut,iLutSym,NIfD)
        if(i.eq.1) then
            !swap
            iLutTemp(:)=iLut(:)
            iLut(:)=iLutSym(:)
            iLutSym(:)=iLutTemp(:)
            tSwapped=.true.
        elseif(i.eq.0) then
            call stop_all("ReturnMomAllowedBitDet","Shouldn't have self-inverse in here")
        else
            tswapped=.false.
        endif

    end subroutine ReturnMomAllowedBitDet

    subroutine ReturnMomAllowedDet(nI,nISym,iLut,iLutSym,tSwapped)
        implicit none
        integer(n_int), intent(inout) :: iLut(0:NIfTot),iLutSym(0:NIfTot)
        integer(n_int) :: iLutTemp(0:NIfTot)
        integer, intent(inout) :: nI(nEl),nISym(nEl)
        integer :: nTemp(NEl),i
        logical, intent(out) :: tSwapped

        i=DetBitLT(iLut,iLutSym,NIfD)
        if(i.eq.1) then
            !swap
            iLutTemp(:)=iLut(:)
            iLut(:)=iLutSym(:)
            iLutSym(:)=iLutTemp(:)
            nTemp(:)=nI(:)
            nI(:)=nISym(:)
            nISym(:)=nTemp(:)
            tSwapped=.true.
        elseif(i.eq.0) then
            call stop_all("ReturnMomAllowedDet","Shouldn't have self-inverse in here")
        else
            tswapped=.false.
        endif

    end subroutine ReturnMomAllowedDet

    !Returns true if the function is an 'allowed' MI function
    !(i.e. the determinant which would be stored in an FCIMC run) 
    logical function IsAllowedMI(nI,iLutnI)
        implicit none
        integer, intent(in) :: nI(NEl)
        integer(n_int), intent(in) :: iLutnI(0:NIfTot)
        integer(n_int) :: MomSymiLut(0:NIfTot),iLutTmp(0:NIfTot)
        logical :: tSwapped

        if(IsMomSelfInv(nI,iLutnI)) then
            IsAllowedMI=.true.
            return
        else
            call InvertMomBitDet(iLutnI,MomSymiLut, nI)
            iLutTmp(:)=iLutnI
            call ReturnMomAllowedBitDet(iLutTmp,MomSymiLut,tSwapped)
            if(tSwapped) then
                IsAllowedMI=.false.
            else
                IsAllowedMI=.true.
            endif
        endif

    end function IsAllowedMI


    !Routine to return whether a determinant is a momentum self inverse or not,
    !using both bit and natural ordered representations.
    !Faster than IsBitMomSelfInv
    pure logical function IsMomSelfInv(nI,iLutnI)
        implicit none
        integer, intent(in) :: nI(NEl)
        integer(n_int), intent(in) :: iLutnI(0:NIfTot)
        integer :: i

        do i=1,NEl
            if(IsNotOcc(iLutnI,MomInvSymOrb(nI(i)))) then
                IsMomSelfInv=.false.
                return
            endif
        enddo

        IsMomSelfInv=.true.

    end function IsMomSelfInv

    pure logical function IsBitMomSelfInv(iLut)
        implicit none
        integer(n_int) , intent(in) :: iLut(0:NIfTot)
        integer :: i,j,elec,Orb

        IsBitMomSelfInv=.true.
        elec=0
        do i=0,NIfD
            do j=0,end_n_int
                if(btest(iLut(i),j)) then
                    !This orbital occupied
                    Orb=MomInvSymOrb((i*bits_n_int)+(j+1))
                    elec=elec+1
                    if(IsNotOcc(iLut,Orb)) then
                        IsBitMomSelfInv=.false.
                        return
                    endif
                    if(elec.eq.NEl) return
                endif
            enddo
        enddo

    end function IsBitMomSelfInv

end module MomInv
