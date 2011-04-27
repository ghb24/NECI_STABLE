#include "macros.h"

Module MomInv
    use constants, only: dp,n_int,end_n_int
    use bit_reps, only: NIfD, NIfTot, NIfDBO
    use DetBitOps, only: DetBitEQ
    use SystemData, only: NEl,G1,nBasis,Brr,Arr,tFixLz,LzTot,G1
    use Parallel
    use sort_mod
    use MomInvData

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
                if(abs(Arr(i,2)-Arr(MomInvSymOrb(i),2)).gt.2.D-4) then
                    call stop_all(t_r,"Ml pairs not degenerate")
                endif
            enddo
        endif

        call MPIBCast(MomInvSymOrb)

        write(6,*) "Sucessfully paired orbitals by momentum."
        write(6,*) "Orbital     Lz_Paired Orbital"
        do i=1,nBasis
            write(6,"(2I10)") i,MomInvSymOrb(i)
        enddo
                
    end subroutine SetupMomInv

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
