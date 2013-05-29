! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
module hphf_integrals
    use constants, only: dp,n_int,sizeof_int
    use SystemData, only: NEl, nBasisMax, G1, nBasis, Brr, tHub, ECore, &
                          ALat, NMSH, tOddS_HPHF, modk_offdiag
    use IntegralsData, only: UMat,FCK,NMAX
    use HPHFRandExcitMod, only: FindDetSpinSym, FindExcitBitDetSym
    use DetBitOps, only: DetBitEQ, FindExcitBitDet, FindBitExcitLevel, &
                         TestClosedShellDet, CalcOpenOrbs
    use sltcnd_mod, only: sltcnd, sltcnd_excit
    use bit_reps, only: NIfD, NIfTot, NIfDBO
    implicit none

    interface hphf_off_diag_helement
        module procedure hphf_off_diag_helement_norm
        module procedure hphf_off_diag_helement_spawn
    end interface

    contains

! NB AJWT has also cannibalized this for CCMC where he just passes in the HElGen.
    function hphf_spawn_sign (nI, nJ, iLutI, iLutJ, ic, ex, &
                                  tParity, HElGen) result (hel)
        integer, intent(in) :: nI(nel), nJ(nel), ic, ex(2,2)
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        logical, intent(in) :: tParity
        integer :: iUnused
        logical :: lUnused
        HElement_t :: hel
        HElement_t, intent(in) :: HElGen

        hel=HElGen

        ! Avoid warnings
        iUnused = IC; iUnused = ex(1,1); iUnused = nI(1); iUnused = nJ(1)
        iUnused = int(iLutI(0),sizeof_int); iUnused = int(iLutJ(0),sizeof_int)
        lUnused = tParity

    end function

    ! TODO: comment as to why!
    function hphf_off_diag_helement_spawn (nI, nJ, iLutI, iLutJ, ic, ex, &
                                           tParity, HElGen) result (hel)
        integer, intent(in) :: nI(nel), nJ(nel), ic, ex(2,2)
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        logical, intent(in) :: tParity
        integer :: iUnused
        logical :: lUnused
        HElement_t :: hel
        HElement_t , intent(in) :: HElGen

        hel = hphf_off_diag_helement_norm (nI, nJ, iLutI, iLutJ)

        if (IC /= 0 .and. modk_offdiag) &
            hel = -abs(hel)

        ! Avoid warnings
        iUnused = IC; iUnused = ex(1,1); lUnused = tParity

    end function

    function hphf_off_diag_helement_norm (nI, nJ, iLutnI, iLutnJ) result(hel)

        ! Find the  between two half-projected hartree-fock 
        ! determinants (different ones). NI and nJ have to be uniquely 
        ! chosen, so that their spin-coupled determinant will not arise.
        !
        ! In:  nI, nJ         - Determinants to consider
        !      iLutnI, iLutnJ - Bit representations of I,J
        ! Ret: hel          - The calculated matrix element

        integer, intent(in) :: nI(nel), nJ(nel)
        integer(kind=n_int), intent(in) :: iLutnI(0:NIfTot), iLutnJ(0:NIfTot)
        HElement_t :: hel

        integer :: nI2(nel), iUnused
        integer(kind=n_int) :: iLutnI2(0:NIfTot)
        integer :: ExcitLevel, OpenOrbsI, OpenOrbsJ, Ex(2,2)
        HElement_t :: MatEl2
        logical :: tSign

        ! Avoid warnings
        iUnused = nJ(1)

        if (DetBitEQ(iLutnI, iLutnJ, NIfDBO)) then
            ! Do not allow a 'diagonal' matrix element. The problem is 
            ! that the HPHF excitation generator can generate the same HPHF 
            ! function. We do not want to allow spawns here.
            hel = (0)
            return
        endif

        hel = sltcnd (nI, iLutnI, iLutnJ)
        if (TestClosedShellDet(iLutnI)) then
            if(tOddS_HPHF) then
                !For odd S states, all matrix elements to CS determinants should be 0
                hel = 0.0_dp
            elseif (.not. TestClosedShellDet(iLutnJ)) then
                ! Closed shell --> Open shell, <X|H|Y> = 1/sqrt(2) [Hia + Hib]
                ! or with minus if iLutnJ has an odd number of spin orbitals.
                ! OTHERWISE Closed shell -> closed shell. Both the alpha and 
                ! beta of the same orbital have been moved to the same new 
                ! orbital. The matrix element is the same as before.
                hel = hel * (sqrt(2.0_dp))
            endif
        else
            if (TestClosedShellDet(iLutnJ)) then
                if(tOddS_HPHF) then
                    !For odd S states, all matrix elements to CS determinants should be 0
                    hel = 0.0_dp
                else
                    ! Open shell -> Closed shell. If one of
                    ! the determinants is connected, then the other is connected 
                    ! with the same IC & matrix element
                    hel = hel * sqrt(2.0_dp)
                endif
            else
                ! Open shell -> Open shell. Find the spin pair of nJ.
                call FindExcitBitDetSym(iLutnI, iLutnI2)
                ExcitLevel = FindBitExcitLevel(iLutnI2, ilutnJ, 2)

                if (ExcitLevel.le.2) then
                    ! We need to find out whether the nJ HPHF wavefunction is 
                    ! symmetric or antisymmetric. This is dependant on the 
                    ! number of open shell orbitals and total spin of the wavefunction.
                    call FindDetSpinSym(nI, nI2, nel)
                    call CalcOpenOrbs(iLutnJ, OpenOrbsJ)

                    ! Original HPHF is antisymmetric if OpenOrbs is odd (and S even), 
                    ! or symmetric if it is even.
                    ! If S is odd, then HPHF is Symmetric if OpenOrbs is odd, and 
                    ! antisymmetric if it is even.
                    call CalcOpenOrbs(iLutnI,OpenOrbsI)
                    Ex(1,1)=ExcitLevel
                    call GetBitExcitation(iLutnI2,iLutnJ,Ex,tSign)

                    MatEl2 = sltcnd_excit (nI2, ExcitLevel, Ex, tSign)

                    if(tOddS_HPHF) then
                        if (((mod(OpenOrbsI,2) == 1).and.(mod(OpenOrbsJ,2) == 1))&
                            .or. ((mod(OpenOrbsI,2) == 1) .and. &
                                  (mod(OpenOrbsJ,2) == 0))) then
                            hel = hel + MatEl2
                        else
                            hel = hel - MatEl2
                        endif
                    else
                        if (((mod(OpenOrbsI,2) == 0).and.(mod(OpenOrbsJ,2) == 0))&
                            .or. ((mod(OpenOrbsI,2) == 0) .and. &
                                  (mod(OpenOrbsJ,2) == 1))) then
                            hel = hel + MatEl2
                        else
                            hel = hel - MatEl2
                        endif
                    endif
                endif
            endif
        endif
    end function


    function hphf_diag_helement (nI, iLutnI) result(hel)

        ! Find the diagonal HElment for a half-projected hartree-fock 
        ! determinant.
        !
        ! In:  nI      - Determinant to consider
        !      iLutnI  - Bit representation of I
        ! Ret: hel   - The calculated matrix element

        integer, intent(in) :: nI(nel) 
        integer(kind=n_int), intent(in) :: iLutnI(0:NIfTot)
        HElement_t :: hel

!        integer :: nI2(nel)
        integer(kind=n_int) :: iLutnI2(0:NIfTot)
        integer :: ExcitLevel, OpenOrbs
        HElement_t :: MatEl2

        hel = sltcnd_excit (nI, 0)
        if (.not. TestClosedShellDet(iLutnI)) then
            ! <i|H|i> = <j|H|j>, so no need to calculate both.
            ! <X|H|X> = 1/2 [ <i|H|i> + <j|H|j> ] + <i|H|j> where i and j are
            ! the two spin-coupled dets which make up X. In the case of the 
            ! antisymmetric pair, the cross term is subtracted.

            ! See if there is a cross-term
            call FindExcitBitDetSym(iLutnI, iLutnI2)
            ExcitLevel = FindBitExcitLevel(iLutnI, iLutnI2, 2)
            if (ExcitLevel.le.2) then
                call CalcOpenOrbs (iLutnI, OpenOrbs)
!                call FindDetSpinSym (nI, nI2, nel)
                MatEl2 = sltcnd (nI,  iLutnI, iLutnI2)

                if (tOddS_HPHF) then
                    if (mod(OpenOrbs,2).eq.1) then
                        ! Subtract cross terms if determinant is antisymmetric.
                        hel = hel + MatEl2
                    else
                        hel = hel - MatEl2
                    endif
                else
                    if (mod(OpenOrbs,2).eq.1) then
                        ! Subtract cross terms if determinant is antisymmetric.
                        hel = hel - MatEl2
                    else
                        hel = hel + MatEl2
                    endif
                endif
            endif
        endif

        hel = hel + (ECore)
    end function hphf_diag_helement

    pure function hphf_sign (ilut) result(sgn)

        ! Is this HPHF  1/sqrt(2)*[X + X'], or 1/sqrt(2)*[X - X']
        ! Returns +-1 respectively

        integer :: sgn, open_orbs
        integer(n_int), intent(in) :: ilut(0:NIfTot)

        call CalcOpenOrbs(ilut, open_orbs)

        if ((mod(open_orbs, 2) == 0) .neqv. tOddS_HPHF) then
            sgn = 1
        else
            sgn = -1
        endif

    end function

end module
