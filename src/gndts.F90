module gndts_mod

    use SystemData, only: BasisFN, BasisFNSize, Symmetry, NullBasisFn, LzTot, &
                          tFixLz
    use sort_mod
    use sym_mod
    implicit none

contains

    ! A recursive version of GNDTS which doesn't require #defines
    ! This sets up variables, and calls the recursion 
    subroutine gndts (nel, nhg, brr, nBasisMax, nMrks, tCount, G1, tSpn, LMS, &
                      tParity, SymRestrict, ii, IFDet)

        integer, intent(in) :: nel, nhg, lms, brr(nhg)
        integer, intent(in) :: nBasisMax(5,*)
        type(BasisFN), intent(in) :: G1(:), SymRestrict
        logical, intent(in) :: tCount, tSpn, tParity
        integer, intent(out) :: ii, IFDet
        integer, intent(inout) :: nMrks(:, :)

        type(BasisFN) :: kJ
        ! Temporary arrays to hold the determinants and spins each electron
        ! is on at the moment.
        integer :: iElecs(nel), iSpins(nel)

        ! Total number of determined determinants
        ii = 0
        if (tParity) then
            kJ=SymRestrict
        else
            kJ=NullBasisFn
        endif
        if (tSpn) then
            kJ%Ms=LMS
        else
            kJ%Ms=0 
        endif
        if (tFixLz) then
            kJ%Ml=LzTot
        endif


        ! Start with the first electron
        call gndts_r (nel, nhg, brr, nBasisMax, nMrks, tCount, G1, tSpn, LMS, &
                      tparity, SymRestrict, ii, iElecs, iSpins, 1, kJ)
        IFDet = 1
    end

    ! The recursive routine in GNDTS.  IELEC is the current electron
    recursive subroutine gndts_r (nel, nhg, brr, nBasisMax, nMrks, tCount, &
                                  G1, tSpn, LMS, tParity, SymRestrict, ii, &
                                  iElecs, iSpins, iElec, kJ)

        integer, intent(in) :: nel, nhg, LMS, brr(nhg), nBasisMax(5,*), iElec
        integer, intent(inout) :: ii, nMrks(:, :), iElecs(nel), iSpins(nel)
        type(BasisFn), intent(in) :: G1(:), kJ, SymRestrict
        logical, intent(in) :: tCount, tSpn, tParity

        type(BasisFn) :: kI
        integer :: nI(nel), ist, iel, i

        ! This will hold the parity of the basis fn
        ! Iterate of the the spin basis functions available for this electron
        ! I_(N+1) = I_N +1 .. NHG
        IST=1
        IF(IELEC.GT.1) IST=IELECS(IELEC-1)+1
        DO IEL=IST,NHG
            IELECS(IELEC)=IEL
            ISPINS(IELEC)=G1(BRR(IELECS(IELEC)))%Ms
            ! If we're on the last electron
            IF(NEL.EQ.IELEC) THEN
                DO I=1,NEL
                    NI(I)=BRR(IELECS(I))
                ENDDO
                call sort (nI)
                CALL GETSYM(NI,NEL,G1,NBASISMAX,KI)
                CALL GetLz(NI,NEL,KI%Ml) 
                IF(.NOT.TPARITY) THEN
                    DO I=1,3
                        KI%k(I)=0
                    ENDDO
                    KI%Sym%s=0 ! Ignore symmetry in generation of determinants.
                ENDIF
                IF(.NOT.TSPN) KI%Ms=0
                !        WRITE(6,*) KI
                IF(LCHKSYM(KI,KJ)) THEN
                    IF((.not.tFixLz).or.(KI%Ml.eq.KJ%Ml)) THEN
                        II=II+1
                        ! If we're generating rather than counting
                        if (.not. tCount) &
                            nMrks(:, ii) = nI
                    ENDIF
                ENDIF
            ELSE
                ! If we're not on the last electron
                ! We recurse over the next electrons
                call gndts_r (nel, nhg, brr, nBasisMax, nMrks, tCount, G1, &
                              tSpn, LMS, tParity, SymRestrict, ii, iElecs, &
                              iSpins, iElec+1, kJ)
            ENDIF
        ENDDO    

    end

end module

