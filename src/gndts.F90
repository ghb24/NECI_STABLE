module gndts_mod

    use constants, only: n_int
    use SystemData, only: BasisFN, BasisFNSize, Symmetry, NullBasisFn, LzTot, &
                          tFixLz
    use sort_mod
    use sym_mod
    use DetBitOps, only: count_open_orbs,encodebitdet
    use CalcData, only: tTruncNOpen,trunc_nopen_max
    use bit_reps, only: niftot,encode_bit_rep
    use constants, only: n_int
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
    end subroutine

    ! The recursive routine in GNDTS.  IELEC is the current electron
    recursive subroutine gndts_r (nel, nhg, brr, nBasisMax, nMrks, tCount, &
                                  G1, tSpn, LMS, tParity, SymRestrict, ii, &
                                  iElecs, iSpins, iElec, kJ)

        integer, intent(in) :: nel, nhg, LMS, brr(nhg), nBasisMax(5,*)
        integer, value :: iElec
        integer, intent(inout) :: ii, nMrks(:, :), iElecs(nel), iSpins(nel)
        type(BasisFn), intent(in) :: G1(:), kJ, SymRestrict
        logical, intent(in) :: tCount, tSpn, tParity
        logical :: tSkip
        integer(n_int) :: ilut(0:NIfTot)


        type(BasisFn) :: kI
        integer :: nI(nel), ist, iel, i

        ! This will hold the parity of the basis fn
        ! Iterate of the the spin basis functions available for this electron
        ! I_(N+1) = I_N +1 .. NHG
        IST=1
        tSkip = .false.
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
                if(tTruncNOpen) then
                    call EncodeBitDet(nI,ilut)
                    !write(6,*) 'ilut:',ilut
                    !write(6,*) "Number of open orbitals: ",count_open_orbs(ilut),trunc_nopen_max
                    if(count_open_orbs(ilut)>trunc_nopen_max) then
                        tSkip = .true.
                    else
                        tSkip = .false.
                    endif
                endif

                IF(LCHKSYM(KI,KJ)) THEN
                    IF((.not.tFixLz).or.(KI%Ml.eq.KJ%Ml)) THEN
                        IF (.not.tSkip) THEN
                            II=II+1
                            ! If we're generating rather than counting
                            if (.not. tCount) &
                                nMrks(:, ii) = nI
                        ENDIF
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

    end subroutine

    subroutine gndts_all_sym_this_proc(ilut_list, tCount, ndets)

        use SystemData, only: nel

        integer(n_int), intent(out) :: ilut_list(0:,:)
        logical, intent(in) :: tCount
        integer, intent(out) :: ndets
        integer :: nI(nel)

        ndets = 0
        nI = 0
        call gndts_all_sym_this_proc_r(ilut_list, tCount, ndets, nI, ielec=1)

    end subroutine gndts_all_sym_this_proc

    recursive subroutine gndts_all_sym_this_proc_r(ilut_list, tCount, ndets, nI, ielec)

        use bit_rep_data, only: NIfD, NIfTot
        use DetBitOps, only: EncodeBitDet
        use load_balance_calcnodes, only: DetermineDetNode
        use Parallel_neci, only: iProcIndex
        use SystemData, only: nel, nbasis

        integer(n_int), intent(inout) :: ilut_list(0:,:)
        logical, intent(in) :: tCount
        integer, intent(inout) :: ndets
        integer, intent(inout) :: nI(nel)
        integer, value :: ielec
        integer :: starting_orb, iorb, proc
        integer(n_int) :: ilut(0:NIfTot)

        ilut = 0_n_int
        starting_orb = 1
        if (ielec > 1) starting_orb = nI(ielec-1) + 1

        do iorb = starting_orb, nbasis
            nI(ielec) = iorb
            ! If we're on the last spin.
            if (ielec == nel) then
                proc = DetermineDetNode(nel,nI,0)
                if (proc == iProcIndex) then
                    ndets = ndets + 1
                    if (.not. tCount) then
                        call EncodeBitDet(nI, ilut)
                        ilut_list(0:NIfD, ndets) = ilut(0:NIfD)
                    end if
                end if
            else
                call gndts_all_sym_this_proc_r(ilut_list, tCount, ndets, nI, ielec+1)
            end if
        end do

    end subroutine gndts_all_sym_this_proc_r

end module
