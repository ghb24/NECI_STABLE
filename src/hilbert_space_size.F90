#include "macros.h"

module hilbert_space_size

    use constants, only: dp, int64, n_int, bits_n_int, sizeof_int, stdout
    use util_mod, only: get_free_unit, operator(.div.), choose_i64
    implicit none

contains

!This routine stochastically finds the size of a given excitation level (iExcitLevelTest), or
!lower, with equal probably given to all determinants.
    SUBROUTINE FindSymMCSizeExcitLevel(IUNIT)
        use SymData, only: TwoCycleSymGens
        use SystemData, only: nEl, G1, nBasis, iMCCalcTruncLev
        use SystemData, only: tUEG, tHPHF, tHub
        use SystemData, only: CalcDetCycles, CalcDetPrint, tFixLz
        use DeterminantData, only: FDet
        use dSFMT_interface
        use soft_exit, only: ChangeVars, tSoftExitFound
        use Parallel_neci
        use DetBitops, only: EncodeBitDet
        use bit_rep_data, only: NIfTot
        IMPLICIT NONE
        INTEGER, intent(in) :: IUNIT
        INTEGER :: j, SpatOrbs, FDetMom, Attempts, iExcitLevTest, i, ExcitLev
        INTEGER(KIND=n_int) :: FDetiLut(0:NIfTot), iLut(0:NIfTot)
        INTEGER :: FDetSym, TotalSym, TotalMom, alpha, beta, ierr, Momx, Momy
        integer(int64) :: Accept, AcceptAll, TotalAttempts, TotalAttemptsAll
        integer(int64) :: ExcitBin(0:iMCCalcTruncLev), ExcitBinAll(0:iMCCalcTruncLev)
        real(dp) :: ExcitLevBias(0:iMCCalcTruncLev)
        real(dp) :: FullSpace, r, Frac, SymSpace
        real(dp) :: SizeLevel(0:iMCCalcTruncLev)
        LOGICAL :: tDummy, tDummy2

        iExcitLevTest = iMCCalcTruncLev

        write(IUNIT, "(A,I6)") "Calculating MC size of symmetry-allowed " &
            //"space of excitation levels up to level: ", iExcitLevTest

        write(IUNIT, "(I18,A,I18,A)") CalcDetCycles, " MC cycles will be used, and " &
            //"statistics printed out every ", CalcDetPrint, " cycles."
        FDetSym = 0
        FDetMom = 0
        ExcitBin(:) = 0
        ExcitBinAll(:) = 0

        do i = 1, NEl
            FDetSym = IEOR(FDetSym, INT(G1(FDet(i))%Sym%S, sizeof_int))
            IF (tFixLz) FDetMom = FDetMom + G1(FDet(i))%Ml
        end do
        CALL EncodeBitDet(FDet, FDetiLut)

        write(IUNIT, *) "Symmetry of HF determinant is: ", FDetSym
        IF (tFixLz) THEN
            write(IUNIT, *) "Imposing momentum sym on size calculation"
            write(IUNIT, *) "Momentum of HF determinant is: ", FDetMom
        end if
        IF (tHPHF) THEN
            write(stdout, *) "Imposing time-reversal symmetry (HPHF) on " &
                //"size of space calculation"
        end if

        Accept = 0
        AcceptAll = 0
        TotalAttemptsAll = 0
        TotalAttempts = 0    !This is the total number of attempts at creating a sym-allowed det (including successful ones)

        !Sz symmetry could be put in here to make it more efficient
        !(would be a little fiddly for OS systems though)
        FullSpace = choose_i64(NEl, iExcitLevTest)    !Pick 4 holes
        FullSpace = FullSpace * choose_i64(nBasis - NEl + iExcitLevTest, iExcitLevTest) !Pick total unoccupied space

        !Calculate excitation level bias due to the way the determinants are constructed.
        do i = 0, iExcitLevTest
            ExcitLevBias(i) = choose_i64(NEl - i, iExcitLevTest - i)
!             write(stdout,*) ExcitLevBias(i)
        end do

        write(IUNIT, *) "Size of excitation level neglecting all symmetry: " &
            , FullSpace

        CALL neci_flush(IUNIT)

        IF (iProcIndex == 0) THEN
            OPEN (14, file="TruncSpaceMCStats", status='unknown', &
                  form='formatted')
        end if

        ! With MerTwistRan the default seed was being used.
        ! dSFMT does not initialise itself if not already initialised.
        call dSFMT_init(5489)

        do i = 1, int(CalcDetCycles, sizeof_int)

            !Create a random determinant up to excitation level iExcitLevTest from FDetiLut
            !Returns det (iLut) and its excitation level, ExcitLevel, and the number of attempts
            !needed to generate the symmetry allowed determinant.
            CALL CreateRandomExcitLevDet(iExcitLevTest, FDet, FDetiLut, iLut, ExcitLev, Attempts)
            TotalAttempts = TotalAttempts + Attempts
            Accept = Accept + 1

!Add to correct bin for the excitation level
            ExcitBin(ExcitLev) = ExcitBin(ExcitLev) + 1

            IF (mod(i, int(CalcDetPrint, sizeof_int)) == 0) THEN
                !Write out statistics
                call MPIReduce(Accept, MPI_SUM, AcceptAll)
                call MPIReduce(TotalAttempts, MPI_SUM, TotalAttemptsAll)
                call MPIReduce(ExcitBin(0:iExcitLevTest), MPI_SUM, ExcitBinAll(0:iExcitLevTest))

                SymSpace = 0.0_dp
                Frac = REAL(AcceptAll, dp) / REAL(TotalAttemptsAll, dp)  !Fraction of the 'full' space which is symmetry allowed
                do j = 0, iExcitLevTest
!                     write(stdout,*) REAL(ExcitBinAll(j),dp),REAL(AcceptAll,dp),Frac,FullSpace,ExcitLevBias(j)
                    SizeLevel(j) = ((REAL(ExcitBinAll(j), dp) / REAL(AcceptAll, dp)) * Frac * FullSpace) / ExcitLevBias(j)
                    SymSpace = SymSpace + SizeLevel(j)
                end do
                IF (iProcIndex == 0) THEN
                    write(14, "(2I16,2G35.15)", advance='no') i, AcceptAll, Frac, SymSpace
                    do j = 0, iExcitLevTest
                        write(14, "(F30.5)", advance='no') SizeLevel(j)
                    end do
                    write(14, "(A)") ""
                end if

                AcceptAll = 0
                ExcitBinAll(0:iExcitLevTest) = 0
                TotalAttemptsAll = 0

                call ChangeVars(tDummy, tDummy2)
                IF (tSoftExitFound) EXIT

            end if

        end do

        call MPIReduce(Accept, MPI_SUM, AcceptAll)
        call MPIReduce(TotalAttempts, MPI_SUM, TotalAttemptsAll)
        call MPIReduce(ExcitBin(0:iExcitLevTest), MPI_SUM, ExcitBinAll(0:iExcitLevTest))

        SymSpace = 0.0_dp
        Frac = REAL(AcceptAll, dp) / REAL(TotalAttemptsAll, dp)  !Fraction of the 'full' space which is symmetry allowed
        do j = 0, iExcitLevTest
!             write(stdout,*) REAL(ExcitBinAll(j),dp),REAL(AcceptAll,dp),Frac,FullSpace,ExcitLevBias(j)
            SizeLevel(j) = ((REAL(ExcitBinAll(j), dp) / REAL(AcceptAll, dp)) * Frac * FullSpace) / ExcitLevBias(j)
            SymSpace = SymSpace + SizeLevel(j)
        end do
        IF (iProcIndex == 0) THEN
            write(14, "(2I16,2G35.15)", advance='no') i, AcceptAll, Frac, SymSpace
            do j = 0, iExcitLevTest
                write(14, "(F30.5)", advance='no') SizeLevel(j)
            end do
            write(14, "(A)") ""
            close(14)
        end if

        write(IUNIT, *) "MC size of truncated space: ", SymSpace
        write(IUNIT, *) "Individual excitation level contributions: "
        do j = 0, iExcitLevTest
            write(IUNIT, "(I5,F30.5)") j, SizeLevel(j)
        end do
        CALL neci_flush(IUNIT)

    END SUBROUTINE FindSymMCSizeExcitLevel

!This routine calls CreateRandomExcitLevDet, but it returns an *unbiased* determinant from the excitation
!levels 0 -> iExcitLevTest. It does this by rejecting determinants, so that the resultant excitations are
!unbiased.
    SUBROUTINE CreateRandomExcitLevDetUnbias(iExcitLevTest, FDet, FDetiLut, iLut, ExcitLev, Attempts)
        use SystemData, only: nEl
        use bit_rep_data, only: NIfTot
        use dSFMT_interface
        INTEGER :: iExcitLevTest, FDet(NEl), ExcitLev, Attempts
        INTEGER(n_int) :: FDetiLut(0:NIfTot), iLut(0:NIfTot)
        real(dp) :: pAcc, r

        do while (.true.)

            call CreateRandomExcitLevDet(iExcitLevTest, FDet, FDetiLut, iLut, ExcitLev, Attempts)

            IF (ExcitLev == iExcitLevTest) then
                RETURN   !Prob of accepting = 1
            else
                pAcc = 1.0_dp / (choose_i64(NEl - ExcitLev, iExcitLevTest - ExcitLev))
                r = genrand_real2_dSFMT()
                if (r <= pAcc) exit
            end if

        end do

    END SUBROUTINE CreateRandomExcitLevDetUnbias

    subroutine create_rand_heisenberg_det(ilut)

        use bit_rep_data, only: NIfTot
        use dSFMT_interface, only: genrand_real2_dSFMT
        use SystemData, only: nbasis, lms

        integer(n_int), intent(out) :: ilut(0:NIfTot)
        integer :: i, nsites, n_up, n_flipped, site_ind, bit_ind, elem
        integer :: beta_ind, alpha_ind
        real(dp) :: r
        logical :: is_alpha, is_beta

        ! Start from all spins down and pick random spins to flip up.
        ilut = 0_n_int
        n_flipped = 0

        nsites = nbasis / 2
        n_up = (lms + nsites) / 2

        do
            r = genrand_real2_dSFMT()
            site_ind = int(r * nsites) + 1
            bit_ind = 2 * site_ind
            elem = (bit_ind - 1) / bits_n_int
            ! If this spin has already been flipped.
            if (btest(ilut(elem), mod(bit_ind - 1, bits_n_int))) cycle
            ! Flip the spin up.
            ilut(elem) = ibset(ilut(elem), mod(bit_ind - 1, bits_n_int))
            n_flipped = n_flipped + 1
            if (n_flipped == n_up) exit
        end do

        do i = 1, nsites
            ! If both alpha and beta bits are down, set the beta one up, to
            ! represent a down spin.
            beta_ind = 2 * i - 1
            alpha_ind = 2 * i
            elem = (alpha_ind - 1) / bits_n_int
            is_alpha = btest(ilut(elem), mod(alpha_ind - 1, bits_n_int))
            is_beta = btest(ilut(elem), mod(beta_ind - 1, bits_n_int))
            if ((.not. is_alpha) .and. (.not. is_beta)) then
                ilut(elem) = ibset(ilut(elem), mod(beta_ind - 1, bits_n_int))
            end if
        end do

    end subroutine create_rand_heisenberg_det

!Create stochastically a random symmetry-allowed determinant from excitation level iExcitLevTest or less, with respect to i
!the bit representation determinant iLutFDet, which is passed in.
!The determinant is returned in bit-form in iLut, and the excitation level of the determinant in ExcitLev.
!This routine will only work when TwoCycleSymGens is true.
!***WARNING*** The determinants created in this way are *NOT* uniform.
!Determinants of a certain excitation level are more likely to be generated than others.
!The bias towards a given determinant is given by:
!(NEl-ExcitLev) Choose (iExcitLevTest-ExcitLev)
    SUBROUTINE CreateRandomExcitLevDet(iExcitLevTest, FDet, FDetiLut, iLut, ExcitLev, Attempts)
        use SystemData, only: nEl, G1, nBasis
        use SystemData, only: tUEG, tHPHF, tHub
        use SystemData, only: tFixLz
        use dSFMT_interface
        use bit_rep_data, only: NIfTot
        use DetBitOps, only: IsAllowedHPHF
        integer, intent(in) :: iExcitLevTest, FDet(NEl)
        integer, intent(out) :: ExcitLev, Attempts
        integer(n_int), intent(out) :: iLut(0:NIfTot)
        integer(n_int), intent(in) :: FDetiLut(0:NIfTot)
        logical :: tSymAllowedDet, tNotAllowed
        integer :: TotalSym, TotalMom, TotalMs, Momx, Momy, Momz, j, Elec, Orb, Hole
        real(dp) :: r

        Attempts = 0 !Count the number of attempts needed to generate the sym-allowed determinant.

        tSymAllowedDet = .false.
        do while (.not. tSymAllowedDet)

            TotalSym = 0
            TotalMom = 0
            TotalMs = 0
            Momx = 0
            Momy = 0
            Momz = 0
            ExcitLev = 0
            iLut(:) = FDetiLut(:)

            !Create random determinant
            !Loop over holes in occupied space
            do j = 1, iExcitLevTest

                tNotAllowed = .true.
                do while (tNotAllowed)  !Loop until we have created an allowed hole.
                    r = genrand_real2_dSFMT()
                    Elec = int(NEl * r) + 1
                    Orb = FDet(Elec)

                    !Electron picked must not be one which has been picked before
                    !i.e. it must be occupied in iLut
                    if (IsOcc(iLut, Orb)) then
                        !Clear orbital to indicate it is gone.
                        clr_orb(iLut, Orb)
                        tNotAllowed = .false.
                        !Deal with totting up the symmetry for the now unocc orbital
                        TotalSym = IEOR(TotalSym, INT((G1(Orb)%Sym%S), sizeof_int))
                        TotalMom = TotalMom + G1(Orb)%Ml
                        TotalMs = TotalMs + G1(Orb)%Ms
                        IF (tUEG .or. tHub) THEN
                            Momx = Momx + G1(Orb)%k(1)
                            Momy = Momy + G1(Orb)%k(2)
                            Momz = Momz + G1(Orb)%k(3)
                        end if
                    end if
                end do
            end do

            !Loop over electrons in the unoccupied space
            do j = 1, iExcitLevTest

                tNotAllowed = .true.
                do while (tNotAllowed)  !Loop until we have created an allowed electron
                    r = genrand_real2_dSFMT()
                    Hole = int(nBasis * r) + 1

                    if (IsNotOcc(iLut, Hole)) then
                        !Set orbital to indicate it is now occupied
                        set_orb(iLut, Hole)
                        tNotAllowed = .false.
                        !Increase excitation level
                        if (IsNotOcc(FDetiLut, Hole)) ExcitLev = ExcitLev + 1
                        !Deal with totting up the symmetry for the now occ orbital
                        TotalSym = IEOR(TotalSym, INT((G1(Hole)%Sym%S), sizeof_int))
                        TotalMom = TotalMom - G1(Hole)%Ml
                        TotalMs = TotalMs - G1(Hole)%Ms
                        if (tUEG .or. tHub) then
                            Momx = Momx - G1(Hole)%k(1)
                            Momy = Momy - G1(Hole)%k(2)
                            Momz = Momz - G1(Hole)%k(3)
                        end if
                    end if
                end do
            end do

            if ((TotalSym == 0) .and. (TotalMom == 0) .and. (Momx == 0) .and. (Momy == 0) .and. (Momz == 0) .and. (TotalMs == 0)) then
                !Created determinant is symmetry allowed.
                if (tHPHF) then
                    if (IsAllowedHPHF(iLut)) tSymAllowedDet = .true.
                else
                    tSymAllowedDet = .true.
                end if
            end if

            Attempts = Attempts + 1

        end do

    END SUBROUTINE CreateRandomExcitLevDet

    subroutine create_rand_det_no_sym(ilut)

        ! Create a determinant with no symmetry imposed whatsoever, apart from using nel electrons.
        ! This routine simply picks nel orbitals uniformly and returns a determinant with these occupied.

        use bit_rep_data, only: NIfTot
        USE dSFMT_interface, only: genrand_real2_dSFMT
        use SystemData, only: nbasis, nel

        integer(n_int), intent(out) :: ilut(0:NIfTot)
        integer :: i, n_occ, orb_ind, elem
        real(dp) :: r

        ! Start from all orbitals unoccupied.
        ilut = 0_n_int
        n_occ = 0

        do
            r = genrand_real2_dSFMT()
            orb_ind = int(r * nbasis) + 1
            elem = (orb_ind - 1) / bits_n_int
            ! If this orbital has already been occupied then choose another.
            if (btest(ilut(elem), mod(orb_ind - 1, bits_n_int))) cycle
            ! Occupy the orbital.
            ilut(elem) = ibset(ilut(elem), mod(orb_ind - 1, bits_n_int))
            n_occ = n_occ + 1
            if (n_occ == nel) exit
        end do

    end subroutine create_rand_det_no_sym

!This routine *stochastically* finds the size of the determinant space. For certain symmetries, its hard to find the
!allowed size of the determinant space. However, it can be simply found using a MC technique.
    SUBROUTINE FindSymMCSizeofSpace(IUNIT)
        use SymData, only: TwoCycleSymGens
        use SystemData, only: nEl, G1, nBasis, nOccAlpha, nOccBeta, tMolpro
        use SystemData, only: tUEG, tHPHF, tHub, tKPntSym, Symmetry
        use SystemData, only: CalcDetCycles, CalcDetPrint, tFixLz
        use CalcData, only: tTruncNOpen, trunc_nopen_max, &
                            t_trunc_nopen_diff, trunc_nopen_diff
        use DeterminantData, only: FDet
        use DetCalcData, only: ICILevel
        use dSFMT_interface
        use soft_exit, only: ChangeVars, tSoftExitFound
        use Parallel_neci
        use DetBitops, only: EncodeBitDet, IsAllowedHPHF, count_open_orbs
        use bit_rep_data, only: NIfTot
        use sym_mod, only: SymProd
        use fcimcdata, only: ilutRef
        IMPLICIT NONE
        INTEGER :: IUNIT, j, SpatOrbs, FDetMom, ExcitLev
        INTEGER(KIND=n_int) :: FDetiLut(0:NIfTot), iLut(0:NIfTot)
        INTEGER :: FDetSym, TotalSym, TotalMom, alpha, beta, ierr, Momx, Momy
        INTEGER :: Momz, nopenorbs, Space_unit
        integer(int64) :: Accept, AcceptAll, i
        integer(int64) :: ExcitBin(0:NEl), ExcitBinAll(0:NEl)
        real(dp) :: FullSpace, r, Frac
        real(dp) :: SizeLevel(0:NEl)
        LOGICAL :: truncate_space, tDummy, tDummy2
        LOGICAL :: tNotAllowed, tAcc
        type(Symmetry) :: FDetKPntMom, KPntMom

        write(IUNIT, *) "Calculating exact size of symmetry-allowed " &
            //"determinant space using MC..."
        write(IUNIT, *) CalcDetCycles, " MC cycles will be used, and " &
            //"statistics printed out every ", CalcDetPrint, " cycles."
        FDetSym = 0
        FDetMom = 0
        ExcitBin(:) = 0
        ExcitBinAll(:) = 0
        Momx = 0
        Momy = 0
        Momz = 0
        FDetKPntMom%S = 0

        do i = 1, NEl
            IF (tFixLz) FDetMom = FDetMom + G1(FDet(i))%Ml
            if (tKPntSym) FDetKPntMom = Symprod(FDetKPntMom, G1(FDet(i))%Sym)
            IF (tUEG .or. tHub) THEN
                Momx = Momx + G1(FDet(i))%k(1)
                Momy = Momy + G1(FDet(i))%k(2)
                Momz = Momz + G1(FDet(i))%k(3)
            else
                FDetSym = IEOR(FDetSym, INT(G1(FDet(i))%Sym%S, sizeof_int))
            end if
        end do

        IF (.not. ((Momx == 0) .and. (Momy == 0) .and. (Momz == 0))) THEN
            write(stdout, *) "Momx: ", Momx, "Momy: ", Momy, "Momz: ", Momz
            call stop_all("FindSymMCSizeofSpace", "Cannot calculate MC size of space with non-zero momentum")
        end if

        IF (ICILevel > 0) THEN
            truncate_space = .true.
        ELSE
            truncate_space = .false.
        end if

        CALL EncodeBitDet(FDet, FDetiLut)

        write(IUNIT, *) "Symmetry of HF determinant is: ", FDetSym
        IF (tFixLz) THEN
            write(IUNIT, *) "Imposing momentum sym on size calculation"
            write(IUNIT, *) "Momentum of HF determinant is: ", FDetMom
        end if
        if (tKPntSym) then
            write(iunit, *) "Using k-point symmetry."
            write(iunit, *) "K-point sym of HF determinant is: ", FDetKPntMom%S
        end if
        IF (tHPHF) THEN
            write(stdout, *) "Imposing time-reversal symmetry (HPHF) on " &
                //"size of space calculation"
        end if

        SpatOrbs = nBasis / 2

        Accept = 0

        FullSpace = choose_i64(SpatOrbs, nOccAlpha)
        FullSpace = FullSpace * choose_i64(SpatOrbs, nOccBeta)

        write(IUNIT, *) "Size of space neglecting all but Sz symmetry: " &
            , FullSpace

        CALL neci_flush(IUNIT)

        IF (iProcIndex == 0) THEN
            Space_unit = get_free_unit()
            OPEN (Space_unit, file="SpaceMCStats", status='unknown', &
                  form='formatted')
        end if

        ! With MerTwistRan the default seed was being used.
        ! dSFMT does not initialise itself if not already initialised.
        call dSFMT_init(5489)

        do i = 1, CalcDetCycles

            KPntMom%S = 0
            TotalSym = 0
            TotalMom = 0
            ExcitLev = 0
            Momx = 0
            Momy = 0
            Momz = 0
            iLut(:) = 0

            !Create random determinant (Correct Sz symmetry)
            !Loop over alpha electrons
            do j = 1, nOccAlpha

                tNotAllowed = .true.
                do while (tNotAllowed)

                    r = genrand_real2_dSFMT()
                    alpha = 2 * (INT(SpatOrbs * r) + 1)
                    IF (.not. BTEST(iLut((alpha - 1) / bits_n_int), mod((alpha - 1), bits_n_int))) THEN
                        !Has *not* been picked before
                        iLut((alpha - 1) / bits_n_int) = &
                            IBSET(iLut((alpha - 1) / bits_n_int), mod(alpha - 1, bits_n_int))
                        tNotAllowed = .false.
                    end if
                end do

                IF (tFixLz) THEN
                    TotalMom = TotalMom + G1(alpha)%Ml
                end if
                IF (tUEG .or. tHub) THEN
                    Momx = Momx + G1(alpha)%k(1)
                    Momy = Momy + G1(alpha)%k(2)
                    Momz = Momz + G1(alpha)%k(3)
                else if (tKPntSym) then
                    KPntMom = Symprod(KPntMom, G1(alpha)%Sym)
                else
                    TotalSym = IEOR(TotalSym, INT((G1(alpha)%Sym%S), sizeof_int))
                end if
                IF (.not. BTEST(FDetiLut((alpha - 1) / bits_n_int), mod((alpha - 1), bits_n_int))) THEN
                    !orbital chosen is *not* in the reference determinant
                    ExcitLev = ExcitLev + 1
                end if

            end do

            !Loop over beta electrons
            do j = 1, nOccBeta

                tNotAllowed = .true.
                do while (tNotAllowed)
                    r = genrand_real2_dSFMT()
                    beta = 2 * (INT(SpatOrbs * r) + 1) - 1
                    IF (.not. BTEST(iLut((beta - 1) / bits_n_int), mod((beta - 1), bits_n_int))) THEN
                        !Has *not* been picked before
                        iLut((beta - 1) / bits_n_int) = &
                            IBSET(iLut((beta - 1) / bits_n_int), mod(beta - 1, bits_n_int))
                        tNotAllowed = .false.
                    end if
                end do

                IF (tFixLz) THEN
                    TotalMom = TotalMom + G1(beta)%Ml
                end if
                IF (tUEG .or. tHub) THEN
                    Momx = Momx + G1(beta)%k(1)
                    Momy = Momy + G1(beta)%k(2)
                    Momz = Momz + G1(beta)%k(3)
                else if (tKPntSym) then
                    KPntMom = Symprod(KPntMom, G1(beta)%Sym)
                ELSE
                    TotalSym = IEOR(TotalSym, INT((G1(beta)%Sym%S), sizeof_int))
                end if
                IF (.not. BTEST(FDetiLut((beta - 1) / bits_n_int), mod((beta - 1), bits_n_int))) THEN
                    !orbital chosen is *not* in the reference determinant
                    ExcitLev = ExcitLev + 1
                end if

            end do

            if (tTruncNOpen .or. t_trunc_nopen_diff) nOpenOrbs = count_open_orbs(iLut)

            tAcc = .true.
            if (TotalSym /= FDetSym) then
                tAcc = .false.
            end if

            if (tAcc .and. tHPHF) then
                if (.not. IsAllowedHPHF(iLut)) then
                    tAcc = .false.
                end if
            end if

            if (tAcc .and. tFixLz) then
                if (TotalMom /= FDetMom) then
                    tAcc = .false.
                end if
            end if

            if (tAcc .and. truncate_space) then
                if (ExcitLev > ICILevel) then
                    tAcc = .false.
                end if
            end if

            if (tAcc .and. (tUEG .or. tHub)) then
                if ((Momx /= 0) .or. (Momy /= 0) .or. (Momz /= 0)) then
                    tAcc = .false.
                end if
            end if

            if (tAcc .and. (tKPntSym)) then
                if (KPntMom%S /= FDetKPntMom%S) then
                    tAcc = .false.
                end if
            end if

            if (tAcc .and. tTruncNOpen) then
                if (nOpenOrbs > trunc_nopen_max) then
                    tAcc = .false.
                end if
            end if

            if (tAcc .and. t_trunc_nopen_diff) then
                if (abs(nOpenOrbs - count_open_orbs(ilutRef(:, 1))) > trunc_nopen_diff) then
                    tAcc = .false.
                end if
            end if

            if (tAcc) Accept = Accept + 1

            IF (tAcc) THEN
!Add to correct bin for the excitation level
                ExcitBin(ExcitLev) = ExcitBin(ExcitLev) + 1
            end if

            IF (mod(i, CalcDetPrint) == 0) THEN
                !Write out statistics
                call MPIAllReduce(Accept, MPI_SUM, AcceptAll)
                call MPIAllReduce(ExcitBin(0:NEl), MPI_SUM, ExcitBinAll(0:NEl))

                Frac = REAL(AcceptAll, dp) / REAL(i * nProcessors, dp)
                do j = 0, NEl
                    SizeLevel(j) = (REAL(ExcitBinAll(j), dp) / REAL(AcceptAll, dp)) * Frac * FullSpace
                end do
                IF (iProcIndex == 0) THEN
                    if (tMolpro) then
                        write(Space_unit, "(I16,G35.15)") i, Frac * FullSpace
                    else
                        write(Space_unit, "(2I16,2G35.15)", advance='no') i, AcceptAll, Frac, Frac * FullSpace
                        do j = 0, NEl
                            write(Space_unit, "(F30.5)", advance='no') SizeLevel(j)
                        end do
                        write(Space_unit, "(A)") ""
                    end if
                end if

                AcceptAll = 0
                ExcitBinAll(0:NEl) = 0

                CALL ChangeVars(tDummy, tDummy2)
                IF (tSoftExitFound) EXIT

            end if

        end do

        call MPIAllReduce(Accept, MPI_SUM, AcceptAll)
        call MPIAllReduce(ExcitBin(0:NEl), MPI_SUM, ExcitBinAll(0:NEl))

        Frac = REAL(AcceptAll, dp) / REAL(i * nProcessors, dp)
        do j = 0, NEl
            SizeLevel(j) = (REAL(ExcitBinAll(j), dp) / REAL(AcceptAll, dp)) * Frac * FullSpace
        end do

        IF (iProcIndex == 0) THEN
            if (tMolpro) then
                write(Space_unit, "(I16,G35.15)") i, Frac * FullSpace
            else
                write(Space_unit, "(2I16,2G35.15)", advance='no') i, AcceptAll, Frac, Frac * FullSpace
                do j = 0, NEl
                    write(Space_unit, "(F30.5)", advance='no') SizeLevel(j)
                end do
                write(Space_unit, "(A)") ""
            end if
            close(Space_unit)
        end if

        write(IUNIT, *) "MC size of space: ", Frac * FullSpace
        write(IUNIT, *) "Individual excitation level contributions: "
        do j = 0, NEl
            write(IUNIT, "(I5,F30.5)") j, SizeLevel(j)
        end do
        CALL neci_flush(IUNIT)

    END SUBROUTINE FindSymMCSizeofSpace

!!This routine finds the size of the determinant space in terms, including all symmetry allowed determinants.
!!This is written to IUNIT. This is only available for molecular (i.e. abelian) systems with a maximum of eigth irreps.
!!This is done in a very crude way. Feel free to optimise it!
    SUBROUTINE FindSymSizeofSpace(IUNIT)
        use SymData, only: TwoCycleSymGens
        use SystemData, only: nEl, G1, nBasis, nOccAlpha, nOccBeta
        use DeterminantData, only: FDet
        IMPLICIT NONE
        INTEGER :: ClassCounts(2, 0:7)
        INTEGER :: Lima(0:7), Limb(0:7), a0, a1, a2, a3, a4, a5, a6, a7, NAlph
        INTEGER :: b0, b1, b2, b3, b4, b5, b6, b7, NBet, i, IUNIT, OverallSym
        INTEGER :: FDetSym
        real(dp) :: Space, SpaceGrow
        LOGICAL :: Sym(0:7)

        write(IUNIT, *) "Calculating exact size of symmetry-allowed determinant space..."
        FDetSym = 0
        do i = 1, NEl
            FDetSym = IEOR(FDetSym, INT(G1(FDet(i))%Sym%S, sizeof_int))
        end do
        write(stdout, *) "Symmetry of HF determinant is: ", FDetSym
        CALL neci_flush(IUNIT)
        ClassCounts(:, :) = 0
!First, we need to find the number of spatial orbitals in each symmetry irrep.
        do i = 1, nBasis, 1
            IF (G1(i)%Ms == 1) THEN
                ClassCounts(1, INT(G1(i)%Sym%S, sizeof_int)) = &
                    ClassCounts(1, INT(G1(i)%Sym%S, sizeof_int)) + 1
            ELSE

                ClassCounts(2, INT(G1(i)%Sym%S, sizeof_int)) = &
                    ClassCounts(2, INT(G1(i)%Sym%S, sizeof_int)) + 1
            end if
        end do
        do i = 0, 7
            IF (mod((ClassCounts(1, i) + ClassCounts(2, i)), 2) /= 0) THEN
                write(stdout, *) 'WARNING: Different number of symmetries between the alpha and beta orbitals.'
            end if
        end do

        Lima(0) = min(nOccAlpha, ClassCounts(1, 0))
        Limb(0) = min(nOccBeta, ClassCounts(2, 0))
        Lima(1) = min(nOccAlpha, ClassCounts(1, 1))
        Limb(1) = min(nOccBeta, ClassCounts(2, 1))
        Lima(2) = min(nOccAlpha, ClassCounts(1, 2))
        Limb(2) = min(nOccBeta, ClassCounts(2, 2))
        Lima(3) = min(nOccAlpha, ClassCounts(1, 3))
        Limb(3) = min(nOccBeta, ClassCounts(2, 3))
        Lima(4) = min(nOccAlpha, ClassCounts(1, 4))
        Limb(4) = min(nOccBeta, ClassCounts(2, 4))
        Lima(5) = min(nOccAlpha, ClassCounts(1, 5))
        Limb(5) = min(nOccBeta, ClassCounts(2, 5))
        Lima(6) = min(nOccAlpha, ClassCounts(1, 6))
        Limb(6) = min(nOccBeta, ClassCounts(2, 6))
        Lima(7) = min(nOccAlpha, ClassCounts(1, 7))
        Limb(7) = min(nOccBeta, ClassCounts(2, 7))
        Space = 0.0_dp


!Loop over each irrep twice, once for alpha electrons and once for beta.
        do a0 = 0, Lima(0)
        do b0 = 0, Limb(0)
            IF (mod(a0 + b0, 2) == 1) THEN
                Sym(0) = .true.
            ELSE
                Sym(0) = .false.
            end if
            do a1 = 0, Lima(1)
            do b1 = 0, Limb(1)
                IF (mod(a1 + b1, 2) == 1) THEN
                    Sym(1) = .true.
                ELSE
                    Sym(1) = .false.
                end if
                do a2 = 0, Lima(2)
                do b2 = 0, Limb(2)
                    IF (mod(a2 + b2, 2) == 1) THEN
                        Sym(2) = .true.
                    ELSE
                        Sym(2) = .false.
                    end if
                    do a3 = 0, Lima(3)
                    do b3 = 0, Limb(3)
                        IF (mod(a3 + b3, 2) == 1) THEN
                            Sym(3) = .true.
                        ELSE
                            Sym(3) = .false.
                        end if
                        do a4 = 0, Lima(4)
                        do b4 = 0, Limb(4)
                            IF (mod(a4 + b4, 2) == 1) THEN
                                Sym(4) = .true.
                            ELSE
                                Sym(4) = .false.
                            end if
                            do a5 = 0, Lima(5)
                            do b5 = 0, Limb(5)
                                IF (mod(a5 + b5, 2) == 1) THEN
                                    Sym(5) = .true.
                                ELSE
                                    Sym(5) = .false.
                                end if
                                do a6 = 0, Lima(6)
                                do b6 = 0, Limb(6)
                                    IF (mod(a6 + b6, 2) == 1) THEN
                                        Sym(6) = .true.
                                    ELSE
                                        Sym(6) = .false.
                                    end if
                                    do a7 = 0, Lima(7)
                                    do b7 = 0, Limb(7)
                                        IF (mod(a7 + b7, 2) == 1) THEN
                                            Sym(7) = .true.
                                        ELSE
                                            Sym(7) = .false.
                                        end if

                                        OverallSym = 0
                                        do i = 0, 7
                                            IF (Sym(i)) THEN
                                                OverallSym = IEOR(OverallSym, i)
                                            end if
                                        end do
                                        IF (OverallSym == FDetSym) THEN
                                            NAlph = a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7
                                            NBet = b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7

                                            IF ((NAlph == NOccAlpha) .and. (NBet == NOccBeta)) THEN

                                                SpaceGrow = 1.0_dp
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(1, 0), a0)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(2, 0), b0)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(1, 1), a1)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(2, 1), b1)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(1, 2), a2)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(2, 2), b2)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(1, 3), a3)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(2, 3), b3)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(1, 4), a4)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(2, 4), b4)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(1, 5), a5)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(2, 5), b5)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(1, 6), a6)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(2, 6), b6)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(1, 7), a7)
                                                SpaceGrow = SpaceGrow * choose_i64(ClassCounts(2, 7), b7)
                                                Space = Space + SpaceGrow
                                            end if
                                        end if

                                    end do
                                    end do
                                end do
                                end do
                            end do
                            end do
                        end do
                        end do
                    end do
                    end do
                end do
                end do
            end do
            end do
        end do
        end do

        write(IUNIT, "(A,G25.16)") " *EXACT* size of symmetry allowed space of determinants is: ", Space
        CALL neci_flush(IUNIT)

    END SUBROUTINE FindSymSizeofSpace

!This routine finds the size of the determinant space in terms, including all symmetry allowed determinants.
!This is written to IUNIT. This is only available for molecular (i.e. abelian) systems with a maximum of eigth irreps.
!This is done in a very crude way. Feel free to optimise it!
    SUBROUTINE FindSymSizeofTruncSpace(IUNIT)
        use SymData, only: TwoCycleSymGens
        use SystemData, only: nEl, G1, nBasis, nOccAlpha, nOccBeta, Brr
        use DeterminantData, only: FDet
        use DetCalcData, only: ICILevel
        IMPLICIT NONE
        INTEGER :: ClassCountsOcc(0:7)
        INTEGER :: ClassCountsVirt(0:7), NAlphOcc, NAlphVirt
        INTEGER :: ClassCountsOccMax(0:7), ClassCountsVirtMax(0:7)
        INTEGER :: LimaOcc(0:7), LimbOcc(0:7), LimaVirt(0:7)
        INTEGER :: LimbVirt(0:7)
        INTEGER :: a0o, a1o, a2o, a3o, a4o, a5o, a6o, a7o
        INTEGER :: a0v, a1v, a2v, a3v, a4v, a5v, a6v, a7v
        INTEGER :: b0o, b1o, b2o, b3o, b4o, b5o, b6o, b7o, OverallSym
        INTEGER :: b0v, b1v, b2v, b3v, b4v, b5v, b6v, b7v, NBetOcc, i, IUNIT
        INTEGER :: FDetSym, NBetVirt
        real(dp) :: Space, SpaceGrow
        LOGICAL :: Sym(0:7)
        character(*), parameter :: this_routine = 'FindSymSizeofTruncSpace'

        IF (.not. TwoCycleSymGens) THEN
            write(IUNIT, *) "Only for molecular abelian symmetry " &
                //" calculations can the exact size of the determinant " &
                //" space be calculated currently..."
            write(IUNIT, *) "Skipping size of space calculation..."
            RETURN
        end if

        write(IUNIT, *) "Calculating exact size of symmetry-allowed " &
            //"determinant space..."
        FDetSym = 0
        do i = 1, NEl
            FDetSym = IEOR(FDetSym, INT(G1(FDet(i))%Sym%S, sizeof_int))
        end do
        write(stdout, *) "Symmetry of HF determinant is: ", FDetSym
        CALL neci_flush(IUNIT)
        ClassCountsOcc(:) = 0
        ClassCountsVirt(:) = 0
!First, we need to find the number of spatial orbitals in each symmetry irrep.
!We need to separate this into occupied and virtual.
        do i = 1, NEl, 1
            ClassCountsOcc(INT(G1(BRR(i))%Sym%S, sizeof_int)) = &
                ClassCountsOcc(INT(G1(BRR(i))%Sym%S, sizeof_int)) + 1
        end do

        do i = NEL + 1, nBasis, 1
            ClassCountsVirt(INT(G1(BRR(i))%Sym%S, sizeof_int)) = &
                ClassCountsVirt(INT(G1(BRR(i))%Sym%S, sizeof_int)) + 1
        end do

!These are still in spin orbitals, so check there are multiple of 2 values in
!each symmetry irrep and then divide by two because we deal with alpha and beta separately.
        do i = 0, 7
        IF (mod((ClassCountsOcc(i) + ClassCountsVirt(i)), 2) /= 0) THEN
            call stop_all(this_routine, 'Error counting determinants')
        end if
        ClassCountsOccMax(i) = CEILING(REAL(ClassCountsOcc(i), dp) / 2.0_dp)
        ClassCountsVirtMax(i) = CEILING(REAL(ClassCountsVirt(i), dp) / 2.0_dp)
        ClassCountsOcc(i) = FLOOR(REAL(ClassCountsOcc(i), dp) / 2.0_dp)
        ClassCountsVirt(i) = FLOOR(REAL(ClassCountsVirt(i), dp) / 2.0_dp)

        end do

        IF (nOccAlpha > nOccBeta) THEN
            do i = 0, 7
                LimaOcc(i) = min(nOccAlpha, ClassCountsOccMax(i))
                LimbOcc(i) = min(nOccBeta, ClassCountsOcc(i))
                LimaVirt(i) = min(ICILevel, ClassCountsVirtMax(i))
                LimbVirt(i) = min(ICILevel, ClassCountsVirt(i))
            end do
        ELSE
            do i = 0, 7
                LimaOcc(i) = min(nOccAlpha, ClassCountsOcc(i))
                LimbOcc(i) = min(nOccBeta, ClassCountsOccMax(i))
                LimaVirt(i) = min(ICILevel, ClassCountsVirt(i))
                LimbVirt(i) = min(ICILevel, ClassCountsVirtMax(i))
            end do
        end if

        Space = 0.0_dp

!Loop over each irrep twice, once for alpha electrons and once for beta.
!a0 is the number of alpha electrons in symmetry 0.
!b0 is the number of beta electrons in symmetry 0.
        do a0o = 0, LimaOcc(0)
        do b0o = 0, LimbOcc(0)
        do a0v = 0, LimaVirt(0)
        do b0v = 0, LimbVirt(0)
            IF (mod(a0o + b0o + a0v + b0v, 2) == 1) THEN
                Sym(0) = .true.
            ELSE
                Sym(0) = .false.
            end if
            do a1o = 0, LimaOcc(1)
            do b1o = 0, LimbOcc(1)
            do a1v = 0, LimaVirt(1)
            do b1v = 0, LimbVirt(1)
                IF (mod(a1o + b1o + a1v + b1v, 2) == 1) THEN
                    Sym(1) = .true.
                ELSE
                    Sym(1) = .false.
                end if
                do a2o = 0, LimaOcc(2)
                do b2o = 0, LimbOcc(2)
                do a2v = 0, LimaVirt(2)
                do b2v = 0, LimbVirt(2)
                    IF (mod(a2o + b2o + a2v + b2v, 2) == 1) THEN
                        Sym(2) = .true.
                    ELSE
                        Sym(2) = .false.
                    end if

                    do a3o = 0, LimaOcc(3)
                    do b3o = 0, LimbOcc(3)
                    do a3v = 0, LimaVirt(3)
                    do b3v = 0, LimbVirt(3)
                        IF (mod(a3o + b3o + a3v + b3v, 2) == 1) THEN
                            Sym(3) = .true.
                        ELSE
                            Sym(3) = .false.
                        end if
                        do a4o = 0, LimaOcc(4)
                        do b4o = 0, LimbOcc(4)
                        do a4v = 0, LimaVirt(4)
                        do b4v = 0, LimbVirt(4)
                            IF (mod(a4o + b4o + a4v + b4v, 2) == 1) THEN
                                Sym(4) = .true.
                            ELSE
                                Sym(4) = .false.
                            end if
                            do a5o = 0, LimaOcc(5)
                            do b5o = 0, LimbOcc(5)
                            do a5v = 0, LimaVirt(5)
                            do b5v = 0, LimbVirt(5)
                                IF (mod(a5o + b5o + a5v + b5v, 2) == 1) THEN
                                    Sym(5) = .true.
                                ELSE
                                    Sym(5) = .false.
                                end if
                                do a6o = 0, LimaOcc(6)
                                do b6o = 0, LimbOcc(6)
                                do a6v = 0, LimaVirt(6)
                                do b6v = 0, LimbVirt(6)
                                    IF (mod(a6o + b6o + a6v + b6v, 2) == 1) THEN
                                        Sym(6) = .true.
                                    ELSE
                                        Sym(6) = .false.
                                    end if
                                    do a7o = 0, LimaOcc(7)
                                    do b7o = 0, LimbOcc(7)
                                    do a7v = 0, LimaVirt(7)
                                    do b7v = 0, LimbVirt(7)
                                        IF (mod(a7o + b7o + a7v + b7v, 2) == 1) THEN
                                            Sym(7) = .true.
                                        ELSE
                                            Sym(7) = .false.
                                        end if

                                        OverallSym = 0
                                        do i = 0, 7
                                            IF (Sym(i)) THEN
                                                OverallSym = IEOR(OverallSym, i)
                                            end if
                                        end do
                                        IF (OverallSym == FDetSym) THEN
                                            NAlphOcc = a0o + a1o + a2o + a3o + a4o + a5o + a6o + a7o
                                            NBetOcc = b0o + b1o + b2o + b3o + b4o + b5o + b6o + b7o
                                            NAlphVirt = a0v + a1v + a2v + a3v + a4v + a5v + a6v + a7v
                                            NBetVirt = b0v + b1v + b2v + b3v + b4v + b5v + b6v + b7v

                                            IF (((NAlphOcc + NAlphVirt) == NOccAlpha) .and. ((NBetOcc + NBetVirt) == NOccBeta)) THEN
                                            IF ((NAlphVirt + NBetVirt) <= ICILevel) THEN

                                                IF (nOccAlpha > nOccBeta) THEN

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(0), a0o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(0), b0o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(0), a0v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(0), b0v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(1), a1o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(1), b1o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(1), a1v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(1), b1v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(2), a2o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(2), b2o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(2), a2v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(2), b2v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(3), a3o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(3), b3o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(3), a3v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(3), b3v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(4), a4o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(4), b4o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(4), a4v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(4), b4v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(5), a5o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(5), b5o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(5), a5v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(5), b5v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(6), a6o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(6), b6o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(6), a6v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(6), b6v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(7), a7o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(7), b7o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(7), a7v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(7), b7v)

                                                    Space = Space + SpaceGrow

                                                ELSE

                                                    SpaceGrow = 1.0_dp
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(0), a0o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(0), b0o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(0), a0v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(0), b0v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(1), a1o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(1), b1o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(1), a1v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(1), b1v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(2), a2o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(2), b2o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(2), a2v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(2), b2v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(3), a3o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(3), b3o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(3), a3v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(3), b3v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(4), a4o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(4), b4o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(4), a4v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(4), b4v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(5), a5o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(5), b5o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(5), a5v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(5), b5v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(6), a6o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(6), b6o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(6), a6v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(6), b6v)

                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOcc(7), a7o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsOccMax(7), b7o)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirt(7), a7v)
                                                    SpaceGrow = SpaceGrow * choose_i64(ClassCountsVirtMax(7), b7v)

                                                    Space = Space + SpaceGrow
                                                end if
                                            end if
                                            end if
                                        end if

                                    end do
                                    end do
                                    end do
                                    end do
                                end do
                                end do
                                end do
                                end do
                            end do
                            end do
                            end do
                            end do
                        end do
                        end do
                        end do
                        end do
                    end do
                    end do
                    end do
                    end do
                end do
                end do
                end do
                end do
            end do
            end do
            end do
            end do
        end do
        end do
        end do
        end do

        write(IUNIT, "(A,G25.16)") " *EXACT* size of symmetry allowed " &
            //"space of determinants is: ", Space
        CALL neci_flush(IUNIT)

    END SUBROUTINE FindSymSizeofTruncSpace

end module hilbert_space_size
