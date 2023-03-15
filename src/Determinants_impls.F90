#include "macros.h"

submodule (Determinants) Determinants_impls
    use SymExcitDataMod, only: ScratchSize, ScratchSize3
    use GenRandSymExcitNUMod, only: SpinOrbSymSetup
    use sym_mod, only: writesymtable, gensymstatepairs
    use SystemData, only: nEl

    better_implicit_none

contains

    subroutine virt_uniform_sym_setup()
        ! We use the third scratch array to store data for single
        ! excitations

        call SpinOrbSymSetup()

        ScratchSize3 = ScratchSize

    end subroutine


    module subroutine DetInit()
        real(dp) :: DNDET
        integer :: i, j
        integer(int64) :: nDet
        integer :: alpha, beta, symalpha, symbeta, endsymstate
        LOGICAL :: tSuccess, tFoundOrbs(nBasis)

        write (stdout, *) "SYMMETRY MULTIPLICATION TABLE"
        CALL WRITESYMTABLE(6)

        CALL GENSymStatePairs(NBASIS / 2, .false.)

!iActiveBasis is a copy of nPaths
        IF (iActiveBasis == -2) then
!  PATHS ACTIVE SETS
            Call GenActiveBasis(ARR, nBasis, nEl, nActiveBasis, nActiveSpace(1), nActiveSpace(2))
        else if (iActiveBasis == -3) then
!  PATHS ACTIVE ORBITALS
            nActiveBasis(1) = nEl + 1 - nActiveSpace(1)
            nActiveBasis(2) = nEl + nActiveSpace(2)
            write (stdout, *) "Active space:", nActiveBasis(1), " TO ", nActiveBasis(2), " (ordered labels)."
            write (stdout, *) "Active space electrons:", nEl - nActiveBasis(1) + 1
        else
            nActiveBasis(1) = 1
            nActiveBasis(2) = nBasis
        end if
!C.. Work out a preliminary Fermi det
!      IF(FDET(1).EQ.0) THEN

!C.. Check if we're blocking the hamiltonian
!C      IF(THFBASIS.AND.TBLOCK) THEN
!C         write(stdout,*) "THFBASIS set and NBLK=0.  ",
!C     &         "Cannot block diagonalize in HF Basis."
!C         STOP
!C      end if
!C      CALL SYMGENEXCITS(FDET,NEL,NBASIS)
!C      CALL LeaveMemoryManager
!C      STOP

!C.. in order to calculate the H matrix, we need to work out all the determinants
!C.. beware with NPATHS - it calcs the list of dets even if we don't calc H
!C.. Could be big.
!C..Now we see how many determinants we need
!C      IF(nBasis.GT.170) THEN
!C..This fix is to stop floating overflow as taking the factorial of (nBasis.GT.170) crashes
!C  using the old FACTRL routine.
        NDET = 1
        DNDET = 1.0_dp
        DO I = 0, NEL - 1
            NDET = (NDET * (nBasis - I)) / (I + 1)
            DNDET = (DNDET * real(nBasis - I, dp)) / real(I + 1, dp)
        end do

        IF (abs(real(NDET) - dndet) > 1.0e-6) THEN
!         write(stdout,*) ' NUMBER OF DETERMINANTS : ' , DNDET
            NDET = -1
        ELSE
!         write(stdout,*) ' NUMBER OF DETERMINANTS : ' , NDET
        end if

!C      CALL TC(I_HMAX,I_P,NWHTAY)

!Check that the symmetry routines have set the symmetry up correctly...
        tSuccess = .true.
        tFoundOrbs(:) = .false.

        IF ((.not. tHub) .and. (.not. tUEG) .and. TwoCycleSymGens) THEN
            do i = 1, nSymLabels
                EndSymState = SymLabelCounts(1, i) + SymLabelCounts(2, i) - 1
                do j = SymLabelCounts(1, i), EndSymState

                    Beta = (2 * SymLabelList(j)) - 1
                    Alpha = (2 * SymLabelList(j))
                    SymAlpha = INT((G1(Alpha)%Sym%S), 4)
                    SymBeta = INT((G1(Beta)%Sym%S), 4)

                    IF (.not. tFoundOrbs(Beta)) THEN
                        tFoundOrbs(Beta) = .true.
                    ELSE
                        CALL Stop_All("SetupParameters", "Orbital specified twice")
                    end if
                    IF (.not. tFoundOrbs(Alpha)) THEN
                        tFoundOrbs(Alpha) = .true.
                    ELSE
                        CALL Stop_All("SetupParameters", "Orbital specified twice")
                    end if

                    IF (G1(Beta)%Ms /= -1) THEN
                        tSuccess = .false.
                    else if (G1(Alpha)%Ms /= 1) THEN
                        tSuccess = .false.
                    else if ((SymAlpha /= (i - 1)) .or. (SymBeta /= (i - 1))) THEN
                        tSuccess = .false.
                    end if
                end do
            end do
            do i = 1, nBasis
                IF (.not. tFoundOrbs(i)) THEN
                    write (stdout, *) "Orbital: ", i, " not found."
                    CALL Stop_All("SetupParameters", "Orbital not found")
                end if
            end do
        end if
        ! SpinOrbSymSetup currently sets up the symmetry arrays for use with
        ! symrandexcit2 excitation routines. These are not currently
        ! compatible with non-abelian symmetry groups, which CPMD jobs
        ! invariably used. To avoid this complication, this symmetry
        ! setup will not be used with CPMD, and thus these excitation
        ! generators won't work.
        IF (.not. tCPMD) THEN
            IF (.not. tSuccess) THEN
                write (stdout, *) "************************************************"
                write (stdout, *) "**                 WARNING!!!                 **"
                write (stdout, *) "************************************************"
                write (stdout, *) "Symmetry information not set up correctly in NECI initialisation"
                write (stdout, *) "Will attempt to set up the symmetry again, but now in terms of spin orbitals"
                write (stdout, *) "Old excitation generators will not work"
                write (stdout, *) "I strongly suggest you check that the reference energy is correct."
                !CALL SpinOrbSymSetup() !.true.)
            ELSE
                write (stdout, *) "Symmetry and spin of orbitals correctly set up for excitation generators."
                write (stdout, *) "Simply transferring this into a spin orbital representation."
                !CALL SpinOrbSymSetup() !.false.)
            end if

            if (tPickVirtUniform) then
                call virt_uniform_sym_setup()
            else
                ! Includes normal & HPHF
                call SpinOrbSymSetup()
            end if

        end if
! From now on, the orbitals are also contained in symlabellist2 and symlabelcounts2.
! These are stored using spin orbitals.
    end subroutine DetInit


    ! Generate the active space from a basis.
    ! The Active basis can be used to in PATHS calculations and (?) as a CASCI

    ! nActiveBasis(1:2) contains (First Active Basis Fn, Last Active Basis Fn)
    ! nDown is the number of orbital sets  below the Fermi level
    ! nUp is the number of orbital sets  above the Fermi level
    subroutine GenActiveBasis(ARR, nBasis, nEl, nActiveBasis, nDown, nUp)
        integer nEl, nActiveBasis(2), nBasis
        real(dp) ARR(nBasis)
        integer I, nDown, nUp, nLeft
        I = nEl + 1
        nLeft = 1 + nUp
        IF (nDown /= 0 .AND. nUp /= 0) write (stdout, *) "Including ", -nDown, ",", nUp, " extra degenerate sets in active space."
        DO WHILE (nLeft > 0 .AND. I < nBasis)
            DO WHILE (I < nBasis .AND. ABS(ARR(I) - ARR(I - 1)) < 1.0e-5_dp)
                I = I + 1
            end do
            nLeft = nLeft - 1
            IF (nLeft == nUp .AND. I /= nEl + 1) write (stdout, *) "Fermi determinant degenerate.  "
            IF (nLeft /= 0) I = I + 2
        end do
        IF (I == nEl + 1 .and. nDown == 0) THEN
    !We have no degeneracies at the Fermi Energy
            write (stdout, *) "Fermi determinant non-degenerate.  "
            IF (nDown == 0) THEN
                write (stdout, *) "Active space empty."
                nActiveBasis(1) = nEl + 1
                nActiveBasis(2) = nEl
                RETURN
            end if
        end if
        nActiveBasis(2) = I - 1
        I = nEl - 1
        nLeft = nDown
        Do WHILE (nLeft > 0 .AND. I > 0)
            DO WHILE (I > 0 .AND. ABS(ARR(I) - ARR(I + 1)) < 1.0e-5_dp)
                I = I - 1
            end do
            nLeft = nLeft - 1
        end do
        nActiveBasis(1) = I + 1
        write (stdout, *) "Active space:", nActiveBasis(1), " TO ", nActiveBasis(2), " (ordered labels)."
        write (stdout, *) "Active space electrons:", nEl - nActiveBasis(1) + 1
    end subroutine


end submodule Determinants_impls
