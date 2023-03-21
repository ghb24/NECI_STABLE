SUBROUTINE DETHAM(NDET, NEL, NMRKS, HAMIL, LAB, NROW, TCOUNT, ICMAX, GC, TMC)
    use global_utilities, only: timer, halt_timer, set_timer
    use constants, only: dp, n_int
    Use Determinants, only: get_helement
    use SystemData, only: BasisFN, tGAS, t_lattice_model
    use SystemData, only: tGUGA
    use CalcData, only: TStar
    use bit_rep_data, only: NIfTot
    use SystemData, only: tHPHF
    use hphf_integrals, only: hphf_diag_helement
    use hphf_integrals, only: hphf_off_diag_helement
    use DetBitOps, only: EncodeBitDet
    use orb_idx_mod, only: SpinOrbIdx_t
    use gasci, only: GAS_spec => GAS_specification

    use guga_matrixElements, only: calcDiagMatEleGuga_nI, calc_guga_matrix_element
    use guga_data, only: ExcitationInformation_t
    use guga_bitrepops, only: CSF_Info_t

    use lattice_mod, only: get_helement_lattice

    use util_mod, only: near_zero

    IMPLICIT NONE
    HElement_t(dp) HAMIL(*)
    INTEGER NDET, NEL
    INTEGER LAB(*)
    INTEGER NMRKS(NEL, *)
    INTEGER NROW(NDET), GC
    INTEGER ICMAX, KI, IBEG, IBEGJ, KJ, IMAX, IDAMAX
    LOGICAL TCOUNT, TMC
    INTEGER STEP, IMAXJ
    integer(n_int) :: ilutI(0:NIfTot), ilutJ(0:NIfTot)
    HElement_t(dp) totSUM
    type(timer), save :: proc_timer
    type(ExcitationInformation_t) :: excitInfo
    integer :: ic = 0
    type(CSF_Info_t) :: csf_i, csf_j
!      LOGICAL TSTAR
! ==-------------------------------------------------------------------==
    proc_timer%timer_name = '    DETHAM'
    call set_timer(proc_timer)
! ==-------------------------------------------------------------------==
!..Global counter
    GC = 0
    NROW(1:NDET) = 0
    IBEG = 0
!..   Now we need to match up any two determinants
    DO KI = 1, NDET
        call EncodeBitDet(NMRKS(:, KI), ilutI)
        IF (mod(KI, 1000) == 0) WRITE(6, *) KI
        IF (KI == 1) THEN
            IBEG = 0
        ELSE
            IBEG = IBEG + NROW(KI - 1)
        END IF
        IF (TMC) THEN
            IBEGJ = 1
            STEP = 1
            IMAXJ = NDET
        ELSE
            IF (TSTAR) THEN
            IF (KI == 1) THEN
                IBEGJ = 1
                IMAXJ = NDET
                STEP = 1
            ELSE
                IBEGJ = KI
                IMAXJ = KI
                STEP = 1
            END IF
            ELSE
            IBEGJ = KI
            IMAXJ = NDET
            STEP = 1
            END IF
        END IF
        IF (STEP == 0) STEP = 1
        if (tGUGA) csf_i = CSF_Info_t(ilutI)

        DO KJ = IBEGJ, IMAXJ, STEP
            call EncodeBitDet(NMRKS(:, KJ), ilutJ)
            if (tGUGA) csf_j = CSF_Info_t(ilutJ)
            if (tHPHF) then
            if (KI == KJ) then
                totsum = hphf_diag_helement(NMRKS(:, KI), ilutI)
            else
                totsum = hphf_off_diag_helement(NMRKS(:, KI), NMRKS(:, KJ), ilutI, ilutJ)
            end if

            else if (tGUGA) then
            if (KI == KJ) then
                totsum = calcDiagMatEleGuga_nI(NMRKS(:, KI))
            else
                call calc_guga_matrix_element(ilutI, csf_i, ilutJ, csf_j, excitInfo, totsum, .true.)
            end if

            else if (t_lattice_model) then
            if (KI == KJ) then
                totsum = get_helement_lattice(NMRKS(:, KI), NMRKS(:, KJ), ic)
            else
                totsum = get_helement_lattice(NMRKS(:, KI), NMRKS(:, KJ))
            end if
            else
            if (KI == KJ) then
                totsum = get_helement(NMRKS(:, KI), NMRKS(:, KJ), 0)
            else
                totsum = get_helement(NMRKS(:, KI), NMRKS(:, KJ), ilutI, ilutJ)
            end if
            end if
            if (tGAS) then
            if (.not. GAS_spec%contains_conf(NMRKS(:, KI)) .or. .not. GAS_spec%contains_conf(NMRKS(:, KJ))) then
                totsum = 0.0_dp
            end if
            end if
            IF (ABS(TOTSUM) < 1.0e-10_dp) TOTSUM = 0.0_dp
            IF (.not. near_zero(TOTSUM) .OR. KI == KJ) THEN
                GC = GC + 1
!..   Stores the number of non-zero elements in each row
                NROW(KI) = NROW(KI) + 1
                IF (.NOT. TCOUNT) THEN
                    LAB(IBEG + NROW(KI)) = KJ
                    HAMIL(IBEG + NROW(KI)) = TOTSUM
                END IF
            END IF
        END DO
    END DO
!..No. of columns
    IF (TCOUNT) THEN
!        IMAX=MAXLOC(NROW)
        IMAX = IDAMAX(NDET, real(NROW, dp), 1)
        ICMAX = NROW(IMAX)
        WRITE(6, *) ' MAXIMUM WIDTH OF HAMIL : ', ICMAX
        WRITE(6, *) ' TOTAL NUMBER OF NON-ZERO ELEMENTS : ', GC
    END IF
! ==-------------------------------------------------------------------==
    call halt_timer(proc_timer)
! ==-------------------------------------------------------------------==
    RETURN
END
! ==---------------------------------------------------------------==
