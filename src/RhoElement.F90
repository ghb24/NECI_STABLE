
FUNCTION Rho2OrderND2(NI, NJ, NEL, NBASISMAX, G1, NBASIS, IC2)
!.. We use a crude method and generate all possible 0th, 1st, and 2nd
!.. excitations of I and of J.  The intersection of these lists is the
!.. selection of dets we want.
    Use Determinants, only: get_helement
    use constants, only: dp
    use SystemData, only: BasisFN, tGUGA
    IMPLICIT NONE
    TYPE(BasisFN) G1(*)
    HElement_t(dp) Rho2OrderND2
    INTEGER NEL, NBASIS, nBasisMax(5, *)
    INTEGER NI(NEL), NJ(NEL), IC2
    INTEGER LSTI(NEL, NBASIS * NBASIS * NEL * NEL)
    INTEGER LSTJ(NEL, NBASIS * NBASIS * NEL * NEL)
    INTEGER NLISTI, NLISTJ, IC, I, J
    INTEGER ICI(NBASIS * NBASIS * NEL * NEL)
    INTEGER ICJ(NBASIS * NBASIS * NEL * NEL), NLISTMAX
    INTEGER CMP, IGETEXCITLEVEL, ICMPDETS
    HElement_t(dp) SUM1
    SUM1 = 0.0_dp
    NLISTMAX = NBASIS * NBASIS * NEL * NEL
    IC = IC2
    IF (IC < 0) IC = IGETEXCITLEVEL(NI, NJ, NEL)

!.. the 1 at the ends ensures K.NE. I or J, as this would make
!.. <I|U'|K> or <K|U'|J> zero (U' has a zero diag)
    CALL GENEXCIT(NI, 2, NBASIS, NEL, LSTI, ICI, NLISTI, 1, G1, .TRUE.,     &
&         NBASISMAX, .FALSE.)
    CALL GENEXCIT(NJ, 2, NBASIS, NEL, LSTJ, ICJ, NLISTJ, 1, G1, .TRUE.,     &
&         NBASISMAX, .FALSE.)
    I = 1
    J = 1
    if (tGUGA) then
        call stop_all("RHO2ORDERND2", "modify get_helement for GUGA")
    end if
!.. Now iterate over K, going along row I
    DO WHILE ((I <= NLISTI) .AND. (J <= NLISTJ))
        CMP = ICMPDETS(LSTI(1, I), LSTJ(1, J), NEL)
!.. While I>J, we increase J
        DO WHILE ((CMP > 0) .AND. (J < NLISTJ))
            J = J + 1
            CMP = ICMPDETS(LSTI(1, I), LSTJ(1, J), NEL)
        end do
        IF (CMP == 0) THEN
            SUM1 = SUM1 + get_helement(nI, lstI(:, I)) * &
                   get_helement(lstJ(:, J), nJ)
        end if
        I = I + 1

    end do
    RHO2ORDERND2 = SUM1
    RETURN
END

!  Get a matrix element of the double-counting corrected unperturbed Hamiltonian.
!  This is just the sum of the Hartree-Fock eigenvalues
!   with the double counting subtracted, Sum_i eps_i - 1/2 Sum_i,j <ij|ij>-<ij|ji>.  (i in HF det, j in excited det)
subroutine GetH0ElementDCCorr(nHFDet, nJ, nEl, G1, ECore, hEl)
    use constants, only: dp
    use Integrals_neci, only: get_umat_el
    use UMatCache
    use SystemData, only: BasisFN, Arr
    implicit none
    integer nEl
    integer nHFDet(nEl), nJ(nEl)
    type(BasisFN) G1(*)
    HElement_t(dp) hEl
    real(dp) ECore
    integer i, j
    integer IDHF(nEl), IDJ(nEl)
    hEl = (ECore)
    do i = 1, nEl
        hEl = hEl + (Arr(nJ(i), 2))
        IDHF(i) = gtID(nHFDet(i))
        IDJ(i) = gtID(nJ(i))
    end do
    do i = 1, nEl
        do j = 1, nEl
!Coulomb term
            hEl = hEl - (0.5_dp) * get_umat_el(IDHF(i), IDJ(j), IDHF(i), IDJ(j))
            if (G1(nHFDet(i))%Ms == G1(nJ(j))%Ms) then
!Exchange term
                hEl = hEl + (0.5_dp) * get_umat_el(IDHF(i), IDJ(j), IDJ(j), IDHF(i))
            end if
        end do
    end do
!         call writedet(77,nj,nel,.false.)
!         write(77,*) "H0DC",hEl
!         call neci_flush(77)
end
