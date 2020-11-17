subroutine InitRIBasis(nBasisMax, Len)
    use constants, only: dp, int64
    Use SymData, only: tAbelian
! lenrec is the number of auxiliary basis functions
    use UMatCache
    use util_mod, only: record_length
    implicit none
    integer nBasisMax(5, *), Len
    integer nBasis
    integer(int64) nAb, nB
    tAbelian = .true.
    open(29, file='RIINTDUMP', status='old', FORM='UNFORMATTED', access='DIRECT', recl=record_length(8))
!.. The first element is the number of aux basis fns.
!.. The second element is the number of basisfunctions.
    read(29, rec=1) nAb
    read(29, rec=2) nB
    write(6, *) nAb, nB
    nAuxBasis = int(nAb, sizeof_int)
    nBasis = int(nB, sizeof_int)
    write(6, *) "Q-Chem auxiliary basis", nAuxBasis, " basis functions:", nBasis
    nBasisMax(1:5, 1:3) = 0
    Len = 2 * nBasis
!.. Note that it's a read in basis.
    nBasisMax(3, 3) = 1
    nBasisMax(4, 1) = -1
    nBasisMax(4, 2) = 1
!.. Correspond to ISS=0
    nBasisMax(1, 3) = 0
!.. Setup Max Sym
    nBasisMax(5, 2) = 0
    close(29)
END

SUBROUTINE GetRI2EInt(a, b, c, d, res)
    use constants, only: dp
    use UMatCache
    implicit none
    integer a, b, c, d
    integer i, GetDFIndex
    integer x, y
    real(dp) res
    res = 0.0_dp
    x = GetDFIndex(a, c)
    y = GetDFIndex(b, d)
! DFOVERLAP        1 - (ij|u|ab)= (ij|u|P)(P|ab)
    do i = 1, nAuxBasis
        res = res + DFCoeffs(i, x) * DFInts(i, y)
    end do
end

SUBROUTINE ReadRI2EIntegrals(nBasis, nOrbUsed)
    use UMatCache
    IMPLICIT NONE
    INTEGER nBasis, nOrbUsed
    IF (nBasis /= nOrbUsed) THEN
! We allocate a small preliminary cache before freezing.
        write(6, *) "Setting up pre-freezing UMatCache"
        call SetupUMatCache(nOrbUsed / 2, .TRUE.)
    ELSE
        call SetupUMatCache(nOrbUsed / 2, .FALSE.)
    end if
    call SetupUMat2D_df()
!Integrals are actually read in READFCIINT.  later on
END

SUBROUTINE ReadRIIntegrals(nBasis, nOrbUsed)
    use UMatCache
    use global_utilities
    use util_mod, only: record_length
    use constants, only: dp, int64
    IMPLICIT NONE
    character(*), parameter :: t_r = 'ReadRIIntegrals'
    INTEGER nBasis, nOrbUsed
    integer i, j, onints, nints, ierr, Q
    real(dp) val
    integer(int64) nA, nB
    integer GetDFIndex
    write(6, *) "Reading QChem C Matrices"
    open(29, file='RIINTDUMP', status='old', FORM='UNFORMATTED', access='DIRECT', recl=record_length(8))
    read(29, rec=1) nA
    read(29, rec=2) nB
    onints = 2
    nints = (nBasis / 2 * (nBasis / 2 + 1)) / 2  !/2 to convert to states
    allocate(DFInts(nAuxBasis, nints), STAT=ierr)
    call LogMemAlloc("DFInts", nints * nAuxBasis, 8, t_r, tagDFInts, ierr)
    write(6, *) nints, nAuxBasis
    DO I = 1, nBasis / 2
        do Q = 1, nAuxBasis
            do J = 1, nBasis / 2
                onints = onints + 1
                read(29, rec=onints) val
                DFInts(Q, GetDFIndex(j, i)) = val
            end do
        end do
    end do
    write(6, *) "Over"
    close(29)
    iDFMethod = 3
    CALL ReadRI2EIntegrals(nBasis, nOrbUsed)
    write(6, *) "DFM:", iDFMethod
end
