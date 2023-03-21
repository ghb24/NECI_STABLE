! lenrec is the number of auxiliary basis functions
SUBROUTINE InitDFBasis(nBasisMax, Len)
    use SystemData, only: tStoreSpinOrbs
    use constants, only: dp, sp, stdout
    use UMatCache
    use util_mod, only: stop_all
    implicit none
    integer nBasisMax(5, *), Len
    character(*), parameter :: C_file = 'SAV_D____a'
    character(3) file_status
    integer(sp) info
    integer lenrec, nrec
    integer nBasis
    character(*), parameter :: this_routine = 'InitDFBasis'

    file_status = 'ADD'
    call stop_all(this_routine, "Reading in of SITUS DF files deprecated")
!         call init_record_handler(C_file,file_status,info)
!         call query_record_handler(C_file,info,file_status,lenrec,nrec,.TRUE.)
!.. lenrec is the number of auxiliary basis functions
!.. nrec is the number of pairs of orbitals.
    nrec = 0
    lenrec = 0
    nBasisPairs = nrec
    nBasis = int(sqrt(nBasisPairs * 2.0_dp))
    nAuxBasis = lenrec
    write(stdout, *) "DALTON/SITUS basis.", nBasis, " basis functions."
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
    tStoreSpinOrbs = .false.     !DF cannot cope currently with UHF/ROHF
END
!  A file to read in density fitted two-electron integrals from SITUS/Dalton.
!.. Also will read in one-electron integrals
!.. Requires various modules from SITUS for the file reading
SUBROUTINE ReadDF2EIntegrals(nBasis, nOrbUsed)
!         use precision
!         use record_handler
    use global_utilities
    use UMatCache
    use constants, only: dp, sp, stdout
    use util_mod, only: stop_all
    implicit none
    character(*), parameter :: C_file = 'SAV_D____a'
    character(*), parameter :: I_file = 'SAV_T____a'
    character(*), parameter :: S_file = 'SAV_S____a'
    character(*), parameter :: nolabel = '        '
    character(3) file_status
    character(*), parameter :: t_r = 'ReadDF2EIntegrals'
    integer(sp) info
    integer i, j, k
    real(dp) r
    integer nBasis, nOrbUsed, ierr
    character(*), parameter :: this_routine = 'ReadDF2EIntegrals'

    call stop_all(this_routine, "Reading in of SITUS DF files deprecated")
    write(stdout, *) "Opening Density fitting matrix files"
    file_status = 'ADD'
!.. We've already got C_file open
!         call init_record_handler(I_file,file_status,info,printinfo=.TRUE.)
!         call init_record_handler(S_file,file_status,info,printinfo=.TRUE.)
!.. lenrec is the number of auxiliary basis functions
!.. nrec is the number of pairs of orbitals.
    write(stdout, *) "Basis Size:", nBasis / 2
    write(stdout, *) "Auxiliary Basis Size:", nAuxBasis

    tDFInts = .TRUE.
    allocate(DFCoeffs(nAuxBasis, nBasisPairs), STAT=ierr)
    call LogMemAlloc("DFCoeffs", nBasisPairs * nAuxBasis, 8, t_r, tagDFCoeffs, ierr)
    allocate(DFInts(nAuxBasis, nBasisPairs), STAT=ierr)
    call LogMemAlloc("DFInts", nBasisPairs * nAuxBasis, 8, t_r, tagDFInts, ierr)
    allocate(DFFitInts(nAuxBasis, nAuxBasis), STAT=ierr)
    call LogMemAlloc("DFFitInts", nAuxBasis * nAuxBasis, 8, t_r, tagDFFitInts, ierr)
    do i = 1, nBasisPairs
!            call read_record(C_file,i,nolabel,DFCoeffs(:,i),info)
!            call read_record(I_file,i,nolabel,DFInts(:,i),info)
    end do
    do i = 1, nAuxBasis
!            call read_record(S_file,i,nolabel,DFFitInts(:,i),info)
    end do
!         call leave_record_handler(C_file,info)
!         call leave_record_handler(I_file,info)
!         call leave_record_handler(S_file,info)

    select case (iDFMethod)
    case (0)
        call stop_all(this_routine, "No DF, but ReadDF2EIntegrals called.")
    case (1)
        write(stdout, *) "DFMETHOD DFOVERLAP        1 - (ij|u|ab)= (ij|u|P)(P|ab)"
    case (2)
        write(stdout, *) "DFMETHOD DFOVERLAP2NDORD  2 - (ij|u|ab)= (ij|u|P)(P|ab)+(ij|P)(P|u|ab)-(ij|P)(P|u|Q)(Q|ab)"
    case (3)
        write(stdout, *) "DFMETHOD DFOVERLAP2       3 - (ij|u|ab)= (ij|P)(P|u|Q)(Q|ab)"
    case (4)
        write(stdout, *) "DFMETHOD DFCOULOMB        4 - (ij|u|ab)= (ij|u|P)[(P|u|Q)^-1](Q|u|ab)"
    case default
        write(stdout, *) "Unknown DF Method: ", iDFMethod
        call stop_all(this_routine, "Unknown DF Method")
    end select

    if (iDFMethod > 2) then
        allocate(DFInvFitInts(nAuxBasis, nAuxBasis), STAT=ierr)
        call LogMemAlloc("DFInvFitInts", nAuxBasis * nAuxBasis, 8, t_r, tagDFInvFitInts, ierr)
    end if
    select case (iDFMethod)
    case (3)
        call DFCalcInvFitInts(0.5_dp)
        write(stdout, *) "Calculating B matrix"
        do i = 1, nBasisPairs
            do j = 1, nAuxBasis
                r = 0
                do k = 1, nAuxBasis
                    r = r + DFInvFitInts(j, k) * DFCoeffs(k, i)
                end do
                DFInts(j, i) = r
            end do
        end do
!DFInts now contains B_ij,P = sum_Q (ij|Q)[(Q|u|P)^1/2]
    case (4)
        call DFCalcInvFitInts(-0.5_dp)
        write(stdout, *) "Calculating B matrix"
        do i = 1, nBasisPairs
            do j = 1, nAuxBasis
                r = 0
                do k = 1, nAuxBasis
                    r = r + DFInvFitInts(j, k) * DFInts(k, i)
                end do
                DFCoeffs(j, i) = r
            end do
        end do
!DFCoeffs now contains B_ij,P = sum_Q (ij|u|Q)[(Q|u|P)^-1/2]
    end select
    IF (nBasis /= nOrbUsed) THEN
! We allocate a small preliminary cache before freezing.
        call SetupUMatCache(nOrbUsed / 2, .TRUE.)
    ELSE
        call SetupUMatCache(nOrbUsed / 2, .FALSE.)
    end if
    call SetupUMat2D_df
END

!.. Get a 2-el integral.  a,b,c,d are indices. <ab|1/r12|cd>
!DFCoeffs(x,yz) is (x|yz)
!DFInts(x,yz) is (x|u|yz)
!DFFitInts(x,y) is (x|u|y)
SUBROUTINE GetDF2EInt(a, b, c, d, res)
    use constants, only: dp
    use UMatCache
    use util_mod, only: stop_all
    implicit none
    integer a, b, c, d
    integer i, j, GetDFIndex
    integer x, y
    real(dp) res
    character(*), parameter :: this_routine = 'GetDF2EInt'

    res = 0.0_dp
    x = GetDFIndex(a, c)
    y = GetDFIndex(b, d)
!         write(stdout,"(7I4)") a,b,c,d,x,y,iDFMethod
    select case (iDFMethod)
! 0 - no DF
    case (0)
        call stop_all(this_routine, "DF Method 0 specified - no DF, but DF called.")
    case (1)
! DFOVERLAP        1 - (ij|u|ab)= (ij|u|P)(P|ab)
        do i = 1, nAuxBasis
            res = res + DFCoeffs(i, x) * DFInts(i, y)
        end do
    case (2)
! DFOVERLAP2NDORD  2 - (ij|u|ab)= (ij|u|P)(P|ab)+(ij|P)(P|u|ab)-(ij|P)(P|u|Q)(Q|ab)
        do i = 1, nAuxBasis
            res = res + DFCoeffs(i, x) * DFInts(i, y) + DFCoeffs(i, y) * DFInts(i, x)
            do j = 1, nAuxBasis
                res = res - DFCoeffs(i, x) * DFCoeffs(j, y) * DFFitInts(i, j)
            end do
        end do
    case (3)
! DFOVERLAP2       3 - (ij|u|ab)= (ij|P)(P|u|Q)(Q|ab)

!DFInts actually contains B_ij,P = sum_Q (ij|Q)[(Q|u|P)^1/2]
        do i = 1, nAuxBasis
            res = res + DFInts(i, x) * DFInts(i, y)
        end do
    case (4)
! DFCOULOMB        4 - (ij|u|ab)= (ij|u|P)[(P|u|Q)^-1](Q|u|ab)
!DFCoeffs actually contains B_ij,P = sum_Q (ij|u|Q)[(Q|u|P)^-1/2]
        do i = 1, nAuxBasis
            res = res + DFCoeffs(i, x) * DFCoeffs(i, y)
        end do
    end select
END

!.. return a DF pair index - i<j (although the pairs are ordered 11 21 22 31 32 33 41 42 ...
INTEGER FUNCTION GetDFIndex(i, j)
    IMPLICIT NONE
    INTEGER I, J
    if (i < j) then
        GetDFIndex = i + j * (j - 1) / 2
    else
        GetDFIndex = j + i * (i - 1) / 2
    end if
END
SUBROUTINE ReadDalton2EIntegrals(nBasis, UMat2D, tUMat2D)
    use constants, only: dp, sp, stdout
    implicit none
    integer nBasis, i, j, k, ilast
    real(dp) val, UMat2D(nBasis, nBasis)
    logical tUMat2D
    tUMat2D = .false.
    open(11, file='HONEEL', status='unknown')
    i = 1
    do while (i /= 0)
        read(11, *) i, j, val
    end do
    i = 1
    ilast = 0
    do while (i /= 0 .and. i <= nBasis .and. j <= nBasis)
        read(11, *, end=20) i, j, k, val
        if (i /= 0 .and. i <= nBasis .and. j <= nBasis) then
            UMat2D(i, j) = val
            ilast = i
        end if
    end do
    tUMat2D = .true.
20  close(11)
    IF (tUMat2D) THEN
        write(stdout, *) "Read in 2-index 2-electron integrals up to orbital ", ilast * 2
    ELSE
        write(stdout, *) "Have not read in 2-index 2-electron integrals."
    end if
    IF (ilast < nBasis) THEN
        write(stdout, *) "Calculating remaining 2-index 2-electron integrals with density fitting."
        open(78, file='UMAT', status='UNKNOWN')
        do i = ilast + 1, nBasis
            do j = 1, i
                if (i < j) then
                    call GetDF2EInt(i, j, i, j, UMat2D(i, j))
                    call GetDF2EInt(i, j, j, i, UMat2D(j, i))
                    write(78, "(3I5,G25.16)") i, j, 0, UMat2D(i, j)
                    IF (i /= j) write(78, "(3I5,G25.16)") j, i, 0, UMat2D(j, i)
                else
                    call GetDF2EInt(i, j, i, j, UMat2D(j, i))
                    call GetDF2EInt(i, j, j, i, UMat2D(i, j))
                    write(78, "(3I5,G25.16)") i, j, 0, UMat2D(j, i)
                    IF (i /= j) write(78, "(3I5,G25.16)") j, i, 0, UMat2D(i, j)
                end if
            end do
        end do
        close(78)
    end if
END
SUBROUTINE ReadDalton1EIntegrals(G1, nBasis, ECore)
    use constants, only: dp
    use SystemData, only: BasisFN, BasisFNSize, Symmetry, NullBasisFn
    USE OneEInts, only: TMATind, TMAT2D, TMATSYM
    use sym_mod, only: TotSymRep
    implicit none
    integer nBasis, i, j
    real(dp) val, ECore
    type(BasisFN) G1(nBasis)
    open(11, file='HONEEL', status='unknown')
    i = 1
    !TMat=0.0_dp
    G1(1:nBasis) = NullBasisFn
    do while (i /= 0)
        read(11, *) i, j, val
!"(2I5,F)"
        if (i == 0) then
            ECore = val
        else if (j /= 0) then
            TMat2D(i * 2 - 1, j * 2 - 1) = (val)
            TMat2D(i * 2, j * 2) = (val)
            TMat2D(j * 2 - 1, i * 2 - 1) = (val)
            TMat2D(j * 2, i * 2) = (val)
        end if
    end do
    close(11)
    do i = 1, nBasis
        G1(i)%Ms = 1 - 2 * iand(i, 1)
        G1(i)%Sym = TotSymRep()
!  We've already read in and ordered the Energies
!            Arr(i,1)=TMat(i,i)
!            Arr(i,2)=TMat(i,i)
!            Brr(i)=i
    end do
END
SUBROUTINE InitDaltonBasis(Arr, Brr, G1, nBasis)
    use constants, only: dp
    use SystemData, only: Symmetry, BasisFN, BasisFNSize, NullBasisFn
    use SymData, only: tAbelian
    use sym_mod, only: TotSymRep
    implicit none
    integer nBasis, Brr(nBasis), i, j
    real(dp) Arr(nBasis, 2), val
    type(BasisFN) G1(nBasis)
    tAbelian = .true.
    open(11, file='HONEEL', status='unknown')
    i = 1
    G1(1:nBasis) = NullBasisFn
    do while (i /= 0)
        read(11, *) i, j, val
        if (j == 0 .and. i /= 0) then
            Arr(i * 2 - 1, 1) = val
            Arr(i * 2 - 1, 2) = val
            Arr(i * 2, 1) = val
            Arr(i * 2, 2) = val
        end if
    end do
    close(11)
    do i = 1, nBasis
        G1(i)%Ms = 1 - 2 * iand(i, 1)
        G1(i)%Sym = TotSymRep()
        Brr(i) = i
    end do
END

!.. Get a 2-el integral.  a,b,c,d are indices. <ab|1/r12|cd>
!DFCoeffs(x,yz) is (x|yz)
!DFInts(x,yz) is (x|u|yz)
!DFFitInts(x,y) is (x|u|y)
!This is slower but calculates more accurately.
SUBROUTINE DFCalcInvFitInts(dPower)
    use constants, only: dp, sp, stdout
    use UMatCache
    use global_utilities
    use MemoryManager, only: TagIntType
    use util_mod, only: stop_all
    implicit none
    real(dp), Pointer :: M(:, :) !(nAuxBasis,nAuxBasis)
    real(dp) Eigenvalues(nAuxBasis), r, dPower
    real(dp) Work(3 * nAuxBasis)
    integer Workl
    integer(sp) info
    integer(TagIntType), save :: tagM = 0
    type(timer), save :: proc_timer
    character(*), parameter :: t_r = 'DFCalcInvFitInts'
    Integer i, j, ierr, k, iMinEigv
    proc_timer%timer_name = 'DFInvFitIn'
    call set_timer(proc_timer)
    allocate(M(nAuxBasis, nAuxBasis), STAT=ierr)
    call LogMemAlloc("M-DFInvFitInts", nAuxBasis * nAuxBasis, 8, t_r, tagM, ierr)
    M = 0.0_dp
    do i = 1, nAuxBasis
        do j = 1, i
            M(i, j) = DFFitInts(i, j)
        end do
    end do
    Workl = 3 * nAuxBasis
    write(stdout, *) "Diagonalizing (P|u|Q)"
    CALL DSYEV('V', 'L', nAuxBasis, M, nAuxBasis, Eigenvalues, WORK, WORKL, INFO)
    IF (INFO /= 0) THEN
        write(stdout, *) 'DYSEV error: ', INFO
        call stop_all(t_r, "DSYEV error")
    end if
    iMinEigv = 1
    if (dPower < 0) then
        do i = 1, nAuxBasis
            if (Eigenvalues(i) < 1d-10) iMinEigv = i + 1
        end do
        write(stdout, *) "Ignoring ", iMinEigv - 1, " eigenvalues <1d-10"
    end if
!M now contains eigenvectors.  Eigenvector i is in M(:,I)
! A=U^T L U (L is the matrix of eigenvalues on the diagonal)
! A^-1 = U^T L^-1 U
    write(stdout, *) "Calculating (P|u|Q)^", dPower
    do i = 1, nAuxBasis
        do j = 1, i
            r = 0
            do k = iMinEigv, nAuxBasis
                r = r + (Eigenvalues(k)**dPower) * M(j, k) * M(i, k)
            end do
            DFInvFitInts(i, j) = r
            DFInvFitInts(j, i) = r
        end do
    end do
    call LogMemDealloc(t_r, tagM)
    Deallocate(M)
    call halt_timer(proc_timer)
END
