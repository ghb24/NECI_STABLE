#include "macros.h"
module Determinants
    use DetBitOps, only: EncodeBitDet, count_open_orbs, spatial_bit_det, GetBitExcitation
    use DeterminantData, only: Fdet, calculated_ms, tagfdet, write_det, write_det_len
    use IntegralsData, only: UMat, FCK, NMAX, nFrozen, nFrozenIn
    use MemoryManager, only: TagIntType
    use OneEInts, only: GetTMatEl
    use SymData, only: nSymLabels, SymLabelList, SymLabelCounts, TwoCycleSymGens
    use SystemData, only: alat, arr, arr, basisfn, basisfnsize, basisfnsizeb, &
        brr, ecore, g1, irreporboffset, k_lattice_constant, k_lattice_vectors, &
        k_momentum, lms, lnosymmetry, nbasis, nbasismax, &
        nclosedorbs, nel, nirreps, nmsh, noccorbs, stot, symmetry, symmetrysize, &
        symmetrysizeb, symrestrict, t_guga_noreorder, t_lattice_model, tcpmd, &
        tfixlz, tguga, thphfints, thub, tmolpro, tpickvirtuniform, &
        tref_not_hf, tspn, tstoreasexcitations, tsymset, tueg, tueg2, &
        tuegspecifymomentum, tusebrillouin, t_k_space_hubbard, &
        tStoquastize
    use bit_rep_data, only: nIfTot
    use blas_interface_mod, only: dcopy
    use constants, only: Pi, Pi2, THIRD, dp, n_int, bits_n_int, int64, maxExcit, stdout
    use excitation_types, only: Excitation_t, Excite_2_t, get_excitation
    use global_utilities, only: timer, set_timer, halt_timer, LogMemAlloc
    use guga_data, only: ExcitationInformation_t
    use guga_matrixElements, only: calcDiagMatEleGUGA_nI
    use lattice_mod, only: get_helement_lattice, lat
    use sltcnd_mod, only: sltcnd, dyn_sltcnd_excit_old, sltcnd_compat, &
        sltcnd_excit, sltcnd_knowIC, SumFock, CalcFockOrbEnergy
    use sort_mod, only: sort
    use sym_mod, only: writesym, GENMOLPSYMTABLE, getsym, MomPbcSym, GetLz, checkMomentumInvalidity
    use util_mod, only: NECI_ICOPY, operator(.div.)
    use error_handling_neci, only: stop_all
    use fmt_utils, only: int_fmt

    better_implicit_none
    private
    public :: get_helement, modifymomentum, get_helement_excit, &
        orderbasis, writebasis, isuhfdet, calcT, NUHFDet, &
        GetUEGKE, DetFreezeBasis, GetH0Element3, tagSpecDet, specdet, tSpecDet, &
        geth0element4, tDefineDet, tagDefDet, DefDet, WriteDetBit, get_helement_det_only, &
        iactivebasis, nActiveBasis, nActiveSpace, DetInit, DetPreFreezeInit, &
        DetPreFreezeInit_old, DetCleanup

    ! TODO: Add an interface for getting a diagonal helement with an ordered
    !       list, or with only a bit-det
    interface get_helement
        module procedure get_helement_compat
        module procedure get_helement_excit
        module procedure get_helement_normal
    end interface

    interface
        module subroutine DetInit()
        end subroutine
    end interface

    save
! Set by Calc on input
    integer :: nActiveSpace(2)
    integer, DIMENSION(:), POINTER :: SPECDET => null()
    integer(TagIntType) :: tagSPECDET = 0
    logical TSPECDET

!nActiveBasis(1) is the lowest non-active orbital
!nActiveBasis(2) is the highest active orbital.  There can be virtuals above this.
!  Active orbitals are used for generating the determinants whose energy/weight is to be found
    integer nActiveBasis(2)
!  Set by input to indicate which type of active basis we need
    integer iActiveBasis
!Used to be from uhfdet.inc
    integer nUHFDet(5000)
    real(dp) E0HFDet

    integer, allocatable :: DefDet(:)
    logical :: tDefineDet
    integer(TagIntType) :: tagDefDet = 0

contains

    Subroutine DetPreFreezeInit_old()
        integer :: ierr, ms, iEl, flagAlpha, msTmp
        integer :: i, j, Lz, OrbOrder(8, 2), FDetTemp(NEl)
        type(BasisFn) s
        logical :: tGenFDet
        HElement_t(dp) :: OrbE
        character(25), parameter :: this_routine = 'DetPreFreezeInit_old'
        allocate (FDet(nEl), stat=ierr)
        LogAlloc(ierr, 'FDet', nEl, 4, tagFDet)
        IF (tDefineDet) THEN
            write (stdout, *) 'Defining FDet according to input'
            do i = 1, NEl
                FDet(i) = DefDet(i)
            end do
            call assignOccOrbs()

            ! A quick check that we have defined a reasonable det.
            ms = sum(get_spin_pn(fdet(1:nel)))
            if (abs(ms) /= abs(lms)) then
                write (stdout, *) 'LMS', lms
                write (stdout, *) 'Calculated Ms', ms
                call stop_all(this_routine, "Defined determinant has the &
                     &wrong Ms value. Change DEFINEDET or &
                     &SPIN-RESTRICT")
            end if

            if (tGUGA) then
                ms = abs(ms)
                if (ms < 0 .or. ms /= STOT) then
                    write (stdout, *) "S: ", STOT
                    write (stdout, *) "calculated S of inputted CSF: ", ms
                    call stop_all(this_routine, " Defined CSF has the &
                        &wrong total spin quantum number. Change DEFINEDET or &
                        &S quantum numnber!")
                end if
            end if

            tRef_Not_HF = .true.
        else
            if ((sum(nOccOrbs) + sum(nClosedOrbs)) == nel) then
                tGenFDet = .false.
                iEl = 1
                msTmp = -1 * lms
                do i = 1, nIrreps
                    ! doubly occupy the closed orbs
                    do j = 1, nClosedOrbs(i)
                        FDet(iEl) = irrepOrbOffset(i) + 2 * j - 1
                        iEl = iEl + 1
                        FDet(iEl) = irrepOrbOffset(i) + 2 * j
                        iEl = iEl + 1
                    end do
                    ! now distribute electrons to the open orbs
                    ! such that the total sz matches the requested
                    ! only consider occ orbs which are not closed
                    do j = nClosedOrbs(i) + 1, nOccOrbs(i)
                        if (msTmp < 0) then
                            flagAlpha = 0
                        else
                            flagAlpha = 1
                        end if
                        FDet(iEl) = irrepOrbOffset(i) + 2 * j - flagAlpha
                        iEl = iEl + 1
                        msTmp = msTmp + (1 - 2 * flagAlpha)
                    end do
                end do
                call sort(FDet)
            else
                CALL GENFDET(FDET)
                call assignOccOrbs()
                IF (tUEGSpecifyMomentum) THEN
                    write (stdout, *) 'Defining FDet according to a momentum input'
                    CALL ModifyMomentum(FDET)
                end if
                tRef_Not_HF = .false.
            end if
        end if
        write (stdout, "(A)", advance='no') " Fermi det (D0):"
        call write_det(stdout, FDET, .true.)
        Call GetSym(FDet, nEl, G1, nBasisMax, s)
        write (stdout, "(A)", advance='no') " Symmetry: "
        Call WriteSym(6, s%Sym, .true.)
        IF (tFixLz) THEN
            Call GetLz(FDet, nEl, Lz)
            write (stdout, "(A,I5)") "Lz of Fermi det:", Lz
        end if
        CALL NECI_ICOPY(NEL, FDET, 1, NUHFDET, 1)
        if (tMolpro) then
            !Orbitals are ordered by occupation number from MOLPRO, and not reordered in NECI
            !Therefore, we know HF determinant is first four occupied orbitals.
            write (stdout, "(A)") "NECI called from MOLPRO, so assuming orbitals ordered by occupation number."
            if (.not. tDefineDet) then
                FDetTemp(:) = FDet(:)
            else
                !We have defined our own reference determinant, but still use the first orbitals for the calculation
                !of 'orbital energies'
                CALL GENFDET(FDETTEMP)
            end if
            write (stdout, "(A)") "Calculating orbital energies..."
            do i = 1, nBasis
                OrbE = CalcFockOrbEnergy(i, FDetTemp)
                Arr(i, 1) = real(OrbE, dp)
                Brr(i) = i
            end do
            write (stdout, "(A)") "Reordering basis by orbital energies..."
            OrbOrder(:, :) = 0
            call ORDERBASIS(NBASIS, Arr, Brr, OrbOrder, nBasisMax, G1)
            !However, we reorder them here
            call writebasis(6, G1, nBasis, Arr, Brr)
        end if
        E0HFDET = ECORE
        DO I = 1, NEL
            E0HFDET = E0HFDET + ARR(NUHFDET(i), 2)
        end do
        write (stdout, *) "Fock operator energy:", E0HFDET

        ! Store the value of Ms for use in other areas
        calculated_ms = sum(get_spin_pn(fdet(1:nel)))

        if (tGUGA) then
            calculated_ms = abs(calculated_ms)
        end if

    contains

        subroutine assignOccOrbs
            integer :: k
            ! assign the occ/closed orbs
            nOccOrbs = 0
            nClosedOrbs = 0
            do k = 1, nel
                if (k > 1) then
                    if (is_alpha(FDet(k)) .and. FDet(k - 1) == FDet(k) - 1) then
                        ! we do not need to resolve the irrep anymore (not relevant
                        ! for further usage)
                        nClosedOrbs(1) = nClosedOrbs(1) + 1
                        cycle
                    end if
                end if
                nOccOrbs(1) = nOccOrbs(1) + 1
            end do
        end subroutine assignOccOrbs
    End Subroutine DetPreFreezeInit_old

    Subroutine DetPreFreezeInit()

        integer ierr, ms, iEl, flagAlpha
        integer i, j, Lz, OrbOrder(8, 2), FDetTemp(NEl), lmsMax
        type(BasisFn) s
        logical :: tGenFDet
        HElement_t(dp) :: OrbE
        character(25), parameter :: this_routine = 'DetPreFreezeInit'
        allocate (FDet(nEl), stat=ierr)
        LogAlloc(ierr, 'FDet', nEl, 4, tagFDet)
        IF (tDefineDet) THEN
            write (stdout, *) 'Defining FDet according to input'
            do i = 1, NEl
                FDet(i) = DefDet(i)
            end do

            ! A quick check that we have defined a reasonable det.
            ms = sum(get_spin_pn(fdet(1:nel)))
            tRef_Not_HF = .true.
        else
            tGenFDet = .true.
            lmsMax = sum(nOccOrbs) - sum(nClosedOrbs)
            if ((sum(nOccOrbs) + sum(nClosedOrbs)) == nel &
                .and. (.not. TSPN .or. abs(LMS) == lmsMax)) then
                tGenFDet = .false.
                if (LMS < 0) then
                    flagAlpha = 0
                else
                    flagAlpha = 1
                end if
                iEl = 1
                do i = 1, nIrreps
                    ! nClosedOrbs is the number of alpha/beta orbs (the ones with minority spin)
                    do j = 1, nClosedOrbs(i)
                        FDet(iEl) = irrepOrbOffset(i) + 2 * j - flagAlpha
                        iEl = iEl + 1
                    end do
                    ! nOccOrbs is the number of majority spin orbs (per irrep)
                    do j = 1, nOccOrbs(i)
                        FDet(iEl) = irrepOrbOffset(i) + 2 * j - (1 - flagAlpha)
                        iEl = iEl + 1
                    end do
                end do
                call sort(FDet)
            end if
            if (tGenFDet) then
                CALL GENFDET(FDET)
                IF (tUEGSpecifyMomentum) THEN
                    write (stdout, *) 'Defining FDet according to a momentum input'
                    CALL ModifyMomentum(FDET)
                end if
                tRef_Not_HF = .false.
            end if
        end if
        write (stdout, "(A)", advance='no') " Fermi det (D0):"
        call write_det(stdout, FDET, .true.)
        Call GetSym(FDet, nEl, G1, nBasisMax, s)
        write (stdout, "(A)", advance='no') " Symmetry: "
        Call WriteSym(6, s%Sym, .true.)
        IF (tFixLz) THEN
            Call GetLz(FDet, nEl, Lz)
            write (stdout, "(A,I5)") "Lz of Fermi det:", Lz
        end if
        CALL NECI_ICOPY(NEL, FDET, 1, NUHFDET, 1)
        if (tMolpro) then
            !Orbitals are ordered by occupation number from MOLPRO, and not reordered in NECI
            !Therefore, we know HF determinant is first four occupied orbitals.
            write (stdout, "(A)") "NECI called from MOLPRO, so assuming orbitals ordered by occupation number."
            if (.not. tDefineDet) then
                FDetTemp(:) = FDet(:)
            else
                !We have defined our own reference determinant, but still use the first orbitals for the calculation
                !of 'orbital energies'
                CALL GENFDET(FDETTEMP)
            end if
            write (stdout, "(A)") "Calculating orbital energies..."
            do i = 1, nBasis
                OrbE = CalcFockOrbEnergy(i, FDetTemp)
                Arr(i, 1) = real(OrbE, dp)
                Brr(i) = i
            end do
            write (stdout, "(A)") "Reordering basis by orbital energies..."
            OrbOrder(:, :) = 0
            call ORDERBASIS(NBASIS, Arr, Brr, OrbOrder, nBasisMax, G1)
            !However, we reorder them here
            call writebasis(6, G1, nBasis, Arr, Brr)
        end if
        E0HFDET = ECORE
        DO I = 1, NEL
            E0HFDET = E0HFDET + ARR(NUHFDET(i), 2)
        end do
        write (stdout, *) "Fock operator energy:", E0HFDET

        ! Store the value of Ms for use in other areas
        calculated_ms = sum(get_spin_pn(fdet(1:nel)))

    End Subroutine DetPreFreezeInit


    function get_helement_compat(nI, nJ, IC, iLutI, iLutJ) result(hel)
        ! Get the matrix element of the hamiltonian. This assumes that we
        ! already know IC. We do not need to know iLutI, iLutJ (although
        ! they are helpful). This better fits the requirements of existing
        ! code than get_helement_normal.
        !
        ! In:  nI, nJ       - The determinants to consider
        !      iLutI, iLutJ - Bit representations of I,J (optional, helpful)
        !      IC           - The number of orbitals I,J differ by
        ! Ret: hel          - The desired matrix element.

        integer, intent(in) :: nI(nel), nJ(nel)
        integer(kind=n_int), intent(in), optional :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(in) :: IC
        HElement_t(dp) :: hel

        character(*), parameter :: this_routine = 'get_helement_compat'

        integer :: temp_ic
        ! GUGA implementation:

        if (tGUGA) then
            if (all(nI == nJ)) then
                hel = calcDiagMatEleGUGA_nI(nI)
            else
                call stop_all(this_routine, "TODO: refactor guga matrix elements!")
            end if
            return
        end if

        if (tHPHFInts) &
            call stop_all(this_routine, "Should not be calling HPHF &
                          &integrals from here.")

        ! nobody actually uses Simons old CSF implementations..
        if (t_lattice_model) then
            temp_ic = ic
            hel = get_helement_lattice(nI, nJ, temp_ic)
            return
        end if

        if (tStoreAsExcitations) &
            call stop_all(this_routine, "tStoreExcitations not supported")

        if (present(iLutJ)) then
            hel = sltcnd_knowIC(nI, iLutI, iLutJ, IC)
        else
            hel = sltcnd_compat(nI, nJ, IC)
        end if

        ! Add in ECore if for a diagonal element
        if (IC == 0) then
            hel = hel + (ECore)
        else if (tStoquastize) then
            hel = -abs(hel)
        end if

    end function

    function get_helement_normal(nI, nJ, iLutI, iLutJ, ICret) result(hel)

        ! Get the matrix element of the hamiltonian.
        !
        ! In:  nI, nJ       - The determinants to consider
        !      iLutI, iLutJ - Bit representations of I,J (optional, helpful)
        ! Out: ICret        - The number of orbitals I,J differ by
        ! Ret: hel          - The desired matrix element.

        integer, intent(in) :: nI(nel), nJ(nel)
        integer(kind=n_int), intent(in), optional :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(out), optional :: ICret
        HElement_t(dp) :: hel

        character(*), parameter :: this_routine = 'get_helement_normal'
        integer :: ex(2, 2), IC
        integer(kind=n_int) :: ilut(0:NIfTot, 2)

        if (tGUGA) then

            if (all(nI == nJ)) then
                hel = calcDiagMatEleGUGA_nI(nI)
            else
                call stop_all(this_routine, "TODO: refactor guga matrix elements!")
            end if
            return
        end if

        if (tHPHFInts) &
            call stop_all(this_routine, "Should not be calling HPHF &
                          &integrals from here.")

        if (t_lattice_model) then
            if (present(ICret)) then
                ic = -1
                hel = get_helement_lattice(nI, nJ, ic)
                ICret = ic
            else
                hel = get_helement_lattice(nI, nJ)
            end if
            return
        end if

        if (tStoreAsExcitations .and. nI(1) == -1 .and. nJ(1) == -1) then

            ex(1, :) = nJ(4:5)
            ex(2, :) = nJ(6:7)
            hel = sltcnd_excit(nI, Excite_2_t(ex), .false.)
        end if

        if (present(iLutJ)) then
            hel = sltcnd(nI, iLutI, iLutJ, IC)
        else
            call EncodeBitDet(nI, iLut(:, 1))
            call EncodeBitdet(nJ, iLut(:, 2))
            ! TODO: This is not an ideal place to end up...
            hel = sltcnd(nI, iLut(:, 1), ilut(:, 2), IC)
        end if

        ! Add in ECore for a diagonal element
        if (IC == 0) then
            hel = hel + (ECore)
        else if (tStoquastize) then
            hel = -abs(hel)
        end if

        ! If requested, return IC
        if (present(ICret)) then
            ICret = IC
        end if

    end function get_helement_normal

    function get_helement_excit(NI, NJ, IC, ExcitMat, TParity) result(hel)

        ! Calculate the Hamiltonian Matrix Element for a given determinant (or
        ! csf), when we have the excitation matrix and parity of the
        ! excitation.
        !
        ! In:  nI, nJ       - The determinants to evaluate
        !      IC           - The number of orbitals I,J differ by
        !      ex           - The excitation matrix
        !      tParity      - The parity of the excitation
        ! Ret: hel          - The H matrix element

        integer, intent(in) :: nI(nel), nJ(nel), IC
        integer, intent(in) :: ExcitMat(2, ic)
        logical, intent(in) :: tParity
        HElement_t(dp) :: hel

        character(*), parameter :: this_routine = 'get_helement_excit'

        ! intermediately put the special call to the hubbard matrix elements
        ! here. Although I want to change that in the whole code to have
        ! procedure pointers similar to the excitation generator, which gets
        ! intialized to the correct function at the beginning of the
        ! excitations
        ! store all the lattice model matrix elements in one call.
        if (t_lattice_model) then
            hel = get_helement_lattice(nI, ic, ExcitMat, tParity)
            return
        end if

        ! GUGA implementation:
        if (tGUGA) then
            if (all(nI == nJ)) then
                hel = calcDiagMatEleGUGA_nI(nI)
                return
            end if

        end if

        if (IC < 0) &
            call stop_all(this_routine, "get_helement_excit should only be &
                         &used if we know the number of excitations and the &
                         &excitation matrix")

        hel = dyn_sltcnd_excit_old(nI, IC, ExcitMat, tParity)

        if (IC == 0) then
            hel = hel + (ECore)
        else if (tStoquastize) then
            hel = -abs(hel)
        end if

    end function get_helement_excit

    function get_helement_det_only(nI, nJ, iLutI, iLutJ, ic, ex, tParity, &
                                   HElGen) result(hel)

        ! Calculate the Hamiltonian Matrix Element for a determinant as above.
        ! This function assumes that we have got it correct for determinants
        ! (i.e. no error checking), and no conditionals. It also has extra
        ! arguments for compatibility with the function pointer methods.
        !
        ! In:  nI, nJ       - The determinants to evaluate
        !      iLutI, iLutJ - Bit representations (unused)
        !      ic           - The number of orbitals i,j differ by
        !      ex           - Excitation matrix
        !      tParity      - Parity of the excitation
        ! Ret: hel          - The H matrix element

        integer, intent(in) :: nI(nel), nJ(nel), ic, ex(2, ic)
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        logical, intent(in) :: tParity
        HElement_t(dp) :: hel
        HElement_t(dp), intent(in) :: HElGen    !Not used - here for compatibility with other interfaces.
        character(*), parameter :: this_routine = "get_helement_det_only"

        unused_var(ilutJ); unused_var(ilutI); unused_var(nJ); unused_var(hElgen);
        ! GUGA implementation:
        if (tGUGA) then
            if (all(nI == nJ)) then
                hel = calcDiagMatEleGUGA_nI(nI)
            else
                call stop_all(this_routine, "TODO: refactor guga matrix elements!")
            end if
            return
        end if

        ! switch to lattice matrix element:
        if (t_lattice_model) then
            hel = get_helement_lattice(nI, ic, ex, tParity)
            return
        end if

        hel = dyn_sltcnd_excit_old(nI, IC, ex, tParity)

        if (IC == 0) then
            hel = hel + ECore
        else if (tStoquastize) then
            hel = -abs(hel)
        end if
    end function

    HElement_t(dp) function GetH0Element4(nI, HFDet)
        ! Returns the matrix element of the unperturbed Hamiltonian,
        ! which is just the sum of the eigenvalues of the occupied
        ! orbitals and the core energy.
        ! HOWEVER, this routine is *SLOWER* than the GetH0Element3
        ! routine, and should only be used if you require this
        ! without reference to the fock eigenvalues.
        ! This is calculated by subtracting the required two electron terms
        ! from the diagonal matrix elements.
        integer, intent(in) :: nI(NEl), HFDet(NEl)
        HElement_t(dp) :: hel

        hel = SumFock(nI, HFDet)
        GetH0Element4 = hel + ECore

    end function GetH0Element4

    HElement_t(dp) function GetH0Element3(nI)
        ! Wrapper for GetH0Element.
        ! Returns the matrix element of the unperturbed Hamiltonian, which is
        ! just the sum of the eigenvalues of the occupied orbitals and the core
        ! energy.
        !  Note that GetH0Element{1,2} don't exist. The name is to be
        !  consistent with GetHElement3, i.e. offer the most abstraction possible.
        ! In:
        !    nI(nEl)  list of occupied spin orbitals in the determinant.
        integer nI(nEl)
        HElement_t(dp) hEl
        call GetH0Element(nI, nEl, Arr(1:nBasis, 1:2), nBasis, ECore, hEl)
        GetH0Element3 = hEl
    end function

    Subroutine DetCleanup()
    End Subroutine DetCleanup


    SUBROUTINE GENFDET(FDET)
        integer, intent(out) :: fdet(nel)
        logical :: invalid

! Working space
        integer NELS(2)
        integer NELR, IREAL, IS, NSBASIS, I, totSpin, nclosed
        integer :: nspins(2)
        integer, dimension(3) :: totK
        logical :: sr

! eigenvector N is in FMAT(i,N,IS), where i is the component of the
! vector unless transposed
        NSBASIS = NBASIS / 2
        NELS(2) = (LMS + NEL) / 2
        NELS(1) = NEL - NELS(2)

        if (tGUGA) then
            write(stdout, '(4(A,'//int_fmt(nel)//'))') "N_neg:", NELS(1), " ; N_pos:", NELS(2), " ; S:", LMS, " ; nEl:", nEL
        else
            WRITE(stdout, '(4(A,'//int_fmt(nel)//'))') " N_alpha:", NELS(1), " ; N_beta:", NELS(2), " ; LMS:", LMS, " ; NEl:", NEL
        end if
        nclosed = minval(nels)

        if (tMolpro) then
!Assume that orbitals are ordered by occupation number
!Assuming LMS positive, NELS(2) > NELS(1)
            do i = 1, nclosed * 2
                FDET(i) = i  !Fill up closed shells
            end do
            do i = 1, abs(LMS) !Fill up open shells in odd orbitals
                FDET(nclosed * 2 + i) = nclosed * 2 + (i * 2 - (1 - sign(1, LMS)) / 2)
            end do
            return
        end if
        NEL = 0
        totK = 0
        nspins = 0
        totSpin = 0
        sr = .false.
        if (tSymSet) then
! here, a momentum eigenstate for a given total momentum is constructed
            do ireal = 1, nbasis
! invalid is true if the desired total momentum cannot be reached anymore
                invalid = checkMomentumInvalidity(ireal, totK, SymRestrict%k, nels(1) - nspins(1), nels(2) - nspins(2))
                if (invalid) cycle
! is ==1 if the spin is -1 and is==2 if it is 1
                is = (3 + (G1(brr(ireal))%Ms)) / 2
                if (nspins(is) <= nels(is)) then
                    nel = nel + 1
                    nspins(is) = nspins(is) + 1
                    fdet(nel) = brr(ireal)
                    if (t_k_space_hubbard) then
                        totK = lat%add_k_vec(totK, G1(brr(ireal))%k)
                    else
                        totK = totK + G1(brr(ireal))%k
                    end if
                    totSpin = totSpin + G1(brr(ireal))%Ms
                    if (tHub .and. .not. t_k_space_hubbard) then
                        call MomPbcSym(totK, nBasisMax)
                    end if
                end if
            end do
        else
            DO IS = 1, 2
                IREAL = 1
                DO I = 1, NSBASIS
                    if (ireal > nbasis) exit
                    DO WHILE (G1(BRR(IREAL))%Ms /= (-3 + 2 * IS))
                        IREAL = IREAL + 1
                        if (ireal > nbasis) exit
                    END DO
                    NELR = (BRR(IREAL) - 1) / 2 + 1
                    IF (I <= NELS(IS)) THEN
                        NEL = NEL + 1
                        FDET(NEL) = BRR(IREAL)
                    END IF
                    IREAL = IREAL + 1
                END DO
            END DO
        end if
        call sort(fDet(1:nel))
#ifdef DEBUG_
        print *, "Total momentum", totK
        print *, "Requested momentum", SymRestrict%k
        print *, "Total spin", totSpin, "Requested spin", lms
        print *, "fDet", fdet
#endif
    END


    subroutine GetH0Element(nI, nEl, Arr, nBasis, ECore, hEl)
        !  Get a matrix element of the unperturbed Hamiltonian.  This is just
        !  the sum of the Hartree-Fock eigenvalues and the core energy.
        !  In:
        !     nI(nEl)  list of occupied spin orbitals in the determinant.
        !     nEl      # of electrons.
        !     Arr      array containing the eigenvalues of the spin-orbitals.
        !              (See System for how it's defined/stored).
        !     nBasis   # spin orbitals.
        !     ECore    Core energy.
        !  Out:
        !     hEl      <D_i|H_0|D_i>, the unperturbed Hamiltonian matrix element.
        integer nEl, nI(nEl), nBasis
        HElement_t(dp) hEl
        real(dp) Arr(nBasis, 2), ECore
        integer i
        if (tStoreAsExcitations .and. nI(1) == -1) then
    !The excitation storage starts with -1.  The next number is the excitation level,L .
    !Next is the parity of the permutation required to lineup occupied->excited.  Then follows
    !a list of the indexes of the L occupied orbitals within the HFDET, and then L virtual spinorbitals.
            hEl = 0.0_dp
            do i = 4, nI(2) + 4 - 1
                hEl = hEl - (Arr(nI(i), 2))
            end do
            do i = i, i + nI(2) - 1
                hEl = hEl + (Arr(nI(i), 2))
            end do
        else
            hEl = ECore
            do i = 1, nEl
                hEl = hEl + (Arr(nI(i), 2))
            end do
        end if
    end subroutine


    subroutine DetFreezeBasis(GG)
        integer i, j
        integer GG(*), Lz
        type(BasisFn) s
        character(*), parameter :: this_routine = 'DetFreezeBasis'
    !C.. Deal with FDET
    !C.. GG(I) is the new position in G of the (old) orb I
        IF (FDET(1) /= 0) THEN
            J = 0
            DO I = 1, NEL
                FDET(I) = GG(FDET(I))
    !C.. any orbitals which no longer exist, we move outside the basis
                IF (FDET(I) == 0) THEN
                    FDET(I) = nBasis + 1
                ELSE
                    J = J + 1
                end if
            end do
            call sort(fdet)
            IF (J /= NEL - NFROZEN - NFROZENIN) THEN
                write (stdout, *) "Failed Freezing Det:"
                call write_det(stdout, FDET, .true.)
                call stop_all(this_routine, "After Freezing, FDET has wrong number of electrons")
            end if
        end if
        IF (nUHFDet(1) /= 0) THEN
            J = 0
            DO I = 1, NEL
                nUHFDET(I) = GG(nUHFDET(I))
    !C.. any orbitals which no longer exist, we move outside the basis
                IF (nUHFDET(I) == 0) THEN
                    nUHFDET(I) = nBasis + 1
                ELSE
                    J = J + 1
                end if
            end do
            call sort(nUHFDet(1:nel))
            IF (J /= NEL - NFROZEN - NFROZENIN) THEN
                write (stdout, *) "Failed Freezing Det:"
                call write_det(stdout, nUHFDET, .true.)
                call stop_all(this_routine, "After Freezing, UHFDET has wrong number of electrons")
            end if
        end if
        write (stdout, "(A)", advance='no') " Post-Freeze Fermi det (D0):"
        call write_det_len(6, fDet, nel - nfrozen - nfrozenin, .true.)
        write (stdout, "(A)", advance='no') " Symmetry: "
        Call GetSym(FDet, nEl - nFrozen - nFrozenIn, G1, nBasisMax, s)
        Call WriteSym(6, s%Sym, .true.)
        IF (tFixLz) THEN
            Call GetLz(FDet, nEl - nFrozen - nFrozenIn, Lz)
            write (stdout, "(A,I5)") " Lz of Fermi det:", Lz
        end if

    end subroutine

    logical FUNCTION ISUHFDET(NI, NEL)
        integer NEL, NI(NEL)
        integer I
        IF (.NOT. TUSEBRILLOUIN) THEN
            ISUHFDET = .FALSE.
            RETURN
        end if
        ISUHFDET = .TRUE.
        DO I = NEL, 1, -1
            IF (NI(I) /= NUHFDET(I)) THEN
                ISUHFDET = .FALSE.
                EXIT
            end if
        end do
    !         ISUHFDET=.FALSE.
        RETURN
    END Function


    ! Calculate the one-electron part of the energy of a det
    FUNCTION CALCT(NI, NEL)
        integer, intent(in) :: nI(nEl), NEL
        HElement_t(dp) :: CALCT
        calct = sum(gettmatel(nI, nI))
    END function

    !Write out the determinant in bit notation
    SUBROUTINE WriteDetBit(nUnit, iLutnI, lTerm)
        integer, intent(in) :: nUnit
        integer(kind=n_int), intent(in) :: iLutnI(0:NIfTot)
        logical, intent(in) :: lTerm
        integer :: i

        do i = 1, nBasis - 1
            if (IsOcc(iLutnI, i)) then
                write (nUnit, "(A1)", advance='no') "1"
            else
                write (nUnit, "(A1)", advance='no') "0"
            end if
        end do
        if (IsOcc(iLutnI, nBasis)) then
            if (lTerm) then
                write (nUnit, "(A1)") "1"
            else
                write (nUnit, "(A1)", advance='no') "1"
            end if
        else
            if (lTerm) then
                write (nUnit, "(A1)") "0"
            else
                write (nUnit, "(A1)", advance='no') "0"
            end if
        end if

    END SUBROUTINE WriteDetBit


    ! This takes a the ground state FDet generated for the UEG and changes its total momentum
    ! According to input options
    subroutine ModifyMomentum(FDet)
        integer :: i, j ! Loop variables
        integer, intent(inout) :: FDet(NEl)
        integer :: k_total(3) ! Stores the total momentum of FDet
        integer :: delta_k(3) ! Stores the difference between the current FDet and the FDet we're aiming for
        integer, allocatable :: kPointToBasisFn(:, :, :, :) ! Look up table for kPoints to basis functions
        integer :: kmaxX, kminX, kmaxY, kminY, kmaxZ, kminZ, iSpinIndex ! Stores the limits of kPointToBasisFn
        integer :: det_sorted(NEl), e_store ! Storage for the sorting routine
        logical :: sorted ! As above
        integer :: wrapped_index
        integer :: k_new

        IF (.not. tUEG) call stop_all("ModifyMomentum", "Only works for UEG")

        ! Finds current momentum, and finds the difference between this and the target momentum
        ! most commonly will be zero to start with
        k_total(1) = 0
        k_total(2) = 0
        k_total(3) = 0
        do j = 1, NEl
            k_total(1) = k_total(1) + G1(FDet(j))%k(1)
            k_total(2) = k_total(2) + G1(FDet(j))%k(2)
            k_total(3) = k_total(3) + G1(FDet(j))%k(3)
        end do
        delta_k = k_momentum - k_total

        if (delta_k(1) == 0 .and. delta_k(2) == 0 .and. delta_k(3) == 0) write (stdout, *) "WARNING: specified momentum is ground state"

        ! Creates a look-up table for k-points (this was taken from symrandexcit2.F90)
        kmaxX = 0
        kminX = 0
        kmaxY = 0
        kminY = 0
        kminZ = 0
        kmaxZ = 0
        do i = 1, nBasis
            IF (G1(i)%k(1) > kmaxX) kmaxX = G1(i)%k(1)
            IF (G1(i)%k(1) < kminX) kminX = G1(i)%k(1)
            IF (G1(i)%k(2) > kmaxY) kmaxY = G1(i)%k(2)
            IF (G1(i)%k(2) < kminY) kminY = G1(i)%k(2)
            IF (G1(i)%k(3) > kmaxZ) kmaxZ = G1(i)%k(3)
            IF (G1(i)%k(3) < kminZ) kminZ = G1(i)%k(3)
        end do
        allocate (kPointToBasisFn(kminX:kmaxX, kminY:kmaxY, kminZ:kmaxZ, 2))
        do i = 1, nBasis
            iSpinIndex = (G1(i)%Ms + 1) / 2 + 1 ! iSpinIndex equals 1 for a beta spin (ms=-1), and 2 for an alpha spin (ms=1)
            kPointToBasisFn(G1(i)%k(1), G1(i)%k(2), G1(i)%k(3), iSpinIndex) = i
        end do

        ! For each of the three dimensions, nudge electrons one at a time by one momentum unit until delta_k is reached
        det_sorted = FDet

        ! Bubble sort to order det_sorted in order of kx of the corresponding electron
        do
            sorted = .true.
            do i = 1, NEl - 1
                j = i + 1
                if (G1(det_sorted(j))%k(1) > G1(det_sorted(i))%k(1)) then
                    sorted = .false.
                    e_store = det_sorted(i)
                    det_sorted(i) = det_sorted(j)
                    det_sorted(j) = e_store
                end if
            end do
            if (sorted) exit
        end do

        ! Nudge momenta one at a time
        if (delta_k(1) > 0) then
            do i = 1, delta_k(1)
                wrapped_index = mod(i, NEl)        ! Take the modulus to know which electron to nudge
                if (wrapped_index == 0) then    ! Deal with the i=NEl case
                    wrapped_index = NEl
                end if
                j = wrapped_index                 ! For convenience asign this to j
                k_new = G1(det_sorted(j))%k(1) + 1  ! Find the new momentum of this electron
                if (k_new > kmaxX) then        ! Check that this momentum isn't outside the cell
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                end if
                iSpinIndex = (G1(j)%Ms + 1) / 2 + 1     ! Spin of the new orbital is the same as the old
    ! Finds basis number for the new momentum
                det_sorted(j) = kPointToBasisFn(k_new, G1(det_sorted(j))%k(2), G1(det_sorted(j))%k(3), iSpinIndex)
            end do
        else if (delta_k(1) < 0) then ! For the negative case, i must run through negative numbers
            do i = -1, delta_k(1), -1
                wrapped_index = mod(i, NEl)
                if (wrapped_index == 0) then
                    wrapped_index = -NEl
                end if
                j = NEl + wrapped_index + 1 ! Now this goes through the list backward (wrapped_index is negative)
                k_new = G1(det_sorted(j))%k(1) - 1 ! Find the new momentum of this electron, this time in the opposite direction
                if (k_new < kminX) then ! Check the limits of the cell again
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                end if
                iSpinIndex = (G1(j)%Ms + 1) / 2 + 1 ! Spin of the new orbital is the same as the old
    ! Finds basis number for the new momentum
                det_sorted(j) = kPointToBasisFn(k_new, G1(det_sorted(j))%k(2), G1(det_sorted(j))%k(3), iSpinIndex)
            end do
        end if

        FDet = det_sorted

        !====ky treated as kx above
        do
            sorted = .true.
            do i = 1, NEl - 1
                j = i + 1
                if (G1(det_sorted(j))%k(2) > G1(det_sorted(i))%k(2)) then
                    sorted = .false.
                    e_store = det_sorted(i)
                    det_sorted(i) = det_sorted(j)
                    det_sorted(j) = e_store
                end if
            end do
            if (sorted) exit
        end do

        ! Nudge momenta one at a time
        if (delta_k(2) > 0) then
            do i = 1, delta_k(2)
                wrapped_index = mod(i, NEl)        ! Take the modulus to know which electron to nudge
                if (wrapped_index == 0) then    ! Deal with the i=NEl case
                    wrapped_index = NEl
                end if
                j = wrapped_index                 ! For convenience asign this to j
                k_new = G1(det_sorted(j))%k(2) + 1  ! Find the new momentum of this electron
                if (k_new > kmaxY) then        ! Check that this momentum isn't outside the cell
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                end if
                iSpinIndex = (G1(j)%Ms + 1) / 2 + 1     ! Spin of the new orbital is the same as the old
    ! Finds basis number for the new momentum
                det_sorted(j) = kPointToBasisFn(G1(det_sorted(j))%k(1), k_new, G1(det_sorted(j))%k(3), iSpinIndex)
            end do
        else if (delta_k(2) < 0) then ! For the negative case, i must run through negative numbers
            do i = -1, delta_k(2), -1
                wrapped_index = mod(i, NEl)
                if (wrapped_index == 0) then
                    wrapped_index = -NEl
                end if
                j = NEl + wrapped_index + 1 ! Now this goes through the list backward (wrapped_index is negative)
                k_new = G1(det_sorted(j))%k(2) - 1 ! Find the new momentum of this electron, this time in the opposite direction
                if (k_new < kminY) then ! Check the limits of the cell again
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                end if
                iSpinIndex = (G1(j)%Ms + 1) / 2 + 1 ! Spin of the new orbital is the same as the old
    ! Finds basis number for the new momentum
                det_sorted(j) = kPointToBasisFn(G1(det_sorted(j))%k(1), k_new, G1(det_sorted(j))%k(3), iSpinIndex)
            end do
        end if

        FDet = det_sorted

        !====kz treated as kx and ky above
        do
            sorted = .true.
            do i = 1, NEl - 1
                j = i + 1
                if (G1(det_sorted(j))%k(3) > G1(det_sorted(i))%k(3)) then
                    sorted = .false.
                    e_store = det_sorted(i)
                    det_sorted(i) = det_sorted(j)
                    det_sorted(j) = e_store
                end if
            end do
            if (sorted) exit
        end do

        ! Nudge momenta one at a time
        if (delta_k(3) > 0) then
            do i = 1, delta_k(3)
                wrapped_index = mod(i, NEl)        ! Take the modulus to know which electron to nudge
                if (wrapped_index == 0) then    ! Deal with the i=NEl case
                    wrapped_index = NEl
                end if
                j = wrapped_index                 ! For convenience asign this to j
                k_new = G1(det_sorted(j))%k(3) + 1  ! Find the new momentum of this electron
                if (k_new > kmaxZ) then        ! Check that this momentum isn't outside the cell
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                end if
                iSpinIndex = (G1(j)%Ms + 1) / 2 + 1     ! Spin of the new orbital is the same as the old
    ! Finds basis number for the new momentum
                det_sorted(j) = kPointToBasisFn(G1(det_sorted(j))%k(1), G1(det_sorted(j))%k(2), k_new, iSpinIndex)
            end do
        else if (delta_k(3) < 0) then ! For the negative case, i must run through negative numbers
            do i = -1, delta_k(3), -1
                wrapped_index = mod(i, NEl)
                if (wrapped_index == 0) then
                    wrapped_index = -NEl
                end if
                j = NEl + wrapped_index + 1 ! Now this goes through the list backward (wrapped_index is negative)
                k_new = G1(det_sorted(j))%k(3) - 1 ! Find the new momentum of this electron, this time in the opposite direction
                if (k_new < kminZ) then ! Check the limits of the cell again
                    call stop_all("ModifyMomentum", "Electron moved outside of the cell limits")
                end if
                iSpinIndex = (G1(j)%Ms + 1) / 2 + 1 ! Spin of the new orbital is the same as the old
                ! Finds basis number for the new momentum
                det_sorted(j) = kPointToBasisFn(G1(det_sorted(j))%k(1), G1(det_sorted(j))%k(2), k_new, iSpinIndex)
            end do
        end if

        FDet = det_sorted

        ! Bubble sort to order FDet back into increasing order by number
        do
            sorted = .true.
            do i = 1, NEl - 1
                j = i + 1
                if (FDet(j) < FDet(i)) then
                    sorted = .false.
                    e_store = FDet(i)
                    FDet(i) = FDet(j)
                    FDet(j) = e_store
                end if
            end do
            if (sorted) exit
        end do

        write (stdout, *) "Total momentum set to", k_momentum

    end subroutine ModifyMomentum


    SUBROUTINE ORDERBASIS(NBASIS, ARR, BRR, ORBORDER, NBASISMAX, G1)
        INTEGER NBASIS, BRR(NBASIS), ORBORDER(8, 2), nBasisMax(5, *)
        INTEGER BRR2(NBASIS)
        TYPE(BASISFN) G1(NBASIS)
        real(dp) ARR(NBASIS, 2), ARR2(NBASIS, 2)
        INTEGER IDONE, I, J, IBFN, ITOT, ITYPE, ISPIN
        real(dp) OEN
        character(*), parameter :: this_routine = 'ORDERBASIS'
        IDONE = 0
        ITOT = 0
    !.. copy the default ordered energies.
        CALL DCOPY(NBASIS, ARR(1, 1), 1, ARR(1, 2), 1)
        CALL DCOPY(NBASIS, ARR(1, 1), 1, ARR2(1, 2), 1)
        write(stdout, *) ''
        write(stdout, "(A,8I4)") "Ordering Basis (Closed): ", (ORBORDER(I, 1), I=1, 8)
        write(stdout, "(A,8I4)") "Ordering Basis (Open  ): ", (ORBORDER(I, 2), I=1, 8)
        IF (NBASISMAX(3, 3) == 1) THEN
    !.. we use the symmetries of the spatial orbitals
    ! actually this is never really used below here it seems.. since orborder
    ! is only zeros, according to output. check that!
    ! and that is independent of the GUGA implementation TODO: check orborder!
            DO ITYPE = 1, 2
                IBFN = 1
                DO I = 1, 8
                    ! 8 probably because at most D2h symmetry giovanni told me about.
                    DO J = 1, ORBORDER(I, ITYPE)
                        DO WHILE (IBFN <= NBASIS .AND. (G1(IBFN)%SYM%s < I - 1 .OR. BRR(IBFN) == 0))
                            IBFN = IBFN + 1
                        end do
                        IF (IBFN > NBASIS) THEN
                            call stop_all(this_routine, "Cannot find enough basis fns of correct symmetry")
                        end if
                        IDONE = IDONE + 1
                        BRR2(IDONE) = IBFN
                        BRR(IBFN) = 0
                        ARR2(IDONE, 1) = ARR(IBFN, 1)
                        IBFN = IBFN + 1
                    end do
                end do
                ! Beta sort
                call sort(arr2(itot + 1:idone, 1), brr2(itot + 1:idone), nskip=2)
                ! Alpha sort
                call sort(arr2(itot + 2:idone, 1), brr2(itot + 2:idone), nskip=2)
                ITOT = IDONE
            end do
            DO I = 1, NBASIS
                IF (BRR(I) /= 0) THEN
                    ITOT = ITOT + 1
                    BRR2(ITOT) = BRR(I)
                    ARR2(ITOT, 1) = ARR(I, 1)
                end if
            end do
            ! what are those doing?
            ! ok those are copying the newly obtained arr2 and brr2 into arr and brr
            CALL NECI_ICOPY(NBASIS, BRR2, 1, BRR, 1)
            CALL DCOPY(NBASIS, ARR2(1, 1), 1, ARR(1, 1), 1)
        end if
        ! i think this is the only reached point: and this means i can make it
        ! similar to the Hubbard implementation to not reorder!
    ! beta sort
        call sort(arr(idone + 1:nbasis, 1), brr(idone + 1:nbasis), nskip=2)
    ! alpha sort
        call sort(arr(idone + 2:nbasis, 1), brr(idone + 2:nbasis), nskip=2)
    !.. We need to now go through each set of degenerate orbitals, and make
    !.. the correct ones are paired together in BRR otherwise bad things
    !.. happen in FREEZEBASIS
    !.. We do this by ensuring that within a degenerate set, the BRR are in
    !.. ascending order
    !         IF(NBASISMAX(3,3).EQ.1) G1(3,BRR(1))=J
        DO ISPIN = 0, 1
            OEN = ARR(1 + ISPIN, 1)
            J = 1 + ISPIN
            ITOT = 2
            DO I = 3 + ISPIN, NBASIS, 2
                IF (ABS(ARR(I, 1) - OEN) > 1.0e-4_dp) THEN
    !.. We don't have degenerate orbitals
    !.. First deal with the last set of degenerate orbitals
    !.. We sort them into order of BRR
                    call sort(brr(i - itot:i - 1), arr(i - itot:i - 1, 1), nskip=2)
    !.. now setup the new degenerate set.
                    J = J + 2
                    ITOT = 2
                ELSE
                    ITOT = ITOT + 2
                end if
                OEN = ARR(I, 1)
                IF (NBASISMAX(3, 3) == 1) THEN
    !.. If we've got a generic spatial sym or hf we mark degeneracies
    !               G(3,BRR(I))=J
                end if
            end do
    ! i is now nBasis+2
            call sort(brr(i - itot:i - 2), arr(i - itot:i - 2, 1), nskip=2)
        end do
    !   if (t_guga_noreorder) then
    !       ! this probably does not work so easy:
    !       allocate(temp_sym(nBasis))
    !       do i = 1, nBasis
    !         temp_sym(i) = G1(i)
    !       end do
    !       do i = 1, nBasis
    !           G1(i) = temp_sym(brr(i))
    !           brr(i) = i
    !       end do
    !       ! could i just do a new molpsymtable here??
    !       ! but only do it if symmetry is not turned off explicetyl!
    !       if (.not. lNoSymmetry) CALL GENMOLPSYMTABLE(NBASISMAX(5,2)+1,G1,NBASIS)
    !   end if

    END subroutine ORDERBASIS

    !dUnscaledEnergy gives the energy without reference to box size and without any offset.
    SUBROUTINE GetUEGKE(I, J, K, ALAT, tUEGTrueEnergies, tUEGOffset, k_offset, Energy, dUnscaledEnergy)
        INTEGER I, J, K
        real(dp) ALat(3), k_offset(3), Energy, E
        LOGICAL tUEGOffset, tUEGTrueEnergies
        real(dp) ::  dUnscaledEnergy
        integer :: kvecX, kvecY, kvecZ
        !==================================
        ! initialize unscaled energy for the case of not using tUEGTrueEnergies
        dunscaledEnergy = 0.0_dp
        if (tUEG2) then
            ! kvectors in cartesian coordinates
            kvecX = k_lattice_vectors(1, 1) * I + k_lattice_vectors(2, 1) * J + k_lattice_vectors(3, 1) * K
            kvecY = k_lattice_vectors(1, 2) * I + k_lattice_vectors(2, 2) * J + k_lattice_vectors(3, 2) * K
            kvecZ = k_lattice_vectors(1, 3) * I + k_lattice_vectors(2, 3) * J + k_lattice_vectors(3, 3) * K

            IF (tUEGTrueEnergies) then
                if (tUEGOffset) then
                    E = (kvecX + k_offset(1))**2 + (kvecY + k_offset(2))**2 + (kvecZ + k_offset(3))**2
                else
                    E = (kvecX)**2 + (kvecY)**2 + (kvecZ)**2
                end if
                Energy = 0.5_dp * E * k_lattice_constant**2
                dUnscaledEnergy = ((kvecX)**2 + (kvecY)**2 + (kvecZ)**2)
            ELSE
                Energy = ((kvecX)**2 + (kvecY)**2 + (kvecZ)**2)
            end if

            return
        end if
        !==================================
        IF (tUEGTrueEnergies) then
            IF (tUEGOffset) then
                E = ((I + k_offset(1))**2 / ALAT(1)**2)
                E = E + ((J + k_offset(2))**2 / ALAT(2)**2)
                E = E + ((K + k_offset(3))**2 / ALAT(3)**2)
            else
                E = (I * I / ALAT(1)**2)
                E = E + (J * J / ALAT(2)**2)
                E = E + (K * K / ALAT(3)**2)
            end if
            Energy = 0.5 * 4 * PI * PI * E
            dUnscaledEnergy = (I * I)
            dUnscaledEnergy = dUnscaledEnergy + (J * J)
            dUnscaledEnergy = dUnscaledEnergy + (K * K)
        ELSE
            E = (I * I)
            E = E + (J * J)
            E = E + (K * K)
            Energy = E
        end if
    END SUBROUTINE GetUEGKE


    SUBROUTINE WRITEBASIS(NUNIT, G1, NHG, ARR, BRR)
        ! Write out the current basis to unit nUnit
        integer, intent(in) :: nunit
        type(basisfn), intent(in) :: g1(nhg)
        integer, intent(in) :: nhg, brr(nhg)

        integer :: pos, i
        real(dp) ARR(NHG, 2), unscaled_energy, kvecX, kvecY, kvecZ

        ! nb. Cannot use EncodeBitDet as would be easy, as nifd, niftot etc are not
        !     filled in yet. --> track pos.
        if (.not. associated(fdet)) &
            write(nunit, '("HF determinant not yet defined.")')
        pos = 1
    !=============================================
        if (tUEG2) then

            DO I = 1, NHG
    !     kvectors in cartesian coordinates
                kvecX = k_lattice_vectors(1, 1) * G1(BRR(I))%K(1) &
                        + k_lattice_vectors(2, 1) * G1(BRR(I))%K(2) &
                        + k_lattice_vectors(3, 1) * G1(BRR(I))%K(3)
                kvecY = k_lattice_vectors(1, 2) * G1(BRR(I))%K(1) &
                        + k_lattice_vectors(2, 2) * G1(BRR(I))%K(2) &
                        + k_lattice_vectors(3, 2) * G1(BRR(I))%K(3)
                kvecZ = k_lattice_vectors(1, 3) * G1(BRR(I))%K(1) &
                        + k_lattice_vectors(2, 3) * G1(BRR(I))%K(2) &
                        + k_lattice_vectors(3, 3) * G1(BRR(I))%K(3)

                unscaled_energy = ((kvecX)**2 + (kvecY)**2 + (kvecZ)**2)

                write(NUNIT, '(6I7)', advance='no') I, BRR(I), G1(BRR(I))%K(1), G1(BRR(I))%K(2), G1(BRR(I))%K(3), G1(BRR(I))%MS
                CALL WRITESYM(NUNIT, G1(BRR(I))%SYM, .FALSE.)
                write(NUNIT, '(I4)', advance='no') G1(BRR(I))%Ml
                write(NUNIT, '(3F19.9)', advance='no') ARR(I, 1), ARR(BRR(I), 2), unscaled_energy

                if (associated(fdet)) then
                    pos = 1
                    do while (pos < nel .and. fdet(pos) < brr(i))
                        pos = pos + 1
                    end do
                    if (brr(i) == fdet(pos)) write(nunit, '(" #")', advance='no')
                end if
                write(nunit, *)
            end do
            RETURN
        end if !UEG2
    !=============================================
        DO I = 1, NHG
            write(NUNIT, '(6I7)', advance='no') I, BRR(I), G1(BRR(I))%K(1), G1(BRR(I))%K(2), G1(BRR(I))%K(3), G1(BRR(I))%MS
            CALL WRITESYM(NUNIT, G1(BRR(I))%SYM, .FALSE.)
            write(NUNIT, '(I4)', advance='no') G1(BRR(I))%Ml
            write(NUNIT, '(2F19.9)', advance='no') ARR(I, 1), ARR(BRR(I), 2)
            if (associated(fdet)) then
                pos = 1
                do while (pos < nel .and. fdet(pos) < brr(i))
                    pos = pos + 1
                end do
                if (brr(i) == fdet(pos)) write(nunit, '(" #")', advance='no')
            end if
            write(nunit, *)
        end do
        RETURN
    END SUBROUTINE WRITEBASIS

end module Determinants
