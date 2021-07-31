#include "macros.h"

module Parallel_Calc

!= Algorithms for calculations in parallel.

!= These routines are still in heavy development (i.e. may well contain numerous
!= bugs) and should not be used in production work blindly. :-)

    use Parallel_neci
    use util_mod, only: near_zero, operator(.isclose.)

    implicit none

contains

    subroutine ParMP2(nI)
        != Calculate the MP2 energy in parallel.
        != Designed for CPMD calculations.
        != Parallel used only for distribution of the wavefunctions and to
        != compute the integrals.
        !=
        != In:
        !=    nI(nEl) : list of occupied spin orbitals in the reference
        !=              determinant.
        != Prints out the <D_0|H|D_0>+MP2 energy.

        != To do:
        != 1. Write a custom routine to calculate the <D_0|H|D_1> elements.
        != This would then enable the excitations ij->ab,iJ->aB,Ij-Ab and IJ->AB
        != to be handled at once, removing the need to any integrals to be stored.
        != 2. Implement for standalone NECI calculations (molecular, Hubbard etc).

        use constants, only: dp
        use System, only: AreSameSpatialOrb
        use SystemData, only: nBasisMax, nEl, Beta, ARR, nBasis, ECore, G1, tCPMD, Symmetry, &
                              t_3_body_excits, tGUGA
        use CalcData, only: NWHTAY
        use Integrals_neci, only: GetUMatEl2
        use UMatCache, only: GTID
        use OneEInts, only: GetTMatEl
        Use Determinants, only: get_helement, GetH0Element3
        use global_utilities
        use SymData, only: SymLabels
        use CPMDData, only: KPntInd
        use sym_mod
        use MemoryManager, only: TagIntType
        use neci_intfce
        IMPLICIT NONE
        integer :: nI(nEl)
        integer :: iMinElec, iMaxElec
        integer :: i, j
        integer :: IA, JA, AA, BA, JJ
        integer :: store(6), Excit(2, 2)
        integer :: ic, exlen(1), iC0, ExLevel
        integer, pointer :: Ex(:)
        integer :: nJ(nEl), weight
        HElement_t(dp) dU(2)
        real(dp) :: dE1, dE2
        HElement_t(dp) :: dE, dEtot(2), dE0
        ! MPIHelSum requires arrays as input/output.
        HElement_t(dp) :: dEarr(2)
        type(Symmetry) :: iSym1, iSym2
        type(Symmetry) :: iSym1Conj, iSym2Conj
        logical :: tSign
        integer :: ierr
        integer(TagIntType) :: tag_Ex
        type(timer), save :: proc_timer
        character(*), parameter :: this_routine = 'ParMP2'
        logical :: dbg

        dbg = .false.

        proc_timer%timer_name = 'ParMP2    '
        call set_timer(proc_timer)

        select case (IAND(nWHTay(1, 1), 24))
        case (0)
            ! Generate both single and double excitations.  CPMD works in a Kohn--Sham
            ! basis, and so Brillouin's theorem does not apply.
            write(stdout, *) 'ParMP2: using single and double excitation.'
            ExLevel = 3
        case (8)
            write(stdout, *) 'ParMP2: using only single excitations.'
            ExLevel = 1
        case (16)
            write(stdout, *) 'ParMP2: using only double excitations.'
            ExLevel = 2
        case (24)
            call stop_all('ParMP2', 'Invalid combination of flags in nWHTay.  Invalid EXCITAIONS specification?')
        end select

        iC0 = 0

        write(stdout, *) "Proc ", iProcIndex + 1, "/", nProcessors
        if (tCPMD) then
            ! For CPMD jobs, we actually want each processor to do the full sum, as each
            ! integral is divided across the processors.
            iMinElec = 1
            iMaxElec = nEl
        else
            ! For other calculations, the sum is split over processors.
            call GetProcElectrons(iProcIndex + 1, iMinElec, iMaxElec)
        end if
        write(stdout, *) "Electrons ", iMinElec, " TO ", iMaxElec

!  The root's "energy"---sum over the eigenvalues of occ. spin orbitals and the
!  core energy. Only used it the H0Element formulation is used (see below).
        dE1 = GetH0Element3(nI)

!  Initialize: get the contribution from the reference determinant.
        dE0 = get_helement(nI, nI, 0)

!  Now enumerate all 2v graphs
!  Setup the spin excit generator
        STORE(1) = 0
!  IC is the excitation level (relative to the reverence det).
        CALL GENSYMEXCITIT3Par(NI, .TRUE., EXLEN, nJ, IC, STORE, ExLevel, iMinElec, iMaxElec)
        allocate(Ex(exLen(1)), stat=ierr)
        call LogMemAlloc('Ex', Exlen(1), 4, this_routine, tag_Ex, ierr)
        EX(1) = 0
        CALL GENSYMEXCITIT3Par(NI, .TRUE., EX, nJ, IC, STORE, ExLevel, iMinElec, iMaxElec)

!  Generate the first excitation
        CALL GENSYMEXCITIT3Par(NI, .False., EX, nJ, IC, STORE, ExLevel, iMinElec, iMaxElec)
        i = 0
        j = 0
        dETot = (0.0_dp)

        DO WHILE (NJ(1) /= 0)
! NJ(1) is zero when there are no more excitations.

            i = i + 1

            if (dbg) then
                ! Quickest attempt (and useful for debugging).
                ! Not as efficient as the code below, but clearer and useful for debugging.
                ! Also, this should be used for unrestricted calculations.

                ! MP2 theory refers to the unperturbed excited determinant
                ! => use GetH0Element3 rather than GetHElement3.
                Excit(1, 1) = 2
                !call GetExcitation(nI,nJ,nEl,Excit,tSign)
                dE2 = GetH0Element3(nJ)
                if (tGUGA) then
                    call stop_all("ParMP2", "modify for GUGA")
                end if
                dU(1) = get_helement(nI, nJ, IC)
                !dU(1) = get_helement (nI, nJ, IC, Excit, tSign)
                call getMP2E(dE1, dE2, dU(1), dE)
                dETot(2) = dETot(2) + dE

            else
                ! Efficient (but more complicated) approach.
                ! Reformulate sum to avoid needing integrals more than once.
                ! This doesn't quite succeed: there are approximately N integrals (#
                ! electrons) which are needed twice in the single excitations if you
                ! take k-point symmetry into account.  It is impossible to avoid this
                ! without storing integrals.

                ! Alternatively, calculate the energy of the excited determinant
                ! in reference to that of the reference determinant (i.e. setting dE1=0).
                ASSERT(.not. t_3_body_excits)
                Excit(1, 1) = 2
                call GetExcitation(nI, nJ, nEl, Excit, tSign)

                ! Assuming a restricted calculation.

                ! Multiply the contribuition of the current excitation by weight.
                ! This allows equivalent excitations to be treated by explicitly
                ! calculating only one of them.
                ! If weight remains zero, we don't explicitly calculate the contribution
                ! of the current excitation.
                weight = 0

                if (Excit(1, 2) == 0) then
                    ! Single excitation
                    if (G1(Excit(1, 1))%Ms == -1) then
                        ! alpha -> alpha single excitation.
                        ! count for beta->beta as well.
                        weight = 2
                    end if
                else
                    ! Double excitation.
                    if (G1(Excit(1, 1))%Ms == -1 .and. G1(Excit(1, 2))%Ms == -1) then
                        ! alpha,alpha -> alpha,alpha double excitation.
                        ! count for beta,beta -> beta,beta as well.
                        weight = 2
                    else if (G1(Excit(1, 1))%Ms == -1 .and. G1(Excit(1, 2))%Ms == 1) then
                        ! alpha,beta -> alpha,beta double excitation.
                        ! We consider just these, and bring in the beta,alpha -> beta,alpha
                        ! excitations via spin symmetry.
                        if (AreSameSpatialOrb(Excit(1, 1), Excit(1, 2))) then
                            ! Excitations from the spin orbitals of the same spatial
                            ! orbital.
                            if (AreSameSpatialOrb(Excit(2, 1), Excit(2, 2))) then
                                ! e.g (1a,1b) -> (2a,2b)
                                ! Unique (occurs only once).
                                weight = 1
                            else if (G1(Excit(2, 1))%Ms == -1) then
                                ! e.g (1a,1b) -> (2a,3b)
                                ! Count also for the identical contribution for (1a,1b) -> (2b,3a).
                                weight = 2
                            end if
                        else
                            ! Excitations from different spatial orbitals.  Use spin
                            ! symmetry:
                            ! (1a,2b) -> (3a,3b) has an identical contribution to the MP2
                            ! as (1b,2a) -> (3a,3b).
                            ! Similarly (1a,2b) -> (3a,4b) and (1b,2a) -> (3b,4a).
                            !
                            ! In a restricted calculation, the integrals needed for
                            ! (1a,2b) -> (3a,4b) and (1b,2a) -> (3b,4a) are also
                            ! used for (1a,2a) -> (3a,4a), so we include these
                            ! the contributions from (1a,2b) -> (3a,4b) and
                            ! (1b,2a) -> (3b,4a) when we evaluate the (1a,2a) -> (3a,4a)
                            ! excitation.
                            !
                            ! This leaves only excitations to the same spatial orbital.
                            ! Count also for (1b,2a) -> (3a,3b).
                            if (AreSameSpatialOrb(Excit(2, 1), Excit(2, 2))) weight = 2
                        end if
                    end if
                end if

                IA = GTID(Excit(1, 1))
                AA = GTID(Excit(2, 1))
                if (Excit(1, 2) /= 0) then
                    JA = GTID(Excit(1, 2))
                    BA = GTID(Excit(2, 2))
                end if

                if (tCPMD) then
                    ! Further reduce the number of symmetry-unique excitations by using
                    ! k-point symmetry.
                    iSym1 = SymLabels(KPntInd(IA))
                    iSym1Conj = SymConj(iSym1)
                    if (Excit(1, 2) == 0) then
                        ! Single excitation, i,k_i -> a,k_i.  k_i==k_a.
                        ! If the k_i > -k_i (comparing indices, rather than the k-point),
                        ! then count the current excitation also for  i,-k_i -> a,-k_i,
                        ! which has the same contribution to the MP2 energy.
                        ! If k_i=-k_i (i.e. self-inverse) then weight is unchanged.
                        if (iSym1%s > iSym1Conj%s) then
                            weight = weight * 2
                        else if (iSym1%s < iSym1Conj%s) then
                            weight = 0
                        end if
                    else
                        ! Double excitation.
                        iSym2 = SymLabels(KPntInd(JA))
                        iSym2Conj = SymConj(iSym2)
                        if (iSym1%s == iSym2Conj%s .and. (Arr(Excit(1, 1), 2) .isclose.Arr(Excit(1, 2), 2))) then
                            ! Excitation from, e.g. 1,k1 1,-k1.  Unique: leave weight
                            ! unchanged.
                        else if (iSym1%s > iSym1Conj%s) then
                            ! Count excitation also for the -k equivalent excitation.
                            weight = weight * 2
                        else if (iSym1%s < iSym1Conj%s) then
                            ! Counted for above.
                            weight = 0
                        else if (iSym2%s > iSym2Conj%s) then
                            ! Excitation from (i,k_i j,k_j), where k_i=-k_i.
                            ! If k_j > -k_j, then count also for (i,k_i j,-k_j).
                            weight = weight * 2
                        else if (iSym2%s < iSym2Conj%s) then
                            ! Counted for above.
                            weight = 0
                        end if
                    end if
                end if

                if (weight /= 0) then

                    j = j + 1
                    dE2 = (Arr(Excit(2, 1), 2) - Arr(Excit(1, 1), 2))
                    dU = (0.0_dp)
                    if (Excit(2, 2) == 0) then
                        ! Single excitation.
                        ! dU=\sum_J 2<IJ|AJ>-<IJ|JA> (in terms of spatial orbitals).
                        IA = GTID(Excit(1, 1))
                        AA = GTID(Excit(2, 1))
                        do JJ = 1, nEl, 2 ! Assuming closed shell.  But we have already assumed restricted. ;-)
                            ! Spatial orbital of the j-th element of the reference determinant.
                            JA = GTID(nI(JJ))
                            ! Try to be as efficient as possible with the integrals...
                            ! Want to ask for each integral only once (we don't *quite*
                            ! succeed), so that the sum is efficient even without a cache.
                            ! \sum_j 2<ij|aj> - <ij|ja>
                            if (JA == IA) then
                                dU(1) = dU(1) + GetUMatEl2(IA, JA, AA, JA)
                            else if (tCPMD) then
                                ! Take advantage of k-point symmetry.
                                iSym2 = SymLabels(KPntInd(JA))
                                iSym2Conj = SymConj(iSym2)
                                if (iSym2%s == iSym1%s .or. iSym2Conj%s == iSym1%s) then
                                    dU(1) = dU(1) + (2) * GetUMatEl2(IA, JA, AA, JA) - GetUMatEl2(IA, JA, JA, AA)
                                else if (iSym2%s > iSym2Conj%s) then
                                    ! <i j,k_j | a j,k_j> = <i j,-k_j | a j,-k_j>
                                    ! Count it here to reduce integrals to be evaluated.
                                    dU(1) = dU(1) + (4) * GetUMatEl2(IA, JA, AA, JA) - GetUMatEl2(IA, JA, JA, AA)
                                else if (iSym2%s < iSym2Conj%s) then
                                    ! Already added <ij|ji>.
                                    dU(1) = dU(1) - GetUMatEl2(IA, JA, JA, AA)
                                end if
                            else
                                dU(1) = dU(1) + (2) * GetUMatEl2(IA, JA, AA, JA) - GetUMatEl2(IA, JA, JA, AA)
                            end if
                        end do
                        dU(1) = dU(1) + GetTMATEl(Excit(1, 1), Excit(2, 1))
                    else
                        ! Double excitation
                        dE2 = dE2 + (Arr(Excit(2, 2), 2) - Arr(Excit(1, 2), 2))
                        ! Obtain the <ij|ab> and <ij|ba> integrals as required.  These
                        ! are combined evaulated below for the various contributions to
                        ! the MP2 energy.
                        if (G1(Excit(1, 1))%Ms == G1(Excit(2, 1))%Ms .AND. G1(Excit(1, 2))%Ms == G1(Excit(2, 2))%Ms) then
                            dU(1) = GetUMatEl2(IA, JA, AA, BA)
                        end if
                        if (G1(Excit(1, 1))%Ms == G1(Excit(2, 2))%Ms .AND. G1(Excit(1, 2))%Ms == G1(Excit(2, 1))%Ms) then
                            dU(2) = GetUMatEl2(IA, JA, BA, AA)
                        end if
                    end if

                    if (Excit(1, 2) == 0) then
                        ! Singles contribution.
                        call getMP2E(0.0_dp, dE2, dU(1), dE)
                        dETot(1) = dETot(1) + (weight) * dE
                    else
                        ! Doubles contributions.
                        if (abs(dU(2)) > 0.0_dp) then
                            ! Get e.g. (1a,2b)->(3a,4b) and (1a,2b)->(3b,4a) for "free"
                            ! when we evaluate (1a,2a)->(3a,4a).
                            call getMP2E(0.0_dp, dE2, dU(1), dE)
                            dETot(2) = dETot(2) + (weight) * dE
                            call getMP2E(0.0_dp, dE2, dU(2), dE)
                            dETot(2) = dETot(2) + (weight) * dE
                        end if
                        dU(1) = dU(1) - dU(2)
                        call getMP2E(0.0_dp, dE2, dU(1), dE)
                        dETot(2) = dETot(2) + (weight) * dE
                    end if

                    ! END of more efficient approach.
                end if

            end if

            !write(stdout,'(2i3,a2,2i3,2f17.8)') Excit(1,:),'->',Excit(2,:),dE

            ! Get next excitation.
            CALL GENSYMEXCITIT3Par(NI, .false., EX, nJ, IC, STORE, ExLevel, iMinElec, iMaxElec)

        end do

        write(stdout, *) 'No. of excitations=', I
        write(stdout, *) 'No. of spin and symmetry unique excitations=', J
        if (.not. tCPMD) then
            write(stdout, '(a28,i3,a1,2f15.8)') 'Contribution from processor', iProcIndex + 1, ':', dEtot
            dEarr = dETot
            call MPISumAll(dEArr, 2, dETot)
        end if
        if (iand(ExLevel, 1) == 1) write(stdout, *) 'MP2 SINGLES=', dETot(1) + dE0
        if (iand(ExLevel, 2) == 2) write(stdout, *) 'MP2 DOUBLES=', dETot(2) + dE0
        write(stdout, *) 'MP2 ENERGY =', dETot(1) + dETot(2) + dE0

        deallocate(Ex)
        call LogMemDealloc(this_routine, tag_Ex)

        call halt_timer(proc_timer)

    end subroutine ParMP2

    Subroutine Par2vSum(nI)
        !=  A parallel version of the 2-vertex sum.
        !=
        !=   In:
        !=     nI(nEl)        The root determinant of the 2v sum.
        !=
        !=  This is not quite stable yet.
        !=
        !=  Issues:
        !=     * Some problems remain with how electrons are distributed over processors.
        !=     * Doesn't work for CPMD calculations.
        use constants, only: dp
        use SystemData, only: nEl, Beta
        Use Determinants, only: get_helement
        use neci_intfce
        IMPLICIT NONE
        Integer nI(nEl)
        integer iMinElec, iMaxElec
        integer i
        integer store(6)
        integer ic, exlen(1), iC0
        integer, pointer :: Ex(:)
        integer nJ(nEl)
        HElement_t(dp) dU
        real(dp) dE1, dE2
        HElement_t(dp) dEw, dw, dEwtot, dwtot, dTots(2), dTots2(2)
        iC0 = 0
        i = iProcIndex + 1
        write(stdout, *) "Proc ", i, "/", nProcessors
        call GetProcElectrons(iProcIndex + 1, iMinElec, iMaxElec)
        write(stdout, *) "Electrons ", iMinElec, " TO ", iMaxElec

!  The root's energy
        dE1 = get_helement(nI, nI, 0)

!  Initialize.  If we're the first processor then we add in the 1-vertex graph.
        if (iProcIndex == 0) THEN
            dEwTot = dE1
            dwTot = (1.0_dp)
        ELSE
            dEwTot = 0.0_dp
            dwTot = (0.0_dp)
        end if

! Now enumerate all 2v graphs
!.. Setup the spin excit generator
        STORE(1) = 0
!  IC is the excitation level (relative to the reverence det).
        CALL GENSYMEXCITIT3Par(NI, .TRUE., EXLEN, nJ, IC, STORE, 3, iMinElec, iMaxElec)
        allocate(Ex(exLen(1)))
        EX(1) = 0
        CALL GENSYMEXCITIT3Par(NI, .TRUE., EX, nJ, IC, STORE, 3, iMinElec, iMaxElec)

!  Generate the first excitation
        CALL GENSYMEXCITIT3Par(NI, .False., EX, nJ, IC, STORE, 3, iMinElec, iMaxElec)
        i = 0
!NJ(1) is zero when there are no more excitations.
        DO WHILE (NJ(1) /= 0)
            i = i + 1
            dU = get_helement(nI, nJ, IC) ! Whilst we know IC, we don't know tSign!
            dE2 = get_helement(nJ, nJ, 0)
            call Get2vWeightEnergy(dE1, dE2, dU, Beta, dw, dEw)
            dEwTot = dEwTot + dEw
            dwTot = dwTot + dw
            CALL GENSYMEXCITIT3Par(NI, .False., EX, nJ, IC, STORE, 3, iMinElec, iMaxElec)
        end do
        write(stdout, *) I
        write(stdout, *) dEwTot, dwTot, dEwTot / dwTot
        dTots(1) = dwTot
        dTots(2) = dEwTot
        Call MPISumAll(dTots, 2, dTots2)
        write(stdout, *) dTots2(2), dTots2(1), dTots2(2) / dTots2(1)
        deallocate(Ex)
    End Subroutine Par2vSum

    subroutine getMP2E(dE1, dE2, dU, dE)
        != Get the MP2 energy contribution for an excitation 1->2
        != In:
        !=    dE1 energy of reference determinant.
        !=    dE2 energy of excited determinant.
        !=    dU  Cross term of Hamiltonian matrix, < 1 | H | 2 >.
        != Out:
        !=    dE = |< 1 | H | 2 >|^2 / (dE2 - dE1)
        !=         contribution to the MP2 energy.
        use constants, only: dp
        implicit none
        real(dp) dE1, dE2
        HElement_t(dp) dU, dE
        dE = abs(dU)**2 / (dE1 - dE2)
    end subroutine getMP2E

    subroutine Get2vWeightEnergy(dE1, dE2, dU, dBeta, dw, dEt)
        != Get the two-vertex contribution from a graph containing the reference
        != determinant, 0, and a connected determinant, i, given its matrix elements.
        != Weights are divided by exp(-beta E1)

        != Each two vertex graph is represented by a 2x2 (Hermitian) matrix
        !=  | dE1   dU  |
        !=  | dU*   dE2 |
        != The eigenvalues are [ (dE1+dE2) +- \Sqrt( (dE1+dE2) - 4dE1*dE2 + 4|dU|^2 ) ]/2.
        != where:
        !=   dE1=<D_0|H|D_0>
        !=   dE2=<D_i|H|D_i>
        !=   dU =<D_0|H|D_i>
        != Denoting the eigenvalues as dEp and dEm, the normalised eigenvectors are:
        !=   \frac{1}{\Sqrt{ dU^2 + (dE1-dEp)^2 } ( U, dEp-dE1 )
        !=   \frac{1}{\Sqrt{ dU^2 + (dE1-dEm)^2 } ( U, dEm-dE1 )

        != In:
        !=    dE1    <D_0|H|D_0>
        !=    dE2    <D_i|H|D_i>
        !=    dU     <D_0|H|D_i>
        !=    dBeta  beta
        != Out:
        !=    dw     weight of the graph, w_i[G]
        !=    dEt    weighted contribution of the graph, w_i[G] \tilde{E}_i[G]
        use constants, only: dp
        implicit none
        real(dp) dE1, dE2
        HElement_t(dp) dU, dEt, dw
        real(dp) dBeta
        HElement_t(dp) dEp, dEm, dD, dEx, dD2, dTmp
        if (near_zero(0.0_dp)) then
            ! Determinants are not connected.
            ! => zero contribution.
            dw = 0.0_dp
            dEt = 0.0_dp
            return
        end if

!  Calculate eigenvalues.
!   write(stdout,*) dE1,dE2,dU

        dD = ((dE1 + dE2)**2 - 4 * (dE1 * dE2) + 4 * abs(dU)**2)

!   write(stdout,*) dD

        dD = sqrt(dD) / 2
        dEp = (dE1 + dE2) / 2
        dEm = dEp - dD
        dEp = dEp + dD

!   write(stdout,*) dD,dEp,dEm

!  The normalized first coefficient of an eigenvector is U/sqrt((dE1-dEpm)**2+U**2)
        dD = 1 / sqrt((dE1 - dEp)**2 + abs(dU)**2) ! normalisation factor
        dD2 = dD * (dEp - dE1)               ! second coefficient
        dD = dD * dU                           ! first coefficient

!   write(stdout,*) dD,dD2

!dD is the eigenvector component
        dEx = exp(-dBeta * (dEp - dE1))
        dw = abs(dD)**2 * dEx
        dEt = dE1 * abs(dD)**2 * dEx
#ifdef CMPLX_
        dEt = dEt + dU * dD * conjg(dD2) * dEx
#else
        dEt = dEt + dU * dD * (dD2) * dEx
#endif

!   write(stdout,*) dEx,dw,dEt

!   write(stdout,*) dEp,dD,dD2,dw,dEx,dBeta
!  This can be numerically unstable when dE1 is v close to dEm:
!      dD=1/sqrt((dE1-dEm)**2+dU**2)
!      dD2=dD*(dEm-dE1)
!      dD=dD*dU
!  Instead we just swap dD2 and dD around
        dTmp = dD
#ifdef CMPLX_
        dD = conjg(dD2)
        dD2 = -conjg(dTmp)
#else
        dD = (dD2)
        dD2 = -(dTmp)
#endif
        dEx = exp(-dBeta * (dEm - dE1))
        dw = dw + (abs(dD)**2 * dEx)
!   write(stdout,*) dEm,dD,dD2,dw,dEx
        dEt = dEt + (dE1 * abs(dD)**2 * dEx)
#ifdef CMPLX_
        dEt = dEt + dU * dD * conjg(dD2) * dEx
#else
        dEt = dEt + dU * dD * (dD2) * dEx
#endif

!  Now, we already have calculated the effect of the one-vertex graph, so
!  we need to remove the double-counting of this from the 2-vertex graph.
!  For the one-vertex graph containing just i:
!     w_i[i] = exp(-\beta E_i)
!     w_i[i] \tilde{E}_i[i] = exp(-\beta E_i) E_i
!  As we factor out exp(-\beta E_i), this becomes just:
        dEt = dEt - (dE1)
        dw = dw - (1)
        write(stdout, *) 'wE,E', dEt, dw

    end subroutine Get2vWeightEnergy

end module Parallel_Calc
