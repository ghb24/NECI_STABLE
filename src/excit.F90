#include "macros.h"

!This is an O[2*N] operation for doubles or O[N] operation for single excitations.
!It takes a determinant nI, and a 2x2 matrix indicating the excited 'to' and 'from'
!orbitals, and finds the determinant nJ, while keeping the final determinant ordered.
!This will save having to order the final determinant.
!We also calculate the parity of the excitation. This is needed for the scrules as
!i.e. <1 2 3 4 | H | 1 3 5 6> is - <2 4 || 5 6> but <1 2 3 4 | H | 1 2 5 6 > = +<3 4 || 5 6>

!1 2 3 4 -> 1 2 5 6 , it would give 3->5 and 4->6 with a +ve parity
!because no interchanges were required to line up the dets
!1 2 3 4 -> 1 3 5 6 we would swap 2 and 3 around to give a negative parity and have  2->5 and 4->6
!so I just need to work out how much the common orbitals between the determinants have moved between
!the ordered list of orbitals. If its an even number, parity is positive (false), otherwise negative (true).

!ExcitMat(1,*) are the **indices** in the determinant to vacate from nI (the i,j pair)
!ExcitMat(2,*) are the orbitals to occupy in nJ (the a,b pair) (not the index, but the actual orbital)
!**However, when the excitmat is returned, all elements refer to ORBITALS, NOT INDICES.**
!IC should be 1 or 2, depending on whether it is a double or single excitation
!Single excitations should just have ExcitMat(1,1) and ExcitMat(2,1) with orbital
!information.
!The algorithm could be improved for double excitations by only searching through the
!determinant once, reducing it from an O[2N] to O[N] operation, though would be a little
!more fiddly...
! [W.D. 11.12.2017]:
! never thought to touch this part of the code ever.. but we need triple excitations
! now too.. so adapt this functionality here.., since this function is
! unfortunately called in too many parts of the code..
module excit_mod
    use CalcData, only: TNEWEXCITATIONS
    use SystemData, only: BasisFN, tFixLz, nEl
    use constants, only: int64
    use lattice_mod, only: sort_unique
    use sym_mod, only: lChkSym, GetLz, getsym
    use util_mod, only: NECI_ICOPY, stop_all
    better_implicit_none
    private
    public :: FindExcitDet, GETEXCITATION, GENEXCIT, isvaliddet

contains
    SUBROUTINE FindExcitDet(ExcitMat, nI, IC, TParity)
        integer, intent(in) :: ic
        integer, intent(inout) :: ExcitMat(2, ic)
        integer, intent(inout) :: nI(nel)
        logical, intent(out) :: tParity
#ifdef DEBUG_
        character(*), parameter :: this_routine = "FindExcitDet"
#endif
        integer :: i, k, pos_moved
        integer :: src(ic), sort_orbs(ic), sort_elecs(ic)

    ! adapt this for double excitations now..
#ifdef DEBUG_
        if (IC < 0 .or. ic > 3) then
            call stop_all(this_routine, "wrong ic!")
        end if
#endif

        if (ic == 3 .or. ic == 2) then
    !Ensure i<j and a<b
            ExcitMat(1, :) = sort_unique(ExcitMat(1, :))
            ExcitMat(2, :) = sort_unique(ExcitMat(2, :))
        end if

        sort_elecs = ExcitMat(1, :)

    ! Return Excitmat with orbitals, rather than indices being specified.
        ExcitMat(1, :) = nI(sort_elecs)

    ! this is the same functionality as in make_double.. but much more
    ! involved here.. so just take the core from make double and others..
    ! look there for the code!
        src = nI(sort_elecs)
        sort_orbs = ExcitMat(2, :)

    ! reshuffle the orbitals..
        if (ic > 1) then
    ! or check it individually:
            if (src(2) < sort_orbs(1)) then
    ! then i hops over j:
                sort_elecs(2) = sort_elecs(2) - 1
            end if
            if (ic == 3) then
            if (src(3) < sort_orbs(1)) then
    ! then i hops over k
    ! (note: this also implies that j hops over k, but treat that
    ! seperately below, to also cover the case, where this if here
    ! is not fullfilled!)
                sort_elecs(3) = sort_elecs(3) - 1
            end if
            if (src(3) < sort_orbs(2)) then
    ! then j hops over k
                sort_elecs(3) = sort_elecs(3) - 1
            end if
            end if
        end if

        pos_moved = 0 !This keeps track of how far the common orbitals have moved.

        do k = 1, ic
        if (src(k) < sort_orbs(k)) then
        if (sort_elecs(k) == nel) then
    ! this can only happen for k == 3
            i = nel + 1
            nI(nel) = sort_orbs(k)
        else
            do i = sort_elecs(k) + 1, nel
            if (sort_orbs(k) < nI(i)) then
                nI(i - 1) = sort_orbs(k)
                exit
            else
                nI(i - 1) = nI(i)
            end if
            end do
            if (i == nel + 1) then
                nI(nel) = sort_orbs(k)
            end if
        end if
        else
        if (sort_elecs(k) == 1) then
            i = 0
            nI(1) = sort_orbs(k)
        else
            do i = sort_elecs(k) - 1, 1, -1
            if (sort_orbs(k) > nI(i)) then
                nI(i + 1) = sort_orbs(k)
                exit
            else
                nI(i + 1) = nI(i)
            end if
            end do
            if (i == 0) then
                nI(1) = sort_orbs(k)
            end if
        end if
        end if

        pos_moved = pos_moved + sort_elecs(k) - i + 1

        end do

        tParity = btest(pos_moved, 0)
    END SUBROUTINE FindExcitDet

    logical pure function isvaliddet(det, nel)
        integer, intent(in) :: nel
        integer, intent(in) :: det(nel)
        integer i
        if (det(1) < 1) then
            isvaliddet = .false.
            return
        end if
        do i = 2, nel
        if (det(i - 1) >= det(i)) then
            isvaliddet = .false.
            return
        end if
        end do
        isvaliddet = .true.
    end

    !.. Get the orbitals which are excited in going from I to J
    !.. EX(1,*) are in I, and EX(2,*) are in J
    !.. TSIGN is set to the parity of the permutations required to line up the orbitals
    !..   TRUE means ODD.
    !.. If there are too many excitations to fit, then we put -excitlevel in EX(1,1) and EX(2,1)
    !  EX(1,1) is the max number of excitations (passed in as a parameter)
    SUBROUTINE GETEXCITATION(NI, NJ, NEL, EX, TSIGN)
        INTEGER NEL, NI(NEL), NJ(NEL), EX(2, *)
        INTEGER I, J, IPAR
        INTEGER IC1, IC2
        INTEGER iMaxExcit
        LOGICAL TSIGN
        iMaxExcit = EX(1, 1)
        EX(1:2, 1:iMaxExcit) = 0
        IC1 = 0
        IC2 = 0
        I = 1
        J = 1
        IPAR = 0
    !         CALL WRITEDET(6,NI,NEL,.TRUE.)
    !         CALL WRITEDET(6,NJ,NEL,.TRUE.)
        DO WHILE (I <= NEL .AND. J <= NEL)
    !.. Differences from I to J
    !            WRITE(stdout,*) "GE",I,J
            DO WHILE (I <= NEL)
                if (NI(I) >= NJ(J)) exit
                IC1 = IC1 + 1
                IF (IC1 <= iMaxExcit) THEN
                    EX(1, IC1) = NI(I)
                    IPAR = IPAR + I
                END IF
                I = I + 1
            END DO
    !.. Differences from J to I
            DO WHILE (I <= NEL .AND. J <= NEL)
                if (NI(I) <= NJ(J)) exit
                IC2 = IC2 + 1
                IF (IC2 <= iMaxExcit) THEN
                    EX(2, IC2) = NJ(J)
                    IPAR = IPAR + J
                END IF
                J = J + 1
            END DO
            IF (I <= NEL .AND. J <= NEL) then
            if (NI(I) == NJ(J)) then
                I = I + 1
                J = J + 1
            end if
            END IF
        END DO
    !.. Deal with remaining I
        DO WHILE (I <= NEL)
            IC1 = IC1 + 1
            IF (IC1 <= iMaxExcit) THEN
                IPAR = IPAR + I
                EX(1, IC1) = NI(I)
            END IF
            I = I + 1
        END DO
    !.. Deal with remaining J
        DO WHILE (J <= NEL)
            IC2 = IC2 + 1
            IF (IC2 <= iMaxExcit) THEN
                EX(2, IC2) = NJ(J)
                IPAR = IPAR + J
            END IF
            J = J + 1
        END DO
        IF (iC1 > iMaxExcit) THEN
    !.. we actually needed more space.  Just list the excitation counts (-ve)
            DO i = 1, iMaxExcit
            IF (i == 1) THEN
                EX(1, 1) = -iC1
                EX(2, 1) = -iC2
            ELSE
                EX(1, i) = 0
                EX(2, i) = 0
            END IF
            END DO
        ELSEIF (iC1 == 0) THEN
            EX(1, 1) = 0
            EX(2, 1) = 0
        END IF
        TSIGN = BTEST(IPAR, 0)
    END

    ! !.. Generate excitations of order at most NORDER from NI (excluding NI
    ! !itself) in order
    ! !.. ICLIST is the list of orders of the excitations
    ! !.. NMIN is the minimum level of excitation
    ! !.. IF NORDER=0, we SET NORDER=NEL
    !
    !.. THis does not use the excitation generators, but enumerates all possible
    !.. determinants within double excitations, which will end up fairly inefficient.
    SUBROUTINE GENEXCIT(NI, NORDER, NBASIS, NEL, LIST, ICLIST, NLIST, NMIN, G1, TSYM, NBASISMAX, TCOUNT)
        integer, intent(in) :: NEL, NI(NEL), NORDER, NBASIS, LIST(NEL, *), ICLIST(*)
        INTEGER NEXCIT(NEL), NLIST, nBasisMax(5, *)
        TYPE(BASISFN) ISYM
        INTEGER NO
        LOGICAL TSYM, TCOUNT
        INTEGER NMIN
        type(BasisFN) G1(*)
        integer(int64) STORE(6)
        INTEGER ICOUNT, ILEVEL
        character(*), parameter :: t_r = "GENEXCIT"

        NO = NORDER
        IF (NO == 0) NO = NEL
        IF (TSYM) CALL GETSYM(NI, NEL, G1, NBASISMAX, ISYM)
        IF (TNEWEXCITATIONS .AND. NORDER <= 2) THEN
            ! Here the non-working routine `SYMSETUPEXCITS` was called.
            call stop_all(t_r, "Should never be here")
        ELSE
            NLIST = 1
            CALL GENEXCIT_R(NI, NO, NEL, 1, NBASIS, LIST, ICLIST, NLIST, NEXCIT, NO, 1, NMIN, G1, ISYM, TSYM, NBASISMAX, TCOUNT)
            NLIST = NLIST - 1
        END IF
    END

    !.. AT 29/1/04 - Looks like this won't scale very well to lots of
    !electrons
    !..      as there will be very few which are within NORDER excitations
    !of
    !..      our original det. (Or for that matter a large basis)
    !.. Recursively go through each electron ( we're on NELEC), filling
    !LIST.
    !.. with up to NORDER excitations.
    !.. NEXCIT is the determinant being constructed.
    !.. NLEFT is the number of possible excitations left
    !.. NSTARTFN is the basis fn to start with for this electron
    RECURSIVE SUBROUTINE GENEXCIT_R(NI, NORDER, NEL, NELEC, NBASIS, LIST, ICLIST, &
                                    NLISTPOS, NEXCIT, NLEFT, NSTARTFN, NMIN, G1, ISYM, TSYM, NBASISMAX, TCOUNT)
        INTEGER NEL, NELEC, NORDER, NLISTPOS, NLEFT, NBASIS, NSTARTFN, NNLEFT
        INTEGER NI(NEL), NEXCIT(NEL), LIST(NEL, *), ICLIST(*)
        INTEGER I, J, NMIN, NMAXEX, NExcitMl
        LOGICAL LISINOLD, LSYM, TSYM, TCOUNT, LZSYM
        TYPE(BASISFN) G1(NBASIS), ISYM, ISYM2
        INTEGER nBasisMax(5, *)
    !.. NMAXEX is the maximum number of excitations we're allowed to  have
    !left
        NMAXEX = NORDER - NMIN
    !..  I is the new basis fn for electron NELEC
    !.. We see if it's in the original determinant
    !.. we need to leave at least NEL-NELEC basis fns for the following
    !electrons
        DO I = NSTARTFN, NBASIS - (NEL - NELEC)
            LISINOLD = .FALSE.
            DO J = 1, NEL
                IF (NI(J) == I) LISINOLD = .TRUE.
            END DO
            NEXCIT(NELEC) = I
            NNLEFT = NLEFT
            IF (.NOT. LISINOLD) NNLEFT = NNLEFT - 1
    !.. NNLEFT is the number of excitations left
    !.. check to see if we're actually exciting an electron, and if we've
    !any
    !.. spare excitations left
    !.. if we've allocated the last electron, we need to see if the det
    !.. is allowed, and if so, store it
            IF (NELEC == NEL) THEN
    !.. if we've excited at least one elec
                IF (NNLEFT <= NMAXEX .AND. NNLEFT >= 0) THEN
    !.. we check whether we're allowed this excitation owing to sym
    !                  CALL WRITEDET(6,NEXCIT,NEL,.TRUE.)
                    IF (TSYM) THEN
                        CALL GETSYM(NEXCIT, NEL, G1, NBASISMAX, ISYM2)
                        LSYM = LCHKSYM(ISYM, ISYM2)
                    ELSE
                        LSYM = .TRUE.
                    END IF

    !We also want to check Lz symmetry
                    IF (tFixLz) THEN
                        CALL GetLz(NEXCIT, NEL, NExcitMl)
                        IF (NExcitMl == (ISYM%Ml)) THEN
                            LZSYM = .true.
                        ELSE
                            LZSYM = .false.
                        END IF
                    ELSE
                        LZSYM = .true.
                    END IF

    !                  WRITE(stdout,*) "ISym:",ISYM
    !                  WRITE(stdout,*) "ISym2",ISYM2
    !                  WRITE(stdout,*) "LSYM",LSYM
                    IF (LSYM .and. LZSYM) THEN
                    IF (.NOT. TCOUNT) THEN
                        CALL NECI_ICOPY(NEL, NEXCIT, 1, LIST(1:NEL, NLISTPOS), 1)
                        ICLIST(NLISTPOS) = NORDER - NNLEFT
                    END IF
                    NLISTPOS = NLISTPOS + 1
                    END IF
                END IF
            ELSEIF (NNLEFT >= 0) THEN
    !.. we start the next electron on the basis fn after this one
                CALL GENEXCIT_R(NI, NORDER, NEL, NELEC + 1, NBASIS, LIST, ICLIST, &
                                NLISTPOS, NEXCIT, NNLEFT, I + 1, NMIN, G1, ISYM, TSYM, NBASISMAX, TCOUNT)
            END IF
        END DO
    END
end module

