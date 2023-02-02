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
      SUBROUTINE FindExcitDet(ExcitMat,nI,IC,TParity)
         use SystemData , only : NEl
         use lattice_mod, only: sort_unique
         use util_mod, only: stop_all
         implicit none
         integer, intent(in) :: ic
         integer, intent(inout) :: ExcitMat(2,ic)
         integer, intent(inout) :: nI(nel)
         logical, intent(out) :: tParity
#ifdef DEBUG_
         character(*), parameter :: this_routine = "FindExcitDet"
#endif
         integer :: i,j,k, pos_moved
         integer :: src(ic), sort_orbs(ic), sort_elecs(ic)

         ! adapt this for double excitations now..
#ifdef DEBUG_
        if (IC < 0 .or. ic > 3) then
            call stop_all(this_routine, "wrong ic!")
        end if
#endif

        if (ic == 3 .or. ic == 2) then
        !Ensure i<j and a<b
            ExcitMat(1,:) = sort_unique(ExcitMat(1,:))
            ExcitMat(2,:) = sort_unique(ExcitMat(2,:))
        end if

        sort_elecs = ExcitMat(1,:)

        ! Return Excitmat with orbitals, rather than indices being specified.
        ExcitMat(1,:) = nI(sort_elecs)

        ! this is the same functionality as in make_double.. but much more
        ! involved here.. so just take the core from make double and others..
        ! look there for the code!
        src = nI(sort_elecs)
        sort_orbs = ExcitMat(2,:)

        ! reshuffle the orbitals..
        if (ic > 1 ) then
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
                            nI(i-1) = sort_orbs(k)
                            exit
                        else
                            nI(i-1) = nI(i)
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
                    do i = sort_elecs(k)-1, 1, -1
                        if (sort_orbs(k) > nI(i)) then
                            nI(i+1) = sort_orbs(k)
                            exit
                        else
                            nI(i+1) = nI(i)
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
         RETURN

      END SUBROUTINE FindExcitDet


      LOGICAL FUNCTION ISVALIDDET(DET,NEL)
         IMPLICIT NONE
         INTEGER, intent(in) :: NEL
         integer, intent(in) :: DET(NEL)
         INTEGER I
         IF(DET(1).LT.1) THEN
            ISVALIDDET=.FALSE.
            RETURN
         ENDIF
         DO I=2,NEL
            IF(DET(I-1).GE.DET(I)) THEN
               ISVALIDDET=.FALSE.
               RETURN
            ENDIF
         ENDDO
         ISVALIDDET=.TRUE.
         RETURN
      END

C.. Get the orbitals which are excited in going from I to J
C.. EX(1,*) are in I, and EX(2,*) are in J
C.. TSIGN is set to the parity of the permutations required to line up the orbitals
C..   TRUE means ODD.
C.. If there are too many excitations to fit, then we put -excitlevel in EX(1,1) and EX(2,1)
!  EX(1,1) is the max number of excitations (passed in as a parameter)
      SUBROUTINE GETEXCITATION(NI,NJ,NEL,EX,TSIGN)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NJ(NEL),EX(2,*)
         INTEGER I,J,IPAR
         INTEGER IC1,IC2
         INTEGER iMaxExcit
         LOGICAL TSIGN
         iMaxExcit=EX(1,1)
         EX(1:2,1:iMaxExcit)=0
         IC1=0
         IC2=0
         I=1
         J=1
         IPAR=0
C         CALL WRITEDET(6,NI,NEL,.TRUE.)
C         CALL WRITEDET(6,NJ,NEL,.TRUE.)
         DO WHILE(I.LE.NEL.AND.J.LE.NEL)
C.. Differences from I to J
C            WRITE(6,*) "GE",I,J
            DO WHILE(I.LE.NEL)
               if (NI(I) >= NJ(J)) exit
               IC1=IC1+1
               IF(IC1.LE.iMaxExcit) THEN
                  EX(1,IC1)=NI(I)
                  IPAR=IPAR+I
               ENDIF
               I=I+1
            ENDDO
C.. Differences from J to I
            DO WHILE(I.LE.NEL.AND.J.LE.NEL)
               if (NI(I) <= NJ(J)) exit
               IC2=IC2+1
               IF(IC2.LE.iMaxExcit) THEN
                  EX(2,IC2)=NJ(J)
                  IPAR=IPAR+J
               ENDIF
               J=J+1
            ENDDO
            IF(I.LE.NEL.AND.J.LE.NEL) then
               if (NI(I) == NJ(J)) then
                   I=I+1
                   J=J+1
               endif
            ENDIF
         ENDDO
C.. Deal with remaining I
         DO WHILE(I.LE.NEL)
            IC1=IC1+1
            IF(IC1.LE.iMaxExcit) THEN
               IPAR=IPAR+I
               EX(1,IC1)=NI(I)
            ENDIF
            I=I+1
         ENDDO
C.. Deal with remaining J
         DO WHILE(J.LE.NEL)
            IC2=IC2+1
            IF(IC2.LE.iMaxExcit) THEN
               EX(2,IC2)=NJ(J)
               IPAR=IPAR+J
            ENDIF
            J=J+1
         ENDDO
         IF(iC1.GT.iMaxExcit) THEN
C.. we actually needed more space.  Just list the excitation counts (-ve)
            DO i=1,iMaxExcit
               IF (i.EQ.1) THEN
                  EX(1,1)=-iC1
                  EX(2,1)=-iC2
               ELSE
                  EX(1,i)=0
                  EX(2,i)=0
               ENDIF
            ENDDO
         ELSEIF(iC1.EQ.0) THEN
            EX(1,1)=0
            EX(2,1)=0
         ENDIF
         TSIGN=BTEST(IPAR,0)
         RETURN
      END


C.. Generate excitations of order at most NORDER from NI (excluding NI
Citself) in order
C.. ICLIST is the list of orders of the excitations
C.. NMIN is the minimum level of excitation
C.. IF NORDER=0, we SET NORDER=NEL

C.. THis does not use the excitation generators, but enumerates all possible
C.. determinants within double excitations, which will end up fairly inefficient.
        SUBROUTINE GENEXCIT(NI,NORDER,NBASIS,NEL,LIST,ICLIST,NLIST,NMIN,
     &            G1,TSYM,NBASISMAX,TCOUNT)
            use CalcData , only : TNEWEXCITATIONS
            use SystemData, only: BasisFN
            use sym_mod, only: getsym
            use constants, only: int64
            use util_mod, only: stop_all
            IMPLICIT NONE
            INTEGER NEL,NI(NEL),NORDER,NBASIS,LIST(NEL,*),ICLIST(*)
            INTEGER NEXCIT(NEL),NLIST,nBasisMax(5,*)
            TYPE(BASISFN) ISYM
            INTEGER NO,IFDET
            LOGICAL TSYM,TCOUNT
            INTEGER NMIN
            type(BasisFN) G1(*)
            INTEGER nBRR(1)
            integer(int64) STORE(6)
            INTEGER ICOUNT,ILEVEL
            character(*), parameter :: t_r="GENEXCIT"


            NO=NORDER
            IF(NO.EQ.0) NO=NEL
            IF(TSYM) CALL GETSYM(NI,NEL,G1,NBASISMAX,ISYM)
            IF(TNEWEXCITATIONS.AND.NORDER.LE.2) THEN
               IF(NORDER.EQ.1) ILEVEL=1
               IF(NORDER.EQ.2) THEN
                  IF(NMIN.EQ.1) THEN
                     ILEVEL=3
                  ELSE
                     ILEVEL=2
                  ENDIF
               ENDIF
C.. use a new excitation routine - THIS DOES NOT PRODUCE AN ORDERED LIST
               CALL SYMSETUPEXCITS(NI,NEL,NBASIS,STORE,
     &            TCOUNT,ICOUNT,ILEVEL,1,nEl)
               IF(.NOT.TCOUNT) THEN
                  IF(ICOUNT.GT.NLIST) THEN
                     WRITE(6,*) "Exctiation list size",NLIST
                     WRITE(6,*) "Requires size",ICOUNT
                     call stop_all(t_r,
     &                          "EXCITATIONS LIST NOT LARGE ENOUGH")
                  ENDIF
                  call stop_all(t_r,"Call to symgenexcits removed due
     &                  to incorrect arguments passes")
!                  CALL SYMGENEXCITS(NI,NEL,STORE,
!     &               LIST,ICLIST,ICOUNT)
               ENDIF
               NLIST=ICOUNT
            ELSE
               NLIST=1
               CALL GENEXCIT_R(NI,NO,NEL,1,NBASIS,LIST,ICLIST,
     &            NLIST,NEXCIT,NO,1,NMIN,G1,ISYM,TSYM,NBASISMAX,TCOUNT)
               NLIST=NLIST-1
            ENDIF
            RETURN
        END

C.. AT 29/1/04 - Looks like this won't scale very well to lots of
Celectrons
C..      as there will be very few which are within NORDER excitations
Cof
C..      our original det. (Or for that matter a large basis)
C.. Recursively go through each electron ( we're on NELEC), filling
CLIST.
C.. with up to NORDER excitations.
C.. NEXCIT is the determinant being constructed.
C.. NLEFT is the number of possible excitations left
C.. NSTARTFN is the basis fn to start with for this electron
      RECURSIVE SUBROUTINE GENEXCIT_R(NI,NORDER,NEL,NELEC,NBASIS,LIST,
     &         ICLIST,NLISTPOS,NEXCIT,NLEFT,NSTARTFN,NMIN,G1,ISYM,TSYM,
     &         NBASISMAX,TCOUNT)
         use SystemData, only: BasisFN,tFixLz
         use sym_mod, only: lChkSym, GetLz, getsym
         use util_mod, only: NECI_ICOPY
         IMPLICIT NONE
         INTEGER NEL,NELEC,NORDER,NLISTPOS,NLEFT,NBASIS,NSTARTFN,NNLEFT
         INTEGER NI(NEL),NEXCIT(NEL),LIST(NEL,*),ICLIST(*)
         INTEGER I,J,NMIN,NMAXEX,NExcitMl
         LOGICAL LISINOLD,LSYM,TSYM,TCOUNT,LZSYM
         TYPE(BASISFN) G1(NBASIS),ISYM,ISYM2
         INTEGER nBasisMax(5,*)
C.. NMAXEX is the maximum number of excitations we're allowed to  have
Cleft
         NMAXEX=NORDER-NMIN
C..  I is the new basis fn for electron NELEC
C.. We see if it's in the original determinant
C.. we need to leave at least NEL-NELEC basis fns for the following
Celectrons
         DO I=NSTARTFN,NBASIS-(NEL-NELEC)
            LISINOLD=.FALSE.
            DO J=1,NEL
               IF(NI(J).EQ.I) LISINOLD=.TRUE.
            ENDDO
            NEXCIT(NELEC)=I
            NNLEFT=NLEFT
            IF(.NOT.LISINOLD) NNLEFT=NNLEFT-1
C.. NNLEFT is the number of excitations left
C.. check to see if we're actually exciting an electron, and if we've
Cany
C.. spare excitations left
C.. if we've allocated the last electron, we need to see if the det
C.. is allowed, and if so, store it
            IF(NELEC.EQ.NEL) THEN
C.. if we've excited at least one elec
               IF(NNLEFT.LE.NMAXEX.AND.NNLEFT.GE.0) THEN
C.. we check whether we're allowed this excitation owing to sym
C                  CALL WRITEDET(6,NEXCIT,NEL,.TRUE.)
                  IF(TSYM) THEN
                     CALL GETSYM(NEXCIT,NEL,G1,NBASISMAX,ISYM2)
                     LSYM=LCHKSYM(ISYM,ISYM2)
                  ELSE
                     LSYM=.TRUE.
                  ENDIF

!We also want to check Lz symmetry
                  IF(tFixLz) THEN
                      CALL GetLz(NEXCIT,NEL,NExcitMl)
                      IF(NExcitMl.eq.(ISYM%Ml)) THEN
                          LZSYM=.true.
                      ELSE
                          LZSYM=.false.
                      ENDIF
                  ELSE
                      LZSYM=.true.
                  ENDIF

!                  WRITE(6,*) "ISym:",ISYM
!                  WRITE(6,*) "ISym2",ISYM2
!                  WRITE(6,*) "LSYM",LSYM
                  IF(LSYM.and.LZSYM) THEN
                     IF(.NOT.TCOUNT) THEN
                        CALL NECI_ICOPY(NEL,NEXCIT,1,
     &                                           LIST(1:NEL,NLISTPOS),1)
                        ICLIST(NLISTPOS)=NORDER-NNLEFT
                     ENDIF
                     NLISTPOS=NLISTPOS+1
                  ENDIF
               ENDIF
            ELSEIF(NNLEFT.GE.0) THEN
C.. we start the next electron on the basis fn after this one
               CALL GENEXCIT_R(NI,NORDER,NEL,NELEC+1,NBASIS,LIST,ICLIST,
     &            NLISTPOS,NEXCIT,NNLEFT,I+1,NMIN,G1,ISYM,TSYM,
     &            NBASISMAX,TCOUNT)
            ENDIF
         ENDDO
         RETURN
      END


