      LOGICAL FUNCTION AREDETSEXCITS(D1,D2,NEL,NBASISMAX,ICE)
         IMPLICIT NONE
         INTEGER ICE,NEL,D1(NEL),D2(NEL),NBASISMAX(*)
         LOGICAL ISUHFDET
         INTEGER IGETEXCITLEVEL
         ICE=IGETEXCITLEVEL(D1,D2,NEL)
         IF(ICE.LE.2) THEN
C.. The new node is connected to this node
            IF(ISUHFDET(D1,NEL,NBASISMAX).OR.
     &         ISUHFDET(D2,NEL,NBASISMAX)) THEN
               IF(ICE.NE.1) THEN
                  AREDETSEXCITS=.TRUE.
               ELSE
                  AREDETSEXCITS=.FALSE.
               ENDIF
            ELSE
               AREDETSEXCITS=.TRUE.
            ENDIF
         ELSE
            AREDETSEXCITS=.FALSE.
         ENDIF
         RETURN
      END

      LOGICAL FUNCTION ISVALIDDET(DET,NEL)
         IMPLICIT NONE
         INTEGER NEL,DET(NEL)
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
C.. IEX(1,*) are in I, and IEX(2,*) are in J
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
         CALL IAZZERO(EX,2*iMaxExcit)
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
            DO WHILE(I.LE.NEL.AND.NI(I).LT.NJ(J))
               IC1=IC1+1
               IF(IC1.LE.iMaxExcit) THEN
                  EX(1,IC1)=NI(I)
                  IPAR=IPAR+I
               ENDIF
               I=I+1
            ENDDO
C.. Differences from J to I
            DO WHILE(I.LE.NEL.AND.J.LE.NEL.AND.NI(I).GT.NJ(J))
               IC2=IC2+1
               IF(IC2.LE.iMaxExcit) THEN
                  EX(2,IC2)=NJ(J)
                  IPAR=IPAR+J
               ENDIF
               J=J+1
            ENDDO
            IF(I.LE.NEL.AND.J.LE.NEL.AND.NI(I).EQ.NJ(J)) THEN
               I=I+1
               J=J+1
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
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NORDER,NBASIS,LIST(NEL,*),ICLIST(*)
         INTEGER NEXCIT(NEL),NLIST,NBASISMAX(5,7)
         INCLUDE 'sym.inc'
         TYPE(BASISFN) ISYM
         INTEGER NO,IFDET
         LOGICAL TSYM,TCOUNT
         INTEGER G1(*),NMIN

         INTEGER STORE(6),ICOUNT,ILEVEL


         NO=NORDER
         IF(NO.EQ.0) NO=NEL
         IF(TSYM) CALL GETSYM(NI,NEL,G1,NBASISMAX,ISYM)
         IF(NBASISMAX(4,7).EQ.1) THEN
            NLIST=1
C.. If we're actually using CSFs
            CALL GNCSFS(NEL,NBASIS,0,NBASISMAX,
     &         LIST,TCOUNT,G1,TSYM,ISYM%MS,TSYM,ISYM,NLIST,IFDET,
     &         .TRUE.,NI,ICLIST,.TRUE.,NO)
         ELSE
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
               CALL SYMSETUPEXCITS(NI,NEL,G1,NBASIS,NBASISMAX,STORE,
     &            TCOUNT,ICOUNT,ILEVEL)
               IF(.NOT.TCOUNT) THEN
                  IF(ICOUNT.GT.NLIST) THEN
                     WRITE(6,*) "Exctiation list size",NLIST
                     WRITE(6,*) "Requires size",ICOUNT
                     STOP "EXCITATIONS LIST NOT LARGE ENOUGH"
                  ENDIF
                  CALL SYMGENEXCITS(NI,NEL,G1,NBASIS,NBASISMAX,STORE,
     &               LIST,ICLIST,ICOUNT)
               ENDIF
               NLIST=ICOUNT
            ELSE
               NLIST=1
               CALL GENEXCIT_R(NI,NO,NEL,1,NBASIS,LIST,ICLIST,
     &            NLIST,NEXCIT,NO,1,NMIN,G1,ISYM,TSYM,NBASISMAX,TCOUNT)
               NLIST=NLIST-1
            ENDIF
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
         IMPLICIT NONE
         INTEGER NEL,NELEC,NORDER,NLISTPOS,NLEFT,NBASIS,NSTARTFN,NNLEFT
         INTEGER NI(NEL),NEXCIT(NEL),LIST(NEL,*),ICLIST(*)
         INTEGER I,J,I1,I2,NMIN,NMAXEX,K
         LOGICAL LISINOLD,LSYM,TSYM,LCHKSYM,TCOUNT
         INCLUDE 'basis.inc'
         TYPE(BASISFN) G1(NBASIS),ISYM,ISYM2
         INTEGER NBASISMAX(5,3)
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
C                  WRITE(6,*) ISYM,ISYM2,LSYM
                  IF(LSYM) THEN
                     IF(.NOT.TCOUNT) THEN
                        CALL ICOPY(NEL,NEXCIT,1,LIST(1:NEL,NLISTPOS),1)
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


