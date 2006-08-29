      INTEGER FUNCTION ICMPCSFS(C1,C2,NEL)
         IMPLICIT NONE
         INTEGER NEL,I
         INTEGER C1(1:NEL),C2(1:NEL)
         INTEGER J,K,L
         INCLUDE 'csf.inc'
         L=0
         DO I=1,NEL
            J=C1(I)
C.. Remove the last two bits
            IF(J.GE.CSF_NBSTART) THEN
               IF(L.EQ.0) L=I
              J=CSF_NBSTART+IAND(J-CSF_NBSTART,65532)
            ENDIF
            K=C2(I)
            IF(K.GE.CSF_NBSTART)
     &         K=CSF_NBSTART+IAND(K-CSF_NBSTART,65532)
            IF(J.LT.K) THEN
               ICMPCSFS=-1
               RETURN
            ELSE IF(J.GT.K) THEN
               ICMPCSFS=1
               RETURN
            ENDIF
         ENDDO
C.. Same spatial, so compare spins
         DO I=L,NEL
            J=C1(I)
C.. Use the last two bits
            IF(J.GE.CSF_NBSTART) THEN
              J=IAND(J-CSF_NBSTART,3)
            ENDIF
            K=C2(I)
            IF(K.GE.CSF_NBSTART)
     &         K=IAND(K-CSF_NBSTART,3)
            IF(J.LT.K) THEN
               ICMPCSFS=-1
               RETURN
            ELSE IF(J.GT.K) THEN
               ICMPCSFS=1
               RETURN
            ENDIF
         ENDDO
C.. Same
         ICMPCSFS=0
         RETURN
      END

      SUBROUTINE GETCSFFROMDET(NI,NCJ,NEL,STOT,SZ)
         IMPLICIT NONE
         INTEGER NEL,NI(NEL),NCJ(NEL),STOT,SZ,IOPEN(NEL)
         INTEGER I,J,K
         INCLUDE 'csf.inc'
         J=0
         I=1
         K=0
         DO WHILE (I.LE.NEL)
            IF(I.LT.NEL.AND.
     &         (IAND(NI(I)-1,65534).EQ.IAND(NI(I+1)-1,65534))) THEN
C.. we've found a closed orbital - use it.
               J=J+1
               NCJ(J)=NI(I)
               I=I+1
               J=J+1
               NCJ(J)=NI(I)
               I=I+1
            ELSE
               K=K+1
C.. open shell - store elsewhere
               IOPEN(K)=IAND(NI(I)-1,65534)
               I=I+1
            ENDIF
         ENDDO
C.. Now deal with the open shell
         I=0
         DO WHILE(I.LT.K)
            J=J+1
            I=I+1
            NCJ(J)=CSF_NBSTART
            NCJ(J)=NCJ(J)+IOPEN(I)*4
            IF(I.GT.(K-SZ)/2) THEN
C.. betas first then alphas
               NCJ(J)=NCJ(J)+1
            ENDIF
            IF(I.LE.(K+STOT)/2) THEN
C.. pluses first then minuses
               NCJ(J)=NCJ(J)+2
            ENDIF
         ENDDO 
         RETURN
      END
      
C.. Get an (S, Sz) encoded LMS from a CSF
      SUBROUTINE CSFGETSPIN(NCI,NEL,LMS) 
         IMPLICIT NONE
         INTEGER NEL,NCI(NEL)
         INTEGER LMS
         INCLUDE 'csf.inc'

         INTEGER ICLOSEDI(NEL),NCLOSEDI,IOPENI(NEL),NOPENI,SZI,SI,
     &      SORDERI(NEL)
         
         CALL CSFGETCSFDATA(NCI,NEL,ICLOSEDI,NCLOSEDI,IOPENI,NOPENI,SZI,
     &      SI,SORDERI)
         LMS=SI*CSF_NBSTART+SZI
         RETURN
      END        
C.. Return true if two csfs have the same spin sym.
      LOGICAL FUNCTION CSFCHKSPINSYM(NCI,NCJ,NEL)
         IMPLICIT NONE
         INTEGER NEL,NCI(NEL),NCJ(NEL)

         INTEGER ICLOSEDI(NEL),NCLOSEDI,IOPENI(NEL),NOPENI,SZI,SI,
     &      SORDERI(NEL)
         INTEGER ICLOSEDJ(NEL),NCLOSEDJ,IOPENJ(NEL),NOPENJ,SZJ,SJ,
     &      SORDERJ(NEL)
         INTEGER I,J,K
         INTEGER ICMPDETS
         
C         IEXLEVEL=1
C         RETURN
         CALL CSFGETCSFDATA(NCI,NEL,ICLOSEDI,NCLOSEDI,IOPENI,NOPENI,SZI,
     &      SI,SORDERI)
         CALL CSFGETCSFDATA(NCJ,NEL,ICLOSEDJ,NCLOSEDJ,IOPENJ,NOPENJ,SZJ,
     &      SJ,SORDERJ)
         CSFCHKSPINSYM=.TRUE.
         IF(SZI.NE.SZJ) CSFCHKSPINSYM=.FALSE.
         IF(SI.NE.SJ) CSFCHKSPINSYM=.FALSE.
         RETURN
      END

      SUBROUTINE GETCSFEXCITLEVEL(NCI,NCJ,NEL,IEXLEVEL)
         IMPLICIT NONE
         INTEGER NEL,NCI(NEL),NCJ(NEL),IEXLEVEL

         INTEGER ICLOSEDI(NEL),NCLOSEDI,IOPENI(NEL),NOPENI,SZI,SI,
     &      SORDERI(NEL)
         INTEGER ICLOSEDJ(NEL),NCLOSEDJ,IOPENJ(NEL),NOPENJ,SZJ,SJ,
     &      SORDERJ(NEL)
         INTEGER I,J,K
         INTEGER ICMPDETS
         
C         IEXLEVEL=1
C         RETURN
         CALL CSFGETCSFDATA(NCI,NEL,ICLOSEDI,NCLOSEDI,IOPENI,NOPENI,SZI,
     &      SI,SORDERI)
         CALL CSFGETCSFDATA(NCJ,NEL,ICLOSEDJ,NCLOSEDJ,IOPENJ,NOPENJ,SZJ,
     &      SJ,SORDERJ)
C.. we need to step through the lists of closed and open orbitals, seeing
C.. what has happened to each type of orbital.
C.. a closed orbital going to an open det is a single excitation.
C.. a closed orbital going to empty orbital is a double
C.. an open going to an empty is a single

         IEXLEVEL=0
C.. start with the closed
         I=1
         J=1
         K=1
         DO WHILE(I.LE.NCLOSEDI)
            DO WHILE(J.LE.NCLOSEDJ.AND.ICLOSEDJ(J).LT.ICLOSEDI(I))
               J=J+2
            ENDDO
            DO WHILE(K.LE.NOPENJ.AND.IOPENJ(K).LT.ICLOSEDI(I))
               K=K+1
            ENDDO   
            IF(J.LE.NCLOSEDJ.AND.ICLOSEDI(I).EQ.ICLOSEDJ(J)) THEN
C.. closed has remained closed
            ELSEIF(K.LE.NOPENJ.AND.ICLOSEDI(I).EQ.IOPENJ(K)) THEN
C.. closed to open
               IEXLEVEL=IEXLEVEL+1
            ELSE
C.. it must have gone to an empty then
               IEXLEVEL=IEXLEVEL+2
            ENDIF
            I=I+2
         ENDDO
C.. now do the open
         I=1
         J=1
         K=1
         DO WHILE(I.LE.NOPENI)
            DO WHILE(J.LE.NCLOSEDJ.AND.ICLOSEDJ(J).LT.IOPENI(I))
               J=J+2
            ENDDO
            DO WHILE(K.LE.NOPENJ.AND.IOPENJ(K).LT.IOPENI(I))
               K=K+1
            ENDDO   
            IF(J.LE.NCLOSEDJ.AND.IOPENI(I).EQ.ICLOSEDJ(J)) THEN
C.. open has gone to closed
C               IEXLEVEL=IEXLEVEL+1
C.. Not an excitation, because the open orb will definitely be filled
            ELSEIF(K.LE.NOPENJ.AND.IOPENI(I).EQ.IOPENJ(K)) THEN
C.. open to open
            ELSE
C.. it must have gone to an empty then
               IEXLEVEL=IEXLEVEL+1
            ENDIF
            I=I+1
         ENDDO
         IF(IEXLEVEL.EQ.0) THEN
            IF(ICMPDETS(NCI,NCJ,NEL).NE.0) IEXLEVEL=1
         ENDIF
      END
      SUBROUTINE SYMGETSSZ(SSYM,S,SZ)
         IMPLICIT NONE
         INTEGER SSYM,S,SZ
         INCLUDE 'csf.inc'
         S=(SSYM+CSF_NBSTART/3)/CSF_NBSTART
         SZ=SSYM-S*CSF_NBSTART
         RETURN
      END
C.. Get a "normal basis fn" corresponding to a 
C.. CSF encoded electron.
C.. This gives the spatial basis with the relevant
C.. spin projection.  The total S contribution (*2*CSF_NBSTART)
C.. is put in SSYM
      SUBROUTINE GETUNCSFELEC(ICEL,IEL,SSYM)
         IMPLICIT NONE
         INTEGER ICEL,IEL,SSYM
         INCLUDE 'csf.inc'
         IF(ICEL.GE.CSF_NBSTART) THEN
            IEL=IAND((ICEL-CSF_NBSTART)/4,65534)
            IEL=IEL+1+IAND(ICEL-CSF_NBSTART,1)
            SSYM=CSF_NBSTART*(IAND(ICEL-CSF_NBSTART,2)-1)
         ELSE
            IEL=ICEL
            SSYM=0
         ENDIF
         RETURN
            
      END

C.. test if NCI is a CSF
      LOGICAL FUNCTION ISCSF(NCI,NEL)
         IMPLICIT NONE
         INTEGER NEL,NCI(NEL)
         INCLUDE 'csf.inc'
         ISCSF=NCI(NEL).GE.CSF_NBSTART
      END

      SUBROUTINE GETNEXTSYMCSF(NEL,NBASIS,BRR,NBASISMAX,
     &  G1,ISYM,NI)
         IMPLICIT NONE
         INTEGER NEL,NBASIS,BRR(NBASIS),NBASISMAX(*)
         INTEGER NMRKS(NEL,1),NDET,G1(*)
         LOGICAL TCOUNT,TSPN,TPARITY
         INTEGER LMS,IFDET
         INCLUDE 'basis.inc'
         TYPE(BASISFN) IPARITY

         INTEGER ORBS(NBASIS/2,2)
         INTEGER NORBS

         INTEGER I,J,IELEC,IORB,IORBTYPE,IORBTYPE2
         INTEGER NI(NEL)
         LOGICAL TBR,BR2,TDO
         INTEGER IBACK
         INTEGER NCLOSED
         INTEGER NOPEN
         TYPE(BASISFN) ISYM,KSYM
         INTEGER SFORCE,SZFORCE
         INTEGER IEXCITS
         INCLUDE 'csf.inc'
         LOGICAL LCHKSYM
         CALL SYMGETSSZ(LMS,SFORCE,SZFORCE)
         CALL ICOPY(BASISFNSIZE,IPARITY,1,ISYM,1)
C.. make sure Sz=0 as this is dealt with in the csfgen part
         ISYM%MS=0
         IFDET=1
         NDET=0
         NORBS=NBASIS/2
         IEXCITS=0
C.. Set all orbitals to be available
         DO I=1,2
            DO J=1,NORBS
               ORBS(J,I)=1
            ENDDO
         ENDDO
         TBR=.TRUE.
C.. Start at the first orbital
         IORB=1
C.. and make it closed shell if we can
         IF(NEL.GT.1) THEN
            IORBTYPE=1
         ELSE
            IORBTYPE=2
         ENDIF
C.. and start puttincg it for the 1st electron
         IELEC=1
         NCLOSED=0
         NOPEN=0
         DO WHILE (TBR)
            IF(IORBTYPE.EQ.1) THEN
C.. Closed shell
               NI(IELEC)=IORB*2-1
               NI(IELEC+1)=IORB*2

C.. mark both closed and open shell orbs unavailable
               ORBS(IORB,1)=0
               ORBS(IORB,2)=0
C.. move up two electrons
               IELEC=IELEC+2
               NCLOSED=NCLOSED+2
            ELSE
C.. we have an open shell electron
               NI(IELEC)=CSF_NBSTART+8*(IORB-1)
               IELEC=IELEC+1
               ORBS(IORB,1)=0
               ORBS(IORB,2)=0
               IORB=IORB+1
               NOPEN=NOPEN+1
            ENDIF
C.. First check if we have enough electrons for a CSF
            IF(IELEC.GT.NEL) THEN
C.. Now check the symmetry's right
C.. we'll deal with spin in the CSF generation
               CALL GETSYM(NI,NEL,G1,NBASISMAX,KSYM)
C.. remove Ms value, as this is dealt with in csfgencsfs
               KSYM%MS=0
               IF(.NOT.TPARITY.OR.LCHKSYM(ISYM,KSYM))
C.. go and create a set of CSFs from these electrons
     &           CALL CSFGENCSFS(NI,NEL,NOPEN,NCLOSED,NMRKS,NDET,TCOUNT,
     &            TSPN,SFORCE,SZFORCE)
            ENDIF
C.. We've now added some electrons, so we need to now check that the new
C.. places our orbitals are pointing at have not been used up already.
            BR2=.TRUE.
            IBACK=0
            IF(IELEC.GT.NEL) IBACK=1
            DO WHILE (BR2)
               IF(IBACK.EQ.0) THEN
C.. if we're on the last one of an odd number of electrons, it must be
C.. open
                  IF(IELEC.EQ.NEL.AND.IAND(NEL,1).EQ.1
     &              .AND.IORBTYPE.EQ.1) THEN
                     IORBTYPE=2
                     IORB=1
                  ENDIF 
                  DO WHILE(ORBS(IORB,IORBTYPE).EQ.0.AND.IORB.LE.NORBS)
                     IORB=IORB+1
                  ENDDO
C.. IORB now either points to NORB+1 or an available orbital
C.. If there are no more orbitals available, we have to go back one or
C.. two electrons and try there, or try the next orbtype
                  IF(IORB.GT.NORBS) THEN
                     IF(IORBTYPE.EQ.1) THEN
C.. Try the next orbtype
                        IORBTYPE=2
                        IORB=1
                     ELSE
                        IBACK=IBACK+1
                    ENDIF
                  ELSE
                     BR2=.FALSE.
                  ENDIF
               ELSE
                  IBACK=IBACK-1
C.. go back an electron or two
C.. See what sort of electron the last one was
                  IF(IELEC.EQ.1) THEN
                     IELEC=0
                  ELSEIF(NI(IELEC-1).GE.CSF_NBSTART) THEN
C.. an open shell
                     IELEC=IELEC-1
                     NOPEN=NOPEN-1
                     IORB=(NI(IELEC)-CSF_NBSTART)/8+1
                     ORBS(IORB,1)=1
                     ORBS(IORB,2)=1
                     IORB=IORB+1
                     IORBTYPE=2
                  ELSE
C.. a closed shell
                     IELEC=IELEC-2
                     NCLOSED=NCLOSED-2
                     IORB=(NI(IELEC)+1)/2
                     ORBS(IORB,1)=1
                     ORBS(IORB,2)=1
                     IORB=IORB+1
                     IORBTYPE=1
                  ENDIF
               ENDIF
               IF(IELEC.EQ.0) THEN
C.. we've run out of electrons to go back entirely, so we're done
                  BR2=.FALSE.
                  TBR=.FALSE.
               ENDIF
            ENDDO
         ENDDO
         RETURN 
      END

C.. find I in sorted LST of size N, returning index or 0 if not there
      INTEGER FUNCTION IFIND(LST,N,I)
         IMPLICIT NONE
         INTEGER N,LST(N),I
         IFIND=1
         DO WHILE (IFIND.LE.N.AND.LST(IFIND).LT.I)
            IFIND=IFIND+1
         ENDDO
         IF(IFIND.GT.N.OR.LST(IFIND).NE.I) IFIND=0
         RETURN
      END
C.. Generate all the CSFs of a system
C.. TCHKEX set if checking whether excitations are at most of order
C.. IEXLEV from NDCH
      SUBROUTINE GNCSFS(NEL,NBASIS,BRR,NBASISMAX,
     &  NMRKS,TCOUNT,G1,TSPN,LMS,TPARITY,
     &  IPARITY,NDET,IFDET,TCHKEX,NDCH,ICLIST,TORDER,IEXLEV)
         IMPLICIT NONE
         INTEGER NEL,NBASIS,BRR(NBASIS),NBASISMAX(*)
         INTEGER NMRKS(NEL,*),NDET,G1(*)
         LOGICAL TCOUNT,TSPN,TPARITY
         INTEGER LMS,IFDET
         INCLUDE 'basis.inc'
         TYPE(BASISFN) IPARITY        
 
         LOGICAL TCHKEX,TORDER
         INTEGER NDCH(NEL),ICLIST(*)

         INTEGER BRR2(NBASIS)
         INTEGER IC
         INTEGER ORBS(NBASIS/2,2)
         INTEGER NORBS

         INTEGER I,J,IELEC,IORB,IORBTYPE,IORBTYPE2
         INTEGER NI(NEL),IND(NEL),NJ(NEL)
         LOGICAL TBR,BR2,TDO,TBR2
         INTEGER IBACK
         INTEGER NCLOSED
         INTEGER NOPEN
         TYPE(BASISFN) ISYM,KSYM
         INTEGER SFORCE,SZFORCE
         INTEGER IEXCITS(NEL+1),OIEX,IEXCITS2(NEL+1),OIEX2
         INTEGER ICLOSEDD(NEL),NCLOSEDD,IOPEND(NEL),NOPEND
         INTEGER SZD,SD,SORDERD(NEL)

         INCLUDE 'csf.inc'
         LOGICAL LCHKSYM
         INTEGER IGETEXCITLEVEL,IFIND
         INTEGER IEXLEV


         IF(TCHKEX) CALL CSFGETCSFDATA(NDCH,NEL,ICLOSEDD,NCLOSEDD,
     &      IOPEND,NOPEND,SZD,SD,SORDERD)

         
         DO I=1,NBASIS
            IF(TORDER) THEN
               BRR2(I)=I
            ELSE
               BRR2(I)=BRR(I)
            ENDIF
         ENDDO
         CALL SYMGETSSZ(LMS,SFORCE,SZFORCE)
         CALL ICOPY(BASISFNSIZE,IPARITY,1,ISYM,1)
C.. ensure Ms=0 as this is dealt with later
         ISYM%MS=0
         IFDET=1
         NDET=0
         NORBS=NBASIS/2
         IEXCITS(1)=0
         IEXCITS2(1)=0
C.. Set all orbitals to be available
         DO I=1,2
            DO J=1,NORBS
               ORBS(J,I)=1
            ENDDO
         ENDDO
         TBR=.TRUE.
C.. Start at the first orbital
         IORB=1
C.. and make it closed shell if we can
         IF(NEL.GT.1) THEN
            IORBTYPE=1
         ELSE
            IORBTYPE=2
         ENDIF
C.. and start putting it for the 1st electron
         IELEC=1
         NCLOSED=0
         NOPEN=0
         DO WHILE (TBR.AND..TRUE.)
            OIEX2=IEXCITS2(IELEC)
            OIEX=IEXCITS(IELEC)
            IF(IORBTYPE.EQ.1) THEN
C.. Closed shell
               IND(IELEC)=IORB*2-1
               IND(IELEC+1)=IORB*2
               NI(IELEC)=BRR2(IND(IELEC))
               NI(IELEC+1)=BRR2(IND(IELEC+1))
C.. mark both closed and open shell orbs unavailable
               ORBS(IORB,1)=0
               ORBS(IORB,2)=0
C.. move up two electrons
               IELEC=IELEC+2
               NCLOSED=NCLOSED+2
            ELSE
C.. we have an open shell electron
               IND(IELEC)=CSF_NBSTART+IORB*2-1
               NI(IELEC)=CSF_NBSTART+4*(BRR2(IORB*2-1)-1)
               IELEC=IELEC+1
               ORBS(IORB,1)=0
               ORBS(IORB,2)=0
               NOPEN=NOPEN+1
            ENDIF
            IF(TCHKEX) THEN
C.. see what excitation level we've gone to
C.. IEXCITS is the level of excitation from NDCH to NI, and
C.. IEXCITS2 is the level of excitation from NI to NDCH
               IF(IFIND(IOPEND,NOPEND,IORB*2-1).NE.0) THEN
C.. it is an open shell in NCHD
                  IF(IORBTYPE.EQ.1) THEN
C.. open->closed Not an excitation because the open will definitely
C.. still be filled.
                     IEXCITS(IELEC)=OIEX
                     IEXCITS2(IELEC)=OIEX2+1
                  ELSE 
                     IEXCITS(IELEC)=OIEX
                     IEXCITS2(IELEC)=OIEX2
                  ENDIF
               ELSEIF(IFIND(ICLOSEDD,NCLOSEDD,IORB*2-1).NE.0) THEN
C.. it's a closed shell in NCHD
                  IF(IORBTYPE.EQ.2) THEN
                     IEXCITS(IELEC)=OIEX+1
                     IEXCITS2(IELEC)=OIEX2
                  ELSE 
                     IEXCITS(IELEC)=OIEX
                     IEXCITS2(IELEC)=OIEX2
                  ENDIF
               ELSE
C.. empty in NDCH
                  IF(IORBTYPE.EQ.1) THEN
                     IEXCITS(IELEC)=OIEX
                     IEXCITS2(IELEC)=OIEX2+2
                  ELSE 
                     IEXCITS(IELEC)=OIEX
                     IEXCITS2(IELEC)=OIEX2+1
                  ENDIF
               ENDIF
            ENDIF
C.. move up one basis fn
               IORB=IORB+1
C.. First check if we have enough electrons for a CSF
            IF(IELEC.GT.NEL) THEN
C.. Now check the symmetry's right
C.. we'll deal with spin in the CSF generation
               CALL GETSYM(NI,NEL,G1,NBASISMAX,KSYM)
               KSYM%MS=0
               IF(.NOT.TPARITY.OR.LCHKSYM(ISYM,KSYM)) THEN
C.. go and create a set of CSFs from these electrons
                CALL ICOPY(NEL,NI,1,NJ,1)
                IF(.NOT.TORDER) CALL SORTI(NEL,NJ)
                IF(TCHKEX) IC=IGETEXCITLEVEL(NDCH,NJ,NEL)
                IF(.NOT.TCHKEX.OR.IC.LE.IEXLEV)
     &           CALL CSFGENCSFS(NJ,NEL,NOPEN,NCLOSED,NMRKS,NDET,TCOUNT,
     &            TSPN,SFORCE,SZFORCE,TCHKEX,ICLIST,IC,NDCH)
               ENDIF
            ENDIF
C.. We've now added some electrons, so we need to now check that the new
C.. places our orbitals are pointing at have not been used up already.
            BR2=.TRUE.
            IBACK=0
            IF(IELEC.GT.NEL) IBACK=1
C.. if we've done too many excitations, go back
C            IF(TCHKEX.AND.(IEXCITS(IELEC).GT.2.OR.IEXCITS2(IELEC).GT.2))
C     &          IBACK=1
            DO WHILE (BR2)
               IF(IBACK.EQ.0) THEN
C.. if we're on the last one of an odd number of electrons, it must be
C.. open
                  IF(IELEC.EQ.NEL.AND.IAND(NEL,1).EQ.1
     &              .AND.IORBTYPE.EQ.1) THEN
                     IORBTYPE=2
                     IORB=1
                  ENDIF 
                  DO WHILE(ORBS(IORB,IORBTYPE).EQ.0.AND.IORB.LE.NORBS)
                     IORB=IORB+1
                  ENDDO
C.. IORB now either points to NORB+1 or an available orbital
C.. If there are no more orbitals available, we have to go back one or
C.. two electrons and try there, or try the next orbtype
                  IF(IORB.GT.NORBS) THEN
                     IF(IORBTYPE.EQ.1) THEN
C.. Try the next orbtype
                        IORBTYPE=2
                        IORB=1
                     ELSE
                        IBACK=IBACK+1
                    ENDIF
                  ELSE
                     BR2=.FALSE.
                  ENDIF
               ELSE
                  IBACK=IBACK-1
C.. go back an electron or two
C.. See what sort of electron the last one was
                  IF(IELEC.EQ.1) THEN
                     IELEC=0
                  ELSEIF(NI(IELEC-1).GE.CSF_NBSTART) THEN
C.. an open shell
C.. remove the spatial sym of this electron - no need, as we roll back the sym
C                     CALL ADDELECSYM(NI(IELEC-1),G1,
C     &                  NBASISMAX,JSYM,.TRUE.)
                     IELEC=IELEC-1
                     NOPEN=NOPEN-1
                     IORB=(IND(IELEC)-CSF_NBSTART)/2+1
                     ORBS(IORB,1)=1
                     ORBS(IORB,2)=1
                     IORB=IORB+1
                     IORBTYPE=2
                  ELSE
C.. a closed shell
C.. remove the spatial symmetries of these electrons - no need as we roll back sym
C               CALL ADDELECSYM(NI(IELEC-1),G1,NBASISMAX,JSYM,.TRUE.)
C               CALL ADDELECSYM(NI(IELEC-2),G1,NBASISMAX,JSYM,.TRUE.)
                     IELEC=IELEC-2
                     NCLOSED=NCLOSED-2
                     IORB=(IND(IELEC)+1)/2
                     ORBS(IORB,1)=1
                     ORBS(IORB,2)=1
                     IORB=IORB+1
                     IORBTYPE=1
                  ENDIF
               ENDIF
               IF(IELEC.EQ.0) THEN
C.. we've run out of electrons to go back entirely, so we're done
                  BR2=.FALSE.
                  TBR=.FALSE.
               ENDIF
            ENDDO
         ENDDO
         RETURN 
      END

C.. From a given closed/open orbital configuration, generate all allowed
C.. CSFs in order
      SUBROUTINE CSFGENCSFS(NI,NEL,NOPEN,NCLOSED,NMRKS,NDET,TCOUNT,
     &            TSPN,SFORCE,SZFORCE,TCHKEX,ICLIST,IC,NDCH)
         IMPLICIT NONE
         INTEGER NI(NEL),NEL,NOPEN,NCLOSED
         INTEGER NMRKS(NEL,*),NDET
         LOGICAL TCOUNT,TSPN
         INTEGER SFORCE,SZFORCE
         LOGICAL TCHKEX
         INTEGER ICLIST(*),IC,NDCH(NEL)
         
C.. the spin order of the current det, and
C.. the DORDER (list of betas then alphas which gives the total M
         INTEGER SORDER(NOPEN),DORDER(NOPEN)
         INTEGER SORDERC(0:NOPEN),DORDERC(0:NOPEN)
         INTEGER S,M,SC,MC
         INTEGER I,J         
         INTEGER IELEC
         LOGICAL TBR,TGEN,TADD,TBR2,TCLOSED

         INTEGER ICMPDETS
         
         TCLOSED=(NOPEN.EQ.0)
         IF(TCLOSED.AND.TSPN.AND.(SZFORCE.NE.0.OR.SFORCE.NE.0)) RETURN
         TBR=.TRUE.
C         IF(TCHKEX) THEN
C            CALL WRITEDET(6,NDCH,NEL,.FALSE.)
C            CALL WRITEDET(6,NI,NEL,.FALSE.)
C            WRITE(6,*) IC
C         ENDIF
C         IF(TCHKEX.AND.IC.EQ.0.AND.NOPEN.EQ.4) THEN
C            TBR=.TRUE.
C         ENDIF
         IELEC=1
         TGEN=.TRUE.
         DORDERC(0)=0
         SORDERC(0)=0
         DO WHILE (IELEC.GT.0)
            IF(TGEN) THEN
C.. we have to reinitialize the SORDER and DORDER arrays from electron IELEC
C.. we traverse down the list, finding the lowest possible combination
C.. each stage.
               DO WHILE(IELEC.LE.NOPEN)
                  IF(IELEC.GT.1) THEN
                     S=SORDERC(IELEC-1)
                     M=DORDERC(IELEC-1)
                  ELSE
                     S=0
                     M=0
                  ENDIF
                  IF(S.GT.0) THEN
                     SORDER(IELEC)=0
                  ELSE
C.. S cannot go below 0
                     SORDER(IELEC)=1
                  ENDIF
                  S=S+2*SORDER(IELEC)-1
                  SORDERC(IELEC)=S
                  DORDER(IELEC)=0
C.. We cannot have a beta after we've had an alpha
C.. (the beta/alpha string must be strictly ascending
C.. to ensure it is unique)
                  IF(IELEC.GT.1) DORDER(IELEC)=DORDER(IELEC-1)
C.. We can force the spin projection at this point.
                  IF(TSPN) THEN
C.. if we're past our allocated betas, force an alpha
                     IF(IELEC.GT.(NOPEN-SZFORCE)/2)
     &                  DORDER(IELEC)=1
                  ENDIF
C.. We cannot have |M|>S
C                  IF(M.LE.-S) DORDER(IELEC)=1
                  M=M+2*DORDER(IELEC)-1
                  DORDERC(IELEC)=M
                  IELEC=IELEC+1
               ENDDO
C.. We've initialised to the end, so set the current electron to the end
               IELEC=NOPEN
               TGEN=.FALSE.
               TADD=.TRUE.
               IF(ABS(M).GT.S) TADD=.FALSE.
C.. We're also guaranteed to have found an allowed csf (not quite)
            ELSE
C.. we haven't reinitialised, so we just increment our current place by 1
               J=SORDER(IELEC)*2+DORDER(IELEC)
               TBR2=.TRUE.
               DO WHILE(TBR2)
                  J=J+1
C.. Test to see if our current +/- alpha/beta setting is allowed
                  SC=IAND(J,2)-1
                  MC=IAND(J,1)*2-1
C.. We cannot have a beta after we've had an alpha
                  IF(IELEC.GT.1.AND.IAND(J,1).LT.DORDER(IELEC-1)) THEN
                     MC=MC+2
                     J=J+1
                  ENDIF
               IF(TSPN.AND.MC.EQ.1.AND.IELEC.LE.(NOPEN-SZFORCE)/2) THEN
C.. If we're trying to put an alpha when there's only a beta allowed
C.. Flip round to the next SC value with a beta
                     J=J+1
                     SC=SC+2
                     MC=MC-2
                  ENDIF
C.. if we've rolled over, we need to go back
                  IF(J.GE.4) THEN
                     TADD=.FALSE.
                     IELEC=IELEC-1
                     TBR2=.FALSE.
C.. If we go back to IELEC=0 our uberloop ends
                  ELSE
                     S=SORDERC(IELEC-1)+SC
                     M=DORDERC(IELEC-1)+MC
C                     IF(ABS(M).LE.S) 
                     TBR2=.FALSE.
                  ENDIF
               ENDDO
               IF(J.LT.4) THEN
C.. We've found an allowed state, so we reinitialize the further
C.. electrons
                  SORDER(IELEC)=(SC+1)/2
                  DORDER(IELEC)=(MC+1)/2
                  SORDERC(IELEC)=S
                  DORDERC(IELEC)=M
                  IF(IELEC.LT.NOPEN) THEN
                     TADD=.FALSE.
                     TGEN=.TRUE.
                     IELEC=IELEC+1
                  ELSE
                     IF(ABS(M).LE.S) THEN
                        TADD=.TRUE.
                     ELSE
                        TADD=.FALSE.
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            IF(TADD) THEN
C.. Check the spin is correct if we're forcing it
               IF(TCLOSED.OR..NOT.TSPN.OR.(DORDERC(IELEC).EQ.SZFORCE
     &            .AND.SORDERC(IELEC).EQ.SFORCE)) THEN
C.. We have found a good combination.  Add it
                  NDET=NDET+1
                  IF(.NOT.TCOUNT) THEN
C.. Copy the closed orbitals
                     CALL ICOPY(NCLOSED,NI,1,NMRKS(1,NDET),1)
C.. now put in the open orbitals
                     DO I=1,NOPEN
                        NMRKS(NCLOSED+I,NDET)=NI(NCLOSED+I)+SORDER(I)*2
     &                                       +DORDER(I)
                     ENDDO
C                     CALL WRITEDET(6,NMRKS(1,NDET),NEL,.TRUE.)
C                     IF(NDET.GT.1.AND.ICMPDETS(NMRKS(1,NDET),
C     &                  NMRKS(1,NDET-1),NEL).LT.1) THEN
C                        CALL WRITEDET(6,NMRKS(1,NDET-1),NEL,.TRUE.)
C                        CALL WRITEDET(6,NMRKS(1,NDET),NEL,.TRUE.)
C                        WRITE(6,*) NDET,ICMPDETS(NMRKS(1,NDET),
C     &                     NMRKS(1,NDET-1),NEL)
C                        STOP "Sort Error"
C                     ENDIF
                     IF(TCHKEX) THEN
                        ICLIST(NDET)=IC
C.. remove NDCH if necessary
                        IF(IC.LE.1
     &                     .AND.ICMPDETS(NMRKS(1,NDET),NDCH,NEL).EQ.0)
     &                        NDET=NDET-1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
C.. If we're only counting, we haven't had a change to exclude NI, so we
C.. do it now.
         IF(TCHKEX.AND.TCOUNT.AND.IC.EQ.0) NDET=NDET-1
         RETURN
      END


C.. Get the Hamiltonian Element between two CSFs
      SUBROUTINE CSFGETHELEMENT(NCI,NCJ,NEL,NBASISMAX,
     &   G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,ECORE,SUM)

         IMPLICIT NONE
         INTEGER G1(*),NBASIS,BRR(*),NMSH,NMAX
         REAL*8 FCK(*),ZIA(*),ALAT(*),UMAT(*)
         INTEGER NEL,NCI(NEL),NCJ(NEL),NBASISMAX(5,2)
         REAL*8 SUM,SUM2,ECORE
         REAL*8 GETHELEMENT2
         
         POINTER (IP_LISTD1,LISTD1),(IP_DOLISTD1,DOLISTD1)
         POINTER (IP_LISTD2,LISTD2),(IP_DOLISTD2,DOLISTD2)
         INTEGER LISTD1(NEL,*),LISTD2(NEL,*)
         INTEGER DOLISTD1(*),DOLISTD2(*)
         INTEGER NOPEN1,NOPEN2,NDETS1,NDETS2
         INTEGER SO1(NEL),SO2(NEL)
         INTEGER I,J,ISUB
         REAL*8 COEFF1,COEFF2,HEL
         REAL*8 COEFFLIST2(*)
         POINTER (IP_COEFFLIST2,COEFFLIST2)
         
         LOGICAL CSFCHKSPINSYM
         IF(.NOT.CSFCHKSPINSYM(NCI,NCJ,NEL)) THEN
C.. If the Ms or S values differ, we return 0.
            SUM=0.D0
            RETURN
         ENDIF
C.. The SORDERS for dets 1 and 2
         CALL TISET('CSFGETHELE',ISUB)

         NOPEN1=0
         NOPEN2=0
         SO1(1)=1
         SO2(1)=1




         CALL CSFGENOPENDETS(NCI,NEL,LISTD1,DOLISTD1,SO1,NOPEN1,NDETS1)
         CALL CSFGENOPENDETS(NCJ,NEL,LISTD2,DOLISTD2,SO2,NOPEN2,NDETS2)

C.. NOPEN has number of open electrons, and NDETS has #open dets
         CALL MEMORY(IP_LISTD1,NDETS1*NEL,'LISTD1')
         CALL MEMORY(IP_DOLISTD1,NDETS1*NOPEN1,'DOLISTD1')
         CALL MEMORY(IP_LISTD2,NDETS2*NEL,'LISTD2')
         CALL MEMORY(IP_DOLISTD2,NDETS2*NOPEN2,'DOLISTD2')
         CALL MEMORY(IP_COEFFLIST2,NDETS2,'COEFFLIST2')
         CALL CSFGENOPENDETS(NCI,NEL,LISTD1,DOLISTD1,SO1,1,NDETS1)
         CALL CSFGENOPENDETS(NCJ,NEL,LISTD2,DOLISTD2,SO2,1,NDETS2)

!         WRITE(6,*) "GHEL"
!         CALL WRITEDET(6,NCI,NEL,.TRUE.)
!         DO I=1,NDETS1
!            CALL CSFGETDETCOEFF(SO1,DOLISTD1((I-1)*NOPEN1+1),NOPEN1,
!     &         COEFF1)
!            CALL WRITEDET(6,LISTD1(1,I),NEL,.FALSE.)
!            WRITE(6,*) COEFF1
!         ENDDO
         
!         CALL WRITEDET(6,NCJ,NEL,.TRUE.)
!         DO I=1,NDETS2
!            CALL CSFGETDETCOEFF(SO1,DOLISTD2((I-1)*NOPEN2+1),NOPEN2,
!     &         COEFF2)
!            CALL WRITEDET(6,LISTD2(1,I),NEL,.FALSE.)
!            WRITE(6,*) COEFF2
!         ENDDO
         SUM=0.D0
         DO I=1,NDETS1
            CALL CSFGETDETCOEFF(SO1,DOLISTD1((I-1)*NOPEN1+1),NOPEN1,
     &         COEFF1)
            SUM2=0.D0
            DO J=1,NDETS2
               IF(I.EQ.1) THEN
                  CALL CSFGETDETCOEFF(SO2,DOLISTD2((J-1)*NOPEN2+1),
     &               NOPEN2,COEFF2)
!                  CALL WRITEDET(6,LISTD2(1,J),NEL,.FALSE.)
!                  WRITE(6,*) COEFF2
                  COEFFLIST2(J)=COEFF2
               ELSE
                  COEFF2=COEFFLIST2(J)
               ENDIF
               HEL=GETHELEMENT2(LISTD1(1,I),LISTD2(1,J),NEL,NBASISMAX,
     &            G1,NBASIS,BRR,NMSH,FCK,ZIA,NMAX,ALAT,UMAT,-1,ECORE)
               SUM2=SUM2+COEFF2*HEL
            ENDDO
            SUM=SUM+SUM2*COEFF1
         ENDDO
         CALL FREEM(IP_LISTD1)
         CALL FREEM(IP_DOLISTD1)
         CALL FREEM(IP_LISTD2)
         CALL FREEM(IP_DOLISTD2)
         CALL FREEM(IP_COEFFLIST2)
         CALL TIHALT('CSFGETHELE',ISUB)
               
      END

C.. Get the coefficient of the det specified by SORDER/DORDER.
C.. This is a product of Clebsch-Gordon coefficients
C.. (and is independent of which spatial orbitals SORDER and DORDER
C..  refer to as long as they are always ordered)
      SUBROUTINE CSFGETDETCOEFF(SORDER,DORDER,NOPEN,COEFF)
         IMPLICIT NONE
         REAL*8 COEFF,CG
         INTEGER NOPEN,SORDER(NOPEN),DORDER(NOPEN)
         INTEGER I,SCUR,MCUR
         REAL*8 S,M
         COEFF=1.D0
         S=0.D0
         M=0.D0
         DO I=1,NOPEN
            SCUR=SORDER(I)
            MCUR=DORDER(I)
            S=S+SCUR-0.5D0
            M=M+MCUR-0.5D0
            CALL GETCLEBSCHGORDON(S,M,SCUR,MCUR,CG)
!            WRITE(6,"(F6.3,A,$)") CG,"*"
            COEFF=COEFF*CG
            IF(COEFF.EQ.0.D0) EXIT
         ENDDO
         RETURN
      END

C.. Return in COEFF, the Clebsch-Gordon coefficient coupling state
C.. with total spin S, projected component of spin M, with state
C.. with spin change SN_=(+/-)1/2 and 
C.. projected spin change MN_=(+/-)1/2
C.. For SN_=+1/2 set SN=1, and for SN_=-1/2 set SN=0
C.. For MN_=+1/2 set MN=1, and for MN_=-1/2 set MN=0
      SUBROUTINE GETCLEBSCHGORDON(S,M,SN,MN,COEFF)
         IMPLICIT NONE
         REAL*8 S,M,COEFF
         INTEGER SN,MN

         REAL*8 MN_
         MN_=MN-0.5D0
         IF(SN.EQ.0) THEN
            COEFF=-2.D0*MN_*SQRT((S+1-2*MN_*M)/(2*S+2))
         ELSEIF(SN.EQ.1) THEN
            COEFF=SQRT((S+2*MN_*M)/(2*S))
         ELSE
            COEFF=0.D0
         ENDIF
      END

C.. Given a CSF in NCI, generate the list of dets which make it up
C.. if ISTART=0, return in IEND the number of dets, and ISTART as NOPEN
C.. otherwise, place dets ISTART to IEND inclusive in LIST
C.. (the first det is numbered 1)
C.. DORDERLIST is a list of DORDER vectors corresponding to the 
C.. spin-projection composition of each determinant
C.. if SORDER2(1) is non-zero, we store SORDER in there on return
      SUBROUTINE CSFGENOPENDETS(NCI,NEL,DETLIST,DORDERLIST,SORDER2,
     &   ISTART,IEND)
         IMPLICIT NONE
         INTEGER NEL,ISTART,IEND
         INTEGER NCI(NEL),DETLIST(NEL,*)
         INTEGER DORDERLIST(*)

         INTEGER NOPEN,IOPEN(NEL),NCLOSED,ICLOSED(NEL)
         INTEGER SZ,S,SORDER(NEL),SORDER2(NEL),ICOUNT,IEN
         INTEGER INDICES(NEL),NBETA,I,J
         LOGICAL TDONE

         INTEGER BINOMIAL
         CALL CSFGETCSFDATA(NCI,NEL,ICLOSED,NCLOSED,IOPEN,NOPEN,SZ,S,
     &      SORDER)
         IF(ISTART.EQ.0) THEN
C.. return the number of dets.
C.. SZ is twice the spin proj (Sz), and we need (nopen \choose nopen/2(+/-)Sz
            IEND=BINOMIAL(NOPEN,(NOPEN+SZ)/2)
            ISTART=NOPEN
            RETURN
         ENDIF
C.. Otherwise we need to just list the dets.  We do this by picking the
C.. NOPEN Beta orbitals in order, and then the remaining orbitals must
C.. be alpha.
         ICOUNT=1
         IEN=IEND-ISTART+1
         NBETA=(NOPEN-SZ)/2
C.. Determine the set of indices corresponding to the start det
         CALL CSFGETCHOOSEINDICES(ISTART,INDICES,NOPEN,NOPEN)
         TDONE=.FALSE.
         DO WHILE(ICOUNT.LE.IEN.AND..NOT.TDONE)
C.. Copy the closed shell orbitals, and then the open shell ones we've
C.. indexed (they being beta), interleaving the alphas of those
C.. which have not been indexed
            CALL ICOPY(NCLOSED,ICLOSED,1,DETLIST(1,ICOUNT),1)
            J=1
C.. J is the next indexed orbital to choose
C.. Also fill out the DORDERLIST
            DO I=1,NOPEN
               IF(J.LE.NBETA.AND.INDICES(J).EQ.I) THEN
C.. we've found a beta
                  DETLIST(NCLOSED+I,ICOUNT)=IOPEN(I)
                  DORDERLIST((ICOUNT-1)*NOPEN+I)=0
                  J=J+1
               ELSE
C.. an alpha      
                  DORDERLIST((ICOUNT-1)*NOPEN+I)=1
                  DETLIST(NCLOSED+I,ICOUNT)=IOPEN(I)+1
               ENDIF
            ENDDO
            CALL SORTI(NEL,DETLIST(1,ICOUNT))
            ICOUNT=ICOUNT+1
            CALL CSFINCCHOOSEINDICES(INDICES,NOPEN,NBETA,TDONE)
         ENDDO
         IF(ICOUNT.LT.IEN) STOP "NOT ENOUGH CSF determinants found"
         IF(SORDER2(1).NE.0) CALL ICOPY(NOPEN,SORDER,1,SORDER2,1)
         RETURN
      END

C..   Given a NITEMS from which to choose NCHOOSE, we list all possible
C..   choices in order, indexing them by ICHOOSE.
C..   This generates a list of the items (in NINDICES) correposnding to
C..   the choice indexed by ICHOOSE.
C..   We can use a method based on traversing a tree (21/4/05)
      SUBROUTINE CSFGETCHOOSEINDICES(ICHOOSE,INDICES,NITEMS,NCHOOSE)
         IMPLICIT NONE
         INTEGER ICHOOSE, NITEMS,NCHOOSE,INDICES(NCHOOSE)

         INTEGER BINOMIAL

         INTEGER N,R,INDX,BIN

         IF(ICHOOSE.GT.NCHOOSE) RETURN
         IF(ICHOOSE.GT.BINOMIAL(NITEMS,NCHOOSE))
     &     STOP "ERROR in CSFGETCHOOISEINDICES - IND too large"
         R=NCHOOSE
         N=NITEMS
         INDX=ICHOOSE
         DO WHILE (R.GT.0)
            BIN=BINOMIAL(N-1,R)
            IF(INDX.GT.BIN) THEN
               INDX=INDX-BIN
               INDICES(R)=N
               R=R-1
            ENDIF
            N=N-1
         ENDDO
      END

C.. Take an INDICES list, and move to the next set in order.
C.. TDONE will be set to TRUE if we've reached the end
      SUBROUTINE CSFINCCHOOSEINDICES(INDICES,NITEMS,NCHOOSE,TDONE)
         IMPLICIT NONE
         INTEGER NITEMS,NCHOOSE,INDICES(NCHOOSE)
         LOGICAL TDONE

         LOGICAL TD2
         INTEGER R,I
         
         R=1
         TD2=.TRUE.
         DO WHILE(TD2)
            IF(R.EQ.NCHOOSE.AND.INDICES(R).EQ.NITEMS) THEN
C.. if we're required to move the last choice past the last item, we 
C.. go off and cry
               TDONE=.TRUE.
               TD2=.FALSE.
            ELSEIF(R.EQ.NCHOOSE.OR.INDICES(R).LT.INDICES(R+1)-1) THEN
C.. Otherwise we see if we can move this choice to the next item.
C.. Do this move, and then reset the lower R to the first R items
               INDICES(R)=INDICES(R)+1
               DO I=1,R-1
                  INDICES(I)=I
               ENDDO
               TD2=.FALSE.
               TDONE=.FALSE.
            ELSE
C.. Failing that, we move to the next choice
               R=R+1
            ENDIF
         ENDDO
         RETURN
      END
C.. take a CSF and from it determine the closed and open orbitals
C.. 
C.. NCLOSED is the number of closed spinorbitals.
C.. ICLOSED will be returned containing a list of close spinorbitals, listed in 
C..   beta/alpha/beta/alpha... order (ascending in basis set numbering).
C.. NOPEN is the number of open orbitals
C.. IOPEN contains the beta-versions of the open orbitals
C.. SZ is the z-projection of the total spin (of both the CSF and the
C..   dets in it) (multiplied by 2)
C.. S  is the total S value (spin) of the CSF (multiplied by 2)
C.. SORDER is a string of 1's and 0's corresponding to the genealogical
C..   method of creating this CSF - 1 is a + and 0 is a -.

      SUBROUTINE CSFGETCSFDATA(NCI,NEL,ICLOSED,NCLOSED,IOPEN,NOPEN,SZ,S,
     &   SORDER)
         IMPLICIT NONE
         INTEGER NCI(NEL),NEL,ICLOSED(NEL),NCLOSED
         INTEGER IOPEN(NEL),NOPEN,SZ,S,SORDER(NEL)
         INCLUDE 'csf.inc'
         
         INTEGER I,J,K
         NCLOSED=0
         NOPEN=0
         SZ=0
         S=0
         DO I=1,NEL
            J=NCI(I)
            IF(J.LT.CSF_NBSTART) THEN
C.. we have a normal orbital
               NCLOSED=NCLOSED+1
               ICLOSED(NCLOSED)=NCI(I)
C.. check and alpha comes after a beta
               IF(NCLOSED.GT.1.AND.IAND(J,1).EQ.0
     &            .AND.ICLOSED(NCLOSED-1).NE.J-1) THEN
                     CALL WRITEDET(6,NCI,NEL,.TRUE.)
                     STOP "CSF has unpaired 'closed orbitals'"
               ENDIF
            ELSE
C.. an open orbital
               NOPEN=NOPEN+1
               K=J-CSF_NBSTART
               IOPEN(NOPEN)=K/4+1
C.. encode bit 1 of K as SORDER (-/+)=(0/1)
               SORDER(NOPEN)=IAND(K,2)/2
C.. add the +/- to total spin
               S=S+SORDER(NOPEN)*2-1
C.. encode bit 0 of K as BETA/Alpha=(0/1)
               SZ=SZ+IAND(K,1)*2-1
            ENDIF
         ENDDO
      END

C.. Generate an SINDEX (a unique ID corresponding to the SORDER list)
C.. from SORDER (the list of 0/1 = -/+ which specifies a particular 
C.. CSF.
C.. NOPEN is the number of open-shell electrons and the length of SORDER
C.. S is 2* the total spin
C.. we loop over increasing N in the genealogical tree, only including
C.. paths which go below or along the one we have (as paths above us 
C.. will have a higher index than us)

      SUBROUTINE CSFGETSINDEX(SORDER,NOPEN,SINDEX,S)
         IMPLICIT NONE
         INTEGER NOPEN,SORDER(NOPEN),S
         INTEGER SLEVEL(0:NOPEN,0:S)
         INTEGER SINDEX
C.. SLEVEL contains the relevant part of the genealogical tree
         INTEGER I,J,CUR
         SLEVEL(0,0)=0
         DO I=1,NOPEN
C.. I is the number of electrons we've added so far.
            DO J=0,MIN(S,I)
C.. J corresponds to the S reached so far
               CUR=0
               IF(J.GT.0) CUR=CUR+SLEVEL(I-1,J-1)
               IF(J.LT.S) CUR=CUR+SLEVEL(I-1,J+1)
               SLEVEL(I,J)=CUR
            ENDDO
         ENDDO
         SINDEX=SLEVEL(NOPEN,S)
      END
