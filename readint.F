      SUBROUTINE INITFROMFCID(NEL,NBASISMAX,LEN,LMS)
         IMPLICIT NONE
         INTEGER NEL,NBASISMAX(5,3),LEN,LMS
         INTEGER NORB,NELEC,MS2,ORBSYM(1000),ISYM
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM
         OPEN(8,FILE='FCIDUMP',STATUS='OLD' )
         READ(8,FCI)
         CLOSE(8)
         IF(NELEC.NE.NEL)  WRITE(6,*) 
     &      '*** WARNING: NEL in FCIDUMP differs from input file ***'
C         NEL=NELEC
         IF(LMS.NE.MS2)   WRITE(6,*)
     &      '*** WARNING: LMS in FCIDUMP differs from input file ***'
         LMS=MS2
         WRITE(6,*) ' NUMBER OF ELECTRONS : ' , NEL
         WRITE(6,*) ' BASIS MS : ' , LMS
         LEN=2*NORB
         CALL IAZZERO(NBASISMAX,12)
         NBASISMAX(1,1)=1
         NBASISMAX(1,2)=NORB
         NBASISMAX(1,3)=2
         NBASISMAX(4,1)=-1
         NBASISMAX(4,2)=1

C.. SHow that there's no momentum conservation
         NBASISMAX(3,3)=1
         RETURN
      END
      SUBROUTINE GETFCIBASIS(NBASISMAX,ARR,BRR,G1,LEN)
         IMPLICIT NONE
         INTEGER NBASISMAX(5,3),BRR(LEN),G1(5,LEN),LEN
         REAL*8 ARR(LEN)
         REAL*8 Z
         INTEGER I,J,K,L
         INTEGER ISYMNUM,ISNMAX
         INTEGER NORB,NELEC,MS2,ORBSYM(1000),ISYM
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM
         OPEN(8,FILE='FCIDUMP',STATUS='OLD')
         READ(8,FCI)
         ISYMNUM=0
         ISNMAX=0
         CALL IAZZERO(G1,5*LEN)
1        READ(8,'(1X,G20.12,4I3)',END=99) Z,I,J,K,L
C.. Each orbital in the file corresponds to alpha and beta spinorbitals
         IF(I.NE.0.AND.K.EQ.0.AND.I.EQ.J) THEN
            IF(I.GT.1) THEN
               IF(ORBSYM(I).NE.ORBSYM(I-1)) THEN
                  IF(ISYMNUM.GT.ISNMAX) ISNMAX=ISYMNUM
                  ISYMNUM=0
               ENDIF
            ENDIF
            ISYMNUM=ISYMNUM+1
            ARR(2*I-1)=Z
            ARR(2*I)=Z
         ENDIF
C.. At the moment we're ignoring the core energy
         GOTO 1
99       CONTINUE
         DO I=1,LEN
            BRR(2*I-1)=2*I-1
            BRR(2*I)=2*I
C            G1(1,2*I-1)=ISYMNUM-1
C            G1(1,2*I)=ISYMNUM-1
            G1(5,2*I-1)=IBSET(0,ORBSYM(I)-1)
            G1(5,2*I)=IBSET(0,ORBSYM(I)-1)
C            G1(3,2*I-1)=I
C            G1(3,2*I)=I
C.. set momentum to 0
            G1(1,2*I-1)=0
            G1(2,2*I)=0
            G1(3,2*I-1)=0
            G1(1,2*I)=0
            G1(2,2*I-1)=0
            G1(3,2*I)=0
            G1(4,2*I-1)=-1
            G1(4,2*I)=1
         ENDDO
         IF(ISNMAX.EQ.1) THEN
         ELSEIF(ISNMAX.LE.2) THEN
         ELSEIF(ISNMAX.LE.4) THEN
            ISNMAX=4
         ELSEIF(ISNMAX.LE.8) THEN
            ISNMAX=8
         ENDIF
         NBASISMAX(1,1)=0
         NBASISMAX(1,2)=ISNMAX-1
         NBASISMAX(2,1)=0
         NBASISMAX(2,2)=ORBSYM(NORB)-1
         CLOSE(8)
         RETURN
      END
      SUBROUTINE READFCIINT(UMAT,TMAT,NBASIS,ECORE)
         IMPLICIT NONE
         INTEGER NBASIS
         REAL*8 TMAT(NBASIS,NBASIS),ECORE
         REAL*8 UMAT(NBASIS/2,NBASIS/2,NBASIS/2,NBASIS/2)
         REAL*8 Z
         INTEGER I,J,K,L
         INTEGER NORB,NELEC,MS2,ORBSYM(1000),ISYM
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM
         OPEN(8,FILE='FCIDUMP',STATUS='OLD')
         READ(8,FCI)
101      READ(8,'(1X,G20.12,4I3)',END=199) Z,I,J,K,L
C.. Each orbital in the file corresponds to alpha and beta spinorbitals
         IF(K.NE.0) THEN
C.. 2-e integrals
C.. UMAT is stored as just spatial orbitals (not spinorbitals)
C..  we're reading in (IJ|KL), but we store <..|..> which is <IK|JL>
CC            UMAT(I,K,J,L)=Z
CC            UMAT(J,L,I,K)=Z
CC            UMAT(K,I,L,J)=Z
CC            UMAT(L,J,K,I)=Z
CC            UMAT(J,K,I,L)=Z
C            UMAT(I,K,J,L)=Z
CC            UMAT(L,I,K,J)=Z
CC            UMAT(K,J,L,I)=Z

            UMAT(I,K,J,L)=Z
            UMAT(J,L,I,K)=Z
            UMAT(K,I,L,J)=Z
            UMAT(L,J,K,I)=Z
            UMAT(J,K,I,L)=Z
            UMAT(I,L,J,K)=Z
            UMAT(L,I,K,J)=Z
            UMAT(K,J,L,I)=Z
         ELSEIF(I.NE.0) THEN
C.. 1-e integrals
C.. These are stored as spinorbitals
            TMAT(2*I-1,2*J-1)=Z
C            TMAT(2*I-1,2*J)=Z
C            TMAT(2*I,2*J-1)=Z
            TMAT(2*I,2*J)=Z
            TMAT(2*J-1,2*I-1)=Z
C            TMAT(2*J-1,2*I)=Z
C            TMAT(2*J,2*I-1)=Z
            TMAT(2*J,2*I)=Z
         ELSE
            ECORE=Z
         ENDIF
         GOTO 101
199      CONTINUE
         CLOSE(8)
         RETURN
      END

