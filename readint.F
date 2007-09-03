      SUBROUTINE INITFROMFCID(NEL,NBASISMAX,LEN,LMS)
         IMPLICIT NONE
         INTEGER NEL,NBASISMAX(5,3),LEN,LMS
         INTEGER NORB,NELEC,MS2,ORBSYM(1000),ISYM
!         USE UmatCache, only : TSTARBIN
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM
         !Will need to read in binary if TSTARBIN is set
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
         CALL IAZZERO(NBASISMAX,15)
C         NBASISMAX(1,1)=1
C         NBASISMAX(1,2)=NORB
         NBASISMAX(1,1)=0
         NBASISMAX(1,2)=0
         NBASISMAX(1,3)=2
         NBASISMAX(4,1)=-1
         NBASISMAX(4,2)=1

C.. SHow that there's no momentum conservation
         NBASISMAX(3,3)=1
         RETURN
      END
      SUBROUTINE GETFCIBASIS(NBASISMAX,ARR,BRR,G1,LEN)
         IMPLICIT NONE
         INTEGER NBASISMAX(5,3),BRR(LEN),LEN
         LOGICAL TSTARBIN
         include 'basis.inc'
         TYPE(BasisFN) G1(LEN)
         REAL*8 ARR(LEN,2)
         REAL*8 Z
         INTEGER I,J,K,L
         INTEGER ISYMNUM,ISNMAX,SYMMAX
         INTEGER NORB,NELEC,MS2,ORBSYM(1000),ISYM
!         USE UMatCache, only : TSTARBIN
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM
         !Will need to read in binary if TSTARBIN is set
         OPEN(8,FILE='FCIDUMP',STATUS='OLD')
         READ(8,FCI)
         ISYMNUM=0
         ISNMAX=0
         IF(LEN.NE.2*NORB) STOP 'LEN .NE. NORB in GETFCIBASIS'
         CALL IAZZERO(G1,BasisFNSize*LEN)
         CALL AZZERO(ARR,2*LEN)
1        READ(8,'(1X,G20.12,4I3)',END=99) Z,I,J,K,L
C.. Each orbital in the file corresponds to alpha and beta spinorbitals
         !Fill ARR with the energy levels
         IF(I.NE.0.AND.K.EQ.0.AND.I.EQ.J) THEN
            IF(I.GT.1) THEN
               IF(ORBSYM(I).NE.ORBSYM(I-1)) THEN
                  IF(ISYMNUM.GT.ISNMAX) ISNMAX=ISYMNUM
                  ISYMNUM=0
               ENDIF
            ENDIF
            ISYMNUM=ISYMNUM+1
            ARR(2*I-1,1)=Z
            ARR(2*I,1)=Z
         ELSEIF(I.NE.0.AND.K.EQ.0.AND.J.EQ.0) THEN
            ARR(2*I-1,1)=Z
            ARR(2*I,1)=Z
         ENDIF
C.. At the moment we're ignoring the core energy
         IF(I.NE.0) GOTO 1
99       CONTINUE
         DO I=1,NORB
            BRR(2*I-1)=2*I-1
            BRR(2*I)=2*I
C            G1(1,2*I-1)=ISYMNUM-1
C            G1(1,2*I)=ISYMNUM-1
            G1(2*I-1)%Sym%s=IBSET(0,ORBSYM(I)-1)
            G1(2*I)%Sym%s=IBSET(0,ORBSYM(I)-1)
C            G1(3,2*I-1)=I
C            G1(3,2*I)=I
C.. set momentum to 0
            G1(2*I-1)%k(1)=0
            G1(2*I)%k(2)=0
            G1(2*I-1)%k(3)=0
            G1(2*I)%k(1)=0
            G1(2*I-1)%k(2)=0
            G1(2*I)%k(3)=0
            G1(2*I-1)%Ms=-1
            G1(2*I)%Ms=1
         ENDDO
         SYMMAX=ORBSYM(NORB)
         IF(SYMMAX.EQ.1) THEN
         ELSEIF(SYMMAX.LE.2) THEN
         ELSEIF(SYMMAX.LE.4) THEN
            SYMMAX=4
         ELSEIF(SYMMAX.LE.8) THEN
            SYMMAX=8
         ENDIF
        WRITE(6,"(A,I)") "Maximum number of symmetries: ",SYMMAX
         NBASISMAX(1,1)=0
         NBASISMAX(1,2)=0
         NBASISMAX(5,1)=0
         NBASISMAX(5,2)=SYMMAX-1
         NBASISMAX(2,1)=0
         NBASISMAX(2,2)=0
         CLOSE(8)
         RETURN
      END
      SUBROUTINE READFCIINT(UMAT,NBASIS,ECORE,ARR,BRR,G1)
         USE HElement
         USE UMatCache , only : UMatInd,TMatind,TMat2D,TMATSYM,
     &                          TSTARBIN
         IMPLICIT NONE
         INTEGER NBASIS
         REAL*8 ECORE,ARR(NBASIS,2)
         REAL*8 UMAT(*)
         REAL*8 Z
         INTEGER I,J,K,L,BRR(NBASIS),G1(*)
         INTEGER NORB,NELEC,MS2,ORBSYM(1000),ISYM
         LOGICAL LWRITE
         include 'sym.inc'
         NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM
         LWRITE=.FALSE.
         OPEN(8,FILE='FCIDUMP',STATUS='OLD')
         READ(8,FCI)
101      READ(8,'(1X,G20.12,4I3)',END=199) Z,I,J,K,L
C.. Each orbital in the file corresponds to alpha and beta spinorbitalsa
         IF(I.EQ.0) THEN
C.. Core energy
            ECORE=Z
         ELSEIF(J.EQ.0) THEN
!C.. HF Eigenvalues
!            ARR(I*2-1,2)=Z
!            ARR(I*2,2)=Z
!            ARR(BRR(I*2-1),1)=Z
!            ARR(BRR(I*2),1)=Z
!            LWRITE=.TRUE.
         ELSEIF(K.EQ.0) THEN
C.. 1-e integrals
            IF(TSTARBIN) THEN
! If TSTARBIN, the one-el integrals are stored in symmetry classes, as spatial orbitals
                TMATSYM(TMatInd(2*J,2*I))=Z
            ELSE
C.. These are stored as spinorbitals (with elements between different spins being 0
                TMAT2D(2*I-1,2*J-1)=Z
                TMAT2D(2*I,2*J)=Z
                TMAT2D(2*J-1,2*I-1)=Z
                TMAT2D(2*J,2*I)=Z
            ENDIF
         ELSE
C.. 2-e integrals
C.. UMAT is stored as just spatial orbitals (not spinorbitals)
C..  we're reading in (IJ|KL), but we store <..|..> which is <IK|JL>
            UMAT(UMatInd(I,K,J,L))=Z
!i            UMAT(I,K,J,L)=Z
!            UMAT(J,L,I,K)=Z
!            UMAT(K,I,L,J)=Z
!            UMAT(L,J,K,I)=Z
!            UMAT(J,K,I,L)=Z
!            UMAT(I,L,J,K)=Z
!            UMAT(L,I,K,J)=Z
!            UMAT(K,J,L,I)=Z
         ENDIF
         IF(I.NE.0) GOTO 101
199      CONTINUE
         CLOSE(8)
C.. If we've changed the eigenvalues, we write out the basis again
!         IF(LWRITE) THEN
!            WRITE(6,*) "1-electron energies have been read in."
!            CALL WRITEBASIS(6,G1,NBASIS,ARR,BRR)
!         ENDIF
         RETURN
      END

