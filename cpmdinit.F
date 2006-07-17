      SUBROUTINE CPMDBASISINIT(NBASISMAX,ARR,BRR,G1,LEN)
         IMPLICIT NONE
         INCLUDE 'cpmddata.inc'
         INCLUDE 'sym.inc'
         INTEGER NBASISMAX(5,2), BRR(LEN),LEN,G1(5,LEN)
         REAL*8 ARR(LEN,2)
         INTEGER I,IDECOMP,J
         LOGICAL ALL1D
C         NBASISMAX(1,1)=1
C         NBASISMAX(1,2)=NSTATES
C         NBASISMAX(1,3)=2
C         NBASISMAX(4,1)=-1
C         NBASISMAX(4,2)=1
         NBASISMAX(1,1)=0
         NBASISMAX(2,1)=0
         NBASISMAX(3,1)=0
         IF(BIGCELL.AND.PROJECT) THEN
            NBASISMAX(1,2)=CELLSIZES(1)-1
            NBASISMAX(2,2)=CELLSIZES(2)-1
            NBASISMAX(3,2)=CELLSIZES(3)-1
         ELSE
            NBASISMAX(1,2)=0
            NBASISMAX(2,2)=0
            NBASISMAX(3,2)=0
         ENDIF
         NBASISMAX(1,3)=2
         NBASISMAX(4,1)=-1
         NBASISMAX(4,2)=1

C.. set ISPINSKIP=0, to tell the SCRs that there's no UMAT
         NBASISMAX(2,3)=0
         CALL EXTRACTSYMLABELS
         CALL GENIRREPS(IMPROPER_OP)
         CALL GENSYMTABLE

C         CALL FAKEEIGENVALUES

         CALL IAZZERO(G1,LEN*5)
         ALL1D=.TRUE.
         DO I=1,NSTATES
           CALL DECOMPOSEREP(SYMLABELCHARS(1,SYMLABELS(1,I)),IDECOMP)
            IF(SYMLABELCHARS(1,SYMLABELS(1,I)).NE.1) ALL1D=.FALSE.
            G1(5,I*2-1)=IDECOMP
            G1(5,I*2)=IDECOMP
            G1(4,I*2-1)=-1
            G1(4,I*2)=1
            IF(BIGCELL.AND.PROJECT) THEN
               DO J=1,3
                  G1(J,I*2-1)=K_VECTORS(J,TRANS_LABEL(I))
                  G1(J,I*2)=K_VECTORS(J,TRANS_LABEL(I))
               ENDDO
            ELSE
               DO J=1,3
                  G1(J,I*2-1)=0
                  G1(J,I*2)=0
               ENDDO
            ENDIF
C.. Fake the eigenvalues
            ARR(2*I-1,1)=EIGENVALUES(I)+XI
            ARR(2*I,1)=EIGENVALUES(I)+XI
            ARR(2*I-1,2)=EIGENVALUES(I)+XI
            ARR(2*I,2)=EIGENVALUES(I)+XI
            BRR(2*I-1)=2*I-1
            BRR(2*I)=2*I
         ENDDO
         IF(ALL1D) THEN
C.. If all reps are 1D, we can block determinants in different symmetries
            NBASISMAX(5,1)=0
            NBASISMAX(5,2)=NSYM-1
            WRITE(6,*)
     &   "All orbitals have 1D symmetry. Using blocking symmetry."
         ELSE
            WRITE(6,*)
     &   "Multidimensional symmetries detected. ",
     &   "Not using symmetry blocking."
         ENDIF
C.. show it's a generic spatial basis
         NBASISMAX(3,3)=1
      END
C.. Called to retrieve number of basis functions from CPMD
      SUBROUTINE CPMDSYSTEMINIT(LEN)
         IMPLICIT NONE
         INCLUDE 'cpmddata.inc'
         INCLUDE 'sym.inc'
         INTEGER LEN
C         CALL FAKEREADCHARTABLE(NSTATES)
         LEN=NSTATES*2
         NROT=NROTOP
      END
      
      SUBROUTINE GENCPMDSYMREPS(G1,NBASIS,ARR,BRR)
         IMPLICIT NONE
         INCLUDE 'cpmddata.inc'
         INTEGER BRR(NBASIS),NBASIS,G1(5,NBASIS)
         REAL*8 ARR(NBASIS,2)
         CALL GENSYMREPS(G1,NBASIS,ARR,BRR,DEGENTOL)
         CALL GENSYMLABELPAIRS(NBASIS/2)
      END

      SUBROUTINE EXTRACTSYMLABELS
         IMPLICIT NONE
         INCLUDE 'cpmddata.inc'
         INCLUDE 'sym.inc'
         INTEGER I,J
         NSYMLABELS=0
         DO I=1,NSTATES
            IF(ROT_LABEL(I).GT.NSYMLABELS) NSYMLABELS=ROT_LABEL(I)
         ENDDO
         CALL MEMORY(IP_SYMLABELCHARS,NSYMLABELS*NROT,'SYMLABELCH')
         CALL IAZZERO(SYMLABELCHARS,NSYMLABELS*NROT)
         CALL MEMORY(IP_SYMLABELS,NSTATES*2,'SYMLABELS')
         DO I=1,NSTATES
            IF(SYMLABELCHARS(1,ROT_LABEL(I)).EQ.0) THEN
C.. we've found an uninitialized symlabel - fill it
               DO J=1,NROT
C.. Indices need to be inverted so cannot straight ICOPY
                  SYMLABELCHARS(J,ROT_LABEL(I))=NINT(ROT_CHAR(I,J))
               ENDDO
            ENDIF
            SYMLABELS(1,I)=ROT_LABEL(I)
         ENDDO
         WRITE(6,*) "SYMMETRY CLASSES"
         CALL WRITEIRREPTAB(6, SYMLABELCHARS,NROT,NSYMLABELS)
C.. Allocate memory gor irreps.
C.. Assume there will be no more than 16 irreps
         CALL MEMORY(IP_IRREPCHARS,NROT*16,"IRREPCH")
      END

      SUBROUTINE FAKEEIGENVALUES
         IMPLICIT NONE
         INCLUDE 'cpmddata.inc'
         INTEGER I
         CALL MEMORY(IP_EIGENVALUES,NSTATES,"EIGENVALUES")
         DO I=1,NSTATES
            EIGENVALUES(I)=I
         ENDDO
      END 
      SUBROUTINE FAKEREADCHARTABLE(NSTATES)
         IMPLICIT NONE
         INCLUDE 'sym.inc'
C.. First set up the labels and chars
         INTEGER CHARS(16,7)
         DATA ((CHARS(J,I),J=1,16),I=1,7)
     &               /1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     &               2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2,
     &               2, 0, 0,-2, 0, 0, 0, 0,-2, 0, 0, 2, 0, 0, 0, 0,
     &               1,-1,-1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1, 1,-1,-1,
     &               4, 0, 0,-4, 0, 0, 0, 0,-4, 0, 0, 4, 0, 0, 0, 0,
     &               1,-1,-1, 1,-1,-1, 1, 1, 1,-1,-1, 1,-1,-1, 1, 1,
     &               2,-2,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,-2,-2/
         INTEGER LABS(19)
         DATA LABS /1,2,2,1,3,3,4,2,2,1,5,5,5,5,
     &              3,3,6,7,7/
         INTEGER NSTATES,I,J
         INTEGER ARRAY(*)
         POINTER (IP_ARRAY,ARRAY)
         NSYMLABELS=7
         NROT=16
         NSTATES=19
         CALL MEMORY(IP_SYMLABELS,NSTATES*2,"SYMLABELS")
         CALL MEMORY(IP_SYMLABELCHARS,NROT*NSYMLABELS,"SYMLABELCH")
         IP_ARRAY=IP_SYMLABELCHARS
         CALL MEMORY(IP_IRREPCHARS,NROT*NSYMLABELS,"IRREPCH")
         I=IP_SYMLABELS
         DO I=1,7
            DO J=1,16
               ARRAY((I-1)*16+J)=CHARS(J,I)
            ENDDO
         ENDDO
         DO I=1,NSTATES
            SYMLABELS(1,I)=LABS(I)
         ENDDO
         WRITE(6,*) "SYMMETRY CLASSES", IP_SYMLABELCHARS,
     &      LOC(SYMLABELCHARS)
         CALL WRITEIRREPTAB(6, SYMLABELCHARS,NROT,NSYMLABELS)
      END

C.. Calculate the h_ij elements for CPMD orbitals.
C.. E_i is the eigenvalue for state i.
C.. f_i is the occupation number for state i.
C.. delta_ij is the Dirac delta.
C.. (ab|cd) is the 4-index Coulomb integral between states, i,j,k,l
C.. h_ij=(E_i+xi) delta_ij - Sum_k f_k (ik|jk)
      SUBROUTINE CPMDINIT2INDINT(TMAT,NHG,NORBUSED,NBASISMAX,ISS,G1,NEL,
     &      ECORE)
         IMPLICIT NONE
         INTEGER NHG,NTFROZEN
         INCLUDE 'cpmddata.inc'
         REAL*8 TMAT(NHG,NHG)
         REAL*8 CPMD1EINT(NSTATES,NSTATES)
         INTEGER I,J,K,NSTATESUSED,ISUB,NORBUSED,II,JJ
         REAL*8 TOT,F,U 
         REAL*8 GETUMATEL
         INTEGER ISS(*),G1(5,*),NBASISMAX(*)
         REAL*8 ECORE
         INTEGER NEL,SYMPROD
         LOGICAL LCHKSYM
         CALL TISET('CPMDINIT2I',ISUB)

         ECORE=EIONION-NEL*NEL*XI/2.0
         WRITE(6,*) "Core Energy:", ECORE 
         NSTATESUSED=NORBUSED/2
         CALL SETUPUMATCACHE(NSTATESUSED)
         CALL CPMDGEN1EINTS(CPMD1EINT)
         OPEN(10,FILE='TMAT',STATUS='UNKNOWN')
C.. Set the UMAT cache to cache things in the lowest possible position
         CALL SETUMATCACHEFLAG(1)
C   JSS  Calculate the elements of the type <ij|ij> and <ij|ji> efficiently.
         CALL CPMDANTISYMINTEL(G1)
         DO I=1,NSTATESUSED
           DO J=I,NSTATESUSED
             II=I*2-1
             JJ=J*2-1
             IF(LCHKSYM(G1(1,II),G1(1,JJ))) THEN
              TOT=0.D0
                 IF(I.EQ.J) TOT=TOT+EIGENVALUES(I)
                 TOT=TOT-CPMD1EINT(I,J)
                   DO K=1,NSTATES
                      CALL CPMDGETOCC(K,F)
C.. F now has the occupation number of state K
                      IF(F.GT.1.D-5) THEN
C.. Get 4-index integral (ik|jk)
                         IF(K.GT.NSTATESUSED) THEN
                            WRITE(6,"(A,I2,A,F10.6)")"TOP FROZEN STATE "
     &                          ,K," has occupation number ",F
                            STOP "TOP FROZEN STATE OCCUPIED"
                         ENDIF
C                         WRITE(6,*) "U",I,K,J,K
                         U=GETUMATEL(NBASISMAX,0.D0,0.D0,NHG,ISS,G1,
     &          I,K,J,K)
C                         CALL INITFINDXI(I,K,J,K,U)
C                         WRITE(6,*) F,U
                         TOT=TOT-F*U
                      ENDIF
                   ENDDO
                 TMAT(2*I-1,2*J-1)=TOT
                 TMAT(2*I,2*J)=TOT
                 TMAT(2*J-1,2*I-1)=TOT
                 TMAT(2*J,2*I)=TOT
                 WRITE(10,*) I,J,TOT
              ENDIF
           END DO
         END DO
C.. Set the UMAT cache to cache normally again
         CALL SETUMATCACHEFLAG(0)
         CLOSE(10)
         CALL TIHALT('CPMDINIT2I',ISUB)
      END   
