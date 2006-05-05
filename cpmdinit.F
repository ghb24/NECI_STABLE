      SUBROUTINE CPMDBASISINIT(NBASISMAX,ARR,BRR,G1,LEN)
         IMPLICIT NONE
         INCLUDE 'cpmddata.inc'
         INCLUDE 'sym.inc'
         INTEGER NBASISMAX(5,2), BRR(LEN),LEN,G1(5,LEN)
         REAL*8 ARR(LEN)
         INTEGER I,IDECOMP
C         NBASISMAX(1,1)=1
C         NBASISMAX(1,2)=NSTATES
C         NBASISMAX(1,3)=2
C         NBASISMAX(4,1)=-1
C         NBASISMAX(4,2)=1
         NBASISMAX(1,1)=0
         NBASISMAX(1,2)=0
         NBASISMAX(1,3)=2
         NBASISMAX(4,1)=-1
         NBASISMAX(4,2)=1

C.. set ISPINSKIP=0, to tell the SCRs that there's no UMAT
         NBASISMAX(2,3)=0

         CALL GENIRREPS
         CALL GENSYMTABLE
         CALL GENSYMLABELPAIRS

         CALL FAKEEIGENVALUES

         CALL IAZZERO(G1,LEN*5)
         DO I=1,NSTATES
           CALL DECOMPOSEREP(SYMLABELCHARS(1,SYMLABELS(1,I)),IDECOMP)
            G1(5,I*2-1)=IDECOMP
            G1(5,I*2)=IDECOMP
            G1(4,I*2-1)=-1
            G1(4,I*2)=1
C.. Fake the eigenvalues
            ARR(2*I-1)=EIGENVALUES(I)
            ARR(2*I)=EIGENVALUES(I)
            BRR(2*I-1)=2*I-1
            BRR(2*I)=2*I
         ENDDO
C.. show it's a generic spatial basis
         NBASISMAX(3,3)=1
      END
C.. Called to retrieve number of basis functions from CPMD
      SUBROUTINE CPMDSYSTEMINIT(LEN)
         IMPLICIT NONE
         INCLUDE 'cpmddata.inc'
         INTEGER LEN
         CALL FAKEREADCHARTABLE(NSTATES)
         LEN=NSTATES*2
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
      SUBROUTINE CPMDINIT2INDINT(TMAT,NHG)
         IMPLICIT NONE
         INTEGER NHG
         INCLUDE 'cpmddata.inc'
         REAL*8 TMAT(NHG,NHG)
         REAL*8 CPMD1EINT(NSTATES,NSTATES)
         INTEGER I,J,K
         REAL*8 TOT,F,U
           
         CALL CPMDGEN1EINTS(CPMD1EINT)
         DO I=1,NSTATES
            DO J=I,NSTATES
               TOT=0.D0
               IF(I.EQ.J) TOT=TOT+EIGENVALUES(I)+XI
               TOT=TOT-CPMD1EINT(I,J)
               DO K=1,NSTATES
                  CALL CPMDGETOCC(K,F)
C.. F now has the occupation number of state K
                  IF(F.GT.0.D0) THEN
C.. Get 4-index integral (ik|jk)
                     CALL INITFINDXI(I,K,J,K,U)
                     TOT=TOT-F*U
                  ENDIF
               ENDDO
               TMAT(2*I-1,2*J-1)=TOT
               TMAT(2*I,2*J)=TOT
               TMAT(2*J-1,2*I-1)=TOT
               TMAT(2*J,2*I)=TOT
            ENDDO
         ENDDO
         
      END   
