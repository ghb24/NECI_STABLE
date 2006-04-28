      SUBROUTINE CPMDBASISINIT(NBASISMAX,ARR,BRR,G1,LEN)
         IMPLICIT NONE
         INCLUDE 'cpmddata.inc'
         INCLUDE 'sym.inc'
         INTEGER NBASISMAX(5,2), BRR(LEN),LEN,G1(5,LEN)
         REAL*8 ARR(LEN)
         INTEGER I,IDECOMP,CNORM
         LOGICAL L
         LOGICAL GETIRREPDECOMP
         NBASISMAX(1,1)=1
         NBASISMAX(1,2)=NSTATES
         NBASISMAX(1,3)=2
         NBASISMAX(4,1)=-1
         NBASISMAX(4,2)=1
         CALL GENIRREPS
         CALL GENSYMTABLE
         CALL GENSYMLABELPAIRS
         CALL IAZZERO(G1,LEN*5)
         DO I=1,NSTATES
            L=GETIRREPDECOMP(SYMLABELCHARS(1,SYMLABELS(I)),IRREPCHARS,
     &         NSYM,NROT,IDECOMP,CNORM)
            G1(5,I*2-1)=IDECOMP
            G1(5,I)=IDECOMP
            G1(4,I*2-1)=-1
            G1(4,I*2)=1
C.. Fake the eigenvalues
            ARR(I)=I*1.D0
            BRR(I)=I
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
         WRITE(6,*) I
         DO I=1,7
            DO J=1,16
               ARRAY((I-1)*16+J)=CHARS(J,I)
            ENDDO
         ENDDO
         DO I=1,NSTATES
            SYMLABELS(1,I)=LABS(I)
         ENDDO
         WRITE(6,*) "SYMMETRY CLASSES"
         CALL WRITEIRREPTAB(6, SYMLABELCHARS,NROT,NSYMLABELS)
      END

