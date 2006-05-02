C.. A stub containing dummy functions in CPMD, which will appear on linking appropriately.

      SUBROUTINE INITFINDXI(I,J,K,L,UMATEL)
         IMPLICIT NONE
         INTEGER I,J,K,L
         REAL*8 UMATEL
         UMATEL=0.D0
      END

      SUBROUTINE CPMDGETOCC(ISTATE,OCC)
         IMPLICIT NONE
C         INCLUDE 'elct.inc'
         REAL*8 F(1),OCC
         INTEGER ISTATE
C         COMMON /ELCT/ F
C         OCC=F(ISTATE)
         OCC=0.D0
      END
      
      SUBROUTINE CPMDGEN1EINTS(INTS)
         IMPLICIT NONE
         REAL*8 INTS(NSTATES,NSTATES)
         INCLUDE 'cpmddata.inc'
         CALL AZZERO(INTS,NSTATES*NSTATES)
      END
