C.. A stub containing dummy functions in CPMD, which will appear on linking appropriately.

      SUBROUTINE INITFINDXI(I,J,K,L,UMATEL,gen2ints)
         IMPLICIT NONE
         INTEGER I,J,K,L
         REAL*8 UMATEL
         logical gen2ints
         UMATEL=0.D0
      END

      SUBROUTINE CPMDANTISYMINTEL 
         IMPLICIT NONE
      END SUBROUTINE

      SUBROUTINE CPMDGETOCC(ISTATE,OCC)
         IMPLICIT NONE
C         INCLUDE 'elct.inc'
         REAL*8 F(1),OCC
         INTEGER ISTATE
C         COMMON /ELCT/ F
C         OCC=F(ISTATE)
         OCC=0.D0
      END
      
      SUBROUTINE CPMDGEN1EINTS(INTS,XIS)
         use CPMDData, only: nStates
         IMPLICIT NONE
         REAL*8 INTS(NSTATES,NSTATES)
         REAL*8 XIS
         CALL AZZERO(INTS,NSTATES*NSTATES)
      END

      subroutine GetKPInd(IndA,IndB,IndC,IndD)
         implicit none
         integer  IndA,IndB,IndC,IndD
         return
      end subroutine GetKPInd

      subroutine CalcHarPIInts(HarInt,nStates)
         implicit none
         integer  nStates
         complex*16 HarInt(nStates,nStates)
         return
      end subroutine CalcHarPIInts



      subroutine KPntSymInt(i,j,k,l,a,b,c,d)
          implicit none
          integer i,j,k,l,a,b,c,d
          return
      end subroutine KPntSymInt
