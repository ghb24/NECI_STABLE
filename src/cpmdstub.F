! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
C.. A stub containing dummy functions in CPMD, which will appear on linking appropriately.

      SUBROUTINE INITFINDXI(I,J,K,L,UMATEL,gen2ints)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER I,J,K,L
         real(dp) UMATEL
         logical gen2ints

         ! Make warnings go away
         i=i; j=j; k=k; l=l; gen2ints=gen2ints
         UMATEL=0.0_dp
      END

      SUBROUTINE CPMDANTISYMINTEL(G1,UMat2D,HarInt,NStates) 
         use constants, only: dp
         use SystemData, only: BasisFN
         IMPLICIT NONE
         type(BasisFN) :: G1(*)
         complex(dp) :: HarInt(*)
         integer :: nStates
         HElement_t :: UMat2D(*)
         !eliminate warnings
         G1(1)=G1(1)
         nStates=nStates
         UMat2D(1)=UMat2D(1)
         HarInt(1)=HarInt(1)
         call stop_all("CPMDANTISYMINTEL","Entering stub routine.")
      END SUBROUTINE

      SUBROUTINE CPMDGETOCC(ISTATE,OCC)
         use constants, only: dp
         IMPLICIT NONE
C         INCLUDE 'elct.inc'
         real(dp) OCC
         INTEGER ISTATE
C         COMMON /ELCT/ F
C         OCC=F(ISTATE)
         OCC=0.0_dp
         istate=istate
      END
      
      SUBROUTINE CPMDGEN1EINTS(INTS,XIS)
         use constants, only: dp
         use CPMDData, only: nStates
         IMPLICIT NONE
         real(dp) INTS(NSTATES,NSTATES)
         real(dp) XIS
         INTS=0.0_dp
         xis=xis
      END

      subroutine GetKPInd(IndA,IndB,IndC,IndD)
         implicit none
         integer  IndA,IndB,IndC,IndD

         ! Eliminate warnings
         IndA=IndA; IndB=IndB; IndC=IndC; IndD=IndD
      end subroutine GetKPInd

      subroutine CalcHarPIInts(HarInt,nStates)
         use constants, only: dp
         implicit none
         integer  nStates
         complex(dp) HarInt(nStates,nStates)

         ! Eliminate Warnings
         HarInt=HarInt; nStates=nStates
      end subroutine CalcHarPIInts



      subroutine KPntSymInt(i,j,k,l,a,b,c,d)
          implicit none
          integer i,j,k,l,a,b,c,d
          ! Make warnings go away
          i=i; j=j; k=k; l=l; a=a; b=b; c=c; d=d
      end subroutine KPntSymInt
