!A simple module with some interfaces to avoid compilation warnings. Ignore.
module neci_intfce
 interface
  SUBROUTINE GENSYMEXCITIT2(NI,NEL,G1,NBASIS,TSETUP,NMEM,NJ,IC,STORE,ILEVEL)
   use SystemData, only: BasisFN
   IMPLICIT NONE
   INTEGER NEL,NI(NEL),NBASIS
   TYPE(BasisFN) G1(nBasis)
   INTEGER STORE(6)
   INTEGER, target :: NMEM(*)
   INTEGER NJ(NEL),IC
   LOGICAL TSETUP
   INTEGER ILEVEL
  END
  SUBROUTINE GENRANDSYMEXCITIT2(NI,NEL,NMEM,NJ,ISEED,ICOUNT,PGEN)
   use constants, only: dp
   IMPLICIT NONE
   INTEGER NEL,NI(NEL)
   INTEGER,target :: NMEM(*)
   integer NJ(NEL),ICOUNT
   INTEGER ISEED
   real(dp) PGEN
  END
  SUBROUTINE GENSYMEXCITIT3Par(NI,TSETUP,NMEM,NJ,IC,STORE,ILEVEL,iMinElec1,iMaxElec1)
   use SystemData, only: nEl
   IMPLICIT NONE
   INTEGER NI(NEL)
   INTEGER, pointer :: DSTORE(:)
   INTEGER STORE(6)
   INTEGER, target ::  NMEM(*)
   INTEGER NJ(NEL),IC
   LOGICAL TSETUP
   INTEGER ILEVEL
   INTEGER iMinElec1, iMaxElec1
  END
  SUBROUTINE GenExcitProb(nI,nJ,nEl,nIExcitor,G1,nBasisMax,Arr,nBasis,pGen)
   use constants, only: dp
   use SystemData, only: BasisFN
   IMPLICIT NONE
   INTEGER nEl,nI(nEl),nJ(nEl),nBasis,nBasisMax(*)
   INTEGER, target :: nIExcitor(*)
   TYPE(BasisFN) G1(nBasis)
   real(dp) pGen
   real(dp) Arr(nBasis,2)
  END
 end interface
end module neci_intfce
