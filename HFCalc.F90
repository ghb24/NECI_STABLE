#include "macros.h"
MODULE HFCalc
   use HElem
   implicit none
   save
   contains
      subroutine HFDoCalc()
      Use global_utilities
      use SystemData, only: tStarStore, BasisFN
      use IntegralsData, only: tHFBasis, tHFCalc, iHFMethod, tReadHF, nHFIt, HFMix, HFCDelta, HFEDelta
      use IntegralsData, only: HFRand, tRHF, ntFrozen, tReadTUMat
      use SystemData, only : tCPMD,  tHFOrder,nBasisMax, G1, Arr, Brr, ECore, nEl, nBasis, iSpinSkip, LMS
      use SystemData, only : tHub, lmsbasis
      Use Logging, only: iLogging
      Use Determinants, only: FDet, nUHFDet
      use IntegralsData, only: UMat, tagUMat
      Use UMatCache, only: GetUMatSize
      Use OneEInts, only: TMat2D, SetupTMat2, DestroyTMat
      character(25), parameter :: this_routine='HFDoCalc'
      REAL*8 HFBASIS(*),HFE(*)
      Type(HElement),pointer :: UMat2(:)
      Type(HElement),pointer :: TMat2D2(:,:)
      POINTER (IP_HFBASIS,HFBASIS),(IP_HFE,HFE)
      integer i, ierr
      integer nOrbUsed
      integer UMatInt, TMatInt
      integer tagUMat2, tagTMat2
         
!C.. If we are using an HF basis instead of our primitive basis, we need
!C.. to load in the coeffs of the HF eigenfunctions in terms of the
!C.. primitive basis.
!C.. We load the coeffs from a file HFBASIS
         IF(THFBASIS.OR.THFCALC.OR.(THFORDER.AND..NOT.TCPMD)) THEN
            CALL MEMORY(IP_HFBASIS,nBasis*nBasis*HElementSize,'HFBASIS')
            CALL AZZERO(HFBASIS,nBasis*nBasis*HElementSize)
!C.. Allocate an array to store the HF Energies
            CALL MEMORY(IP_HFE,nBasis,'HFE')
            CALL AZZERO(HFE,nBasis)
            IF(THFORDER.AND..NOT.THFBASIS) THEN
!C.. If we're not using HF, but just calculating the HF order
!C.. We generate the HF energies (this has no mixing or randomisation, so should jsut
!C.. re-order the orbitals and give us some energy)
!C.. HF basis is NOT using the LMS value set in the input
              CALL CALCHFBASIS(NBASIS,ISPINSKIP,NBASISMAX,G1,ARR,BRR,ECORE,UMAT,HFE,HFBASIS,1,NEL,LMSBASIS,1.D0,HFEDELTA,HFCDELTA,.TRUE.,0,TREADHF,0.D0,FDET,ILOGGING)
               CALL ORDERBASISHF(ARR,BRR,HFE,HFBASIS,G1,NBASIS,FDET,NEL)
            ELSEIF(THFCALC) THEN
              CALL CALCHFBASIS(NBASIS,ISPINSKIP,NBASISMAX,G1,ARR,BRR,ECORE,UMAT,HFE,HFBASIS,NHFIT,NEL,LMS,HFMIX,HFEDELTA,HFCDELTA,TRHF,IHFMETHOD,TREADHF,HFRAND,FDET,ILOGGING)
               CALL SETUPHFBASIS(NBASISMAX,G1,NBASIS,HFE,ARR,BRR)
            ELSEIF(THFBASIS) THEN
               CALL READHFBASIS(HFBASIS,HFE,G1,NBASIS)
               CALL SETUPHFBASIS(NBASISMAX,G1,NBASIS,HFE,ARR,BRR)
            ENDIF
            WRITE(6,*) "FINAL HF BASIS"
            CALL WRITEBASIS(6,G1,nBasis,ARR,BRR)

            WRITE(6,"(A,$)") "Fermi det (D0):"
            CALL WRITEDET(6,FDET,NEL,.TRUE.)
            CALL FLUSH(6)
!C.. If in Hubbard, we generate site-spin occupations
            IF(THUB) THEN
!  Don't think this works
!               CALL GENSITESPINOCC(NBASIS,NBASIS/ISPINSKIP,ISPINSKIP,NBASISMAX,G1,NEL,LMS,BRR,HFBASIS)
            ENDIF
!C.. We now generate a new U matrix corresponding to the HF basis fns
!C.. This requires a new matrix UMAT2 to be set up, as our HF basis is
!C.. unrestricted
            !THIS ROUTINE NO LONGER WORKS WITH NEW TMAT/UMAT MODULARISATION
            IF(THFBASIS) THEN
               WRITE(6,*) "Allocating TMAT2"
!               CALL MEMORY(IP_TMAT2,nBasis**2,'TMAT2')
!               CALL AZZERO(TMAT2,nBasis**2)
               IF(TSTARSTORE) STOP 'TSTARSTORE WITH HFBASIS?!'
               CALL SetupTMAT2(nBasis,2,TMATINT)
               NORBUSED=NBASIS-NTFROZEN
               IF(TREADTUMAT) THEN
                  CALL READHFTMAT(NBASIS,HFBASIS)
               ELSE
                  CALL CALCHFTMAT(NBASIS,HFBASIS,NORBUSED)
               ENDIF
                CALL DestroyTMAT(.false.)
                TMAT2D => TMAT2D2
                NULLIFY(TMAT2D2)
!               CALL FREEM(IP_TMAT)
!               IP_TMAT=IP_TMAT2
!               IP_TMAT2=NULL
!C.. Allocate the new matrix
               CALL GetUMatSize(nBasis,nEl,1,UMATINT)
               Allocate(UMat2(UMatInt), stat=ierr)
               LogAlloc(ierr,'UMAT2', UMatInt, HElementSizeB, tagUMat2)
               CALL AZZERO(UMAT2,HElementSize*UMATINT)
!C.. We need to pass the TMAT to CALCHFUMAT as TMAT is no longer diagona
!C.. This also modified G1, ARR, BRR
               IF(TREADTUMAT) THEN
                 CALL READHFUMAT(UMAT,UMAT2,NBASIS,NBASISMAX,G1,HFBASIS,ISPINSKIP,HFE,ARR,BRR)
               ELSE
                 CALL CALCHFUMAT(UMAT,UMAT2,NBASIS,NBASISMAX,G1,HFBASIS,ISPINSKIP,HFE,ARR,BRR,NORBUSED)
               ENDIF
!C.. Now we can remove the old UMATRIX, and set the pointer UMAT to point
!C.. to UMAT2
               LogDealloc(tagUMat)
               Deallocate(UMat)
               UMat=>UMat2
               nullify(UMat2)
               tagUMat=tagUMat2
               tagUMat2=0
!C.. spinskip values
               ISPINSKIP=1
               NBASISMAX(2,3)=1
!C.. Indicate that we're in UHF and that <D0|H|D1>=0 for D1 being a
!C.. single excitation
               IF(LMS.EQ.0) THEN
                  NBASISMAX(4,5)=1
                  DO I=1,NEL
                     NUHFDET(I)=BRR(I)
                  ENDDO
                  CALL SORTI(NEL,NUHFDET)
               ELSE
                  NBASISMAX(4,5)=2
               ENDIF
            ENDIF
!C.. Now deallocate the HF arrays
            CALL FREEM(IP_HFBASIS)
            CALL FREEM(IP_HFE)
            CALL FLUSH(6)
         ENDIF
      End Subroutine HFDoCalc
End Module HFCalc
