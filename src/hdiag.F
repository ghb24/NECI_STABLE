#include "macros.h"

!Diagonalize a compressed hamiltonian (in HAMIL), returning eigenvalues in W and eigenvectors in CK.
       SUBROUTINE HDIAG_neci(NDET,HAMIL,LAB,NROW,CK,W,WORK2,WORK,
     &                     NBLOCKSTARTS,NBLOCKS)
         use constants, only: dp,int32!,n_int
         use SystemData, only: BasisFN
         use global_utilities
         use HElem, only: HElement_t_size
         use util_mod, only: stop_all
         IMPLICIT NONE
!         INTEGER(n_int) :: BitRep(0:NIfTot),BitRepSym(0:NIfTot)
!         INTEGER(n_int) :: TempBitRep(0:NIfTot),BitRepHF(0:NIfTot)
         INTEGER NDET
         INTEGER NROW(NDET),LAB(*)
         HElement_t(dp) HAMIL(*),CK(NDET,NDET),WORK2(3*NDET),WORK
!         HElement_t(dp) CopyCK(NDET,NDET)
         real(dp) W(NDET)
         INTEGER I,J,IND,INDZ,NBS!,nopen,k,Excit
         INTEGER(int32) INFO
         INTEGER NBLOCKS,NBLOCKSTARTS(NBLOCKS+1)
         type(timer), save :: proc_timer
         real(dp) GSEN
         character(len=*), parameter :: t_r="HDIAG_neci"
#ifndef CMPLX_
         unused_var(WORK)
#endif

         proc_timer%timer_name='HDIAG     '
         call set_timer(proc_timer)
         GSEN=1.D100
         DO I=1,NDET
            DO J=1,NDET
               CK(I,J)=0.0_dp
            ENDDO
         ENDDO
C.. Now we fill the RIJ array
         IND=1
         INDZ=1
         DO I=1,NDET
            INDZ=INDZ+NROW(I)
            DO WHILE (IND.LT.INDZ)
               CK(I,LAB(IND))=HAMIL(IND)
               IND=IND+1
            ENDDO
         ENDDO
         IF(HElement_t_size.ne.1) THEN
             !ensure hermitian
             do i=1,ndet
                 do j=1,i-1
#ifdef CMPLX_
                     !To avoid compiler warnings
                     CK(i,j)=CONJG(CMPLX(CK(j,i),kind=kind(0.0_dp)))
#endif
                 enddo
             enddo
         ENDIF

!****************************
C.. Diagonalize
         DO I=1,NBLOCKS
            NBS=NBLOCKSTARTS(I)
            IF(HElement_t_size.EQ.1) THEN
#ifndef CMPLX_
               CALL DSYEV('V','U',NBLOCKSTARTS(I+1)-NBS,CK(NBS,NBS),
     &            NDET,W(NBS),WORK2,3*NDET,INFO)
#else
               call stop_all(t_r,"This block for real calcs only")
#endif
            ELSE
#ifdef CMPLX_
                 CALL ZHEEV('V','U',NBLOCKSTARTS(I+1)-NBS,CK(NBS,NBS),
     &            NDET,W(NBS),WORK,4*NDET,WORK2,INFO)
#else
                 call stop_all(t_r,"This block for complex calcs only")
#endif
            ENDIF
            IF(INFO.NE.0) THEN
               WRITE(6,*) 'DYSEV error: ',INFO
               call stop_all(t_r, "DSYEV error")
            ENDIF
            IF(W(NBS).LT.GSEN) GSEN=W(NBS)
         ENDDO

C.. CK now contains the eigenvectors, and W the eigenvalues
         WRITE(6,"(A,F19.11,I4)") "GROUND E=",GSEN
         call halt_timer(proc_timer)
         RETURN

      END
