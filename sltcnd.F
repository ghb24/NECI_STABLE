!.. Use the Slater-Condon Rules to evaluate the H matrix element between two dets
      SUBROUTINE SLTCND(NEL,NBASISMAX,
     &	NHG,NDET1,NDET2,GBS,IC,NMSH,FCK,
     &	NMAX,ALAT,UMat,SUM)
      USE HElem
      use SystemData , only : FCOUL,TEXCH
!      use global_utilities
      IMPLICIT NONE
C..
C.. 
#ifdef POINTER8
      INTEGER*8 LOC
#else
      INTEGER LOC
#endif
      INTEGER NDET1(NEL),NDET2(NEL)
      INTEGER N1E(NHG),N2E(NHG)
      INTEGER GBS(4,*)
      TYPE(HElement) UMAT(*),SUM
      REAL*8 KP
      INTEGER nBasisMax(5,*)
C..Cube Arrays
      COMPLEX*16 FCK(-NMSH/2:NMSH/2-1,-NMSH/2:NMSH/2-1,-NMSH/2:NMSH/2-1)
      REAL*8 ALAT(3)
      INTEGER ISPINSKIP,NEL,NHG,NMSH
      INTEGER IC,NMAX
      INTEGER IFLAG,IPTOT
      REAL*8 COULCOUPLE
!      type(timer), save :: proc_timer
!      proc_timer%timer_name='SLTCND    '
!      call set_timer(proc_timer)
      ISPINSKIP=NBASISMAX(2,3)
C..
C..Routine to read in UMAT
      SUM=0.D0
C..We now check to see which Slater-Condon rule to apply
      COULCOUPLE=FCOUL
      IF(IC.EQ.(NEL-2)) THEN
C..Need to get the determinants in the form of max. coincidence
!        CALL IAZZERO(N1E,NHG)
!        CALL IAZZERO(N2E,NHG)
C..Also returns the sign change IPTOT
!        CALL LINEUP(NEL,NDET1,NDET2,NHG,N1E,N2E,IFLAG,IPTOT)
C..Differ by two 
        CALL SCR2(NEL,NBASISMAX,NDET1,NDET2,GBS,
     &         NHG,UMAT,ALAT,ISPINSKIP,COULCOUPLE,TEXCH,SUM)
      ELSEIF(IC.EQ.(NEL-1)) THEN
C..Need to get the determinants in the form of max. coincidence
!        CALL IAZZERO(N1E,NHG)
!        CALL IAZZERO(N2E,NHG)
C..Also returns the sign change IPTOT
!        CALL LINEUP(NEL,NDET1,NDET2,NHG,N1E,N2E,IFLAG,IPTOT)
C..Differ by one 
        CALL SCR1(NEL,NBASISMAX,NDET1,NDET2,GBS,
     &         NHG,UMAT,ALAT,ISPINSKIP,COULCOUPLE,TEXCH,SUM)
      ELSEIF(IC.EQ.NEL) THEN
C..Differ by none 
!        IPTOT=1
        CALL SCR0(NEL,NBASISMAX,NDET1,NDET2,GBS,
     &         NHG,UMAT,ALAT,ISPINSKIP,COULCOUPLE,TEXCH,SUM)
C.. If we're in the UEG, we need to add in the part of the background
C.. contribution which doesn't cancel
         IF(NBASISMAX(1,3).EQ.-1) THEN
C            WRITE(6,*) FCK(0,0,0)
CC removed by AJWT 6/3/06
C            SUM=SUM-NEL*REAL(FCK(0,0,0))/2

         ENDIF
      ELSE
        SUM=0.D0
        RETURN
      ENDIF
!      SUM=SUM*DFLOAT(IPTOT)
!      call halt_timer(proc_timer)
      RETURN
      END

      REAL*8 FUNCTION COULDAMP(I,NHG)
         use SystemData , only : FCOULDAMPBETA,FCOULDAMPMU
         IMPLICIT NONE
         INTEGER I,NHG,J
         REAL*8 DERFC
         IF(FCOULDAMPBETA.LT.0.D0) THEN
           COULDAMP=1.D0
         ELSE
            STOP 'COULDAMP non-functional'
!            COULDAMP=DERFC(FCOULDAMPBETA*(ARRSTORE(NHG+I)-FCOULDAMPMU))
!&                  /2.D0
         ENDIF
         RETURN
      END
