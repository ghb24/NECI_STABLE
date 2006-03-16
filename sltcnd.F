C.. In fact ZIA is not used by this routine at all, so we've instead
C.. hidden the T matrix in it.
C.. ZIA(1:NBASIS)
      SUBROUTINE SLTCND(NEL,NBASISMAX,
     &	NHG,NDET1,NDET2,GBS,IC,NMSH,FCK,ZIA,
     &	NMAX,ALAT,UMAT,SUM)
      IMPLICIT NONE
C..
C.. 
#ifdef POINTER8
      INTEGER*8 GET_ADDR
#else
      INTEGER GET_ADDR
#endif
      INTEGER NDET1(NEL),NDET2(NEL)
      INTEGER N1E(NHG),N2E(NHG)
      INTEGER GBS(4,*)
      REAL*8 UMAT(*)
      REAL*8 KP
      INTEGER NBASISMAX(4,3)
C..Cube Arrays
      COMPLEX*16 FCK(-NMSH/2:NMSH/2-1,-NMSH/2:NMSH/2-1,-NMSH/2:NMSH/2-1)
C      COMPLEX*16 ZIA(-NMSH/2:NMSH/2,NMAX,NMAX)
      REAL*8 ALAT(3),ZIA(*),SUM
      INTEGER ISPINSKIP,NEL,NHG,NMSH
      INTEGER IC,NMAX
      INTEGER IFLAG,IPTOT
      REAL*8 COULCOUPLE
      INCLUDE 'fcoul.inc'
      ISPINSKIP=NBASISMAX(2,3)
C..
C..Routine to read in UMAT
C..Need to get the determinants in the form of max. coincidence
      CALL IAZZERO(N1E,NHG)
      CALL IAZZERO(N2E,NHG)
C..Also returns the sign change IPTOT
      CALL LINEUP(NEL,NDET1,NDET2,NHG,N1E,N2E,IFLAG,IPTOT)
      SUM=0.D0
C..We now check to see which Slater-Condon rule to apply
      COULCOUPLE=FCOUL
      IF(IC.EQ.(NEL-2)) THEN
C..Differ by two 
        CALL SCR2(NEL,NBASISMAX,NDET1,NDET2,GBS,
     &         NHG,UMAT,ALAT,ISPINSKIP,COULCOUPLE,SUM)
      ELSEIF(IC.EQ.(NEL-1)) THEN
C..Differ by one 
        CALL SCR1(NEL,NBASISMAX,NDET1,NDET2,GBS,
     &         NHG,UMAT,ALAT,ISPINSKIP,ZIA,COULCOUPLE,SUM)
      ELSEIF(IC.EQ.NEL) THEN
C..Differ by none 
        CALL SCR0(NEL,NBASISMAX,NDET1,NDET2,GBS,
     &         NHG,UMAT,ALAT,ISPINSKIP,ZIA,COULCOUPLE,SUM)
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
      SUM=SUM*DFLOAT(IPTOT)
      RETURN
      END
C=================================================================
C      SUBROUTINE GETIDM(NBASISMAX,N,G,ID)
C      IMPLICIT REAL*8 (A-H,O-Z)
C      INTEGER G(4,*),NBASISMAX(4,2),N
C      ID=((G(1,N)-NBASISMAX(1,1))*(NBASISMAX(2,2)-NBASISMAX(2,1)+1)+
C     & (G(2,N)-NBASISMAX(2,1)))*(NBASISMAX(3,2)-NBASISMAX(3,1)+1)+
C     & (G(3,N)-NBASISMAX(3,1))+1
C      RETURN
C      END
