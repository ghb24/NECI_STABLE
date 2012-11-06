      SUBROUTINE GEN_XCHOLE(CIJKL,PSIR,IOBS,JOBS,KOBS,G1,SITAB,
     &    NMAX,NMSH,
     &    NHG,IXD,IYD,IZD,RHO,TCOMPLETE,XCHOLE,SPAC,ALAT,OMEGA,
     &    NMRKS,NDET,NEVAL,NEL)
      use global_utilities
      use MemoryManager, only: TagIntType
      use constants, only: dp,pi
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
C..
      integer :: IOBS,JOBS,KOBS,KI,KJ,ICC,IC,IL,IX,IY,IZ,NEL,I,J
      INTEGER :: NDET,NEVAL,NMSH,NHG,NMAX
      real(dp),ALLOCATABLE,SAVE :: CFXCS(:)
      integer(TagIntType), save :: tagCFXCS=0
      real(dp) :: CIJKL(NDET,NEVAL)
      integer :: NMRKS(NEL,*)
      real(dp) :: RHO(NMSH,NMSH,NMSH),XCHOLE(NMSH,NMSH,NMSH)
      real(dp) :: AUXC(NMSH,NMSH,NMSH)
      INTEGER G1(5,NHG)
C..Cube arrays
      real(dp) :: PSIR(-NMSH:NMSH),OMEGA,CST,CF
      real(dp) :: SITAB(NMSH,NMAX),ALAT(3)
      LOGICAL TCOMPLETE
      real(dp) :: eps,zspacing,yspacing,xspacing,spac,SUM,CHKRHO
      integer :: IXD,IYD,IZD
C..
!      PI=ACOS(-1.0_dp)
      EPS=1.0e-6_dp_dp
C..
      IF(RHO(IOBS,JOBS,KOBS).EQ.0.0_dp) THEN
        STOP ' !!! DENSITY AT REFERENCE POSITION IS ZERO !!! '
      ENDIF
C..
      allocate(CFXCS(NDET))
      call LogMemAlloc('CFXCS',NDET,8,'Gen_XCHole',tagCFXCS)
      CFXCS=0.0_dp
      CALL DCOPY(NDET,CIJKL(1,1),1,CFXCS,1)
      XCHOLE=0.0_dp
C..Evaluates second-order spinless density matrix
      DO KI=1,NDET
        IF(ABS(CFXCS(KI)).LT.EPS) GOTO 300
        DO KJ=KI,NDET
          IF(ABS(CFXCS(KJ)).LT.EPS) GOTO 500
          ICC=2
          IF(KI.EQ.KJ) ICC=1
C..IC is the number of pairs which are the *SAME* in both determinants
          IC=0
          DO I=1,NEL
            DO J=1,NEL
              IF(NMRKS(I,KJ).EQ.NMRKS(J,KI)) IC=IC+1
            ENDDO
          ENDDO
C..If determinants differ by more than 2 spin orbitals IC < NEL-2 
          IF(IC.LT.(NEL-2)) GOTO 500
          AUXC=0.0_dp
          CALL RLSXC(NEL,NMSH,SITAB,NMAX,IOBS,JOBS,KOBS,
     &      NMRKS(1:NEL,KI),NMRKS(1:NEL,KJ),G1,NHG,IC,AUXC)
C..
          CF=CFXCS(KI)*CFXCS(KJ)*ICC
C..          
          DO IZ=1,NMSH
            DO IY=1,NMSH
              DO IX=1,NMSH
                XCHOLE(IX,IY,IZ)=XCHOLE(IX,IY,IZ)+CF*AUXC(IX,IY,IZ)
              ENDDO
            ENDDO
          ENDDO
500   CONTINUE
        ENDDO
300   CONTINUE
      ENDDO
C..Evaluates exchange-correlation hole
      CST=1.0_dp/(OMEGA*OMEGA)
      CALL DSCAL(NMSH*NMSH*NMSH,0.5_dp*CST,XCHOLE,1)
      IF(TCOMPLETE) THEN
        CHKRHO=0.0_dp
        DO IZ=1,NMSH
          DO IY=1,NMSH
            DO IX=1,NMSH
C..This is an integral over the second order density matrix
              CHKRHO=CHKRHO+2.0_dp*XCHOLE(IX,IY,IZ) 
C..This is the EXCHANGE-CORRELATION hole n_xc(r,r')
              XCHOLE(IX,IY,IZ)=
     &          2.0_dp*XCHOLE(IX,IY,IZ)/RHO(IOBS,JOBS,KOBS)-RHO(IX,IY,IZ)
            ENDDO
          ENDDO
        ENDDO
C..Integral of exchange-correlation hole
        SUM=0.0_dp
        DO IZ=1,NMSH
          DO IY=1,NMSH
            DO IX=1,NMSH
              SUM=SUM+XCHOLE(IX,IY,IZ)
        ENDDO
          ENDDO
        ENDDO
        SUM=SUM*OMEGA/real(NMSH*NMSH*NMSH,dp)
        CHKRHO=CHKRHO*OMEGA/real(NMSH*NMSH*NMSH,dp)
        WRITE(6,'(A,F19.6)') ' RHO AT OBS. POINT : '
     &     ,real(NEL-1,dp)*RHO(IOBS,JOBS,KOBS)
        WRITE(6,'(A,F19.6)')  
     &      ' INTEGRAL OVER 2ND ORDER DENSITY MATRIX : ' 
     &       ,CHKRHO 
        WRITE(6,'(A,F19.6)') ' INTEGRAL OVER XC HOLE : ' , SUM
        WRITE(6,'(A,F19.6)') ' ERROR IN XCHOLE (PERCENT) : ' 
     &          ,(-1.0_dp-SUM)*100.0_dp
      ELSE
        XSPACING=ALAT(1)/NMSH*IXD
        YSPACING=ALAT(2)/NMSH*IYD
        ZSPACING=ALAT(3)/NMSH*IZD
        SPAC=SQRT(XSPACING**2+YSPACING**2+ZSPACING**2)
        DO IL=-NMSH,NMSH
          IX=IOBS+IL*IXD
          IY=JOBS+IL*IYD
          IZ=KOBS+IL*IZD
          IF(IX.LE.0.OR.IX.GT.NMSH.OR.IY.LE.0.OR.IY.GT.NMSH.OR.
     &        IZ.LE.0.OR.IZ.GT.NMSH) THEN
            PSIR(IL)=0.0_dp
          ELSE
            PSIR(IL)=2.0_dp*XCHOLE(IX,IY,IZ)/RHO(IOBS,JOBS,KOBS)
     &        -RHO(IX,IY,IZ)
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END
C ==-------------------------------------------------------------------==
      SUBROUTINE RLSXC(NEL,NMSH,SITAB,NMAX,I,J,K,
     &    NDET1,NDET2,G1,NHG,IC,AUXC)
      use constants, only: dp
!     IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      INTEGER NHG,G1(5,NHG),NEL,NMSH,NMAX
      INTEGER :: NDET1(NEL),NDET2(NEL)
      real(dp) :: AUXC(NMSH,NMSH,NMSH)
      INTEGER N1E(NHG),N2E(NHG)
C..Cube arrays
      real(dp) :: SITAB(NMSH,NMAX)
      integer :: IPTOT,IFLAG,IC,I,J,K,LEN

      N1E(1:NHG)=0
      N2E(1:NHG)=0
      CALL LINEUP(NEL,NDET1,NDET2,NHG,N1E,N2E,IFLAG,IPTOT)
C..
      IF(IC.EQ.(NEL-2)) THEN
C..Determinants differ by two spin orbitals
        CALL RLSXC2_SC(NEL,NMSH,NMAX,SITAB,I,J,K,
     &        NDET1,NDET2,G1,NHG,AUXC)
      ELSEIF(IC.EQ.(NEL-1)) THEN
C..Determinants differ by three spin orbitals
        CALL RLSXC1_SC(NEL,NMSH,NMAX,SITAB,I,J,K,
     &        NDET1,NDET2,G1,NHG,AUXC)
      ELSEIF(IC.EQ.NEL) THEN
C..Determinants are the same
        CALL RLSXC_SC(NEL,NMSH,NMAX,SITAB,I,J,K,
     &        NDET1,G1,NHG,AUXC)
      ELSE
        RETURN
      ENDIF
C..
      LEN=NMSH*NMSH*NMSH
      CALL DSCAL(LEN,real(IPTOT,dp),AUXC,1)
      RETURN
      END
