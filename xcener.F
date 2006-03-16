C     ==================================================================
      SUBROUTINE XC(MFXCX,MFXCC,RHO,EX,EC,VX,VC)
C     ==--------------------------------------------------------------==
C     ==  LDA EXCHANGE AND CORRELATION FUNCTIONALS                    ==
C     ==                                                              ==
C     ==  EXCHANGE  :  SLATER alpha                                   ==
C     ==  CORRELATION : CEPERLEY & ALDER (PERDEW-ZUNGER PARAMETERS)   ==
C     ==                VOSKO, WILK & NUSSAIR                         ==
C     ==                LEE, YANG & PARR                              ==
C     ==                PERDEW & WANG                                 ==
C     ==                WIGNER                                        ==
C     ==                HEDIN & LUNDQVIST                             ==
C     ==                ORTIZ & BALLONE (PERDEW-ZUNGER FORMULA)       ==
C     ==                ORTIZ & BALLONE (PERDEW-WANG FORMULA)         ==
C     ==                HCTH                                          ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (SMALL=1.D-10)
      PARAMETER (PI34= 0.75D0 / 3.1415 92653 58979 3D+00
     &         ,THIRD=1.D0/3.D0)
C     ==--------------------------------------------------------------==
C..Exchange
      IF(MFXCX.EQ.1) THEN
        SALPHA=2.D0/3.D0
        CALL SLATERX(RHO,EX,VX,SALPHA)
      ELSE
        EX=0.0D0
        VX=0.0D0
      ENDIF
      IF(RHO.LE.SMALL) THEN
        EC = 0.0D0
        VC = 0.0D0
        EX = 0.0D0
        VX = 0.0D0
      ELSE IF(MFXCC.EQ.1) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.1.0D0) IFLG=1
        CALL PZ(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.2) THEN
        RS = (PI34/RHO)**THIRD
        CALL VWN(RS,EC,VC)
      ELSEIF(MFXCC.EQ.3) THEN
        CALL LYP(RHO,EC,VC)
      ELSEIF(MFXCC.EQ.4) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.0.5D0) IFLG=1
        IF(RS.GT.100.D0) IFLG=3
        CALL PW(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.5) THEN
        CALL WIGNER(RHO,EC,VC)
      ELSEIF(MFXCC.EQ.6) THEN
        CALL HEDIN(RHO,EC,VC)
      ELSEIF(MFXCC.EQ.7) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.1.0D0) IFLG=1
        CALL OBPZ(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.8) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.0.5D0) IFLG=1
        IF(RS.GT.100.D0) IFLG=3
        CALL OBPW(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.9) THEN
        RS=(PI34/RHO)**THIRD
        CALL PADE(RS,EC,VC)
      ELSE
        EC=0.0D0
        VC=0.0D0
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SLATERX(RHO,EX,VX,ALPHA)
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (SMALL=1.D-10)
C..This F1 = -9/8*(3/PI)**THIRD
      PARAMETER (F1 = -1.10783814957303361D0)
      PARAMETER (THIRD=1.D0/3.D0,F43=4.D0/3.D0)
C     ==--------------------------------------------------------------==
      IF(RHO.LE.SMALL) THEN
        EX = 0.0D0
        VX = 0.0D0
      ELSE
        RS = RHO**THIRD
C..If ALPHA = 2/3 then this reduces to the Dirac Exchange formula
        EX = F1*ALPHA*RS
        VX = F43*F1*ALPHA*RS
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PZ(RS,EPZ,VPZ,IFLG)
C     ==--------------------------------------------------------------==
C     ==  J.P. PERDEW AND ALEX ZUNGER PRB 23, 5048 (1981)             ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.0311D0,B=-0.048D0,C=0.0020D0,D=-0.0116D0,
     &           GC=-0.1423D0,B1=1.0529D0,B2=0.3334D0)
C     ==--------------------------------------------------------------==
      EPWC=0.0D0
      VPWC=0.0D0
      IF(IFLG.EQ.1) THEN
C..High density formula
        XLN=LOG(RS)
        EPZ=A*XLN+B+C*RS*XLN+D*RS
        VPZ=A*XLN+(B-A/3.D0)+2.D0/3.D0*C*RS*XLN+
     &              (2.D0*D-C)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
C..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        OX=1.D0+B1*RS1+B2*RS2
        DOX=1.D0+7./6.*B1*RS1+4./3.*B2*RS2
        EPZ=GC/OX
        VPZ=EPZ*DOX/OX
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE VWN(RS,EVWN,VVWN)
C     ==--------------------------------------------------------------==
C     ==  S.H VOSKO, L.WILK, AND M. NUSAIR,                           ==
C     ==                 CAN. J. PHYS. 58 1200  (1980)                ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.0310907,B=3.72744,C=12.9352,X0=-0.10498)
      PARAMETER (TWO=2.0D0)
C     ==--------------------------------------------------------------==
      Q  = SQRT(4.D0*C-B*B)
      F1 = TWO*B/Q
      F2 = B*X0/(X0*X0+B*X0+C)
      F3 = TWO*(TWO*X0+B)/Q
      X  = SQRT(RS)
      FX = X*X+B*X+C
      QX = ATAN(Q/(TWO*X+B))
      EVWN=A*(LOG(RS/FX)+F1*QX-F2*(LOG((X-X0)**2/FX)+F3*QX))
      TXPB=TWO*X+B
      TTQQ=TXPB*TXPB+Q*Q
      VVWN=EVWN - X*A/6.*(TWO/X-TXPB/FX-4.*B/TTQQ-F2*(TWO/(X-X0)
     &          -TXPB/FX-4.*(TWO*X0+B)/TTQQ))
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE LYP(RHO,ELYP,VLYP)
C     ==--------------------------------------------------------------==
C     ==  C. LEE, W. YANG, AND R.G. PARR, PRB 37, 785 (1988)          ==
C     ==  THIS IS ONLY THE LDA PART                                   ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.04918,B=0.132,C=0.2533,D=0.349)
      PARAMETER (CF=2.87123400018819108)
C     ==--------------------------------------------------------------==
      RS=RHO**(-1./3.)
      ECRS=B*CF*EXP(-C*RS)
      OX=1.D0/(1.D0+D*RS)
      ELYP=-A*OX*(1.D0+ECRS)
      VLYP=ELYP-RS/3.*A*OX*(D*OX+ECRS*(D*OX+C))
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PW(RS,EPWC,VPWC,IFLG)
C     ==--------------------------------------------------------------==
C     ==  J.P. PERDEW AND YUE WANG PRB 45, 13244 (1992)               ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.031091,A1=0.21370,B1=7.5957,B2=3.5876,B3=1.6382,
     &           B4=0.49294,C0=A,C1=0.046644,C2=0.00664,C3=0.01043,
     &           D0=0.4335,D1=1.4408)
C     ==--------------------------------------------------------------==
      EPWC=0.0D0
      VPWC=0.0D0
      IF(IFLG.EQ.1) THEN
C..High density formula
        XLN=LOG(RS)
        EPWC=C0*XLN-C1+C2*RS*XLN-C3*RS
        VPWC=C0*XLN-(C1+C0/3.D0)+2.D0/3.D0*C2*RS*XLN-
     &              (2.D0*C3+C2)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
C..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        RS3=RS2*RS1
        RS4=RS2*RS2
        OM=2.D0*A*(B1*RS1+B2*RS2+B3*RS3+B4*RS4)
        DOM=2.D0*A*(0.5*B1*RS1+B2*RS2+1.5D0*B3*RS3+2.D0*B4*RS4)
        OLOG=LOG(1.D0+1.0D0/OM)
        EPWC=-2.D0*A*(1+A1*RS)*OLOG
        VPWC=-2.*A*(1.D0+2./3.*A1*RS)*OLOG
     &       -2./3.*A*(1.+A1*RS)*DOM/(OM*(OM+1.))
      ELSEIF(IFLG.EQ.3) THEN
C..Low density formula
        EPWC=-D0/RS+D1/RS**1.5D0
        VPWC=-4.D0/3.D0*D0/RS+1.5D0*D1/RS**1.5
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE WIGNER(RHO,EXC,FXC)
      IMPLICIT REAL*8 (A-H,O-Z)
      RH=RHO
      X=RH**0.33333333333333333D0
      FXC=-X*((0.943656+8.8963*X)/(1.0+12.57*X)**2)
      EXC=-0.738*X*(0.959/(1.0+12.57*X))
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE HEDIN(RHO,ECP,FCP)
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE RH
      IF(RH .EQ. 0.0D0) RETURN
      RH=RHO
      RSM1=0.62035049D0*RH**(0.3333333333333333D0)
      ALN=DLOG(1.0D0 + 21.0D0*RSM1)
      X=21.0D0/RSM1
      ECP = ALN+(X**3*ALN-X*X)+X/2-1.0D0/3.0D0
      ECP = -0.0225*ECP
      FCP = -0.0225*ALN
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE OBPZ(RS,EPZ,VPZ,IFLG)
C     ==--------------------------------------------------------------==
C     ==  G.ORTIZ AND P. BALLONE PRB 50, 1391 (1994)                  ==
C     ==  PERDEW-ZUNGER FORMULA                                       ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.031091D0,B=-0.046644D0,C=0.00419D0,D=-0.00983D0,
     &           GC=-0.103756D0,B1=0.56371D0,B2=0.27358D0)
C     ==--------------------------------------------------------------==
      EPWC=0.0D0
      VPWC=0.0D0
      IF(IFLG.EQ.1) THEN
C..High density formula
        XLN=LOG(RS)
        EPZ=A*XLN+B+C*RS*XLN+D*RS
        VPZ=A*XLN+(B-A/3.D0)+2.D0/3.D0*C*RS*XLN+
     &              (2.D0*D-C)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
C..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        OX=1.D0+B1*RS1+B2*RS2
        DOX=1.D0+7./6.*B1*RS1+4./3.*B2*RS2
        EPZ=GC/OX
        VPZ=EPZ*DOX/OX
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE OBPW(RS,EPWC,VPWC,IFLG)
C     ==--------------------------------------------------------------==
C     ==  G.ORTIZ AND P. BALLONE PRB 50, 1391 (1994)                  ==
C     ==  PERDEW-WANG FORMULA                                         ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.031091,A1=0.026481,B1=7.5957,B2=3.5876,B3=-0.46647,
     &           B4=0.13354,C0=A,C1=0.046644,C2=0.00664,C3=0.01043,
     &           D0=0.4335,D1=1.4408)
C     ==--------------------------------------------------------------==
      EPWC=0.0D0
      VPWC=0.0D0
      IF(IFLG.EQ.1) THEN
C..High density formula
        XLN=LOG(RS)
        EPWC=C0*XLN-C1+C2*RS*XLN-C3*RS
        VPWC=C0*XLN-(C1+C0/3.D0)+2.D0/3.D0*C2*RS*XLN-
     &              (2.D0*C3+C2)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
C..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        RS3=RS2*RS1
        RS4=RS2*RS2
        OM=2.D0*A*(B1*RS1+B2*RS2+B3*RS3+B4*RS4)
        DOM=2.D0*A*(0.5*B1*RS1+B2*RS2+1.5D0*B3*RS3+2.D0*B4*RS4)
        OLOG=LOG(1.D0+1.0D0/OM)
        EPWC=-2.D0*A*(1+A1*RS)*OLOG
        VPWC=-2.*A*(1.D0+2./3.*A1*RS)*OLOG
     &       -2./3.*A*(1.+A1*RS)*DOM/(OM*(OM+1.))
      ELSEIF(IFLG.EQ.3) THEN
C..Low density formula
        EPWC=-D0/RS+D1/RS**1.5D0
        VPWC=-4.D0/3.D0*D0/RS+1.5D0*D1/RS**1.5
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PADE(RS,EC,VC)
C     ==--------------------------------------------------------------==
C     ==  PADE APPROXIMATION                                          ==
C     ==  S. GOEDECKER, M. TETER, J. HUTTER, PRB in press             ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A0=0.4581652932831429D0,A1=2.217058676663745D0,
     &           A2=0.7405551735357053D0,A3=0.01968227878617998D0)
      PARAMETER (B1=1.0000000000000000D0,B2=4.504130959426697D0,
     &           B3=1.110667363742916D0,B4=0.02359291751427506D0)
      PARAMETER (O3=1.D0/3.D0)
C     ==--------------------------------------------------------------==
      TOP=A0+RS*(A1+RS*(A2+RS*A3))
      DTOP=A1+RS*(2.D0*A2+3.D0*A3*RS)
      BOT=RS*(B1+RS*(B2+RS*(B3+RS*B4)))
      DBOT=B1+RS*(2.D0*B2+RS*(3.D0*B3+RS*4.D0*B4))
      EC=-TOP/BOT
      VC=EC+RS*O3*(DTOP/BOT-TOP*DBOT/(BOT*BOT))
C     ==--------------------------------------------------------------==
      RETURN
      END
