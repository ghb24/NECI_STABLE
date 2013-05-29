! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
C.. Perform the contour integration to calculate the Z-sum
      SUBROUTINE ZSUM(XIS,I_H,I_P,I_POINTS,R,SUM,SUMI)
         use constants, only: dp,pi
         complex(dp) Z,ZI,ZRET,ZS,DZ
         INTEGER I_H,I_P
         real(dp) THETA,R,SUM,SUMI,DTHETA,P
         real(dp) XIS(1:I_H)
C.. Simpson's Rule: Integral(f(x)) = (1/3)h(f0+4(f1+f3+...+f(2n-1))
C..                                  +2(f2+f4+...+f(2n-2)+f2n
C.. We use I_POINTS+1, but have both THETA=0 and 2PI, so we divide the 
C.. unit circle into I_POINTS 
         ZI=(0,1)
         P=I_P
         DTHETA=2*PI/I_POINTS
         DZ=EXP(ZI*DTHETA)
         THETA=0.0_dp
         Z=(0.0_dp,0.0_dp)
         Z=Z+R
C         WRITE(6,*),DTHETA,DZ
         ZS=(0.0_dp,0.0_dp)
         DO I=1,I_POINTS/2
C.. Calculate an even point
C            WRITE(6,*),Z
            CALL ZPRIMITIVE(XIS,I_H,P,Z,R,THETA,ZRET)
            ZS=ZS+2*ZRET
            IF(I.EQ.1) ZS=ZS-ZRET
C.. Move round DTHETA in the complex plane
            THETA=THETA+DTHETA
            Z=Z*DZ
C            WRITE(6,*),Z
C.. Calculate an odd point
            CALL ZPRIMITIVE(XIS,I_H,P,Z,R,THETA,ZRET)
            ZS=ZS+4*ZRET
C.. Move round DTHETA in the complex plane
            THETA=THETA+DTHETA
            Z=Z*DZ
         ENDDO
C         WRITE(6,*),Z
C.. Now do the last point, f2n
         CALL ZPRIMITIVE(XIS,I_H,P,Z,R,THETA,ZRET)
         ZS=ZS+ZRET
         ZS=ZS*DTHETA/3.0_dp
         ZS=ZS/(2.0_dp*PI*ZI)
         SUM=DBLE(ZS)
         SUMI=AIMAG(ZS)
      END

C.. Calculate a primitive element in the contour integration
C.. NB P is no longer integer, but is REAL
      SUBROUTINE ZPRIMITIVE(XIS,I_H,P,Z,R,THETA,ZRET)
         use constants, only: dp
         INTEGER I_H,I
         real(dp) P
         complex(dp) Z,ZI,ZREIT,ZRET
         real(dp) XIS(1:I_H)
         real(dp) R,THETA
         ZI=(0,1)
         ZREIT=R*EXP(ZI*THETA)
         ZRET=R*ZI*EXP(ZI*THETA)*(((ZREIT)**P)-1.0_dp)
         ZRET=ZRET/(ZREIT-1)
         DO I=1,I_H
            ZRET=ZRET/(Z-XIS(I))
         ENDDO
         RETURN
      END
C.. Calculate a primitive element in the contour integration - handle
C.. the XIS
C.. NB P is no longer integer, but is REAL
      SUBROUTINE ZPRIM1(XIS,I_H,Z,ZRET)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER I_H,I
         complex(dp) Z,ZRET
         real(dp) XIS(1:I_H)
         ZRET=(1,0)
         DO I=1,I_H
            ZRET=ZRET/(Z-XIS(I))
         ENDDO
         RETURN
      END
C.. Calculate a primitive element in the contour integration - handle
C.. the Z^P section
C.. NB P is no longer integer, but is REAL
      SUBROUTINE ZPRIM2(P,Z,ZMZ0,ZRET)
         use constants, only: dp
         IMPLICIT NONE
         real(dp) P
         complex(dp) Z,ZI,ZMZ0,ZRET,ZR
         LOGICAL LSAME2
         ZI=(0,1)
         IF(LSAME2(Z,1.0_dp)) THEN
            WRITE(25,*) P,(Z**P-1)/(Z-1)
            ZRET=ZRET*ZI*ZMZ0*P
         ELSE
            ZR=ZI*ZMZ0*(((Z)**P)-1.0_dp)
            ZRET=ZRET*ZR/(Z-1)
         ENDIF
         RETURN
      END

C.. ZSFUN(ZSD1,ZSD2,I_H,P,Z,R,THETA,ZRET) is a user-defined function
C.. ZSD1 and ZSD2 are user defined data
C.. Perform the contour integration to calculate the Z-sum
      SUBROUTINE ZSUM2(ZSD1,ZSD2,ZSFUN,Z0,I_P,I_POINTS,R,SUM,SUMI)
         use constants, only: dp,pi
         IMPLICIT NONE
         real(dp) ZSD1(*),ZSD2(*),ZSFUN
         complex(dp) Z,ZI,ZS,DZ,Z0,Z2
         INTEGER I_P,I_POINTS,I
         real(dp) THETA,R,SUM,SUMI,DTHETA,P
         character(*), parameter :: t_r="ZSUM2"
C.. Simpson's Rule: Integral(f(x)) = (1/3)h(f0+4(f1+f3+...+f(2n-1))
C..                                  +2(f2+f4+...+f(2n-2)+f2n
C.. We use I_POINTS+1, but have both THETA=0 and 2PI, so we divide the 
C.. unit circle into I_POINTS 
         ZI=(0,1)
         P=I_P
         DTHETA=2*PI/I_POINTS
         DZ=EXP(ZI*DTHETA)
         THETA=0.0_dp
         Z=(0.0_dp,0.0_dp)
         Z=Z+R
         Z2=Z0+Z
C         WRITE(6,*),DTHETA,DZ
         ZS=(0.0_dp,0.0_dp)
         DO I=1,I_POINTS/2
C.. Calculate an even point
C            WRITE(6,*),Z
             call stop_all(t_r,"ZSFUN should be a function!")
!            CALL ZSFUN(ZSD1,ZSD2,P,Z2,R,THETA,Z,ZRET)
!            ZS=ZS+2*ZRET
!            IF(I.EQ.1) ZS=ZS-ZRET
C.. Move round DTHETA in the complex plane
            THETA=THETA+DTHETA
            Z=Z*DZ
            Z2=Z0+Z
C            WRITE(6,*),Z
C.. Calculate an odd point
!            CALL ZSFUN(ZSD1,ZSD2,P,Z2,R,THETA,Z,ZRET)
!            ZS=ZS+4*ZRET
C.. Move round DTHETA in the complex plane
            THETA=THETA+DTHETA
            Z=Z*DZ
            Z2=Z0+Z
         ENDDO
C         WRITE(6,*),Z
C.. Now do the last point, f2n
!         CALL ZSFUN(ZSD1,ZSD2,P,Z2,R,THETA,Z,ZRET)
!         ZS=ZS+ZRET
         ZS=ZS*DTHETA/3.0_dp
         ZS=ZS/(2.0_dp*PI*ZI)
         SUM=DBLE(ZS)
         SUMI=AIMAG(ZS)
C         WRITE(6,*) SUM,SUMI
      END

C.. Calculate the contour integral around ZSFUN by integrating around
C.. each of the ZPOLES(I).  ZPOLES(I) is not necessarily unique.
C.. The radius of the contour is |ZPOLES(I)-ZPOLES(ICLOSEST)|/P
      SUBROUTINE ZSUM3(ZSD1,ZSD2,ZSFUN,I_POINTS,I_P,IPOLES,
     &            ZPOLES,TAUXPOLY,SUM,SUMI)
C.. Go through the list of ZPOLES.  For each pole, find the closest
C.. neighbour distance (removing poles which are within a distance EPS)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER I_P,IPOLES,I_POINTS
         real(dp) ZSD1(*),ZSD2(0:IPOLES+1),ZSFUN
         complex(dp) ZPOLES(IPOLES)
         LOGICAL LPOLES(IPOLES),TAUXPOLY
         real(dp) SUM,SUMI,S,SI
         complex(dp) ZAUXPOLY(0:IPOLES+1)
      
         INTEGER I,J,ICOUNT
         real(dp) DIST,DISTMIN,EPS,DISTMAX,DM
         complex(dp) ZD,Z0
         character(*), parameter :: t_r="ZSUM3"
      
         SUM=0.0_dp
         SUMI=0.0_dp

C.. two poles are the same if they are within EPS of each other
         EPS=1.0e-6_dp
C.. in fact, all that needs to be done is to take a contour around all
C.. poles which are similarly close to each other - we do this by
C.. computing their mean, and taking a radius which is at least as large
C.. the max distance between them.

         DO I=1,IPOLES
C.. mark all the poles as to be processed
            LPOLES(I)=.TRUE.
         ENDDO

         DO I=1,IPOLES
            IF(LPOLES(I)) THEN
               Z0=ZPOLES(I)
               DISTMIN=1.D100
               DISTMAX=0.0_dp
               ICOUNT=1
               DO J=1,IPOLES
                  IF(I.NE.J) THEN
                     ZD=ZPOLES(J)-ZPOLES(I)
                     DIST=SQRT(DBLE(ZD)*DBLE(ZD)+AIMAG(ZD)*AIMAG(ZD))
                     IF(DIST.LE.EPS) THEN
                        LPOLES(J)=.FALSE.
                        IF(DIST.GT.DISTMAX) DISTMAX=DIST
                        Z0=Z0+ZPOLES(J)
                        ICOUNT=ICOUNT+1
                     ELSEIF(DIST.LT.DISTMIN) THEN
                        DISTMIN=DIST
                     ENDIF
                  ENDIF
               ENDDO
               Z0=Z0/ICOUNT
               IF(IPOLES.EQ.1) DISTMIN=1e-5
C.. we have the distance to the closest pole.  Divide it by P
               DM=DISTMIN/(10*I_P)
C.. if we go too small, we might not be on our pole, so we just take
C.. half the distance to the closest pole
               IF(DM.LT.1D-8) DM=DISTMIN/2.0_dp
C/(2*I_P)
C.. If this calculated distance is closer than the furthest pole which
C.. wish to count in this integral, we increase the distance
               IF(DM.LT.DISTMAX) DM=DISTMAX*1.1_dp
C.. do the contour integral
               IF(TAUXPOLY) THEN
                  CALL CREATEAUXPOLY(ZSD2(1),IPOLES,
     &                     ZPOLES(I),ZAUXPOLY(1))
C                  WRITE(30,"(8G25.16e3)",advance='no') (ZAUXPOLY(J),J=5,1,-1)
C                  WRITE(30,*)
                   ZAUXPOLY(0)=ZSD2(0)
                   call stop_all(t_r,"ZAUXPOLY inconsistent data type")
!                  CALL ZSUM2(ZSD1,ZAUXPOLY,ZSFUN,ZPOLES(I),I_P,I_POINTS,
!     &               DM,S,SI)
               ELSE
                  CALL ZSUM2(ZSD1,ZSD2,ZSFUN,ZPOLES(I),I_P,I_POINTS,
     &               DM,S,SI)
               ENDIF
C               WRITE(31,"(A2,5G25.16e3,A2)",advance='no') "(",ZPOLES(I),DM,S,SI,")"
               SUM=SUM+S
               SUMI=SUMI+SI
            ENDIF
         ENDDO
C         WRITE(31,*)
         RETURN
      END
