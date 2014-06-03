C.. Factor a polynomial of order ORDER into complex roots ZROOTS
      SUBROUTINE FACTORPOLY(COEFFS,ORDER,ZROOTS)
         use constants, only : dp
C.. Use Newton-Raphson to find roots and then divide them out one by one
C         IMPLICIT NONE
         INTEGER ORDER,IROOTS,WORDER,I
         real(dp) COEFFS(0:ORDER)
         real(dp) C(0:ORDER),POLY(0:2),REM(0:2)
         complex(dp) ZROOTS(ORDER),ZROOT
         LOGICAL LCOMPLEX
         LCOMPLEX=.FALSE.
         DO I=0,ORDER
            C(I)=COEFFS(I)
         ENDDO
         IROOTS=0
         WORDER=ORDER
C         WRITE(6,*) "POLY=",(WCOEFFS(I),I=0,WORDER)
         DO WHILE (IROOTS.LT.ORDER)
C            IF(MOD(WORDER,2).EQ.0) THEN
               CALL FINDROOT(C,WORDER,ZROOT,LCOMPLEX)
C            ELSE
C.. we've an odd order poly, so we can use false descent
C               CALL FINDROOTFD(C,WORDER,XROOT)
C               ZROOT=DCMPLX(XROOT,0.0_dp)
C            ENDIF
            IROOTS=IROOTS+1
            IF(AIMAG(ZROOT).EQ.0.0_dp) THEN
C.. a real root
               ZROOTS(IROOTS)=ZROOT
               POLY(0)=-DBLE(ZROOT)
               POLY(1)=1.0_dp
               CALL POLYDIV(C,WORDER,POLY,1,REM)
            ELSE
C.. we have a complex root, so we use its conjugate as well
               ZROOTS(IROOTS)=ZROOT
               IROOTS=IROOTS+1
               ZROOTS(IROOTS)=CMPLX(DBLE(ZROOT),-AIMAG(ZROOT),kind=dp)
C.. we divide out the product of ZROOT and ZROOT*
               POLY(0)=DBLE(ZROOT)*DBLE(ZROOT)+AIMAG(ZROOT)*AIMAG(ZROOT)
               POLY(1)=-2.0_dp*DBLE(ZROOT)
               POLY(2)=1.0_dp
               CALL POLYDIV(C,WORDER,POLY,2,REM)
            ENDIF
            IF(ABS(REM(0)).GT.1.0e-5_dp) THEN
               WRITE(6,*) 'POLYNOMIAL FACTORIZATION FAILED'
               WRITE(6,*) "POLY=",(COEFFS(I),I=0,ORDER)
   
               WRITE(6,*) "REM=",REM(0)
               WRITE(6,*) "Roots=",(ZROOTS(I),I=1,IROOTS)
               WRITE(6,*) "Quotient=",(C(I),I=0,WORDER-1)
               STOP 'POLYNOMIAL FACTORIZATION FAILED'
            ENDIF
         ENDDO
C.. All should be well now
         RETURN      
      END

C.. Divide POLY1 of order ORDER1 by POLY2 of order ORDER2, placing the
C,, result in POLY1, and the remainder in REM (of order ORDER2-1)
      SUBROUTINE POLYDIV(POLY1,ORDER1,POLY2,ORDER2,REM)
         use constants, only : dp
         IMPLICIT NONE
         INTEGER ORDER1,ORDER2
         real(dp) POLY1(0:ORDER1),POLY2(0:ORDER2),REM(0:ORDER2-1)
         real(dp) T
         INTEGER I,J
         
         DO I=ORDER1,ORDER2,-1
C.. Place the result of this stage of the division in POLY1(I) 
            T=POLY1(I)/POLY2(ORDER2)
            POLY1(I)=T
            DO J=I-1,I-ORDER2,-1
               POLY1(J)=POLY1(J)-T*POLY2(ORDER2-(I-J))
            ENDDO
         ENDDO
C.. Now just move the remainder to REM
         DO I=0,ORDER2-1
            REM(I)=POLY1(I)
         ENDDO
C.. the result is of length ORDER1-ORDER2
         ORDER1=ORDER1-ORDER2
         DO I=0,ORDER1
            POLY1(I)=POLY1(I+ORDER2)
         ENDDO
         RETURN
      END

C.. Use Newton Raphson to find a root
C.. if LCOMPLEX we start in the complex plane.
C.. if we need to go to the complex plane, we set LCOMPLEX
      SUBROUTINE FINDROOT(WCOEFFS,WORDER,ZROOT,LCOMPLEX)
         use constants, only : dp
         IMPLICIT NONE
         INTEGER WORDER
         real(dp) WCOEFFS(0:WORDER),MAG,MAG2,REPS
         complex(dp) ZROOT,ZEPS,Z1,Z2,Z3,ZOR,Z4
         INTEGER ICOUNT,I
         LOGICAL LCOMPLEX,TBIS
         IF(LCOMPLEX) THEN
            ZROOT=(1.0_dp,0.1_dp)
         ELSE
C.. As our roots are close to 1, so are possible minima which will throw 
C.. the NR method.  we start quite far away.
            ZROOT=(1.01_dp,0.0_dp)
         ENDIF
         REPS=1.0e-14_dp
         ICOUNT=0
         CALL POLYVALUE(WCOEFFS,WORDER,ZROOT,ZEPS)
         MAG=SQRT(DBLE(ZEPS)*DBLE(ZEPS)+AIMAG(ZEPS)*AIMAG(ZEPS))
         Z2=(0.1_dp,0.0_dp)
C         WRITE(6,*) "POLY=",(WCOEFFS(I),I=0,WORDER)
         DO WHILE (MAG.GT.REPS)
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GE.100.AND.ICOUNT.LE.145.AND.MOD(ICOUNT,13).EQ.0)
     &          REPS=REPS*10.0_dp
            IF(ICOUNT.EQ.150) THEN
C.. if we can't converge on the real axis, we try the complex plane
               LCOMPLEX=.TRUE.
               ZROOT=(1.0_dp,1.0_dp)
            ENDIF
            IF(ICOUNT.EQ.550) THEN
               WRITE(6,*) "POLY=",(WCOEFFS(I),I=0,WORDER)
               STOP 'Newton-Raphson did not converge'
            ENDIF
            CALL POLYDERIV(WCOEFFS,WORDER,ZROOT,Z1)
            IF(ICOUNT.GT.260) THEN
               WRITE(6,*) "I,DERIV,Z,ZROOT",ICOUNT,Z1,ZEPS,ZROOT
            ENDIF
            TBIS=.FALSE.
C.. save the deriv for next time
               Z2=Z1
C            WRITE(6,*) "VALUE=",ZEPS
               ZEPS=ZEPS/Z1
               Z4=ZEPS/Z1            
            IF (SQRT(DBLE(ZEPS)**2+AIMAG(ZEPS)**2).GT.1.D2) THEN
C.. the deriv is too close to zero.  use bisection
               ZEPS=ZROOT-ZOR
               TBIS=.TRUE.
            ENDIF
C            IF (SQRT(DBLE(Z1)**2+AIMAG(Z1)**2).LT.1.0e-2_dp
C     &            .AND.ICOUNT.GE.60) Z3=Z3*Z1 
C            WRITE(6,*) ICOUNT,ZROOT,ZEPS
            ZOR=ZROOT
            ZROOT=ZROOT-ZEPS
C.. every 10 moves after 50, we assume we're bouncing around between two points,
C.. so we aid convergence by only going by ZEPS/2
            IF(ICOUNT.GE.30.AND.ICOUNT.LT.150.AND.MOD(ICOUNT,6).EQ.0
     &         .OR.ICOUNT.GE.290.AND.MOD(ICOUNT,6).EQ.0) TBIS=.TRUE.
            IF(TBIS) ZROOT=ZROOT+ZEPS/2.0_dp
            CALL POLYVALUE(WCOEFFS,WORDER,ZROOT,ZEPS)
            Z3=ZEPS
            IF(AIMAG(Z3).NE.0.0_dp) THEN
               MAG=SQRT(DBLE(Z3)*DBLE(Z3)+AIMAG(Z3)*AIMAG(Z3))
            ELSE
               MAG=ABS(DBLE(Z3))
            ENDIF
C            WRITE(6,*) MAG,REPS,ZEPS
         ENDDO
C.. Check to see if just the real version is better than the complex
         Z3=DBLE(ZROOT)
         CALL POLYVALUE(WCOEFFS,WORDER,Z3,Z1)
         CALL POLYVALUE(WCOEFFS,WORDER,ZROOT,Z2)
         MAG=ABS(DBLE(Z1))
         MAG2=SQRT(DBLE(Z2)*DBLE(Z2)+AIMAG(Z2)*AIMAG(Z2))
C         WRITE(6,*) MAG,MAG2
         IF(MAG.LT.MAG2) ZROOT=Z3         
C         WRITE(6,*) "ROOT=",ZROOT,"ICOUNT=",ICOUNT,MAG,MAG2         
         RETURN
      END

C.. return the value of the polynomial at Z
      SUBROUTINE POLYVALUE(WCOEFFS,WORDER,Z,ZRET)
         use constants, only : dp
         IMPLICIT NONE
         INTEGER WORDER,I
         real(dp) WCOEFFS(0:WORDER)
         complex(dp) Z,ZRET
         ZRET=(0.0_dp,0.0_dp)
         DO I=WORDER,0,-1
            ZRET=ZRET*Z+WCOEFFS(I)
         ENDDO
         RETURN
      END

C.. return the value of the polynomial at Z
      SUBROUTINE POLYVALUER(WCOEFFS,WORDER,X,XRET)
         use constants, only : dp
         IMPLICIT NONE
         INTEGER WORDER,I
         real(dp) WCOEFFS(0:WORDER)
         real(dp) X,XRET
         XRET=0.0_dp
         DO I=WORDER,0,-1
            XRET=XRET*X+WCOEFFS(I)
         ENDDO
         RETURN
      END

C.. return the derivative of the polynomial at Z
      SUBROUTINE POLYDERIV(WCOEFFS,WORDER,Z,ZRET)
         use constants, only : dp
         IMPLICIT NONE
         INTEGER WORDER,I
         real(dp) WCOEFFS(0:WORDER)
         complex(dp) Z,ZRET
         ZRET=(0.0_dp,0.0_dp)
         DO I=WORDER,1,-1
            ZRET=ZRET*Z+WCOEFFS(I)*I
         ENDDO
         RETURN
      END
C.. return the derivative of the polynomial at Z
      SUBROUTINE POLYNDERIV(WCOEFFS,WORDER,Z,ZRET,N)
         use constants, only : dp
         IMPLICIT NONE
         INTEGER WORDER,I,N,K
         real(dp) WCOEFFS(0:WORDER),D
         complex(dp) Z,ZRET
         ZRET=(0.0_dp,0.0_dp)
         DO I=WORDER,1,-1
            D=1.0_dp
            DO K=0,N-1
               D=(D*(I-K))/(K+1)
            ENDDO
            ZRET=ZRET*Z+WCOEFFS(I)*D
         ENDDO
         RETURN
      END
C.. False descent root finding - only suitable for odd order polys
C.. if we need to go to the complex plane, we set LCOMPLEX
      SUBROUTINE FINDROOTFD(WCOEFFS,WORDER,XROOT)
         use constants, only : dp
         IMPLICIT NONE
         INTEGER WORDER
         real(dp) WCOEFFS(0:WORDER),MAG
         real(dp) XROOT,XEPS,X1,X2,X3,X4,X5,A,B
         INTEGER ICOUNT,I
C.. Choose the bracket to be big
         A=-2
         B=1 
         ICOUNT=0
         CALL POLYVALUER(WCOEFFS,WORDER,A,X1)
         CALL POLYVALUER(WCOEFFS,WORDER,B,X2)
         DO WHILE(X1*X2.GT.0.0_dp)
            A=A*2
            B=B*2
            CALL POLYVALUER(WCOEFFS,WORDER,A,X1)
            CALL POLYVALUER(WCOEFFS,WORDER,B,X2)
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.40) THEN
               WRITE(6,*) "POLY=",(WCOEFFS(I),I=0,WORDER)
               STOP 'False descent could not find suitable bracket.'
            ENDIF 
         ENDDO
         ICOUNT=0
         MAG=1.0_dp;
C         WRITE(6,*) "POLY=",(WCOEFFS(I),I=0,WORDER)
         DO WHILE (MAG.GT.1.0e-9_dp)
            ICOUNT=ICOUNT+1
            IF(ICOUNT.EQ.200) THEN
               WRITE(6,*) "POLY=",(WCOEFFS(I),I=0,WORDER)
               STOP 'False Descent did not converge'
            ENDIF
C.. Save old posn
            X5=X3
            X3=(X2*A-X1*B)/(X2-X1)
            CALL POLYVALUER(WCOEFFS,WORDER,X3,X4)
            WRITE(6,*) A,B,X1,X2,X3,X4
            IF(X4*X1.GT.0.0_dp) THEN
C.. if X4 and X1 have the same sign
               A=X3
               X1=X4
            ELSE
               B=X3
               X2=X4
            ENDIF
            XEPS=X4
            MAG=ABS(XEPS)
         ENDDO
C.. Check to see if just the real version is better than the complex
         XROOT=X3
         RETURN
      END

C.. take x_n z^n+x_(n-1) z^(n-1)+...+x_1 z + x_0
C.. and substitute = z= z0+z'
C.. to create a new polynomial in z' 
C.. coeffs x'_n = Sum_i(x_i z0^i Binomial[i,n])/z0^n
      SUBROUTINE CREATEAUXPOLY(WCOEFFS,WORDER,Z0,ZACOEFFS)
         use constants, only : dp
         INTEGER WORDER
         real(dp) WCOEFFS(0:WORDER)
         complex(dp) ZACOEFFS(0:WORDER),Z0
         INTEGER I,J
         real(dp) CFACT
         DO I=0,WORDER
            ZACOEFFS(I)=(0.0_dp,0.0_dp)
         ENDDO
         DO I=WORDER,0,-1
            CFACT=1.0_dp
            DO J=0,I
               ZACOEFFS(J)=ZACOEFFS(J)*Z0+WCOEFFS(I)*CFACT
               CFACT=(CFACT*(I-J))
               CFACT=CFACT/(J+1)
            ENDDO
         ENDDO
         RETURN
      END

      SUBROUTINE POLYVALUEC(ZCOEFFS,WORDER,Z,ZRET)
         use constants, only : dp
         IMPLICIT NONE
         INTEGER WORDER,I
         complex(dp) ZCOEFFS(0:WORDER)
         complex(dp) Z,ZRET
         ZRET=(0.0_dp,0.0_dp)
         DO I=WORDER,0,-1
            ZRET=ZRET*Z+ZCOEFFS(I)
         ENDDO
         RETURN
      END

