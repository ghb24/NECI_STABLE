! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi

      SUBROUTINE PLOTROOTSTAR(N,A,B)
         use constants, only: dp
         use legacy_data, only: CALCP_N, CALCP_NT, CALCP_SUM, CALCP_HSUM
         IMPLICIT NONE
         INTEGER N,I,J
         HElement_t A(0:N),B(0:N)
         real(dp) HIGHZ,DELTA,Z
         real(dp) STAR,STARPROD,PROD,PROD2,STARPROD2

         !Calculate the poles of the new product formulation


         Z=0.9994
         HIGHZ=0.99941
         DELTA=0.0000001
!         WRITE(97,*) N," excitations" 
!         WRITE(97,*) "rhoii   rhoij"
!         DO I=1,N
!            WRITE(97,*) (A(I)), (B(I))
!         ENDDO
!         CALL neci_flush(97)
         
         WRITE(66,*) "Z     STARVAL     PRODVAL    PRODVALALL"
         DO WHILE(Z.LE.HIGHZ)
            STAR=0.0_dp
            DO I=1,N
                STAR=STAR+abs(B(I))**2/(Z-(A(I)))
            ENDDO
            STARPROD=0.0_dp
            DO I=1,N
                PROD=0.0_dp
                DO J=1,N
                    IF(J.eq.I) CYCLE
!                    WRITE(97,*) I,J
                    PROD=PROD+abs(B(J))**2/(Z-((A(I))*(A(J))))
                ENDDO
                STARPROD=STARPROD+abs(B(I))**2/(Z-(A(I))-PROD)
            ENDDO
            STARPROD2=0.0_dp
            DO I=1,N
                PROD2=0.0_dp
                DO J=1,N
                    PROD2=PROD2+ 
     &                        abs(B(J))**2/(Z-((A(I))*(A(J))))
                ENDDO
                STARPROD2=STARPROD2+abs(B(I))**2/(Z-(A(I))-PROD2)
            ENDDO
            Z=Z+DELTA
            WRITE(66,"(4F14.10)") Z,STAR,STARPROD,STARPROD2
         ENDDO
         CALL neci_flush(66)
         RETURN
      END
                
            

C.. A procedure to use the regular falsi (false descent) method
C.. to find the roots of a polynomial in the form

C.. a0-x + b1^2/(x-a1) +b2^2/(x-a2) + ...

C.. which is the characteristic polynomial of a matrix corresponding to
C.. a star-graph, with a as the diagonal terms, and b as the off
C.. diagonal (first row) terms.  a0 is the diagonal of the pivot of the
C.. star.  THere are N spokes and 1 pivot, and thus N+1 roots

C.. A and B are sorted in ascending order of A

      SUBROUTINE FINDROOTSTAR(N,A,B,ROOTS,NROOTS)
         use constants, only: dp
         use IntegralsData , only : NRCONV,RFCONV,NRSTEPSMAX
         IMPLICIT NONE
         INTEGER N,I,J,NROOTS,NMIN,K
         HElement_t A(0:N),B(0:N)
         real(dp) ROOTS(0:NROOTS),INITHIGH,INITLOW,DELTA
         real(dp) XHIGH,XLOW,X,VAL,DERIV,OX,VALLOW,VALHIGH
         LOGICAL LOWPOLE,HIGHPOLE
         NMIN=N-NROOTS
!         write(76,*) n,nRoots
!         do i=0,n
!            write(76,*) i,a(i),b(i)
!         enddo
!         write(76,*)
C.. We find the first root by starting very slightly above the highest
C.. pole
         roots(0)=0
         I=NROOTS
!         WRITE(76,*) "Total roots = ",NROOTS
!  This will be bounded from below by A(I+NMIN)
         XLOW=(A(I+NMIN))
         XHIGH=2*(A(I+NMIN))
C.. The value at XHIGH is always bigger than at XLOW

!  When I=0, there is no min.
!I points to the root we're finding, starting from the highest, NROOTS, and working down.
         DO WHILE(I.GT.0)
             
!            WRITE(76,*),I,A(I+NMIN),B(I+NMIN)
            VAL=1.0_dp
            J=0
            K=0
            LOWPOLE=.true.
            HIGHPOLE=.true.
            DELTA=1.0e-13_dp
            INITHIGH=XHIGH
            INITLOW=XLOW
!            WRITE(55,*) "***NEW ROOT***",I
!            WRITE(55,*) "INITIAL XHIGH/XLOW: ", XHIGH, XLOW
            !Do some Bisection first, to ensure that the poly is not evaluated at the poles...
            DO WHILE(LOWPOLE.or.HIGHPOLE)
                IF(J.gt.2) THEN
                    IF(LOWPOLE) THEN
                        IF(((XHIGH-XLOW).lt.NRCONV).or.
     &                      (J.GT.NRSTEPSMAX)) GOTO 123
!                        WRITE(55,*) "LOW BOUND NOT MOVED"
                        X=XLOW+DELTA
!                        WRITE(55,*) "ATTEMPTED XLOW = ",X
                        CALL CALCSTARPOLYVALUE(N,A,B,X,VAL)
                        IF(VAL.gt.0.0_dp) THEN
                            XLOW=X
                            VALLOW=VAL
                            LOWPOLE=.false.
                        ELSE
                            DELTA=DELTA/10.0_dp
                            XHIGH=X
                            VALHIGH=VAL
                            HIGHPOLE=.false.
                        ENDIF
                    ENDIF
                    IF(HIGHPOLE) THEN
                        IF(((XHIGH-XLOW).lt.NRCONV).or.
     &                      (J.GT.NRSTEPSMAX)) GOTO 123
!                        WRITE(55,*)"HIGHBOUND NOT MOVED, XHIGH= ",XHIGH
!                        WRITE(55,*)"DELTA = ", DELTA
                        X=XHIGH-DELTA
!                        WRITE(55,*) "ATTEMPTED XHIGH = ",X
                        CALL CALCSTARPOLYVALUE(N,A,B,X,VAL)
!                        WRITE(55,*) "VALUE GIVEN IS ", VAL
                        IF(VAL.lt.0.0_dp) THEN
                            XHIGH=X
                            VALHIGH=VAL
                            HIGHPOLE=.false.
                        ELSE
                            DELTA=DELTA/10.0_dp
                            XLOW=X
                            VALLOW=VAL
                            LOWPOLE=.false.
                        ENDIF
                    ENDIF
                ELSE
!                    X=(XLOW+XHIGH)/2.0_dp
                    !Experience says root normally closer to higher pole, so do decasector search...
                    X=XHIGH-(0.01*(XHIGH-XLOW))
!                    WRITE(55,*) "Decasection X = ", X
                    CALL CALCSTARPOLYVALUE(N,A,B,X,VAL)
                    IF(VAL.GT.0.0_dp) THEN
                        XLOW=X
                        VALLOW=VAL
                        LOWPOLE=.false.
!                        WRITE(55,*) "TOO LOW"
                    ELSE
                        XHIGH=X
                        VALHIGH=VAL
                        HIGHPOLE=.false.
!                        WRITE(55,*) "TOO HIGH"
                    ENDIF
                ENDIF
                J=J+1
            ENDDO

            !If looking at the highest root, it should converge quickly by NR from the lowest bound,
            !and vice versa for the lowest root
            IF(I.eq.NROOTS) THEN
                X=XLOW
                VAL=VALLOW
            ELSEIF(I.eq.1) THEN
                X=XHIGH
                VAL=VALHIGH
            ENDIF
            
!            WRITE(55,*) "END OF BISECTION, XHIGH/XLOW MOVED BY", 
!     &            (INITHIGH-XHIGH),(XLOW-INITLOW)
!            WRITE(55,*) "XHIGH/XLOW NOW: ", XHIGH,XLOW
!             WRITE(55,*) "END OF BISECTION - ", J," STEPS NEEDED"
!            WRITE(55,*) "No. steps needed: ",J
!            CALL neci_flush(55)

            !No longer need to do extended RF since have finite boundary
!            WRITE(55,*) "***RF starting..."
!  Now, do some regula falsi, which has faster convergence than bisection, but still stable...
!            J=0
!            DO WHILE(ABS(VAL).GT.1.0e-9_dp.AND.ABS(XLOW-XHIGH).GT.RFCONV)
!                X=(-(VALHIGH*(XLOW-XHIGH))/(VALLOW-VALHIGH))+XHIGH
!                WRITE(55,*) J,VALHIGH,VALLOW,X
!                WRITE(55,*) XHIGH,XLOW

!                CALL CALCSTARPOLYVALUE(N,A,B,X,VAL)
!                J=J+1
!               WRITE(76,"(4G23.15)") XLOW,XHIGH,X,VAL
!               IF(VAL.GT.0.0_dp) THEN
!                  XLOW=X
!                  VALLOW=VAL
!               ELSEIF(VAL.LT.0.0_dp) THEN
!                  XHIGH=X
!                  VALHIGH=VAL
!               ELSE
!                  XLOW=X
!                  XHIGH=X
!               ENDIF
!            ENDDO
C.. Now polish off with some Newton Raphson
            J=0
            OX=0.0_dp
!            WRITE(76,*) "PreRoot found: ",X,A(NMIN+I),B(NMIN+I)
            !VAL is value of polynomial whose root we want,ie VAL=0 - X is root, OX old root

!       WRITE(55,*) "RF finished, giving limits of ",XHIGH," and ",XLOW
!       WRITE(55,*) "NR STARTING..."
            DO WHILE((ABS(VAL).GT.NRCONV).AND.(J.LT.NRSTEPSMAX)
     &      .AND.(ABS(X-OX).GT.NRCONV))
!                WRITE(55,*) I,J,VAL,X,OX,X-OX,DERIV
!                WRITE(55,*) "X = ",X
                CALL CALCSTARPOLYDERIV(N,A,B,X,DERIV)
               
                
                !If DERIV too small, then revert to regula falsi
                !Also if DERIV too large, z~=rho_jj and error large
                !Also errors if DERIV is inf
!          WRITE(55,*) "**DERIV TOO SMALL/LARGE - calculating using RF**"
!                K=0
                 IF((I.ne.1).and.(I.ne.NROOTS)) THEN
                     DO WHILE((ABS(DERIV).lt.1.0e-13_dp).or.
     &                       (ABS(DERIV).gt.1.D+18))
!                       WRITE(55,*) "DERIVATIVE UOB",K,DERIV
!                       WRITE(55,*) "VALHIGH/LOW = ",VALHIGH,VALLOW
                        OX=X
                        J=J+1
                        K=K+1
                      X=(-(VALHIGH*(XLOW-XHIGH))/(VALLOW-VALHIGH))+XHIGH
!                       X=XHIGH-(0.5*(XHIGH-XLOW))
!                       WRITE(55,*) "NEW X = ", X
                        CALL CALCSTARPOLYVALUE(N,A,B,X,VAL)
                        IF(VAL.GT.0.0_dp) THEN
                            XLOW=X
                            VALLOW=VAL
                        ENDIF
                        IF(VAL.LT.0.0_dp) THEN
                            XHIGH=X
                            VALHIGH=VAL
                        ENDIF
                        !test convergence
                        IF((ABS(VAL).LT.NRCONV).OR.(ABS(OX-X)
     &                   .LT.NRCONV).OR.(J.GT.NRSTEPSMAX)) THEN
                            GOTO 123
                        ENDIF
                        CALL CALCSTARPOLYDERIV(N,A,B,X,DERIV)
!                       WRITE(55,*) "RFNEW XHIGH/LOW= ",XHIGH,", ",XLOW
                     ENDDO
                 ENDIF
                      
                 !Perform NR calculation
!                 WRITE(55,*) "VAL, DERIV = ", VAL,DERIV
!                WRITE(76,"(4G23.15)") X,VAL,DERIV
                 OX=X
                 X=X-VAL/DERIV
                 CALL CALCSTARPOLYVALUE(N,A,B,X,VAL)
!                 WRITE(55,*) "NEW NR X = ", X
                 J=J+1
                 IF (((X.GT.XHIGH).or.(X.LT.XLOW)).and.
     &                   ((I.ne.NROOTS).and.(I.ne.1))) THEN
                     K=K+1
!  The NR has become unstable so we just use a regula Falsi step instead.
!                     WRITE(55,*) "*RFSTEP* - out of bounds"
!                     WRITE(55,*) "XHIGH/LOW = ", XHIGH, XLOW
                      X=(-(VALHIGH*(XLOW-XHIGH))/(VALLOW-VALHIGH))+XHIGH
!                      X=XHIGH-(0.5*(XHIGH-XLOW))
!                      WRITE(55,*) "NEW RF X = ", X
                      CALL CALCSTARPOLYVALUE(N,A,B,X,VAL)
!                      WRITE(55,*) "NEW VAL = ", VAL
                      IF(VAL.GT.0.0_dp) THEN
                          XLOW=X
                          VALLOW=VAL
                      ELSEIF(VAL.LT.0.0_dp) THEN
                          XHIGH=X
                          VALHIGH=VAL
                      ENDIF
                      !test convergence
                      IF((ABS(VAL).LT.NRCONV).OR.(ABS(OX-X)
     &                   .LT.NRCONV).OR.(J.GT.NRSTEPSMAX)) THEN
                          GOTO 123
                      ENDIF
                 
                 ELSE
                    IF(VAL.GT.0.0_dp) THEN
                        XLOW=X
                        VALLOW=VAL
                    ELSEIF(VAL.LT.0.0_dp) THEN
                        XHIGH=X
                        VALHIGH=VAL
                    ENDIF
                 ENDIF
            ENDDO


123         CONTINUE            
!            WRITE(55,*) "FINISHED NR WITH J= ", J, " and K= ",K
!            WRITE(55,*) ""
!            WRITE(55,*) "Number of steps with wrong deriv= ", K
!            CALL neci_flush(55)
!            JTOT=JTOT+J
            
!  If this STOP is reached then we've likely hit a pole with a very small numerator, so the root is very close to it.
!    Try increasing the value for NRSTEPSMAX.  Setting a RhoEpsilon might help 
            IF(J.GT.(NRSTEPSMAX-1)) THEN
                STOP "Newton Raphson failed to converge."
!                WRITE(55,*) "Newton Raphson failed to converge."
            ENDIF

            ROOTS(I)=X
!            WRITE(76,*) "Root found: ",X
            DO WHILE(I.GT.0.AND.
     &          ABS(REAL(A(NMIN+I),dp)-REAL(A(NMIN+I-1),dp))
     &          .LT.1.0e-9_dp)
               I=I-1
               ROOTS(I)=(A(NMIN+I))
C               WRITE(6,*),I,A(I),B(I)
!               WRITE(76,*) "Root found2:",A(I+NMIN),B(I+NMIN)
            ENDDO
            XHIGH=(A(NMIN+I))
            I=I-1
            IF(I.EQ.1) THEN
               XLOW=-1.0_dp
            ELSE
               XLOW=(A(NMIN+I))
            ENDIF
         ENDDO
!         WRITE(76,*) ROOTS(:)
         RETURN
      END

      SUBROUTINE CALCSTARPOLYVALUE(N,A,B,X,VAL)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER N
         HElement_t A(0:N),B(0:N)
         real(dp) X,VAL
         INTEGER I
         VAL=(A(0))-X
         DO I=1,N
            VAL=VAL+abs(B(I))**2/(X-(A(I)))
         ENDDO
         RETURN
      END
      SUBROUTINE CALCSTARPOLYDERIV(N,A,B,X,VAL)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER N
         HElement_t A(0:N),B(0:N)
         real(dp) X,VAL
         INTEGER I
         VAL=-1.0_dp
         DO I=1,N
            VAL=VAL-abs(B(I))**2/((X-(A(I)))**2)
         ENDDO
         RETURN
      END
!  This is used by the READSTARDIAG called from FMCPR3NVSTAR routine.  It allows contributions multiple times from some spokes.
!  It probably doesn't work (AJWT 6/6/07)
      SUBROUTINE FINDROOTSTAR2(N,A,B,NR,ROOTS,IRMIN)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER N,I,J
         real(dp) A(0:N),B(0:N),ROOTS(0:N)
         INTEGER NR(0:N)
         real(dp) XHIGH,XLOW,X,VAL,DERIV,XVLOW,XVHIGH,SMIN,OH,OX
         real(dp) OXLOW,OXHIGH
         INTEGER IMODE,IROOT,IRMIN
C.. We find the first root by starting very slightly above the highest
C.. pole
         XLOW=A(N)
         XHIGH=2*A(N)
C.. The value at XHIGH is always bigger than at XLOW
         I=N
         IROOT=N
         DO WHILE(I.GE.0)
C            WRITE(6,*),I,A(I),B(I)
            
            CALL CALCSTARPOLYVALUE2(N,A,B,XHIGH,NR,XVHIGH,.FALSE.)
            CALL CALCSTARPOLYVALUE2(N,A,B,XLOW,NR,XVLOW,.TRUE.)
C            WRITE(6,*) "HIGH:",XHIGH,VAL
C            WRITE(6,*) "LOW:",XLOW,DERIV
C            WRITE(6,*) VAL*DERIV
            VAL=XVHIGH*XVLOW
            IF(ABS(VAL).EQ.VAL) THEN
C.. No roots to be found here by regular falsi.
C.. However, if we find the minimum in the current section, we can try
C.. using that.
C               WRITE(6,*) I,"Falsi Fail:", XLOW,XHIGH,XVHIGH,XVLOW
               OH=XHIGH
               CALL FINDSTARPOLYMIN(N,A,B,XHIGH,XLOW,NR,SMIN)
               CALL CALCSTARPOLYVALUE2(N,A,B,SMIN,NR,XVLOW,.FALSE.)
C               WRITE(6,*) IROOT,"Min:",SMIN,XVLOW
               XLOW=SMIN
               XHIGH=OH
               VAL=XVLOW*XVHIGH
               IF(ABS(VAL).EQ.VAL) THEN
C.. Definitely no roots to be found.  give up
                  X=0.0_dp
                  IMODE=0
               ELSE
C.. There should be two roots
                  WRITE(6,*) "2 roots"
                  IMODE=2
               ENDIF
            ELSE
               IMODE=1
            ENDIF
            OXLOW=XLOW
            OXHIGH=XHIGH
            DO WHILE(IMODE.GT.0)
             VAL=1.0_dp
C            WRITE(6,*) "HIGH:", XHIGH,XVHIGH
C            WRITE(6,*) "LOW :", XLOW,XVLOW
             DO WHILE(ABS(VAL).GT.1.0e-9_dp.AND.ABS(XLOW-XHIGH)
     &                 .GT.1.0e-12_dp)
               X=(XLOW+XHIGH)/2.0_dp
               CALL CALCSTARPOLYVALUE2(N,A,B,X,NR,VAL,.FALSE.)
C               WRITE(6,"(4G19.12)") XLOW,XHIGH,X,VAL
               IF(VAL.GT.0) THEN
                  IF(XVLOW.GT.0) THEN
                     XLOW=X
                     XVLOW=VAL
                  ELSE
                     XHIGH=X
                     XVHIGH=VAL
                  ENDIF
               ELSEIF(VAL.LT.0) THEN
                  IF(XVHIGH.LT.0) THEN
                     XHIGH=X
                     XVHIGH=VAL
                  ELSE
                     XLOW=X
                     XVLOW=VAL
                  ENDIF
               ELSE
                  XLOW=X
                  XHIGH=X
               ENDIF
             ENDDO
C.. Now polish off with some Newton Raphson
            J=0
            IF(X.EQ.0.0_dp) J=11
            OX=0.0_dp
            DO WHILE(ABS(VAL).GT.3.0e-13_dp.AND.J.LT.10
     &         .AND.ABS(X-OX).GT.1.0e-13_dp)
               CALL CALCSTARPOLYVALUE2(N,A,B,X,NR,VAL,.FALSE.)
               CALL CALCSTARPOLYDERIV2(N,A,B,X,NR,DERIV)
C               WRITE(6,*) X,VAL
               OX=X
               X=X-VAL/DERIV
               J=J+1
            ENDDO
            IF(X.LT.OXLOW.OR.X.GT.OXHIGH) J=11
            IF(J.LT.10) THEN
C             WRITE(6,*) "Root found: ",IROOT,X
             ROOTS(IROOT)=X
             IROOT=IROOT-1
            ELSE
               WRITE(6,*) A(I+2),B(I+2),NR(I+2)
               WRITE(6,*) A(I+1),B(I+1),NR(I+1)
               WRITE(6,*) A(I),B(I),NR(I)
               WRITE(6,*) A(I-1),B(I-1),NR(I-1)
               WRITE(6,*) A(I-2),B(I-2),NR(I-2)
               WRITE(6,*) X,OXLOW,OXHIGH
               WRITE(6,*) I,"NO ROOT FOUND"
               CALL CALCSTARPOLYVALUE2(N,A,B,OXHIGH,NR,XVHIGH,.FALSE.)
               CALL CALCSTARPOLYVALUE2(N,A,B,OXLOW,NR,XVLOW,.TRUE.)
               WRITE(6,*) VAL,XVLOW,XVHIGH
C               STOP "NO ROOT FOUND"
            ENDIF
             IMODE=IMODE-1
             IF(IMODE.EQ.1) THEN
C.. one more root to find
               XHIGH=SMIN
               XLOW=A(I)
               CALL CALCSTARPOLYVALUE2(N,A,B,XHIGH,NR,XVHIGH,.FALSE.)
               CALL CALCSTARPOLYVALUE2(N,A,B,XLOW,NR,XVLOW,.TRUE.)
             ENDIF
            ENDDO
             DO WHILE(I.GT.0.AND.ABS(A(I)-A(I-1)).LT.1.0e-9_dp)
               I=I-1
               ROOTS(IROOT)=0.0_dp
C..A(I)
C               WRITE(6,*) "Root found2:",IROOT,A(I)
C               IROOT=IROOT-1
C               WRITE(6,*),I,A(I),B(I)
             ENDDO
            XHIGH=A(I)
            I=I-1
            IF(I.EQ.0) THEN
               XLOW=-1.0_dp
            ELSE
               XLOW=A(I)
            ENDIF
         ENDDO
         DO I=IROOT,0,-1
            ROOTS(I)=0.0_dp
         ENDDO
         IRMIN=IROOT+1
         RETURN
      END         
         

      SUBROUTINE CALCSTARPOLYVALUE2(N,A,B,X,NR,VAL,TABOVE)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER N
         real(dp) A(0:N),B(0:N),X,VAL,OV,X2,DX
         INTEGER NR(0:N)
         INTEGER I
         LOGICAL TABOVE
         VAL=A(0)-X
         OV=0.0_dp
         DO I=1,N
            VAL=VAL+NR(I)*B(I)*B(I)/(X-A(I))
            IF(X.EQ.A(I)) OV=OV+NR(I)*B(I)*B(I)
         ENDDO
         IF(OV.NE.0.0_dp) THEN 
C.. We're actually on a root.
C. if above, it then x-a >0, so just the sign of OV
C            VAL=1.D10*OV   
C            IF(.NOT.TABOVE) VAL=-VAL
            IF(TABOVE) THEN
               DX=1.0e-16_dp
            ELSE
               DX=-1.0e-16_dp
            ENDIF
            X2=X+DX
            DO WHILE (X2.EQ.X)
               DX=DX*2.0_dp
               X2=X+DX
            ENDDO
            VAL=A(0)-X2
            DO I=1,N
               VAL=VAL+NR(I)*B(I)*B(I)/(X2-A(I))
            ENDDO
         ENDIF
         RETURN
      END
      SUBROUTINE CALCSTARPOLYDERIV2(N,A,B,X,NR,VAL)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER N
         real(dp) A(0:N),B(0:N),X,VAL
         INTEGER NR(0:N)
         INTEGER I
         VAL=-1.0_dp
         DO I=1,N
            VAL=VAL-NR(I)*B(I)*B(I)/((X-A(I))**2)
         ENDDO
         RETURN
      END
      SUBROUTINE CALCSTARPOLYDERIV3(N,A,B,X,NR,VAL)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER N
         real(dp) A(0:N),B(0:N),X,VAL
         INTEGER NR(0:N)
         INTEGER I
         VAL=0.0_dp
         DO I=1,N
            VAL=VAL+2*NR(I)*B(I)*B(I)/((X-A(I))**3)
         ENDDO
         RETURN
      END
      SUBROUTINE FINDSTARPOLYMIN(N,A,B,XH,XL,NR,ROOT)
         use constants, only: dp
         IMPLICIT NONE
         INTEGER N,J
         real(dp) A(0:N),B(0:N)
         INTEGER NR(0:N)
         real(dp) XHIGH,XLOW,X,VAL,DERIV,XVLOW,XVHIGH,ROOT,OX
         real(dp) XH,XL
         XHIGH=XH
         XLOW=XL
         CALL CALCSTARPOLYDERIV2(N,A,B,XHIGH*0.99999999_dp,NR,XVHIGH)
         CALL CALCSTARPOLYDERIV2(N,A,B,XLOW *1.00000001_dp,NR,XVLOW)
         VAL=1.0_dp
         DO WHILE(ABS(VAL).GT.1.0e-9_dp.AND.ABS(XLOW-XHIGH)
     &            .GT.1.0e-12_dp)
            X=(XLOW+XHIGH)/2.0_dp
            CALL CALCSTARPOLYDERIV2(N,A,B,X,NR,VAL)
            IF(VAL.GT.0) THEN
               IF(XVLOW.GT.0) THEN
                  XLOW=X
                  XVLOW=VAL
               ELSE
                  XHIGH=X
                  XVHIGH=VAL
               ENDIF
            ELSEIF(VAL.LT.0) THEN
               IF(XVHIGH.LT.0) THEN
                  XHIGH=X
                  XVHIGH=VAL
               ELSE
                  XLOW=X
                  XVLOW=VAL
               ENDIF
            ELSE
               XLOW=X
               XHIGH=X
            ENDIF
         ENDDO
C.. Now polish off with some Newton Raphson
         J=0
         OX=0.0_dp
         DO WHILE(ABS(VAL).GT.1.0e-13_dp.AND.J.LT.10
     &      .AND.ABS(X-OX).GT.1.0e-13_dp)
            CALL CALCSTARPOLYDERIV2(N,A,B,X,NR,VAL)
            CALL CALCSTARPOLYDERIV3(N,A,B,X,NR,DERIV)
            OX=X
            X=X-VAL/DERIV
            J=J+1
         ENDDO
         ROOT=X
         IF(X.LT.XL.OR.X.GT.XH) ROOT=(XHIGH+XLOW)/2.0_dp
         RETURN
      END
