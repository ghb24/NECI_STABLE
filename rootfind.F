C.. A procedure to use the regular falsi (false descent) method
C.. to find the roots of a polynomial in the form

C.. a0-x + b1^2/(x-a1) +b2^2/(x-a2) + ...

C.. which is the characteristic polynomial of a matrix corresponding to
C.. a star-graph, with a as the diagonal terms, and b as the off
C.. diagonal (first row) terms.  a0 is the diagonal of the pivot of the
C.. star.  THere are N spokes and 1 pivot, and thus N+1 roots

C.. A and B are sorted in ascending order of A

      SUBROUTINE FINDROOTSTAR(N,A,B,ROOTS,NROOTS)
         USE HElement
         IMPLICIT NONE
         INTEGER N,I,J,NROOTS,NMIN
         TYPE(HElement) A(0:N),B(0:N)
         real*8 ROOTS(0:NROOTS)
         REAL*8 XHIGH,XLOW,X,VAL,DERIV,OX!,VALLOW,VALHIGH
         logical tonestep
         include 'calcp.inc'
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
!  This will be bounded from below by A(I+NMIN)
         XLOW=DREAL(A(I+NMIN))
         XHIGH=2*DREAL(A(I+NMIN))
C.. The value at XHIGH is always bigger than at XLOW

!  When I=0, there is no min.
!I points to the root we're finding, starting from the highest, NROOTS, and working down.
         DO WHILE(I.GT.0)
            
!            WRITE(76,*),I,A(I+NMIN),B(I+NMIN)
            VAL=1.D0
!            CALL CALCSTARPOLYVALUE(N,A,B,XLOW,VALLOW)
!            CALL CALCSTARPOLYVALUE(N,A,B,XHIGH,VALHIGH)
!            WRITE(55,*) "*********NEXT ROOT******** - RF starting..."
!            WRITE(55,*),I,A(I+NMIN),B(I+NMIN)
!  Do some regular falsi
!  Ensure we always do at least 1 step irrespective of the convergence thresholds - this gives us  and initial position X, and value VAL.
            tOneStep=.true.
            DO WHILE(tOneStep.or.
     &       (ABS(VAL).GT.1.D-9.AND.ABS(XLOW-XHIGH).GT.RFCONV))
               tOneStep=.false.
!                X=(-(VALHIGH*(XLOW-XHIGH))/(VALLOW-VALHIGH))+XHIGH

                X=(XLOW+XHIGH)/2.D0
                CALL CALCSTARPOLYVALUE(N,A,B,X,VAL)
!               WRITE(76,"(4G23.15)") XLOW,XHIGH,X,VAL
               IF(VAL.GT.0.D0) THEN
                  XLOW=X
!                  VALLOW=VAL
               ELSEIF(VAL.LT.0.D0) THEN
                  XHIGH=X
!                  VALHIGH=VAL
               ELSE
                  XLOW=X
                  XHIGH=X
               ENDIF
            ENDDO
C.. Now polish off with some Newton Raphson
            J=0
            OX=0.D0
!            WRITE(76,*) "PreRoot found: ",X,A(NMIN+I),B(NMIN+I)
            !VAL is value of polynomial whose root we want,ie VAL=0 - X is root, OX old root

!       WRITE(55,*) "RF finished, giving limits of ",XHIGH," and ",XLOW
!       WRITE(55,*) "NR STARTING..."
            DO WHILE(ABS(VAL).GT.NRCONV.AND.(J.LT.NRSTEPSMAX)
     &      .AND.(ABS(X-OX).GT.NRCONV))
!                WRITE(55,*) I,J,VAL,X,OX,X-OX,DERIV
                
                CALL CALCSTARPOLYDERIV(N,A,B,X,DERIV)
                
                !If DERIV too small, then revert to regular falsi
                !Also if DERIV too large, z~=rho_jj and error large
                !Also errors if DERIV is inf
                   IF((ABS(DERIV).lt.1.D-13).or.
     &                  (ABS(DERIV).gt.1.D+08)) THEN
!          WRITE(55,*) "**DERIV TOO SMALL/LARGE - calculating using RF**"
                      OX=X
                      J=J+1
                      IF(X.GT.XHIGH) THEN
                          X=(OX+XHIGH)/2
                      ENDIF
                      IF(X.LT.XLOW) THEN
                          X=(OX+XLOW)/2
                      ENDIF
!                      X=(-(VALHIGH*(XLOW-XHIGH))/(VALLOW-VALHIGH))+XHIGH
                      CALL CALCSTARPOLYVALUE(N,A,B,X,VAL)
                      IF(VAL.GT.0.D0) THEN
                          XLOW=X
!                          VALLOW=VAL
                      ELSEIF(VAL.LT.0.D0) THEN
                          XHIGH=X
!                          VALHIGH=VAL
                      ELSE
                          XLOW=X
                          XHIGH=X
                      ENDIF
!                      WRITE(55,*) "NEW XHIGH/LOW= ",XHIGH,", ",XLOW
                   ELSE
                       !Perform NR calculation
                       CALL CALCSTARPOLYVALUE(N,A,B,X,VAL)
!                      WRITE(76,"(4G23.15)") X,VAL,DERIV
                       OX=X
                       X=X-VAL/DERIV
                       J=J+1
                       IF((X.GT.XHIGH).or.(X.LT.XLOW)) THEN
!  The NR has become unstable so we just use a regular Falsi step instead.
!                          WRITE(55,*) "*RFSTEP* - out of bounds"
                            IF(X.GT.XHIGH) THEN
                                X=(OX+XHIGH)/2
                            ENDIF
                            IF(X.LT.XLOW) THEN
                                X=(OX+XLOW)/2
                            ENDIF

!                          X=(-(VALHIGH*(XLOW-XHIGH))/
!     &                       (VALLOW-VALHIGH))+XHIGH
                          CALL CALCSTARPOLYVALUE(N,A,B,X,VAL)
                          IF(VAL.GT.0.D0) THEN
                              XLOW=X
!                              VALLOW=VAL
                          ELSEIF(VAL.LT.0.D0) THEN
                              XHIGH=X
!                              VALHIGH=VAL
                          ELSE
                              XLOW=X
                              XHIGH=X
                          ENDIF
!                          WRITE(55,*) "NEW XHIGH/LOW= ",XHIGH,", ",XLOW
                          
!                          X=(OX+XHIGH)/2
                       ENDIF
!                       IF(X.LT.XLOW) THEN
!                           WRITE(55,*) "*RFSTEP* - too low"
!!  The NR has become unstable so we just use a regular Falsi step instead.
!                           X=(OX+XLOW)/2
!                       ENDIF
                   ENDIF
            ENDDO

!  If this STOP is reached then we've likely hit a pole with a very small numerator, so the root is very close to it.
!    Try increasing the value for NRSTEPSMAX.  Setting a RhoEpsilon might help 
            IF(J.GT.(NRSTEPSMAX-1)) THEN
                STOP "Newton Raphson failed to converge."
!                WRITE(55,*) "Newton Raphson failed to converge."
            ENDIF

            ROOTS(I)=X
!            WRITE(76,*) "Root found: ",X
            DO WHILE(I.GT.0.AND.
     &          ABS(DREAL(A(NMIN+I))-DREAL(A(NMIN+I-1))).LT.1.D-9)
               I=I-1
               ROOTS(I)=DREAL(A(NMIN+I))
C               WRITE(6,*),I,A(I),B(I)
!               WRITE(76,*) "Root found2:",A(I+NMIN),B(I+NMIN)
            ENDDO
            XHIGH=DREAL(A(NMIN+I))
            I=I-1
            IF(I.EQ.1) THEN
               XLOW=-1.D0
            ELSE
               XLOW=DREAL(A(NMIN+I))
            ENDIF
         ENDDO
         RETURN
      END         
         

      SUBROUTINE CALCSTARPOLYVALUE(N,A,B,X,VAL)
         USE HElement
         IMPLICIT NONE
         INTEGER N
         TYPE(HElement) A(0:N),B(0:N)
         REAL*8 X,VAL
         INTEGER I
         VAL=DREAL(A(0))-X
         DO I=1,N
            VAL=VAL+SQ(B(I))/(X-DREAL(A(I)))
         ENDDO
         RETURN
      END
      SUBROUTINE CALCSTARPOLYDERIV(N,A,B,X,VAL)
         USE HElement
         IMPLICIT NONE
         INTEGER N
         TYPE(HElement) A(0:N),B(0:N)
         REAL*8 X,VAL
         INTEGER I
         VAL=-1.D0
         DO I=1,N
            VAL=VAL-SQ(B(I))/((X-DREAL(A(I)))**2)
         ENDDO
         RETURN
      END
!  This is used by the READSTARDIAG called from FMCPR3NVSTAR routine.  It allows contributions multiple times from some spokes.
!  It probably doesn't work (AJWT 6/6/07)
      SUBROUTINE FINDROOTSTAR2(N,A,B,NR,ROOTS,IRMIN)
         IMPLICIT NONE
         INTEGER N,I,J
         REAL*8 A(0:N),B(0:N),ROOTS(0:N)
         INTEGER NR(0:N)
         REAL*8 XHIGH,XLOW,X,VAL,DERIV,XVLOW,XVHIGH,SMIN,OH,OX
         REAL*8 OXLOW,OXHIGH
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
                  X=0.D0
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
             VAL=1.D0
C            WRITE(6,*) "HIGH:", XHIGH,XVHIGH
C            WRITE(6,*) "LOW :", XLOW,XVLOW
             DO WHILE(ABS(VAL).GT.1.D-9.AND.ABS(XLOW-XHIGH).GT.1.D-12)
               X=(XLOW+XHIGH)/2.D0
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
            IF(X.EQ.0.D0) J=11
            OX=0.D0
            DO WHILE(ABS(VAL).GT.3.D-13.AND.J.LT.10
     &         .AND.ABS(X-OX).GT.1.D-13)
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
             DO WHILE(I.GT.0.AND.ABS(A(I)-A(I-1)).LT.1.D-9)
               I=I-1
               ROOTS(IROOT)=0.D0
C..A(I)
C               WRITE(6,*) "Root found2:",IROOT,A(I)
C               IROOT=IROOT-1
C               WRITE(6,*),I,A(I),B(I)
             ENDDO
            XHIGH=A(I)
            I=I-1
            IF(I.EQ.0) THEN
               XLOW=-1.D0
            ELSE
               XLOW=A(I)
            ENDIF
         ENDDO
         DO I=IROOT,0,-1
            ROOTS(I)=0.D0
         ENDDO
         IRMIN=IROOT+1
         RETURN
      END         
         

      SUBROUTINE CALCSTARPOLYVALUE2(N,A,B,X,NR,VAL,TABOVE)
         IMPLICIT NONE
         INTEGER N
         REAL*8 A(0:N),B(0:N),X,VAL,OV,X2,DX
         INTEGER NR(0:N)
         INTEGER I
         LOGICAL TABOVE
         VAL=A(0)-X
         OV=0.D0
         DO I=1,N
            VAL=VAL+NR(I)*B(I)*B(I)/(X-A(I))
            IF(X.EQ.A(I)) OV=OV+NR(I)*B(I)*B(I)
         ENDDO
         IF(OV.NE.0.D0) THEN 
C.. We're actually on a root.
C. if above, it then x-a >0, so just the sign of OV
C            VAL=1.D10*OV   
C            IF(.NOT.TABOVE) VAL=-VAL
            IF(TABOVE) THEN
               DX=1.D-16
            ELSE
               DX=-1.D-16
            ENDIF
            X2=X+DX
            DO WHILE (X2.EQ.X)
               DX=DX*2.D0
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
         IMPLICIT NONE
         INTEGER N
         REAL*8 A(0:N),B(0:N),X,VAL
         INTEGER NR(0:N)
         INTEGER I
         VAL=-1.D0
         DO I=1,N
            VAL=VAL-NR(I)*B(I)*B(I)/((X-A(I))**2)
         ENDDO
         RETURN
      END
      SUBROUTINE CALCSTARPOLYDERIV3(N,A,B,X,NR,VAL)
         IMPLICIT NONE
         INTEGER N
         REAL*8 A(0:N),B(0:N),X,VAL
         INTEGER NR(0:N)
         INTEGER I
         VAL=0.D0
         DO I=1,N
            VAL=VAL+2*NR(I)*B(I)*B(I)/((X-A(I))**3)
         ENDDO
         RETURN
      END
      SUBROUTINE FINDSTARPOLYMIN(N,A,B,XH,XL,NR,ROOT)
         IMPLICIT NONE
         INTEGER N,I,J
         REAL*8 A(0:N),B(0:N)
         INTEGER NR(0:N)
         REAL*8 XHIGH,XLOW,X,VAL,DERIV,XVLOW,XVHIGH,ROOT,OX
         REAL*8 XH,XL
         XHIGH=XH
         XLOW=XL
         CALL CALCSTARPOLYDERIV2(N,A,B,XHIGH*0.99999999D0,NR,XVHIGH)
         CALL CALCSTARPOLYDERIV2(N,A,B,XLOW *1.00000001D0,NR,XVLOW)
         VAL=1.D0
         DO WHILE(ABS(VAL).GT.1.D-9.AND.ABS(XLOW-XHIGH).GT.1.D-12)
            X=(XLOW+XHIGH)/2.D0
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
         OX=0.D0
         DO WHILE(ABS(VAL).GT.1.D-13.AND.J.LT.10
     &      .AND.ABS(X-OX).GT.1.D-13)
            CALL CALCSTARPOLYDERIV2(N,A,B,X,NR,VAL)
            CALL CALCSTARPOLYDERIV3(N,A,B,X,NR,DERIV)
            OX=X
            X=X-VAL/DERIV
            J=J+1
         ENDDO
         ROOT=X
         IF(X.LT.XL.OR.X.GT.XH) ROOT=(XHIGH+XLOW)/2.D0
         RETURN
      END
