C.. This works for large P (may need to take more terms - 8/12/03)
C.. approximate w_h(1+e_1,1+e_2,...,1+e_h,P-h) = (P-h h)
C..    + (P-h h+1)Sum_i=1,h(e_i) + (P-h h+2) Sum_(i=1,h;j=1,h)(e_i e_j) + ...
C.. = (P h) (1+((P-h)/(h+1))Sum_(i=1,h)(e_i (1+ ((P-h-1)/(h+2))Sum_(j=i,h) e_j (1+ ...)
C.. return (RH_II^-I_H)*RH_IJ*RH_JK*...*RH_LI*W_h(RH_JJ/RH_II,...RH_LL/RH_II,P)
C.. This is recursive and will generate up to order NWHTAY
      REAL*8 RECURSIVE FUNCTION CALCPATH6_2(IPATH,RHOII,RHOCUML,I_H,I_P,
     &      NWHTAY,ISTART,NORDER) RESULT (CALCPATHS6_2RES)
         REAL*8 RHOII(0:I_H),DFACT,CALCPATHS6_2RES,FNCR
         INTEGER NWHTAY
         F1=0.D0
         DFACT=(I_P-I_H-NORDER)/(I_H+NORDER+1.D0)
         DO I=ISTART,I_H-1
            IF(NORDER.LT.(NWHTAY-1)) THEN
              F2=CALCPATH6_2(IPATH,RHOII,RHOCUML,I_H,I_P,
     &               NWHTAY,I,NORDER+1)         
            ELSE
               F2=1.D0
            ENDIF
            F1=F1+(RHOII(I)/RHOII(0)-1)*F2
         ENDDO
         F1=1.D0+F1*DFACT
         IF(NORDER.EQ.0) THEN
            CALCPATHS6_2RES=FNCR(I_P,I_H)*F1
C.CALCPATHS6_2 now contains the weight
            CALCPATHS6_2RES=CALCPATHS6_2RES*(RHOII(0)**(-I_H))*RHOCUML
         ELSE
            CALCPATHS6_2RES=F1
         ENDIF
         RETURN
      END

C.. Generate a real that is (P H) (the binomial coefficient)
C.. Best for small H
      REAL*8 FUNCTION FNCR(I_P,I_H)
         FNCR=1.D0
         DO I=0,I_H-1
            FNCR=FNCR*(I_P-I)/(I+1)
         ENDDO
         RETURN
      END

      LOGICAL FUNCTION LSAME2(A,B)
         IMPLICIT REAL*8(A-H,O-Z)
         IMPLICIT LOGICAL(L)
         LSAME2=.FALSE.
         IF((ABS((A/B)-1.D0).LT.1.D-10).OR.
     &      ABS(A-B).LT.1.D-10) LSAME2=.TRUE.
C         IF(LSAME2) WRITE(10,*) "LSAME2",A,B,LSAME2
         RETURN
      END
C.........................................
C.       CALCPATH8 not used or finished
C.........................................
      REAL*8 FUNCTION CALCPATH8(IPATH,RHOII,RHOCUML,I_HH,
     &               I_P,NWHTAY)
         IMPLICIT REAL*8 (A-H,O-Z)
         IMPLICIT LOGICAL (L)
         INTEGER IPATH(0:I_HH)
         REAL*8 RHOII(0:I_HH-1)
         REAL*8 S(I_HH-1)
         INTEGER IGS(I_HH-1),IG
         LOGICAL LFLIP

         CALCPATH8=0
         I_H=I_HH
         IINCR=0
         IF(NWHTAY.GT.0) THEN
           CALCPATH7=CALCPATH6_2(IPATH,RHOII,RHOCUML,I_H,I_P,NWHTAY,1,0)
            RETURN
         ENDIF
C.. We need to be careful if RHOII(j)=0, in which case we move from 
C.. Z_h(ijk...l,p-h) to Z_(h-1)(ik...l,p-h)
         DO I=1,I_HH-1
            IF(RHOII(I).EQ.0.D0) THEN
               IINCR=IINCR+1
               I_H=I_H-1
            ELSE
               S(I-IINCR)=RHOII(I)/RHOII(0)-1.D0
            ENDIF
         ENDDO
         IF(I_H.LT.I_HH) THEN
            I=1
         ENDIF
C.. To determine which formula to use to take into account degenerecies,
C.. we need to sort S1,S2,S3.  We bubble sort (as I_H will never be too large)
         LFLIP=.TRUE.
         DO WHILE (LFLIP)
            IG=1
            IGS(1)=1
            LFLIP=.FALSE.
            DO I=1,I_H-2
               IF(S(I).NE.S(I+1)) THEN
                  IG=IG+1
                  IGS(IG)=1
               ELSE
                  IGS(IG)=IGS(IG)+1
               ENDIF
               IF(S(I).GT.S(I+1)) THEN
                  ST=S(I)
                  S(I)=S(I+1)
                  S(I+1)=ST
                  LFLIP=.TRUE.
               ENDIF 
            ENDDO
         ENDDO
C.. Here we need to sort the now grouped degenerate values into order of
C.. group size.  The group size was worked out on the last pass of the
C.. bubble sort
         IONE=0
         DO I=1,I_H-1
            IF(LSAME2(S(I),0.D0)) IONE=IONE+1
         ENDDO
C.. Now we go to specifics for each I_H
         IF(I_H.EQ.4) THEN
            IIP=I_P
            IF(LSAME2(S(1),S(2))) THEN
C..      XXY
               IF(LSAME2(S(2),S(3))) THEN
C..      XXX
                  IF(IONE.EQ.3) THEN
C..  111
                     CALCPATH7=WW4111(IIP)
                  ELSE
C..  AAA
                     CALCPATH7=WW4AAA(S(1),IIP)
                  ENDIF
               ELSE
C..      XXY
                  IF(IONE.EQ.2) THEN
C..  11A
                     CALCPATH7=WW4A11(S(3),IIP)
                  ELSEIF(IONE.EQ.1) THEN
C..  AA1
                     CALCPATH7=WW4AA1(S(1),IIP)
                  ELSE
C..  AAB
                     CALCPATH7=WW4AAB(S(1),S(3),IIP)
                  ENDIF
               ENDIF
            ELSEIF(LSAME2(S(2),S(3))) THEN
C..     XYY
               IF(IONE.EQ.1) THEN
C..  1AA       
                  CALCPATH7=WW4AA1(S(2),IIP)
               ELSEIF(IONE.EQ.2) THEN
C..  A11
                  CALCPATH7=WW4A11(S(1),IIP)
               ELSE
C..  ABB
                  CALCPATH7=WW4AAB(S(1),S(3),IIP)
               ENDIF
            ELSEIF(IONE.EQ.1) THEN
C..     XYZ
               IF(LSAME2(S(1),0.D0)) THEN
C..  1AB
                  CALCPATH7=WW4AB1(S(2),S(3),IIP)
               ELSEIF(LSAME2(S(2),0.D0)) THEN
C..  A1B
                  CALCPATH7=WW4AB1(S(1),S(3),IIP)
               ELSE
C..  AB1
                  CALCPATH7=WW4AB1(S(1),S(2),IIP)
               ENDIF
            ELSE
C..  ABC
               CALCPATH7=WW4ABC(S(1),S(2),S(3),IIP)
            ENDIF           

            CALCPATH7=(RHOII(0)**(-I_HH))*RHOCUML*CALCPATH7
 
         ELSEIF(I_H.EQ.3) THEN
            IIP=I_P
            IF(LSAME2(S(1),S(2))) THEN
               IF(LSAME2(S(1),0.D0)) THEN
                  CALCPATH7=(IIP*(2 - 3*IIP + IIP**2))/6.
               ELSE
                  CALCPATH7=
     &  (2 - 2*(1 + S(1))**IIP + S(1)*(1 + (1 + S(1))**(-1 + IIP))*IIP)/
     -   S(1)**3


               ENDIF  
            ELSEIF(LSAME2(S(1),0.D0).OR.LSAME2(S(2),0.D0)) THEN
               IF(LSAME2(S(2),0.D0)) S(2)=S(1)
               CALCPATH7=
     & -((1 - (1 + S(2))**IIP + S(2)*IIP + (S(2)**2*(-1 + IIP)*IIP)/2.)/
     -     S(2)**3)


            ELSE
               CALCPATH7=
     &    ((-1 + (1 + S(1))**IIP)*S(2)**2 - S(1)*S(2)**2*IIP +
     -     S(1)**2*(1 - (1 + S(2))**IIP + S(2)*IIP))/
     -   (S(1)**2*(S(1) - S(2))*S(2)**2)

            ENDIF
            CALCPATH7=(RHOII(0)**(-I_HH))*RHOCUML*CALCPATH7

         ELSEIF(I_H.EQ.2) THEN
            IF(LSAME2(S(1),0.D0)) THEN
               CALCPATH7=(RHOII(0)**(-I_HH))*RHOCUML
     &            *I_P*(I_P-1)/2.D0
            ELSE
               CALCPATH7=(RHOII(0)**(-I_HH))*RHOCUML
     &      *(-1 -I_P*S(1) + (1+S(1))**I_P)/ (S(1)**2)
            ENDIF
         ELSE
            CALCPATH7=0
C            STOP 'I_H > 4'
         ENDIF
         RETURN
      END


      REAL*8 FUNCTION CALCPATH7(IPATH,RHOII,RHOCUML,I_HH,
     &               I_P,NWHTAY)
         IMPLICIT REAL*8 (A-H,O-Z)
         IMPLICIT LOGICAL (L)
         INTEGER IPATH(0:I_HH)
         REAL*8 RHOII(0:I_HH-1)
         REAL*8 S(0:I_HH-1),SUM,SUMI
         LOGICAL LFLIP
         I_H=I_HH
         IINCR=0
         IF(NWHTAY.GT.0) THEN
           CALCPATH7=CALCPATH6_2(IPATH,RHOII,RHOCUML,I_H,I_P,NWHTAY,1,0)
            RETURN
         ENDIF
         IF(NWHTAY.EQ.-1) THEN
C.. We use contour integration instead
C.. We first have to divide all our RHOIIS by RHOII(0).  We do this and calculate the radius of the contour integration as slightly larger.
            S(0)=1.D0
            SMAX=1.D0
            DO I=1,I_HH-1
               S(I)=RHOII(I)/RHOII(0)
               IF(S(I).GT.SMAX) SMAX=S(I)
            ENDDO
C.. Try a radius 10% larger than the max S
            CALL ZSUM(S(0:I_HH-1),I_HH,I_P,1500,SMAX*1.1D0,SUM,SUMI)
            CALCPATH7=SUM
            CALCPATH7=(RHOII(0)**(-I_HH))*RHOCUML*CALCPATH7
            RETURN
         ENDIF
C.. We need to be careful if RHOII(j)=0, in which case we move from 
C.. Z_h(ijk...l,p-h) to Z_(h-1)(ik...l,p-h)
         DO I=1,I_HH-1
            IF(RHOII(I).EQ.0.D0) THEN
               IINCR=IINCR+1
               I_H=I_H-1
            ELSE
               S(I-IINCR)=RHOII(I)/RHOII(0)-1.D0
            ENDIF
         ENDDO
         IF(I_H.LT.I_HH) THEN
            I=1
         ENDIF
C.. To determine which formula to use to take into account degenerecies,
C.. we need to sort S1,S2,S3.  We bubble sort (as I_H will never be too large)
         LFLIP=.TRUE.
         DO WHILE (LFLIP)
            LFLIP=.FALSE.
            DO I=1,I_H-2
               IF(S(I).GT.S(I+1)) THEN
                  ST=S(I)
                  S(I)=S(I+1)
                  S(I+1)=ST
                  LFLIP=.TRUE.
               ENDIF 
            ENDDO
         ENDDO
         IONE=0
         DO I=1,I_H-1
            IF(LSAME2(S(I),0.D0)) IONE=IONE+1
         ENDDO
C.. Now we go to specifics for each I_H
         IF(I_H.EQ.4) THEN
            IIP=I_P
            IF(LSAME2(S(1),S(2))) THEN
C..      XXY
               IF(LSAME2(S(2),S(3))) THEN
C..      XXX
                  IF(IONE.EQ.3) THEN
C..  111
                     CALCPATH7=WW4111(IIP)
                  ELSE
C..  AAA
                     CALCPATH7=WW4AAA(S(1),IIP)
                  ENDIF
               ELSE
C..      XXY
                  IF(IONE.EQ.2) THEN
C..  11A
                     CALCPATH7=WW4A11(S(3),IIP)
                  ELSEIF(IONE.EQ.1) THEN
C..  AA1
                     CALCPATH7=WW4AA1(S(1),IIP)
                  ELSE
C..  AAB
                     CALCPATH7=WW4AAB(S(1),S(3),IIP)
                  ENDIF
               ENDIF
            ELSEIF(LSAME2(S(2),S(3))) THEN
C..     XYY
               IF(IONE.EQ.1) THEN
C..  1AA       
                  CALCPATH7=WW4AA1(S(2),IIP)
               ELSEIF(IONE.EQ.2) THEN
C..  A11
                  CALCPATH7=WW4A11(S(1),IIP)
               ELSE
C..  ABB
                  CALCPATH7=WW4AAB(S(1),S(3),IIP)
               ENDIF
            ELSEIF(IONE.EQ.1) THEN
C..     XYZ
               IF(LSAME2(S(1),0.D0)) THEN
C..  1AB
                  CALCPATH7=WW4AB1(S(2),S(3),IIP)
               ELSEIF(LSAME2(S(2),0.D0)) THEN
C..  A1B
                  CALCPATH7=WW4AB1(S(1),S(3),IIP)
               ELSE
C..  AB1
                  CALCPATH7=WW4AB1(S(1),S(2),IIP)
               ENDIF
            ELSE
C..  ABC
               CALCPATH7=WW4ABC(S(1),S(2),S(3),IIP)
            ENDIF           

            CALCPATH7=(RHOII(0)**(-I_HH))*RHOCUML*CALCPATH7
 
         ELSEIF(I_H.EQ.3) THEN
            IIP=I_P
            IF(LSAME2(S(1),S(2))) THEN
               IF(LSAME2(S(1),0.D0)) THEN
                  CALCPATH7=(IIP*(2 - 3*IIP + IIP**2))/6.
               ELSE
                  CALCPATH7=
     &  (2 - 2*(1 + S(1))**IIP + S(1)*(1 + (1 + S(1))**(-1 + IIP))*IIP)/
     -   S(1)**3


               ENDIF  
            ELSEIF(LSAME2(S(1),0.D0).OR.LSAME2(S(2),0.D0)) THEN
               IF(LSAME2(S(2),0.D0)) S(2)=S(1)
               CALCPATH7=
     & -((1 - (1 + S(2))**IIP + S(2)*IIP + (S(2)**2*(-1 + IIP)*IIP)/2.)/
     -     S(2)**3)


            ELSE
               CALCPATH7=
     &    ((-1 + (1 + S(1))**IIP)*S(2)**2 - S(1)*S(2)**2*IIP +
     -     S(1)**2*(1 - (1 + S(2))**IIP + S(2)*IIP))/
     -   (S(1)**2*(S(1) - S(2))*S(2)**2)

            ENDIF
            CALCPATH7=(RHOII(0)**(-I_HH))*RHOCUML*CALCPATH7

         ELSEIF(I_H.EQ.2) THEN
            IF(LSAME2(S(1),0.D0)) THEN
               CALCPATH7=(RHOII(0)**(-I_HH))*RHOCUML
     &            *I_P*(I_P-1)/2.D0
            ELSE
               CALCPATH7=(RHOII(0)**(-I_HH))*RHOCUML
     &      *(-1  - I_P*S(1) + (1+S(1))**I_P)/ S(1)**2
            ENDIF
         ELSE
            CALCPATH7=0
C            STOP 'I_H > 4'
         ENDIF
         RETURN
      END

      REAL*8 FUNCTION WW4ABC(A,B,C,IPP)
         WW4ABC=
     &((-1 + (1 + A)**IPP)*B**2*(B - C)*C**2 + A*B**2*C**2*(-B + C)*IPP+
     -     A**3*(-((-1 + (1 + B)**IPP)*C**2) + B*C**2*IPP +
     -        B**2*(-1 + (1 + C)**IPP - C*IPP)) +
     -     A**2*((-1 + (1 + B)**IPP)*C**3 - B*C**3*IPP +
     -        B**3*(1 - (1 + C)**IPP + C*IPP)))/
     -   (A**2*(A - B)*B**2*(A - C)*(B - C)*C**2)

         RETURN
      END

      REAL*8 FUNCTION WW4AB1(A,B,IPP)
         WW4AB1=
     &   (2*A**3 - 2*B**3 + 2*(1 + A)**IPP*B**3 - 2*A**3*(1 + B)**IPP +
     -     2*A**3*B*IPP - A**3*B**2*IPP - 2*A*B**3*IPP + A**2*B**3*IPP +
     -     A**3*B**2*IPP**2 - A**2*B**3*IPP**2)/(2.*A**3*(A - B)*B**3)
 
         RETURN
      END

      REAL*8 FUNCTION WW4AA1(A,IPP)
         WW4AA1=
     &    (6 - 6*(1 + A)**IPP + 2*A*(2 + (1 + A)**(-1 + IPP))*IPP +
     -     A**2*(-1 + IPP)*IPP)/(2.*A**4)

         RETURN
      END

      REAL*8 FUNCTION WW4A11(A,IPP)
         WW4A11=
     &    -(6 - 6*(1 + A)**IPP + 6*A*IPP + 3*A**2*(-1 + IPP)*IPP +
     -      A**3*IPP*(2 - 3*IPP + IPP**2))/(6.*A**4)
         RETURN
      END

      REAL*8 FUNCTION WW4111(IPP)
         WW4111=(IPP*(-6 + 11*IPP - 6*IPP**2 + IPP**3))/24
         RETURN
      END

      REAL*8 FUNCTION WW4AAB(A,B,IPP)
         WW4AAB=
     &  (2*(-1 + (1 + A)**IPP)*B**3 + A**4*(-1 + (1 + B)**IPP - B*IPP) +
     -     A**3*(-1 + (1 + B)**IPP - B*IPP + 2*B**2*IPP) -
     - A**2*B**2*(3*(-1 + (1 + A)**IPP) + (-2 - (1 + A)**IPP + B)*IPP) -
     -     A*B**2*(3*(-1 + (1 + A)**IPP) +
     -        B*(2 - 2*(1 + A)**IPP + IPP + (1 + A)**IPP*IPP)))/
     -   (A**3*(1 + A)*(A - B)**2*B**2)

         RETURN
      END

      REAL*8 FUNCTION WW4AAA(A,IPP)
         WW4AAA=
     &    (6*(-1 + (1 + A)**IPP) - 2*A**3*IPP -
     -     2*A*(6 - 6*(1 + A)**IPP + IPP + 2*(1 + A)**IPP*IPP) +
     -     A**2*(-6 + 6*(1 + A)**IPP - 4*IPP - 5*(1 + A)**IPP*IPP +
     -        (1 + A)**IPP*IPP**2))/(2.*A**4*(1 + A)**2)
         RETURN
      END


      REAL*8 FUNCTION CALCPATH5(IPATH,RHOII,RHOCUML,I_H,
     &               I_P,NWHTAY)
         IMPLICIT REAL*8 (A-H,O-Z)
         IMPLICIT LOGICAL (L)
         INTEGER IPATH(*)
         REAL*8 RHOII(0:*)
         REAL*8 S(I_H)
         LOGICAL LFLIP
         IF(NWHTAY.GT.0) THEN
           CALCPATH5=CALCPATH6_2(IPATH,RHOII,RHOCUML,I_H,I_P,NWHTAY,1,0)
            RETURN
         ENDIF
         DO I=1,I_H-1
            S(I)=RHOII(I)/RHOII(0)
         ENDDO
C.. To determine which formula to use to take into account degenerecies,
C.. we need to sort S1,S2,S3.  We bubble sort (as I_H will never be too large)
         LFLIP=.TRUE.
         DO WHILE (LFLIP)
            LFLIP=.FALSE.
            DO I=1,I_H-2
               IF(S(I).GT.S(I+1)) THEN
                  ST=S(I)
                  S(I)=S(I+1)
                  S(I+1)=ST
                  LFLIP=.TRUE.
               ENDIF 
            ENDDO
         ENDDO
         IONE=0
         DO I=1,I_H-1
            IF(LSAME2(S(I),1.D0)) IONE=IONE+1
         ENDDO
C.. Now we go to specifics for each I_H
         IF(I_H.EQ.4) THEN
            IIP=I_P
            IF(LSAME2(S(1),S(2))) THEN
C..      XXY
               IF(LSAME2(S(2),S(3))) THEN
C..      XXX
                  IF(IONE.EQ.3) THEN
C..  111
                     CALCPATH5=W4111(IIP)
                  ELSE
C..  AAA
                     CALCPATH5=W4AAA(S(1),IIP)
                  ENDIF
               ELSE
C..      XXY
                  IF(IONE.EQ.2) THEN
C..  11A
                     CALCPATH5=W4A11(S(3),IIP)
                  ELSEIF(IONE.EQ.1) THEN
C..  AA1
                     CALCPATH5=W4AA1(S(1),IIP)
                  ELSE
C..  AAB
                     CALCPATH5=W4AAB(S(1),S(3),IIP)
                  ENDIF
               ENDIF
            ELSEIF(LSAME2(S(2),S(3))) THEN
C..     XYY
               IF(IONE.EQ.1) THEN
C..  1AA       
                  CALCPATH5=W4AA1(S(2),IIP)
               ELSEIF(IONE.EQ.2) THEN
C..  A11
                  CALCPATH5=W4A11(S(1),IIP)
               ELSE
C..  ABB
                  CALCPATH5=W4AAB(S(1),S(3),IIP)
               ENDIF
            ELSEIF(IONE.EQ.1) THEN
C..     XYZ
               IF(LSAME2(S(1),1.D0)) THEN
C..  1AB
                  CALCPATH5=W4AB1(S(2),S(3),IIP)
               ELSEIF(LSAME2(S(2),1.D0)) THEN
C..  A1B
                  CALCPATH5=W4AB1(S(1),S(3),IIP)
               ELSE
C..  AB1
                  CALCPATH5=W4AB1(S(1),S(2),IIP)
               ENDIF
            ELSE
C..  ABC
               CALCPATH5=W4ABC(S(1),S(2),S(3),IIP)
            ENDIF           

            CALCPATH5=(RHOII(0)**(-I_H))*RHOCUML*CALCPATH5
 
         ELSEIF(I_H.EQ.3) THEN
            IIP=I_P
            IF(LSAME2(S(1),S(2))) THEN
               IF(LSAME2(S(1),1.D0)) THEN
                  CALCPATH5=(IIP*(2 - 3*IIP + IIP**2))/6.
               ELSE
                  CALCPATH5=(-((-2 + IIP)*S(1)) + IIP*S(1)**2 - 
     &               IIP*S(1)**IIP
     &               +(-2 + IIP)*S(1)**(1 + IIP))/((-1 + S(1))**3*S(1)) 
               ENDIF  
            ELSEIF(LSAME2(S(1),1.D0).OR.LSAME2(S(2),1.D0)) THEN
               IF(LSAME2(S(2),1.D0)) S(2)=S(1)
               CALCPATH5=(-(IIP**2*(-1 + S(2))**2) + 
     &               IIP*(3 - 4*S(2) + S(2)**2) +
     -               2*(-1 + S(2)**IIP))/(2.*(-1 + S(2))**3)
            ELSE
               CALCPATH5=(S(1)**IIP*(-1 + S(2))**2 - 
     &               (-2 + IIP)*S(2) + (-1 + IIP)*S(2)**2 - S(2)**IIP -
     -               S(1)**2*(-1 + IIP - IIP*S(2) + S(2)**IIP) +
     -               S(1)*(-2 + IIP - IIP*S(2)**2 + 2*S(2)**IIP))/
     -               ((-1 + S(1))**2*(S(1) - S(2))*(-1 + S(2))**2)
            ENDIF
            CALCPATH5=(RHOII(0)**(-I_H))*RHOCUML*CALCPATH5

         ELSEIF(I_H.EQ.2) THEN
            IF(LSAME2(S(1),1.D0)) THEN
               CALCPATH5=(RHOII(0)**(-I_H))*RHOCUML
     &            *I_P*(I_P-1)/2.D0
            ELSE
               CALCPATH5=(RHOII(0)**(-I_H))*RHOCUML
     &      *(-1 + I_P - I_P*S(1) + S(1)**I_P)/(-1 + S(1))**2
            ENDIF
         ELSE
            STOP 'I_H > 3'
         ENDIF
         RETURN
      END


      REAL*8 FUNCTION W4ABC(A,B,C,IPP)
         W4ABC=
     &  (A**IPP*(-1 + B)**2*(B - C)*(-1 + C)**2 + B**IPP*(-1 + C)**2*C -
     -     B*(C**IPP + C**2*(-3 + IPP) - C**3*(-2 + IPP)) -
     -     B**3*(C**IPP + C*(-2 + IPP) - C**2*(-1 + IPP)) +
     -     B**2*(2*C**IPP + C*(-3 + IPP) - C**3*(-1 + IPP)) +
     -     A**3*(-(B**IPP*(-1 + C)**2) + C**IPP + C*(-2 + IPP) -
     -        C**2*(-1 + IPP) + B**2*(-1 + C**IPP + IPP - C*IPP) +
     -        B*(2 - 2*C**IPP - IPP + C**2*IPP)) +
     -     A**2*(-2*C**IPP + B**IPP*(-1 + C)**2*(2 + C) - C*(-3 + IPP) +
     -        C**3*(-1 + IPP) - B**3*(-1 + C**IPP + IPP - C*IPP) +
     -        B*(-3 + 3*C**IPP + IPP - C**3*IPP)) +
     -     A*(C**IPP - B**IPP*(-1 + C)**2*(1 + 2*C) + C**2*(-3 + IPP) -
     -        C**3*(-2 + IPP) + B**3*(-2 + 2*C**IPP + IPP - C**2*IPP) +
     -        B**2*(3 - 3*C**IPP - IPP + C**3*IPP)))/
     -   ((-1 + A)**2*(A - B)*(-1 + B)**2*(A - C)*(B - C)*(-1 + C)**2)
         RETURN
      END

      REAL*8 FUNCTION W4AB1(A,B,IPP)
         W4AB1=
     &   (2*A**IPP*(-1 + B)**3 + 2*B**IPP - B*(6 - 5*IPP + IPP**2) +
     -     2*B**2*(3 - 4*IPP + IPP**2) - B**3*(2 - 3*IPP + IPP**2) +
     -     A*(6 - 6*B**IPP - 5*IPP - 3*B**2*(-3 + IPP)*IPP +
     -        2*B**3*(-2 + IPP)*IPP + IPP**2) +
     -     A**3*(2 - 2*B**IPP - 3*IPP - 2*B*(-2 + IPP)*IPP +
     -        B**2*(-1 + IPP)*IPP + IPP**2) +
     -     A**2*(6*B**IPP + 3*B*(-3 + IPP)*IPP - B**3*(-1 + IPP)*IPP -
     -     2*(3 - 4*IPP + IPP**2)))/(2.*(-1 + A)**3*(A - B)*(-1 + B)**3)
 
         RETURN
      END

      REAL*8 FUNCTION W4AA1(A,IPP)
         W4AA1=
     & (2*A**(1 + IPP)*(-3 + IPP)-2*A**IPP*IPP - 2*A**2*(-3 + IPP)*IPP +
     -  A**3*(-1 + IPP)*IPP + A*(6 - 5*IPP + IPP**2))/(2.*(-1 + A)**4*A)
 
         RETURN
      END

      REAL*8 FUNCTION W4A11(A,IPP)
         W4A11=
     -      (-6 + 6*A**IPP + 11*IPP - 6*IPP**2 + IPP**3 -
     -  3*A*IPP*(6 - 5*IPP + IPP**2) + 3*A**2*IPP*(3 - 4*IPP + IPP**2) -
     -     A**3*IPP*(2 - 3*IPP + IPP**2))/(6.*(-1 + A)**4)

         RETURN
      END

      REAL*8 FUNCTION W4111(IPP)
         W4111=(IPP*(-6 + 11*IPP - 6*IPP**2 + IPP**3))/24
         RETURN
      END

      REAL*8 FUNCTION W4AAB(A,B,IPP)
         W4AAB=(A**IPP*(-1 + B)**3*B +
     -     B*(-(B**2*(-3 + IPP)) + B**(1 + IPP)*(-3 + IPP) +
     -        B**3*(-1 + IPP) - B**IPP*(-1 + IPP)) +
     -     A**2*(-(B*(-3 + IPP)) - B**2*(-3 + IPP) +
     -        B**(1 + IPP)*(-3 + IPP) + B**(2 + IPP)*(-3 + IPP) +
     -        2*B**3*IPP - 2*B**IPP*IPP) +
     -     A**3*(B*(-2 + IPP) - B**(1 + IPP)*(-2 + IPP) - B**2*IPP +
     -        B**IPP*IPP) + A*(2*B**2*(-3 + IPP) -
     -    2*B**(2 + IPP)*(-3 + IPP) - B**3*IPP - B**4*IPP + B**IPP*IPP +
     -        B**(1 + IPP)*IPP))/((-1 + A)**2*(A - B)**2*(-1 + B)**3*B)
         RETURN
      END

      REAL*8 FUNCTION W4AAA(A,IPP)
         W4AAA=
     - (2*A**2*(-3 + IPP) - 2*A**3*IPP - 2*A**(1 + IPP)*(-3 + IPP)*IPP +
     -     A**IPP*(-1 + IPP)*IPP + A**(2 + IPP)*(6 - 5*IPP + IPP**2))/
     -   (2.*(-1 + A)**4*A**2)
         RETURN
      END


C.. Recursively generate Z-sum for use with MCPATHS
C.. use formula with approximations - not v good (14/11/03)
C.. R_R(R,X,NL) = (1/(1-x))(R_R(r,s_(nl-1)/r,nl-1)-x*R_R(r*x,s_(nl-1)/rx,nl-1))
C.. R_R(R,X,0) = r**L0 (1-x**(L0+1))/(1-x)
      REAL*8 RECURSIVE FUNCTION R_R2(RHO,CHI,NL,RHOII,L0) 
     &    RESULT (R_R2RES)
         IMPLICIT REAL*8 (A-H,O-Z)
         REAL*8 RHOII(0:*)
         R_R=0.D0
         CCHI=CHI
C         WRITE(10,"(A,2E25.16,I2)") "RR(",RHO,CHI,NL
         IF(NL.EQ.0) THEN
            IF(CCHI.EQ.0.D0)  THEN
               CCHI=1D-10
            ENDIF
C               R_R=(R**L0)*(L0+1)
C            ELSE
               R_R2RES=((1+RHO)**L0)*(1-(1+CCHI)**(L0+1))/(-CCHI)
C            ENDIF
         ELSE
C.. if we have a degeneracy, the sums blow up, so we break the degeneracy
C.. very slightly
            SIGMAN=(RHOII(NL-1)/RHOII(0))-1
            IF(CCHI.EQ.0.D0) THEN
               CCHI=1D-10
            ENDIF
C               RD=R_R(R+EPS,SN/(R+EPS),NL-1,RHOII,L0)
C               R2=R_R(R,SN/R,NL-1,RHOII,L0)
C               R_R=(RD-R2)/EPS+RD
C            ELSE
C            IF(CX.EQ.1.D0) RETURN
               R1=R_R2(RHO,SIGMAN-RHO,NL-1,RHOII,L0)
C            IF(R1.EQ.0.D0) RETURN
               R2=(1+CCHI)*R_R2(RHO+CHI+RHO*CHI,SIGMAN-(RHO+CHI+RHO*CHI)
     &               ,NL-1,RHOII,L0)
C            IF(R2.EQ.0.D0) RETURN
               R1=R1-R2
               R_R2RES=-R1/CCHI
C            ENDIF
         ENDIF
         WRITE(10,"(A,2E25.16,I2,A)",advance='no') "RR(",RHO,CHI,NL,")"
         WRITE(10,*) R_R2RES
         RETURN
      END


C.. use formula
C.. R_R(R,X,NL) = (1/(1-x))(R_R(r,s_(nl-1)/r,nl-1)-x*R_R(r*x,s_(nl-1)/rx,nl-1))
C.. R_R(R,X,0) = r**L0 (1-x**(L0+1))/(1-x)
      REAL*8 RECURSIVE FUNCTION R_R(R,X,NL,RHOII,L0) RESULT (R_RRES)
         IMPLICIT REAL*8 (A-H,O-Z)
         REAL*8 RHOII(0:*)
         R_RRES=0.D0
         CX=X
         WRITE(10,"(A,2F19.7,I2)") "RR(",R,X,NL
         IF(NL.EQ.0) THEN
            IF(CX.EQ.1.D0)  THEN
               CX=1.D0+1D-10
            ENDIF
C               R_R=(R**L0)*(L0+1)
C            ELSE
               R_RRES=(R**L0)*(1-CX**(L0+1))/(1-CX)
C            ENDIF
         ELSE
C.. if we have a degeneracy, the sums blow up, so we break the degeneracy
C.. very slightly
            SN=RHOII(NL-1)/RHOII(0)
            IF(CX.EQ.1.D0) THEN
               CX=1.D0+1D-10
            ENDIF
C               RD=R_R(R+EPS,SN/(R+EPS),NL-1,RHOII,L0)
C               R2=R_R(R,SN/R,NL-1,RHOII,L0)
C               R_R=(RD-R2)/EPS+RD
C            ELSE
C            IF(CX.EQ.1.D0) RETURN
               R1=R_R(R,SN/R,NL-1,RHOII,L0)
C            IF(R1.EQ.0.D0) RETURN
               R2=CX*R_R(R*CX,SN/(R*CX),NL-1,RHOII,L0)
C            IF(R2.EQ.0.D0) RETURN
               R1=R1-R2
               R_RRES=R1/(1-CX)
C            ENDIF
         ENDIF
         WRITE(10,"(A,2F19.7,I2,A)",advance='no') "RR(",R,X,NL,")"
         WRITE(10,*) R_RRES
         RETURN
      END

C.. Recursively go through each possible path with nodes IPATH, and
C.. sum.  I_HIND is the recursion variable, specifying the node of interest
C.. When called it places zero to I_BOUND self-hops on node I_HIND, then
C.. recurses to the next node. (5/11/03)
      RECURSIVE REAL*8 FUNCTION CALCPATH(IPATH,RHOII,RHOCUML,I_H,I_HIND,
     &                        I_BOUND) RESULT (CALCPATHRES)
         IMPLICIT REAL*8 (A-H,O-Z)
         INTEGER IPATH(*)
         REAL*8 RHOII(0:*)
         
         TOTAL=0.D0
         RII=1.D0
         IF(I_HIND.NE.I_H) THEN
            DO I=1,I_BOUND
C.. Display
C               DO J=1,I_HIND-1
C                  WRITE(6,"A$"),"   "
C               ENDDO
C               WRITE(10,"I2,A$"),I,","
               TOTAL=TOTAL+
     &            CALCPATH(IPATH,RHOII,RHOCUML*RII,I_H,I_HIND+1,
     &                  I_BOUND-(I-1))
               RII=RII*RHOII(I_HIND)
            ENDDO
            CALCPATHRES=TOTAL
         ELSE
            RC=RHOCUML*(RHOII(I_HIND)**(I_BOUND-1))
C            WRITE(10,"I2,F"),I_BOUND,RC
C.. RHOCUML now contains the appropriate rho for the path
C.. Output the path, and RHOCUML, and continue
C            WRITE(10,*) ((IPATH(I),"*",,I=0,I_H)
            CALCPATHRES=RC
         ENDIF
         RETURN
      END

C.. GETRHOELEMENT
C.. Get matrix element of rho
      REAL*8 FUNCTION GETRHOELEMENT(II,IJ,RHOIJ,RHOIJLAB,RHOIJROWS,NDET)
      REAL*8 RHOIJ(*)
      INTEGER RHOIJLAB(*)
      INTEGER RHOIJROWS(NDET)
C.. We only have half of RHO, so if J<I, return the symmetrical (J,I) element
C.. Or if we have the whole RHO, it's quicker to look closer to its beginning
      IF(IJ<II) THEN
         I=IJ
         J=II
      ELSE
         I=II
         J=IJ
      ENDIF
C      IF(I.GT.0) THEN
      GETRHOELEMENT=0.D0
      INDXROW=RHOIJROWS(I)
      IMAX=RHOIJROWS(I+1)-1
      DO K=INDXROW,IMAX
         IF(RHOIJLAB(K).GT.J) RETURN
         IF(RHOIJLAB(K).EQ.J) THEN
            GETRHOELEMENT=RHOIJ(K)
            RETURN
         ENDIF
      ENDDO
C      ELSE
C         JK=1
C      ENDIF
      GETRHOELEMENT=0.D0
      RETURN
      END

C.. Recursively sum all paths with nodes IPATH, using geometric series
C.. sum techniques. (7/11/03)
C.. RHOCUML contains RHO_i0i1*RHO_i1i2...
C.. This calls a recursive function which recurses over the nodes.
C.. This works for I_H<=3 (exactly)
      REAL*8 FUNCTION CALCPATH2(IPATH,RHOII,RHOCUML,I_H,
     &               I_P)
         IMPLICIT REAL*8 (A-H,O-Z)
         INTEGER IPATH(*)
         REAL*8 RHOII(0:*)
         CALCPATH2=(RHOII(0)**(I_P-I_H))*RHOCUML
     &         *R_R(1.D0,RHOII(I_H-1)/RHOII(0),I_H-1,RHOII,I_P-I_H)
C         CALCPATH2=(RHOII(0)**(I_P-I_H))*RHOCUML
C     &         *R_R2(0.D0,(RHOII(I_H-1)/RHOII(0))-1,I_H-1,RHOII,I_P-I_H)
         RETURN
      END

C.. Sum all paths up to order I_H=4, using 1st order
C.. Taylor expansions (18/11/03)
      REAL*8 FUNCTION CALCPATH3(IPATH,RHOII,RHOCUML,I_H,
     &               I_P)
         IMPLICIT REAL*8 (A-H,O-Z)
         INTEGER IPATH(*)
         REAL*8 RHOII(0:*)
         IF(I_H.EQ.4) THEN
            S1=RHOII(1)/RHOII(0)
            S2=RHOII(2)/RHOII(0)
            S3=RHOII(3)/RHOII(0)
            CALCPATH3=(RHOII(0)**(I_P-I_H))*RHOCUML
     &  *(I_P*(-6 + 11*I_P - 6*I_P**2 + I_P**3)*
     &  (1254 - 568*S2 + I_P**3*(-1 + S1)*(-1 + S2)*(-1 + S3) - 568*S3 +
     &       260*S2*S3 - 4*S1*(142 - 65*S2 - 65*S3 + 30*S2*S3) +
     &       I_P**2*(36 - 29*S2 - 29*S3 + 22*S2*S3 +
     &          S1*(-29 + 22*S2 + 22*S3 - 15*S2*S3)) +
     &       I_P*(-389 + 242*S2 + 242*S3 - 137*S2*S3 +
     &          S1*(242 - 137*S2 - 137*S3 + 74*S2*S3))))/5040.
  
         ELSEIF(I_H.EQ.3) THEN
            S1=RHOII(1)/RHOII(0)
            S2=RHOII(2)/RHOII(0)
            CALCPATH3=(RHOII(0)**(I_P-I_H))*RHOCUML
     &*(I_P*(2 -3*I_P + I_P**2)*(62 + I_P**2*(-1 +S1)*(-1 +S2) - 27*S2 +
     &       3*S1*(-9 + 4*S2) + I_P*(-17 + S1*(12 - 7*S2) + 12*S2)))/120
         ELSEIF(I_H.EQ.2) THEN
            S1=RHOII(1)/RHOII(0)
            CALCPATH3=(RHOII(0)**(I_P-I_H))*RHOCUML
     &      *((-1 + I_P)*I_P*(3 + (-2 + I_P)*(-1 + S1)))/6
 
         ELSE
            STOP 'I_H > 4'
         ENDIF
         RETURN
      END
C.. Sum all paths up to order I_H=4, using 2nd order
C.. Taylor expansions (18/11/03)
      REAL*8 FUNCTION CALCPATH4(IPATH,RHOII,RHOCUML,I_H,
     &               I_P)
         IMPLICIT REAL*8 (A-H,O-Z)
         INTEGER IPATH(*)
         REAL*8 RHOII(0:*)
         IF(I_H.EQ.4) THEN
            S1=RHOII(1)/RHOII(0)
            S2=RHOII(2)/RHOII(0)
            S3=RHOII(3)/RHOII(0)
            CALCPATH4=(RHOII(0)**(I_P-I_H))*RHOCUML
     &*     (I_P*(-6 + 11*I_P - 6*I_P**2 + I_P**3)*
     -     (I_P**6*(-1 + S1)**2*(-1 + S2)**2*(-1 + S3)**2 -
     -       I_P**5*(-1 + S1)*(-1 + S2)*(-1 + S3)*
     -        (-69 + S2*(59 - 49*S3) + 59*S3 +
     -          S1*(59 - 49*S3 + S2*(-49 + 39*S3))) +
     -       5*I_P**4*(413 - 694*S3 + 299*S3**2 +
     -          S2**2*(299 - 484*S3 + 203*S3**2) -
     -          2*S2*(347 - 571*S3 + 242*S3**2) +
     -          S1**2*(299 - 484*S3 + 203*S3**2 +
     -             S2*(-484 + 758*S3 - 310*S3**2) +
     -             S2**2*(203 - 310*S3 + 125*S3**2)) -
     -          2*S1*(347 - 571*S3 + 242*S3**2 +
     -             S2*(-571 + 914*S3 - 379*S3**2) +
     -             S2**2*(242 - 379*S3 + 155*S3**2))) -
     -       5*I_P**3*(6567 - 9832*S3 + 3949*S3**2 -
     -          8*S2*(1229 - 1772*S3 + 696*S3**2) +
     -          S2**2*(3949 - 5568*S3 + 2159*S3**2) -
     -          8*S1*(1229 - 1772*S3 + 696*S3**2 +
     -             S2*(-1772 + 2433*S3 - 931*S3**2) +
     -             S2**2*(696 - 931*S3 + 352*S3**2)) +
     -          S1**2*(3949 - 5568*S3 + 2159*S3**2 -
     -             8*S2*(696 - 931*S3 + 352*S3**2) +
     -             S2**2*(2159 - 2816*S3 + 1053*S3**2))) +
     -       240*(10164 - 9128*S3 + 2897*S3**2 +
     -          S2*(-9128 + 8303*S3 - 2649*S3**2) +
     -          S2**2*(2897 - 2649*S3 + 847*S3**2) +
     -          S1**2*(2897 - 2649*S3 + 847*S3**2 +
     -             S2*(-2649 + 2443*S3 - 784*S3**2) +
     -             7*S2**2*(121 - 112*S3 + 36*S3**2)) +
     -          S1*(-9128 + 8303*S3 - 2649*S3**2 +
     -             S2**2*(-2649 + 2443*S3 - 784*S3**2) +
     -             S2*(8303 - 7626*S3 + 2443*S3**2))) +
     -       2*I_P**2*(144587 - 188644*S3 + 70472*S3**2 +
     -          S2*(-188644 + 233963*S3 - 85189*S3**2) +
     -          S2**2*(70472 - 85189*S3 + 30692*S3**2) +
     -          S1**2*(70472 - 85189*S3 + 30692*S3**2 +
     -             S2*(-85189 + 98603*S3 - 34924*S3**2) +
     -             S2**2*(30692 - 34924*S3 + 12287*S3**2)) +
     -          S1*(-188644 + 233963*S3 - 85189*S3**2 +
     -             S2**2*(-85189 + 98603*S3 - 34924*S3**2) +
     -             S2*(233963 - 276226*S3 + 98603*S3**2))) -
     -       4*I_P*(329724 - 365408*S3 + 125819*S3**2 +
     -          S2*(-365408 + 386231*S3 - 131073*S3**2) +
     -          S2**2*(125819 - 131073*S3 + 44269*S3**2) +
     -          S1**2*(125819 - 131073*S3 + 44269*S3**2 +
     -             S2*(-131073 + 134011*S3 - 44968*S3**2) +
     -             S2**2*(44269 - 44968*S3 + 15054*S3**2)) +
     -          S1*(-365408 + 386231*S3 - 131073*S3**2 +
     -             S2**2*(-131073 + 134011*S3 - 44968*S3**2) +
     -             S2*(386231 - 398202*S3 + 134011*S3**2)))))/3.6288e6
 
 
         ELSEIF(I_H.EQ.3) THEN
            S1=RHOII(1)/RHOII(0)
            S2=RHOII(2)/RHOII(0)
            CALCPATH4=(RHOII(0)**(I_P-I_H))*RHOCUML
     &*         (I_P*(2 - 3*I_P + I_P**2)*
     -     (I_P**4*(-1 + S1)**2*(-1 + S2)**2 -
     -     I_P**3*(-1 + S1)*(-1 + S2)*(32 - 25*S2 + S1*(-25 + 18*S2)) +
     -       7*I_P**2*(59 - 88*S2 + 35*S2**2 +
     -    S1**2*(35 - 46*S2 + 17*S2**2) - 2*S1*(44 - 61*S2 + 23*S2**2))
     -        + 6*(802 - 687*S2 + 214*S2**2 +
     -          S1*(-687 + 604*S2 - 190*S2**2) +
     -          2*S1**2*(107 - 95*S2 + 30*S2**2)) +
     -       I_P*(-2302 + 2763*S2 - 965*S2**2 +
     -          S1**2*(-965 + 1013*S2 - 342*S2**2) +
     -          S1*(2763 - 2978*S2 + 1013*S2**2))))/5040.

         ELSEIF(I_H.EQ.2) THEN
            S1=RHOII(1)/RHOII(0)
            CALCPATH4=(RHOII(0)**(I_P-I_H))*RHOCUML
     &*         ((-1 + I_P)*I_P*(12 + 4*(-2 + I_P)*(-1 + S1) +
     -       (-3 + I_P)*(-2 + I_P)*(-1 + S1)**2))/24.
 
         ELSE
            STOP 'I_H > 4'
         ENDIF
         RETURN
      END

