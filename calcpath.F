
      FUNCTION CALCPATHS(IPATH,RHOII,RHOIJ,I_V,
     &      I_HMAX,I_P,NWHTAY,NEL,I_VMAX,ILOGGING)
         USE HElem
         use global_utilities
         IMPLICIT NONE
         TYPE(HDElement) CALCPATHS
         INTEGER NEL,I_V,I_P,NWHTAY,I_HMAX,I_VMAX
C.. IPATH contains the vertices we can choose from
C.. IPATH2 contains the path we are constructing (listing the dets)
C.. IP2 contains the path we are constructing as an index into IPATH
         INTEGER IPATH(1:NEL,0:I_V),IPATH2(1:NEL,0:I_HMAX),IP2(0:I_HMAX)
C.. RHOII contains the RHO_II for each vertex.
C.. 12/7/05 - this should now just be the number 0.D0 - it is not used


C.. RHOII2 is constructed to contain the RHOII for each point in the path
         TYPE(HElement) RHOIJ(0:I_V,0:I_V),RH
         TYPE(HDElement) TOTAL,RHOCUML,RHOII(0:I_V),RHOII2(0:I_HMAX)
         TYPE(HDElement) DLWDB
         INTEGER ILOGGING
         type(timer), save :: proc_timer
         TYPE(HDElement) RhoDiag_CPP, CalcPaths_2, CalcPaths_3,
     &      CalcPaths_4, CalcPaths_R
         INTEGER I_HIND

         proc_timer%timer_name='CALCPATHS '
         call set_timer(proc_timer,50)
C.. If we're using the contour integral trick to sum all possible graphs with the given vertices
C.. REMEMBER TO MODIFY CALCPATHS_N too!
         IF(I_HMAX.EQ.-1.OR.(I_HMAX.LE.-3.AND.I_HMAX.GE.-8)) THEN
C.. we use a diagonalization version rather than contour etc.
C.. subtracting out lower vertex paths at this level
C            WRITE(6,*) "-----"
            CALCPATHS=RHODIAG_CPP(RHOIJ,I_P,I_V,I_V+1,.TRUE.,
     &         0.D0, DLWDB)
         ELSEIF(I_HMAX.EQ.-2.OR.I_HMAX.EQ.-6) THEN
C.. we use a diagonalization version rather than contour etc.
C.. we deal with removing lower vertex paths later
C            WRITE(6,*) "-----"
            CALCPATHS=RHODIAG_CPP(RHOIJ,I_P,I_V,I_V+1,.FALSE.,
     &         0.D0, DLWDB)
         ELSEIF(I_HMAX.EQ.0) THEN
            WRITE(6,*) "I_HMAX=0"
            STOP 'WARNING: Running old disfunctional CALCPATHS I_HMAX=0'
            CALCPATHS=0.D0
            IF(I_V.EQ.2) CALCPATHS=CALCPATHS_2(RHOII,RHOIJ,I_P)
            IF(I_V.EQ.3) CALCPATHS=CALCPATHS_3(RHOII,RHOIJ,I_P)
            IF(I_V.EQ.4) CALCPATHS=CALCPATHS_4(RHOII,RHOIJ,I_P)
C            IF(I_V.EQ.3) CALCPATHS=CP3CI2(RHOII,RHOIJ,I_P)
         ELSE
            WRITE(6,*) "I_HMAX>0"
            STOP 'WARNING: Running old disfunctional CALCPATHS I_HMAX>0'
      
            CALCPATHS=CALCPATHS_R(IPATH,RHOII,RHOIJ,I_V,I_HMAX,I_P,
     &               NWHTAY,I_HIND,NEL,IPATH2,IP2,1.D0,RHOII2,I_VMAX,
     &               ILOGGING)
         ENDIF
         call halt_timer(proc_timer)
         RETURN
      END 

      FUNCTION CALCPATHS_N(IPATH,RHOII,RHOIJ,I_V,
     &      I_HMAX,I_P,WEIGHT,NEL,I_VMAX,ILOGGING,DBETA,DLWDB,HIJS,ICLS)
         USE HElem
         use CalcData , only : CALCP_SUB2VSTAR,CALCP_LOGWEIGHT
         use global_utilities
         IMPLICIT NONE
         TYPE(HDElement) CalcPaths_N
         INTEGER NEL,I_V,I_P,NWHTAY,I_HMAX,I_VMAX,ICLS
C.. IPATH contains the vertices we can choose from
C.. IPATH2 contains the path we are constructing (listing the dets)
C.. IP2 contains the path we are constructing as an index into IPATH
         INTEGER IPATH(1:NEL,0:I_V)
!Removed as no longer used
!  ,IPATH2(1:NEL,0:i_hmax),IP2(0:i_hmax)
C.. RHOII contains the RHO_II for each vertex.
C.. 12/7/05 - this should now just be the number 0.D0 - it is not used


C.. RHOII2 is constructed to contain the RHOII for each point in the path
         TYPE(HDElement) RHOII(0:I_V),TOTAL,RHOCUML
! removed as no longer used
!         TYPE(HDElement) RHOII2(0:i_hmax)
         TYPE(HDElement) CALCPATHS
         TYPE(HElement) HIJS(1:I_V),RHOIJ(0:I_V,0:I_V),RH
         INTEGER ICLASS,ILG,ISGN
         REAL*8 WEIGHT
!DBETA has two parts
!DBETA(1) is DBETA
!DBETA(2) is GRAPHEPSILON
         REAL*8 DBETA(2)
         TYPE(HDElement) DLWDB
         INTEGER ILOGGING
         type(timer), save :: proc_timer
         TYPE(HDElement) RhoDiag_CPP, CalcPaths_2, CalcPaths_3,
     &      CalcPaths_4, CalcPaths_R,HDiag_CPP
         TYPE(HDElement) RhoDiag_CPPS2VS,HDiag_CPPS2VS
         INTEGER I_HIND
         proc_timer%timer_name='CALCPATHS '
         call set_timer(proc_timer,50)
C.. If we're using the contour integral trick to sum all possible graphs with the given vertices
         IF(I_HMAX.EQ.-1.OR.(I_HMAX.LE.-3.AND.I_HMAX.GE.-8).OR.
     &      I_HMAX.LE.-12) THEN
C.. we use a diagonalization version rather than contour etc.
C.. subtracting out lower vertex paths at this level
C            WRITE(6,*) "-----"
            IF(I_HMAX.EQ.-19.OR.I_HMAX.EQ.-20) THEN
!.NOT.(RHOII(0).AGT.0.D0)) THEN
C.. We actually diagonalize a matrix of Hij elements rather than a matrix of
C.. RHOIJ elements.  We send beta in from RHOII(1)
               if(calcp_sub2vstar) then
                  CALCPATHS=HDIAG_CPP(RHOIJ,I_P,I_V,I_V+1,.FALSE.,
     &                RHOII(1),DLWDB,HIJS)
                 CALCPATHS=CALCPATHS+HDIAG_CPPS2VS(RHOIJ,I_P,I_V,I_V+1,
     &               .TRUE.,RHOII(1),DLWDB,HIJS)
               else
                  CALCPATHS=HDIAG_CPP(RHOIJ,I_P,I_V,I_V+1,.TRUE.,
     &            RHOII(1),DLWDB,HIJS)
               endif
            ELSE
               if(calcp_sub2vstar) then
                  CALCPATHS=RHODIAG_CPP(RHOIJ,I_P,I_V,I_V+1,
     &            .FALSE.,DBETA,DLWDB,HIJS)
               CALCPATHS=CALCPATHS+RHODIAG_CPPS2VS(RHOIJ,I_P,I_V,I_V+1,
     &            .TRUE.,DBETA,DLWDB,HIJS)
               else
                  CALCPATHS=RHODIAG_CPP(RHOIJ,I_P,I_V,I_V+1,
     &            .TRUE.,DBETA,DLWDB,HIJS,calcp_logweight)
               endif
            ENDIF
         ELSEIF(I_HMAX.EQ.-2.OR.I_HMAX.EQ.-6) THEN
C.. we use a diagonalization version rather than contour etc.
C.. we deal with removing lower vertex paths later
C            WRITE(6,*) "-----"
            CALCPATHS=RHODIAG_CPP(RHOIJ,I_P,I_V,I_V+1,.FALSE.,
     &         DBETA,DLWDB,HIJS,calcp_logweight)
         ELSEIF(I_HMAX.EQ.0) THEN
            WRITE(6,*) "I_HMAX=0"
            STOP 'WARNING: Running old disfunctional CALCPATHS I_HMAX=0'
            CALCPATHS=0.D0
            IF(I_V.EQ.2) CALCPATHS=CALCPATHS_2(RHOII,RHOIJ,I_P)
            IF(I_V.EQ.3) CALCPATHS=CALCPATHS_3(RHOII,RHOIJ,I_P)
            IF(I_V.EQ.4) CALCPATHS=CALCPATHS_4(RHOII,RHOIJ,I_P)
C            IF(I_V.EQ.3) CALCPATHS=CP3CI2(RHOII,RHOIJ,I_P)
         ELSE
            WRITE(6,*) "I_HMAX>0"
            STOP 'WARNING: Running old disfunctional CALCPATHS I_HMAX>0'
!            CALCPATHS=CALCPATHS_R(IPATH,RHOII,RHOIJ,I_V,I_HMAX,I_P,
!     &               NWHTAY,I_HIND,NEL,IPATH2,IP2,1.D0,RHOII2,I_VMAX,
!     &               ILOGGING)
         ENDIF
!DBETA(2) is GRAPHEPSILON
         IF(.NOT.(CALCPATHS.AGE.DBETA(2))) THEN
            CALCPATHS=0.D0
            DLWDB=0.D0
         ENDIF
         CALCPATHS_N=CALCPATHS
         IF(I_V.LE.5) THEN
            IF(WEIGHT.GT.0.D0)
     &       CALL CLASSPATHS(DREAL(CALCPATHS)/WEIGHT,
     &            DREAL(DLWDB)/WEIGHT,WEIGHT,RHOIJ,
     &            I_V,ICLS)
         ENDIF
         call halt_timer(proc_timer)
         RETURN
      END
      SUBROUTINE CLASSPATHS(RHOC,HC,WEIGHT,RHOIJ,I_V,ICLASS)
         USE HElem
         IMPLICIT NONE
         INCLUDE 'calcp.inc'
         TYPE(HElement) RHOIJ(0:I_V,0:I_V)
         REAL*8 RHOC,HC,WEIGHT
         INTEGER I_V,ILG,ISGN,ICLASS
         INTEGER INCR
C.. allow for subtractions
         IF(I_V.GT.5) RETURN
         INCR=1
         IF(WEIGHT.LT.0.D0) INCR=-1
         IF(ICLASS.EQ.0)
     &      CALL CLASSPATH(RHOIJ,I_V,ICLASS)
C            WRITE(56,*) ICLASS,RHOC
C            WRITE(32,*) I_V,ICLASS,RHOC
            ILG=LOG(ABS(RHOC))/LOG(10.D0)
            IF(ILG.GE.0) THEN
               ILG=1
            ELSEIF(ILG.LT.-15) THEN
               ILG=15
            ELSE
               ILG=-ILG
            ENDIF
            IF(RHOC.GT.0.D0) THEN
               ISGN=1
            ELSEIF(RHOC.LT.0.D0) THEN
               ISGN=2
            ELSE
               ISGN=0
            ENDIF
            IF(ISGN.NE.0)
     &         CALCP_NT(ICLASS,ILG,ISGN)=CALCP_NT(ICLASS,ILG,ISGN)+INCR
            CALCP_N(ICLASS)=CALCP_N(ICLASS)+INCR
            CALCP_SUM(ICLASS,1)=CALCP_SUM(ICLASS,1)+RHOC
            IF(ISGN.NE.0) THEN
               CALCP_SUM(ICLASS,ISGN+1)=CALCP_SUM(ICLASS,ISGN+1)+RHOC
            ENDIF
            CALCP_HSUM(ICLASS)=CALCP_HSUM(ICLASS)+HC
      END
      SUBROUTINE CLASSPATH(RHOIJ,N,ICLASS)
         INTEGER I,J,N,ICLASS
         REAL*8 RHOIJ(0:N,0:N)
         ICLASS=0
         DO I=0,N-1
            DO J=I+1,N-1
               ICLASS=ICLASS*2
               IF(RHOIJ(I,J).NE.0.D0) THEN
                  ICLASS=ICLASS+1
               ENDIF
            ENDDO
         ENDDO
      END

      SUBROUTINE WRITECLASSPATHS()
         INCLUDE 'calcp.inc'
         INTEGER I,J
         OPEN(45,FILE='CLASSPATHS',STATUS='UNKNOWN')
         OPEN(46,FILE='CLASSPATHS2',STATUS='UNKNOWN')
         DO I=0,1023
            IF(CALCP_N(I).GT.0) THEN
               WRITE(45,"(I5,I11,4G19.11)") I,CALCP_N(I),CALCP_SUM(I,1),
     &            CALCP_HSUM(I),CALCP_SUM(I,2),CALCP_SUM(I,3)
               WRITE(46,"(I5,I11,2G19.11)") I,CALCP_N(I),CALCP_SUM(I,2),
     &            CALCP_SUM(I,3)
               DO J=1,15
                  WRITE(46,*) -J,CALCP_NT(I,J,1),CALCP_NT(I,J,2)
               ENDDO
            ENDIF
         ENDDO
         CLOSE(45)
         CLOSE(46)
      END


      SUBROUTINE CLEARCLASSPATHS()
         INCLUDE 'calcp.inc'
         CALCP_N=0
         CALCP_N=0
         CALCP_SUM=0.d0
         CALCP_HSUM=0.d0
      END

      RECURSIVE REAL*8 FUNCTION CALCPATHS_R(IPATH,RHOII,RHOIJ,I_V,
     &      I_HMAX,I_P,NWHTAY,I_HIND,NEL,IPATH2,IP2,RHOCUML,RHOII2,
     &      I_VMAX,ILOGGING) RESULT (CALCPATHS_RRES)
         IMPLICIT REAL*8 (A-H,O-Z)
         INTEGER NEL,I_V,I_P,NWHTAY,I_HMAX,I_VMAX
C.. IPATH contains the vertices we can choose from
C.. IPATH2 contains the path we are constructing (listing the dets)
C.. IP2 contains the path we are constructing as an index into IPATH
         INTEGER IPATH(1:NEL,0:I_V),IPATH2(1:NEL,0:I_HMAX),IP2(0:I_HMAX)
         INTEGER IP1(0:I_V)
C.. RHOII contains the RHO_II for each vertex.
C.. RHOII2 is constructed to contain the RHOII for each point in the path
         REAL*8 RHOII(0:I_V),RHOIJ(0:I_V,0:I_V),TOTAL,RHOCUML,RH
         REAL*8 RHOII2(0:I_HMAX),T2
         INTEGER I_H,I_HIND,I,J,TOT
         TOTAL=0.D0
c           WRITE(10,"A$") "{{"
c           DO J=0,I_HIND-1
c              WRITE(10,"A$") "("
c              DO K=1,NEL
c                 WRITE(10,"I3,A$") IPATH2(K,J),","
c              ENDDO
c              WRITE(10,"A$") "),"
c           ENDDO
c           WRITE(10,*) "}}"
C.. We need to sum over all permutations of IJK...L which are valid, up to a
C.. length I_HMAX.  We do this crudely here.  Slight optimizations are possible
         IF(I_HIND.EQ.0) THEN
C.. first node must be I
            CALL NECI_ICOPY(NEL,IPATH(1:NEL,0),1,
     &                     IPATH2(1:NEL,0),1)
            IP2(0)=0
            RHOII2(0)=RHOII(0)
            TOTAL=TOTAL+CALCPATHS_R(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &               I_P,NWHTAY,1,NEL,IPATH2,IP2,RHOCUML,RHOII2,I_VMAX,
     &               ILOGGING)
         ELSE
            DO I=0,I_V-1
C.. we can't have two successive nodes the same
               IF(IP2(I_HIND-1).NE.I) THEN
                  RH=RHOIJ(IP2(I_HIND-1),I)
                  IF(RH.NE.0.D0) THEN
C.. we try all possible vertices in this position.  If one happens to be I
C.. then we calculate the Z-sum for that, and add it into the total.  We then
C.. recurse to the next node if we're below I_HMAX)
                     CALL NECI_ICOPY(NEL,IPATH(1:NEL,I),1,
     &                     IPATH2(1:NEL,I_HIND),1)
                     IP2(I_HIND)=I
                     RHOII2(I_HIND)=RHOII(I)
                     IF(I.EQ.0) THEN
C.. need to check we have used all the vertices at least once
                        DO J=0,I_V
                           IP1(J)=0
                        ENDDO
                        TOT=0
                        DO J=0,I_HIND
                           IF(IP1(IP2(J)).EQ.0) THEN
                              IP1(IP2(J))=1
                              TOT=TOT+1
                           ENDIF
                        ENDDO
                        IF(TOT.EQ.I_V) THEN
C.. do the Z-sum
                           T2=CALCPATH7(IPATH2,RHOII2,RHOCUML*RH,I_HIND,
     &                        I_P,NWHTAY)
                           TOTAL=TOTAL+T2
         IF(BTEST(ILOGGING,4)) THEN
            WRITE(10,"(A)",advance='no') "{"
            DO J=0,I_HIND
               WRITE(10,"(A)",advance='no') "("
               DO K=1,NEL
                  WRITE(10,"(I3,A)",advance='no') IPATH2(K,J),","
               ENDDO
               WRITE(10,"(A)",advance='no') "),"
            ENDDO
            WRITE(10,"(A)",advance='no') "}"
            IF(BTEST(ILOGGING,5)) THEN
               WRITE(10,"(A)",advance='no') "("
               WRITE(10,"(E25.16,A)",advance='no') 1.D0,","
               DO J=1,I_HIND-1
                  WRITE(10,"(E25.16,A)",advance='no')
     &                                       RHOIJ(IP2(J),IP2(J-1)),","
               ENDDO
               WRITE(10,"(E25.16)",advance='no')
     &                                      RHOIJ(IP2(0),IP2(I_HIND-1))
               WRITE(10,*) ")"
            ENDIF
            WRITE(10,*) T2
         ENDIF
                        ENDIF
                     ENDIF
                     IF(I_HIND.LT.I_HMAX) THEN
                        TOTAL=TOTAL+
     &                     CALCPATHS_R(IPATH,RHOII,RHOIJ,I_V,I_HMAX,
     &                    I_P,NWHTAY,I_HIND+1,NEL,IPATH2,IP2,RHOCUML*RH,
     &                     RHOII2,I_VMAX,ILOGGING)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
         CALCPATHS_RRES=TOTAL
         RETURN
      END


C.. Calculate all possible paths of two vertices, using a contour integral trick
C.. and summing the resulant geometric progression.
C.. We first calculate the roots of (z-rii)(z-rjj)-Rij=0 as z1 and z2
C.. The result is  Rij (z1^P-rii^P)/((z1-rii)(z1-z2))
C..               +Rij (z2^P-rii^P)/((z2-rii)(z2-z1))
C..              = (Rij/(z1-z2))((z1^P-rii^P)/(z1-rii)-(z2^P-rii^P)/(z2-rii))

      REAL*8 FUNCTION CALCPATHS_2(RHOII,RHOIJ,P)
         INTEGER P
         REAL*8 RHOII(0:2),RHOIJ(0:2,0:2)
         REAL*8 Z1,Z2,ZD,RIIP
         LOGICAL LSAME2
         IF(RHOIJ(0,1).EQ.0.D0) THEN
            CALCPATHS_2=0.D0
            RETURN
         ENDIF
C.. We've got a quadratic.  Solve easily (Numerical Recipes)
         A=-(RHOII(0)/RHOII(0)+RHOII(1)/RHOII(0))
         B=(RHOII(0)*RHOII(1)-RHOIJ(0,1)**2)/(RHOII(0)**2)
         Q=SQRT(A*A-4*B)
         IF(A.LT.0.D0) Q=-Q
         Q=-0.5*(A+Q)
         Z1=Q
         Z2=B/Q
C.. We must be careful.  If z1 or z2 = 1, then this expression
C.. blows up, so we use the limited version.

         IF(LSAME2(Z1,1.D0)) THEN
            CALCPATHS_2=P
         ELSE
            CALCPATHS_2=(Z1**P-1)/(Z1-1)
         ENDIF
         IF(LSAME2(Z2,1.D0)) THEN
            IF(LSAME2(Z1,1.D0)) THEN
C.. Z1=Z2=1, so we use the residue as pole is order 2
               CALCPATHS_2=P*(P-1)*(RHOIJ(0,1)/RHOII(0)**2)/(2*RIIP)
               RETURN
            ELSE
               CALCPATHS_2=CALCPATHS_2-P
            ENDIF
         ELSE
            CALCPATHS_2=CALCPATHS_2
     &               -(Z2**P-1)/(Z2-1)
         ENDIF
         CALCPATHS_2=CALCPATHS_2*((RHOIJ(0,1)/RHOII(0))**2)
     &                  /(Z1-Z2)
         RETURN
      END

      REAL*8 FUNCTION CP2(RII,RJJ,RIJ,P)
         INTEGER P
         REAL*8 RII,RJJ,RIJ,RTII(0:2),RTIJ(0:2,0:2)
         RTII(0)=RII
         RTII(2)=RII
         RTII(1)=RJJ
         RTIJ(0,1)=RIJ
         RTIJ(1,0)=RIJ
         CP2=CALCPATHS_2(RTII,RTIJ,P)
         RETURN
      END

C.. For f[z]=(z^P-rii^P)(2Rijk+Rij(z-rkk)+Rik(z-rjj)/(z-rii)
      REAL*8 FUNCTION CP3F(Z,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
         REAL*8 Z,RII,RJJ,RKK,RIJK,RIJ,RIK
         INTEGER P
         LOGICAL LSAME2
         IF(LSAME2(Z,1.D0)) THEN
C.. deal with limiting expression
            CP3F=P*(2*RIJK+RIJ*(Z-RKK)+RIK*(Z-RJJ))/Z
         ELSE
            CP3F=(Z**P-1)*(2*RIJK+RIJ*(Z-RKK)+RIK*(Z-RJJ))/(Z-1)
         ENDIF
C         WRITE(10,*) "CP3F=",CP3F
         RETURN
      END
C            CP3F=((Z/RII)**P-1)*(2*RIJK+RIJ*(Z-RKK)+RIK*(Z-RJJ))/(Z-RII)
 

      REAL*8 FUNCTION DCP3F(Z1,Z2,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
         REAL*8 Z1,Z2,RII,RJJ,RKK,RIJK,RIJ,RIK
         INTEGER P
         LOGICAL LSAME2
         REAL*8 A,B,C,AP,BP,CP
C.. d(AB/C)/dz = (A'BC+AB'C-ABC')/CC
         A=(Z1**P-1)/(Z1-1)
         B=RIJK+RIJ*(Z1-RKK)+RIK*(Z1-RJJ)
         C=Z1-Z2
         AP=-A/(Z1-1)+P*A/Z1
         BP=RIJ+RIK
         CP=1
         DCP3F=(AP*B*C+A*BP*C-A*B*CP)/(C*C)
         DCP3F=DCP3F
         RETURN
      END


C.. For f[z]=(z^P-rii^P)(2Rijk+Rij(z-rkk)+Rik(z-rjj)/(z-rii)
      COMPLEX*16 FUNCTION CP3FC(Z,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
         REAL*8 RII,RJJ,RKK,RIJK,RIJ,RIK
         COMPLEX*16 Z
         INTEGER P
         CP3FC=(Z**P-1)*(2*RIJK+RIJ*(Z-RKK)+RIK*(Z-RJJ))/(Z-1)
         RETURN
      END

C.. Calculate all possible paths of 3 vertices, using a contour integral trick
C.. and summing the resulant geometric progression.
C.. We first calculate the roots of 
C.. (z-rii)(z-rjj)(z-rkk)-(2Rijk+Rij(z-rkk)+Rik(z-rjj)+Rjk(z-rii))=0 
C.. as z1, z2 and z3
C.. For f[z]=(z^P-rii^P)(2Rijk+Rij(z-rkk)+Rik(z-rjj)/(z-rii)
C.. The result is
C.. f[z1]/((z1-z2)(z1-z3))+f[z2]/((z2-z1)(z2-z3))+f[z3]/((z3-z1)(z3-z2))
C.. - CALCPATHS_2(i,j) - CALCPATHS_2(i,k)
C.. Although we want the answer /rii^P, so we divide that out
      REAL*8 FUNCTION CP3_SOLVE(A,B,C,CP3TF,DCP3TF,CP3TFC,P,
     &   RII,RJJ,RKK,RIJ,RIK,RIJK)
         IMPLICIT NONE
         REAL*8 CP3TF,DCP3TF 
         INTEGER P
         COMPLEX*16 CP3TFC
         REAL*8 A,B,C,RA
         REAL*8 RII,RJJ,RKK,RIJ,RIK,RJK,RIJK
         
         COMPLEX*16 ZZ2,ZZ3,ZZ1,ZZ,ZA
         REAL*8 Z1,Z2,Z3
         LOGICAL LCOMPLEX
         LOGICAL LSAME2

         REAL*8 POLY(0:3)
         COMPLEX*16 ZR(3)
         REAL*8 CP2
C.. Deal with a special case
         IF(RJK.EQ.0.D0.AND.ABS((RIJ-RIK)/RIJ).LE.1.D-10) THEN
C..      The situation reduces to 2Aij/(1-2Aij)
            CP3_SOLVE=CP2(RII,RJJ,2.D0*RIJ,P)
            RETURN
         ENDIF
         POLY(3)=1.D0
         POLY(2)=A
         POLY(1)=B
         POLY(0)=C
C         CALL FACTORPOLY(POLY,3,ZR)
         CALL SOLVE_CUBIC(A,B,C,ZZ1,ZZ2,ZZ3,LCOMPLEX)
C            WRITE(23,
C     & "(3G25.16e3)")
C     &   A,B,C
C            WRITE(23,
C     & "((G25.16e3,G25.16e3),(G25.16e3, G25.16e3),(G25.16e3,G25.16e3))")
C     &   ZZ1,ZZ2,ZZ3
C            WRITE(23,
C     & "((G25.16e3,G25.16e3),(G25.16e3, G25.16e3),(G25.16e3,G25.16e3))")
C     &   ZR(1),ZR(2),ZR(3) 
         IF(LCOMPLEX) THEN
C            WRITE(22,
C     &   "(G25.16e3,(G25.16e3, G25.16e3),(G25.16e3,G25.16e3)")
C     &   DBLE(ZZ1),ZZ2,ZZ3
         ELSE
C            WRITE(22,"3G25.16e3")
C     &         DBLE(ZZ1),DBLE(ZZ2),DBLE(ZZ3)
         ENDIF
C.. ZZ1,ZZ2,ZZ3 now contain the roots to the cubic. LCOMPLEX is set if
C.. ZZ2 and ZZ3 are complex
         IF(LCOMPLEX) THEN
            Z1=DBLE(ZZ1)
            ZA=CP3TFC(ZZ2,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((ZZ2-Z1)*(ZZ2-ZZ3))
            ZZ=ZA
C            WRITE(31,"((G25.16e3, G25.16e3))",advance='no') ZA
         
            ZA=CP3TFC(ZZ3,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((ZZ3-Z1)*(ZZ3-ZZ2))
            ZZ=ZZ+ZA
C            WRITE(31,"((G25.16e3, G25.16e3))",advance='no') ZA
            Z2=DBLE(ZZ2)
            Z3=DIMAG(ZZ2)
            CP3_SOLVE=
     &         CP3TF(Z1,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /(Z1*Z1-2*Z2*Z1+Z2*Z2+Z3*Z3)
C            WRITE(31,"(G25.16e3)") CP3_SOLVE 
            CP3_SOLVE=CP3_SOLVE+DBLE(ZZ)
         ELSE
C.. all roots are real
C.. Only one degeneracy can occur, when QRRAT=1, Z1=Z3, so we handle 
C.. that separately, as it has a different set of residues
C               WRITE(10,*) Z1,Z2,Z3
            Z1=DBLE(ZZ1)
            Z2=DBLE(ZZ2)
            Z3=DBLE(ZZ3)
            IF(LSAME2(Z1,Z3)) THEN
               RA=CP3TF(Z2,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((Z2-Z1)*(Z2-Z3))
               CP3_SOLVE=RA
C               WRITE(31,"(G25.16e3)",advance='no') RA
               RA=DCP3TF(Z1,Z2,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
C               WRITE(31,"(G25.16e3)") RA
               CP3_SOLVE=CP3_SOLVE+RA
            ELSE
               RA=CP3TF(Z1,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((Z1-Z2)*(Z1-Z3))
               CP3_SOLVE=RA
C               WRITE(31,"(G25.16e3)",advance='no') RA
               RA=CP3TF(Z2,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((Z2-Z1)*(Z2-Z3))
               CP3_SOLVE=CP3_SOLVE+RA
C               WRITE(31,"(G25.16e3)",advance='no') RA
               RA=CP3TF(Z3,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((Z3-Z1)*(Z3-Z2))
               CP3_SOLVE=CP3_SOLVE+RA
C               WRITE(31,"(G25.16e3)") RA
            ENDIF
         ENDIF
         RETURN
      END
C.. Solve the cubic Z^3 + A Z^2 + B Z + C=0
C.. If the ZZ2 and ZZ3 are complex, LCOMPLEX is set
      SUBROUTINE SOLVE_CUBIC(A,B,C,ZZ1,ZZ2,ZZ3,LCOMPLEX)
         IMPLICIT NONE
         COMPLEX*16 ZZ1,ZZ2,ZZ3
         REAL*8 A,B,C
         LOGICAL LCOMPLEX         
         REAL*8 Q,R,RT,QRRAT,THETAB3,TRQ,AA,BB,Z1,Z2,Z3,TPB3
         TPB3=3.14159265358979323846464D0*2.D0/3.D0
         IF(C.EQ.0.D0) THEN
C.. We've actually got a quadratic as some of our vertices have no self
C..   loops, so C is zero.  Solve easily (Numerical Recipes)
            Q=DSQRT(A*A-4.D0*B)
            IF(A.LT.0.D0) Q=-Q
            Q=-0.5D0*(A+Q)
            ZZ1=0.D0
            ZZ2=Q
            ZZ3=B/Q
            LCOMPLEX=.FALSE.
         ELSE
            Q=(A*A-3.D0*B)/9.D0
            R=(A*(2.D0*A*A-9.D0*B)+27.D0*C)/54.D0
            RT=R*R-Q*Q*Q
            IF(RT.LT.0.D0) THEN
C.. All roots are real
               QRRAT=R/DSQRT(Q*Q*Q)
               THETAB3=DACOS(QRRAT)/3.D0
               TRQ=-2.D0*DSQRT(Q)
               ZZ1=TRQ*DCOS(THETAB3)-A/3.D0
               ZZ2=TRQ*DCOS(THETAB3+TPB3)-A/3.D0
               ZZ3=TRQ*DCOS(THETAB3-TPB3)-A/3.D0 
               LCOMPLEX=.FALSE.
           ELSE
C.. We have a real root, and a complex conjugate pair of roots.
               AA=-((ABS(R)+DSQRT(RT))**(1.D0/3.D0))
               IF(R.LT.0.D0) AA=-AA
               BB=Q/AA
               Z2=AA+BB
               ZZ1=Z2-A/3.D0
C.. Re(ZZ2)
               Z2=-Z2/2.D0-A/3.D0
C.. Im(ZZ2)
               Z3=SQRT(3.D0)*(AA-BB)/2.D0
               ZZ2=DCMPLX(Z2,Z3)
               ZZ3=DCMPLX(Z2,-Z3)
               LCOMPLEX=.TRUE.
            ENDIF
         ENDIF
         RETURN
      END

CC.. This is the old version of CP3_SOLVE to be removed at some point 
      REAL*8 FUNCTION CP3_SOLVE2(A,B,C,CP3TF,DCP3TF,CP3TFC,P,
     &   RII,RJJ,RKK,RIJ,RIK,RIJK)
         REAL*8 A,B,C
         REAL*8 CP3TF,DCP3TF 
         INTEGER P
         COMPLEX*16 CP3TFC
         COMPLEX*16 ZZ2,ZZ3,ZZ
         REAL*8 RII,RJJ,RKK,RIJ,RIK,RJK,RIJK
         REAL*8 Q,Z1,Z2,Z3,R,RT,QRRAT,THETAB3,TPB3,TRQ
         LOGICAL LSAME2
         TPB3=3.14159265358979323846464D0*2.D0/3.D0
         IF(C.EQ.0.D0) THEN
C.. We've actually got a quadratic as some of our vertices have no self
C..   loops, so C is zero.  Solve easily (Numerical Recipes)
            Q=SQRT(A*A-4*B)
            IF(A.LT.0.D0) THEN Q=-Q
            Q=-0.5*(A+Q)
            Z1=0.D0
            Z2=Q
            Z3=B/Q
            CP3_SOLVE2=
     &      CP3TF(Z1,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((Z1-Z2)*(Z1-Z3))
     &     +CP3TF(Z2,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((Z2-Z1)*(Z2-Z3))
     &     +CP3TF(Z3,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((Z3-Z1)*(Z3-Z2))
         ELSE
            Q=(A*A-3.D0*B)/9.D0
            R=(A*(2.D0*A*A-9.D0*B)+27.D0*C)/54.D0
            RT=R*R-Q*Q*Q
C            WRITE(10,*) R*R,RT,RT.LT.0.D0,SQRT(RT),SQRT(R*R-Q*Q*Q)
            IF(RT.LT.0.D0) THEN
C.. All roots are real
               QRRAT=R/SQRT(Q*Q*Q)
               THETAB3=ACOS(QRRAT)/3.D0
               TRQ=-2.D0*SQRT(Q)
               Z1=TRQ*COS(THETAB3)-A/3.D0
               Z2=TRQ*COS(THETAB3+TPB3)-A/3.D0
               Z3=TRQ*COS(THETAB3-TPB3)-A/3.D0
C.. Only one degeneracy can occur, when QRRAT=1, Z1=Z3, so we handle 
C.. that separately, as it has a different set of residues
C               WRITE(10,*) Z1,Z2,Z3
               IF(LSAME2(Z1,Z3)) THEN
                  CP3_SOLVE2=
     &      CP3TF(Z2,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((Z2-Z1)*(Z2-Z3))
     &     +DCP3TF(Z1,Z2,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
 
               ELSE
C.. Now just calculate the answer
                  CP3_SOLVE2=
     &      CP3TF(Z1,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((Z1-Z2)*(Z1-Z3))
     &     +CP3TF(Z2,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((Z2-Z1)*(Z2-Z3))
     &     +CP3TF(Z3,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((Z3-Z1)*(Z3-Z2))
               ENDIF
           ELSE

C.. We have a real root, and a complex conjugate pair of roots.
C.. Thankfully, we can just take the real part of these, and the answer
C.. will come out nicely
               AA=-((ABS(R)+SQRT(RT))**(1/3.D0))
               IF(R.LT.0.D0) AA=-AA
               BB=Q/AA
               Z2=AA+BB
               Z1=Z2-A/3.D0
C.. Re(ZZ2)
               Z2=-Z2/2.D0-A/3.D0
C.. Im(ZZ2)
               Z3=SQRT(3.D0)*(AA-BB)/2.D0
               ZZ2=DCMPLX(Z2,Z3)
               ZZ3=DCMPLX(Z2,-Z3)
C.. Now just calculate the answer
               ZZ=CP3TFC(ZZ2,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((ZZ2-Z1)*(ZZ2-ZZ3))
     &           +CP3TFC(ZZ3,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /((ZZ3-Z1)*(ZZ3-ZZ2))
               
               CP3_SOLVE2=
     &      CP3TF(Z1,P,RII,RJJ,RKK,RIJK,RIJ,RIK)
     &               /(Z1*Z1-2*Z2*Z1+Z2*Z2+Z3*Z3)
     &     +DBLE(ZZ)
C               WRITE(10,*) AA,BB,Z1,ZZ2,ZZ3
            ENDIF
         ENDIF
         RETURN
      END

      REAL*8 FUNCTION CALCPATHS_3(RHOII,RHOIJ,P)
         INTEGER P
         REAL*8 RHOII(0:3),RHOIJ(0:3,0:3)
         REAL*8 CP3_SOLVE,CP3CI
         REAL*8 RIJ,RIK,RIJK,RJK,A,B,C,T,TT
         REAL*8 RII,RJJ,RKK
         EXTERNAL CP3F,DCP3F
         EXTERNAL CP3FC
         LOGICAL LSAME2
         RII=RHOII(0)/RHOII(0)
         RJJ=RHOII(1)/RHOII(0)
         RKK=RHOII(2)/RHOII(0)
         RIJ=(RHOIJ(0,1)/RHOII(0))**2
         RIK=(RHOIJ(0,2)/RHOII(0))**2
         RIJK=RHOIJ(0,1)*RHOIJ(1,2)*RHOIJ(2,0)/(RHOII(0)**3)
         RJK=(RHOIJ(1,2)/RHOII(0))**2
         A=-(RII+RJJ+RKK)
         B=(RII*RJJ+RII*RKK+RJJ*RKK)-(RIJ+RIK+RJK)
         C=-RII*RJJ*RKK-2*RIJK+RIJ*RKK+RIK*RJJ+RJK*RII
         T=CP3_SOLVE(A,B,C,CP3F,DCP3F,CP3FC,P,
     &               RII,RJJ,RKK,
     &               RIJ,RIK,RIJK)
      IF(.FALSE.) THEN
C.. A breakdown of the subloops
         B=(RII*RJJ+RII*RKK+RJJ*RKK)
         C=-RII*RJJ*RKK-2*RIJK
         T1=CP3_SOLVE(A,B,C,CP3F,DCP3F,CP3FC,P,
     &               RII,RJJ,RKK,
     &               0.D0,0.D0,RIJK)
         B=(RII*RJJ+RII*RKK+RJJ*RKK)-(RIJ+RJK)
         C=-RII*RJJ*RKK+RIJ*RKK+RJK*RII
         T2=CP3_SOLVE(A,B,C,CP3F,DCP3F,CP3FC,P,
     &               RII,RJJ,RKK,
     &               RIJ,0.D0,0.D0)
         B=(RII*RJJ+RII*RKK+RJJ*RKK)-(RIK+RJK)
         C=-RII*RJJ*RKK+RIK*RJJ+RJK*RII
         T3=CP3_SOLVE(A,B,C,CP3F,DCP3F,CP3FC,P,
     &               RII,RJJ,RKK,
     &               0.D0,RIK,0.D0)
         B=(RII*RJJ+RII*RKK+RJJ*RKK)-(RIJ+RIK)
         C=-RII*RJJ*RKK-2*RIJK+RIJ*RKK+RIK*RJJ
         T4=CP3_SOLVE(A,B,C,CP3F,DCP3F,CP3FC,P,
     &               RII,RJJ,RKK,
     &               RIJ,RIK,0.D0)
     
      ENDIF    
         T5=0.D0
         T6=0.D0
         IF(RHOII(0).NE.0.D0) THEN
            IF(RHOIJ(0,1).NE.0.D0)
     &            T5=CP2(RHOII(0),RHOII(1),RHOIJ(0,1),P)
            IF(RHOIJ(0,2).NE.0.D0)
     &            T6=CP2(RHOII(0),RHOII(2),RHOIJ(0,2),P)
         ENDIF
         CALCPATHS_3=T-T5-T6
C         T0=CALCPATHS_3_(RHOII,RHOIJ,P)
C         T7=CP3CI3(RHOII,RHOIJ,P)
C         T8=CP3CI2(RHOII,RHOIJ,P)
C         WRITE(21,"(13E25.16)")
C     &     T,T1,T2,T3,T4,T2-T5,T3-T6,T4-T5-T6,T5,T6,CALCPATHS_3,T0,T7
C       WRITE(21,"(8G25.16e3)")
C    & CALCPATHS_3,T8-T6-T5,(T-T8)/T,T,T8,T7,T5,T6
C         CALCPATHS_3=T8-T6-T5

C         IF(ABS((T-T8)/T).GE.1.D-3) STOP "CI and CU differ"
C         
C            T0=CALCPATHS_3_(RHOII,RHOIJ,P)
C            STOP 'CP3'
C         ENDIF
         RETURN
      END

      REAL*8 FUNCTION CALCPATHS_3_(RHOII,RHOIJ,P)
         INTEGER P
         REAL*8 RHOII(0:3),RHOIJ(0:3,0:3)
         REAL*8 Z1,Z2,Z3,Q,R,THETAB3,A,B,C,QRRAT,TRQ
         REAL*8 TPB3,AA,BB,RT
         REAL*8 RIJ,RIK,RIJK,RTII(0:2),RTIJ(0:2,0:2)
         COMPLEX*16 ZZ2,ZZ3,ZZ,CP3FC
         LOGICAL LSAME2
         TPB3=3.14159265358979323846464D0*2.D0/3.D0
C.. Solve the cubic Z^3+A Z^2+B Z+C = 0
C.. (Numerical Recipes in FORTRAN 2nd Ed p179
         RIJ=RHOIJ(0,1)**2
         RIK=RHOIJ(0,2)**2
         RIJK=RHOIJ(0,1)*RHOIJ(1,2)*RHOIJ(2,0)
         A=-(RHOII(0)+RHOII(1)+RHOII(2))
         B=RHOII(0)*RHOII(1)+RHOII(0)*RHOII(2)+RHOII(1)*RHOII(2)
     &    -(RIJ+RIK+RHOIJ(1,2)**2)
         C=-RHOII(0)*RHOII(1)*RHOII(2)-2*RIJK
     &         +RIJ*RHOII(2)+RIK*RHOII(1)+(RHOIJ(1,2)**2)*RHOII(0)
         IF(C.EQ.0.D0) THEN
C.. We've actually got a quadratic as some of our vertices have no self
C..   loops, so C is zero.  Solve easily (Numerical Recipes)
            Q=SQRT(A*A-4*B)
            IF(A.LT.0.D0) THEN Q=-Q
            Q=-0.5*(A+Q)
            Z1=0.D0
            Z2=Q
            Z3=B/Q
               CALCPATHS_3_=
     &      CP3F(Z1,P,RHOII(0),RHOII(1),RHOII(2),RIJK,RIJ,RIK)
     &               /((Z1-Z2)*(Z1-Z3))
     &     +CP3F(Z2,P,RHOII(0),RHOII(1),RHOII(2),RIJK,RIJ,RIK)
     &               /((Z2-Z1)*(Z2-Z3))
     &     +CP3F(Z3,P,RHOII(0),RHOII(1),RHOII(2),RIJK,RIJ,RIK)
     &               /((Z3-Z1)*(Z3-Z2))
         ELSE
            Q=(A*A-3.D0*B)/9.D0
            R=(A*(2.D0*A*A-9.D0*B)+27.D0*C)/54.D0
            RT=R*R-Q*Q*Q
C            WRITE(10,*) R*R,RT,RT.LT.0.D0,SQRT(RT),SQRT(R*R-Q*Q*Q)
            IF(RT.LT.0.D0) THEN
C.. All roots are real
               QRRAT=R/SQRT(Q*Q*Q)
               THETAB3=ACOS(QRRAT)/3.D0
               TRQ=-2.D0*SQRT(Q)
               Z1=TRQ*COS(THETAB3)-A/3.D0
               Z2=TRQ*COS(THETAB3+TPB3)-A/3.D0
               Z3=TRQ*COS(THETAB3-TPB3)-A/3.D0
C.. Only one degeneracy can occur, when QRRAT=1, Z1=Z3, so we handle 
C.. that separately, as it has a different set of residues
C               WRITE(10,*) Z1,Z2,Z3
               IF(LSAME2(Z1,Z3)) THEN
                  CALCPATHS_3_=
     &      CP3F(Z2,P,RHOII(0),RHOII(1),RHOII(2),RIJK,RIJ,RIK)
     &               /((Z2-Z1)*(Z2-Z3))
     &     +DCP3F(Z1,Z2,P,RHOII(0),RHOII(1),RHOII(2),RIJK,RIJ,RIK)
 
               ELSE
C.. Now just calculate the answer
                  CALCPATHS_3_=
     &      CP3F(Z1,P,RHOII(0),RHOII(1),RHOII(2),RIJK,RIJ,RIK)
     &               /((Z1-Z2)*(Z1-Z3))
     &     +CP3F(Z2,P,RHOII(0),RHOII(1),RHOII(2),RIJK,RIJ,RIK)
     &               /((Z2-Z1)*(Z2-Z3))
     &     +CP3F(Z3,P,RHOII(0),RHOII(1),RHOII(2),RIJK,RIJ,RIK)
     &               /((Z3-Z1)*(Z3-Z2))
               ENDIF
           ELSE

C.. We have a real root, and a complex conjugate pair of roots.
C.. Thankfully, we can just take the real part of these, and the answer
C.. will come out nicely
               AA=-((ABS(R)+SQRT(RT))**(1/3.D0))
               IF(R.LT.0.D0) AA=-AA
               BB=Q/AA
               Z2=AA+BB
               Z1=Z2-A/3.D0
C.. Re(ZZ2)
               Z2=-Z2/2.D0-A/3.D0
C.. Im(ZZ2)
               Z3=SQRT(3.D0)*(AA-BB)/2.D0
               ZZ2=DCMPLX(Z2,Z3)
               ZZ3=DCMPLX(Z2,-Z3)
C.. Now just calculate the answer
               ZZ=CP3FC(ZZ2,P,RHOII(0),RHOII(1),RHOII(2),RIJK,RIJ,RIK)
     &               /((ZZ2-Z1)*(ZZ2-ZZ3))
     &           +CP3FC(ZZ3,P,RHOII(0),RHOII(1),RHOII(2),RIJK,RIJ,RIK)
     &               /((ZZ3-Z1)*(ZZ3-ZZ2))
               
               CALCPATHS_3_=
     &      CP3F(Z1,P,RHOII(0),RHOII(1),RHOII(2),RIJK,RIJ,RIK)
     &               /(Z1*Z1-2*Z2*Z1+Z2*Z2+Z3*Z3)
     &     +DBLE(ZZ)
C               WRITE(10,*) AA,BB,Z1,ZZ2,ZZ3
            ENDIF
         ENDIF
         RTII(0)=RHOII(0)
         RTII(1)=RHOII(1)
         RTII(2)=RHOII(0)
         RTIJ(0,1)=RHOIJ(0,1)
         RTIJ(1,0)=RHOIJ(1,0)
C         WRITE(10,*) CALCPATHS_3
         IF(RTII(0).NE.0.D0) 
     &         CALCPATHS_3_=CALCPATHS_3_-CALCPATHS_2(RTII,RTIJ,P)
         RTII(1)=RHOII(2)
         RTIJ(0,1)=RHOIJ(0,2)
         RTIJ(1,0)=RHOIJ(2,0)
C         WRITE(10,*) CALCPATHS_3
         IF(RTII(0).NE.0.D0) 
     &         CALCPATHS_3_=CALCPATHS_3_-CALCPATHS_2(RTII,RTIJ,P)
C         WRITE(10,*) CALCPATHS_3
      END

      REAL*8 FUNCTION CP3CI3(RHOII,RHOIJ,I_P)
         IMPLICIT NONE
         INTEGER I_P
         REAL*8 RHOII(0:3),RHOIJ(0:3,0:3)
         EXTERNAL CP3CIZSFUN2
         REAL*8 SUM,SUMI
         REAL*8 XIS(3,4),RIS(4)
         REAL*8 RII,RJJ,RKK,RIJ,RIK,RIJK,RJK
         REAL*8 POLYNUM(0:4),POLYDEN(0:4)
         COMPLEX*16 ZROOTS(3)
         RII=RHOII(0)/RHOII(0)
         RJJ=RHOII(1)/RHOII(0)
         RKK=RHOII(2)/RHOII(0)
         RIJ=(RHOIJ(0,1)/RHOII(0))**2
         RIK=(RHOIJ(0,2)/RHOII(0))**2
         RIJK=RHOIJ(0,1)*RHOIJ(1,2)*RHOIJ(2,0)/(RHOII(0)**3)
         RJK=(RHOIJ(1,2)/RHOII(0))**2
         XIS(1,1)=RII
         XIS(2,1)=RJJ
         XIS(3,1)=RKK
         XIS(1,2)=RII
         XIS(2,2)=RJJ
         XIS(1,3)=RII
         XIS(2,3)=RKK
         XIS(1,4)=RJJ
         XIS(2,4)=RKK
         RIS(1)=RIJK
         RIS(2)=RIJ
         RIS(3)=RIK
         RIS(4)=RJK
         POLYDEN(4)=1.D0
         POLYDEN(3)=-(RII+RJJ+RKK)
         POLYDEN(2)=(RII*RJJ+RII*RKK+RJJ*RKK)-(RIJ+RIK+RJK)
         POLYDEN(1)=-RII*RJJ*RKK-2*RIJK+RIJ*RKK+RIK*RJJ+RJK*RII
         POLYDEN(0)=DBLE(3)
         POLYNUM(2)=RIJ+RIK
         POLYNUM(1)=2*RIJK-RIJ*RKK-RIK*RJJ
         POLYNUM(0)=DBLE(1)        
 
C         OPEN(30,FILE="CIPATH3",STATUS="UNKNOWN")

C         WRITE(30,"(3G25.16e3)",advance='no') POLYDEN(3),POLYDEN(2),POLYDEN(1)
C         WRITE(30,"(G25.16e3)",advance='no') POLYNUM(1)/POLYNUM(2)
         CALL FACTORPOLY(POLYDEN(1),3,ZROOTS)
C         WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZROOTS(1)
C         WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZROOTS(2)
C         WRITE(30,"((G25.16e3, G25.16e3))") ZROOTS(3)
C         WRITE(30,"(4G25.16e3)") RIS(1),RIS(2),RIS(3),RIS(4)
         CALL ZSUM3(POLYNUM,POLYDEN,CP3CIZSFUN2,1000,I_P,3,ZROOTS,
     &            .TRUE.,SUM,SUMI)
C         CLOSE(30)
         CP3CI3=SUM
      END

      REAL*8 FUNCTION CP3CI2(RHOII,RHOIJ,I_P)
         IMPLICIT NONE
         INTEGER I_P
         REAL*8 RHOII(0:3),RHOIJ(0:3,0:3)
         EXTERNAL CP3CIZSFUN
         REAL*8 SUM,SUMI
         REAL*8 XIS(3,4),RIS(4)
         REAL*8 RII,RJJ,RKK,RIJ,RIK,RIJK,RJK
         REAL*8 POLY(0:3)
         COMPLEX*16 ZROOTS(3)
         RII=RHOII(0)/RHOII(0)
         RJJ=RHOII(1)/RHOII(0)
         RKK=RHOII(2)/RHOII(0)
         RIJ=(RHOIJ(0,1)/RHOII(0))**2
         RIK=(RHOIJ(0,2)/RHOII(0))**2
         RIJK=RHOIJ(0,1)*RHOIJ(1,2)*RHOIJ(2,0)/(RHOII(0)**3)
         RJK=(RHOIJ(1,2)/RHOII(0))**2
         XIS(1,1)=RII
         XIS(2,1)=RJJ
         XIS(3,1)=RKK
         XIS(1,2)=RII
         XIS(2,2)=RJJ
         XIS(1,3)=RII
         XIS(2,3)=RKK
         XIS(1,4)=RJJ
         XIS(2,4)=RKK
         RIS(1)=RIJK
         RIS(2)=RIJ
         RIS(3)=RIK
         RIS(4)=RJK
         POLY(3)=1.D0
         POLY(2)=-(RII+RJJ+RKK)
         POLY(1)=(RII*RJJ+RII*RKK+RJJ*RKK)-(RIJ+RIK+RJK)
         POLY(0)=-RII*RJJ*RKK-2*RIJK+RIJ*RKK+RIK*RJJ+RJK*RII
         
C         OPEN(30,FILE="CIPATH2",STATUS="UNKNOWN")

C         WRITE(30,"(3G25.16e3)",advance='no') POLY(2),POLY(1),POLY(0)
         CALL FACTORPOLY(POLY,3,ZROOTS)
C         WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZROOTS(1)
C         WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZROOTS(2)
C         WRITE(30,"((G25.16e3, G25.16e3))") ZROOTS(3)
C         WRITE(30,"(4G25.16e3)") RIS(1),RIS(2),RIS(3),RIS(4)
C         ZROOTS(1)=0.5D0
C         CALL ZSUM3(XIS,RIS,CP3CIZSFUN,1000,I_P,1,ZROOTS,
C     &            SUM,SUMI)
         CALL ZSUM3(XIS,RIS,CP3CIZSFUN,1000,I_P,3,ZROOTS,
     &            .FALSE.,SUM,SUMI)
C         CLOSE(30)
         CP3CI2=SUM
      END

      REAL*8 FUNCTION CP3CI(RHOII,RHOIJ,I_P)
         IMPLICIT NONE
         INTEGER I_P
         REAL*8 RHOII(0:3),RHOIJ(0:3,0:3)
         EXTERNAL CP3CIZSFUN
         REAL*8 SUM,SUMI
         REAL*8 XIS(3,4),RIS(4)
         REAL*8 RII,RJJ,RKK,RIJ,RIK,RIJK,RJK
         RII=RHOII(0)/RHOII(0)
         RJJ=RHOII(1)/RHOII(0)
         RKK=RHOII(2)/RHOII(0)
         RIJ=(RHOIJ(0,1)/RHOII(0))**2
         RIK=(RHOIJ(0,2)/RHOII(0))**2
         RIJK=RHOIJ(0,1)*RHOIJ(1,2)*RHOIJ(2,0)/(RHOII(0)**3)
         RJK=(RHOIJ(1,2)/RHOII(0))**2
         XIS(1,1)=RII
         XIS(2,1)=RJJ
         XIS(3,1)=RKK
         XIS(1,2)=RII
         XIS(2,2)=RJJ
         XIS(1,3)=RII
         XIS(2,3)=RKK
         XIS(1,4)=RJJ
         XIS(2,4)=RKK
         RIS(1)=RIJK
         RIS(2)=RIJ
         RIS(3)=RIK
         RIS(4)=RJK
C         OPEN(30,FILE="CIPATH",STATUS="UNKNOWN")
         
          CALL ZSUM2(XIS,RIS,CP3CIZSFUN,(1.D0,0.D0),I_P,10000,1D-4,
     &           .FALSE.,SUM,SUMI)
C         CLOSE(30)
         CP3CI=SUM
      END

      SUBROUTINE CP3CIZSFUN(XIS,RIS,P,Z,R,THETA,ZMZ0,ZRET)
         IMPLICIT NONE
         REAL*8 XIS(3,4), RIS(4)
         COMPLEX*16 Z,ZRET,ZT,ZS,ZZT,ZMZ0
         REAL*8 R,THETA,P
C         WRITE(30,"((G25.16e3, G25.16e3),2G25.16e3)",advance='no') Z,R,THETA
C.. AIJK
         CALL ZPRIM1(XIS(1,1),3,Z,ZS)
         ZZT=2*ZS*RIS(1)
         ZT=ZZT
C         WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZZT
C..AIJ
         CALL ZPRIM1(XIS(1,2),2,Z,ZS)
         ZZT=ZS*RIS(2)
         ZT=ZT+ZZT
C         WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZZT
C..AIK
         CALL ZPRIM1(XIS(1,3),2,Z,ZS)
         ZZT=ZS*RIS(3)
         ZT=ZT+ZZT
C         WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZZT
         ZRET=ZT
C.. AJK
         CALL ZPRIM1(XIS(1,4),2,Z,ZS)
         ZZT=ZS*RIS(4)
         ZT=ZT+ZZT
C         WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZZT
         ZRET=ZRET/(1-ZT)
         CALL ZPRIM2(P,Z,R,THETA,ZMZ0,ZRET)
C         ZRET=Z*(R*EXP((0.D0,1.D0)*THETA))*(0.D0,1.D0)/(Z-0.5D0)
C          WRITE(30,"((G25.16e3, G25.16e3))") ZRET
      END

      SUBROUTINE CP3CIZSFUN2(POLYNUM,POLYDEN,P,Z,R,THETA,ZMZ0,ZRET)
         IMPLICIT NONE
         REAL*8 POLYNUM(0:*)
         COMPLEX*16 POLYDEN(0:*)
         COMPLEX*16 Z,ZRET,ZT,ZS,ZZT,ZMZ0
         REAL*8 R,THETA,P
C        WRITE(30,"((G25.16e3, G25.16e3),2G25.16e3)",advance='no') Z,R,THETA

         CALL POLYVALUE(POLYNUM(1),NINT(POLYNUM(0)),Z,ZRET)
C.. This is the auxiliary taylor expanded poly around z0
         CALL POLYVALUEC(POLYDEN(1),NINT(DBLE(POLYDEN(0))),ZMZ0,ZT)
C         WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZRET
C         WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZT
         ZRET=ZRET/ZT
C         WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZRET
         CALL ZPRIM2(P,Z,R,THETA,ZMZ0,ZRET)
C         ZRET=Z*(R*EXP((0.D0,1.D0)*THETA))*(0.D0,1.D0)/(Z-0.5D0)
C          WRITE(30,"((G25.16e3, G25.16e3))") ZRET
      END

      REAL*8 FUNCTION CALCPATHS_4(RHOII,RHOIJ,I_P)
         IMPLICIT NONE
         INTEGER I_P,I
         REAL*8 RHOII(0:4),RHOIJ(0:4,0:4)
         REAL*8 RHOII2(0:2),RHOIJ2(0:2,0:2),A,B,C
         EXTERNAL CP3CIZSFUN2,DCP3F,CP3FC,CP3F
         REAL*8 SUMR,SUMI
         REAL*8 RII,RJJ,RKK,RLL,RIJ,RIK,RIL,RJK,RJL,RKL
         REAL*8 RIJK,RIJL,RIKL,RJKL,RIJKL,RIJLK,RIKJL
         REAL*8 POLYNUM(0:5),POLYDEN(0:5)
         REAL*8 T,T2,CP3_SOLVE,CP2
         COMPLEX*16 ZROOTS(4)
         LOGICAL TLOG
         TLOG=.FALSE.
         RII=RHOII(0)/RHOII(0)
         RJJ=RHOII(1)/RHOII(0)
         RKK=RHOII(2)/RHOII(0)
         RLL=RHOII(3)/RHOII(0)
         RIJ=(RHOIJ(0,1)/RHOII(0))**2
         RIK=(RHOIJ(0,2)/RHOII(0))**2
         RIL=(RHOIJ(0,3)/RHOII(0))**2
         RJK=(RHOIJ(1,2)/RHOII(0))**2
         RJL=(RHOIJ(1,3)/RHOII(0))**2
         RKL=(RHOIJ(2,3)/RHOII(0))**2
         RIJK=RHOIJ(0,1)*RHOIJ(1,2)*RHOIJ(2,0)/(RHOII(0)**3)
         RIJL=RHOIJ(0,1)*RHOIJ(1,3)*RHOIJ(3,0)/(RHOII(0)**3)
         RIKL=RHOIJ(0,2)*RHOIJ(2,3)*RHOIJ(3,0)/(RHOII(0)**3)
         RJKL=RHOIJ(1,2)*RHOIJ(2,3)*RHOIJ(3,1)/(RHOII(0)**3)
         RIJKL=RHOIJ(0,1)*RHOIJ(1,2)*RHOIJ(2,3)*RHOIJ(3,0)
     &         /(RHOII(0)**4)
         RIJLK=RHOIJ(0,1)*RHOIJ(1,3)*RHOIJ(3,2)*RHOIJ(2,0)
     &         /(RHOII(0)**4)
         RIKJL=RHOIJ(0,2)*RHOIJ(2,1)*RHOIJ(1,3)*RHOIJ(3,0)
     &         /(RHOII(0)**4)
         POLYDEN(5)=1.D0
         POLYDEN(4)=-(RII+RJJ+RKK+RLL)
         POLYDEN(3)=(RII*(RJJ+RKK+RLL)+RJJ*(RKK+RLL)+RKK*RLL)
     &               -(RIJ+RIK+RIL+RJK+RJL+RKL)
         POLYDEN(2)=-2*(RIJK+RIJL+RIKL+RJKL)+RII*(RJK+RJL+RKL)
     &         +RJJ*(RIK+RIL+RKL)+RKK*(RIJ+RIL+RJL)+RLL*(RIJ+RIK+RJK)
     &         -RII*(RJJ*RKK+RJJ*RLL+RKK*RLL)-RJJ*RKK*RLL
         POLYDEN(1)=-2*(RIJKL+RIJLK+RIKJL)
     &               +2*(RIJK*RLL+RIJL*RKK+RIKL*RKK+RJKL*RII)
     &               +(RIJ*RKL+RIK*RJL+RIL*RJK)
     &               -(RIJ*RKK*RLL+RIK*RJJ*RLL+RIL*RJJ*RKK
     &                +RJK*RII*RLL+RJL*RII*RKK+RKL*RII*RJJ)
     &               +RII*RJJ*RKK*RLL
         POLYDEN(0)=4.D0
         POLYNUM(5)=1.D0
         POLYNUM(4)=POLYDEN(4)
         POLYNUM(3)=(RII*(RJJ+RKK+RLL)+RJJ*(RKK+RLL)+RKK*RLL)
     &               -(RJK+RJL+RKL)
         POLYNUM(2)=-2*RJKL+RJK*(RII+RLL)+RJL*(RII+RKK)+RKL*(RII+RJJ)
     &            -RII*RJJ*(RKK+RLL)-(RII+RJJ)*RKK*RLL
         POLYNUM(1)=RII*((2.D0*RJKL-(RJK*RLL+RJL*RKK+RKL*RJJ))
     &                     +RJJ*RKK*RLL)
         POLYNUM(0)=4.D0 
         WRITE(34,"(8G25.16e3)") RII,RJJ,RKK,RLL,RJK,RJL,RKL,RJKL
         CALL FACTORPOLY(POLYDEN(1),4,ZROOTS)
 
         IF(TLOG) THEN
            OPEN(30,FILE="CIPATH4",STATUS="UNKNOWN")
            WRITE(30,"(4G25.16e3)",advance='no') (POLYDEN(I),I=4,1,-1)
            WRITE(30,"(4G25.16e3)",advance='no') (POLYNUM(I),I=4,1,-1)
            WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZROOTS(1)
            WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZROOTS(2)
            WRITE(30,"((G25.16e3, G25.16e3))",advance='no') ZROOTS(3)
            WRITE(30,"((G25.16e3, G25.16e3))") ZROOTS(4)
         ENDIF
         CALL ZSUM3(POLYNUM,POLYDEN,CP3CIZSFUN2,2000,I_P,4,ZROOTS,
     &            .TRUE.,SUMR,SUMI)
         CALCPATHS_4=SUMR
         IF(TLOG) WRITE(30,"(2G25.16e3)",advance='no') SUMR,SUMI 
C.. Now subtract out the 3 and 2 combinations
C.. ijkl=ijkl'-ijk-ijl-ikl-ij-ik-il
C..     =ijkl'-ijk'-ijl'-ikl'+ij'+ik'+il'
C.. (primed include smaller subgraphs)
         IF(RIJ.EQ.0.D0.AND.RIK.EQ.0.D0) THEN
            T=0.D0
         ELSE
            CALL GENCP3ABC(RII,RJJ,RKK,RIJ,RIK,RJK,RIJK,A,B,C)
            T=CP3_SOLVE(A,B,C,CP3F,DCP3F,CP3FC,I_P,
     &               RII,RJJ,RKK,
     &               RIJ,RIK,RIJK)
         ENDIF
         IF(TLOG) WRITE(30,"(G25.16e3)",advance='no') T 
         IF(RIJ.EQ.0.D0.AND.RIL.EQ.0.D0) THEN
            T2=0.D0
         ELSE
            CALL GENCP3ABC(RII,RJJ,RLL,RIJ,RIL,RJL,RIJL,A,B,C)
            T2=CP3_SOLVE(A,B,C,CP3F,DCP3F,CP3FC,I_P,
     &               RII,RJJ,RLL,
     &               RIJ,RIL,RIJL)
         ENDIF
         T=T+T2
         IF(TLOG) WRITE(30,"(G25.16e3)",advance='no') T2 
         IF(RIK.EQ.0.D0.AND.RIL.EQ.0.D0) THEN
            T2=0.D0
         ELSE
            CALL GENCP3ABC(RII,RKK,RLL,RIK,RIL,RKL,RIKL,A,B,C)
            T2=CP3_SOLVE(A,B,C,CP3F,DCP3F,CP3FC,I_P,
     &               RII,RKK,RLL,
     &               RIK,RIL,RIKL)
         ENDIF
         T=T+T2
         IF(TLOG) WRITE(30,"(G25.16e3)",advance='no') T2 
         T2=CP2(RHOII(0),RHOII(1),RHOIJ(0,1),I_P)
         T=T-T2
         IF(TLOG) WRITE(30,"(G25.16e3)",advance='no') T2 
         T2=CP2(RHOII(0),RHOII(2),RHOIJ(0,2),I_P)
         T=T-T2
         IF(TLOG) WRITE(30,"(G25.16e3)",advance='no') T2 
         T2=CP2(RHOII(0),RHOII(3),RHOIJ(0,3),I_P)
         T=T-T2
         IF(TLOG) WRITE(30,"(G25.16e3)",advance='no') T2 

         CALCPATHS_4=CALCPATHS_4-T
         IF(TLOG) WRITE(30,"(G25.16e3)",advance='no') CALCPATHS_4
         IF(TLOG) CLOSE(30)
         RETURN
      END

      SUBROUTINE GENCP3ABC(RII,RJJ,RKK,RIJ,RIK,RJK,RIJK,A,B,C)
         IMPLICIT NONE
         REAL*8 RII,RJJ,RKK,RIJ,RIK,RJK,RIJK,A,B,C
         A=-(RII+RJJ+RKK)
         B=(RII*RJJ+RII*RKK+RJJ*RKK)-(RIJ+RIK+RJK)
         C=-RII*RJJ*RKK-2*RIJK+RIJ*RKK+RIK*RJJ+RJK*RII
         RETURN
      END
