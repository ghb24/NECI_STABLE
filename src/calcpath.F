! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 

      FUNCTION CALCPATHS(IPATH,RHOII,RHOIJ,I_V,
     &      I_HMAX,I_P,NWHTAY,NEL,I_VMAX,ILOGGING)
         use constants
         use global_utilities
         IMPLICIT NONE
         real(dp) CALCPATHS
         INTEGER NEL,I_V,I_P,NWHTAY,I_HMAX,I_VMAX
C.. IPATH contains the vertices we can choose from
C.. IPATH2 contains the path we are constructing (listing the dets)
C.. IP2 contains the path we are constructing as an index into IPATH
         INTEGER IPATH(1:NEL,0:I_V),IPATH2(1:NEL,0:I_HMAX),IP2(0:I_HMAX)
C.. RHOII contains the RHO_II for each vertex.
C.. 12/7/05 - this should now just be the number 0.0_dp - it is not used


C.. RHOII2 is constructed to contain the RHOII for each point in the path
         HElement_t RHOIJ(0:I_V,0:I_V)
         real(dp) RHOII(0:I_V),RHOII2(0:I_HMAX)
         real(dp) DLWDB
         INTEGER ILOGGING
         type(timer), save :: proc_timer
         real(dp) RhoDiag_CPP, CalcPaths_2, CalcPaths_3,
     &      CalcPaths_4, CalcPaths_R
         INTEGER I_HIND
         CHARACTER(*), PARAMETER :: t_r="CALCPATHS"

         CALCPATHS = 0.0_dp !initialisation
         proc_timer%timer_name='CALCPATHS '
         call set_timer(proc_timer,50)
C.. If we're using the contour integral trick to sum all possible graphs with the given vertices
C.. REMEMBER TO MODIFY CALCPATHS_N too!
         IF(I_HMAX.EQ.-1.OR.(I_HMAX.LE.-3.AND.I_HMAX.GE.-8)) THEN
C.. we use a diagonalization version rather than contour etc.
C.. subtracting out lower vertex paths at this level
C            WRITE(6,*) "-----"
                   call stop_all(t_r,"RhoDiag_CPP does not have right   &
     &                   // number of arguments")
!            CALCPATHS=RHODIAG_CPP(RHOIJ,I_P,I_V,I_V+1,.TRUE.,
!     &         0.0_dp, DLWDB)
         ELSEIF(I_HMAX.EQ.-2.OR.I_HMAX.EQ.-6) THEN
C.. we use a diagonalization version rather than contour etc.
C.. we deal with removing lower vertex paths later
C            WRITE(6,*) "-----"
                   call stop_all(t_r,"RhoDiag_CPP does not have right   &
     &                   // number of arguments")
!            CALCPATHS=RHODIAG_CPP(RHOIJ,I_P,I_V,I_V+1,.FALSE.,
!     &         0.0_dp, DLWDB)
         ELSEIF(I_HMAX.EQ.0) THEN
            WRITE(6,*) "I_HMAX=0"
            STOP 'WARNING: Running old disfunctional CALCPATHS I_HMAX=0'
            CALCPATHS=0.0_dp
C            IF(I_V.EQ.2) CALCPATHS=CALCPATHS_2(RHOII,RHOIJ,I_P)
C            IF(I_V.EQ.3) CALCPATHS=CALCPATHS_3(RHOII,RHOIJ,I_P)
C            IF(I_V.EQ.4) CALCPATHS=CALCPATHS_4(RHOII,RHOIJ,I_P)
C            IF(I_V.EQ.3) CALCPATHS=CP3CI2(RHOII,RHOIJ,I_P)
         ELSE
            WRITE(6,*) "I_HMAX>0"
            STOP 'WARNING: Running old disfunctional CALCPATHS I_HMAX>0'
      
C            CALCPATHS=CALCPATHS_R(IPATH,RHOII,RHOIJ,I_V,I_HMAX,I_P,
C     &               NWHTAY,I_HIND,NEL,IPATH2,IP2,1.0_dp,RHOII2,I_VMAX,
C     &               ILOGGING)
         ENDIF
         call halt_timer(proc_timer)
         RETURN
      END 

      FUNCTION CALCPATHS_N(RHOII,RHOIJ,I_V,
     &      I_HMAX,I_P,WEIGHT,DBETA,DLWDB,HIJS,ICLS)
         use constants, only: dp
         use CalcData , only : CALCP_SUB2VSTAR,CALCP_LOGWEIGHT
         use CalcData, only: GraphEpsilon
         use global_utilities
         IMPLICIT NONE
         real(dp) CalcPaths_N
         INTEGER I_V,I_P,I_HMAX,ICLS
C.. IPATH contains the vertices we can choose from
C.. IPATH2 contains the path we are constructing (listing the dets)
C.. IP2 contains the path we are constructing as an index into IPATH
!Removed as no longer used
!  ,IPATH2(1:NEL,0:i_hmax),IP2(0:i_hmax)
C.. RHOII contains the RHO_II for each vertex.
C.. 12/7/05 - this should now just be the number 0.0_dp - it is not used


C.. RHOII2 is constructed to contain the RHOII for each point in the path
         real(dp) RHOII(0:I_V)
! removed as no longer used
!         real(dp) RHOII2(0:i_hmax)
         real(dp) CALCPATHS
         HElement_t HIJS(1:I_V),RHOIJ(0:I_V,0:I_V)
         real(dp) WEIGHT
         real(dp) DBETA
         real(dp) DLWDB
         type(timer), save :: proc_timer
         real(dp) RhoDiag_CPP, CalcPaths_2, CalcPaths_3,
     &      CalcPaths_4, HDiag_CPP
         real(dp) RhoDiag_CPPS2VS,HDiag_CPPS2VS
         character(*), parameter :: t_r="CALCPATHS_N"
         proc_timer%timer_name='CALCPATHS '
         call set_timer(proc_timer,50)
C.. If we're using the contour integral trick to sum all possible graphs with the given vertices
         IF(I_HMAX.EQ.-1.OR.(I_HMAX.LE.-3.AND.I_HMAX.GE.-8).OR.
     &      I_HMAX.LE.-12) THEN
C.. we use a diagonalization version rather than contour etc.
C.. subtracting out lower vertex paths at this level
C            WRITE(6,*) "-----"
            IF(I_HMAX.EQ.-19.OR.I_HMAX.EQ.-20) THEN
!.NOT.(abs(RHOII(0)).gt.0.0_dp)) THEN
C.. We actually diagonalize a matrix of Hij elements rather than a matrix of
C.. RHOIJ elements.  We send beta in from RHOII(1)
               if(calcp_sub2vstar) then
                  CALCPATHS=HDIAG_CPP(RHOIJ,I_P,I_V,I_V+1,.FALSE.,
     &                RHOII(1),DLWDB,HIJS)
                 CALCPATHS=CALCPATHS+HDIAG_CPPS2VS(RHOIJ,I_V,
     &               RHOII(1),DLWDB,HIJS)
               else
                  CALCPATHS=HDIAG_CPP(RHOIJ,I_P,I_V,I_V+1,.TRUE.,
     &            RHOII(1),DLWDB,HIJS)
               endif
            ELSE
               if(calcp_sub2vstar) then
                   call stop_all(t_r,"RhoDiag_CPP does not have right   &
     &                   // number of arguments")
!                  CALCPATHS=RHODIAG_CPP(RHOIJ,I_P,I_V,I_V+1,
!     &            .FALSE.,DBETA,DLWDB,HIJS)
!               CALCPATHS=CALCPATHS+RHODIAG_CPPS2VS(RHOIJ,I_P,I_V,
!     &            DBETA,DLWDB,HIJS)
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
            CALCPATHS=0.0_dp
C            IF(I_V.EQ.2) CALCPATHS=CALCPATHS_2(RHOII,RHOIJ,I_P)
C            IF(I_V.EQ.3) CALCPATHS=CALCPATHS_3(RHOII,RHOIJ,I_P)
C            IF(I_V.EQ.4) CALCPATHS=CALCPATHS_4(RHOII,RHOIJ,I_P)
C            IF(I_V.EQ.3) CALCPATHS=CP3CI2(RHOII,RHOIJ,I_P)
         ELSE
            WRITE(6,*) "I_HMAX>0"
            STOP 'WARNING: Running old disfunctional CALCPATHS I_HMAX>0'
!            CALCPATHS=CALCPATHS_R(IPATH,RHOII,RHOIJ,I_V,I_HMAX,I_P,
!     &               NWHTAY,I_HIND,NEL,IPATH2,IP2,1.0_dp,RHOII2,I_VMAX,
!     &               ILOGGING)
         ENDIF
         IF(.NOT.(abs(CALCPATHS).GE.abs(GraphEpsilon))) THEN
            CALCPATHS=0.0_dp
            DLWDB=0.0_dp
         ENDIF
         CALCPATHS_N=CALCPATHS
         IF(I_V.LE.5) THEN
            IF(WEIGHT.GT.0.0_dp)
     &       CALL CLASSPATHS((CALCPATHS)/WEIGHT,
     &            (DLWDB)/WEIGHT,WEIGHT,RHOIJ,
     &            I_V,ICLS)
         ENDIF
         call halt_timer(proc_timer)
         RETURN
      END
      SUBROUTINE CLASSPATHS(RHOC,HC,WEIGHT,RHOIJ,I_V,ICLASS)
         use constants
         use legacy_data, only: CALCP_N, CALCP_NT, CALCP_SUM, CALCP_HSUM
         IMPLICIT NONE
         INTEGER I_V
         HElement_t RHOIJ(0:I_V,0:I_V)
         real(dp) RHOC,HC,WEIGHT
         INTEGER ILG,ISGN,ICLASS
         INTEGER INCR
C.. allow for subtractions
         IF(I_V.GT.5) RETURN
         INCR=1
         IF(WEIGHT.LT.0.0_dp) INCR=-1
         IF(ICLASS.EQ.0)
     &      CALL CLASSPATH(RHOIJ,I_V,ICLASS)
C            WRITE(56,*) ICLASS,RHOC
C            WRITE(32,*) I_V,ICLASS,RHOC
            ILG=int(LOG(ABS(RHOC))/LOG(10.0_dp),sizeof_int)
            IF(ILG.GE.0) THEN
               ILG=1
            ELSEIF(ILG.LT.-15) THEN
               ILG=15
            ELSE
               ILG=-ILG
            ENDIF
            IF(RHOC.GT.0.0_dp) THEN
               ISGN=1
            ELSEIF(RHOC.LT.0.0_dp) THEN
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
        use constants
         INTEGER I,J,N,ICLASS
         HElement_t RHOIJ(0:N,0:N)
         ICLASS=0
         DO I=0,N-1
            DO J=I+1,N-1
               ICLASS=ICLASS*2
               IF(abs(RHOIJ(I,J)).NE.0.0_dp) THEN
                  ICLASS=ICLASS+1
               ENDIF
            ENDDO
         ENDDO
      END

      SUBROUTINE WRITECLASSPATHS()
         USE Parallel_neci , only : iProcIndex
         use legacy_data, only: CALCP_N, CALCP_NT, CALCP_SUM, CALCP_HSUM
         INTEGER I,J
         IF(iProcIndex.eq.0) THEN
             OPEN(45,FILE='CLASSPATHS',STATUS='UNKNOWN')
             OPEN(46,FILE='CLASSPATHS2',STATUS='UNKNOWN')
             DO I=0,1023
                IF(CALCP_N(I).GT.0) THEN
               WRITE(45,"(I5,I11,4G19.11)") I,CALCP_N(I),CALCP_SUM(I,1),
     &              CALCP_HSUM(I),CALCP_SUM(I,2),CALCP_SUM(I,3)
               WRITE(46,"(I5,I11,2G19.11)") I,CALCP_N(I),CALCP_SUM(I,2),
     &              CALCP_SUM(I,3)
                    DO J=1,15
                        WRITE(46,*) -J,CALCP_NT(I,J,1),CALCP_NT(I,J,2)
                    ENDDO
                ENDIF
            ENDDO
            CLOSE(45)
            CLOSE(46)
         ENDIF
      END


      SUBROUTINE CLEARCLASSPATHS()
         use legacy_data, only: CALCP_N, CALCP_NT, CALCP_SUM, CALCP_HSUM
         CALCP_N=0
         CALCP_N=0
         CALCP_SUM=0
         CALCP_HSUM=0
      END

