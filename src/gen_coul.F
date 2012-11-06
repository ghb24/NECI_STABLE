      SUBROUTINE GEN_COUL(NBASISMAX,
     &   NHG,G1,NMSH,NMAX,FCK,UMAT,ZIA)
      use SystemData, only: BasisFN
      use UMatCache, only: GTID, UMatInd
      use global_utilities
      use constants, only: dp
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      INTEGER :: NMSH,NMAX,NHG,I,J,K,L,ID1,ID2,ID3,ID4,ISPINSKIP,II
      TYPE(BASISFN) G1(NHG)
      INTEGER nBasisMax(5,*)
      type(timer), save :: proc_timer
      real(dp) :: UMAT(*)
C..Cube arrays
      complex(dp) FCK(NMSH,NMSH,NMSH)
      complex(dp) ZIA(-NMSH/2:NMSH/2,NMAX,NMAX)
      real(dp) :: SUM
C..This routine generates *ALL* possible combinations of Coulomb integrals
C..stored in the form (u1 u2 | U | u1' u2') = UMAT(n1 n2 n3 n4)
C..
C..This first call calculates the inner integral 
C..The call to SCOUL calculates the outer integral 
C ==--------------------------------------------------------------------==
      proc_timer%timer_name='GEN_COUL  '
      call set_timer(proc_timer)
C ==--------------------------------------------------------------------==
      OPEN(10,FILE='UMAT',STATUS='UNKNOWN')
      II=0
      ISPINSKIP=NBASISMAX(2,3)
      DO I=1,NHG,ISPINSKIP
        DO J=1,NHG,ISPINSKIP
          DO K=1,NHG,ISPINSKIP
            DO L=1,NHG,ISPINSKIP
              SUM=0.0_dp
C..Original call
C              CALL SLATCOULFOU(G1(1,I),G1(1,J),
C     &        G1(1,K),G1(1,L),NMSH,FCK,NMAX,SUM)
              ID1 = GTID(I)
              ID2 = GTID(J)
              ID3 = GTID(K)
              ID4 = GTID(L)
              IF(UMAT(UMatInd(ID1,ID2,ID3,ID4,0,0)).NE.0.0_dp) GOTO 100
              CALL SLATCOULFOU(G1(I)%k,G1(J)%k,
     &        G1(K)%k,G1(L)%k,NMSH,FCK,NMAX,ZIA,SUM)
c              SUM=SUM*COUPLE 
C..Original call
              UMAT(UMatInd(ID1,ID2,ID3,ID4,0,0))=SUM
!              UMAT(ID1,ID2,ID3,ID4)=SUM
C..Symmetries
!              UMAT(ID3,ID4,ID1,ID2)=SUM
!              UMAT(ID2,ID1,ID4,ID3)=SUM
!              UMAT(ID4,ID3,ID2,ID1)=SUM
              IF(SUM.LE.1.D-10) GOTO 100 
              WRITE(10,'(4I7,F19.9)') ID1,ID2,ID3,ID4,SUM
100           CONTINUE
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CLOSE(10)
      WRITE(6,*) ' !!! FINISHED CALCULATING ALL 2E INTEGRALS !!! '
C ==--------------------------------------------------------------------==
      call halt_timer(proc_timer)
C ==--------------------------------------------------------------------==
      RETURN
      END
C ==------------------------------------------------------------------==
