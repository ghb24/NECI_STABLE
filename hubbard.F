C ===================================================================
      SUBROUTINE CALCUMATHUBREAL(NEL,NBASIS,NBASISMAX,G1,UHUB,
     &            UMAT)
         IMPLICIT NONE
         INTEGER NEL,NBASIS,NBASISMAX(4,3),G1(4,NBASIS)
         REAL*8 UMAT(NBASIS,NBASIS,NBASIS,NBASIS)
         REAL*8 UHUB
         INTEGER I,J,K,IC
         OPEN(10,FILE='UMAT',STATUS='UNKNOWN')
         IC=1
         DO I=1,NBASIS-1,2
            UMAT(I,I+1,I,I+1)=UHUB
            UMAT(I+1,I,I+1,I)=UHUB
            WRITE(10,'(4I7,F19.9)') I,I+1,I,I+1,UHUB
         ENDDO
         CLOSE(10)
         RETURN
      END


      SUBROUTINE HUBKIN(I,J,K,NBASISMAX,BHUB,TTILT,SUM,TREAL)
C returns the non-interacting energy of state with 
C..quatum numbers (i,j,k) for Hubbard model of 
C..lengths (LX,LY,LZ) PBC.
C..Returned in sum
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NBASISMAX(4,5),AX,AY
      INTEGER LX,LY,LZ,DX,DY,DZ
      LOGICAL TTILT,TREAL
      INCLUDE 'cons.inc'
      LX=NBASISMAX(1,2)-NBASISMAX(1,1)+1
      LY=NBASISMAX(2,2)-NBASISMAX(2,1)+1
      LZ=NBASISMAX(3,2)-NBASISMAX(3,1)+1
      AX=NBASISMAX(1,4)
      AY=NBASISMAX(2,4)
      IF(TREAL) THEN
         SUM=0.D0
      ELSE
       IF(TTILT) THEN
CCCC.. NBASISMAX goes from -NMAXX+1 to MAXX so LX=2MAXX
         LX=NBASISMAX(1,5)
         LX=(LX*(AX*AX+AY*AY))
         SUM=COS(2*PI*(I*AX+J*AY)/LX)+COS(2*PI*(I*AY-J*AX)/LX)
       ELSE
         IF(LX.GT.1) SUM=COS(2*PI*I/(LX))
         IF(LY.GT.1) SUM=SUM+COS(2*PI*J/(LY))
         IF(LZ.GT.1) SUM=SUM+COS(2*PI*K/(LZ))
       ENDIF
       SUM=SUM*2.D0*BHUB
      ENDIF
      RETURN
      END

      SUBROUTINE HUBKINN(I,J,K,NBASISMAX,BHUB,TTILT,SUM,TREAL)
C returns the non-interacting energy of state with 
C..quatum numbers (i,j,k) for Hubbard model of 
C..lengths (LX,LY,LZ).  NON-PBC
C..Returned in sum
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NBASISMAX(4,5),AX,AY
      LOGICAL TTILT,TREAL
      INTEGER II,JJ,KK
      INCLUDE 'cons.inc'
      LX=NBASISMAX(1,2)-NBASISMAX(1,1)+1
      LY=NBASISMAX(2,2)-NBASISMAX(2,1)+1
      LZ=NBASISMAX(3,2)-NBASISMAX(3,1)+1
      AX=NBASISMAX(1,4)
      AY=NBASISMAX(2,4)
      IF(TTILT) THEN
CCCC.. NBASISMAX goes from -NMAXX+1 to MAXX so LX=2MAXX
         LX=NBASISMAX(1,5)
         LX=(LX*(AX*AX+AY*AY))
         SUM=COS(PI*(I*AX-J*AY)/(LX+1))+COS(PI*(I*AY+J*AX)/(LX+1))
      ELSE
         II=I-NBASISMAX(1,1)+1
         JJ=J-NBASISMAX(2,1)+1
         KK=K-NBASISMAX(3,1)+1
         IF(LX.GT.1) SUM=COS(PI*II/(LX+1))
         IF(LY.GT.1) SUM=SUM+COS(PI*JJ/(LY+1))
         IF(LZ.GT.1) SUM=SUM+COS(PI*KK/(LZ+1))
C         IF(LY.GT.1) SUM=SUM+COS(PI*J/(LY+1))
C         IF(LZ.GT.1) SUM=SUM+COS(PI*K/(LZ+1))
      ENDIF
      SUM=SUM*2.D0*BHUB
      IF(TREAL) SUM=0.D0
      RETURN
      END

      SUBROUTINE CALCTMATHUB(TMAT,NBASIS,NBASISMAX,BHUB,TTILT,G,TREAL,
     &         TPBC)
         IMPLICIT NONE
         INTEGER NBASIS,G(4,NBASIS),NBASISMAX(4,2)
         REAL*8 TMAT(NBASIS,NBASIS),BHUB
         INTEGER I,J
         INTEGER DX,DY,DZ,LX,LY,LZ
         LOGICAL TTILT,TREAL,TPBC
         REAL*8 SUM,S1
         OPEN(10,FILE='TMAT',STATUS='UNKNOWN')
         SUM=0.D0
         LX=NBASISMAX(1,2)-NBASISMAX(1,1)
         LY=NBASISMAX(2,2)-NBASISMAX(2,1)
         LZ=NBASISMAX(3,2)-NBASISMAX(3,1)
         IF(LY.EQ.0) LY=-1
         IF(LZ.EQ.0) LZ=-1
         IF(TREAL) THEN
            DO I=1,NBASIS
               DO J=1,NBASIS
                  DX=ABS(G(1,I)-G(1,J))
                  DY=ABS(G(2,I)-G(2,J))
                  DZ=ABS(G(3,I)-G(3,J))
                  SUM=BHUB
                  IF(TPBC) THEN
                     IF(DX.EQ.LX) THEN
                        DX=1
                        IF(LX.EQ.1) SUM=SUM+BHUB
                     ENDIF
                     IF(DY.EQ.LY) THEN
                        DY=1
                        IF(LY.EQ.1) SUM=SUM+BHUB
                     ENDIF
                     IF(DZ.EQ.LZ) THEN
                        DZ=1
                        IF(LZ.EQ.1) SUM=SUM+BHUB
                     ENDIF
                  ENDIF
                  IF(DX+DY+DZ.EQ.1.AND.G(4,I).EQ.G(4,J)) THEN
                     TMAT(I,J)=SUM
                     WRITE(10,*) I,J,SUM
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            IF(TPBC) THEN
               DO I=1,NBASIS
                  CALL HUBKIN(G(1,I),G(2,I),G(3,I),NBASISMAX,BHUB,
     &               TTILT,SUM)
                  TMAT(I,I)=SUM
                  WRITE(10,*) I,I,TMAT(I,I)
               ENDDO
            ELSE
               DO I=1,NBASIS
                  CALL HUBKINN(G(1,I),G(2,I),G(3,I),NBASISMAX,BHUB,
     &               TTILT,SUM)
                  TMAT(I,I)=SUM
                  WRITE(10,*) I,I,TMAT(I,I)
               ENDDO
            ENDIF
         ENDIF
         CLOSE(10)
         RETURN
      END
C.. NBASISMAX descriptor (1,3)
C
C HUBBARD:
C 0 Non-Tilted Lattice - pbc
C 1 Tilted Lattice - pbc
C 2 Non-Tilted lattice - no pbc
C 3 Tilted Lattice - no pbc
C.. four following are REAL
C 4 Non-Tilted Lattice - pbc
C 5 Tilted Lattice - pbc
C 6 Non-Tilted lattice - no pbc
C 7 Tilted Lattice - no pbc
C
      SUBROUTINE SETBASISLIM_HUB(NBASISMAX,NMAXX,NMAXY,NMAXZ,LEN,
     &            TPBC,TREAL)
         IMPLICIT NONE
         INTEGER NBASISMAX(4,3),NMAXX,NMAXY,NMAXZ,LEN
         LOGICAL TPBC,TREAL
         IF(TPBC) THEN
            NBASISMAX(1,3)=0
         ELSE
C.. Non-tilted, not pbc
            NBASISMAX(1,3)=2
         ENDIF
         IF(TREAL) NBASISMAX(1,3)=NBASISMAX(1,3)+4
C         IF(.NOT.TPBC.AND..NOT.TREAL) THEN
C.. non-pbc has Huckel MOs starting k from 1
C            NBASISMAX(1,1)=1 
C            NBASISMAX(2,1)=1 
C            NBASISMAX(3,1)=1
C            NBASISMAX(1,2)=NMAXX
C            NBASISMAX(2,2)=NMAXY
C            NBASISMAX(3,2)=NMAXZ
C         ELSE
         IF(MOD(NMAXX,2).EQ.0) THEN
            NBASISMAX(1,2)=NMAXX/2
            NBASISMAX(1,1)=-NMAXX/2+1
         ELSE
            NBASISMAX(1,2)=NMAXX/2
            NBASISMAX(1,1)=-NMAXX/2
         ENDIF
         IF(MOD(NMAXY,2).EQ.0) THEN
            NBASISMAX(2,2)=NMAXY/2
            NBASISMAX(2,1)=-NMAXY/2+1
         ELSE
            NBASISMAX(2,2)=NMAXY/2
            NBASISMAX(2,1)=-NMAXY/2
         ENDIF
         IF(MOD(NMAXZ,2).EQ.0) THEN
            NBASISMAX(3,2)=NMAXZ/2
            NBASISMAX(3,1)=-NMAXZ/2+1
         ELSE
            NBASISMAX(3,2)=NMAXZ/2
            NBASISMAX(3,1)=-NMAXZ/2
         ENDIF
C         ENDIF
         NBASISMAX(1,4)=0
         NBASISMAX(2,4)=1
         NBASISMAX(1,5)=NMAXX
         NBASISMAX(2,5)=NMAXY
         NBASISMAX(3,5)=NMAXZ
         LEN=NMAXX*NMAXY*NMAXZ*((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
      END

      SUBROUTINE SETBASISLIM_HUBTILT(NBASISMAX,NMAXX,NMAXY,NMAXZ,LEN,
     &         TPBC,ITILTX,ITILTY)
         INTEGER NBASISMAX(4,5),NMAXX,NMAXY,NMAXZ,LEN
         LOGICAL TPBC
         INTEGER ITILTX,ITILTY
         IF(TPBC) THEN
C.. Indicate tilted
            NBASISMAX(1,3)=1
         ELSE
C.. Not periodic boundaries
            NBASISMAX(1,3)=3
         ENDIF
         NBASISMAX(1,2)=(NMAXX*((ITILTX+ITILTY)/2.D0))
         NBASISMAX(1,1)=-NBASISMAX(1,2)
C+1-MOD(NMAXX,2)
         IF(NMAXY.NE.NMAXX)
     &      STOP 'CANNOT HANDLE NON-SQUARE TILTED HUBBARD'
         NBASISMAX(2,2)=(NMAXX*((ITILTY+ITILTX)/2.D0))
         NBASISMAX(2,1)=-NBASISMAX(2,2)
C+1-MOD(NMAXX,2)
         IF(NMAXZ.GT.1) STOP 'CANNOT HANDLE TILTED 3D HUBBARD'
         NBASISMAX(3,2)=0
         NBASISMAX(3,1)=0
         NBASISMAX(1,4)=ITILTX
         NBASISMAX(2,4)=ITILTY
         NBASISMAX(1,5)=NMAXX
         NBASISMAX(2,5)=NMAXY
         NBASISMAX(3,5)=NMAXZ
         LEN=NMAXX*NMAXY*(ITILTX*ITILTX+ITILTY*ITILTY)
     &         *((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
C         LEN=NMAXX*NMAXY*2
C     &         *((NBASISMAX(4,2)-NBASISMAX(4,1))/2+1)
      END

