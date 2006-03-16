C ==--------------------------------------------------------------------==
      SUBROUTINE SCATTER_SR(NSO,I1,I2,I3,I4,V)
      IMPLICIT REAL*8 (A-H,O-Z)
C..This routine takes the four integer labels of the determinant 
C..and turns them into a vector containing solely 1's and 0's
      INTEGER LB(4),V(NSO)
C..
      CALL IAZZERO(LB,4)
C..
      LB(1)=I1
      LB(2)=I2
      LB(3)=I3
      LB(4)=I4
C..        
      CALL IAZZERO(V,NSO)
      DO I=1,4
        V(LB(I))=1
      ENDDO
C..
      RETURN
      END 
C ==--------------------------------------------------------------------==
      SUBROUTINE FLIP(NSO,V,VOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
C..This routine flips ALPHA spin (10) to BETA spin (01) for each electron
      INTEGER V(NSO),VOUT(NSO)
      CALL IAZZERO(VOUT,NSO)
C..
      DO I=1,NSO/2
        II=2*I
C..
 	VOUT(II)=V(II-1)
        VOUT(II-1)=V(II)
      ENDDO
C..
      RETURN
      END
C ==--------------------------------------------------------------------==
      SUBROUTINE GATHER(NSO,V,I1,I2,I3,I4)
      IMPLICIT REAL*8 (A-H,O-Z)
C..This routine *gathers* a scattered set of integers
      INTEGER LB(4),V(NSO)
      II=0
      CALL IAZZERO(LB,4)
C..
      DO I=1,NSO
        IF(V(I).EQ.1) THEN
          II=II+1
          LB(II)=I
        ENDIF
      ENDDO 
C..Now ensure they're in ascending order
100   CONTINUE
      DO I=1,3
        IF(LB(I).GT.LB(I+1)) THEN
          IT=LB(I)
          LB(I)=LB(I+1)
          LB(I+1)=IT
          GOTO 100
        ENDIF
      ENDDO
C..
      I1=LB(1)
      I2=LB(2)
      I3=LB(3)
      I4=LB(4)
C..
      END 
C ==--------------------------------------------------------------------==
      SUBROUTINE CYCLIC_PERMUTE(NSO,V,V1)
      IMPLICIT REAL*8 (A-H,O-Z)
C..This routine cyclically permutes the electrons *throughout the orbitals
C..they occupy ONLY*
      INTEGER V(NSO),V1(NSO),NPOS(4)
C..Zero V1
      CALL IAZZERO(V1,NSO)
C..Need to find out which orbitals the electrons are in
      II=0
      DO I=1,NSO
        IF(V(I).EQ.1) THEN
          II=II+1
          IF(MOD(I,2).EQ.1) IP1=(I+1)/2
          IF(MOD(I,2).EQ.0) IP1=I/2
          NPOS(II)=IP1
          IF(II.GE.4) GOTO 100
        ENDIF
      ENDDO  
 100  CONTINUE
C..
      V1(2*NPOS(1)-1)=V(2*NPOS(4)-1)
      V1(2*NPOS(1))=V(2*NPOS(4))
      DO I=1,3
        V1(2*NPOS(I+1)-1)=V(2*NPOS(I)-1)
        V1(2*NPOS(I+1))=V(2*NPOS(I))
      ENDDO
C      IS1=V(NSO-1)
C      IS2=V(NSO)
C..
C      DO I=1,NSO-2
C        V1(I+2)=V(I)
C      ENDDO
C      V1(1)=IS1
C      V1(2)=IS2
C..
      RETURN
      END
C ==--------------------------------------------------------------------==
      SUBROUTINE FLIP_MIDDLE(NSO,V,VOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
C..This routine finds the *middle* two electrons and flips there spins
C..ALPHA (10) -> BETA (01)
      INTEGER V(NSO),VOUT(NSO),IP1,IP2
C..
      CALL IAZZERO(VOUT,NSO)
C..
      DO I=1,NSO
        VOUT(I)=V(I)
      ENDDO
C..Find orbital which contains electron 2
      II=0
      DO I=1,NSO
        IF(VOUT(I).EQ.1) II=II+1
        IF(II.EQ.2) THEN
          IF(MOD(I,2).EQ.1) IP1=(I+1)/2
          IF(MOD(I,2).EQ.0) IP1=I/2
          GOTO 100
        ENDIF
      ENDDO  
 100  CONTINUE
C..Find orbital which contains electron 3
      II=0
      DO I=1,NSO
        IF(VOUT(I).EQ.1) II=II+1
        IF(II.EQ.3) THEN
          IF(MOD(I,2).EQ.1) IP2=(I+1)/2
          IF(MOD(I,2).EQ.0) IP2=I/2
          GOTO 200
        ENDIF
      ENDDO
 200  CONTINUE
C..Flip spins of second orbital
      ITMP=VOUT(2*IP1-1)
      VOUT(2*IP1-1)=VOUT(2*IP1)
      VOUT(2*IP1)=ITMP
C..Flip spins of third orbital
      ITMP=VOUT(2*IP2-1)
      VOUT(2*IP2-1)=VOUT(2*IP2)
      VOUT(2*IP2)=ITMP
C..
C      DO I=2,3
C        II=2*I
C..
C        VOUT(II)=V(II-1)
C        VOUT(II-1)=V(II)
C      ENDDO
C..
C      VOUT(1)=V(1)
C      VOUT(2)=V(2)
C      VOUT(7)=V(7)
C      VOUT(8)=V(8)
C      VOUT(9)=V(9)
C      VOUT(10)=V(10)
C..
      RETURN
      END
C ==--------------------------------------------------------------------==
      SUBROUTINE FLIP_END(NSO,V,VOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
C..This routine flips the *end* two electrons
C..ALPHA (10) -> BETA (01)
      INTEGER V(NSO),VOUT(NSO)
C..
      CALL IAZZERO(VOUT,NSO)
C..
      DO I=1,NSO
        VOUT(I)=V(I)
      ENDDO
C..
      II=0
      DO I=1,NSO
        IF(VOUT(I).EQ.1) II=II+1
        IF(II.EQ.3) THEN
          IF(MOD(I,2).EQ.1) IP1=(I+1)/2
          IF(MOD(I,2).EQ.0) IP1=I/2
          GOTO 100
        ENDIF
      ENDDO
 100  CONTINUE
      II=0
      DO I=1,NSO
        IF(VOUT(I).EQ.1) II=II+1
        IF(II.EQ.4) THEN
          IF(MOD(I,2).EQ.1) IP2=(I+1)/2
          IF(MOD(I,2).EQ.0) IP2=I/2
          GOTO 200
        ENDIF
      ENDDO
 200  CONTINUE
C..
      ITMP=VOUT(2*IP1-1)
      VOUT(2*IP1-1)=VOUT(2*IP1)
      VOUT(2*IP1)=ITMP
C..
      ITMP=VOUT(2*IP2-1)
      VOUT(2*IP2-1)=VOUT(2*IP2)
      VOUT(2*IP2)=ITMP
C..
C      DO I=3,4
C        II=2*I
C..
C 	VOUT(II)=V(II-1)
C        VOUT(II-1)=V(II)
C      ENDDO
C..
C      VOUT(1)=V(1)
C      VOUT(2)=V(2)
C      VOUT(3)=V(3)
C      VOUT(4)=V(4)
C..
      RETURN
      END
C ==--------------------------------------------------------------------==
      SUBROUTINE NODES(I1,I2,I3,I4,INDE)
      IMPLICIT REAL*8 (A-H,O-Z)
C..This routine counts the number of nodes in the determinant i.e.
C..ALPHA BETA ALPHA BETA => NNODE=3 
      L1=MOD(I1,2)
      L2=MOD(I2,2)
      L3=MOD(I3,2)
      L4=MOD(I4,2)
C..
      INDE=0
      IF(L1.NE.L2) INDE=INDE+1
      IF(L2.NE.L3) INDE=INDE+1
      IF(L3.NE.L4) INDE=INDE+1
C..
      RETURN
      END
C ==--------------------------------------------------------------------==
      SUBROUTINE STORE(NEL,NSO,BRR,I1,I2,I3,I4,LAB,IKILL)
      IMPLICIT REAL*8 (A-H,O-Z)
C..This routine stores the set of integers as an accepted determinant
      DIMENSION LAB(NEL),IKILL(NSO,NSO,NSO,NSO)
      INTEGER BRR(NSO)
C..
      LAB(1)=BRR(I1)
      LAB(2)=BRR(I2)
      LAB(3)=BRR(I3)
      LAB(4)=BRR(I4)
C..
      IKILL(LAB(4),LAB(3),LAB(2),LAB(1))=1
C..
      RETURN
      END
C ==--------------------------------------------------------------------==
