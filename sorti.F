C.. This would appear to be a standard quicksort from Numerical Recipes [AJWT 20040407]
C.. Based on SORT, SORTI sorts integers rather than reals
      SUBROUTINE SORTI(N,RA)
      IMPLICIT REAL*8 (A-H,O-Z) 
      INTEGER N,I,L,IR,J
      INTEGER RA(N)
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END

C.. Based on SORTI, SORT2I sorts two arrays of integers (ordering according to the first)
      SUBROUTINE SORT2I(N,RA,RB)
      IMPLICIT REAL*8 (A-H,O-Z) 
      INTEGER N,I,L,IR,J
      INTEGER RA(N)
      INTEGER RB(N)
      INTEGER RRA,RRB
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10
      END

C.. Based on SORTI, SORTIW sorts two arrays. It sorts lists of type(Walker). These consist of an array of integers, and a corresponding array of logicals. The array is returned ordered according to the integers - ghb24 - 31/5/08.
      SUBROUTINE SORTIW(N,RA)
      IMPLICIT NONE
      TYPE Walker
          INTEGER :: Det
          LOGICAL :: WSign
      END TYPE
      TYPE(Walker) RA(N)
      INTEGER N,I,L,IR,J
      INTEGER RRA
      LOGICAL RRB
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)%Det
          RRB=RA(L)%WSign
        ELSE
          RRA=RA(IR)%Det
          RA(IR)%Det=RA(1)%Det
          RRB=RA(IR)%WSign
          RA(IR)%WSign=RA(1)%WSign
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)%Det=RRA
            RA(1)%WSign=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF((RA(J)%Det).LT.(RA(J+1)%Det)) J=J+1
          ENDIF
          IF(RRA.LT.(RA(J)%Det))THEN
            RA(I)%Det=RA(J)%Det
            RA(I)%WSign=RA(J)%WSign
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)%Det=RRA
        RA(I)%WSign=RRB
      GO TO 10
      END


C.. Based on SORTI, SORTIIARR sorts an arrays of integers and a 2D array of integers (ordering according to the first)
      SUBROUTINE SORTIIARR(N,RA,RB,NRB)
      IMPLICIT REAL*8 (A-H,O-Z) 
      INTEGER N,I,L,IR,J
      INTEGER RA(N)
      INTEGER RB(NRB,N)
      INTEGER RRA,RRB(NRB),NRB
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          CALL ICOPY(NRB,RB(1,L),1,RRB,1)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          CALL ICOPY(NRB,RB(1,IR),1,RRB,1)
          CALL ICOPY(NRB,RB(1,1),1,RB(1,IR),1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            CALL ICOPY(NRB,RRB,1,RB(1,1),1)
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            CALL ICOPY(NRB,RB(1,J),1,RB(1,I),1)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        CALL ICOPY(NRB,RRB,1,RB(1,I),1)
      GO TO 10
      END
C.. Based on SORTI, SORTSYMIARR sorts aan arrays of types whose first elements are symmetries and a 2D array of integers (ordering according to the symmetries)
C.. RA is the array of types whose first element is a symmetry.
C.. NRA is the length (in numbers of integers) of each element of NRA
C.. RB is the 2D array of integers
C.. NRB is the length of each row in RB
      SUBROUTINE SORTSYMIARR(N,RA,NRA,RB,NRB)
      IMPLICIT REAL*8 (A-H,O-Z) 
      INTEGER N,I,L,IR,J
      INTEGER RA(NRA,N)
      INTEGER RB(NRB,N)
      INTEGER RRA(NRA),RRB(NRB),NRB

      LOGICAL SYMLT
 
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          CALL ICOPY(NRA,RA(1,L),1,RRA,1)
          CALL ICOPY(NRB,RB(1,L),1,RRB,1)
        ELSE
          CALL ICOPY(NRA,RA(1,IR),1,RRA,1)
          CALL ICOPY(NRA,RA(1,1),1,RA(1,IR),1)
          CALL ICOPY(NRB,RB(1,IR),1,RRB,1)
          CALL ICOPY(NRB,RB(1,1),1,RB(1,IR),1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            CALL ICOPY(NRA,RRA,1,RA(1,1),1)
            CALL ICOPY(NRB,RRB,1,RB(1,1),1)
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(SYMLT(RA(1,J),RA(1,J+1))) J=J+1
          ENDIF
          IF(SYMLT(RRA,RA(1,J)))THEN
            CALL ICOPY(NRA,RA(1,J),1,RA(1,I),1)
            CALL ICOPY(NRB,RB(1,J),1,RB(1,I),1)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        CALL ICOPY(NRA,RRA,1,RA(1,I),1)
        CALL ICOPY(NRB,RRB,1,RB(1,I),1)
      GO TO 10
      END
C.. Based on SORTSYMIARR, SORTSYMARR sorts an arrays of types whose first elements are symmetries ordering according to the symmetries
C.. RA is the array of types whose first element is a symmetry.
C.. NRA is the length (in numbers of integers) of each element of NRA
      SUBROUTINE SORTSYMARR(N,RA,NRA)
         IMPLICIT NONE
         INTEGER N,NRA,RA(NRA)
         CALL SORTSYMIARR(N,RA,NRA,0,0)
      END

C.. Based on SORTI, SORTDETS sorts arrays of integers, and takes the corresponding element from RB with it
C.. RA is the array of integers of length N to sort
C.. NRA is the length (in numbers of integers) of each element of NRA
C.. RB is the 2D array of integers
C.. NRB is the length of each row in RB
      SUBROUTINE SORTDETS(N,RA,NRA,RB,NRB)
      IMPLICIT REAL*8 (A-H,O-Z) 
      INTEGER N,I,L,IR,J
      INTEGER RA(NRA,N)
      INTEGER RB(NRB,N)
      INTEGER RRA(NRA),RRB(NRB),NRB

      INTEGER DETLT
 
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          CALL ICOPY(NRA,RA(1,L),1,RRA,1)
          CALL ICOPY(NRB,RB(1,L),1,RRB,1)
        ELSE
          CALL ICOPY(NRA,RA(1,IR),1,RRA,1)
          CALL ICOPY(NRA,RA(1,1),1,RA(1,IR),1)
          CALL ICOPY(NRB,RB(1,IR),1,RRB,1)
          CALL ICOPY(NRB,RB(1,1),1,RB(1,IR),1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            CALL ICOPY(NRA,RRA,1,RA(1,1),1)
            CALL ICOPY(NRB,RRB,1,RB(1,1),1)
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF((DETLT(RA(1,J),RA(1,J+1),NRA)).eq.-1) J=J+1
          ENDIF
          IF((DETLT(RRA,RA(1,J),NRA)).eq.-1)THEN
            CALL ICOPY(NRA,RA(1,J),1,RA(1,I),1)
            CALL ICOPY(NRB,RB(1,J),1,RB(1,I),1)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        CALL ICOPY(NRA,RRA,1,RA(1,I),1)
        CALL ICOPY(NRB,RRB,1,RB(1,I),1)
      GO TO 10
      END
