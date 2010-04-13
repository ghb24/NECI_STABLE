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
