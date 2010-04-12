! Based on SortBitSign SortEvecbyEval, however in SortEvecbyEval, RA is an array of reals (evalues) to be sorted in descending order of absolute size, and
! RB is an array of elements going from 1:NVec (evectors) to be taken with the element of RA. 
! RA has length N.
      SUBROUTINE SortEvecbyEval(N,RA,NVec,RB)
      INTEGER N,NVec,I,L,IR,J
      REAL*8 RA(N)
      REAL*8 RB(1:NVec,N)
      REAL*8 RRA,RRB(1:NVec)

      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB(1:NVec)=RB(1:NVec,L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          RRB(1:NVec)=RB(1:NVec,IR)
          RB(1:NVec,IR)=RB(1:NVec,1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1:NVec,1)=RRB(1:NVec)
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ABS(RA(J)).gt.ABS(RA(J+1))) J=J+1
          ENDIF
          IF(ABS(RRA).gt.ABS(RA(J))) THEN
            RA(I)=RA(J)
            RB(1:NVec,I)=RB(1:NVec,J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(1:NVec,I)=RRB(1:NVec)

      GO TO 10

      END SUBROUTINE SortEvecbyEval
 

! Based on SortBitSign SortEvecbyEvalPlus1, however in SortEvecbyEvalPlus1, RA is an array of reals (evalues) to be sorted in descending order of absolute size, and
! RB is an array of elements going from 0:NVec (evectors) to be taken with the element of RA.  RC is a third array of integers also to be taken with the element of RA.
! RA has length N.
      SUBROUTINE SortEvecbyEvalPlus1(N,RA,NVec,RB,RC)
      INTEGER N,NVec,I,L,IR,J
      REAL*8 RA(N)
      REAL*8 RB(1:NVec,1:N)
      INTEGER RC(N)
      REAL*8 RRA,RRB(1:NVec)
      INTEGER RRC
 
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB(1:NVec)=RB(1:NVec,L)
          RRC=RC(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          RRB(1:NVec)=RB(1:NVec,IR)
          RB(1:NVec,IR)=RB(1:NVec,1)
          RRC=RC(IR)
          RC(IR)=RC(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1:NVec,1)=RRB(1:NVec)
            RC(1)=RRC
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ABS(RA(J)).gt.ABS(RA(J+1))) J=J+1
          ENDIF
          IF(ABS(RRA).gt.ABS(RA(J))) THEN
            RA(I)=RA(J)
            RB(1:NVec,I)=RB(1:NVec,J)
            RC(I)=RC(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(1:NVec,I)=RRB(1:NVec)
        RC(I)=RRC

      GO TO 10

      END SUBROUTINE SortEvecbyEvalPlus1
 
