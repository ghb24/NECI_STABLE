
C.. This would seem to sort both ARR and BRR into ascending ARR order [AJWT20031001]
C.. Further investigation shows it is the heapsort from Numerical Recipes, taking two arrays to be sorted, ARR, BRR according to the order of ARR [JSS 03/03/2010]
C.. We only sort every NSKIP item in the array.
      SUBROUTINE SORT2SKIP(N,ARR,BRR,NSKIP)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER, PARAMETER :: M=7,NSTACK=50
      DIMENSION ARR(NSKIP,N/NSKIP),ISTACK(NSTACK)
      INTEGER BRR(NSKIP,N/NSKIP),NSKIP
      JSTACK=0
      L=1
      IR=N/NSKIP
10    IF(IR-L.LT.M) THEN
        DO J=L+1,IR
          A=ARR(1,J)
          B=BRR(1,J)
          DO I=J-1,1,-1
            IF(ARR(1,I).LE.A) GOTO 20
            ARR(1,I+1)=ARR(1,I)
            BRR(1,I+1)=BRR(1,I)
          ENDDO
          I=0
20        ARR(1,I+1)=A
          BRR(1,I+1)=B
        ENDDO
        IF(JSTACK.EQ.0) RETURN
        IR=ISTACK(JSTACK)
        L=ISTACK(JSTACK-1)
        JSTACK=JSTACK-2
      ELSE
        K=(L+IR)/2
        TEMP=ARR(1,K)
        ARR(1,K)=ARR(1,L+1)
        ARR(1,L+1)=TEMP
        TEMP=BRR(1,K)
        BRR(1,K)=BRR(1,L+1)
        BRR(1,L+1)=TEMP
        IF(ARR(1,L+1).GT.ARR(1,IR)) THEN
          TEMP=ARR(1,L+1)
          ARR(1,L+1)=ARR(1,IR)
          ARR(1,IR)=TEMP
          TEMP=BRR(1,L+1)
          BRR(1,L+1)=BRR(1,IR)
          BRR(1,IR)=TEMP
        ENDIF
        IF(ARR(1,L).GT.ARR(1,IR)) THEN
          TEMP=ARR(1,L)
          ARR(1,L)=ARR(1,IR)
          ARR(1,IR)=TEMP
          TEMP=BRR(1,L)
          BRR(1,L)=BRR(1,IR)
          BRR(1,IR)=TEMP
        ENDIF 
        IF(ARR(1,L+1).GT.ARR(1,L)) THEN
          TEMP=ARR(1,L+1)
          ARR(1,L+1)=ARR(1,L)
          ARR(1,L)=TEMP
          TEMP=BRR(1,L+1)
          BRR(1,L+1)=BRR(1,L)
          BRR(1,L)=TEMP
        ENDIF
        I=L+1
        J=IR
        A=ARR(1,L)
        B=BRR(1,L)
30      CONTINUE
        I=I+1
        IF(ARR(1,I).LT.A) GOTO 30
40      CONTINUE
        J=J-1
        IF(ARR(1,J).GT.A) GOTO 40
        IF(J.LT.I) GOTO 50
        TEMP=ARR(1,I)
        ARR(1,I)=ARR(1,J)
        ARR(1,J)=TEMP
        TEMP=BRR(1,I)
        BRR(1,I)=BRR(1,J)
        BRR(1,J)=TEMP
        GOTO 30
50      ARR(1,L)=ARR(1,J)
        ARR(1,J)=A
        BRR(1,L)=BRR(1,J)
        BRR(1,J)=B
        JSTACK=JSTACK+2
        IF(JSTACK.GT.NSTACK) STOP ' NSTACK TOO SMALL INSORT2 '
        IF(IR-I+1.GE.J-L) THEN
          ISTACK(JSTACK)=IR
          ISTACK(JSTACK-1)=I
          IR=J-1
        ELSE
          ISTACK(JSTACK)=J-1
          ISTACK(JSTACK-1)=L
          L=I
        ENDIF
      ENDIF
      GOTO 10
      RETURN
      END 

C.. This would seem to sort both ARR and BRR into ascending ARR order [AJWT20031001]
C.. Further investigation shows it is the heapsort from Numerical Recipes, taking two arrays to be sorted, ARR, BRR according to the order of ARR [JSS 03/03/2010]
C.. We only sort every NSKIP item in the array.
C.. This version has ARR as an int and BRR as a real 
      SUBROUTINE SORT2SKIP_(N,ARR,BRR,NSKIP)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER, PARAMETER :: M=7,NSTACK=50
      DIMENSION ISTACK(NSTACK)
      INTEGER ARR(NSKIP,N/NSKIP)
      REAL*8 BRR(NSKIP,N/NSKIP)
      INTEGER NSKIP
      JSTACK=0
      L=1
      IR=N/NSKIP
10    IF(IR-L.LT.M) THEN
        DO J=L+1,IR
          A=ARR(1,J)
          B=BRR(1,J)
          DO I=J-1,1,-1
            IF(ARR(1,I).LE.A) GOTO 20
            ARR(1,I+1)=ARR(1,I)
            BRR(1,I+1)=BRR(1,I)
          ENDDO
          I=0
20        ARR(1,I+1)=A
          BRR(1,I+1)=B
        ENDDO
        IF(JSTACK.EQ.0) RETURN
        IR=ISTACK(JSTACK)
        L=ISTACK(JSTACK-1)
        JSTACK=JSTACK-2
      ELSE
        K=(L+IR)/2
        TEMP=ARR(1,K)
        ARR(1,K)=ARR(1,L+1)
        ARR(1,L+1)=TEMP
        TEMP=BRR(1,K)
        BRR(1,K)=BRR(1,L+1)
        BRR(1,L+1)=TEMP
        IF(ARR(1,L+1).GT.ARR(1,IR)) THEN
          TEMP=ARR(1,L+1)
          ARR(1,L+1)=ARR(1,IR)
          ARR(1,IR)=TEMP
          TEMP=BRR(1,L+1)
          BRR(1,L+1)=BRR(1,IR)
          BRR(1,IR)=TEMP
        ENDIF
        IF(ARR(1,L).GT.ARR(1,IR)) THEN
          TEMP=ARR(1,L)
          ARR(1,L)=ARR(1,IR)
          ARR(1,IR)=TEMP
          TEMP=BRR(1,L)
          BRR(1,L)=BRR(1,IR)
          BRR(1,IR)=TEMP
        ENDIF 
        IF(ARR(1,L+1).GT.ARR(1,L)) THEN
          TEMP=ARR(1,L+1)
          ARR(1,L+1)=ARR(1,L)
          ARR(1,L)=TEMP
          TEMP=BRR(1,L+1)
          BRR(1,L+1)=BRR(1,L)
          BRR(1,L)=TEMP
        ENDIF
        I=L+1
        J=IR
        A=ARR(1,L)
        B=BRR(1,L)
30      CONTINUE
        I=I+1
        IF(ARR(1,I).LT.A) GOTO 30
40      CONTINUE
        J=J-1
        IF(ARR(1,J).GT.A) GOTO 40
        IF(J.LT.I) GOTO 50
        TEMP=ARR(1,I)
        ARR(1,I)=ARR(1,J)
        ARR(1,J)=TEMP
        TEMP=BRR(1,I)
        BRR(1,I)=BRR(1,J)
        BRR(1,J)=TEMP
        GOTO 30
50      ARR(1,L)=ARR(1,J)
        ARR(1,J)=A
        BRR(1,L)=BRR(1,J)
        BRR(1,J)=B
        JSTACK=JSTACK+2
        IF(JSTACK.GT.NSTACK) STOP ' NSTACK TOO SMALL INSORT2 '
        IF(IR-I+1.GE.J-L) THEN
          ISTACK(JSTACK)=IR
          ISTACK(JSTACK-1)=I
          IR=J-1
        ELSE
          ISTACK(JSTACK)=J-1
          ISTACK(JSTACK-1)=L
          L=I
        ENDIF
      ENDIF
      GOTO 10
      RETURN
      END 

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
 
! RA and RB have length N.
! Based on all the other sorts - (but in this case modified from SortEvecbyEval) - this sorts 2 REAL arrays based on descending order of RA.
! No absolutes are taken so -ve's appear before +ve values.
      SUBROUTINE Sort2Real(N,RA,RB)
      INTEGER N,I,L,IR,J
      REAL*8 RA(N)
      REAL*8 RB(N)
      REAL*8 RRA,RRB

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
            IF(RA(J).gt.RA(J+1)) J=J+1
          ENDIF
          IF(RRA.gt.RA(J)) THEN
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

      END SUBROUTINE Sort2Real
 

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
 
