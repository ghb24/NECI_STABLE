C.. This would seem to sort both ARR and BRR into ascending ARR order [AJWT20031001]
C.. Further investigation shows it is the quicksort from Numerical Recipes, taking two arrays to be sorted, ARR, BRR according to the order of ARR [AJWT20040407]
      SUBROUTINE NECI_SORT2(N,ARR,BRR)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER, PARAMETER :: M=7,NSTACK=50
      DIMENSION ARR(N),ISTACK(NSTACK)
      INTEGER BRR(N)
      JSTACK=0
      L=1
      IR=N
10    IF(IR-L.LT.M) THEN
        DO J=L+1,IR
          A=ARR(J)
          B=BRR(J)
          DO I=J-1,1,-1
            IF(ARR(I).LE.A) GOTO 20
            ARR(I+1)=ARR(I)
            BRR(I+1)=BRR(I)
          ENDDO
          I=0
20        ARR(I+1)=A
          BRR(I+1)=B
        ENDDO
        IF(JSTACK.EQ.0) RETURN
        IR=ISTACK(JSTACK)
        L=ISTACK(JSTACK-1)
        JSTACK=JSTACK-2
      ELSE
        K=(L+IR)/2
        TEMP=ARR(K)
        ARR(K)=ARR(L+1)
        ARR(L+1)=TEMP
        TEMP=BRR(K)
        BRR(K)=BRR(L+1)
        BRR(L+1)=TEMP
        IF(ARR(L+1).GT.ARR(IR)) THEN
          TEMP=ARR(L+1)
          ARR(L+1)=ARR(IR)
          ARR(IR)=TEMP
          TEMP=BRR(L+1)
          BRR(L+1)=BRR(IR)
          BRR(IR)=TEMP
        ENDIF
        IF(ARR(L).GT.ARR(IR)) THEN
          TEMP=ARR(L)
          ARR(L)=ARR(IR)
          ARR(IR)=TEMP
          TEMP=BRR(L)
          BRR(L)=BRR(IR)
          BRR(IR)=TEMP
        ENDIF 
        IF(ARR(L+1).GT.ARR(L)) THEN
          TEMP=ARR(L+1)
          ARR(L+1)=ARR(L)
          ARR(L)=TEMP
          TEMP=BRR(L+1)
          BRR(L+1)=BRR(L)
          BRR(L)=TEMP
        ENDIF
        I=L+1
        J=IR
        A=ARR(L)
        B=BRR(L)
30      CONTINUE
        I=I+1
        IF(ARR(I).LT.A) GOTO 30
40      CONTINUE
        J=J-1
        IF(ARR(J).GT.A) GOTO 40
        IF(J.LT.I) GOTO 50
        TEMP=ARR(I)
        ARR(I)=ARR(J)
        ARR(J)=TEMP
        TEMP=BRR(I)
        BRR(I)=BRR(J)
        BRR(J)=TEMP
        GOTO 30
50      ARR(L)=ARR(J)
        ARR(J)=A
        BRR(L)=BRR(J)
        BRR(J)=B
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
C.. Further investigation shows it is the quicksort from Numerical Recipes, taking two arrays to be sorted, ARR, BRR according to the order of ARR [AJWT20040407]
C.. This version has ARR as an int and BRR as a real 
      SUBROUTINE SORT2_(N,ARR,BRR)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, PARAMETER :: M=7,NSTACK=50
      INTEGER ARR(N)
      DIMENSION ISTACK(NSTACK)
      REAL*8 BRR(N)
      JSTACK=0
      L=1
      IR=N
10    IF(IR-L.LT.M) THEN
        DO J=L+1,IR
          A=ARR(J)
          B=BRR(J)
          DO I=J-1,1,-1
            IF(ARR(I).LE.A) GOTO 20
            ARR(I+1)=ARR(I)
            BRR(I+1)=BRR(I)
          ENDDO
          I=0
20        ARR(I+1)=A
          BRR(I+1)=B
        ENDDO
        IF(JSTACK.EQ.0) RETURN
        IR=ISTACK(JSTACK)
        L=ISTACK(JSTACK-1)
        JSTACK=JSTACK-2
      ELSE
        K=(L+IR)/2
        TEMP=ARR(K)
        ARR(K)=ARR(L+1)
        ARR(L+1)=TEMP
        TEMP=BRR(K)
        BRR(K)=BRR(L+1)
        BRR(L+1)=TEMP
        IF(ARR(L+1).GT.ARR(IR)) THEN
          TEMP=ARR(L+1)
          ARR(L+1)=ARR(IR)
          ARR(IR)=TEMP
          TEMP=BRR(L+1)
          BRR(L+1)=BRR(IR)
          BRR(IR)=TEMP
        ENDIF
        IF(ARR(L).GT.ARR(IR)) THEN
          TEMP=ARR(L)
          ARR(L)=ARR(IR)
          ARR(IR)=TEMP
          TEMP=BRR(L)
          BRR(L)=BRR(IR)
          BRR(IR)=TEMP
        ENDIF 
        IF(ARR(L+1).GT.ARR(L)) THEN
          TEMP=ARR(L+1)
          ARR(L+1)=ARR(L)
          ARR(L)=TEMP
          TEMP=BRR(L+1)
          BRR(L+1)=BRR(L)
          BRR(L)=TEMP
        ENDIF
        I=L+1
        J=IR
        A=ARR(L)
        B=BRR(L)
30      CONTINUE
        I=I+1
        IF(ARR(I).LT.A) GOTO 30
40      CONTINUE
        J=J-1
        IF(ARR(J).GT.A) GOTO 40
        IF(J.LT.I) GOTO 50
        TEMP=ARR(I)
        ARR(I)=ARR(J)
        ARR(J)=TEMP
        TEMP=BRR(I)
        BRR(I)=BRR(J)
        BRR(J)=TEMP
        GOTO 30
50      ARR(L)=ARR(J)
        ARR(J)=A
        BRR(L)=BRR(J)
        BRR(J)=B
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
C.. Further investigation shows it is the quicksort from Numerical Recipes, taking two arrays to be sorted, ARR, BRR according to the order of ARR [AJWT20040407]
C.. This version has ARR as an int and BRR as an array of length BD of reals
      SUBROUTINE SORTIRN(N,ARR,BRR,BD)
      IMPLICIT NONE
      INTEGER, PARAMETER :: M=7,NSTACK=50
      INTEGER ARR(N),A,L,IR,N,JSTACK,ISTACK(NSTACK),I,J,K
      INTEGER ITEMP,BD
      REAL*8 BRR(BD,N),B(BD),TEMP(BD)
      JSTACK=0
      L=1
      IR=N
10    IF(IR-L.LT.M) THEN
        DO J=L+1,IR
          A=ARR(J)
          B=BRR(:,J)
          DO I=J-1,1,-1
            IF(ARR(I).LE.A) GOTO 20
            ARR(I+1)=ARR(I)
            BRR(:,I+1)=BRR(:,I)
          ENDDO
          I=0
20        ARR(I+1)=A
          BRR(:,I+1)=B
        ENDDO
        IF(JSTACK.EQ.0) RETURN
        IR=ISTACK(JSTACK)
        L=ISTACK(JSTACK-1)
        JSTACK=JSTACK-2
      ELSE
        K=(L+IR)/2
        ITEMP=ARR(K)
        ARR(K)=ARR(L+1)
        ARR(L+1)=ITEMP
        TEMP=BRR(:,K)
        BRR(:,K)=BRR(:,L+1)
        BRR(:,L+1)=TEMP
        IF(ARR(L+1).GT.ARR(IR)) THEN
          ITEMP=ARR(L+1)
          ARR(L+1)=ARR(IR)
          ARR(IR)=ITEMP
          TEMP=BRR(:,L+1)
          BRR(:,L+1)=BRR(:,IR)
          BRR(:,IR)=TEMP
        ENDIF
        IF(ARR(L).GT.ARR(IR)) THEN
          ITEMP=ARR(L)
          ARR(L)=ARR(IR)
          ARR(IR)=ITEMP
          TEMP=BRR(:,L)
          BRR(:,L)=BRR(:,IR)
          BRR(:,IR)=TEMP
        ENDIF 
        IF(ARR(L+1).GT.ARR(L)) THEN
          ITEMP=ARR(L+1)
          ARR(L+1)=ARR(L)
          ARR(L)=ITEMP
          TEMP=BRR(:,L+1)
          BRR(:,L+1)=BRR(:,L)
          BRR(:,L)=TEMP
        ENDIF
        I=L+1
        J=IR
        A=ARR(L)
        B=BRR(:,L)
30      CONTINUE
        I=I+1
        IF(ARR(I).LT.A) GOTO 30
40      CONTINUE
        J=J-1
        IF(ARR(J).GT.A) GOTO 40
        IF(J.LT.I) GOTO 50
        ITEMP=ARR(I)
        ARR(I)=ARR(J)
        ARR(J)=ITEMP
        TEMP=BRR(:,I)
        BRR(:,I)=BRR(:,J)
        BRR(:,J)=TEMP
        GOTO 30
50      ARR(L)=ARR(J)
        ARR(J)=A
        BRR(:,L)=BRR(:,J)
        BRR(:,J)=B
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
C.. This would seem to sort both ABS(ARR) and BRR into ascending ARR order [AJWT20031001]
C.. Further investigation shows it is the quicksort from Numerical Recipes, taking two arrays to be sorted, ARR, BRR according to the order of ARR [AJWT20040407]
C.. This version has ARR as a real and BRR as an array of length BD of integers
      SUBROUTINE SORTRNI(N,ARR,BRR,BD)
      IMPLICIT NONE
      INTEGER, PARAMETER :: M=7,NSTACK=50
      INTEGER L,IR,N,JSTACK,ISTACK(NSTACK),I,J,K
      INTEGER BD
      INTEGER BRR(BD,N),B(BD),TEMP(BD)
      REAL*8 ARR(N),A,ITEMP
      JSTACK=0
      L=1
      IR=N
10    IF(IR-L.LT.M) THEN
        DO J=L+1,IR
          A=ARR(J)
          B=BRR(:,J)
          DO I=J-1,1,-1
            IF(ABS(ARR(I)).GE.ABS(A)) GOTO 20
            ARR(I+1)=ARR(I)
            BRR(:,I+1)=BRR(:,I)
          ENDDO
          I=0
20        ARR(I+1)=A
          BRR(:,I+1)=B
        ENDDO
        IF(JSTACK.EQ.0) RETURN
        IR=ISTACK(JSTACK)
        L=ISTACK(JSTACK-1)
        JSTACK=JSTACK-2
      ELSE
        K=(L+IR)/2
        ITEMP=ARR(K)
        ARR(K)=ARR(L+1)
        ARR(L+1)=ITEMP
        TEMP(:)=BRR(:,K)
        BRR(:,K)=BRR(:,L+1)
        BRR(:,L+1)=TEMP
        IF(ABS(ARR(L+1)).LT.ABS(ARR(IR))) THEN
          ITEMP=ARR(L+1)
          ARR(L+1)=ARR(IR)
          ARR(IR)=ITEMP
          TEMP=BRR(:,L+1)
          BRR(:,L+1)=BRR(:,IR)
          BRR(:,IR)=TEMP(:)
        ENDIF
        IF(ABS(ARR(L)).LT.ABS(ARR(IR))) THEN
          ITEMP=ARR(L)
          ARR(L)=ARR(IR)
          ARR(IR)=ITEMP
          TEMP=BRR(:,L)
          BRR(:,L)=BRR(:,IR)
          BRR(:,IR)=TEMP
        ENDIF 
        IF(ABS(ARR(L+1)).LT.ABS(ARR(L))) THEN
          ITEMP=ARR(L+1)
          ARR(L+1)=ARR(L)
          ARR(L)=ITEMP
          TEMP(:)=BRR(:,L+1)
          BRR(:,L+1)=BRR(:,L)
          BRR(:,L)=TEMP(:)
        ENDIF
        I=L+1
        J=IR
        A=ARR(L)
        B=BRR(:,L)
30      CONTINUE
        I=I+1
        IF(ABS(ARR(I)).GT.ABS(A)) GOTO 30
40      CONTINUE
        J=J-1
        IF(ABS(ARR(J)).LT.ABS(A)) GOTO 40
        IF(J.LT.I) GOTO 50
        ITEMP=ARR(I)
        ARR(I)=ARR(J)
        ARR(J)=ITEMP
        TEMP(:)=BRR(:,I)
        BRR(:,I)=BRR(:,J)
        BRR(:,J)=TEMP(:)
        GOTO 30
50      ARR(L)=ARR(J)
        ARR(J)=A
        BRR(:,L)=BRR(:,J)
        BRR(:,J)=B
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

      SUBROUTINE SORT3R(N,ARR,BRR,CRR)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER, PARAMETER :: M=7,NSTACK=50
      DIMENSION ARR(N),ISTACK(NSTACK),BRR(N),CRR(N)
      JSTACK=0
      L=1
      IR=N
10    IF(IR-L.LT.M) THEN
        DO J=L+1,IR
          A=ARR(J)
          B=BRR(J)
          C=CRR(J)
          DO I=J-1,1,-1
            IF(ARR(I).LE.A) GOTO 20
            ARR(I+1)=ARR(I)
            BRR(I+1)=BRR(I)
            CRR(I+1)=CRR(I)
          ENDDO
          I=0
20        ARR(I+1)=A
          BRR(I+1)=B
          CRR(I+1)=C
        ENDDO
        IF(JSTACK.EQ.0) RETURN
        IR=ISTACK(JSTACK)
        L=ISTACK(JSTACK-1)
        JSTACK=JSTACK-2
      ELSE
        K=(L+IR)/2
        TEMP=ARR(K)
        ARR(K)=ARR(L+1)
        ARR(L+1)=TEMP
        TEMP=BRR(K)
        BRR(K)=BRR(L+1)
        BRR(L+1)=TEMP
        TEMP=CRR(K)
        CRR(K)=CRR(L+1)
        CRR(L+1)=TEMP
        IF(ARR(L+1).GT.ARR(IR)) THEN
          TEMP=ARR(L+1)
          ARR(L+1)=ARR(IR)
          ARR(IR)=TEMP
          TEMP=BRR(L+1)
          BRR(L+1)=BRR(IR)
          BRR(IR)=TEMP
          TEMP=CRR(L+1)
          CRR(L+1)=CRR(IR)
          CRR(IR)=TEMP
        ENDIF
        IF(ARR(L).GT.ARR(IR)) THEN
          TEMP=ARR(L)
          ARR(L)=ARR(IR)
          ARR(IR)=TEMP
          TEMP=BRR(L)
          BRR(L)=BRR(IR)
          BRR(IR)=TEMP
          TEMP=CRR(L)
          CRR(L)=CRR(IR)
          CRR(IR)=TEMP
        ENDIF 
        IF(ARR(L+1).GT.ARR(L)) THEN
          TEMP=ARR(L+1)
          ARR(L+1)=ARR(L)
          ARR(L)=TEMP
          TEMP=BRR(L+1)
          BRR(L+1)=BRR(L)
          BRR(L)=TEMP
          TEMP=CRR(L+1)
          CRR(L+1)=CRR(L)
          CRR(L)=TEMP
        ENDIF
        I=L+1
        J=IR
        A=ARR(L)
        B=BRR(L)
        C=CRR(L)
30      CONTINUE
        I=I+1
        IF(ARR(I).LT.A) GOTO 30
40      CONTINUE
        J=J-1
        IF(ARR(J).GT.A) GOTO 40
        IF(J.LT.I) GOTO 50
        TEMP=ARR(I)
        ARR(I)=ARR(J)
        ARR(J)=TEMP
        TEMP=BRR(I)
        BRR(I)=BRR(J)
        BRR(J)=TEMP
        TEMP=CRR(I)
        CRR(I)=CRR(J)
        CRR(J)=TEMP
        GOTO 30
50      ARR(L)=ARR(J)
        ARR(J)=A
        BRR(L)=BRR(J)
        BRR(J)=B
        CRR(L)=CRR(J)
        CRR(J)=C
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


      SUBROUTINE SORT3R1I(N,ARR,BRR,CRR,DRR)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER, PARAMETER :: M=7,NSTACK=50
      DIMENSION ARR(N),ISTACK(NSTACK),BRR(N),CRR(N)
      INTEGER DRR(N),D,ITMP
      JSTACK=0
      L=1
      IR=N
10    IF(IR-L.LT.M) THEN
        DO J=L+1,IR
          A=ARR(J)
          B=BRR(J)
          C=CRR(J)
          D=DRR(J)
          DO I=J-1,1,-1
            IF(ARR(I).LE.A) GOTO 20
            ARR(I+1)=ARR(I)
            BRR(I+1)=BRR(I)
            CRR(I+1)=CRR(I)
            DRR(I+1)=DRR(I) 
         ENDDO
          I=0
20        ARR(I+1)=A
          BRR(I+1)=B
          CRR(I+1)=C
          DRR(I+1)=D
        ENDDO
        IF(JSTACK.EQ.0) RETURN
        IR=ISTACK(JSTACK)
        L=ISTACK(JSTACK-1)
        JSTACK=JSTACK-2
      ELSE
        K=(L+IR)/2
        TEMP=ARR(K)
        ARR(K)=ARR(L+1)
        ARR(L+1)=TEMP
        TEMP=BRR(K)
        BRR(K)=BRR(L+1)
        BRR(L+1)=TEMP
        TEMP=CRR(K)
        CRR(K)=CRR(L+1)
        CRR(L+1)=TEMP
        ITMP=DRR(K)
        DRR(K)=DRR(L+1)
        DRR(L+1)=ITMP
        IF(ARR(L+1).GT.ARR(IR)) THEN
          TEMP=ARR(L+1)
          ARR(L+1)=ARR(IR)
          ARR(IR)=TEMP
          TEMP=BRR(L+1)
          BRR(L+1)=BRR(IR)
          BRR(IR)=TEMP
          TEMP=CRR(L+1)
          CRR(L+1)=CRR(IR)
          CRR(IR)=TEMP
          ITMP=DRR(L+1)
          DRR(L+1)=DRR(IR)
          DRR(IR)=ITMP
        ENDIF
        IF(ARR(L).GT.ARR(IR)) THEN
          TEMP=ARR(L)
          ARR(L)=ARR(IR)
          ARR(IR)=TEMP
          TEMP=BRR(L)
          BRR(L)=BRR(IR)
          BRR(IR)=TEMP
          TEMP=CRR(L)
          CRR(L)=CRR(IR)
          CRR(IR)=TEMP
          ITMP=DRR(L)
          DRR(L)=DRR(IR)
          DRR(IR)=ITMP
        ENDIF 
        IF(ARR(L+1).GT.ARR(L)) THEN
          TEMP=ARR(L+1)
          ARR(L+1)=ARR(L)
          ARR(L)=TEMP
          TEMP=BRR(L+1)
          BRR(L+1)=BRR(L)
          BRR(L)=TEMP
          TEMP=CRR(L+1)
          CRR(L+1)=CRR(L)
          CRR(L)=TEMP
          ITMP=DRR(L+1)
          DRR(L+1)=DRR(L)
          DRR(L)=ITMP
        ENDIF
        I=L+1
        J=IR
        A=ARR(L)
        B=BRR(L)
        C=CRR(L)
        D=DRR(L)
30      CONTINUE
        I=I+1
        IF(ARR(I).LT.A) GOTO 30
40      CONTINUE
        J=J-1
        IF(ARR(J).GT.A) GOTO 40
        IF(J.LT.I) GOTO 50
        TEMP=ARR(I)
        ARR(I)=ARR(J)
        ARR(J)=TEMP
        TEMP=BRR(I)
        BRR(I)=BRR(J)
        BRR(J)=TEMP
        TEMP=CRR(I)
        CRR(I)=CRR(J)
        CRR(J)=TEMP
        ITMP=DRR(I)
        DRR(I)=DRR(J)
        DRR(J)=ITMP
        GOTO 30
50      ARR(L)=ARR(J)
        ARR(J)=A
        BRR(L)=BRR(J)
        BRR(J)=B
        CRR(L)=CRR(J)
        CRR(J)=C
        DRR(L)=DRR(J)
        DRR(J)=D
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
C.. Further investigation shows it is the quicksort from Numerical Recipes, taking two arrays to be sorted, ARR, BRR according to the order of ARR [AJWT20040407]
C.. This version has ARR BRR and CRR as arrays of length BD of reals.  Sort is based on Arr(1,:)
      SUBROUTINE SORT3RN(N,ARR,BRR,CRR,BD)
      IMPLICIT NONE
      INTEGER, PARAMETER  :: M=7,NSTACK=50
      INTEGER L,IR,N,JSTACK,I,J,K,BD
      REAL*8 ARR(BD,N),A(BD),ATEMP(BD)
      REAL*8 BRR(BD,N),B(BD),BTEMP(BD)
      REAL*8 CRR(BD,N),C(BD),CTEMP(BD)
      integer ISTACK(NSTACK)
      JSTACK=0
      L=1
      IR=N
10    IF(IR-L.LT.M) THEN
        DO J=L+1,IR
          A=ARR(:,J)
          B=BRR(:,J)
          C=CRR(:,J)
          DO I=J-1,1,-1
            IF(ARR(1,I).LE.A(1)) GOTO 20
            ARR(:,I+1)=ARR(:,I)
            BRR(:,I+1)=BRR(:,I)
            CRR(:,I+1)=CRR(:,I)
          ENDDO
          I=0
20        ARR(:,I+1)=A
          BRR(:,I+1)=B
          CRR(:,I+1)=C
        ENDDO
        IF(JSTACK.EQ.0) RETURN
        IR=ISTACK(JSTACK)
        L=ISTACK(JSTACK-1)
        JSTACK=JSTACK-2
      ELSE
        K=(L+IR)/2
        ATEMP=ARR(:,K)
        ARR(:,K)=ARR(:,L+1)
        ARR(:,L+1)=ATEMP
        BTEMP=BRR(:,K)
        BRR(:,K)=BRR(:,L+1)
        BRR(:,L+1)=BTEMP
        CTEMP=CRR(:,K)
        CRR(:,K)=CRR(:,L+1)
        CRR(:,L+1)=CTEMP
        IF(ARR(1,L+1).GT.ARR(1,IR)) THEN
          ATEMP=ARR(:,L+1)
          ARR(:,L+1)=ARR(:,IR)
          ARR(:,IR)=ATEMP
          BTEMP=BRR(:,L+1)
          BRR(:,L+1)=BRR(:,IR)
          BRR(:,IR)=BTEMP
          CTEMP=CRR(:,L+1)
          CRR(:,L+1)=CRR(:,IR)
          CRR(:,IR)=CTEMP
        ENDIF
        IF(ARR(1,L).GT.ARR(1,IR)) THEN
          ATEMP=ARR(:,L)
          ARR(:,L)=ARR(:,IR)
          ARR(:,IR)=ATEMP
          BTEMP=BRR(:,L)
          BRR(:,L)=BRR(:,IR)
          BRR(:,IR)=BTEMP
          CTEMP=CRR(:,L)
          CRR(:,L)=CRR(:,IR)
          CRR(:,IR)=CTEMP
        ENDIF 
        IF(ARR(1,L+1).GT.ARR(1,L)) THEN
          ATEMP=ARR(:,L+1)
          ARR(:,L+1)=ARR(:,L)
          ARR(:,L)=ATEMP
          BTEMP=BRR(:,L+1)
          BRR(:,L+1)=BRR(:,L)
          BRR(:,L)=BTEMP
          CTEMP=CRR(:,L+1)
          CRR(:,L+1)=CRR(:,L)
          CRR(:,L)=CTEMP
        ENDIF
        I=L+1
        J=IR
        A=ARR(:,L)
        B=BRR(:,L)
        C=CRR(:,L)
30      CONTINUE
        I=I+1
        IF(ARR(1,I).LT.A(1)) GOTO 30
40      CONTINUE
        J=J-1
        IF(ARR(1,J).GT.A(1)) GOTO 40
        IF(J.LT.I) GOTO 50
        ATEMP=ARR(:,I)
        ARR(:,I)=ARR(:,J)
        ARR(:,J)=ATEMP
        BTEMP=BRR(:,I)
        BRR(:,I)=BRR(:,J)
        BRR(:,J)=BTEMP
        CTEMP=CRR(:,I)
        CRR(:,I)=CRR(:,J)
        CRR(:,J)=CTEMP
        GOTO 30
50      ARR(:,L)=ARR(:,J)
        ARR(:,J)=A
        BRR(:,L)=BRR(:,J)
        BRR(:,J)=B
        CRR(:,L)=CRR(:,J)
        CRR(:,J)=C
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
C.. Further investigation shows it is the quicksort from Numerical Recipes, taking two arrays to be sorted, ARR, BRR according to the order of ARR [AJWT20040407]
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
C.. Further investigation shows it is the quicksort from Numerical Recipes, taking two arrays to be sorted, ARR, BRR according to the order of ARR [AJWT20040407]
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
 
