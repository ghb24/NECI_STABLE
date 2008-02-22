C.. This would seem to sort both ARR and BRR into ascending ARR order [AJWT20031001]
C.. Further investigation shows it is the quicksort from Numerical Recipes, taking two arrays to be sorted, ARR, BRR according to the order of ARR [AJWT20040407]
      SUBROUTINE SORT2(N,ARR,BRR)
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
