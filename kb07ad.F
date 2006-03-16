C     ==================================================================
      SUBROUTINE KB07AD(COUNT,N,INDEX)
C     ==--------------------------------------------------------------==
C     == Sorting routine for the reciprocal space vectors (g)         ==
C     == KB07AD HANDLES DOUBLE PRECISION VARIABLES                    ==
C     == STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)            ==
C     == THE WORK-SPACE 'MARK' OF LENGTH 50 PERMITS UP TO 2**(50/2)   ==
C     == NUMBERS TO BE SORTED. THIS IS MORE THAN THE IBM VIRTUAL      ==
C     == MEMORY SPACE WILL HOLD .                                     ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N
      REAL*8  COUNT(N)
      INTEGER INDEX(N)
C     Variables
      INTEGER MARK(50),I,K,K1,LNGTH,IP,IFK,INTEST,INT,IY,J,IS1,IFKA,
     &        MLOOP,IF,IS,LA,M
      REAL*8  AV,X
C     ==--------------------------------------------------------------==
C     ==  SET INDEX ARRAY TO ORIGINAL ORDER .                         ==
C     ==--------------------------------------------------------------==
      DO I=1,N
        INDEX(I)=I
      ENDDO
C     ==--------------------------------------------------------------==
C     == CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED.              ==
C     ==--------------------------------------------------------------==
      IF(N.EQ.1)GOTO 200
      IF(N.GE.1)GOTO 30
      WRITE(*,20)
   20 FORMAT(/,20X,
     *'***KB07AD***NO NUMBERS TO BE SORTED ** RETURN TO CALLING PROGRAM'
     *)
      GOTO 200
C     ==--------------------------------------------------------------==
C     == 'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER  ==
C     == THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.         ==
C     ==--------------------------------------------------------------==
   30 M=12
C     ==--------------------------------------------------------------==
C     == SET UP INITIAL VALUES.                                       ==
C     ==--------------------------------------------------------------==
      LA=2
      IS=1
      IF=N
      DO 190 MLOOP=1,N
C     ==--------------------------------------------------------------==
C     ==  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE. ==
C     ==--------------------------------------------------------------==
        IFKA=IF-IS
        IF((IFKA+1).GT.M)GOTO 70
C     ==--------------------------------------------------------------==
C     == FINAL SORTING  ( A SIMPLE BUBBLE SORT )                      ==
C     ==--------------------------------------------------------------==
        IS1=IS+1
        DO 60 J=IS1,IF
          I=J
   40     IF(COUNT(I-1).LT.COUNT(I))GOTO 60
          IF(COUNT(I-1).GT.COUNT(I))GOTO 50
          IF(INDEX(I-1).LT.INDEX(I))GOTO 60
   50     AV=COUNT(I-1)
          COUNT(I-1)=COUNT(I)
          COUNT(I)=AV
          INT=INDEX(I-1)
          INDEX(I-1)=INDEX(I)
          INDEX(I)=INT
          I=I-1
          IF(I.GT.IS)GOTO  40
   60   CONTINUE
        LA=LA-2
        GOTO 170
C     ==--------------------------------------------------------------==
C     ==                *******  QUICKSORT  ********                  ==
C     == SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS  ==
C     == THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENTS ==
C     == HIGHEST ADDRESS.                                             ==
C     ==--------------------------------------------------------------==
   70   IY=(IS+IF)/2
        X=COUNT(IY)
        INTEST=INDEX(IY)
        COUNT(IY)=COUNT(IF)
        INDEX(IY)=INDEX(IF)
C     ==--------------------------------------------------------------==
C     == THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END ==
C     == OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE   ==
C     == OF X .                                                       ==
C     ==--------------------------------------------------------------==
        K=1
        IFK=IF
C     ==--------------------------------------------------------------==
C     == WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE ==
C     == INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AND INDICES AS   ==
C     == NECESSARY, UNTIL THEY MEET .                                 ==
C     ==--------------------------------------------------------------==
        DO 110 I=IS,IF
          IF(X.GT.COUNT(I))GOTO 110
          IF(X.LT.COUNT(I))GOTO 80
          IF(INTEST.GT.INDEX(I))GOTO 110
   80     IF(I.GE.IFK)GOTO 120
          COUNT(IFK)=COUNT(I)
          INDEX(IFK)=INDEX(I)
          K1=K
          DO 100 K=K1,IFKA
            IFK=IF-K
            IF(COUNT(IFK).GT.X)GOTO 100
            IF(COUNT(IFK).LT.X)GOTO 90
            IF(INTEST.LE.INDEX(IFK))GOTO 100
   90       IF(I.GE.IFK)GOTO 130
            COUNT(I)=COUNT(IFK)
            INDEX(I)=INDEX(IFK)
            GO TO 110
  100     CONTINUE
          GOTO 120
  110   CONTINUE
C     ==--------------------------------------------------------------==
C     == RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER  ==
C     == WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO ==
C     == 2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN OR EQUAL ==
C     == TO ANY ELEMENT IN THE SECOND PART, AND THEY MAY NOW BE SORTED==
C     == INDEPENDENTLY .                                              ==
C     ==--------------------------------------------------------------==
  120   COUNT(IFK)=X
        INDEX(IFK)=INTEST
        IP=IFK
        GOTO 140
  130   COUNT(I)=X
        INDEX(I)=INTEST
        IP=I
C     ==--------------------------------------------------------------==
C     ==  STORE THE LONGER SUBDIVISION IN WORKSPACE.                  ==
C     ==--------------------------------------------------------------==
  140   IF((IP-IS).GT.(IF-IP))GOTO 150
        MARK(LA)=IF
        MARK(LA-1)=IP+1
        IF=IP-1
        GOTO 160
  150   MARK(LA)=IP-1
        MARK(LA-1)=IS
        IS=IP+1
C     ==--------------------------------------------------------------==
C     == FIND THE LENGTH OF THE SHORTER SUBDIVISION.                  ==
C     ==--------------------------------------------------------------==
  160   LNGTH=IF-IS
        IF(LNGTH.LE.0)GOTO 180
C     ==--------------------------------------------------------------==
C     == IF IT CONTAINS MORE THAN ONE ELEMENT SUPPLY IT WITH WORKSPACE==
C     ==--------------------------------------------------------------==
        LA=LA+2
        GOTO 190
  170   IF(LA.LE.0)GOTO 200
C     ==--------------------------------------------------------------==
C     == OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT==
C     ==--------------------------------------------------------------==
  180   IF=MARK(LA)
        IS=MARK(LA-1)
  190 CONTINUE
C     ==--------------------------------------------------------------==
  200 RETURN
      END
C     ==================================================================
      SUBROUTINE KB07AI(COUNT,N,INDEX)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER N,COUNT(N),INDEX(N)
C     Variables
      INTEGER MARK(50),I,LNGTH,IP,K,K1,IFK,INTEST,IY,INT,J,IS1,IFKA,
     &        MLOOP,IF,M,LA,IS
      INTEGER AV,X
C     ==--------------------------------------------------------------==
C     ==  SET INDEX ARRAY TO ORIGINAL ORDER .                         ==
C     ==--------------------------------------------------------------==
      DO I=1,N
        INDEX(I)=I
      ENDDO
C     ==--------------------------------------------------------------==
C     == CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED.              ==
C     ==--------------------------------------------------------------==
      IF(N.EQ.1)GOTO 200
      IF(N.GE.1)GOTO 30
      WRITE(*,20)
   20 FORMAT(/,20X,
     *'***KB07AD***NO NUMBERS TO BE SORTED ** RETURN TO CALLING PROGRAM'
     *)
      GOTO 200
C     ==--------------------------------------------------------------==
C     == 'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER  ==
C     == THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.         ==
C     ==--------------------------------------------------------------==
   30 M=12
C     ==--------------------------------------------------------------==
C     == SET UP INITIAL VALUES.                                       ==
C     ==--------------------------------------------------------------==
      LA=2
      IS=1
      IF=N
      DO 190 MLOOP=1,N
C     ==--------------------------------------------------------------==
C     ==  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE. ==
C     ==--------------------------------------------------------------==
        IFKA=IF-IS
        IF((IFKA+1).GT.M)GOTO 70
C     ==--------------------------------------------------------------==
C     == FINAL SORTING  ( A SIMPLE BUBBLE SORT )                      ==
C     ==--------------------------------------------------------------==
        IS1=IS+1
        DO 60 J=IS1,IF
          I=J
   40     IF(COUNT(I-1).LT.COUNT(I))GOTO 60
          IF(COUNT(I-1).GT.COUNT(I))GOTO 50
          IF(INDEX(I-1).LT.INDEX(I))GOTO 60
   50     AV=COUNT(I-1)
          COUNT(I-1)=COUNT(I)
          COUNT(I)=AV
          INT=INDEX(I-1)
          INDEX(I-1)=INDEX(I)
          INDEX(I)=INT
          I=I-1
          IF(I.GT.IS)GOTO  40
   60   CONTINUE
        LA=LA-2
        GOTO 170
C     ==--------------------------------------------------------------==
C     ==                *******  QUICKSORT  ********                  ==
C     == SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS  ==
C     == THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENTS ==
C     == HIGHEST ADDRESS.                                             ==
C     ==--------------------------------------------------------------==
   70   IY=(IS+IF)/2
        X=COUNT(IY)
        INTEST=INDEX(IY)
        COUNT(IY)=COUNT(IF)
        INDEX(IY)=INDEX(IF)
C     ==--------------------------------------------------------------==
C     == THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END ==
C     == OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE   ==
C     == OF X .                                                       ==
C     ==--------------------------------------------------------------==
        K=1
        IFK=IF
C     ==--------------------------------------------------------------==
C     == WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE ==
C     == INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AND INDICES AS   ==
C     == NECESSARY, UNTIL THEY MEET .                                 ==
C     ==--------------------------------------------------------------==
        DO 110 I=IS,IF
          IF(X.GT.COUNT(I))GOTO 110
          IF(X.LT.COUNT(I))GOTO 80
          IF(INTEST.GT.INDEX(I))GOTO 110
   80     IF(I.GE.IFK)GOTO 120
          COUNT(IFK)=COUNT(I)
          INDEX(IFK)=INDEX(I)
          K1=K
          DO 100 K=K1,IFKA
            IFK=IF-K
            IF(COUNT(IFK).GT.X)GOTO 100
            IF(COUNT(IFK).LT.X)GOTO 90
            IF(INTEST.LE.INDEX(IFK))GOTO 100
   90       IF(I.GE.IFK)GOTO 130
            COUNT(I)=COUNT(IFK)
            INDEX(I)=INDEX(IFK)
            GO TO 110
  100     CONTINUE
          GOTO 120
  110   CONTINUE
C     ==--------------------------------------------------------------==
C     == RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER  ==
C     == WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO ==
C     == 2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN OR EQUAL ==
C     == TO ANY ELEMENT IN THE SECOND PART, AND THEY MAY NOW BE SORTED==
C     == INDEPENDENTLY .                                              ==
C     ==--------------------------------------------------------------==
  120   COUNT(IFK)=X
        INDEX(IFK)=INTEST
        IP=IFK
        GOTO 140
  130   COUNT(I)=X
        INDEX(I)=INTEST
        IP=I
C     ==--------------------------------------------------------------==
C     ==  STORE THE LONGER SUBDIVISION IN WORKSPACE.                  ==
C     ==--------------------------------------------------------------==
  140   IF((IP-IS).GT.(IF-IP))GOTO 150
        MARK(LA)=IF
        MARK(LA-1)=IP+1
        IF=IP-1
        GOTO 160
  150   MARK(LA)=IP-1
        MARK(LA-1)=IS
        IS=IP+1
C     ==--------------------------------------------------------------==
C     == FIND THE LENGTH OF THE SHORTER SUBDIVISION.                  ==
C     ==--------------------------------------------------------------==
  160   LNGTH=IF-IS
        IF(LNGTH.LE.0)GOTO 180
C     ==--------------------------------------------------------------==
C     == IF IT CONTAINS MORE THAN ONE ELEMENT SUPPLY IT WITH WORKSPACE==
C     ==--------------------------------------------------------------==
        LA=LA+2
        GOTO 190
  170   IF(LA.LE.0)GOTO 200
C     ==--------------------------------------------------------------==
C     == OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT==
C     ==--------------------------------------------------------------==
  180   IF=MARK(LA)
        IS=MARK(LA-1)
  190 CONTINUE
C     ==--------------------------------------------------------------==
  200 RETURN
      END
C     ==================================================================
#ifdef POINTER8
      SUBROUTINE KB07AI8(COUNT,N,INDEX)
C     ==--------------------------------------------------------------==
C     == VERSION FOR INTEGER*8 (used in PRINT_MEMORY)                 ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      INTEGER   N,INDEX(N)
      INTEGER*8 COUNT(N)
C     Variables
      INTEGER   MARK(50),I,LNGTH,IP,K,K1,IFK,INTEST,IY,INT,J,IS1,IFKA,
     &          MLOOP,IF,M,LA,IS
      INTEGER*8 AV,X
C     ==--------------------------------------------------------------==
C     ==  SET INDEX ARRAY TO ORIGINAL ORDER .                         ==
C     ==--------------------------------------------------------------==
      DO I=1,N
        INDEX(I)=I
      ENDDO
C     ==--------------------------------------------------------------==
C     == CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED.              ==
C     ==--------------------------------------------------------------==
      IF(N.EQ.1)GOTO 200
      IF(N.GE.1)GOTO 30
      WRITE(*,20)
   20 FORMAT(/,20X,
     *'***KB07AD***NO NUMBERS TO BE SORTED ** RETURN TO CALLING PROGRAM'
     *)
      GOTO 200
C     ==--------------------------------------------------------------==
C     == 'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER  ==
C     == THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.         ==
C     ==--------------------------------------------------------------==
   30 M=12
C     ==--------------------------------------------------------------==
C     == SET UP INITIAL VALUES.                                       ==
C     ==--------------------------------------------------------------==
      LA=2
      IS=1
      IF=N
      DO 190 MLOOP=1,N
C     ==--------------------------------------------------------------==
C     ==  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE. ==
C     ==--------------------------------------------------------------==
        IFKA=IF-IS
        IF((IFKA+1).GT.M)GOTO 70
C     ==--------------------------------------------------------------==
C     == FINAL SORTING  ( A SIMPLE BUBBLE SORT )                      ==
C     ==--------------------------------------------------------------==
        IS1=IS+1
        DO 60 J=IS1,IF
          I=J
   40     IF(COUNT(I-1).LT.COUNT(I))GOTO 60
          IF(COUNT(I-1).GT.COUNT(I))GOTO 50
          IF(INDEX(I-1).LT.INDEX(I))GOTO 60
   50     AV=COUNT(I-1)
          COUNT(I-1)=COUNT(I)
          COUNT(I)=AV
          INT=INDEX(I-1)
          INDEX(I-1)=INDEX(I)
          INDEX(I)=INT
          I=I-1
          IF(I.GT.IS)GOTO  40
   60   CONTINUE
        LA=LA-2
        GOTO 170
C     ==--------------------------------------------------------------==
C     ==                *******  QUICKSORT  ********                  ==
C     == SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS  ==
C     == THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENTS ==
C     == HIGHEST ADDRESS.                                             ==
C     ==--------------------------------------------------------------==
   70   IY=(IS+IF)/2
        X=COUNT(IY)
        INTEST=INDEX(IY)
        COUNT(IY)=COUNT(IF)
        INDEX(IY)=INDEX(IF)
C     ==--------------------------------------------------------------==
C     == THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END ==
C     == OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE   ==
C     == OF X .                                                       ==
C     ==--------------------------------------------------------------==
        K=1
        IFK=IF
C     ==--------------------------------------------------------------==
C     == WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE ==
C     == INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AND INDICES AS   ==
C     == NECESSARY, UNTIL THEY MEET .                                 ==
C     ==--------------------------------------------------------------==
        DO 110 I=IS,IF
          IF(X.GT.COUNT(I))GOTO 110
          IF(X.LT.COUNT(I))GOTO 80
          IF(INTEST.GT.INDEX(I))GOTO 110
   80     IF(I.GE.IFK)GOTO 120
          COUNT(IFK)=COUNT(I)
          INDEX(IFK)=INDEX(I)
          K1=K
          DO 100 K=K1,IFKA
            IFK=IF-K
            IF(COUNT(IFK).GT.X)GOTO 100
            IF(COUNT(IFK).LT.X)GOTO 90
            IF(INTEST.LE.INDEX(IFK))GOTO 100
   90       IF(I.GE.IFK)GOTO 130
            COUNT(I)=COUNT(IFK)
            INDEX(I)=INDEX(IFK)
            GO TO 110
  100     CONTINUE
          GOTO 120
  110   CONTINUE
C     ==--------------------------------------------------------------==
C     == RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER  ==
C     == WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO ==
C     == 2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN OR EQUAL ==
C     == TO ANY ELEMENT IN THE SECOND PART, AND THEY MAY NOW BE SORTED==
C     == INDEPENDENTLY .                                              ==
C     ==--------------------------------------------------------------==
  120   COUNT(IFK)=X
        INDEX(IFK)=INTEST
        IP=IFK
        GOTO 140
  130   COUNT(I)=X
        INDEX(I)=INTEST
        IP=I
C     ==--------------------------------------------------------------==
C     ==  STORE THE LONGER SUBDIVISION IN WORKSPACE.                  ==
C     ==--------------------------------------------------------------==
  140   IF((IP-IS).GT.(IF-IP))GOTO 150
        MARK(LA)=IF
        MARK(LA-1)=IP+1
        IF=IP-1
        GOTO 160
  150   MARK(LA)=IP-1
        MARK(LA-1)=IS
        IS=IP+1
C     ==--------------------------------------------------------------==
C     == FIND THE LENGTH OF THE SHORTER SUBDIVISION.                  ==
C     ==--------------------------------------------------------------==
  160   LNGTH=IF-IS
        IF(LNGTH.LE.0)GOTO 180
C     ==--------------------------------------------------------------==
C     == IF IT CONTAINS MORE THAN ONE ELEMENT SUPPLY IT WITH WORKSPACE==
C     ==--------------------------------------------------------------==
        LA=LA+2
        GOTO 190
  170   IF(LA.LE.0)GOTO 200
C     ==--------------------------------------------------------------==
C     == OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT==
C     ==--------------------------------------------------------------==
  180   IF=MARK(LA)
        IS=MARK(LA-1)
  190 CONTINUE
C     ==--------------------------------------------------------------==
  200 RETURN
      END
C     ==================================================================
#endif
