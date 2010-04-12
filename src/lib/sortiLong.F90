! These routines use INTEGER(SELECTED_INT_KIND(18)), and so are in a seperate F90 file.
!.. Based on SortI, this will sort a list of INTEGER(KIND=S_I_K(18)), and order the next three integer arrays and logical array according to the first.
! This is used for when we also have to take the excitation level with the hash.
      SUBROUTINE SORT4I1LLong(N,RA,RB,RC,RD,RE)
      IMPLICIT NONE
      INTEGER N,I,L,IR,J
      INTEGER(KIND=SELECTED_INT_KIND(18)) :: RA(N)
      INTEGER(KIND=SELECTED_INT_KIND(18)) :: RRA
      LOGICAL :: RE(N)
      LOGICAL :: RRE
      INTEGER RB(N)
      INTEGER RC(N),RD(N)
      INTEGER RRB,RRC,RRD
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
          RRC=RC(L)
          RRD=RD(L)
          RRE=RE(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          RRC=RC(IR)
          RC(IR)=RC(1)
          RRD=RD(IR)
          RD(IR)=RD(1)
          RRE=RE(IR)
          RE(IR)=RE(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RC(1)=RRC
            RD(1)=RRD
            RE(1)=RRE
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
            RC(I)=RC(J)
            RD(I)=RD(J)
            RE(I)=RE(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
        RC(I)=RRC
        RD(I)=RRD
        RE(I)=RRE
      GO TO 10
      END


! These routines use INTEGER(SELECTED_INT_KIND(18)), and so are in a seperate F90 file.
!.. Based on SortI, this will sort a list of INTEGER(KIND=S_I_K(18)), and order the next two integer arrays and logical array according to the first.
      SUBROUTINE SORT3I1LLong(N,RA,RB,RC,RD)
      IMPLICIT NONE
      INTEGER N,I,L,IR,J
      INTEGER(KIND=SELECTED_INT_KIND(18)) :: RA(N)
      INTEGER(KIND=SELECTED_INT_KIND(18)) :: RRA
      LOGICAL :: RD(N)
      LOGICAL :: RRD
      INTEGER RB(N)
      INTEGER RC(N)
      INTEGER RRB,RRC
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
          RRC=RC(L)
          RRD=RD(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          RRC=RC(IR)
          RC(IR)=RC(1)
          RRD=RD(IR)
          RD(IR)=RD(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RC(1)=RRC
            RD(1)=RRD
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
            RC(I)=RC(J)
            RD(I)=RD(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
        RC(I)=RRC
        RD(I)=RRD
      GO TO 10
      END

! These routines use INTEGER(SELECTED_INT_KIND(18)), and so are in a seperate F90 file.
!.. Based on SortI, this will sort a list of INTEGER(KIND=S_I_K(18)), and order the next three integer arrays according to the first.
      SUBROUTINE SORT4ILong(N,RA,RB,RC,RD)
      IMPLICIT NONE
      INTEGER N,I,L,IR,J
      INTEGER(KIND=SELECTED_INT_KIND(18)) :: RA(N)
      INTEGER(KIND=SELECTED_INT_KIND(18)) :: RRA
      INTEGER :: RD(N)
      INTEGER :: RRD
      INTEGER RB(N)
      INTEGER RC(N)
      INTEGER RRB,RRC
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
          RRC=RC(L)
          RRD=RD(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          RRC=RC(IR)
          RC(IR)=RC(1)
          RRD=RD(IR)
          RD(IR)=RD(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RC(1)=RRC
            RD(1)=RRD
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
            RC(I)=RC(J)
            RD(I)=RD(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
        RC(I)=RRC
        RD(I)=RRD
      GO TO 10
      END

!.. Based on SortI, this will sort a list of INTEGERS, and order the next integer array, integer(KIND=S_I_K(18)) array and logical array according to the first.
      SUBROUTINE SORT2IILongL(N,RA,RB,RC,RD)
      IMPLICIT NONE
      INTEGER N,I,L,IR,J
      INTEGER(KIND=SELECTED_INT_KIND(18)) :: RC(N)
      INTEGER(KIND=SELECTED_INT_KIND(18)) :: RRC
      LOGICAL :: RD(N)
      LOGICAL :: RRD
      INTEGER RB(N)
      INTEGER RA(N)
      INTEGER RRB,RRA
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
          RRC=RC(L)
          RRD=RD(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          RRC=RC(IR)
          RC(IR)=RC(1)
          RRD=RD(IR)
          RD(IR)=RD(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RC(1)=RRC
            RD(1)=RRD
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
            RC(I)=RC(J)
            RD(I)=RD(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
        RC(I)=RRC
        RD(I)=RRD
      GO TO 10
      END

!.. Based on SortI, this will sort a list of INTEGERS, and order the next integer array, integer(KIND=S_I_K(18)) array and further integer array according to the first.
      SUBROUTINE SORT2IILongI(N,RA,RB,RC,RD)
      IMPLICIT NONE
      INTEGER N,I,L,IR,J
      INTEGER(KIND=SELECTED_INT_KIND(18)) :: RC(N)
      INTEGER(KIND=SELECTED_INT_KIND(18)) :: RRC
      INTEGER :: RD(N)
      INTEGER :: RRD
      INTEGER RB(N)
      INTEGER RA(N)
      INTEGER RRB,RRA

      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
          RRC=RC(L)
          RRD=RD(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          RRC=RC(IR)
          RC(IR)=RC(1)
          RRD=RD(IR)
          RD(IR)=RD(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RC(1)=RRC
            RD(1)=RRD
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
            RC(I)=RC(J)
            RD(I)=RD(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
        RC(I)=RRC
        RD(I)=RRD
      GO TO 10
      END

