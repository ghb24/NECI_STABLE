! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
      SUBROUTINE LINEUP(N,LST1,LST2,NBASIS,L1E,L2E,IFLAG,IPTOT)
      use util_mod, only: NECI_ICOPY
      INTEGER N, LST1(N), LST2(N)
      INTEGER LIST1(N),LIST2(N),L1E(NBASIS),L2E(NBASIS),
     &    LD(NBASIS)

      CALL NECI_ICOPY(N,LST1,1,LIST1,1)
      CALL NECI_ICOPY(N,LST2,1,LIST2,1)
C..
      CALL SCATTER_neci(N,NBASIS,LIST1,L1E)
      CALL SCATTER_neci(N,NBASIS,LIST2,L2E)
      IC=0
      DO I=1,NBASIS
         LD(I)=L1E(I)-L2E(I)
         IF(LD(I).GT.0) IC=IC+1
      ENDDO
      IF(IC.GT.2) THEN 
C        WRITE(6,*) ' MORE THAN TWO DIFFER'
        IFLAG=3
        RETURN
      ENDIF
      IF(IC.EQ.0) THEN 
C        WRITE(6,*) ' THE TWO ARE THE SAME'
        CALL ASCENDING_ORDER(N,LIST1,IP1)
        CALL ASCENDING_ORDER(N,LIST2,IP2)
        IPTOT=IP1*IP2
C        WRITE(6,*) (LIST1(I),I=1,N)
C        WRITE(6,*) (LIST2(I),I=1,N)
C        WRITE(6,*) ' TOTAL SIGN CHANGE:',IPTOT
        IFLAG=0
        RETURN
      ENDIF
      IF(IC.EQ.1) THEN 
         IPTOT=1
C..find out which one in list one differ
        DO I=1,NBASIS
          IF(LD(I).EQ.1) THEN
            DO J=1,N
              IF(LIST1(J).EQ.I) THEN  
                CALL CYCLIC_REORDER(N-J+1,LIST1(J),IP)
                IPTOT=IPTOT*IP
              ENDIF
            ENDDO
          ENDIF
          CONTINUE
          IF(LD(I).EQ.-1) THEN
            DO J=1,N
              IF(LIST2(J).EQ.I) THEN
                CALL CYCLIC_REORDER(N-J+1,LIST2(J),IP)
                IPTOT=IPTOT*IP
              ENDIF
            ENDDO
          ENDIF
       ENDDO
       CONTINUE
       CALL ASCENDING_ORDER(N-1,LIST1,IP1)
       CALL ASCENDING_ORDER(N-1,LIST2,IP2)
       IPTOT=IPTOT*IP1*IP2
C       WRITE(6,*) ' ONE DIFFERS:'
C       WRITE(6,*) (LIST1(I),I=1,N)
C       WRITE(6,*) (LIST2(I),I=1,N)
C       WRITE(6,*) ' TOTAL SIGN CHANGE:',IPTOT
       IFLAG=1
       RETURN
      ENDIF
C..
      IF(IC.EQ.2) THEN 
         IPTOT=1
C         WRITE(6,*) ' TWO DIFFER'
         IFLAG=2
C..find out which two in list one differ
        DO I=1,NBASIS
          IF(LD(I).EQ.1) THEN
            DO J=1,N
              IF(LIST1(J).EQ.I) THEN
                CALL CYCLIC_REORDER(N-J+1,LIST1(J),IP)
                IPTOT=IPTOT*IP
              ENDIF
            ENDDO
          ENDIF
          IF(LD(I).EQ.-1) THEN
            DO J=1,N
              IF(LIST2(J).EQ.I) THEN
                CALL CYCLIC_REORDER(N-J+1,LIST2(J),IP)
                IPTOT=IPTOT*IP
              ENDIF
            ENDDO
          ENDIF
       ENDDO
       CALL ASCENDING_ORDER(N-2,LIST1,IP1)
       CALL ASCENDING_ORDER(N-2,LIST2,IP2)
       IPTOT=IPTOT*IP1*IP2
       CALL ASCENDING_ORDER(2,LIST1(N-1),IP1)
       CALL ASCENDING_ORDER(2,LIST2(N-1),IP2)
       IPTOT=IPTOT*IP1*IP2
C       WRITE(6,*) (LIST1(I),I=1,N)
C       WRITE(6,*) (LIST2(I),I=1,N)
C       WRITE(6,*) ' TOTAL SIGN CHANGE:',IPTOT
      ENDIF
      RETURN
      END
C==============================================================
      SUBROUTINE SCATTER_neci(N,NBASIS,LIST,LE)
      DIMENSION LIST(N),LE(NBASIS)
      LE(1:NBASIS)=0
      DO I=1,N
        LE(LIST(I))=1
      ENDDO
      RETURN
      END
C===============================================================
      SUBROUTINE CYCLIC_REORDER(N,L,IP)
      INTEGER L(N)
      IS=L(1)
      DO I=2,N
         L(I-1)=L(I)
      ENDDO
      L(N)=IS
      IP=(-1)**(N-1)
      RETURN
      END
C===============================================================
      SUBROUTINE ASCENDING_ORDER(N,L,IP)
      INTEGER L(N)
      IP=1
 100  CONTINUE
      DO I=1,N-1 
         IF(L(I).GT.L(I+1)) THEN 
            IT=L(I)
            L(I)=L(I+1)
            L(I+1)=IT
            IP=IP*(-1)
            GOTO 100
         ENDIF
      ENDDO
      RETURN
      END
