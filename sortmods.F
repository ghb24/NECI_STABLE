C ==--------------------------------------------------------------==
      OPEN(10,FILE='PSI',STATUS='UNKNOWN')
      DO J=1,NEVAL
        IF(J.EQ.1) THEN
          WRITE(10,*) ' GROUND STATE COEFFICIENTS  '
        ELSE
          WRITE(10,*) ' COEFFICIENTS FOR EXCITED STATE NUMBER : ' , J
        ENDIF
        S=0.D0
        MS=0
        TMS=0.D0
        DO I=1,NDET
         IF(ABS(CG(I,J)).GT.1.D-15) THEN
            DO IEL=1,NEL
               WRITE(10,"I1,I1,2I2,2X$") (G(L,NM(1,I)),L=1,4)
            ENDDO
            WRITE(10,"F19.9,1X,I7") CG(I,J),I
         ENDIF
         S=S+CG(I,J)*CG(I,J)
        ENDDO
        WRITE(10,'(A,F19.5)') ' SQUARE OF COEFFICIENTS : ' , S
        WRITE(10,*)
      ENDDO
      CLOSE(10)
      RETURN
      END
