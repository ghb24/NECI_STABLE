C     ==================================================================
      SUBROUTINE READSR(STRING,IIN,IOUT,FPTNUM,ERROR)
C     ==--------------------------------------------------------------==
C     == READ A FULL PRECISION (64 BIT) FLOATING POINT NUMBER         ==
C     == FROM A STRING                                                ==
C     == NOTE: ACCEPTED FORM IS                                       ==
C     ==    (SIGN)___.___E(SIGN)__/(SIGN)___.___E(SIGN)__             ==
C     ==    THE EXPONENT PORTION IS OPTIONAL AS IS THE DECIMAL POINT. ==
C     ==    A FIELD STARTED BY E IS ASSUMED TO HAVE MANTISSA OF 1.    ==
C     ==--------------------------------------------------------------==
C     == INPUT:                                                       ==
C     ==   STRING*(*) Character string                                ==
C     ==   IIN  Index of the first read character in STRING           ==
C     == OUTPUT:                                                      ==
C     ==   IOUT Index of the last read character (STRING(IIN:IOUT))   ==
C     ==   FPTNUM Full precision floating point number                ==
C     ==   ERROR .TRUE. if error occured                              ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      CHARACTER STRING*(*)
      INTEGER   IIN,IOUT
      REAL*8    FPTNUM
      LOGICAL   ERROR
C     Variables
      REAL*8    ZERO,ONE,TEN
      PARAMETER (ZERO=0.D0,ONE=1.D0,TEN=10.D0)
      CHARACTER IV*80
      INTEGER   I1,I2,LENGTH,INV,ISEXP,NEXP,J,INDEX
      REAL*8    V1,SIGN,SCALE
      LOGICAL   DIGIT,EXP,POINT,TSIGN,SLASH
C     ==--------------------------------------------------------------==
      CALL XSTRING(STRING(IIN:LEN(STRING)),I1,I2)
      IOUT=IIN+I2
      IV=STRING(IIN+I1-1:IIN+I2-1)
      LENGTH=I2-I1+1
      ERROR=.FALSE.
      INV    = 0
      FPTNUM = ZERO
      IF(LENGTH.EQ.0) THEN
        ERROR=.TRUE.
        RETURN
      ENDIF
      V1=ZERO
      SIGN=ONE
      SCALE=ONE
      DIGIT = .FALSE.
      EXP   = .FALSE.
      POINT = .FALSE.
      TSIGN = .FALSE.
      SLASH = .FALSE.
      ISEXP = 1
      NEXP  = 0
C     ==--------------------------------------------------------------==
      DO 150 J=1,LENGTH
C       
C       --- IDENTIFY THE J-TH CHARACTER ---
C       
        IF(IV(J:J).EQ.'0') THEN
          INDEX=0
        ELSEIF(IV(J:J).EQ.'1') THEN
          INDEX=1
        ELSEIF(IV(J:J).EQ.'2') THEN
          INDEX=2
        ELSEIF(IV(J:J).EQ.'3') THEN
          INDEX=3
        ELSEIF(IV(J:J).EQ.'4') THEN
          INDEX=4
        ELSEIF(IV(J:J).EQ.'5') THEN
          INDEX=5
        ELSEIF(IV(J:J).EQ.'6') THEN
          INDEX=6
        ELSEIF(IV(J:J).EQ.'7') THEN
          INDEX=7
        ELSEIF(IV(J:J).EQ.'8') THEN
          INDEX=8
        ELSEIF(IV(J:J).EQ.'9') THEN
          INDEX=9
        ELSEIF(IV(J:J).EQ.'/') THEN
          INDEX=10
        ELSEIF(IV(J:J).EQ.'.') THEN
          INDEX=11
        ELSEIF(IV(J:J).EQ.'+') THEN
          INDEX=12
        ELSEIF(IV(J:J).EQ.'-') THEN
          INDEX=13
        ELSEIF(IV(J:J).EQ.'E') THEN
          INDEX=14
        ELSEIF(IV(J:J).EQ.'e') THEN
          INDEX=14
        ELSEIF(IV(J:J).EQ.'D') THEN
          INDEX=15
        ELSEIF(IV(J:J).EQ.'d') THEN
          INDEX=15
        ELSE
C         ILLEGAL CHARACTER FOUND
          ERROR=.TRUE.
          RETURN
        ENDIF
C       
C       --- TEST IF THIS CHARACTER IS A NON-NUMERAL ---
C       
        IF(INDEX.LE.9) THEN
C         FIND A DIGIT
          DIGIT=.TRUE.
          TSIGN=.FALSE.
C         TEST IF PART OF EXPONENT
          IF(EXP) THEN
C           ADD DIGIT TO EXPONENT
            NEXP=10*NEXP+INDEX
            GO TO 150
          ENDIF
C         ADD DIGIT TO MANTISSA
          INV=INV+1
          IF(POINT) THEN
            SCALE=SCALE/TEN
            FPTNUM=FPTNUM+DBLE(INDEX)*SCALE
          ELSE
            FPTNUM=TEN*FPTNUM+DBLE(INDEX)
          END IF
          GO TO 150
        ENDIF
C       
C       --- PROCESS NON-NUMERALS CHARACTERS ---
        IF(INDEX.EQ.10) THEN
          GOTO 120
        ELSEIF(INDEX.EQ.11) THEN
          GOTO 125
        ELSEIF(INDEX.EQ.12) THEN
          GOTO 130
        ELSEIF(INDEX.EQ.13) THEN
          GOTO 135
        ELSEIF(INDEX.EQ.14.OR.INDEX.EQ.15) THEN
          GOTO 145
        ENDIF
C       
C       --- SLASH DETECTED (FIELD EXPRESSED AS A FRACTION) NUMERATOR COMP
C       
 120    CONTINUE
        IF(SLASH) THEN
          ERROR=.TRUE.
          RETURN
        ENDIF
        SLASH=.TRUE.
        V1=FPTNUM*SIGN*(TEN**(ISEXP*NEXP))
        FPTNUM=ZERO
        TSIGN=.FALSE.
        DIGIT=.FALSE.
        INV=0
        IF(V1.EQ.ZERO) GO TO 155
        SIGN=ONE
        ISEXP=1
        EXP=.FALSE.
        SCALE=ONE
        POINT=.FALSE.
        NEXP=0
        GO TO 150
C       
C       --- DECIMAL POINT DETECTED ---
C       
 125    CONTINUE
        IF(EXP) THEN
          ERROR=.TRUE.
          RETURN
        ENDIF
        POINT=.TRUE.
        GO TO 150
C       
C       --- PLUS TSIGN DETECTED, TEST IF START OF MANTISSA OR EXPONENT ---
C
 130    CONTINUE
        IF(TSIGN) THEN
          ERROR=.TRUE.
          RETURN
        ELSE
          TSIGN=.TRUE.
        ENDIF
        IF(EXP) THEN
          IF(NEXP.NE.0) THEN
            ERROR=.TRUE.
            RETURN
          ENDIF
          ISEXP=1
        ELSE
          IF(INV.NE.0) THEN
            ERROR=.TRUE.
            RETURN
          ENDIF
          SIGN=ONE
        ENDIF
        GOTO 150
C       
C       --- MINUS TSIGN, TEST IF START OF MANTISSA OF EXPONENT ---
C       
 135    CONTINUE
        IF(TSIGN) THEN
          ERROR=.TRUE.
          RETURN
        ELSE
          TSIGN=.TRUE.
        ENDIF
        IF(EXP) THEN
          IF(NEXP.NE.0) THEN
            ERROR=.TRUE.
            RETURN
          ENDIF
          ISEXP=-1
        ELSE
          IF(INV.NE.0) THEN
            ERROR=.TRUE.
            RETURN
          ENDIF
          SIGN=-ONE
        ENDIF
        GOTO 150
C       
C       --- E, D OR EMBEDDED + OR - STARTS EXPONENT FIELD ---
C       
 145    CONTINUE
        IF(EXP) THEN
          ERROR=.TRUE.
          RETURN
        ENDIF
        EXP=.TRUE.
C       A FIELD STARTED BY E IS ASSUMED TO HAVE MANTISSA OF 1.
        IF(INV.EQ.0) FPTNUM=ONE
        TSIGN=.FALSE.
        GO TO 150
 150  CONTINUE
C     ==--------------------------------------------------------------==
      IF(DIGIT) THEN
        FPTNUM=FPTNUM*SIGN*(TEN**(ISEXP*NEXP))
      ELSE
C       NO DIGIT DETECTED
        ERROR=.TRUE.
        RETURN
      ENDIF
C     
C     --- THE NUMERATOR IS FINISHED, TEST FOR NO DENOMINATOR ---
C     
      IF(V1.EQ.ZERO) GO TO 155
      IF(FPTNUM.EQ.ZERO) FPTNUM=ONE
      FPTNUM=V1/FPTNUM
 155  CONTINUE
      IF(DIGIT) THEN
        V1=FPTNUM
      ELSE
C       NO DIGIT DETECTED
        ERROR=.TRUE.
        RETURN
      ENDIF 
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE READSI(STRING,IIN,IOUT,INTNUM,ERROR)
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      CHARACTER STRING*(*)
      INTEGER INTNUM,IIN,IOUT
      LOGICAL ERROR
C     Variables
      REAL*8 FPTNUM
C     ==--------------------------------------------------------------==
      CALL READSR(STRING,IIN,IOUT,FPTNUM,ERROR)
      INTNUM=NINT(FPTNUM)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE XSTRING(STRING,IA,IE)
C     ==--------------------------------------------------------------==
C     == Give index of the first and last non-blank character         ==
C     == blank characters are " " and CHAR(9) (tab)                   ==
C     == more blank characters are char(13) (cr) and CHAR(10) (newl)  ==
C     == STRING (INPUT)                                               ==
C     == IA,IE  (OUTPUT)                                              ==
C     ==--------------------------------------------------------------==
      IMPLICIT NONE
C     Arguments
      CHARACTER STRING*(*)
      INTEGER IA,IE
C     Variables
      INTEGER SLEN,I
C     ==--------------------------------------------------------------==
      SLEN=LEN(STRING)
      IA=1
      DO I=1,SLEN
         IF(STRING(I:I).NE.' '.AND.
     &       STRING(I:I).NE.CHAR(9).AND.      !tab character
     &       STRING(I:I).NE.CHAR(10).AND.     !newline character
     &       STRING(I:I).NE.CHAR(13) ) THEN   !carriage return character
           IA=I
           GOTO 10
         ENDIF
      ENDDO
 10   CONTINUE
      DO I=IA,SLEN
         IF(STRING(I:I).EQ.' '.OR.
     &       STRING(I:I).EQ.CHAR(0).OR. !\0 character
     &       STRING(I:I).EQ.CHAR(9)) THEN !tab character
           IE=I-1
           GOTO 20
         ENDIF
      ENDDO
      IE=SLEN
 20   CONTINUE
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
