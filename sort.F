      PROGRAM NECICUBE
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IPARITY(3)
      LOGICAL TSPN,TALPHA,TPARITY,TCORR,TFODM,
     &	TREAD,TRHOOFR
C..Read in input parameters
      READ(5,*) NEL
      READ(5,*) TSPN
      READ(5,*) LMS
      READ(5,*) TPARITY,(IPARITY(I),I=1,3)
      READ(5,*) NMAXX,NMAXY,NMAXZ
      READ(5,*) NMSH
      READ(5,*) BOX,BOA,COA
      READ(5,*) TALPHA,ALPHA
      READ(5,*) ISTATE
      READ(5,*) NEVAL
      READ(5,*) B2L
      READ(5,*) NCYCLE
      READ(5,*) NBLK
      READ(5,*) NKRY
      READ(5,*) TRHOOFR
      READ(5,*) TREAD
      READ(5,*) TCORR
      READ(5,*) TFODM
      READ(5,*) IOBS,JOBS,KOBS 
C..
      CALL NECICB(NEL,TSPN,LMS,TPARITY,IPARITY,NMAXX,NMAXY,NMAXZ,
     &	NMSH,BOX,BOA,COA,TALPHA,ALPHA,ISTATE,NEVAL,B2L,NCYCLE,NBLK,
     &	NKRY,TRHOOFR,TREAD,TCORR,TFODM,IOBS,JOBS,KOBS)
C..
      CALL PMAXMEM
      CALL TIPRI 
      END
C ==------------------------------------------------------------------==
      SUBROUTINE NECICB(NEL,TSPN,LMS,TPARITY,IPARITY,NMAXX,NMAXY,
     &	NMAXZ,NMSH,BOX,BOA,COA,TALPHA,ALPHA,ISTATE,NEVAL,B2L,NCYCLE,
     &	NBLK,NKRY,TRHOOFR,TREAD,TCORR,TFODM,IOBS,JOBS,KOBS)
      IMPLICIT REAL*8 (A-H,O-Z) 
      INCLUDE 'cons.inc'
      INCLUDE 'irat.inc'
C..Pointers
      POINTER (IP_NMRKS,NMRKS),(IP_ARR,ARR),(IP_BRR,BRR),
     &	(IP_HAMIL,HAMIL),(IP_CKN,CKN),(IP_G1,G1),
     &	(IP_LAB,LAB),(IP_NROW,NROW),(IP_CK,CK),(IP_A,A),
     &	(IP_W,W),(IP_V,V),(IP_AM,AM),(IP_BM,BM),(IP_T,T),
     &	(IP_WT,WT),(IP_SCR,SCR),(IP_ISCR,ISCR),(IP_INDEX,INDEX),
     &	(IP_WH,WH),(IP_WORK2,WORK2),(IP_V2,V2),(IP_IKILL,IKILL),
     &	(IP_TKE,TKE)
C..These are the CUBE specific pointers
      POINTER (IP_FCK,FCK),(IP_COEFF,COEFF),(IP_ZIA,ZIA),
     &	(IP_RHO,RHO),(IP_SITAB,SITAB),(IP_XCHOLE,XCHOLE),
     &	(IP_PSIR,PSIR),(IP_DLINE,DLINE),(IP_SUMA,SUMA),
     &	(IP_SCRTCH,SCRTCH)
C..
      SAVE IP_NMRKS,IP_ARR,IP_BRR,IP_HAMIL,IP_CKN,
     &	IP_G1,IP_LAB,IP_NROW,IP_CK,IP_A,IP_W,IP_V,IP_AM,IP_BM,
     &	IP_T,IP_WT,IP_SCR,IP_ISCR,IP_INDEX,IP_WH,IP_WORK2,IP_V2,
     &	IP_IKILL,IP_TKE
C..Save the Cube specific pointers
      SAVE IP_FCK,IP_COEFF,IP_ZIA,IP_RHO,IP_SITAB,IP_XCHOLE,
     &	IP_PSIR,IP_DLINE,IP_SUMA,IP_SCRTCH
#ifdef POINTER8
      INTEGER*8 GET_ADDR
#else
      INTEGER GET_ADDR
#endif
      INTEGER NMRKS(NEL,*),G1(4,*),LAB(*),NROW(*),IKILL(*)
      INTEGER IPARITY(3),GC
      REAL*8 ARR(*),BRR(*),HAMIL(*),CKN(*),CK(*)
      REAL*8 A(*),W(*),V(*),AM(*),BM(*),T(*),WT(*),SCR(*)
      REAL*8 WH(*),WORK2(*),V2(*),RHO(*),SITAB(*),TKE(*)
C..Cube arrays
      REAL*8 ALAT(3),XCHOLE(*),SUMA(*)
      REAL*8 PSIR(*),DLINE(*),SCRTCH(*)
      COMPLEX*16 FCK(*)
      COMPLEX*16 COEFF(*)
      COMPLEX*16 ZIA(*)
      CHARACTER CPAR(3)*1,CPARITY*3
C..
      INTEGER INDEX(*),ISCR(*)
      LOGICAL TSPN,TALPHA,TPARITY,TCORR,TFODM,TREAD,TRHOOFR
C ==-------------------------------------------------------------------==
      CALL TISET('NECICUBE',ISUB)
C ==-------------------------------------------------------------------==
C..Input parameters
      WRITE(6,*) ' NUMBER OF ELECTRONS : ' , NEL
      WRITE(6,*) ' TSPN : ' , TSPN
      WRITE(6,*) ' GLOBAL MS : ' , LMS
C..
      IF(TPARITY) THEN
        WRITE(6,*) ' ******************************* '
        WRITE(6,*) ' PARITY IS ON '
        DO I=1,3
          IF(IPARITY(I).EQ.1) THEN
            CPAR(I)='G'
          ELSEIF(IPARITY(I).EQ.-1) THEN
            CPAR(I)='U'
          ELSE 
            STOP ' !!! PROBLEM WITH PARITY !!! '
          ENDIF
        ENDDO
        CPARITY=CPAR(1)//CPAR(2)//CPAR(3)
        WRITE(6,*) ' PARITY : ' , CPARITY
        WRITE(6,*) ' ******************************* '
      ELSE
        WRITE(6,*) ' PARITY IS OFF '
      ENDIF
C..
      NMAX=MAX(NMAXX,NMAXY,NMAXZ)
      NNR=NMSH*NMSH*NMSH
      WRITE(6,*) ' NMAXX : ' , NMAXX
      WRITE(6,*) ' NMAXY : ' , NMAXY
      WRITE(6,*) ' NMAXZ : ' , NMAXZ
      WRITE(6,*) ' NMSH : ' , NMSH 
C..
      WRITE(6,'(1X,A,F19.5)') ' BOX LENGTH : ' , BOX
      WRITE(6,'(1X,A,F19.5)') ' B/A : ' , BOA
      WRITE(6,'(1X,A,F19.5)') ' C/A : ' , COA
      ALAT(1)=BOX
      ALAT(2)=BOX*BOA
      ALAT(3)=BOX*COA
      OMEGA=ALAT(1)*ALAT(2)*ALAT(3)
      RS=(3.D0*OMEGA/(4.D0*PI*NEL))**THIRD
      WRITE(6,'(1X,A,F19.5)') ' VOLUME : ' , OMEGA
      WRITE(6,*) ' TALPHA : ' , TALPHA
      WRITE(6,'(1X,A,F19.5)') ' ALPHA : ' , ALPHA
      ALPHA=MIN(ALAT(1),ALAT(2),ALAT(3))*ALPHA
      WRITE(6,'(1X,A,F19.5)') ' SCALED ALPHA : ' , ALPHA
C..
      WRITE(6,*) ' ISTATE : ' , ISTATE 
      WRITE(6,*) ' NEVAL : ' , NEVAL
      WRITE(6,'(1X,A,E19.3)') ' B2LIMIT : ' , B2L
      WRITE(6,*) ' NCYCLE : ' , NCYCLE
      WRITE(6,*) ' NBLK : ' , NBLK 
      WRITE(6,*) ' NKRY : ' , NKRY
      WRITE(6,*) ' TREAD : ' , TREAD
      WRITE(6,*) ' TCORR : ' , TCORR
      IF(TCORR) THEN
        WRITE(6,*) ,
     &	' *** EXCHANGE-CORRELATION HOLE WILL BE CALCULATED *** ' 
      ENDIF
      WRITE(6,*) ' IOBS : ' , IOBS 
      WRITE(6,*) ' JOBS : ' , JOBS 
      WRITE(6,*) ' KOBS : ' , KOBS 
      IF(IOBS.GT.NMSH.OR.IOBS.LE.0.OR.JOBS.GT.NMSH.OR.JOBS.LE.0.
     &	OR.KOBS.GT.NMSH.OR.KOBS.LE.0) THEN
        STOP ' !!! REFERENCE PARTICLE NOT IN BOX !!! '
      ENDIF
C..
      LEN=2*NMAXX*NMAXY*NMAXZ
      CALL MEMORY(IP_ARR,LEN,'ARR')
      CALL AZZERO(ARR,LEN)
      CALL MEMORY(IP_BRR,LEN,'BRR')
      CALL AZZERO(BRR,LEN)
      CALL MEMORY(IP_G1,4*LEN,'G1')
      CALL IAZZERO(G1,4*LEN)
C..Calculate number of basis functions
      IG=0
      DO I=1,NMAXX
        DO J=1,NMAXY
          DO K=1,NMAXZ
            DO L=-1,1,2
C..
              IG=IG+1
              SUM=(BOX**2)*((I*I/ALAT(1)**2)+(J*J/ALAT(2)**2)+
     &			(K*K/ALAT(3)**2))
              ARR(IG)=SUM
              BRR(IG)=IG
C..These are the quantum numbers: n,l,m and sigma
              G1(1,IG)=I
              G1(2,IG)=J
              G1(3,IG)=K
              G1(4,IG)=L
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C..Check to see if all's well
      WRITE(6,*) ' NUMBER OF BASIS FUNCTIONS : ' , IG 
      IF(LEN.NE.IG) STOP ' LEN NE IG ' 
      NHG=IG
C..Now we sort them using SORT2 and then SORT
      CALL SORT2(NHG,ARR,BRR)
      CALL SORT(NHG,BRR)
      DO I=1,NHG
        WRITE(6,'(5I7,F19.5)') INT(BRR(I)),G1(1,INT(BRR(I))),
     &		G1(2,INT(BRR(I))),G1(3,INT(BRR(I))),
     &		G1(4,INT(BRR(I)))
      ENDDO
C..Now we see how many determinants we need
      IF(NHG.GT.170) THEN
C..This fix is to stop floating overflow as FACTRL(NHG.GT.170) crashes
        WRITE(6,*) ' !!! NHG TOO LARGE ASSUMING NDET EQ. 1 !!! '
        NDET=1
        GOTO 200
      ENDIF
      NDET=FACTRL(NHG)/(FACTRL(NEL)*FACTRL(NHG-NEL)) 
200   CONTINUE
      WRITE(6,*) ' NUMBER OF DETERMINANTS : ' , NDET
C..Need to determine the determinants 
      CALL GNDTS(NEL,NHG,BRR,NMAXX,NMAXY,NMAXZ,NMRKS,.TRUE.,
     &		NDET,G1,TSPN,LMS,TPARITY,IPARITY,II)
C..
      CALL MEMORY(IP_NMRKS,NEL*II,'NMRKS')
      CALL IAZZERO(NMRKS,NEL*II)
C..
      CALL GNDTS(NEL,NHG,BRR,NMAXX,NMAXY,NMAXZ,NMRKS,.FALSE.,
     &		NDET,G1,TSPN,LMS,TPARITY,IPARITY,II)
C..
      NDET=II    
      WRITE(6,*) ' NUMBER OF DETS. EQ NUMBER COUNTED : ' , NDET,II
C ==----------------------------------------------------------------==
C..Set up memory for c's, nrow and the label
      CALL MEMORY(IP_CKN,NDET*NEVAL,'CKN')
      CALL AZZERO(CKN,NDET*NEVAL)
C..
      CALL MEMORY(IP_CK,NDET*NEVAL,'CK')
      CALL AZZERO(CK,NDET*NEVAL)
C..
      IF(TREAD) THEN
        CALL READ_PSI(BOX,BOA,COA,NDET,NEVAL,
     &		NMAXX,NMAXY,NMAXZ,NEL,CK)
        GOTO 100
      ENDIF
      CALL MEMORY(IP_NROW,NDET,'NROW')
      CALL IAZZERO(NROW,NDET)
C..Need to initialise the Fourier arrays
      CALL MEMORY(IP_FCK,2*NMSH**3,'FCK')
      CALL MEMORY(IP_COEFF,2*(3*NMSH+48),'COEFF')
      CALL MEMORY(IP_ZIA,2*(NMSH+1)*NMAX*NMAX,'ZIA')
C..
      CALL INITFOU(NMSH,FCK,COEFF,NMAX,ZIA,ALAT,TALPHA,
     &	ALPHA,OMEGA)
C..We need to measure HAMIL and LAB first 
      ICMAX=1
      CALL DETHAM(NDET,NEL,NMRKS,NMAXX,NMAXY,NMAXZ,
     &	NHG,HAMIL,G1,LAB,NROW,.TRUE.,
     &	NMSH,FCK,ZIA,NMAX,ALAT,ICMAX,GC,.FALSE.)
      WRITE(6,*) ' FINISHED COUNTING '
C..Now we know size, allocate memory to HAMIL and LAB
      LENHAMIL=GC
      CALL MEMORY(IP_HAMIL,LENHAMIL,'HAMIL')
      CALL AZZERO(HAMIL,LENHAMIL)
C..
      CALL MEMORY(IP_LAB,LENHAMIL/IRAT+1,'LAB')
      CALL IAZZERO(LAB,LENHAMIL)
C..Now we store HAMIL and LAB 
      CALL DETHAM(NDET,NEL,NMRKS,NMAXX,NMAXY,NMAXZ,
     &	NHG,HAMIL,G1,LAB,NROW,.FALSE.,
     &	NMSH,FCK,ZIA,NMAX,ALAT,ICMAX,GC,.FALSE.)
C..Things needed for Friesner-Pollard diagonalisation
      NKRY1=NKRY+1
      NBLOCK=MIN(NEVAL,NBLK)
      LSCR=MAX(NDET*NEVAL,8*NBLOCK*NKRY)
      LISCR=6*NBLOCK*NKRY
C..
      WRITE(*,'(/,/,8X,64(1H*))')
      WRITE(*,'(7X," *",62X,"*")')
      WRITE(*,'(7X," *",19X,A,18X,"*")') ' LANCZOS DIAGONALISATION '
      WRITE(*,'(7X," *",62X,"*")')
      WRITE(*,'(7X,1X,64(1H*))')
C..Set up memory for FRSBLKH
      CALL MEMORY(IP_A,NEVAL*NEVAL,'A')
      CALL AZZERO(A,NEVAL*NEVAL)
C..
      CALL MEMORY(IP_W,NEVAL,'W')
      CALL AZZERO(W,NEVAL)
C..
      CALL MEMORY(IP_V,NDET*NBLOCK*NKRY1,'V')
      CALL AZZERO(V,NDET*NBLOCK*NKRY1)
C..
      CALL MEMORY(IP_AM,NBLOCK*NBLOCK*NKRY1,'AM') 
      CALL AZZERO(AM,NBLOCK*NBLOCK*NKRY1)
C..
      CALL MEMORY(IP_BM,NBLOCK*NBLOCK*NKRY,'BM')
      CALL AZZERO(BM,NBLOCK*NBLOCK*NKRY)
C..
      CALL MEMORY(IP_T,3*NBLOCK*NKRY*NBLOCK*NKRY,'T')
      CALL AZZERO(T,3*NBLOCK*NKRY*NBLOCK*NKRY)
C..
      CALL MEMORY(IP_WT,NBLOCK*NKRY,'WT')
      CALL AZZERO(WT,NBLOCK*NKRY)
C..
      CALL MEMORY(IP_SCR,LSCR,'SCR')
      CALL AZZERO(SCR,LSCR)
      CALL MEMORY(IP_ISCR,LISCR,'ISCR')
      CALL IAZZERO(ISCR,LISCR)
      CALL MEMORY(IP_INDEX,NEVAL,'INDEX')
      CALL IAZZERO(INDEX,NEVAL)
C..
      CALL MEMORY(IP_WH,NDET,'WH')
      CALL AZZERO(WH,NDET)
      CALL MEMORY(IP_WORK2,3*NDET,'WORK2')
      CALL AZZERO(WORK2,3*NDET)
      CALL MEMORY(IP_V2,NDET*NEVAL,'V2')
      CALL AZZERO(V2,NDET*NEVAL)
C..Lanczos iterative diagonalising routine
      CALL FRSBLKH(NDET,ICMAX,NEVAL,HAMIL,LAB,CK,CKN,
     &	NKRY,NKRY1,NBLOCK,NROW,LSCR,LISCR,A,W,V,AM,BM,T,WT,
     &	SCR,ISCR,INDEX,WH,WORK2,V2,NCYCLE,B2L)
      CALL DSCAL(NEVAL,-1.D0,W,1) 
C..
      CALL MEMORY(IP_TKE,NEVAL,'TKE')
C..Free HAMIL AND LAB memory as we no longer need them
      CALL FREEM(IP_HAMIL)
      CALL FREEM(IP_LAB)
C..
      CALL CFF_CHCK(NDET,NEVAL,NMRKS,NMAXX,NMAXY,NMAXZ,NEL,G1,CK,
     &	ALAT,TKE)
      CALL WRITE_PSI(BOX,BOA,COA,NDET,NEVAL,
     &	NMAXX,NMAXY,NMAXZ,NEL,CK)
      WRITE(6,*) ,
     &	'	==--------------------------------------------------== '
      WRITE(6,'(A5,5X,A15,1X,A18,1x,A20)') 'STATE',
     &	'KINETIC ENERGY', 'COULOMB ENERGY', 'TOTAL ENERGY'
      DO IN=1,NEVAL
      	WRITE(6,'(I3,2X,3(F19.11,2x))')
     &		IN,TKE(IN),W(IN)-TKE(IN),W(IN)
      ENDDO
      WRITE(6,*) ,
     &	'	==--------------------------------------------------== '
100   CONTINUE
      IF(TRHOOFR) THEN
C..Generate memory for RHO and SITAB
        CALL MEMORY(IP_RHO,NMSH*NMSH*NMSH,'RHO')
        CALL MEMORY(IP_SITAB,NMSH*NMAX,'SITAB')
        CALL MEMORY(IP_SCRTCH,NMSH*NMSH*NMSH,'SCRTCH')
C..Calculate RHOOFR
        CALL RHOOFR(NHG,CK,G1,RHO,NMSH,SITAB,NMAX,NMRKS,NEL,
     &	  NDET,NEVAL,RS,ALAT,OMEGA,SCRTCH)
C..
        CALL MEMORY(IP_DLINE,NMSH,'DLINE')
        CALL AZZERO(DLINE,NMSH)
C..Calculate RHOOFR in certain directions
C..001
        CALL PLANARAV(RHO,NMSH,DLINE,0,0,1,SPAC,ALAT) 
        CALL WRITE_LINE(8,'RHOAV001',DLINE,1,NMSH,-1,-1,-1,SPAC,RS)
C..100
        CALL PLANARAV(RHO,NMSH,DLINE,1,0,0,SPAC,ALAT)
        CALL WRITE_LINE(8,'RHOAV100',DLINE,1,NMSH,-1,-1,-1,SPAC,RS)
C..010
        CALL PLANARAV(RHO,NMSH,DLINE,0,1,0,SPAC,ALAT)
        CALL WRITE_LINE(8,'RHOAV010',DLINE,1,NMSH,-1,-1,-1,SPAC,RS)
        IF(TCORR) THEN
C..Now generate memory for XCHOLE 
          CALL MEMORY(IP_XCHOLE,NMSH*NMSH*NMSH,'XCHOLE')
          CALL MEMORY(IP_PSIR,2*NMSH+1,'PSIR')
          CALL AZZERO(PSIR,2*NMSH+1)
C..
          IXD=1
          IYD=0
          IZD=0
          SPAC=0.D0
          CALL GEN_XCHOLE(CK,PSIR,IOBS,JOBS,KOBS,G1,SITAB,
     &		NMAX,NMSH,NHG,IXD,IYD,IZD,RHO,.TRUE.,
     &	  	XCHOLE,SPAC,ALAT,OMEGA,NMRKS,NDET,NEVAL,NEL)
          CALL WRITE_RHO(10,'COMPXCHOLE',XCHOLE,NMSH,NMSH,NMSH,
     &		ALAT,.FALSE.,.TRUE.,RS)
C..
          CALL XCHOLES(CK,PSIR,IOBS,JOBS,KOBS,G1,SITAB,
     &		NMAX,NMSH,NHG,RHO,XCHOLE,SPAC,RS,ALAT,DLINE,
     &	  	OMEGA,NMRKS,NDET,NEL,NEVAL)
C          OPEN(11,FILE='XCEDENS',STATUS='UNKNOWN')
C          DO IZ=1,NMSH
C            CALL XCENERGY_DENSITY(CK,PSIR,DLINE,IOBS,JOBS,IZ,G1,
C     &	      SITAB,NMAX,NMSH,NHG,RHO,XCHOLE,ALAT,OMEGA,
C     &	      NMRKS,NDET,NEVAL,NEL)
C          ENDDO
C          CLOSE(11)
        ENDIF 
      ENDIF
      IF(TFODM) THEN
        CALL MEMORY(IP_SUMA,NMAXX*NMAXY*NMAXZ,'SUMA')
        CALL AZZERO(SUMA,NMAXX*NMAXY*NMAXZ)
        CALL FODMAT(NEL,NHG,NDET,NEVAL,ISTATE,NMRKS,G1,CK,
     &		NMAXX,NMAXY,NMAXZ,SUMA)
      ENDIF
      CALL MEMORY_CHECK
C ==-------------------------------------------------------------------==
      CALL TIHALT('NECICUBE',ISUB)
C ==-------------------------------------------------------------------==
      RETURN
      END
C ==-------------------------------------------------------------------==
      SUBROUTINE CFF_CHCK(NDET,NEVAL,NM,NMAXX,NMAXY,NMAXZ,NEL,G,CG,
     &	ALAT,TKE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CG(NDET,NEVAL),NM(NEL,*)
      DIMENSION ALAT(3),TKE(NEVAL)
      INTEGER G(4,2*NMAXX*NMAXY*NMAXZ)
C..Calculate the expectation value of the kinetic energy
C..<Psi|T|Psi>
      PI=ACOS(-1.D0)
      DO IN=1,NEVAL
        TKE(IN)=0.D0
        DO I=1,NDET
          SUM=0.D0
          DO J=1,NEL
            AUX=((ALAT(1)**2)*((G(1,NM(J,I))**2)/(ALAT(1)**2)+
     &		(G(2,NM(J,I))**2)/(ALAT(2)**2)+
     &		(G(3,NM(J,I))**2)/(ALAT(3)**2)))
            SUM=SUM+AUX
          ENDDO
C..Cube multiplier
          CST=PI*PI/(2.D0*ALAT(1)*ALAT(1))
          SUM=CST*SUM 
          TKE(IN)=TKE(IN)+SUM*CG(I,IN)*CG(I,IN)
        ENDDO
      ENDDO
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
#ifdef NEL2
           WRITE(10,'(2(I1,I1,2I2,2X),F19.9,1X,I7)')
#endif
#ifdef NEL3
	   WRITE(10,'(3(I1,I1,2I2,2X),F19.9,1X,I7)')
#endif
#ifdef NEL4
           WRITE(10,'(4(I1,I1,2I2,2X),F19.9,1X,I7)')
#endif
#ifdef NEL5
           WRITE(10,'(5(I1,I1,2I2,2X),F19.9,1X,I7)')
#endif
#ifdef NEL6
           WRITE(10,'(6(I1,I1,2I2,2X),F19.9,1X,I7)')
#endif
     &	(G(L,NM(1,I)),L=1,4)
     &	,(G(K,NM(2,I)),K=1,4)
#if NEL3 || NEL4 || NEL5 || NEL6
     &	,(G(K,NM(3,I)),K=1,4)
#endif
#if NEL4 || NEL5 || NEL6
     &	,(G(K,NM(4,I)),K=1,4)
#endif
#if NEL5 || NEL6
     &	,(G(K,NM(5,I)),K=1,4)
#endif
#if NEL6
     &	,(G(K,NM(6,I)),K=1,4)
#endif
     &	,CG(I,J),I
 	  ENDIF
          S=S+CG(I,J)*CG(I,J)
        ENDDO
        WRITE(10,'(A,F19.5)') ' SQUARE OF COEFFICIENTS : ' , S
        WRITE(10,*)
      ENDDO
      CLOSE(10)
C      DO I=1,NDET
C        WRITE(6,*) I,(NM(J,I),J=1,NEL)
C      ENDDO
      RETURN
      END
