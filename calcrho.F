      INTEGER FUNCTION ICMPDETS(D1,D2,NEL)
         IMPLICIT NONE
         INTEGER NEL,I
         INTEGER D1(1:NEL),D2(1:NEL),ICMPCSFS
         LOGICAL ISCSF
         IF(ISCSF(D1,NEL).OR.ISCSF(D2,NEL)) THEN
            ICMPDETS=ICMPCSFS(D1,D2,NEL)
            RETURN
         ENDIF
         DO I=1,NEL
            IF(D1(I).LT.D2(I)) THEN
               ICMPDETS=-1
               RETURN
            ELSE IF(D1(I).GT.D2(I)) THEN
               ICMPDETS=1
               RETURN
            ENDIF
         ENDDO
         ICMPDETS=0
         RETURN
      END

      INTEGER FUNCTION IGETEXCITLEVEL(NI,NJ,NEL)
         INTEGER I,J,NEL,NI(NEL),NJ(NEL),IC
         LOGICAL ISCSF
         IF(ISCSF(NI,NEL).OR.ISCSF(NJ,NEL)) THEN
            CALL GETCSFEXCITLEVEL(NI,NJ,NEL,IGETEXCITLEVEL)
C            IGETEXCITLEVEL=2
            RETURN
         ENDIF
         IC=0
         DO I=1,NEL
            DO J=1,NEL
               IF(NI(I).EQ.NJ(J)) IC=IC+1
            ENDDO
         ENDDO
         IGETEXCITLEVEL=NEL-IC
         RETURN
      END

C.. Very unnecessarily complex, but would be faster for large numbers of
C.. electrons.

      INTEGER FUNCTION IGETEXCITLEVEL_(NI,NJ,NEL)
         INTEGER I,J,NEL,NI(NEL),NJ(NEL),IC
         IC=0
         I=1
         J=1
C.. We only count differences from I to J
         DO WHILE(I.LE.NEL.AND.J.LE.NEL)
            DO WHILE(NI(I).LT.NJ(J).AND.I.LE.NEL)
               I=I+1
               IC=IC+1
            ENDDO
            DO WHILE(NI(I).GT.NJ(J).AND.I.LE.NEL.AND.J.LE.NEL)
               J=J+1
C               IC=IC+1
            ENDDO
            IF(NI(I).EQ.NJ(J)) THEN
               I=I+1
               J=J+1
            ENDIF            
         ENDDO
         IC=IC+(NEL+1-I)
         IGETEXCITLEVEL_=IC 
         RETURN
      END
      
      SUBROUTINE  CALCRHOIIMAX(BETA,I_P,NEL,
     &           NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,NTAY,
     &           RHOMAX,ECORE)
         INTEGER IDET(NEL)
         CALL CALCRHOIIMAXR(BETA,I_P,NEL,
     &           NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,NTAY,
     &           RHOMAX,IDET,1,0,ECORE)

         RETURN
      END

      RECURSIVE SUBROUTINE CALCRHOIIMAXR(BETA,I_P,NEL,
     &           NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,NTAY,
     &           RHOMAX,IDET,NTHISEL,ILASTBASIS,ECORE)
         INTEGER IDET(NEL),NTHISEL,I,ILASTBASIS
         REAL*8 RH,RHOMAX
         DO I=ILASTBASIS+1,NBASIS
            IDET(NTHISEL)=I
            IF(NTHISEL.EQ.NEL) THEN
               CALL CALCRHO2(IDET,IDET,BETA,I_P,NEL,NBASISMAX,
     &            G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,RH,NTAY,0,ECORE)
               IF(RH.GT.RHOMAX) RHOMAX=RH
            ELSE
               CALL CALCRHOIIMAXR(BETA,I_P,NEL,
     &           NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,NTAY,
     &           RHOMAX,IDET,NTHISEL+1,I,ECORE)
            ENDIF
         ENDDO
         RETURN
      END
 
      SUBROUTINE  CALCRHOIIMAX2(BETA,I_P,NEL,
     &           NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,NTAY,
     &           ISTEPS,ECORE,RHOMAX)
         INTEGER NI(NEL)
         REAL*8 RHOMAX
         INTEGER BRR(NBASIS)
         RHOMAX=0.D0
         DO I=1,NEL
            NI(I)=BRR(I)
         ENDDO
         DO I=1,ISTEPS
            CALL GENNEWDET2(NI,NI,NEL,NBASIS,NBASISMAX,G1,ISEED)
            CALL CALCRHO2(NI,NI,BETA,I_P,NEL,NBASISMAX,
     &            G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,RH,NTAY,0,ECORE)
            IF(RH.GT.RHOMAX) RHOMAX=RH
         ENDDO
         RETURN
      END

C.. CALCRHO(I,J,BETA,P,HAMIL,LAB,NROW,NDET,RH)
C..
C.. Calculate RH = exp(-(BETA/P)H) matrix element <I|RH|J> (= RH_IJ)
C.. using Taylor expansion.
C.. Let BETA/P = B
C.. <I|exp(-BH)|J> = delta_IJ -B <I|H|J> + (B^2)/2 <I|H^2|J> + ...
C.. <I|H^2|J> = Sum_K <I|H|K><K|H|J>
C.. <I|H^3|J> = Sum_K <I|H|K><K|H^2|J>
C.. I, J = the indices of the determinants for which to get the matrix element
C.. BETA = 1/kT , P = a parameter for splitting the calc of exp(BETA H) 
C.. HAMIL = HAMIL matrix
C.. LAB = column indices of elements of HAMIL matrix
C.. NROW = # elements in each row of H (so SUM(i=1 to N-1) (NROW(i)) = index of
C..         row N in HAMIL
C.. NDET = total number of determinants K will range over
C.. RH = output RH_IJ

      SUBROUTINE CALCRHO(I,J,BETA,P,HAMIL,LAB,NROW,NDET,NTAY,RH)
      REAL*8 RH
      REAL*8 HAMIL(*)
      INTEGER LAB(*)
      INTEGER NROW(NDET)
      INTEGER IRI,IRJ,INI,INJ,P
      CALL TISET('CALCRHO   ',ISUB)
      B=BETA/P
      RH=0
C.. Deal with the delta_IJ
      IF(I.EQ.J) RH=RH+1
C.. Now do the first order term
      IF(NTAY.GE.1)
     &   RH=RH-B*GETHELEMENT(I,J,HAMIL,LAB,NROW,NDET)
C.. Now the 2nd order term
      IF(NTAY.GE.2)
     &   RH=RH+B*B*RHO2ORDER(I,J,HAMIL,LAB,NROW,NDET)/2
C.. And now the 3rd order term
      IF(NTAY.GE.3)
     &   RH=RH-B*B*B*RHO3ORDER(I,J,HAMIL,LAB,NROW,NDET)/6
      CALL TIHALT('CALCRHO   ',ISUB)
      RETURN
      END

C.. RHO2ORDER(I,J,HAMIL,LAB,NROW,NDET
C.. Calculate the 2nd order term in rho

      REAL*8 FUNCTION RHO2ORDER(I,J,HAMIL,LAB,NROW,NDET)
      REAL*8 HAMIL(*)
      INTEGER LAB(*)
      INTEGER NROW(NDET)
      INTEGER IRI,IRJ,INI,INJ

C.. Now the second order term
      IRI=IGETHROW(I,HAMIL,LAB,NROW,NDET)
      IRJ=IGETHROW(J,HAMIL,LAB,NROW,NDET)
      INI=NROW(I)
      INJ=NROW(J)
C.. If there are more elements in I than J, we choose J to iterate over, by 
C.. swapping I and J
      IF(INI.GT.INJ) THEN
         ITEMP=IRI
         IRI=IRJ
         IRJ=ITEMP
         ITEMP=INI
         INI=INJ
         INJ=ITEMP
      ENDIF
      INI=INI+IRI-1
      INJ=INJ+IRJ-1
      SUM=0.D0
C.. Now iterate over K, going along row I
      DO WHILE ((IRI.LE.INI).AND.(IRJ.LE.INJ))
         ILAB=LAB(IRI)
         DO WHILE ((LAB(IRJ).LT.ILAB).AND.(IRJ.LE.INJ))
            IRJ=IRJ+1
         ENDDO
         IF((ILAB.EQ.LAB(IRJ)).AND.(IRJ.LE.INJ)) 
     &		SUM=SUM+HAMIL(IRI)*HAMIL(IRJ)
         IRI=IRI+1
      ENDDO
      RHO2ORDER=SUM
      RETURN
      END


C.. RHO3ORDER(I,J,HAMIL,LAB,NROW,NDET
C.. Calculate the 3rd order term in rho

      REAL*8 FUNCTION RHO3ORDER(I,J,HAMIL,LAB,NROW,NDET)
      REAL*8 HAMIL(*)
      INTEGER LAB(*)
      INTEGER NROW(NDET)
      INTEGER IRI,IRJ,INI,INJ

C.. Now the second order term
      IRI=IGETHROW(I,HAMIL,LAB,NROW,NDET)
      IRJ=IGETHROW(J,HAMIL,LAB,NROW,NDET)
      INI=NROW(I)
      INJ=NROW(J)
      JNEW=J
C.. If there are more elements in I than J, we choose J to iterate over, by
C.. swapping I and J - Not sure whether this is the best for 3rd order
      IF(INI.GT.INJ) THEN
         IRI=IRJ
         INI=INJ
         JNEW=I
      ENDIF
      INI=INI+IRI-1
      SUM=0.D0
C.. Now iterate over K, going along row I
      DO IK=IRI,INI
         SUM=SUM+HAMIL(IK)*RHO2ORDER(LAB(IK),JNEW,HAMIL,LAB,NROW,NDET)
      ENDDO
      RHO3ORDER=SUM
      RETURN
      END

C.. GETHELEMENT
C.. Get matrix element of the hamiltonian
      REAL*8 FUNCTION GETHELEMENT(II,IJ,HAMIL,LAB,NROW,NDET)
      REAL*8 HAMIL(*)
      INTEGER LAB(*)
      INTEGER NROW(NDET)
C.. We only have half of H, so if J<I, return the symmetrical (J,I) element
C.. Or if we have the whole H, it's quicker to look closer to its beginning
      IF(IJ.LT.II) THEN
         I=IJ
         J=II
      ELSE
         I=II
         J=IJ
      ENDIF
      GETHELEMENT=0.D0
      INDXROW=1
C.. Find the Ith row
      DO K=1,I-1
         INDXROW=INDXROW+NROW(K)
      ENDDO
      IMAX=INDXROW+NROW(I)-1
      DO K=INDXROW,IMAX
         IF(LAB(K).GT.J) RETURN
         IF(LAB(K).EQ.J) THEN
            GETHELEMENT=HAMIL(K)
            RETURN
         ENDIF
      ENDDO
      RETURN
      END

C.. GETHROW
C.. Get the index within HAMIL of a row II of the H matrix
      INTEGER FUNCTION IGETHROW(II,HAMIL,LAB,NROW,NDET)
      REAL*8 HAMIL(*)
      INTEGER LAB(*)
      INTEGER NROW(NDET)
      IGETHROW=1
C.. Find the Ith row
      DO K=1,II-1
         IGETHROW=IGETHROW+NROW(K)
      ENDDO
      RETURN
      END

C.. For the given (I,J) pairs extracted from  RHOIJLAB,RHOIJROWS
C.. calculate approx RHOIJ, and store
      SUBROUTINE CALCSTOREAPRHOIJ(NDET,BETA,I_P,HAMIL,LAB,NROW,
     &      RHOIJROWS,RHOIJ,RHOIJLAB,TTROT,NTAY)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RH
      REAL*8 RHOIJ(*)
      INTEGER RHOIJLAB(*)
      INTEGER RHOIJROWS(NDET+1)
      LOGICAL TTROT
      CALL TISET('CSARHOIJ  ',ISUB)

      NMAX=RHOIJROWS(NDET+1)-1
      I=1
      DO N=1,NMAX
         IF(N.GE.RHOIJROWS(I+1)) I=I+1
         IF(TTROT) THEN
            CALL CALCRHOTROT(I,RHOIJLAB(N),BETA,I_P,HAMIL,LAB,NROW,
     &               NDET,NTAY,RH)
         ELSE
         
            CALL CALCRHO(I,RHOIJLAB(N),BETA,I_P,HAMIL,LAB,NROW,
     &               NDET,NTAY,RH)
         ENDIF
         RHOIJ(N)=RH
      ENDDO
      CALL TIHALT('CSARHOIJ  ',ISUB)
      RETURN
      END

C.. For the given (I,J) pairs extracted from  RHOIJLAB,RHOIJROWS
C.. calculate exact RHOIJ, and store
      SUBROUTINE CALCSTOREEXRHOIJ(NDET,BETA,I_P,NEVAL,CK,W,ECORE
     &      RHOIJROWS,RHOIJ,RHOIJLAB)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RH
      REAL*8 RHOIJ(*)
      INTEGER RHOIJLAB(*)
      INTEGER RHOIJROWS(NDET+1)
      DIMENSION CK(NDET,NEVAL)
      DIMENSION W(NEVAL)
      CALL TISET('CSERHOIJ  ',ISUB)

      NMAX=RHOIJROWS(NDET+1)-1
      I=1
      DO N=1,NMAX
         IF(N.GE.RHOIJROWS(I+1)) I=I+1
         J=RHOIJLAB(N) 
         RH=0.D0
         DO IK=1,NEVAL
            RH=RH+CK(I,IK)*CK(J,IK)*EXP(-(W(IK))*BETA/I_P)
         ENDDO
         RHOIJ(N)=RH
      ENDDO
      CALL TIHALT('CSERHOIJ  ',ISUB)
      RETURN
      END


C.. Calculate all RHOIJ using diag energies, coeffs
C.. RHOEPSILON is now calculated here (21/11/03) - previously it was 
C.. just specified.  Now RHOEPSILON=RHOEPS*RHIIMAX
      SUBROUTINE CALCRHOIJ(NDET,NEVAL,CK,W,BETA,I_P,RHOMIN,RHOMAX,
     &            NRHOS,RHOEPS,ECORE)
         IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION CK(NDET,NEVAL)
         DIMENSION W(NEVAL)
                                                                       
C.. Calculate rho_IJ for I<=J, and put in bins
         CALL TISET('RHOS      ',ISUB2)
         RHOMAX=tiny(1.d0)
         RHOMIN=huge(1.d0)
         OPEN(10,FILE='RHOS',STATUS='UNKNOWN')
         RHOEPSILON=RHOEPS
         WRITE(6,*) "RHOEPSILON: ", RHOEPSILON
         
         DO I=1,NDET
            DO J=1,NDET
               RH=0.D0
               DO IK=1,NEVAL
                  RH=RH+CK(I,IK)*CK(J,IK)*EXP(-(W(IK))*BETA/I_P)
               ENDDO
               IF(RH>RHOMAX) RHOMAX=RH
               IF(RH<RHOMIN) RHOMIN=RH
               IF(ABS(RH).GE.RHOEPSILON) THEN
                  WRITE(10,'(2I12,E25.16)')  I,J,RH
                  NRHOS=NRHOS+1
               ENDIF
            ENDDO
         ENDDO
         CLOSE(10)
         OPEN(10,FILE="NRHOS",STATUS="UNKNOWN")
         WRITE(10,"(I10)") NRHOS
         CLOSE(10)
                                                                      
         WRITE(6,*) "Number of RHOs calculated:",NDET*NDET
         WRITE(6,*) "Number of non-zero RHO:",NRHOS
         CALL TIHALT('RHOS      ',ISUB2)
         RETURN
      END
                                                                      
      SUBROUTINE CALCAPPROXRHOIJ(BETA,I_P,HAMIL,LAB,NROW,NDET,
     &              RHOMIN,RHOMAX,NRHOS,RHOEPS,TTROT,NTAY)
                                                        
         LOGICAL TTROT          
C.. Calculate rho_IJ for I<=J, and put in bins
         CALL TISET('RHOS      ',ISUB2)
                                                                     
C.. First work out RHOEPSILON:
         RHOEPSILON=RHOEPS
         WRITE(6,*) "RHOEPSILON: ", RHOEPSILON
         
         OPEN(10,FILE='RHOS',STATUS='UNKNOWN')
         NRHOS=0
         DO I=1,NDET
            DO J=1,NDET
               IF(TTROT) THEN
                  CALL CALCRHOTROT(I,J,BETA,I_P,HAMIL,LAB,NROW,NDET,
     &                  NTAY,RH)
               ELSE
                  CALL CALCRHO(I,J,BETA,I_P,HAMIL,LAB,NROW,NDET,NTAY,RH)
               ENDIF
               IF(RH>RHOMAX) RHOMAX=RH
               IF(RH<RHOMIN) RHOMIN=RH
               IF(ABS(RH).GE.RHOEPSILON) THEN
                  WRITE(10,'(2I12,E25.16)')  I,J,RH
                  NRHOS=NRHOS+1
               ENDIF
            ENDDO
         ENDDO
         CLOSE(10)
         OPEN(10,FILE="NRHOS",STATUS="UNKNOWN")
         WRITE(10,"(I10)") NRHOS
         CLOSE(10)
         WRITE(6,*) "Number of RHOs calculated:",NDET*NDET
         WRITE(6,*) "Number of non-zero RHO:",NRHOS
         CALL TIHALT('RHOS      ',ISUB2)
         RETURN
      END
                                                                      
      SUBROUTINE STORERHOIJ(NDET,NRHOS,RHOIJROWS,RHOIJ,
     &                  RHOIJLAB,NRHOMAXWIDTH)
         IMPLICIT REAL*8 (A-H,O-Z)
         INTEGER NDET,NRHO
         INTEGER RHOIJLAB(NRHOS),RHOIJROWS(NDET+1)
         REAL*8 RHOIJ(NRHOS)
                                                                        
         WRITE(6,*) "Storing Rho"
         OPEN(10,FILE='RHOS',STATUS="OLD")
         IROW=1
         NRHOMAXWIDTH=0
         RHOIJROWS(1)=1
         DO N=1,NRHOS
            READ(10,"(I12,I12,E25.16)") I,J,RH
C            WRITE(6,*) I,J,RH
C.. see if we need to move to a new row.
            IF(IROW.LT.I) THEN
               DO IROW=IROW+1,I
                  RHOIJROWS(IROW)=N
                 IF((RHOIJROWS(IROW)-RHOIJROWS(IROW-1)).GT.NRHOMAXWIDTH)
     &                  NRHOMAXWIDTH=1+RHOIJROWS(IROW)-RHOIJROWS(IROW-1)
               ENDDO
               IROW=I
            ENDIF
            RHOIJLAB(N)=J
            RHOIJ(N)=RH
         ENDDO
         DO IROW=IROW+1,NDET+1
            RHOIJROWS(IROW)=N
         ENDDO
         CLOSE(10)
         WRITE(6,*) "Max width of rho:",NRHOMAXWIDTH
         RETURN
      END
                                                                        

      SUBROUTINE CALCRHOTROT(I,J,BETA,P,HAMIL,LAB,NROW,NDET,NTAY,RH)
      REAL*8 UEXP
      REAL*8 HAMIL(*)
      INTEGER LAB(*)
      INTEGER NROW(NDET)
      INTEGER IRI,IRJ,INI,INJ,P
      CALL TISET('CALCRHOTRO',ISUB)

      B=BETA/P
      UEXP=0
C  We decompose H=T'+U'  T' being diagonal, and U' non-diag
C.. Deal with the delta_IJ
      IF(I.EQ.J) THEN
         UEXP=UEXP+1
      ELSE
C.. Now do the first order term, which only exists for non-diag
         IF(NTAY.GE.1)
     &      UEXP=UEXP-B*GETHELEMENT(I,J,HAMIL,LAB,NROW,NDET)
      ENDIF
C.. Now the 2nd order term
      IF(NTAY.GE.2)
     &   UEXP=UEXP+B*B*RHO2ORDERND(I,J,HAMIL,LAB,NROW,NDET)/2
C.. And now the 3rd order term
      IF(NTAY.GE.3)
     &   UEXP=UEXP-B*B*B*RHO3ORDERND(I,J,HAMIL,LAB,NROW,NDET)/6
      RH=EXP(-B*( GETHELEMENT(I,I,HAMIL,LAB,NROW,NDET)+
     &            GETHELEMENT(J,J,HAMIL,LAB,NROW,NDET))/2)
     &   *UEXP
      CALL TIHALT('CALCRHOTRO',ISUB)
      RETURN
      END


C.. RHO2ORDERND(I,J,HAMIL,LAB,NROW,NDET
C.. Calculate the 2nd order term in rho for non diag elements

      REAL*8 FUNCTION RHO2ORDERND(I,J,HAMIL,LAB,NROW,NDET)
      REAL*8 HAMIL(*)
      INTEGER LAB(*)
      INTEGER NROW(NDET)
      INTEGER IRI,IRJ,INI,INJ

C.. Now the second order term
      IRI=IGETHROW(I,HAMIL,LAB,NROW,NDET)
      IRJ=IGETHROW(J,HAMIL,LAB,NROW,NDET)
      INI=NROW(I)
      INJ=NROW(J)
C.. If there are more elements in I than J, we choose J to iterate over, by 
C.. swapping I and J
      IF(INI.GT.INJ) THEN
         ITEMP=IRI
         IRI=IRJ
         IRJ=ITEMP
         ITEMP=INI
         INI=INJ
         INJ=ITEMP
      ENDIF
      INI=INI+IRI-1
      INJ=INJ+IRJ-1
      SUM=0.D0
C.. Now iterate over K, going along row I
      DO WHILE ((IRI.LE.INI).AND.(IRJ.LE.INJ))
         ILAB=LAB(IRI)
         DO WHILE ((LAB(IRJ).LT.ILAB).AND.(IRJ.LE.INJ))
            IRJ=IRJ+1
         ENDDO
C.. we must reject diag elements
         IF((ILAB.EQ.LAB(IRJ)).AND.(IRJ.LE.INJ).AND.
     &      (ILAB.NE.I).AND.(LAB(IRJ).NE.J)) 
     &		SUM=SUM+HAMIL(IRI)*HAMIL(IRJ)
         IRI=IRI+1
      ENDDO
      RHO2ORDERND=SUM
      RETURN
      END


C.. RHO3ORDER(I,J,HAMIL,LAB,NROW,NDET
C.. Calculate the 3rd order term in rho for non-dial elements

      REAL*8 FUNCTION RHO3ORDERND(I,J,HAMIL,LAB,NROW,NDET)
      REAL*8 HAMIL(*)
      INTEGER LAB(*)
      INTEGER NROW(NDET)
      INTEGER IRI,IRJ,INI,INJ

C.. Now the second order term
      IRI=IGETHROW(I,HAMIL,LAB,NROW,NDET)
      IRJ=IGETHROW(J,HAMIL,LAB,NROW,NDET)
      INI=NROW(I)
      INJ=NROW(J)
      JNEW=J
C.. If there are more elements in I than J, we choose J to iterate over, by
C.. swapping I and J - Not sure whether this is the best for 3rd order
      IF(INI.GT.INJ) THEN
         IRI=IRJ
         INI=INJ
         JNEW=I
      ENDIF
      INI=INI+IRI-1
      SUM=0.D0
C.. Now iterate over K, going along row I
      DO IK=IRI,INI
C.. reject diag terms
       IF(LAB(IK).NE.I)
     &   SUM=SUM+HAMIL(IK)*RHO2ORDERND(LAB(IK),JNEW,HAMIL,LAB,NROW,NDET)
      ENDDO
      RHO3ORDERND=SUM
      RETURN
      END


C..   FCK->RHOIJ
C..   ZIA->RHOIJLAB
C..   UMAT->RHOIJROWS
C..   NMSH->NMRKS
C..   BETA->NHG
      SUBROUTINE GETRHOEXND(NI,NJ,NEL,NHG,NMRKS,RHOIJ,RHOIJLAB,
     &               RHOIJROWS,RH)
         INTEGER II,IJ,NI(NEL),NJ(NEL)
         REAL*8 RH
         REAL*8 GETRHOELEMENT
         CALL FINDDET(NI,NEL,NHG,NMRKS,II)
         CALL FINDDET(NJ,NEL,NHG,NMRKS,IJ)
         IF(II.EQ.0.OR.IJ.EQ.0) THEN
            WRITE(0,*) 'DETERMINANT NOT FOUND in GETRHOEXND ('
     &         ,(NI(K),K=1,NEL),'),(',(NJ(K),K=1,NEL),')'
            CALL FINDDET(NI,NEL,NHG,NMRKS,II)
            CALL FINDDET(NJ,NEL,NHG,NMRKS,IJ)
            STOP 'DETERMINANT NOT FOUND in GETRHOEXND'
         ENDIF
         RH=GETRHOELEMENT(II,IJ,RHOIJ,RHOIJLAB,RHOIJROWS,NDET)
         
      END

      SUBROUTINE FINDDET(NI,NEL,NDET,NMRKS,II)
         INTEGER NEL
         INTEGER NI(NEL),NDET,NMRKS(NEL,NDET)
         INTEGER II
         LOGICAL T1
         DO I=1,NDET
            T1=.TRUE.
            DO J=1,NEL
               IF(NMRKS(J,I).NE.NI(J)) T1=.FALSE.
            ENDDO
            IF(T1) THEN
               II=I
               RETURN
            ENDIF
         ENDDO
         II=0
         RETURN
      END
      SUBROUTINE CALCRHOEXND(NI,NJ,NEL,BETA,NEVAL,W,CK,NDET,NMRKS,
     &            NHG,I_P,ECORE,RH)
         INTEGER II,IJ,NI(NEL),NJ(NEL),IK,NEVAL,NDET,I_P,NBASIS
         REAL*8 RH,W(NEVAL),CK(NDET,NEVAL),BETA,ECORE
         CALL FINDDET(NI,NEL,NDET,NMRKS,II)
         CALL FINDDET(NJ,NEL,NDET,NMRKS,IJ)
         IF(II.EQ.0.OR.IJ.EQ.0) THEN
            WRITE(0,*) 'DETERMINANT NOT FOUND in CALCRHOEXND ('
     &         ,(NI(K),K=1,NEL),'),(',(NJ(K),K=1,NEL),')'
            CALL FINDDET(NI,NEL,NDET,NMRKS,II)
            CALL FINDDET(NJ,NEL,NDET,NMRKS,IJ)
            STOP 'DETERMINANT NOT FOUND in CALCRHOEXND'
         ENDIF
         RH=0.D0
         DO IK=1,NEVAL
            RH=RH+CK(II,IK)*CK(IJ,IK)*EXP(-(W(IK))*BETA/I_P)
         ENDDO
         RETURN
      END  
