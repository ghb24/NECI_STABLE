!.. Calculate RHO_IJ using a Trotter decomposition without having a stored H
!.. Calculate RHO = exp(-(BETA/P)H) matrix element <I|RHO|J> (= RHO_IJ)
!.. Trotter = RHO ~ exp(-(BETA/P)H'/2)exp(-(BETA/P)U')exp(-(BETA/P)H'/2)
!.. where H' is the diag part of H, and U' is the non-diag
!.. 
!.. NTAY is the order of the Taylor expansion for U'
!.. IC is the number of basis fns by which NI and NJ differ (or -1 if not known)
!.. 
SUBROUTINE CALCRHO2(NI,NJ,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,&
                     NMAX,ALAT,UMAT,RH,NTAY,IC2,ECORE)
      USE HElement
      TYPE(HElement) UMat(*),TMat(*),RH
      INTEGER I_P,I_HMAX,NTAY,NTRUNC,NEL,NBASIS,NBASISMAX(5,2)
      INTEGER NI(NEL),NJ(NEL),NMAX(3),IC,IC2
      REAL*8 BETA
      LOGICAL LSAME      
      INTEGER NMSH,ISUB,I,BRR(NBASIS),J
      INCLUDE 'basis.inc'
      TYPE(BasisFN) G1(*)
      COMPLEX*16 FCK(*)
      REAL*8 ALAT(3)  
      TYPE(HElement) hE,hE2,UExp,B
      IF(NTAY.LT.0) THEN
!.. We've actually hidden a matrix of rhos in the coeffs for calcing RHO
         CALL GETRHOEXND(NI,NJ,NEL,BETA,NMSH,FCK,ZIA,UMAT,RH)
         RETURN
      ELSEIF(NTAY.EQ.0) THEN
!.. NTAY=0 signifying we're going to calculate the RHO values when we
!.. need them from the list of eigenvalues.
!.. Hide NMSH=NEVAL
!..      FCK=W
!..      ZIA=CK
!..      UMAT=NDET
!..      ALAT=NMRKS
          CALL CALCRHOEXND(NI,NJ,NEL,BETA,NMSH,FCK,TMat,UMAT,ALAT,NBASIS,I_P,ECORE,RH)
         RETURN
      ENDIF
      CALL TISET('CALCRHO2  ',ISUB)
      IC=IC2
      B=BETA/I_P
      UEXP=0.D0
      IF(IC.LT.0)   IC=IGETEXCITLEVEL(NI,NJ,NEL)
      IF(IC.EQ.0) THEN
         LSAME=.TRUE.
      ELSE
         LSAME=.FALSE.
      ENDIF
          

      IF(LSAME) THEN
         UExp=UExp+HElement(1.D0)
      ELSE
!.. Now do the first order term, which only exists for non-diag
         IF(NTAY.GE.1) THEN
            UExp=UExp-B*GETHELEMENT2(NI,NJ,NEL,NBASISMAX,G1,&
     &      NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,IC,ECORE)
         ENDIF
      ENDIF
!.. Now the 2nd order term
      IF(NTAY.GE.2) UExp=UExp+B*B*HElement(RHO2ORDERND2(NI,NJ,NEL,NBASISMAX,    &
     &            G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,IC,ECORE)/2.D0)

      hE =(GETHELEMENT2(NI,NI,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,0,ECORE) &
          +GETHELEMENT2(NJ,NJ,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,0,ECORE))&
            /HElement(2.D0)
      RH=EXP(HElement(-BETA/I_P)*hE)*UExp
!WRDET
!      WRITE(6,"(A,$)") "RHO:"
!      CALL WRITEDET(6,NI,NEL,.FALSE.)
!      CALL WRITEDET(6,NJ,NEL,.FALSE.)
!      WRITE(6,*) RH
      CALL TIHALT('CALCRHO2  ',ISUB)
      RETURN
     END
!.. RHO2ORDERND2(I,J,HAMIL,LAB,NROW,NDET
!.. Calculate the 2nd order term in rho for non diag elements
!.. IC2 is the number of basis fns that differ in NI and NJ (or -1 if not known)

     FUNCTION Rho2OrderND2(NI,NJ,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMat,IC2,ECORE)
!.. We use a crude method and generate all possible 0th, 1st, and 2nd
!.. excitations of I and of J.  The intersection of these lists is the
!.. selection of dets we want.
         USE HElement
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(*)
         TYPE(HElement) Rho2OrderND2
         COMPLEX*16 FCK(*)
         REAL*8 ALAT(*),ECORE
         TYPE(HElement) UMat(*),TMat(*)
         INTEGER NEL,NBASIS,NBASISMAX(5,2),BRR(*)
         INTEGER NI(NEL),NJ(NEL),IC2,NMSH,NMAX
         INTEGER LSTI(NEL,NBASIS*NBASIS*NEL*NEL)
         INTEGER LSTJ(NEL,NBASIS*NBASIS*NEL*NEL)
         INTEGER NLISTI, NLISTJ, IC,I,J,ICC
         INTEGER ICI(NBASIS*NBASIS*NEL*NEL)
         INTEGER ICJ(NBASIS*NBASIS*NEL*NEL),NLISTMAX
         INTEGER CMP,IGETEXCITLEVEL,ICMPDETS
         TYPE(HElement) SUM
         SUM=0.D0
         NLISTMAX=NBASIS*NBASIS*NEL*NEL
         IC=IC2
         IF(IC.LT.0) IC=IGETEXCITLEVEL(NI,NJ,NEL)

!.. the 1 at the ends ensures K.NE. I or J, as this would make
!.. <I|U'|K> or <K|U'|J> zero (U' has a zero diag)
         CALL GENEXCIT(NI,2,NBASIS,NEL,LSTI,ICI,NLISTI,1,G1,.TRUE.,     &
     &         NBASISMAX,.FALSE.)                                     
         CALL GENEXCIT(NJ,2,NBASIS,NEL,LSTJ,ICJ,NLISTJ,1,G1,.TRUE.,     &
     &         NBASISMAX,.FALSE.)
         I=1
         J=1
!.. Now iterate over K, going along row I
         DO WHILE ((I.LE.NLISTI).AND.(J.LE.NLISTJ))
            CMP=ICMPDETS(LSTI(1,I),LSTJ(1,J),NEL)
!.. While I>J, we increase J
            DO WHILE ((CMP.GT.0).AND.(J.LT.NLISTJ))
               J=J+1
               CMP=ICMPDETS(LSTI(1,I),LSTJ(1,J),NEL)
            ENDDO
            IF(CMP.EQ.0) THEN 
               SUM=SUM+  GETHELEMENT2(NI,LSTI(1,I),NEL,NBASISMAX,       &
     &      G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMat,ICI(I),ECORE)     &
     &               *GETHELEMENT2(LSTJ(1,J),NJ,NEL,NBASISMAX,          &
     &      G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMat,ICJ(J),ECORE)
            ENDIF
            I=I+1
         
         ENDDO
         RHO2ORDERND2=SUM
         RETURN
      END
