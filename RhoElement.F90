!.. Calculate RHO_IJ using a Trotter decomposition without having a stored H
!.. Calculate RHO = exp(-(BETA/P)H) matrix element <I|RHO|J> (= RHO_IJ)
!.. Trotter = RHO ~ exp(-(BETA/P)H'/2)exp(-(BETA/P)U')exp(-(BETA/P)H'/2)
!.. where H' is the diag part of H, and U' is the non-diag
!.. 
!.. NTAY is the order of the Taylor expansion for U'
!.. IC is the number of basis fns by which NI and NJ differ (or -1 if not known)
!.. 
SUBROUTINE CALCRHO2(NI,NJ,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,&
                     NMAX,ALAT,UMAT,RH,NTAY,IC2,ECORE)
      USE HElem
      USE System , only : TSTOREASEXCITATIONS
      IMPLICIT NONE
      TYPE(HElement) UMat(*),RH
      INTEGER I_P,I_HMAX,NTAY(2),NTRUNC,NEL,NBASIS,NBASISMAX(5,5)
      INTEGER NI(NEL),NJ(NEL),NMAX,IC,IC2
      REAL*8 BETA,ECORE
      LOGICAL LSAME      
      INTEGER NMSH,ISUB,I,BRR(NBASIS),J,IGETEXCITLEVEL
      INCLUDE 'basis.inc'
      INCLUDE 'uhfdet.inc'
      TYPE(BasisFN) G1(*)
      COMPLEX*16 FCK(*)
      REAL*8 ALAT(3)  
      TYPE(HElement) hE,hE2,UExp,B,EDIAG
      IF(NTAY(1).LT.0) THEN
!.. We've actually hidden a matrix of rhos in the coeffs for calcing RHO
         CALL GETRHOEXND(NI,NJ,NEL,BETA,NMSH,FCK,UMAT,RH)
         RETURN
      ELSEIF(NTAY(1).EQ.0) THEN
!.. NTAY=0 signifying we're going to calculate the RHO values when we
!.. need them from the list of eigenvalues.
!.. Hide NMSH=NEVAL
!..      FCK=W
!..      ZIA=CK
!..      UMAT=NDET
!..      ALAT=NMRKS
          CALL CALCRHOEXND(NI,NJ,NEL,BETA,NMSH,FCK,UMAT,ALAT,NBASIS,I_P,ECORE,RH)
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
          
      IF(tStoreAsExcitations.AND.nI(1).eq.-1.AND.nJ(1).eq.-1) THEN
!Store as excitations.
         IF(NTAY(2).NE.3) STOP "Store as Excitations only works for Fock-Partition-Lowdiag"
!Partition with Trotter with H(0) having just the Fock Operators
!Fock-Partition-Lowdiag
         IF(LSAME) THEN
            call GetH0Element(nJ,nEl,nMax,nBasis,ECore,EDiag)
            EDiag=EDiag+HElement(E0HFDET)
            RH=EXP(-B*EDiag)
         ELSE
            call GetH0Element(nI,nEl,nMax,nBasis,ECore,UExp)
            UExp=UExp+HElement(E0HFDET)
            call GetH0Element(nJ,nEl,nMax,nBasis,ECore,EDiag)
            EDiag=(UExp+UExp+EDiag)/HElement(2.D0)
            UExp=GETHELEMENT2(NI,NJ,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,IC,ECORE)
            UExp=-B*UExp
            RH=EXP(-B*EDiag)*UExp
         ENDIF
         RETURN
      ELSEIF(nI(1).eq.-1.or.nJ(1).eq.-1) THEN
         STOP "Store as Excitations used, but not allowed in CALCRHO2"
      ENDIF
      IF(NTAY(2).EQ.1) THEN
!Diag-Partition
!Partition with Trotter using H(0) containing the complete diagonal
         IF(LSAME) THEN
            UExp=UExp+HElement(1.D0)
         ELSE
!.. Now do the first order term, which only exists for non-diag
            UExp=UExp-B*GETHELEMENT2(NI,NJ,NEL,NBASISMAX,G1,&
     &      NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,IC,ECORE)
         ENDIF
!.. Now the 2nd order term
!      IF(NTAY.GE.2) UExp=UExp+B*B*HElement(RHO2ORDERND2(NI,NJ,NEL,NBASISMAX,    &
!     &            G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,IC,ECORE)/2.D0)

         hE =(GETHELEMENT2(NI,NI,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,0,ECORE) &
          +GETHELEMENT2(NJ,NJ,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,0,ECORE))&
            /HElement(2.D0)
         RH=EXP(HElement(-BETA/I_P)*hE)*UExp
      ELSEIF(NTAY(2).EQ.2) THEN
!Partition with Trotter with H(0) having just the Fock Operators
         IF(LSAME) THEN
            call GetH0Element(nI,nEl,nMax,nBasis,ECore,EDiag)
            UExp=1.D0
!Fock-Partition
            UExp=UExp-B*(GETHELEMENT2(NI,NI,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,0,ECORE)-EDiag)
            RH=EXP(-B*EDiag)*UExp
         ELSE
            call GetH0Element(nI,nEl,nMax,nBasis,ECore,UExp)
            call GetH0Element(nJ,nEl,nMax,nBasis,ECore,EDiag)
            EDiag=(UExp+EDiag)/HElement(2.D0)
            UExp=GETHELEMENT2(NI,NJ,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,IC,ECORE)
            UExp=-B*UExp
            RH=EXP(-B*EDiag)*UExp
         ENDIF
      ELSEIF(NTAY(2).EQ.3) THEN
!Partition with Trotter with H(0) having just the Fock Operators
!Fock-Partition-Lowdiag
         IF(LSAME) THEN
            call GetH0Element(nI,nEl,nMax,nBasis,ECore,EDiag)
            RH=EXP(-B*EDiag)
         ELSE
            call GetH0Element(nI,nEl,nMax,nBasis,ECore,UExp)
            call GetH0Element(nJ,nEl,nMax,nBasis,ECore,EDiag)
            EDiag=(UExp+EDiag)/HElement(2.D0)
            UExp=GETHELEMENT2(NI,NJ,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,IC,ECORE)
            UExp=-B*UExp
            RH=EXP(-B*EDiag)*UExp
         ENDIF
      ELSEIF(NTAY(2).EQ.4) THEN
!         Do a simple Taylor expansion on the whole lot
         IF(LSAME) THEN
            UEXP=1.D0
         ELSE
            UEXP=0.D0
         ENDIF
         RH=UEXP-B*GETHELEMENT2(NI,NJ,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,IC,ECORE)
      ELSEIF(NTAY(2).EQ.5) THEN
!Fock-Partition-DCCorrect-LowDiag
!Partition with Trotter with H(0) having just the Fock Operators.  Taylor diagonal to zeroeth order, and off-diag to 1st.
! Instead of 
         IF(LSAME) THEN
            call GetH0ElementDCCorr(nUHFDet,nI,nEl,nBasisMax,G1,nBasis,Brr,NMSH,FCK,NMAX,ALAT,UMat,ECore,EDiag)
            RH=EXP(-B*EDiag)
         ELSE
            call GetH0ElementDCCorr(nUHFDet,nI,nEl,nBasisMax,G1,nBasis,Brr,NMSH,FCK,NMAX,ALAT,UMat,ECore,UExp)
            call GetH0ElementDCCorr(nUHFDet,nJ,nEl,nBasisMax,G1,nBasis,Brr,NMSH,FCK,NMAX,ALAT,UMat,ECore,EDiag)
            EDiag=(UExp+EDiag)/HElement(2.D0)
            UExp=GETHELEMENT2(NI,NJ,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,IC,ECORE)
            UExp=-B*UExp
            RH=EXP(-B*EDiag)*UExp
         ENDIF
      ENDIF

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

     FUNCTION Rho2OrderND2(NI,NJ,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMat,IC2,ECORE)
!.. We use a crude method and generate all possible 0th, 1st, and 2nd
!.. excitations of I and of J.  The intersection of these lists is the
!.. selection of dets we want.
         USE HElem
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(*)
         TYPE(HElement) Rho2OrderND2
         COMPLEX*16 FCK(*)
         REAL*8 ALAT(*),ECORE
         TYPE(HElement) UMat(*)
         INTEGER NEL,NBASIS,NBASISMAX(5,5),BRR(*)
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
     &      G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMat,ICI(I),ECORE)     &
     &               *GETHELEMENT2(LSTJ(1,J),NJ,NEL,NBASISMAX,          &
     &      G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMat,ICJ(J),ECORE)
            ENDIF
            I=I+1
         
         ENDDO
         RHO2ORDERND2=SUM
         RETURN
      END

!  Get a matrix element of the double-counting corrected unperturbed Hamiltonian.
!  This is just the sum of the Hartree-Fock eigenvalues 
!   with the double counting subtracted, Sum_i eps_i - 1/2 Sum_i,j <ij|ij>-<ij|ji>.  (i in HF det, j in excited det)
      subroutine GetH0ElementDCCorr(nHFDet,nJ,nEl,nBasisMax,G1,nBasis,Brr,NMSH,FCK,Arr,ALAT,UMat,ECore,hEl)
         USE HElem
         use UMatCache
         implicit none
         include 'basis.inc'
         integer nHFDet(nEl),nJ(nEl),nEl,nBasis
         type(BasisFN) G1(*)
         integer Brr(nBasis),nBasisMax(5,5)
         type(HElement) UMat(*)  
         type(HElement) hEl
         real*8 Arr(nBasis,2),ECore
         integer i,j
         INTEGER NMSH
         COMPLEX*16 FCK(*)
         REAL*8 ALAT(3)
         integer IDHF(nEl),IDJ(nEl)
         hEl=HElement(ECore)
         do i=1,nEl
            hEl=hEl+HElement(Arr(nJ(i),2))
            call gtID(nBasisMax,nHFDet(i),IDHF(i))
            call gtID(nBasisMax,nJ(i),IDJ(i))
         enddo
         do i=1,nEl
            do j=1,nEl
!Coulomb term
               hEl=hEl-HElement(0.5d0)*GetUMatEl(nBasisMax,UMat,ALat,nBasis,nBasisMax(2,3),G1,IDHF(i),IDJ(j),IDHF(i),IDJ(j))
               if(G1(nHFDet(i))%Ms.eq.G1(nJ(j))%Ms) then
!Exchange term
                  hEl=hEl+HElement(0.5d0)*GetUMatEl(nBasisMax,UMat,ALat,nBasis,nBasisMax(2,3),G1,IDHF(i),IDJ(j),IDJ(j),IDHF(i))
               endif
            enddo
         enddo
!         call writedet(77,nj,nel,.false.)
!         write(77,*) "H0DC",hEl
!         call flush(77)
      end
