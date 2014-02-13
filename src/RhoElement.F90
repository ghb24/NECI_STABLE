! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
!.. Calculate RHO_IJ using a Trotter decomposition without having a stored H
!.. Calculate RHO = exp(-(BETA/P)H) matrix element <I|RHO|J> (= RHO_IJ)
!.. Trotter = RHO ~ exp(-(BETA/P)H'/2)exp(-(BETA/P)U')exp(-(BETA/P)H'/2)
!.. where H' is the diag part of H, and U' is the non-diag
!.. 
!.. NTAY is the order of the Taylor expansion for U'
!.. IC is the number of basis fns by which NI and NJ differ (or -1 if not known)
!.. 
SUBROUTINE CALCRHO2(NI,NJ,BETA,I_P,NEL,G1,NBASIS,NMSH,FCK,&
                     NMAX,ALAT,UMAT,RH,NTAY,IC2,ECORE)
      Use Determinants, only: get_helement, nUHFDet, &
                              E0HFDet
      use constants, only: dp
      use SystemData , only : TSTOREASEXCITATIONS,BasisFN
      use global_utilities
      IMPLICIT NONE
      HElement_t UMat(*),RH
      INTEGER I_P,NTAY(2),NEL,NBASIS
      INTEGER NI(NEL),NJ(NEL),NMAX,IC,IC2
      real(dp) BETA,ECORE
      LOGICAL tSameD      
      INTEGER NMSH,IGETEXCITLEVEL
      type(timer), save :: proc_timer
      TYPE(BasisFN) G1(*)
      complex(dp) FCK(*)
      real(dp) ALAT(3)  
      HElement_t hE,UExp,B,EDIAG
      IF(NTAY(1).LT.0) THEN
!.. We've actually hidden a matrix of rhos in the coeffs for calcing RHOa
          STOP "GETRHOEXND has been removed."
!         CALL GETRHOEXND(NI,NJ,NEL,BETA,NMSH,FCK,UMAT,RH)
         RETURN
      ELSEIF(NTAY(1).EQ.0) THEN
!.. NTAY=0 signifying we're going to calculate the RHO values when we
!.. need them from the list of eigenvalues.
!.. Hide NMSH=NEVAL
!..      FCK=W
!..      ZIA=CK
!..      UMAT=NDET
!..      ALAT=NMRKS
          STOP "Exact RHO calculation broken."
!          CALL CALCRHOEXND(NI,NJ,NEL,BETA,NMSH,FCK,UMAT,ALAT,I_P,RH)
         RETURN
      ENDIF
      proc_timer%timer_name='CALCRHO2  '
      call set_timer(proc_timer,60)
      IC=IC2
      B=BETA/I_P
      UEXP=0.0_dp
      IF(IC.LT.0)   IC=IGETEXCITLEVEL(NI,NJ,NEL)
      IF(IC.EQ.0) THEN
         tSAMED=.TRUE.
      ELSE
         tSAMED=.FALSE.
      ENDIF
          
      IF(tStoreAsExcitations.AND.nI(1).eq.-1.AND.nJ(1).eq.-1) THEN
!Store as excitations.
         IF(NTAY(2).NE.3) STOP "Store as Excitations only works for Fock-Partition-Lowdiag"
!Partition with Trotter with H(0) having just the Fock Operators
!Fock-Partition-Lowdiag
         IF(tSAMED) THEN
            call GetH0Element(nJ,nEl,nMax,nBasis,ECore,EDiag)
            EDiag=EDiag+(E0HFDET)
            RH=EXP(-B*EDiag)
         ELSE
            call GetH0Element(nI,nEl,nMax,nBasis,ECore,UExp)
            UExp=UExp+(E0HFDET)
            call GetH0Element(nJ,nEl,nMax,nBasis,ECore,EDiag)
            EDiag=(UExp+UExp+EDiag)/(2.0_dp)
            UExp = get_helement (nI, nJ, IC)
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
         IF(tSAMED) THEN
            UExp=UExp+(1.0_dp)
         ELSE
!.. Now do the first order term, which only exists for non-diag
            UExp=UExp-B*get_helement(nI, nJ, IC)
         ENDIF
!.. Now the 2nd order term
!      IF(NTAY.GE.2) UExp=UExp+B*B*(RHO2ORDERND2(NI,NJ,NEL,NBASISMAX,    &
!     &            G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,IC,ECORE)/2.0_dp)

         hE = (get_helement(nI, nI, 0) + &
               get_helement (nJ, nJ, 0)) / (2.0_dp)
         RH=EXP((-BETA/I_P)*hE)*UExp
      ELSEIF(NTAY(2).EQ.2) THEN
!Partition with Trotter with H(0) having just the Fock Operators
         IF(tSAMED) THEN
            call GetH0Element(nI,nEl,nMax,nBasis,ECore,EDiag)
            UExp=1.0_dp
!Fock-Partition
            UExp=UExp-B*(get_helement(nI, nI, 0)-EDiag)
            RH=EXP(-B*EDiag)*UExp
         ELSE
            call GetH0Element(nI,nEl,nMax,nBasis,ECore,UExp)
            call GetH0Element(nJ,nEl,nMax,nBasis,ECore,EDiag)
            EDiag=(UExp+EDiag)/(2.0_dp)
            UExp=get_helement(nI, nJ, IC)
            UExp=-B*UExp
            RH=EXP(-B*EDiag)*UExp
         ENDIF
      ELSEIF(NTAY(2).EQ.3) THEN
!Partition with Trotter with H(0) having just the Fock Operators
!Fock-Partition-Lowdiag
         IF(tSAMED) THEN
            call GetH0Element(nI,nEl,nMax,nBasis,ECore,EDiag)
            RH=EXP(-B*EDiag)
         ELSE
            call GetH0Element(nI,nEl,nMax,nBasis,ECore,UExp)
            call GetH0Element(nJ,nEl,nMax,nBasis,ECore,EDiag)
            EDiag=(UExp+EDiag)/(2.0_dp)
            UExp=get_helement(nI, nJ, IC)
            UExp=-B*UExp
            RH=EXP(-B*EDiag)*UExp
         ENDIF
      ELSEIF(NTAY(2).EQ.4) THEN
!         Do a simple Taylor expansion on the whole lot
         IF(tSAMED) THEN
            UEXP=1.0_dp
         ELSE
            UEXP=0.0_dp
         ENDIF
         RH=UEXP-B*get_helement(nI, nJ, IC)
      ELSEIF(NTAY(2).EQ.5) THEN
!Fock-Partition-DCCorrect-LowDiag
!Partition with Trotter with H(0) having just the Fock Operators.  Taylor diagonal to zeroeth order, and off-diag to 1st.
! Instead of 
         IF(tSAMED) THEN
            call GetH0ElementDCCorr(nUHFDet,nI,nEl,G1,nBasis,NMAX,ECore,EDiag)
            RH=EXP(-B*EDiag)
         ELSE
            call GetH0ElementDCCorr(nUHFDet,nI,nEl,G1,nBasis,NMAX,ECore,UExp)
            call GetH0ElementDCCorr(nUHFDet,nJ,nEl,G1,nBasis,NMAX,ECore,EDiag)
            EDiag=(UExp+EDiag)/(2.0_dp)
            UExp=get_helement(nI, nJ, IC)
            UExp=-B*UExp
            RH=EXP(-B*EDiag)*UExp
         ENDIF
      ENDIF

!WRDET
!      WRITE(6,"(A)",advance='no') "RHO:"
!      CALL WRITEDET(6,NI,NEL,.FALSE.)
!      CALL WRITEDET(6,NJ,NEL,.FALSE.)
!      WRITE(6,*) RH
      call halt_timer(proc_timer)
      RETURN
     END
!.. RHO2ORDERND2(I,J,HAMIL,LAB,NROW,NDET
!.. Calculate the 2nd order term in rho for non diag elements
!.. IC2 is the number of basis fns that differ in NI and NJ (or -1 if not known)

     FUNCTION Rho2OrderND2(NI,NJ,NEL,NBASISMAX,G1,NBASIS,IC2)
!.. We use a crude method and generate all possible 0th, 1st, and 2nd
!.. excitations of I and of J.  The intersection of these lists is the
!.. selection of dets we want.
         Use Determinants, only: get_helement
         use constants, only: dp
         use SystemData, only: BasisFN
         IMPLICIT NONE
         TYPE(BasisFN) G1(*)
         HElement_t Rho2OrderND2
         INTEGER NEL,NBASIS,nBasisMax(5,*)
         INTEGER NI(NEL),NJ(NEL),IC2
         INTEGER LSTI(NEL,NBASIS*NBASIS*NEL*NEL)
         INTEGER LSTJ(NEL,NBASIS*NBASIS*NEL*NEL)
         INTEGER NLISTI, NLISTJ, IC,I,J
         INTEGER ICI(NBASIS*NBASIS*NEL*NEL)
         INTEGER ICJ(NBASIS*NBASIS*NEL*NEL),NLISTMAX
         INTEGER CMP,IGETEXCITLEVEL,ICMPDETS
         HElement_t SUM1
         SUM1=0.0_dp
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
               SUM1=SUM1+  get_helement (nI, lstI(:,I)) * &
                           get_helement(lstJ(:,J), nJ)
            ENDIF
            I=I+1
         
         ENDDO
         RHO2ORDERND2=SUM1
         RETURN
      END

!  Get a matrix element of the double-counting corrected unperturbed Hamiltonian.
!  This is just the sum of the Hartree-Fock eigenvalues 
!   with the double counting subtracted, Sum_i eps_i - 1/2 Sum_i,j <ij|ij>-<ij|ji>.  (i in HF det, j in excited det)
      subroutine GetH0ElementDCCorr(nHFDet,nJ,nEl,G1,nBasis,NMAX,ECore,hEl)
         use constants, only: dp
         use Integrals_neci, only: get_umat_el
         use UMatCache
         use SystemData, only: BasisFN,Arr
         implicit none
         integer nEl,nBasis
         integer nHFDet(nEl),nJ(nEl)
         type(BasisFN) G1(*)
         HElement_t hEl
         real(dp) ECore
         integer i,j,NMAX
         integer IDHF(nEl),IDJ(nEl)
         hEl=(ECore)
         do i=1,nEl
            hEl=hEl+(Arr(nJ(i),2))
            IDHF(i) = gtID(nHFDet(i))
            IDJ(i) = gtID(nJ(i))
         enddo
         do i=1,nEl
            do j=1,nEl
!Coulomb term
               hEl=hEl-(0.5_dp)*get_umat_el(IDHF(i),IDJ(j),IDHF(i),IDJ(j))
               if(G1(nHFDet(i))%Ms.eq.G1(nJ(j))%Ms) then
!Exchange term
                  hEl=hEl+(0.5_dp)*get_umat_el(IDHF(i),IDJ(j),IDJ(j),IDHF(i))
               endif
            enddo
         enddo
!         call writedet(77,nj,nel,.false.)
!         write(77,*) "H0DC",hEl
!         call neci_flush(77)
      end
