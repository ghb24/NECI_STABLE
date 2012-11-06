      SUBROUTINE FODMAT(NEL,NHG,NDET,NEVAL,ISTATE,
     &    NMRKS,G1,CK,NMAXX,NMAXY,NMAXZ,SUMA)
      use constants, only: dp
      use SystemData, only: BasisFN
      use global_utilities
      use MemoryManager, only: TagIntType
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      real(dp), ALLOCATABLE :: FODM(:,:),Z(:,:),W(:)
C..Diagonalisation vectors
      real(dp), ALLOCATABLE :: WORK(:)
      real(dp), ALLOCATABLE :: IWORK(:),IFAIL(:)
      integer(TagIntType), save :: tagFODM=0,tagZ=0,tagW=0
      integer(TagIntType), save :: tagWORK=0
      integer(TagIntType), save :: tagIWORK=0,tagIFAIL=0
      integer :: NEl,NHG,NDET,NEVAL,ISTATE,NMAXX,NMAXY,NMAXZ
      REAL(dp) :: CK(NDET,NEVAL),SUM
      INTEGER :: N1E(NHG),N2E(NHG),IG,JG,IN,ID,ID1,ID2,I,J,K,INFO
      INTEGER :: M,IPTOT,IFLAG,ND2,ND1,IQ,NMRKS(NEL,*)
      REAL(dp) :: SUMA(NMAXX,NMAXY,NMAXZ),EPS,D,SUMC,ABSTOL,VL,VU
      TYPE(BASISFN) G1(NHG)
      character(*),parameter :: this_routine='FODMAT'
C..
      WRITE(6,'(1X,A)') '*** NATURAL ORBITAL EVALUATION ***'
      WRITE(6,'(2X,A)') 'FILE ./NATURAL BEING UPDATED'
      EPS=1.D-8 
C..
      allocate(FODM(NHG,NHG))
      call LogMemAlloc('FODM',NHG**2,8,this_routine,tagFODM)
      FODM=0.0_dp
C..
      N1E(1:NHG)=0
      N2E(1:NHG)=0
C..
      DO IG=1,NHG     ! n   (orbitals) 
        DO JG=IG,NHG    ! nprime (orbitals)
          IF(JG.EQ.IG) THEN
C..
            SUM=0.0_dp
            DO ID=1,NDET
              IF(ABS(CK(ID,ISTATE)).LT.EPS) GOTO 100 
              DO IN=1,NEL
                IF(NMRKS(IN,ID).EQ.IG) THEN
                  SUM=SUM+CK(ID,ISTATE)*CK(ID,ISTATE)
                  GOTO 100
                ENDIF
              ENDDO
 100          CONTINUE
            ENDDO
            FODM(JG,IG)=SUM
          ELSE
            SUM=0.0_dp
            DO ID1=1,NDET
              IF(ABS(CK(ID1,ISTATE)).LT.EPS) GOTO 300 
              DO ID2=1,NDET
                IF(ABS(CK(ID2,ISTATE)).LT.EPS) GOTO 200 
C..   Compare NDET1 to NDET2
                 IQ=0
                 K=0
                 DO I=1,NEL
                    DO J=1,NEL
                       IF(NMRKS(I,ID1).EQ.NMRKS(J,ID2)) THEN 
                          GOTO 40
                       ENDIF
                    ENDDO
                    K=K+1
 40                 CONTINUE
                 ENDDO
                 IF(K.NE.1) GOTO 200
C..   we have the single excitation from i   
                 CALL WHCH(NMRKS,NEL,ID1,ID2,ND1,ND2)
                 CALL LINEUP(NEL,NMRKS(1:NEL,ID1),NMRKS(1:NEL,ID2),
     &                NHG,N1E,N2E,IFLAG,IPTOT)
                 IF(ND1.EQ.IG.AND.ND2.EQ.JG) THEN 
                    SUM=SUM + REAL(IPTOT,dp)*
     &                   CK(ID1,ISTATE)*CK(ID2,ISTATE)
                 ENDIF
 200             CONTINUE
              ENDDO             ! id2
 300          CONTINUE
           ENDDO                ! id1 
           FODM(JG,IG)=SUM
        ENDIF  
C..   
      ENDDO                     ! jg
      ENDDO    !ig
C..   
C      OPEN(10,FILE='FODM',STATUS='UNKNOWN')
C      CALL WRITE_MATRIX('   FODM: ',NHG,NHG,FODM)
C      CLOSE(10)
C..Now diagonalise FODM
      allocate(WORK(8*NHG))
      call LogMemAlloc('WORK',8*NHG,8,this_routine,tagWORK)
      allocate(IWORK(5*NHG))
      call LogMemAlloc('IWORK',5*NHG,4,this_routine,tagIWORK)
      allocate(IFAIL(NHG))
      call LogMemAlloc('IFAIL',NHG,4,this_routine,tagIFAIL)
      allocate(Z(NHG,NHG))
      call LogMemAlloc('Z',NHG**2,8,this_routine,tagZ)
      allocate(W(NHG))
      call LogMemAlloc('W',NHG,8,this_routine,tagW)
      WORK=0.0_dp
      IWORK(1:5*NHG)=0
      IFAIL(1:NHG)=0
C..Natural orbitals and occupation numbers
      Z=0.0_dp
      W=0.0_dp
C..
      ABSTOL=1.0e-15_dp_dp
      CALL DSYEVX('V','A','L',NHG,FODM,NHG,VL,VU,1,NHG,ABSTOL,M,
     &    W,Z,NHG,WORK,8*NHG,IWORK,IFAIL,INFO)
      IF(INFO.GT.0) THEN
        WRITE(6,*) ' !!!!! WARNING FROM DIAGONALISOR !!!!!'
        WRITE(6,*) ' INFO:',INFO
        WRITE(6,*) '    M:',M
        WRITE(6,*) ' !!!!! WARNING FROM DIAGONALISOR !!!!!'
      ENDIF
C..
      OPEN(10,FILE='NATURAL',STATUS='UNKNOWN')
      DO I=1,NHG   
        WRITE(10,*) ' NATURAL ORBITAL NO. ' , I
        WRITE(10,'(F19.5)') W(I)
        WRITE(10,*) ' COEFFICIENT EXPANSION : ' 
        SUMC=0.0_dp
        DO J=1,NHG
          SUMC=SUMC+Z(J,I)*Z(J,I)
          IF(ABS(Z(J,I)).GT.1.0e-12_dp_dp) THEN
            WRITE(10,'(I7,1X,5(I7,1X),E19.6)') J,(G1(J)%K(K),K=1,5),
     &         Z(J,I)
          ENDIF
        ENDDO
        WRITE(10,'(1X,A,F19.5)') ' SQUARE OF COEFFICIENTS : ' , SUMC
      ENDDO
C..
      SUM=0.0_dp
      WRITE(6,'(1X,A,F19.5)') ' OCCUPATION NUMBERS : '
      DO I=NHG,1,-1   
        SUM=SUM+W(I)
        WRITE(6,'(I4,2F10.3)') NHG-I+1,W(I),SUM
      ENDDO
C      WRITE(6,'(1X,A,F19.5)') ' TOTAL SUM: ' , SUM
      CLOSE(10)
C.. Population analysis
      WRITE(6,*) ' POPULATION ANALYSIS OF THE NATURAL ORBITALS'
      SUMA=0.0_dp
      DO I=1,NHG
        D=0.0_dp
        DO J=1,NHG
          D=D+Z(I,J)**2*W(J)
        ENDDO
         SUMA(G1(I)%K(1),G1(I)%K(2),G1(I)%K(3))=
     &            SUMA(G1(I)%K(1),G1(I)%K(2),G1(I)%K(3))+D
      ENDDO
      DO I=1,NMAXZ  ! over nmaxx 
        DO J=1,NMAXY  ! over nmaxy
          DO K=1,NMAXX  ! over nmaxx
            IF(SUMA(K,J,I).GT.0) THEN
              WRITE(6,'(3I3,F9.2)') K,J,I,SUMA(K,J,I) 
            ENDIF
          ENDDO
        ENDDO
      ENDDO 
      WRITE(6,*) ' TOTAL ANGULAR-MOMENTUM DECOMPOSED POPULATION:'
      DO I=1,NMAXZ ! over nmaxz
        SUM=0.0_dp
        DO J=1,NMAXY ! over nmaxy
          DO K=1,NMAXX
            SUM=SUM+SUMA(K,J,I)
          ENDDO
        ENDDO
        WRITE(6,'(I3,F9.3)') I,SUM
      ENDDO       
      deallocate(FODM,Z,W,WORK,IWORK,IFAIL)
      call LogMemDealloc(this_routine,tagZ)
      call LogMemDealloc(this_routine,tagW)
      call LogMemDealloc(this_routine,tagFODM)
      call LogMemDealloc(this_routine,tagWORK)
      call LogMemDealloc(this_routine,tagIWORK)
      call LogMemDealloc(this_routine,tagIFAIL)
C..Now we calculate the overlap of the single det. composed of the 
C..natural orbitals with the many-body wavefunction
C..<Phi_k|Psi>
      RETURN
      END
C ==---------------------------------------------------------------------==
      SUBROUTINE WHCH(NMRKS,NEL,ID1,ID2,ND1,ND2)
      use constants, only: dp
!      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT NONE
      INTEGER :: ND1,ND2,ID1,ID2,NEL,IQ,I,J,K
      INTEGER :: NMRKS(NEL,*)
C..ID1 and ID2 are the labels of the determinants under investigation 
C..ND1 and ND2 are the labels of orbitals they differ by
C..Compare NDET1 to NDET2
      ND1=0
      ND2=0
      IQ=0
      DO I=1,NEL
        K=0
        DO J=1,NEL
          IF(NMRKS(I,ID1).EQ.NMRKS(J,ID2)) CONTINUE
          IF(NMRKS(I,ID1).NE.NMRKS(J,ID2)) K=K+1
        ENDDO
        IF(K.EQ.NEL) THEN
          IQ=IQ+1
          ND1=NMRKS(I,ID1)
          IF(IQ.GE.1) GOTO 100
        ENDIF
      ENDDO
 100  CONTINUE
C..   Now we compare NDET2 to NDET1
      IQ=0
      DO I=1,NEL
        K=0
        DO J=1,NEL
          IF(NMRKS(I,ID2).EQ.NMRKS(J,ID1)) CONTINUE
          IF(NMRKS(I,ID2).NE.NMRKS(J,ID1)) K=K+1
        ENDDO
        IF(K.EQ.NEL) THEN
          IQ=IQ+1
          ND2=NMRKS(I,ID2)
          IF(IQ.GE.1) GOTO 200
        ENDIF
      ENDDO
 200  CONTINUE
C..
      RETURN
      END
