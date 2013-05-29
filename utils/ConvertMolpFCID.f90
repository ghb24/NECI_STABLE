PROGRAM ConvertMolpFCID 

    IMPLICIT NONE
    INTEGER :: NORB,NELEC,MS2,ORBSYM(500),ISYM,SYMMAX,i,j,a,k,l,Elec,nISym,UnpairedElec
    INTEGER :: iPairs,iSize,ierr,UMatInd
    INTEGER , ALLOCATABLE :: nI(:),SymCounts(:),AlphaOcc(:),BetaOcc(:),SymNum(:)
    REAL*8 , ALLOCATABLE :: UMAT(:),ARR(:),TMAT(:,:)
    REAL*8 :: Z,ECore
    LOGICAL :: tROHF,exists
    NAMELIST /FCI/ NORB,NELEC,MS2,ORBSYM,ISYM

    ORBSYM(:)=0

    OPEN(8,FILE='FCIDUMP',STATUS='OLD',FORM='FORMATTED')
    READ(8,FCI)

    IF(NORB.gt.500) STOP 'Too Many Orbitals'
    IF(MS2.ne.0) THEN
        WRITE(6,*) "High spin HF determinant - converting orbitals to spin orbitals"
        WRITE(6,*) "I am assuming that this is an ROHF calculation, rather than UHF!"
        tROHF=.true.
    ELSE
        WRITE(6,*) "RESTRICTED HF calculation - orbitals left as spatial orbitals"
        tROHF=.false.
    ENDIF


    SYMMAX=0
    do i=1,NORB
        IF(ORBSYM(i).gt.SYMMAX) SYMMAX=ORBSYM(i)
    enddo
    WRITE(6,*) SYMMAX," symmetries found."

    ALLOCATE(SymCounts(SYMMAX))
    ALLOCATE(SymNum(SYMMAX))
    ALLOCATE(nI(NELEC))
    ALLOCATE(AlphaOcc(SYMMAX))
    ALLOCATE(BetaOcc(SYMMAX))

    nI(:)=0

!Symcounts(1:SYMMAX) gives the cumulative index of orbitals in each symmetry
    SymNum(:)=0
    do i=1,NORB
        SymNum(ORBSYM(i))=SymNum(ORBSYM(i))+1
    enddo

    SymCounts(1)=1

    do i=2,SYMMAX
        SymCounts(i)=SymNum(i-1)+SymCounts(i-1)
    enddo

    INQUIRE(FILE='DetSpec',EXIST=Exists)
    IF(Exists) THEN
        WRITE(6,*) "DetSpec file found. Reading determinant"
        OPEN(11,FILE='DetSpec',STATUS='OLD')
        READ(11,*) nI(:)
        UnpairedElec=0
        do i=1,NELEC
            IF((mod(nI(i),2).eq.1).and.(i.eq.NELEC)) THEN
                UnpairedElec=UnpairedElec+1
            ELSEIF((mod(nI(i),2).eq.1).and.(nI(i+1).ne.(nI(i)+1))) THEN
                UnpairedElec=UnpairedElec+1
            ENDIF
        enddo
        CLOSE(11,status='delete')
        Elec=NELEC+1
        nISym=ISYM-1

    ELSE

        INQUIRE(FILE='MolpStats',EXIST=Exists)
        IF(Exists) THEN
            OPEN(11,FILE='MolpStats',STATUS='OLD')
            IF(tROHF) THEN
                READ(11,*) AlphaOcc(1:SYMMAX)
                READ(11,*) BetaOcc(1:SYMMAX)
            ELSE
                READ(11,*) BetaOcc(1:SYMMAX)
                AlphaOcc(:)=BetaOcc(:)
            ENDIF
            CLOSE(11,status='delete')
        ELSE
            STOP 'Cannot find MolpStats file'
        ENDIF

        Elec=1

        nISym=0
        UnpairedElec=0

        do i=1,SYMMAX
            IF(mod(AlphaOcc(i)+BetaOcc(i),2).eq.1) THEN
                UnpairedElec=UnpairedElec+1
            ENDIF
            do j=1,AlphaOcc(i)
                nI(Elec)=(SymCounts(i)+j-1)*2-1
                nISym=IEOR(nISym,i-1)
                Elec=Elec+1
                IF(Elec.gt.NELEC) GOTO 20
            enddo
            do j=1,BetaOcc(i)
                nI(Elec)=(SymCounts(i)+j-1)*2
                nISym=IEOR(nISym,i-1)
                Elec=Elec+1
                IF(Elec.gt.NELEC) GOTO 20
            enddo
        enddo

20      CONTINUE
        CALL NECI_SORTI(NELEC,nI)
    ENDIF

    WRITE(6,*) "Spin orbital representation of HF det is: "
    do i=1,NELEC
        WRITE(6,"(I5)",advance='no') nI(i)
    enddo
    WRITE(6,*) ""

    IF(Elec.ne.NELEC+1) STOP 'Incorrect number of electrons when finding nI'
    IF(UnpairedElec.ne.MS2) THEN
        WRITE(6,*) UnpairedElec,MS2
        STOP 'Incorrect spin in HF Det'
    ENDIF
    IF((nISym+1).ne.ISYM) STOP 'HF det has wrong sym'

    iPairs=(NORB*(NORB+1))/2
    iSize=(iPairs*(iPairs+1))/2

    ALLOCATE(UMAT(0:iSize),stat=ierr)
    IF(ierr.ne.0) STOP 'Allocate error'
    IF(tROHF) THEN
        ALLOCATE(ARR(NORB*2),stat=ierr)
    ELSE
        ALLOCATE(ARR(NORB),stat=ierr)
    ENDIF
    ALLOCATE(TMAT(0:NORB,0:NORB),stat=ierr)
    IF(ierr.ne.0) STOP 'Allocate error'
    UMAT(:)=0.D0
    TMAT(:,:)=0.D0
    ARR(:)=0.D0

    OPEN(12,FILE='NeciFCIDUMP',STATUS='UNKNOWN')

    IF(tROHF) THEN
        WRITE(12,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',NORB*2,'NELEC=',NELEC,',MS2=',MS2,','
    ELSE
        WRITE(12,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',NORB,'NELEC=',NELEC,',MS2=',MS2,','
    ENDIF
    WRITE(12,'(A9)',advance='no') 'ORBSYM='
    DO i=1,NORB
        IF(tROHF) THEN
            WRITE(12,'(I1,A1,I1,A1)',advance='no') ORBSYM(i),',',ORBSYM(i),','
        ELSE
            WRITE(12,'(I1,A1)',advance='no') ORBSYM(i),','
        ENDIF
    ENDDO
    WRITE(12,*)
    IF(tROHF) THEN
        WRITE(12,'(A7,I1,A12)') 'ISYM=',ISYM,' UHF=.TRUE.'
        WRITE(12,'(A5)') '&END'
    ELSE
        WRITE(12,'(A7,I1,A12)') 'ISYM=',ISYM,' UHF=.FALSE.'
        WRITE(12,'(A5)') '&END'
    ENDIF

    DO WHILE(.true.)

        READ(8,'(1X,G20.12,4I3)',END=99) Z,I,J,K,L    

        IF(I.ne.0.and.J.eq.0.and.K.eq.0) THEN
            !Fock energies
            STOP 'Fock energies should not be being printed!'

        ELSEIF(I.ne.0.and.J.ne.0.and.K.eq.0.and.L.eq.0) THEN
            !One-electron energies

            TMAT(I,J)=Z
            TMAT(J,I)=Z
            IF(tROHF) THEN
                WRITE(12,'(1X,G20.12,4I3)') Z,I*2,J*2,0,0
                WRITE(12,'(1X,G20.12,4I3)') Z,I*2-1,J*2-1,0,0
            ELSE
                WRITE(12,'(1X,G20.12,4I3)') Z,I,J,K,L
            ENDIF
        ELSEIF(I.eq.0) THEN
            !Core energy
            ECore=Z
        ELSE
            !2e int
            UMAT(UMatInd(I,J,K,L))=Z
            IF(abs(Z).gt.1.D-08) THEN
                IF(tROHF) THEN
                    WRITE(12,'(1X,G20.12,4I3)') Z,I*2,J*2,K*2,L*2
                    WRITE(12,'(1X,G20.12,4I3)') Z,I*2-1,J*2-1,K*2-1,L*2-1
                    WRITE(12,'(1X,G20.12,4I3)') Z,I*2,J*2,K*2-1,L*2-1
                    WRITE(12,'(1X,G20.12,4I3)') Z,I*2-1,J*2-1,K*2,L*2
                ELSE
                    WRITE(12,'(1X,G20.12,4I3)') Z,I,J,K,L
                ENDIF
            ENDIF
        ENDIF

    END DO

99  CONTINUE

!Integrals have now been read in
    CLOSE(8)

!Create fock energies (diagonal fock matrix elements)

    IF(tROHF) THEN
        !Work in spin orbitals
        do i=1,NORB*2
            IF(mod(i,2).eq.0) THEN
                !even orbital
                Arr(i)=TMAT(i/2,i/2)
                do a=1,NELEC
                    IF(mod(nI(a),2).eq.0) THEN
                        !even elec
                        Arr(i)=Arr(i)+(UMAT(UMatInd(i/2,i/2,nI(a)/2,nI(a)/2)))!/2.D0
                        Arr(i)=Arr(i)-(UMAT(UMatInd(i/2,nI(a)/2,nI(a)/2,i/2)))!/2.D0
                    ELSE
                        !odd elec - no exchange term
                        Arr(i)=Arr(i)+(UMAT(UMatInd(i/2,i/2,(nI(a)+1)/2,(nI(a)+1)/2)))!/2.D0
                    ENDIF
                enddo
            ELSE
                !odd orbital
                Arr(i)=TMAT((i+1)/2,(i+1)/2)
                do a=1,NELEC
                    IF(mod(nI(a),2).eq.0) THEN
                        !even elec - no exchange term
                        Arr(i)=Arr(i)+(UMAT(UMatInd((i+1)/2,(i+1)/2,nI(a)/2,nI(a)/2)))!/2.D0
                    ELSE
                        !odd elec
                        Arr(i)=Arr(i)+(UMAT(UMatInd((i+1)/2,(i+1)/2,(nI(a)+1)/2,(nI(a)+1)/2)))!/2.D0
                        Arr(i)=Arr(i)-(UMAT(UMatInd((i+1)/2,(nI(a)+1)/2,(nI(a)+1)/2,(i+1)/2)))!/2.D0
                    ENDIF
                enddo
            ENDIF
        enddo

    ELSE
        !Work in spatial orbitals
        do i=1,NORB
            Arr(i)=TMAT(I,I)
            do a=1,NELEC,2
                !Add coulomb term
                Arr(i)=Arr(i)+(UMAT(UMatInd(i,i,(nI(a)+1)/2,(nI(a)+1)/2)))*2.D0
                !Add exchange term
                Arr(i)=Arr(i)-UMAT(UMatInd(i,(nI(a)+1)/2,(nI(a)+1)/2,i))
            enddo
        enddo
    ENDIF

    IF(tROHF) THEN
        WRITE(6,*) "Fock Energies are (spin orbs): "
        do i=1,NORB*2
            WRITE(6,"(i8,G25.17)") i,Arr(i)
            WRITE(12,'(1X,G20.12,4I3)') Arr(i),I,0,0,0
        enddo
    ELSE
        WRITE(6,*) "Fock Energies are (spatial orbs): "
        do i=1,NORB
            WRITE(6,"(i8,G25.17)") i,Arr(i)
            WRITE(12,'(1X,G20.12,4I3)') Arr(i),I,0,0,0
        enddo
    ENDIF


!Print out core energy
    WRITE(12,'(1X,G20.12,4I3)') ECore,0,0,0,0
    CLOSE(12)

    DEALLOCATE(UMAT)
    DEALLOCATE(ARR)
    DEALLOCATE(TMAT)


END PROGRAM



FUNCTION UMatInd(I,J,K,L)
    IMPLICIT NONE
    INTEGER :: I,J,K,L,A,B,UMatInd

    IF(I.eq.0.or.J.eq.0.or.K.eq.0.or.L.eq.0) THEN
        UMatInd=0
        RETURN
    ENDIF

    IF(I.gt.J) THEN
        A=(I*(I-1))/2+J
    ELSE
        A=(J*(J-1))/2+I
    ENDIF

    IF(K.gt.L) THEN
        B=(K*(K-1))/2+L
    ELSE
        B=(L*(L-1))/2+K
    ENDIF
    IF(A.gt.B) THEN
        UMatInd=(A*(A-1))/2+B
    ELSE
        UMatInd=(B*(B-1))/2+A
    ENDIF

END FUNCTION UMatInd

      SUBROUTINE NECI_SORTI(N,RA)
      IMPLICIT NONE 
      INTEGER N,I,L,IR,J
      INTEGER RA(N),RRA
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
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
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END

