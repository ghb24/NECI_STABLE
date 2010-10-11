Program ModelFCIQMC
    IMPLICIT NONE
    INTEGER :: NDet=100
    INTEGER :: LScr
    REAL*8, PARAMETER :: Tau=0.05
    REAL*8, PARAMETER :: SftDamp=0.05
    INTEGER, PARAMETER :: StepsSft=25
    INTEGER, PARAMETER :: NMCyc=1000000
    INTEGER, PARAMETER :: InitialWalk=1 
    INTEGER, PARAMETER :: TargetWalk=1000
    REAL*8, PARAMETER :: InitialShift=0.D0
    REAL*8 :: Norm,GrowRate,GroundShift,rat,r,Norm1
    INTEGER :: GroundTotParts,Iter,OldGroundTotParts
    REAL*8 :: Check
    INTEGER :: ierr,i,j,Die,Create,k
    ! Arrays for storing info.  All dimensions are NDet.
    REAL*8, allocatable :: EigenVec(:,:),EValues(:),KMat(:,:) 
    INTEGER, allocatable :: WalkListGround(:),WalkListGroundSpawn(:)
    INTEGER*8, allocatable :: SumWalkListGround(:)
    ! workspace for lapack.  LScr = 4*NDet is more than sufficient
    REAL*8, allocatable :: Scr(:) 
    LOGICAL :: tFixedShift, taverage = .false.
    integer :: nVaryShiftCycles
    real*8 :: AverageShift, SumShift

    integer :: iargc
    character(255) :: input_file
    logical :: tinput_file = .false.

!Initialise rand
    call random_seed()

    if (iargc() == 1) then
        call getarg(1,input_file)
        write (6,*) 'Reading Hamiltonian from file: ',trim(input_file)
        ! Have been given an input file containing the Hamiltonian.
        inquire(file=input_file,exist=tinput_file)
        if (tinput_file) then
            open(13,file=input_file,status='old',form='formatted')
            read(13,*) NDet
            close(13,status='keep')
        else
            write (6,*) 'Input file does not exist.'
            stop
        end if
    end if

    !Array allocation.
    LScr = 4*NDet
    allocate(KMat(NDet,NDet))
    allocate(EigenVec(NDet,NDet))
    allocate(EValues(NDet))
    allocate(WalkListGround(NDet))
    allocate(WalkListGroundSpawn(NDet))
    allocate(SumWalkListGround(NDet))
    allocate(Scr(LScr))

!Set up KMat
    if (tinput_file) then
        call ReadInKMat(NDet,input_file,KMat)
    else
        WRITE(6,*) "Setting up random K-Matrix..."
        CALL SetUpKMat(NDet,KMat)
    end if

!Diagonalise KMat
    EigenVec(:,:)=KMat(:,:)
    WRITE(6,*) "Diagonalising K-Matrix..."
    CALL DSYEV('V','U',NDet,EigenVec,NDet,Evalues,Scr,LScr,ierr)
    IF(ierr.ne.0) THEN
        STOP 'Error diagonalising matrix'
    ENDIF
    OPEN(9,FILE="Eigenvectors",STATUS='UNKNOWN')
    do i=1,NDet
        WRITE(9,"(I5)",advance='no') i
        do j=1,9
            WRITE(9,"(F20.12)",advance='no') EigenVec(i,j)
        enddo
        WRITE(9,"(F20.12)") EigenVec(i,10)
    enddo
    CLOSE(9)

    WRITE(6,*) "Lowest eigenvalues: "
    OPEN(9,FILE="Eigenvalues",STATUS='UNKNOWN')
    do i=1,min(10,NDet)
        WRITE(9,*) i,EValues(i)
        WRITE(6,*) i,EValues(i)
    enddo
    CLOSE(9)

    WRITE(6,*) "Performing Spawning..."

!Setup spawning
    WalkListGround(:)=0
    GroundShift=InitialShift
    SumShift = 0.0d0
    nVaryShiftCycles = 0
    AverageShift = 0.0d0
    SumWalkListGround(:)=0
    OPEN(12,FILE='ModelFCIMCStats',STATUS='unknown')

    WalkListGround(1)=InitialWalk
    tFixedShift=.true.
    OldGroundTotParts=InitialWalk
    
    do Iter=1,NMCyc

!Every so often, update the shift
        IF((mod(Iter,StepsSft).eq.0).and.(Iter.ne.1)) THEN

            GrowRate=REAL(GroundTotParts,8)/REAL(OldGroundTotParts,8)
            IF(.not.tFixedShift) THEN
                GroundShift=GroundShift-(log(GrowRate)*SftDamp)/(Tau*(StepsSft+0.D0))
                ! Start averaging?
                if (abs((GroundShift-EValues(1))/EValues(1)) < 0.05d0) then
                    taverage = .true.
!                    write (6,*) GroundShift,EValues(1),abs((GroundShift-EValues(1))/EValues(1))
!                    stop
                end if
                if (taverage) then
                    SumShift = SumShift + GroundShift
                    nVaryShiftCycles = nVaryShiftCycles + 1
                    AverageShift = SumShift/nVaryShiftCycles
                end if
            ELSE
                IF(GroundTotParts.ge.TargetWalk) THEN
                    tFixedShift=.false.
                ENDIF
            ENDIF
            OldGroundTotParts=GroundTotParts

            !Write out stats
            WRITE(6,"(I8,2F25.12,I15,I7)") Iter,GroundShift,AverageShift,GroundTotParts,WalkListGround(1)
            WRITE(12,"(I8,F25.12,I15,I7)") Iter,GroundShift,GroundTotParts,WalkListGround(1)

            Norm=0.D0
            Norm1=0.D0
            do i=1,NDet
                Norm=Norm+REAL(WalkListGround(i),8)**2
                Norm1=Norm1+REAL(SumWalkListGround(i),8)**2
            enddo
            Norm=SQRT(Norm)
            Norm1=SQRT(Norm1)
            OPEN(13,FILE='GroundWavevec',STATUS='unknown')
!            Check=0.D0
            do i=1,NDet
!                Check=Check+(REAL(SumWalkListGround(i),8)/Norm1)**2
                WRITE(13,*) i,REAL(WalkListGround(i),8)/Norm,REAL(SumWalkListGround(i),8)/Norm1
            enddo
!            WRITE(6,*) "Check = ",Check
            CLOSE(13)
        ENDIF

!Rezero spawning arrays    
        WalkListGroundSpawn(:)=0

!Simulate dynamic for calculation of GS
        do i=1,NDet     !Run through all determinants
!            WRITE(6,*) "Determinant: ",i

            do j=1,abs(WalkListGround(i))  !Run through all walkers

!Simulate full spawning by running through all connections.
                do k=1,NDet

                    IF(KMat(i,k).eq.0.D0) CYCLE
                    IF(i.eq.k) CYCLE

                    rat=abs(Tau*KMat(i,k))
!                    WRITE(6,*) rat
                    Create=INT(rat)
                    rat=rat-REAL(Create)
                    call random_number(r)
                    IF(rat.gt.r) THEN
!                        WRITE(6,*) "CREATED PARTICLE"
                        Create=Create+1
                    ENDIF

                    !create particles.
                    IF(KMat(i,k).gt.0.D0) THEN
                        !Flip child sign
                        IF(WalkListGround(i).lt.0) THEN
                            !Positive children
                            WalkListGroundSpawn(k)=WalkListGroundSpawn(k)+Create
                        ELSE
                            WalkListGroundSpawn(k)=WalkListGroundSpawn(k)-Create
                        ENDIF
                    ELSE
                        !Same sign as parent
                        IF(WalkListGround(i).gt.0) THEN
                            !Positive children
                            WalkListGroundSpawn(k)=WalkListGroundSpawn(k)+Create
                        ELSE
                            WalkListGroundSpawn(k)=WalkListGroundSpawn(k)-Create
                        ENDIF
                    ENDIF

                enddo

            enddo

!            do j=1,abs(WalkListGround(i))   
!!Each walker has unique probability to die
!                rat=Tau*(KMat(i,i)-GroundShift)
!                Die=INT(abs(rat))
!                IF(rat.lt.0.D0) THEN
!                    Die=-Die
!                    rat=rat+REAL(Die)
!                ELSE
!                    rat=rat-REAL(Die)
!                ENDIF
!                call random_number(r)
!                IF(abs(rat).gt.r) THEN
!                    IF(rat.gt.0.D0) THEN
!                        Die=Die+1
!                    ELSE
!                        Die=Die-1
!                    ENDIF
!                ENDIF
!                IF(WalkListGround(i).gt.0) THEN
!                    WalkListGround(i)=WalkListGround(i)-Die
!                ELSE
!                    WalkListGround(i)=WalkListGround(i)+Die
!                ENDIF
!            enddo

!Attempt to die simultaneously to all particles
            rat=REAL(abs(WalkListGround(i)),8)*Tau*(KMat(i,i)-GroundShift)
            Die=INT(rat)
            rat=rat-REAL(Die)
            call random_number(r)
            IF(abs(rat).gt.r) THEN
                IF(rat.gt.0.D0) THEN
                    Die=Die+1
                ELSE
                    Die=Die-1
                ENDIF
            ENDIF
            IF(Die.gt.abs(WalkListGround(i))) STOP 'Trying to create anti-particles'
            IF(WalkListGround(i).gt.0) THEN
                WalkListGround(i)=WalkListGround(i)-Die
            ELSE
                WalkListGround(i)=WalkListGround(i)+Die
            ENDIF

        enddo

!Combine lists (annihilation)
        GroundTotParts=0
        do i=1,NDet
            WalkListGround(i)=WalkListGround(i)+WalkListGroundSpawn(i)
            SumWalkListGround(i)=SumWalkListGround(i)+WalkListGround(i)
            GroundTotParts=GroundTotParts+abs(WalkListGround(i))
        enddo
        IF(GroundTotParts.eq.0) THEN
            WRITE(6,*) "ALL WALKERS DIED - RESTARTING"
            GroundShift=InitialShift
            WalkListGround(1)=InitialWalk
            tFixedShift=.true.
            OldGroundTotParts=InitialWalk
        ENDIF

    enddo
    CLOSE(12)
        


End Program ModelFCIQMC


SUBROUTINE SetUpKMat(NDet,KMat)
    ! Sets up a random K Matrix.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDet
    REAL*8, INTENT(OUT) :: KMat(NDet,NDet)
    INTEGER :: i,j
    REAL*8 :: StartEl,EndEl,Step,ProbNonZero,OffDiagEl
    REAL*8 :: r

    KMat(:,:)=0.D0

    StartEl=0.5
    EndEl=10.D0
    Step=(EndEl-StartEl)/REAL(NDet-1,8)
    KMat(2,2)=StartEl
    do i=3,NDet
        KMat(i,i)=KMat(i-1,i-1)+Step
    enddo

    WRITE(6,*) "RefDet = ", KMat(1,1)
    WRITE(6,*) "MaxDet = ", KMat(NDet,NDet)

    ProbNonZero=0.2     !This is probability that off-diagonal matrix elements are non-zero
    OffDiagEl=5.D-2     !This is the magnitude of the off-diagonal matrix elements.
    do i=1,NDet
        do j=1,i-1
            call random_number(r)
            IF(r.gt.ProbNonZero) THEN

                IF(r.gt.0.51) THEN
                    !Matrix element is negative with probability 0.49
                    KMat(i,j)=-OffDiagEl
                    KMat(j,i)=-OffDiagEl
                ELSE
                    !Matrix element is positive with probability 0.51
                    KMat(i,j)=OffDiagEl
                    KMat(j,i)=OffDiagEl
                ENDIF
            ENDIF
        enddo
    enddo

END SUBROUTINE SetUpKMat

SUBROUTINE ReadInKMat(NDet,input_file,KMat)
    ! Reads in a Hamiltonian matrix from
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDet
    character(255), intent(in) :: input_file
    REAL*8, INTENT(OUT) :: KMat(NDet,NDet)
    integer :: i,j
    real*8 :: H00

    ! The first line just contains the size of the Hamiltonian matrix.
    ! This has already been read...
    open(13,file=input_file,status='old',form='formatted')
    read(13,*) i

    do i = 1,NDet
        read (13,*) (KMat(i,j),j=1,NDet)
    end do

    close(13,status='keep')

    ! K_ij = H_ij - H00 d_ij
    H00 = KMat(1,1)
    forall (i=1:NDet) KMat(i,i) = KMat(i,i) - H00

end subroutine ReadInKMat
