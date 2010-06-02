Program ModelFCIQMC

    IMPLICIT NONE

    INTEGER :: NDet=25
    INTEGER, PARAMETER :: lenof_sign=2   !Number of integers needed to store a walker
    REAL*8, PARAMETER :: Tau=0.05
    REAL*8, PARAMETER :: SftDamp=0.05
    INTEGER, PARAMETER :: StepsSft=25
    INTEGER, PARAMETER :: NMCyc=100000
    INTEGER, PARAMETER :: InitialWalk=1 
    INTEGER, PARAMETER :: TargetWalk=5000
    REAL*8, PARAMETER :: InitialShift=0.D0
    INTEGER, PARAMETER :: dp=8

    CALL RunModelFCIQMC()

CONTAINS

    SUBROUTINE RunModelFCIQMC()
        IMPLICIT NONE
        INTEGER :: LScr,LWork2
        REAL*8 :: Norm,RealGrowRate,ImagGrowRate,RealShift,ImagShift,GrowRate,GroundShift,rat,r,Norm1
        INTEGER :: Iter,OldBothTotParts,BothTotParts
        INTEGER, DIMENSION(lenof_sign) :: GroundTotParts,OldGroundTotParts
        REAL*8 :: Check
        INTEGER :: ierr,i,j,k
        INTEGER, DIMENSION(lenof_sign) :: Child,nDie
        ! Arrays for storing info.  All dimensions are NDet.
        REAL*8, allocatable :: EValues(:),Work2(:)
        COMPLEX*16, allocatable :: KMat(:,:),EigenVec(:,:)
        INTEGER, allocatable :: WalkListGround(:,:),WalkListGroundSpawn(:,:)
        INTEGER*8, allocatable :: SumWalkListGround(:,:)
        ! workspace for lapack.  LScr = 4*NDet is more than sufficient
        COMPLEX*16, allocatable :: Scr(:) 
        LOGICAL :: tFixedShift, taverage = .false.
        integer :: nVaryShiftCycles
        real*8 :: AverageShift, SumShift,fac
        integer :: iargc,part_type
        character(255) :: input_file
        logical :: tinput_file = .false.
        COMPLEX*16 :: rh

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
        LWork2 = 3*NDet-2
        allocate(KMat(NDet,NDet))
        allocate(EigenVec(NDet,NDet))
        allocate(EValues(NDet))
        allocate(WalkListGround(lenof_sign,NDet))
        allocate(WalkListGroundSpawn(lenof_sign,NDet))
        allocate(SumWalkListGround(lenof_sign,NDet))
        allocate(Scr(LScr))
        allocate(Work2(LWork2))

    !Set up KMat
        if (tinput_file) then
            call ReadInKMat(input_file,KMat)
        else
            WRITE(6,*) "Setting up random K-Matrix..."
            CALL SetUpKMat(KMat)
            CALL WriteOutKMat(KMat)
        end if

    !Diagonalise KMat
        EigenVec(:,:)=KMat(:,:)
        WRITE(6,*) "Diagonalising K-Matrix..."
        CALL ZHEEV('V','U',NDet,EigenVec,NDet,Evalues,Scr,LScr,Work2,ierr)
    !    CALL DSYEV('V','U',NDet,EigenVec,NDet,Evalues,Scr,LScr,ierr)
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
        WalkListGround(:,:)=0
        GroundShift=InitialShift
        SumShift = 0.0d0
        nVaryShiftCycles = 0
        AverageShift = 0.0d0
        SumWalkListGround(:,:)=0
        OPEN(12,FILE='ModelFCIMCStats',STATUS='unknown')
        WRITE(12,*) "#1.Iter  2.GroundShift  3.RealShift  4.ImagShift  5.AverageShift  6.RealParts  7.ImParts  8.ReRoot  9.ImRoot"
        WRITE(6,*) "1.Iter  2.GroundShift  3.RealShift  4.ImagShift  5.AverageShift  6.RealParts  7.ImParts  8.ReRoot  9.ImRoot"

        !Start with just real walkers
        WalkListGround(1,1)=InitialWalk
        tFixedShift=.true.
        OldGroundTotParts(1)=InitialWalk
        OldBothTotParts=InitialWalk
        
        do Iter=1,NMCyc

!Every so often, update the shift
            IF((mod(Iter,StepsSft).eq.0).and.(Iter.ne.1)) THEN

                RealGrowRate=REAL(GroundTotParts(1),8)/REAL(OldGroundTotParts(1),8)
                ImagGrowRate=REAL(GroundTotParts(2),8)/REAL(OldGroundTotParts(2),8)
                GrowRate=REAL(BothTotParts,8)/REAL(OldBothTotParts,8)
                IF(.not.tFixedShift) THEN
                    GroundShift=GroundShift-(log(GrowRate)*SftDamp)/(Tau*(StepsSft+0.D0))
                    RealShift=RealShift-(log(RealGrowRate)*SftDamp)/(Tau*(StepsSft+0.D0))
                    ImagShift=ImagShift-(log(ImagGrowRate)*SftDamp)/(Tau*(StepsSft+0.D0))
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
                    IF(BothTotParts.ge.TargetWalk) THEN
                        tFixedShift=.false.
                    ENDIF
                ENDIF
                OldGroundTotParts=GroundTotParts
                OldBothTotParts=BothTotParts

                !Write out stats
                WRITE(6,"(I8,4F25.12,2I15,2I7)") Iter,GroundShift,RealShift,ImagShift,AverageShift,GroundTotParts(1),GroundTotParts(2),WalkListGround(1,1),WalkListGround(2,1)
                WRITE(12,"(I8,4F25.12,2I15,2I7)") Iter,GroundShift,RealShift,ImagShift,AverageShift,GroundTotParts(1),GroundTotParts(2),WalkListGround(1,1),WalkListGround(2,1)
                CALL FLUSH(12)
                CALL FLUSH(6)

                Norm=0.D0
                Norm1=0.D0
                do i=1,NDet
                    Norm=Norm+(REAL(WalkListGround(1,i),8)**2)+(REAL(WalkListGround(2,i),8)**2)
                    Norm1=Norm1+(REAL(SumWalkListGround(1,i),8)**2)+(REAL(SumWalkListGround(2,i),8)**2)
                enddo
                Norm=SQRT(Norm)
                Norm1=SQRT(Norm1)
                OPEN(13,FILE='GroundWavevec',STATUS='unknown')
    !            Check=0.D0
                do i=1,NDet
    !                Check=Check+(REAL(SumWalkListGround(i),8)/Norm1)**2
                    WRITE(13,*) i,REAL(WalkListGround(1,i),8)/Norm,REAL(WalkListGround(2,i),8)/Norm,REAL(SumWalkListGround(1,i),8)/Norm1,REAL(SumWalkListGround(2,i),8)/Norm1
                enddo
    !            WRITE(6,*) "Check = ",Check
                CLOSE(13)

            ENDIF

!Rezero spawning arrays    
            WalkListGroundSpawn(:,:)=0

!Simulate dynamic for calculation
            do i=1,NDet     !Run through all determinants
    !            WRITE(6,*) "Determinant: ",i

                do part_type=1,lenof_sign   !Loop over real & imag walkers

                    do j=1,abs(WalkListGround(part_type,i))  !Run through all walkers

!Simulate full spawning by running through all connections.
                        do k=1,NDet

                            IF(KMat(i,k).eq.0.D0) CYCLE
                            IF(i.eq.k) CYCLE

                            rh=KMat(i,k)
                            !All spawning attempts here
                            child=Attempt_create(rh,i,k,WalkListGround(:,i),part_type)
            
                            !`Child` returns the number of real & im walkers to create at k
                            !create particles.
                            WalkListGroundSpawn(:,k)=WalkListGroundSpawn(:,k)+child

                        enddo

                    enddo

                enddo   !end looping over real & Im

!Attempt to die simultaneously to all particles
                fac=Tau*(REAL(KMat(i,i),dp)-GroundShift)
                ndie = Attempt_die(fac,WalkListGround(:,i))

                do part_type=1,lenof_sign
                    !Update lists with the number that have died
                    IF(WalkListGround(part_type,i).gt.0) THEN
                        WalkListGround(part_type,i)=WalkListGround(part_type,i)-ndie(part_type)
                    ELSE
                        WalkListGround(part_type,i)=WalkListGround(part_type,i)+ndie(part_type)
                    ENDIF
                enddo

            enddo

!Combine lists (annihilation)
            GroundTotParts=0
            do i=1,NDet
                WalkListGround(:,i)=WalkListGround(:,i)+WalkListGroundSpawn(:,i)
                SumWalkListGround(:,i)=SumWalkListGround(:,i)+WalkListGround(:,i)
                GroundTotParts(:)=GroundTotParts(:)+abs(WalkListGround(:,i))
            enddo

            BothTotParts=GroundTotParts(1)+GroundTotParts(2)

            IF(BothTotParts.eq.0) THEN
                WRITE(6,*) "ALL WALKERS DIED - RESTARTING"
                GroundShift=InitialShift
                WalkListGround(1,1)=InitialWalk
                tFixedShift=.true.
                OldGroundTotParts(1)=InitialWalk
                OldBothTotParts=InitialWalk
            ENDIF

        enddo
        CLOSE(12)

    End SUBROUTINE RunModelFCIQMC


    FUNCTION Attempt_die(fac,wSign) result(nDie)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: fac
        INTEGER, DIMENSION(lenof_sign), INTENT(IN) :: wSign
        INTEGER, DIMENSION(lenof_sign) :: nDie
        REAL*8 :: rat,r
        INTEGER :: part_type

        do part_type=1,lenof_sign

            rat=fac*abs(wSign(part_type))
            nDie(part_type) = int(rat)
            rat = rat - real(ndie(part_type), dp)
            ! Choose to die or not stochastically
            call random_number(r)
            if (abs(rat) > r) ndie(part_type) = ndie(part_type) + nint(sign(1.0_dp, rat))
            IF(ndie(part_type).gt.abs(wSign(part_type))) STOP 'Trying to create anti-particles'
        enddo

    END FUNCTION Attempt_die

    FUNCTION Attempt_create(rh,parent_ind,child_ind,parent_sign,part_type) result(child)
        IMPLICIT NONE
        COMPLEX*16, INTENT(IN) :: rh
        REAL*8 :: MatEl,rat,r 
        INTEGER, INTENT(IN) :: parent_ind,child_ind,part_type
        INTEGER, DIMENSION(lenof_sign), INTENT(IN) :: parent_sign
        INTEGER, DIMENSION(lenof_sign) :: child
        INTEGER :: ExtraCreate,i

        !We are dealing with spawning from real and imaginary elements, and assume
        !that rh is complex
        IF(part_type.eq.1) THEN
            !Real parent particle

            do i=1,lenof_sign
                !Run over spawnings from both the real and imaginary part of the matrix element

                IF(i.eq.1) THEN
                    !We want to use the real part of the matrix element to create real walkers
                    MatEl=REAL(rh,8)
                ELSE
                    !We want to use the imaginary part of the matrix element to create imaginary walkers
                    MatEl=AIMAG(rh)
                ENDIF

                !Attempt spawning
                rat = tau * abs(MatEl)

                ! If probability > 1, then we just create multiple children at the
                ! chosen determinant.
                extraCreate = int(rat)
                rat = rat - real(extraCreate, dp)

                ! Stochastically choose whether to create or not.
                call random_number(r)
                if (rat > r) then
                    !Create child
                    child(i) = -nint(sign(1.0_dp, parent_sign(part_type)*MatEl))  !Will return +- one depending on the desired sign of the stochastically created child.
                    child(i) = child(i) + sign(extraCreate, child(i))
                else
                    !Just return if any extra particles created
                    child(i) = -extraCreate*nint(sign(1.0_dp, parent_sign(part_type)*MatEl))
                endif
            enddo

        ELSE
            !Imaginary parent particle - rules are slightly different...
            !Attempt to spawn REAL walkers with prob +AIMAG(Hij)/P
            !Attempt to spawn IMAG walkers with prob -REAL(Hij)/P
            do i=1,lenof_sign
                !Run over spawnings from both the real and imaginary part of the matrix element

                IF(i.eq.1) THEN
                    !We want to use the imaginary part of the matrix element to create real walkers
                    MatEl=AIMAG(rh)
                ELSE
                    !We want to use the real part of the matrix element to create imaginary walkers
                    MatEl=REAL(rh,dp)
                ENDIF

                !Attempt spawning
                rat = tau * abs(MatEl)

                ! If probability > 1, then we just create multiple children at the
                ! chosen determinant.
                extraCreate = int(rat)
                rat = rat - real(extraCreate, dp)

                ! Stochastically choose whether to create or not.
                call random_number(r)
                IF(i.eq.1) THEN
                    !Prob = +AIMAG(Hij)/P to create real children
                    if (rat > r) then
                        !Create child
                        child(i) = nint(sign(1.0_dp, parent_sign(part_type)*MatEl))  !Will return +- one depending on the desired sign of the stochastically created child.
                        child(i) = child(i) + sign(extraCreate, child(i))
                    else
                        !Just return if any extra particles created
                        child(i) = extraCreate*nint(sign(1.0_dp, parent_sign(part_type)*MatEl))
                    endif
                ELSE
                    !Prob = -REAL(Hij)/P to create imaginary children
                    if (rat > r) then
                        !Create child
                        child(i) = -nint(sign(1.0_dp, parent_sign(part_type)*MatEl))  !Will return +- one depending on the desired sign of the stochastically created child.
                        child(i) = child(i) + sign(extraCreate, child(i))
                    else
                        !Just return if any extra particles created
                        child(i) = -extraCreate*nint(sign(1.0_dp, parent_sign(part_type)*MatEl))
                    endif
                ENDIF
            enddo

        ENDIF   ! Type of parent

    END FUNCTION Attempt_create


    SUBROUTINE SetUpKMat(KMat)
        ! Sets up a random COMPLEX K Matrix.
        IMPLICIT NONE
        COMPLEX*16, INTENT(OUT) :: KMat(NDet,NDet)
        INTEGER :: i,j,l
        REAL*8 :: StartEl,EndEl,Step,ProbNonZero,OffDiagEl
        REAL*8 :: r

        KMat(:,:)=CMPLX(0.D0,0.D0)

        StartEl=0.5     !Initial diagonal matrix element - these must all be real obv.
        EndEl=10.D0
        Step=(EndEl-StartEl)/REAL(NDet-1,8) !Rate of increase of diagonal matrix elements
        KMat(2,2)=CMPLX(StartEl,0.D0)
        do i=3,NDet
            KMat(i,i)=CMPLX(REAL(KMat(i-1,i-1),8)+Step,0.D0)
        enddo

        WRITE(6,*) "RefDet = ", KMat(1,1)
        WRITE(6,*) "MaxDet = ", KMat(NDet,NDet)

        ProbNonZero=0.4     !This is probability that off-diagonal matrix elements are non-zero
        OffDiagEl=5.D-2     !This is the magnitude of the off-diagonal matrix elements.
        do l=1,2            !loop over real and imaginary parts to allow a chance for both the real and imaginary parts to become non-zero
            do i=1,NDet
                do j=1,i-1
                    call random_number(r)
                    IF(r.gt.ProbNonZero) THEN

                        IF(r.gt.0.51) THEN
                            !Matrix element is negative with probability 0.49
                            IF(l.eq.1) THEN
                                !Real part
                                KMat(i,j)=CMPLX(-OffDiagEl,0.D0)
                                KMat(j,i)=CMPLX(-OffDiagEl,0.D0)
                            ELSE
                                !Imaginary part - must ensure hermicity
                                KMat(i,j)=CMPLX(REAL(KMat(i,j),8),-OffDiagEl)
                                KMat(j,i)=CMPLX(REAL(KMat(i,j),8),OffDiagEl)
                            ENDIF
                        ELSE
                            !Matrix element is positive with probability 0.51
                            IF(l.eq.1) THEN
                                !Real part
                                KMat(i,j)=CMPLX(OffDiagEl,0.D0)
                                KMat(j,i)=CMPLX(OffDiagEl,0.D0)
                            ELSE
                                KMat(i,j)=CMPLX(REAL(KMat(i,j),8),OffDiagEl)
                                KMat(j,i)=CMPLX(REAL(KMat(i,j),8),OffDiagEl)
                            ENDIF
                        ENDIF
                    ENDIF
                enddo
            enddo
        enddo

    END SUBROUTINE SetUpKMat

    SUBROUTINE WriteoutKMat(KMat)
        IMPLICIT NONE
        INTEGER :: i,j
        COMPLEX*16, INTENT(IN) :: KMat(NDet,NDet)

        WRITE(6,*) "*** Writing out KMat ***"

        do i=1,NDet
            do j=1,NDet
                write(6,'("(")',advance='no')
                write(6,"(G16.7)",advance='no') KMat(i,j)
                write(6,'(")")',advance='no')
            enddo
            WRITE(6,*) ""
        enddo

        CALL FLUSH(6)

    END SUBROUTINE WriteoutKMat

    SUBROUTINE ReadInKMat(input_file,KMat)
        ! Reads in a Hamiltonian matrix from
        IMPLICIT NONE
        character(255), intent(in) :: input_file
        COMPLEX*16, INTENT(OUT) :: KMat(NDet,NDet)
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

END Program ModelFCIQMC
