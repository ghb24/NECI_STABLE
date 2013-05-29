Program ModelFCIQMC

    IMPLICIT NONE

    !Parameters to control hamiltonian setup
    REAL*8, PARAMETER :: StartEl=0.2     !Initial diagonal matrix element - these must all be real obv.
    REAL*8, PARAMETER :: EndEl=7.D0
    !Off diagonal matrix elements
    REAL*8, PARAMETER :: ProbNonZero=0.3     !This is probability that off-diagonal matrix elements are non-zero
    REAL*8, PARAMETER :: OffDiagEl=8.D-2     !This is the magnitude of the off-diagonal matrix elements.

    INTEGER :: NDet=28 

    INTEGER, PARAMETER :: lenof_sign=2   !Number of integers needed to store a walker
    REAL*8, PARAMETER :: Tau=0.02 
    REAL*8, PARAMETER :: SftDamp=0.2 
    INTEGER, PARAMETER :: StepsSft=50 
    INTEGER, PARAMETER :: NMCyc=100000
    INTEGER, PARAMETER :: InitialWalk=1 
    INTEGER, PARAMETER :: TargetWalk=3000
    REAL*8, PARAMETER :: InitialShift=0.D0
    INTEGER, PARAMETER :: dp=8
    LOGICAL, PARAMETER :: tRotateWavefunction=.false. 
    LOGICAL, PARAMETER :: tSeperateShift=.false. 
    LOGICAL, PARAMETER :: tKeepWalkerFiles=.false. 
    LOGICAL, PARAMETER :: tDumpKMat=.false.
    REAL*8, PARAMETER :: pi=3.14159265358979323846264338327950288419716939937510D0
    REAL*8 :: piby2=pi/2.D0
    REAL*8 :: pi2=pi*2.D0

    !Run-time variables
    REAL*8 :: GroundShift,ImagShift,RealShift

    CALL RunModelFCIQMC()

CONTAINS

    SUBROUTINE RunModelFCIQMC()
        IMPLICIT NONE
        INTEGER :: LScr,LWork2
        REAL*8 :: Norm,RealGrowRate,ImagGrowRate,GrowRate,rat,r,Norm1
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
        real*8 :: AverageShift, SumShift,fac,ProjE
        integer :: iargc,part_type
        integer, allocatable :: Seed(:)
        character(255) :: input_file
        logical :: tinput_file = .false.
        COMPLEX*16 :: rh
        character(len=26) :: abstr

    !Initialise rand
        call random_seed
        call random_seed(SIZE = K)
        WRITE(*,*) ' Number of integers for starting value = ', K
        ALLOCATE(Seed(1:K))
        Seed(1)=73
        CALL random_seed(PUT=Seed(1:K))
        call random_number(r)
        WRITE(6,*) "Random number seed is: ",Seed(1:K)
        WRITE(6,*) "First random number is: ",r
        call random_number(r)
        WRITE(6,*) "Second random number is: ",r

        IF(tRotateWavefunction.and.tSeperateShift) STOP 'Should not have rotatewavefunction AND seperateshift both on'

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
        
        WRITE(6,*) "NDet = ",NDet

    !Set up KMat
        if (tinput_file) then
            call ReadInKMat(input_file,KMat)
            CALL WriteOutKMat(KMat)
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
            WRITE(9,"(2F20.12,A)",advance='no') REAL(EigenVec(i,1),dp), AIMAG(EigenVec(i,1)), "  *  "
            WRITE(9,"(F20.12)",advance='no') REAL(EigenVec(i,2),dp)
            WRITE(9,"(F20.12)") AIMAG(EigenVec(i,2))
        enddo
        CLOSE(9)
        
        WRITE(6,*) "Lowest eigenvalues: "
        OPEN(9,FILE="Eigenvalues",STATUS='UNKNOWN')
        do i=1,min(10,NDet)
            WRITE(9,*) i,EValues(i)
            WRITE(6,*) i,EValues(i)
        enddo
        CLOSE(9)

        if(tdumpkmat) THEN
            call dumpkmat(KMat)
        endif

        WRITE(6,*) "Performing Spawning..."

    !Setup spawning
        WalkListGround(:,:)=0
        GroundShift=InitialShift
        SumShift = 0.0d0
        nVaryShiftCycles = 0
        AverageShift = 0.0d0
        SumWalkListGround(:,:)=0
        OPEN(12,FILE='ModelFCIMCStats',STATUS='unknown')
        WRITE(12,'(A)') "#1.Iter  2.GroundShift  3.RealShift  4.ImagShift  5.AverageShift  6.RealParts  7.ImParts  8.ReRoot  9.ImRoot  10.ProjE"
        WRITE(6,*) "1.Iter  2.GroundShift  3.RealShift  4.ImagShift  5.AverageShift  6.RealParts  7.ImParts  8.ReRoot  9.ImRoot  10. ProjE"

        IF(tRotateWavefunction) THEN
            OPEN(14,FILE='RotateKilledParts',STATUS='unknown')
            WRITE(14,"(A)") "#Real killed Parts     Imag killed parts"
        ENDIF

        !Start with just real walkers
        WalkListGround(1,1)=InitialWalk
        tFixedShift=.true.
        OldGroundTotParts(1)=InitialWalk
        OldBothTotParts=InitialWalk
        
        do Iter=1,NMCyc

!Every so often, update the shift
            IF((mod(Iter,StepsSft).eq.0).and.(Iter.ne.1)) THEN
                
                CALL CalcProjE(KMat,ProjE,WalkListGround)

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
                WRITE(6,"(I8,4F25.12,2I15,2I7,F25.12)") Iter,GroundShift,RealShift,ImagShift,AverageShift,GroundTotParts(1),GroundTotParts(2),WalkListGround(1,1),WalkListGround(2,1),ProjE
                WRITE(12,"(I8,4F25.12,2I15,2I7,F25.12)") Iter,GroundShift,RealShift,ImagShift,AverageShift,GroundTotParts(1),GroundTotParts(2),WalkListGround(1,1),WalkListGround(2,1),ProjE
                CALL FLUSH(12)
                CALL FLUSH(6)

                IF(tKeepWalkerFiles) THEN

                    Norm=0.D0
                    Norm1=0.D0
                    do i=1,NDet
                        Norm=Norm+(REAL(WalkListGround(1,i),8)**2)+(REAL(WalkListGround(2,i),8)**2)
                        Norm1=Norm1+(REAL(SumWalkListGround(1,i),8)**2)+(REAL(SumWalkListGround(2,i),8)**2)
                    enddo
                    Norm=SQRT(Norm)
                    Norm1=SQRT(Norm1)

                    abstr=''
                    write(abstr,'(I12)') Iter
                    abstr='GroundWavevec-'//adjustl(abstr)
                    OPEN(13,FILE=abstr,STATUS='UNKNOWN')
!                Check=0.D0
                    do i=1,NDet
!                    Check=Check+(REAL(SumWalkListGround(1,i),8)/Norm1)**2+(REAL(SumWalkListGround(2,i),8)/Norm1)**2
                        WRITE(13,"(I8,3F25.12,A,3F25.12)") i,REAL(WalkListGround(1,i),8)/Norm,REAL(WalkListGround(2,i),8)/Norm,SQRT((REAL(WalkListGround(1,i),8)/Norm)**2+(REAL(WalkListGround(2,i),8)/Norm)**2),"   *   ",REAL(SumWalkListGround(1,i),8)/Norm1,REAL(SumWalkListGround(2,i),8)/Norm1,SQRT((REAL(SumWalkListGround(1,i),8)/Norm1)**2+(REAL(SumWalkListGround(2,i),8)/Norm1)**2)
                    enddo
!                WRITE(6,*) "Check = ",Check
                    CLOSE(13)

                ENDIF

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

                            IF((REAL(KMat(i,k),dp).eq.0.D0).and.(AIMAG(KMat(i,k)).eq.0.D0)) CYCLE
                            IF(i.eq.k) CYCLE

                            rh=KMat(k,i)
                            !All spawning attempts here
                            child=Attempt_create(rh,WalkListGround(:,i),part_type)
            
                            !`Child` returns the number of real & im walkers to create at k
                            !create particles.
                            WalkListGroundSpawn(:,k)=WalkListGroundSpawn(:,k)+child

                        enddo

                    enddo

                enddo   !end looping over real & Im

!Attempt to die simultaneously to all particles
                ndie = Attempt_die(KMat,i,WalkListGround(:,i))

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

            !GroundTotParts is now calculated in the RotateWavefunction routine, since this can change the number of 
            !particles in the system
!Rotate the phase of each occupied determinant in the entire wavefunction, to ensure that the root determinant has
!only real walkers on it.
            IF(tRotateWavefunction) THEN
                CALL RotateWavefunction(WalkListGround,GroundTotParts)
            ENDIF

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

    SUBROUTINE RotateWavefunction(WalkListGround,GroundTotParts)
        IMPLICIT NONE
        INTEGER, INTENT(INOUT) :: WalkListGround(lenof_sign,NDet)
        INTEGER, DIMENSION(lenof_sign), INTENT(INOUT) :: GroundTotParts
        INTEGER, DIMENSION(lenof_sign) :: PreRotGroundTotParts
        REAL*8 :: HFAngle,HFMag,rat,r,OrigAngle,Mag,RealPart,ImagPart,NewAngle
        INTEGER :: i

        PreRotGroundTotParts=GroundTotParts
        GroundTotParts=0

        HFMag=SQRT((REAL(WalkListGround(1,1))**2)+(REAL(WalkListGround(2,1))**2))
        HFAngle=ATAN(REAL(WalkListGround(2,1),dp)/REAL(WalkListGround(1,1)))

        IF(WalkListGround(1,1).eq.0) THEN
            !No real particles - kill Im ones and return
            WalkListGround(2,1)=0
            !No phase rotation possible - no real particles at root determinant
            !Still need to count total particles
            do i=1,NDet
                GroundTotParts(:)=GroundTotParts(:)+abs(WalkListGround(:,i))
            enddo

            RETURN
        ELSEIF(WalkListGround(2,1).eq.0) THEN
            !No Im particles - wavefunction phase coherent?
            !Still need to count total particles
            do i=1,NDet
                GroundTotParts(:)=GroundTotParts(:)+abs(WalkListGround(:,i))
            enddo
            RETURN
        ENDIF
        
        !Work out angle information
        !ATAN should return angles -piby2 -> piby2
        !My convention is 0 -> pi angles for positive imaginary amplitudes
        !and 0 -> -pi angles for negative imaginary amplitudes
        IF((WalkListGround(1,1).gt.0).and.(WalkListGround(2,1).gt.0)) THEN
            !First quadrant - phase 0 -> pi/2
            !Should be correct - check this
            IF(HFAngle.gt.piby2.or.HFAngle.lt.0.D0) STOP 'ATAN Error'
        ELSEIF((WalkListGround(1,1).lt.0).and.(WalkListGround(2,1).gt.0)) THEN
            !Second quadrant - phase pi/2 -> pi
            !ATAN should return negative angle
            IF(HFAngle.gt.0.D0) STOP 'ATAN Error 2'
            HFAngle=ABS(HFAngle)+piby2
        ELSEIF((WalkListGround(1,1).lt.0).and.(WalkListGround(2,1).lt.0)) THEN
            !Third quadrant - phase -pi -> -piby2
            !ATAN should return positive angle
            IF(HFAngle.lt.0.D0) STOP 'ATAN Error 3'
            HFAngle=-HFAngle-piby2
        ELSEIF((WalkListGround(1,1).gt.0.D0).and.(WalkListGround(2,1).lt.0)) THEN
            !Fourth quadrant - phase 0 -> -piby2
            !Should be correct - check this
            IF(HFAngle.lt.-piby2.or.HFAngle.gt.0.D0) STOP 'ATAN Error 4'
        ELSE
            STOP 'Error - should not get here'
        ENDIF

        !Rotate the root determinant to the correct phase (0)
        WalkListGround(2,1)=0   !Kill all imaginary particles
        WalkListGround(1,1)=INT(HFMag)   !Ensure magnitude of wavefunction stays as close as possible to the same
        rat=HFMag-WalkListGround(1,1)
        !Stochastically try to create the final particle.
        call random_number(r)
        IF(abs(rat).gt.r) THEN
            WalkListGround(1,1)=WalkListGround(1,1)+nint(sign(1.D0,rat))
        ENDIF

        GroundTotParts(:)=GroundTotParts(:)+abs(WalkListGround(:,1))

        do i=2,NDet

            IF(WalkListGround(1,i).eq.0.and.WalkListGround(2,i).eq.0) CYCLE

            !Now we need to rotate the wavefunction on all other determinants by the same amount.
            !Extreme care needs to be undertaken re. phases.
            OrigAngle=ATAN(REAL(WalkListGround(2,i))/REAL(WalkListGround(1,i)))
            !Work out angle information
            !ATAN should return angles -piby2 -> piby2
            !My convention is 0 -> pi angles for positive imaginary amplitudes
            !and 0 -> -pi angles for negative imaginary amplitudes
            IF((WalkListGround(1,i).gt.0).and.(WalkListGround(2,i).gt.0)) THEN
                !First quadrant - phase 0 -> pi/2
                !Should be correct - check this
                IF(OrigAngle.gt.piby2.or.OrigAngle.lt.0.D0) STOP 'ATAN Error'
            ELSEIF((WalkListGround(1,i).lt.0).and.(WalkListGround(2,i).gt.0)) THEN
                !Second quadrant - phase pi/2 -> pi
                !ATAN should return negative angle
                IF(OrigAngle.gt.0.D0) STOP 'ATAN Error 2'
                OrigAngle=ABS(OrigAngle)+piby2
            ELSEIF((WalkListGround(1,i).lt.0).and.(WalkListGround(2,i).lt.0)) THEN
                !Third quadrant - phase -pi -> -piby2
                !ATAN should return positive angle
                IF(OrigAngle.lt.0.D0) STOP 'ATAN Error 3'
                OrigAngle=-OrigAngle-piby2
            ELSEIF((WalkListGround(1,i).gt.0.D0).and.(WalkListGround(2,i).lt.0)) THEN
                !Fourth quadrant - phase 0 -> -piby2
                !Should be correct - check this
                IF(OrigAngle.lt.-piby2.or.OrigAngle.gt.0.D0) STOP 'ATAN Error 4'
            ELSEIF(WalkListGround(1,i).eq.0) THEN
                !No real particles
                IF(WalkListGround(2,i).gt.0) THEN
                    OrigAngle=piby2
                ELSE
                    OrigAngle=-piby2
                ENDIF
            ELSEIF(WalkListGround(2,i).eq.0) THEN
                !No Imag particles
                IF(WalkListGround(1,i).gt.0) THEN
                    OrigAngle=0.D0
                ELSE
                    OrigAngle=pi
                ENDIF
            ENDIF

            !Original angle is now the correct convention.
            NewAngle=OrigAngle+HFAngle
            IF(NewAngle.gt.pi) then
                NewAngle=NewAngle-pi2
            ELSEIF(NewAngle.lt.-pi) THEN
                NewAngle=NewAngle+pi2
            ENDIF
            
            Mag=SQRT((REAL(WalkListGround(1,i))**2)+(REAL(WalkListGround(2,i))**2))

            !Stochastically realise this magnitude and angle
            RealPart=Mag*COS(NewAngle)
            ImagPart=Mag*SIN(NewAngle)
            WalkListGround(1,i)=INT(RealPart)
            rat=RealPart-INT(RealPart)
            call random_number(r)
            IF(abs(rat).gt.r) THEN
                WalkListGround(1,i)=WalkListGround(1,i)+nint(sign(1.D0,rat))
            ENDIF

            WalkListGround(2,i)=INT(ImagPart)
            rat=ImagPart-INT(ImagPart)
            call random_number(r)
            IF(abs(rat).gt.r) THEN
                WalkListGround(2,i)=WalkListGround(2,i)+nint(sign(1.D0,rat))
            ENDIF

            GroundTotParts(:)=GroundTotParts(:)+abs(WalkListGround(:,i))

        enddo

        WRITE(14,*) GroundTotParts(1)-PreRotGroundTotParts(1),GroundTotParts(2)-PreRotGroundTotParts(2)
            
    END SUBROUTINE RotateWavefunction



    SUBROUTINE CalcProjE(KMat,ProjE,WalkListGround)
        IMPLICIT NONE
        COMPLEX*16, INTENT(IN) :: KMat(NDet,NDet)
        INTEGER, INTENT(IN) :: WalkListGround(lenof_sign,NDet)
        REAL*8 :: ProjE
        complex*16 :: calc_proje
        INTEGER :: i

        calc_proje=cmplx(0.D0,0.D0)

        do i=2,NDet

            calc_proje=calc_proje + KMat(1,i)*cmplx(WalkListGround(1,i),WalkListGround(2,i))

        enddo

        calc_proje=calc_proje/cmplx(WalkListGround(1,1),WalkListGround(2,1))

        if (abs(aimag(calc_proje)) > 1.e-6) write (6,*) 'warning: proje not real!', calc_proje

        proje = real(calc_proje)

    END SUBROUTINE CalcProjE

    FUNCTION Attempt_die(KMat,i,wSign) result(nDie)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: i
        COMPLEX*16, INTENT(IN) :: KMat(NDet,NDet)
        INTEGER, DIMENSION(lenof_sign), INTENT(IN) :: wSign
        INTEGER, DIMENSION(lenof_sign) :: nDie
        REAL*8 :: rat,r,fac
        INTEGER :: part_type
                
        IF(.not.tSeperateShift) THEN
            fac=Tau*(REAL(KMat(i,i),dp)-GroundShift)
        ENDIF

        do part_type=1,lenof_sign

            IF(tSeperateShift) THEN
                IF(part_type.eq.1) THEN
                    !Use real shift
                    fac=Tau*(REAL(KMat(i,i),dp)-RealShift)
                ELSE
                    fac=Tau*(REAL(KMat(i,i),dp)-ImagShift)
                ENDIF
            ENDIF


            rat=fac*abs(wSign(part_type))
            nDie(part_type) = int(rat)
            rat = rat - real(ndie(part_type), dp)
            ! Choose to die or not stochastically
            call random_number(r)
            if (abs(rat) > r) ndie(part_type) = ndie(part_type) + nint(sign(1.0_dp, rat))
            IF(ndie(part_type).gt.abs(wSign(part_type))) STOP 'Trying to create anti-particles'
        enddo

    END FUNCTION Attempt_die

    FUNCTION Attempt_create(rh,parent_sign,part_type) result(child)
        IMPLICIT NONE
        COMPLEX*16, INTENT(IN) :: rh
        REAL*8 :: MatEl,rat,r 
        INTEGER, INTENT(IN) :: part_type
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
        REAL*8 :: Step
        REAL*8 :: r

        KMat(:,:)=CMPLX(0.D0,0.D0)

        Step=(EndEl-StartEl)/REAL(NDet-1,8) !Rate of increase of diagonal matrix elements
        KMat(2,2)=CMPLX(StartEl,0.D0)
        do i=3,NDet
            KMat(i,i)=CMPLX(REAL(KMat(i-1,i-1),8)+Step,0.D0)
        enddo

        WRITE(6,*) "RefDet = ", KMat(1,1)
        WRITE(6,*) "MaxDet = ", KMat(NDet,NDet)

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
                                KMat(j,i)=CMPLX(REAL(KMat(i,j),8),-OffDiagEl)
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
        integer :: i,j,indi,indj
        real*8 :: H00

        ! The first line just contains the size of the Hamiltonian matrix.
        ! This has already been read...
        open(13,file=input_file,status='old',form='formatted')
        READ(13,*) i
        IF(i.ne.NDet) STOP 'NDet not correct'

        do j = 1,NDet
            do i = 1,NDet
                read (13,*) indi,indj,KMat(i,j)
                IF(indi.ne.i.or.indj.ne.j) THEN
                    WRITE(6,*) indi,indj,i,j
                    STOP 'Error reading in hamiltonian'
                ENDIF
            enddo
        end do

        close(13,status='keep')

        do i=1,NDet
            IF(abs(AIMAG(KMat(i,i))).gt.1.D-8) THEN
                STOP 'Diagonal HElements not real'
            ENDIF
            do j=1,NDet
                IF(ABS(AIMAG(KMat(i,j))+AIMAG(KMat(j,i))).gt.1.D-8) THEN
                    STOP 'Matrix not hermitian'
                ENDIF
            enddo
        enddo

        ! K_ij = H_ij - H00 d_ij
        H00 = REAL(KMat(1,1),dp)
        forall (i=1:NDet) KMat(i,i) = CMPLX(REAL(KMat(i,i),dp) - H00,0.D0)

    end subroutine ReadInKMat

    subroutine dumpkmat(KMat)
        IMPLICIT NONE
        COMPLEX*16 KMat(NDet,NDet)
        INTEGER :: i,j

        OPEN(19,FILE='DumpedKMat',status='unknown')

        do i=1,NDet
            do j=1,NDet
                WRITE(19,"(2G25.15)") REAL(KMat(i,j),dp),AIMAG(KMat(i,j))
            enddo
        enddo

        CLOSE(19)

    end subroutine dumpkmat

END Program ModelFCIQMC
