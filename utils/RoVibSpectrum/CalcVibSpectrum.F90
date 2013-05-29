PROGRAM CalcVibSpectrum
    IMPLICIT NONE
    REAL*8, PARAMETER :: PI=3.1415926535897932384626433832795029D0
    REAL*8 :: Atom1Mass,Atom2Mass,Rin,Rout,RedMass,Rspacing,Requilib,Rcentre,SumTemp
    INTEGER :: Rpoints,JRotValue,JRotValueMax,ierr,i,j,k,kmax,WorkSize,WorkCheck,v,Maxi,x,y,MaxVibPlot,MaxvForFit
    REAL*8 , ALLOCATABLE :: Hamil(:,:),akvalues(:),EvectorsInit(:),Evectors(:),Work(:),VeqPlusGv(:)
    REAL*8 :: alpha,beta,Const,VibPlotMinR,VibPlotMaxR,Veq,TwiceAtomEnergy,DissV
    LOGICAL :: exists

! First open the 'inputVib' file which has some input parameters.    
        INQUIRE(FILE='inputVib',EXIST=exists)
        IF(.not.exists) THEN
            WRITE(6,*) 'No input file (name inputVib) detected.'
            STOP
        ELSE
            OPEN(17,FILE='inputVib',Status='old')
        ENDIF

        READ(17,*) MaxvForFit           !The maximum v to be used for the fitting to get rotational and vibrational parameters.
        READ(17,*) JRotValueMax         !Maximum J value you want to look at.
        READ(17,*) Atom1Mass            !Masses in amu.
        READ(17,*) Atom2Mass
        READ(17,*) Rin                  !Starting R in bohr.
        READ(17,*) Rout                 !Ending R in bohr
        READ(17,*) Rpoints              !Number of points in grid.
        READ(17,*) MaxVibPlot           !Maximum vibrational level we want to plot.
        READ(17,*) VibPlotMinR,VibPlotMaxR  !The R values between which we want to plot them (and the potential).
        CLOSE(17)

! Second input file contains the data on the potential.        
        INQUIRE(FILE='PotData',EXIST=exists)
        IF(.not.exists) THEN
            WRITE(6,*) 'No potential data file (name PotData) detected.'
            STOP
        ELSE
            OPEN(17,FILE='PotData',Status='old')
        ENDIF
        READ(17,*) kmax
        READ(17,*) Requilib
        READ(17,*) alpha
        READ(17,*) beta
        READ(17,*) TwiceAtomEnergy
        ALLOCATE(akvalues(0:kmax),stat=ierr)
        do i=0,kmax
            READ(17,*) akvalues(i)
        enddo
        CLOSE(17)

! Calculate the reduced mass (in terms of *m_e*!!!!).        
        Const=(1.66053886E-027)/(0.91093826E-30)
        RedMass=((Atom1Mass*Atom2Mass)/(Atom1Mass+Atom2Mass))*Const
! Reduced mass of the diatomic. 

! Rspacing is deltaR, the spacing of the grid for the kinetic energy.
        Rspacing=(Rout-Rin)/(REAL(Rpoints))
! Need the grid to go from -infinity to +infinity, so the grid points start at -(Number of points)/2.        
        Maxi=Rpoints/2
! Rcentre is the centre of the grid.        
        Rcentre=(Rout+Rin)/2.D0

! Allocate matrix for Hamiltonian.
        ALLOCATE(Hamil(1:(Rpoints-1),1:(Rpoints+1)),stat=ierr)

! The EvectorsInit that come out are E_v,J = V_eq + G_v + F_v(J)
! For J=0, F_v(0) = 0 and G_v is the pure vibrational spectrum.
        ALLOCATE(EvectorsInit(1:(Rpoints-1)),stat=ierr)
! For J=0 then, the evectors are originally the V_eq + G_v part.
! Can get the F_v(J) part later (for higher values of J) by taking this away from the eigenvalues obtained with larger J.
        ALLOCATE(VeqPlusGv(1:(Rpoints-1)),stat=ierr)
        VeqPlusGv(:)=0.D0
        
! For J=0, Evectors is G_v, for higher J, Evectors is F_v(J).         
        ALLOCATE(Evectors(1:(Rpoints-1)),stat=ierr)

        WorkCheck=3*(RPoints-1)+1
        WorkSize=WorkCheck
        ALLOCATE(Work(WorkSize),stat=ierr)

! This is the file containing the data at different J's.  
! This is later used to get out the rotational constants etc.
! Not quite sure what we want to be printing at this stage.
        IF(JRotValueMax.gt.0) THEN
            OPEN(17,FILE="ROTDATA",status='unknown')
            OPEN(22,FILE="ROTVIBSPECTRUM",status='unknown')
            WRITE(17,'(A1,A5)',advance='no') '#','J'
            WRITE(22,'(A)') '# These are the theoretical values for what would be experimentally measured: E_v,J - E_0,0 = (G_v - G_0) + F_v(J)'
            WRITE(22,'(A)') '# i.e. the pure vib spectrum (J=0) with the rotational energy levels also included (all in cm-1)'
            WRITE(22,'(A1,A5)',advance='no') '#','J'
            do v=1,(MaxvForFit+1)
                IF(v.eq.(MaxvForFit+1)) THEN
                    WRITE(17,'(A18,I2)') 'F_v(J) v =',(v-1)
                    WRITE(22,'(A18,I2)') 'E_v,J;  v =',(v-1)
                ELSE
                    WRITE(17,'(A18,I2)',advance='no') 'F_v(J) v =',(v-1)
                    WRITE(22,'(A18,I2)',advance='no') 'E_v,J;  v =',(v-1)
                ENDIF
            enddo
        ENDIF

        do JRotValue=0,JRotValueMax
            Hamil(:,:)=0.D0
            EvectorsInit(:)=0.D0
            Evectors(:)=0.D0
            IF((JRotValueMax.gt.0).and.(JRotValue.ne.0)) THEN
                WRITE(17,'(I6)',advance='no') JRotValue
                WRITE(22,'(I6)',advance='no') JRotValue
            ELSEIF(JRotValueMax.gt.0) THEN
                WRITE(22,'(I6)',advance='no') JRotValue
            ENDIF

! x and y are the positions in the matrix, but i and j are the positions on the grid (go from -infinity -> infinity), but in 
! reality they go from Rin -> Rout.
            j=(-1)*Maxi
            do y=1,(Rpoints-1)
                i=(-1)*Maxi
                do x=1,(Rpoints-1)

                    IF(i.eq.j) THEN
                        ! First add in the kinetic energy term for when i=j.
                        Hamil(x,y)=(1.D0/(2.D0*RedMass*(Rspacing**2.D0)))*((-1.D0)**(i-j))*((PI**2.D0)/3.D0)

                        ! Then add in the J term (this is 0 when J=0).
                        Hamil(x,y)=Hamil(x,y) + (REAL(JRotValue*(JRotValue+1))/((((REAL(i)*Rspacing)+Rcentre)**2.D0)*(2.D0*RedMass)))    
                        
                        ! Then add in the value of the potential which is diagonal in this real space grid.
                        do k=0,kmax
                            Hamil(x,y)=Hamil(x,y)+(akvalues(k)*(EXP((-1.D0)*alpha*(beta**k)*(((REAL(i)*Rspacing)+Rcentre)**2.D0))))
                        enddo

                    ELSE

                        ! The off diagonal terms only include the kinetic energy terms.
                        Hamil(x,y)=(1.D0/(2.D0*RedMass*(Rspacing**2.D0)))*((-1.D0)**(i-j))*(2.D0/((i-j)**2.D0))

                    ENDIF
                    i=i+1
                enddo
                j=j+1
            enddo

!            WRITE(6,*) 'Hamil'
!            do j=1,Rpoints+1
!                do i=1,Rpoints+1
!                    WRITE(6,*) Hamil(i,j),Hamil(j,i)
!                enddo
!            enddo
!            CALL FLUSH(6)

! Diagonalise the hamiltonian to get out the eigenvalues and the eigenvectors.
            CALL DSYEV('V','U',(Rpoints-1),Hamil(1:(Rpoints-1),1:(Rpoints-1)),(Rpoints-1),EvectorsInit(1:(Rpoints-1)),Work,WorkSize,ierr)
            ! Hamil goes in as the original Hamil, comes out as the eigenvectors (Coefficients).
            ! TMAT2DBlock comes out as the eigenvalues in ascending order.


            Veq=0.D0
            do k=0,kmax
                Veq=Veq+(akvalues(k)*(EXP((-1.D0)*alpha*(beta**k)*(Requilib**2.D0))))
            enddo

            DissV=0.D0
            do k=0,kmax
                DissV=DissV+(akvalues(k)*(EXP((-1.D0)*alpha*(beta**k)*((5.D0*Requilib)**2.D0))))
            enddo

            OPEN(18,FILE="SPEC_CONSTS",status='unknown')
            WRITE(18,'(A)') "SPECTROSCOPIC PROPERTIES"
            WRITE(18,'(A)') "Note: These are formatted for comparison to papers such as D.Feller and J.A.Sordo,J.Chem.Phys.V113.2000"
            WRITE(18,*) ''
            WRITE(18,'(A,F15.6,A)') 'E_total (in E_h) : ',Veq+TwiceAtomEnergy,'     The energy of the diatomic at the eqm bond length.'
            WRITE(18,'(A,F10.2,A)') 'D_e (in kcal/mol) : ',(Veq-DissV)*(-1.D0)*627.509,'     The dissociation energy.'
            WRITE(18,'(A,F10.5,A)') 'r_e (in Angstroms) : ',Requilib*0.529177249,'     The equilibrium bond length.'
            CLOSE(18)


! Now we just basically print the stats we want.            
            IF(JRotValue.eq.0) THEN

                ! We have the pure vibrational spectrum.
                do i=1,(Rpoints-1)
                    VeqPlusGv(i)=EvectorsInit(i)
                enddo
                do i=1,(Rpoints-1)
                    Evectors(i)=EvectorsInit(i)
                    Evectors(i)=Evectors(i)-Veq
                enddo
                OPEN(18,FILE="PUREVIB",status='unknown')
                WRITE(18,'(A1,A6,3A30)') '#',' v','G_v (E_h)','G_v - G_0 (E_h)','G_v - G_0 (cm-1)'
                do v=1,(Rpoints-1)
! Only print points that have vibrational frequencies less than the depth of the well (well the depth of the ZPE since that where we vibrationally excite from).
                    IF(((Evectors(v)*219474.63)-(Evectors(1)*219474.63)).lt.(ABS(Veq*219474.63)-(Evectors(1)*219474.63))) THEN
                        WRITE(18,'(I7,3F30.10)') (v-1),Evectors(v),(Evectors(v)-Evectors(1)),((Evectors(v)*219474.63)-(Evectors(1)*219474.63))
                    ENDIF
                enddo
                WRITE(18,'(A1,A6,A30,2F30.10)') '#',' ZPE','E_h : cm-1',Evectors(1),Evectors(1)*219474.63
                WRITE(18,'(A1,A6,A30,2F30.10)') '#',' D_0','E_h : cm-1',ABS(Veq)-Evectors(1),((ABS(Veq)*219474.63)-(Evectors(1)*219474.63))
                WRITE(18,'(A1,A6,A30,2F30.10)') '#',' D_e','E_h : cm-1',ABS(Veq),ABS(Veq)*219474.63
                WRITE(18,'(A1,A6,A30,F30.10)') '#',' R_e','bohr',Requilib
                CLOSE(18)

                OPEN(18,FILE="POT.WAVFNS",status='unknown')
                WRITE(18,'(A1,A4,A20,A30,A20)') '#','v','R (bohr)','V(R) (mE_h)','Psi(R)'
                i=(-1)*Maxi
                do v=1,(Rpoints-1)
                    IF((((REAL(i)*Rspacing)+Rcentre).ge.VibPlotMinR).and.(((REAL(i)*Rspacing)+Rcentre).le.VibPlotMaxR)) THEN
                        SumTemp=0.D0
                        do k=0,kmax
                            SumTemp=SumTemp+(akvalues(k)*(EXP((-1.D0)*alpha*(beta**k)*(((REAL(i)*Rspacing)+Rcentre)**2.D0))))
                        enddo
                        WRITE(18,'(I5,F20.10,F30.10)',advance='no') v,(REAL(i)*Rspacing)+Rcentre,(SumTemp*1E3)
                        do j=1,MaxVibPlot+1 
!                            WRITE(18,'(F20.10)',advance='no') ((Hamil(v,j)/((REAL(i)*Rspacing)+Rcentre))+(Evectors(v)+Veq))*1E3)
                            IF(j.eq.(MaxVibPlot+1)) THEN
                                WRITE(18,'(F20.10)') (((Hamil(v,j)/((REAL(i)*Rspacing)+Rcentre))*10)+(Veq*1E3)+(Evectors(j)*1E3))
                            ELSE
                                WRITE(18,'(F20.10)',advance='no') (((Hamil(v,j)/((REAL(i)*Rspacing)+Rcentre))*10)+(Veq*1E3)+(Evectors(j)*1E3))
                            ENDIF
                        enddo
                    ENDIF
                    i=i+1
                enddo
                CLOSE(18)

               
                OPEN(18,FILE="GnuplotVib.gpi",status='unknown')
                WRITE(18,*) 'set ylabel "mE_h"'
                WRITE(18,*) 'set xlabel "R / bohr"'
                WRITE(18,'(A4)',advance='no') "plot"
                do j=1,MaxVibPlot+1
                    WRITE(18,'(F20.10,A25)') (Evectors(j)*(1E3))+(Veq*1E3)," lc rgb 'blue' notitle, \"
                enddo
                WRITE(18,'(A48)') "'POT.WAVFNS' u 2:3 w l lc rgb 'black' notitle, \"
                do j=1,MaxVibPlot+1
                    IF(j.eq.(MaxVibPlot+1)) THEN
                        WRITE(18,'(A17,I2,A25)') "'POT.WAVFNS' u 2:",(j+3)," w l lc rgb 'red' notitle"
                    ELSE
                        IF(j.lt.7) THEN
                            WRITE(18,'(A17,I1,A28)') "'POT.WAVFNS' u 2:",(j+3)," w l lc rgb 'red' notitle, \"
                        ELSE
                            WRITE(18,'(A17,I2,A28)') "'POT.WAVFNS' u 2:",(j+3)," w l lc rgb 'red' notitle, \"
                        ENDIF
                    ENDIF
                enddo
                CLOSE(18)

                IF(JRotValueMax.gt.0) THEN
                    ! Write out value of F_v(J) which for J=0 is 0.
!                    do v=1,(MaxvForFit+1)
!                        IF(v.eq.(MaxvForFit+1)) THEN
!                            WRITE(17,'(F20.10)') 0.D0 
!                        ELSE
!                            WRITE(17,'(F20.10)',advance='no') 0.D0 
!                        ENDIF
!                    enddo
                    do v=1,MaxvForFit+1
                        IF(v.eq.(MaxvForFit+1)) THEN
                            WRITE(22,'(F20.10)') ((Evectors(v)*219474.63)-(Evectors(1)*219474.63))
                        ELSE
                            WRITE(22,'(F20.10)',advance='no') ((Evectors(v)*219474.63)-(Evectors(1)*219474.63))
                        ENDIF
                    enddo
                ENDIF

            ELSE
                do v=1,(Rpoints-1)
                    Evectors(v)=EvectorsInit(v)-VeqPlusGv(v)
                enddo

                do v=1,MaxvForFit+1
                    IF(v.eq.(MaxvForFit+1)) THEN
                        WRITE(17,'(F20.10)') Evectors(v)/(JRotValue*(JRotValue+1))
                        WRITE(22,'(F20.10)') (EvectorsInit(v)-VeqPlusGv(1))*219474.63
                    ELSE
                        WRITE(17,'(F20.10)',advance='no') Evectors(v)/(JRotValue*(JRotValue+1))
                        WRITE(22,'(F20.10)',advance='no') (EvectorsInit(v)-VeqPlusGv(1))*219474.63
                    ENDIF
                enddo

                do v=1,MaxVibPlot+1
                    IF(v.eq.(MaxVibPlot+1)) THEN
                        WRITE(22,'(F20.10)') (EvectorsInit(v)-VeqPlusGv(1))*219474.63
                    ELSE
                        WRITE(22,'(F20.10)',advance='no') (EvectorsInit(v)-VeqPlusGv(1))*219474.63
                    ENDIF
                enddo
            ENDIF
        enddo

        OPEN(18,FILE="GnuplotRotVib.gpi",status='unknown')
        WRITE(18,*) 'set xlabel "E_v,J - E_0,0 / cm-1"'
        WRITE(18,'(A14)',advance='no') "plot [][0:1.5]"
        do j=1,MaxVibPlot+1
            IF(j.eq.(MaxVibPlot+1)) THEN
                WRITE(18,'(A19,I2,A39)') "'ROTVIBSPECTRUM' u ",(j+1),":(($1)/($1)) w i lc rgb 'black' notitle"
            ELSE
                IF(j.lt.9) THEN
                    WRITE(18,'(A19,I1,A42)') "'ROTVIBSPECTRUM' u ",(j+1),":(($1)/($1)) w i lc rgb 'black' notitle, \"
                ELSE
                    WRITE(18,'(A19,I2,A42)') "'ROTVIBSPECTRUM' u ",(j+1),":(($1)/($1)) w i lc rgb 'black' notitle, \"
                ENDIF
            ENDIF
        enddo
        CLOSE(18)



        IF(JRotValueMax.gt.0) THEN
            CLOSE(17)
            CLOSE(22)
        ENDIF
        DEALLOCATE(Hamil)
        DEALLOCATE(EvectorsInit)
        DEALLOCATE(Evectors)
        DEALLOCATE(VeqPlusGv)
        DEALLOCATE(Work)


    END PROGRAM
