! This function is based on the principle that the correlation energy of a system is simply the sum of
! the hamiltonian elements between HF and each excitation, multiplied by the probability of a random
! walker populating the determinant in the long time limit, normalised by the probability of the random
! walker populating the HF determinant. Because of this, all graphs can be calculated from a discrete 
! diffusion problem. In this case, the principle is used to calculate the energy of the star graph, by 
! picking a determinant of the graph at random (including HF), and then pushing electron density onto
! connected determinants according to the rho matrix element connecting them. This can be consdered a
! local application of the rho matrix a random determinant, and so the correct energy should be eventually
! converged upon.

MODULE MCStarMod
    USE HElem
    USE MemoryManager , only : LogMemAlloc,LogMemDealloc 
    USE System , only : NEl
    USE Determinants , only : FDet
    IMPLICIT NONE
    SAVE

!Array to hold the excitation info for the determinants
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ
    TYPE(HElement) , ALLOCATABLE :: ExcitInfo(:,:)
    INTEGER :: ExcitInfoTag=0

    TYPE(HElement) , ALLOCATABLE :: Eigenvector(:)
    INTEGER :: EigenvectorTag=0

!NoExcits is the total number of connected excitations in the full star graph
    INTEGER :: Seed,NoExcits

    REAL*8 :: NormFactor,HFReNorm

    TYPE(HElement) :: rhii,Hii

!The intermediate weights and energys of the star graph, calculated with the updated wavevector are stored
    TYPE(HElement) :: TempEnergyxw,TempWeight
      
    contains

    SUBROUTINE MCStar(Weight,Energyxw)
        USE System, only: Alat,Beta,Brr,ECore,G1,nBasis,nBasisMax,Arr
        USE Calc , only : i_P,G_VMC_Seed
        USE Integrals, only : fck,nMax,nMsh,UMat,nTay
        USE Determinants , only : GetHElement2
        IMPLICIT NONE
        TYPE(HDElement) :: Weight,Energyxw
        CHARACTER(len=*), PARAMETER :: this_routine='MCStar'
        
        OPEN(63,file='MCStarStats',Status='unknown')
        IF(HElementSize.ne.1) STOP 'Only real orbitals allowed in MCStar so far'
        Weight=HDElement(0.D0)
        Energyxw=HDElement(0.D0)

!Initialise random number generator
        Seed=G_VMC_Seed

!Find rho_ii value, which all rho elements will be divided by, and Hii value
        CALL CalcRho2(FDet,FDet,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rhii,nTay,0,ECore)
        Hii=GetHElement2(FDet,FDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)

!First fill excitinfo
        CALL FindStarExcits()

!Then create initial trial normalised wavevector
        CALL CreateInitTrialWavevector()

!Routine to pick determinants at random (inc. HF) according to probability given by eigenvector component.
!The wavevector is then propagated in all possible allowed directions, and the energy updated
        CALL PropagateLocalWavevector()

        CLOSE(63)

!Deallocate info...
        DEALLOCATE(Eigenvector)
        CALL LogMemDealloc(this_routine,EigenvectorTag)
        DEALLOCATE(ExcitInfo)
        CALL LogMemDealloc(this_routine,ExcitInfoTag)

!Return final info
        Weight=HDElement((TempWeight%v)-1.D0)
        Energyxw=HDElement(TempEnergyxw%v)

        RETURN
    END SUBROUTINE MCStar

!This routine finds all Hij, rho_ij and rho_jj elements for all double excitations in the star graph, and fills ExcitInfo with this
    SUBROUTINE FindStarExcits()
        USE System , only : G1,Alat,Beta,Brr,ECore,nBasis,nBasisMax,Arr
        USE Calc , only : i_P,RhoEps,dBeta
        USE Integrals , only : fck,nMax,nMsh,UMat,nTay
        USE Determinants , only : GetHElement2
        USE Logging , only : iLogging
        IMPLICIT NONE
        INTEGER :: nStore(6),exFlag,nExcitMemLen,iMaxExcit,nJ(NEl)
        INTEGER , ALLOCATABLE :: nExcit(:)
        INTEGER :: nExcitTag=0
        TYPE(HElement) :: rh
        REAL*8 :: StarWeight,DLWDB
        INTEGER :: iExcit,i,j,k,iSubFindStar,nRoots
        CHARACTER(len=*), PARAMETER :: this_routine='FindStarExcits'
        INTEGER :: ierr,ExcitCount
        LOGICAL :: TCountExcits 

        CALL TISET('FindStarExcits',iSubFindStar)

!HFReNorm is equal to the increase in total probability when the propagation step is from the HF determinant
!This is simply equal to the increase electron probability on each determinant, i.e. 1 (rhii/rhii) + sum_j (|rhij|/rhii)^x
        HFReNorm=1.D0

!        TCountExcits=BTEST(nWHTay,8)
!.. Allow both singles and doubles
        exFlag=3

!.. Count the excitations. - First call of GenSymExcitIt2 calculates memory needed for internal use in excitation generators
        nStore(1)=0
        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
        Allocate(nExcit(nExcitMemLen))

!Second call calculates size of arrays needed to store all symmetry allowed excitations - further calls will generate excitation on-the-fly(shown by the false in arg(6)
        nExcit(1)=0
        CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)

! iMaxExcit now contains the number of excitations.
!TCountExcits will run through all excitations possible, determine if they are connected, and then only store these.
!Will be twice as expensive, as needs to run through all excitations twice - however, will only store memory needed.
        
!        IF(TCountExcits) THEN
!            Write(6,"(A,I10,A)") "Counting excitations - Running through all ",iMaxExcit," excitations to determine number connected"
!            ExcitCount=0
!            do while(.true.)
!                CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,exFlag)
!                IF(nJ(1).eq.0) exit
!                CALL CalcRho2(FDet,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,iExcit,ECore)
!                IF(rh.agt.RhoEps) ExcitCount=ExcitCount+1
!            enddo
!
!!Set number of excitations to number of connected determinants, and reset generator
!            Deallocate(nExcit)
!            nStore(1)=0
!            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
!            Allocate(nExcit(nExcitMemLen))
!            nExcit(1)=0
!            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
!
!!set number of excitations to exact number there are
!            iMaxExcit=ExcitCount
!                        
!        ENDIF
        
!.. Allocate memory for the lists
        Write(6,*) "Allocating storage for ",iMaxExcit," excitations."
        Allocate(ExcitInfo(0:iMaxExcit,0:2),stat=iErr)
        CALL LogMemAlloc('ExcitInfo',(iMaxExcit+1)*3,HElementSize,this_routine,ExcitInfoTag)
        CALL AZZERO(ExcitInfo,(iMaxExcit+1)*3*HElementSize)

        i=0
        ExcitInfo(i,0)=1.D0
        ExcitInfo(i,1)=1.D0
        ExcitInfo(i,2)=Hii
                    
        do while(.true.)
            CALL GenSymExcitIt2(FDet,NEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,exFlag)
            IF(nJ(1).eq.0) exit

!Calculate and store rhoij element
            CALL CalcRho2(FDet,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,iExcit,ECore)
            
            if(rh .agt. RhoEps) then
               i=i+1
!Divide all elements though by rhoii
               ExcitInfo(i,1)=rh/rhii

!Calculate a value which will be equal to the change in total probability if the HF determinant is chosen in a propagation step
                HFReNorm=HFReNorm+((ABS(ExcitInfo(i,1)%v))**2)
               
!Calculate and store rho_jj elements
               CALL CalcRho2(nJ,nJ,Beta,i_P,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,Arr,ALat,UMat,rh,nTay,0,ECore)
               ExcitInfo(i,0)=rh/rhii
               ExcitInfo(i,2)=GetHElement2(FDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,iExcit,ECore)

            endif
        enddo

!The total number of excitations is now put into NoExcits
        NoExcits=i
        Deallocate(nExcit)

!If we want, we should be able to determine the value from polynomial diagonalisation of the full star matrix to compare
        WRITE(6,*) "Calculating highest eigenvector of star by polynomial diagonalisation for comparison..."
        nRoots=1
        StarWeight=0.D0
        DLWDB=0.D0
        CALL StarDiag2(0,NEl,NoExcits+1,ExcitInfo,iMaxExcit+1,Beta,i_P,StarWeight,dBeta(1),DLWDB,nRoots,iLogging)
        DLWDB=DLWDB+Hii%v
        StarWeight=StarWeight+1.D0
        WRITE(6,*) "Energy of Star Graph calculated to be: ", DLWDB/StarWeight

        IF(ierr.ne.0) STOP 'Problem in allocation somewhere in FindStarExcits'
        CALL TIHALT('FindStarExcits',iSubFindStar)

    END SUBROUTINE FindStarExcits

!This routine is designed to calculate the initial attempt at an eigenvector to the star matrix problem.
!This needs to be normalised (or at least sum of probabilities known)
!It is desireable that the initial wavevector will have zero components for determinants which are not connected to HF
!This will ensure that determinants which are not connected do not try and propagate any wavefunctio from them
!They would not succeed (we put rhjj=0, and rhij=0 by definition), though we could not get an exact result since there would
!be no way to reduce the eigenvector component which is initially given to them
    SUBROUTINE CreateInitTrialWavevector()
        IMPLICIT NONE
        CHARACTER(len=*), PARAMETER :: this_routine='CreateInitTrial'
        INTEGER :: iErr,i
        REAL*8 :: StarEnergy,NormCheck
        TYPE(HElement) :: NormCons

!We first need to allocate memory to hold the trial wavevector
        Allocate(Eigenvector(0:NoExcits),stat=iErr)
        IF(iErr.ne.0) STOP 'Problem in allocation of Eigenvector'
        CALL LogMemAlloc('Eigenvector',NoExcits+1,HElementSize,this_routine,EigenvectorTag)
        CALL AZZERO(Eigenvector,(NoExcits+1)*HElementSize)

!Since NoExcits counts only connected determinants, we only need an eigenvector this length
!For an initial guess, let the components be equal to their rhoij/rhii elements (can be negative)
!For this specific case, we already know the normalisation constant from HFReNorm
        NormCons=HElement(SQRT(HFReNorm))
        Eigenvector(0)=HElement(1.D0)/NormCons
        NormCheck=(Eigenvector(0)%v)**2
        do i=1,NoExcits
            Eigenvector(i)=ExcitInfo(i,1)/NormCons
            NormCheck=NormCheck+(Eigenvector(i)%v)**2
        enddo
        IF((ABS(NormCheck-1.D0)).gt.1.D-08) THEN
            WRITE(6,*) "Initial trial wavevector not correctly normalised"
            WRITE(6,*) NormCheck
            STOP "Initial trial wavevector not correctly normalised"
        ENDIF

!Calculate initial energy of star with trial wavevector
        TempEnergyxw=HElement(0.D0)
        do i=1,NoExcits
            TempEnergyxw=TempEnergyxw+(Eigenvector(i)*ExcitInfo(i,2))
        enddo
        TempWeight=Eigenvector(0)
        StarEnergy=(TempEnergyxw%v/TempWeight%v)+Hii%v

!Set the NormFactor to be the initial normalisation factor of the trial wavevector
!This is simply equal to the sum of the probabilities, or the sum of the squares of the trial eigenvector
        NormFactor=1.D0

        WRITE(63,"(I15,2G22.14)") 1,StarEnergy,NormFactor

    END SUBROUTINE CreateInitTrialWavevector


!In this routine, determinants are picked stochastically according to the magnitude of their component of 
!the wavevector at that time squared. If the root is picked, then the weight of all excitations are increased
!according to the size of their rij/rhii values, and the root component is increased by 1. Normalisation constants
!are also increased accordingly. If an excitation is chosen, then the root increases by rij/rhii, and the
!excitation chosen is increased by rhjj/rhii. Again the results affect normalisation.
    SUBROUTINE PropagateLocalWavevector()
        USE Calc , only : Iters
        IMPLICIT NONE
        INTEGER :: Iterations,iSubProp,i,j
        REAL*8 :: r,RAN2,StarEnergy,OrigRoot,OrigExcit
        CHARACTER(len=*), PARAMETER :: this_routine='PropLocalWaveVec'
        
        CALL TISET('PropLocalWaveVec',iSubProp)

!Cycle over the number of iterations of local application of rho to the graph we want
        do Iterations=1,Iters

!Multiply by NormFactor - the sum of the probabilities
!This renormalises the wavevector each time
            r=RAN2(Seed)*NormFactor

!Set i=-1, to allow all excitation + root
            i=-1

!We first need to choose an excitation or root according to the renormalised trial wavevector
            do while ((r.gt.0.D0).and.(i.lt.NoExcits))
                i=i+1

                r=r-((Eigenvector(i)%v)**2)
            enddo

            IF(r.gt.0.D0) THEN
!Error in normalisation of wavevector
                WRITE(6,*) "Error in normalisation of trial wavevector - exiting..."
                STOP 'Error in normalisation of trial wavevector'
            ENDIF

            IF(i.eq.0) THEN
!Root is selected to propagate from...
                WRITE(6,*) "ROOT PICKED"
!Set normfactor back to zero - all values change, so we have to create a full new normfactor, and energy factors
                TempEnergyxw=HElement(0.D0)
                TempWeight=HElement(0.D0)
                NormFactor=0.D0
                
                do j=1,NoExcits
!Add to the wavevector components of the excitations rhoij/rhii
!However, since the determinant was chosen with probability given by the component squared, 
!one of these factors needs to be divided out to achieve the correct weighting of each determiant
                    Eigenvector(j)=Eigenvector(j)+ExcitInfo(j,1)/Eigenvector(0)
                    NormFactor=NormFactor+((Eigenvector(j)%v)**2)
!Calculate energy again...
                    TempEnergyxw=TempEnergyxw+HElement(ExcitInfo(j,2)%v*Eigenvector(j)%v)
                enddo

!The root also increases its wavevector component by rhii/rhii, and is divided through by the same
!factor for the same reason
                Eigenvector(0)=Eigenvector(0)+HElement(1.D0)/Eigenvector(0)
                NormFactor=NormFactor+Eigenvector(0)%v

                TempWeight=Eigenvector(0)

            ELSE
!Excitation i is selected to propagate from...
                WRITE(6,*) "EXCIT PICKED"
                
!The root is increased by a proportion given by rhij/rhii, divided by the component of the eigenvector chosen
!The original values of the eigenvectors are needed to avoid renormalisation
                OrigRoot=Eigenvector(0)%v
                Eigenvector(0)=Eigenvector(0)+ExcitInfo(i,1)/Eigenvector(i)

!The excitation itself is increased by the diagonal element of the rho-matrix, again divided by the cpt. of eigenvector
                OrigExcit=Eigenvector(i)%v
                Eigenvector(i)=Eigenvector(i)+ExcitInfo(i,0)/Eigenvector(i)

!The full normalisation does not need to be calculated again
                NormFactor=((Eigenvector(0)%v)**2)+((Eigenvector(i)%v)**2)
                NormFactor=NormFactor+(1.D0-(OrigRoot**2)-(OrigExcit**2))

!The same trick can be used to calculate the desired terms in the energy
                TempEnergyxw=TempEnergyxw-(ExcitInfo(i,2)*(HElement(OrigRoot)-Eigenvector(i)))
                TempWeight=Eigenvector(0)

            ENDIF
                
            StarEnergy=(TempEnergyxw%v/TempWeight%v)+Hii%v
            WRITE(63,"(I15,2G22.14)") 1,StarEnergy,NormFactor

        enddo
        
        CALL TIHALT('PropLocalWaveVec',iSubProp)

    END SUBROUTINE PropagateLocalWavevector

END MODULE MCStarMod
