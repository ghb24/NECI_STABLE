!This code takes a trial CI vector/Graph, and then systematically improves it.
!The aim is to apply (exp[-Delta*H])^P, which is expanded out as iterative application of
!Psi(t+1)=Psi(t)-Delta*HPsi(t). This H is applied locally to every determinant in the graph.

MODULE GraphMorph2

    USE System , only : NEl
    USE Determinants , only : FDet
!Iters is the number of interations of morphing the graph, DetsMax is the maximum size of the CI vector
    USE Calc , only : Iters,DetsMax,DeltaH,TStoch
    USE MemoryManager , only : LogMemAlloc,LogMemDealloc
    USE HElem

    IMPLICIT NONE
    SAVE

!Array to hold CI vector
    REAL*8 , ALLOCATABLE :: CIVect(:)
    INTEGER :: CIVectTag=0

!Array to hold configuration of each determinant in the CIvect
!The zeroth element of each determinant holds the excitation level from the HF determinant
    INTEGER , ALLOCATABLE :: CIConfig(:)
    INTEGER :: CIConfigTag=0

    INTEGER :: Seed

!Dets is the current number of excited determinants in the CIvector
    INTEGER :: Dets

!TempEnergyxw and TempWeight are the weight and energyxw of the current CI vector
    REAL*8 :: TempEnergyxw,TempWeight

    REAL*8 :: Hii

    contains

    SUBROUTINE MorphGraph2(Weight,Energyxw)
        USE System, only: Alat,Beta,Brr,ECore,G1,nBasis,nBasisMax
        USE Calc , only : G_VMC_Seed
        USE Integrals, only : fck,nMax,nMsh,UMat
        USE Determinants , only : GetHElement2
        IMPLICIT NONE
        TYPE(HDElement) :: Weight,Energyxw
        TYPE(HElement) :: Helii
        INTEGER :: ierr,i,Iteration
        CHARACTER(len=*), PARAMETER :: this_routine='MorphGraph2'
    
        OPEN(64,file='MCMorphStats',Status='unknown')
        
        IF(HElementSize.ne.1) STOP 'Only real orbitals allowed in GraphMorph2 so far'

        Weight=HDElement(0.D0)
        Energyxw=HDElement(0.D0)
    
!Initialise random number generator
        Seed=G_VMC_Seed
        Helii=GetHElement2(FDet,FDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)
        Hii=Helii%v

!Allocate memory to hold CI vector
        ALLOCATE(CIConfig(0:DetsMax,0:NEl),stat=ierr)
        CALL LogMemAlloc('CIConfig',(DetsMax+1)*(NEl+1),4,this_routine,CIConfigTag)
        CALL IAZZERO(CIConfig,(DetsMax+1)*(NEl+1))
        ALLOCATE(CIVect(0:DetsMax),stat=ierr)
        CALL LogMemAlloc('CIVect',DetsMax+1,8,this_routine,CIVectTag)
        CALL AZZERO(CIVect,DetsMax+1)
        IF(ierr.ne.0) STOP 'Problem in allocating memory for GraphMorph2'

!Create initial normalised trial CI vector
        CALL CreateInitTrialWavevect()

        do Iteration=1,Iters
            
            CALL ApplyHMat()

        enddo

        RETURN
    END SUBROUTINE MorphGraph2

    SUBROUTINE ApplyHMat()
        IMPLICIT NONE
        CHARACTER(len=*), PARAMETER :: 

    SUBROUTINE CreateInitTrialWavevect()
        IMPLICIT NONE
        INTEGER :: iErr,i
        REAL*8 :: StarEnergy,NormCheck
        TYPE(HElement) :: NormCons

!Just let the initial CI vector be the HF determinant
        CIVect(0)=1.D0
!No excitations in current vector, therefore Dets=0
        Dets=0
!Store the HF determinant
        do i=1,NEl
            CIConfig(0,i)=FDet(i)
        enddo
!Set the excitation level for the HF determinant equal to 0
        CIConfig(0,0)=0

!Calculate initial energy of trial wavevector
        CALL CalcCIVectE()

        WRITE(64,"(I15,2G22.14)") 1,TempEnergyxw/TempWeight,TempWeight

        RETURN
    END SUBROUTINE CreateInitTrialWavevect

!This simply calculates the energy of a CI vector
    SUBROUTINE CalcCIVectE()
        USE System, only: Alat,Beta,Brr,ECore,G1,nBasis,nBasisMax
        USE Integrals, only : fck,nMax,nMsh,UMat
        USE Determinants , only : GetHElement2
        IMPLICIT NONE
        TYPE(HElement) :: Hij
        INTEGER :: i

        TempEnergyxw=0.D0
        TempWeight=0.D0

!Search through the CI vector for double excitations
        do i=1,Dets
            IF(CIConfig(i,0).eq.2) THEN
!Double excitation found - calculate Hij value
                Hij=GetHElement2(FDet,CIConfig(i,1:NEl),NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,2,ECore)
                TempEnergyxw=TempEnergyxw+(Hij%v*CIVect(i))
            ENDIF
        enddo
        TempWeight=CIVect(0)

        RETURN
    END SUBROUTINE CalcCIVectE



END MODULE GraphMorph2



