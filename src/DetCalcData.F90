! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Contain data used by other modules after DetCalc has done its stuff
module DetCalcData
      use constants, only: dp,n_int
        use MemoryManager, only: TagIntType

      INTEGER NDET  ! The total number of determinants we have listed
      INTEGER Det   ! The number of determinants with the same sym
                    ! as the reference det.  This is the number of
                    ! dets in FCIDets
!This will contain a list of determinants of the same symmetry as the reference det, with dets 
!in compressed form.  Usually (NIfTot, Det)
      INTEGER(kind=n_int), Allocatable :: FCIDets(:,:)  
!This indicates where the excitation levels start in the FCIDets array(will go from 0->NEl+1).
      INTEGER, Allocatable :: FCIDetIndex(:)
      INTEGER ICILEVEL ! The maximum excitation level up to which to enumerate dets.  
      INTEGER(TagIntType) :: tagNMRKS=0
      INTEGER, pointer :: NMRKS(:,:)=>null() !(NEL-NFROZEN,nDet)  A list of all determinants which have been enumerated.  
      HElement_t, pointer :: CK(:,:)  !  (nDet,nEval) This will store the eventual eigenvectors
      INTEGER(TagIntType) :: tagCK=0
      real(dp), pointer :: W(:)  ! (nEval) This will contain the eigenvalues
      INTEGER(TagIntType) tagW
      HElement_t, pointer :: HAMIL(:)    !The Hamiltonian in compressed form.  
      !Contains only non-zero elements.  The total number of elements is in LenHamil
      INTEGER(TagIntType) :: tagHamil=0
      INTEGER LenHamil                       !The Total number of non-zero elements in the compressed Hamiltonian

      INTEGER , ALLOCATABLE :: LAB(:),NROW(:),ReIndex(:)
      INTEGER(TagIntType) :: LabTag=0,NRowTag=0
      
      INTEGER NCYCLE !The Max number of Lanczos cycles
      real(dp) B2L  ! From Calc
      INTEGER NEVAL  !The number of eigenvectors requested
      INTEGER NBLK   !The number of Lanczos Blocks
      INTEGER NKRY   !The number of Lanczos Krylov vectors

end module
