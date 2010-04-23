! Contain data used by other modules after DetCalc has done its stuff
module DetCalcData
      use constants, only: dp
      INTEGER NDET  ! The total number of determinants we have listed
      INTEGER Det   ! The number of determinants with the same sym
                    ! as the reference det.  This is the number of
                    ! dets in FCIDets
      INTEGER(kind=n_int), Allocatable :: FCIDets(:,:)  !This will contain a list of determinants of the same symmetry as the reference det, with dets in compressed form.  Usually (NIfTot, Det)
      INTEGER, Allocatable :: FCIDetIndex(:)!This indicates where the excitation levels start in the FCIDets array(will go from 0->NEl+1).
      INTEGER ICILEVEL ! The maximum excitation level up to which to enumerate dets.  
      INTEGER :: tagNMRKS=0
      INTEGER, pointer :: NMRKS(:,:)=>null() !(NEL-NFROZEN,nDet)  A list of all determinants which have been enumerated.  
      HElement_t, pointer :: CK(:,:)  !  (nDet,nEval) This will store the eventual eigenvectors
      INTEGER :: tagCK=0
      REAL*8, pointer :: W(:)  ! (nEval) This will contain the eigenvalues
      INTEGER tagW
      
      REAL*8 B2L  ! From Calc
      INTEGER NEVAL  !The number of eigenvectors requested
      INTEGER NBLK   !The number of Lanczos Blocks
      INTEGER NKRY   !The number of Lanczos Krylov vectors

end module
