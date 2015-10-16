! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
!This function takes all {a,b} excitations from each i,j pair, and prediagonalises them, fully connecting all the 
!determinants to each other. This is what is meant by a 'node' - a set of double excitations, 
!all of which are excited from the same {i,j}. 
!The transformed objects from this diagonalisation of each node are then attached back to root in a star and solved.
!This should approximate CID, although O[N^2 M] connections are missing throughout the structure, correponding to 
!the few connections there are between nodes.
!Current scaling is N^12 if the final star is explicitly diagonalised, or N^8 if it is 
!done solving the polynomial (to do!).

!!BEWARE!! Bugs expected in code - should NOT reattach the original determinants - have already 
!included them in the prediagonalisation when
!considering the projection of the original determinant. Needs fixing

!Also to do : consider the linear approximation of eigenvectors, which would mean that only a few of the nodes would 
!need to be explicitly diagonalised, before a linear approximation
!could be taken to approximate the eigensystems from the other nodes. This would rely on the matrices begin the same, 
!and simply the diagonal elements being multiplied by a constant
!This would reduce the scaling to M^6 - same as CID

    MODULE NODEDIAG
      use constants, only: dp,int32
      use SystemData, only: BasisFN
      use IntegralsData, only: tDiscoNodes
      use Determinants, only: get_helement, get_helement_excit
      use global_utilities
      use MemoryManager, only: TagIntType
      implicit none

!Stores the excitations by their {a,b} value, according to which {i,j} family they are under.
      INTEGER, ALLOCATABLE :: EXCITSTORE(:,:,:)

!Counts the {a,b} excitations in each {i,j}
      INTEGER, ALLOCATABLE :: ABCOUNTER(:)

!Gives the excited {i,j} orbitals in each node
      INTEGER, ALLOCATABLE :: ijorbs(:,:)

!Stores the final star matrix in the same way as in stardiag.F90
      HElement_t(dp), ALLOCATABLE :: EXCITINFO(:,:)

!The rho_ii rho matrix element
      HElement_t(dp) :: rhii
      real(dp) :: totlinks
      INTEGER :: crosslinks

! Memory tags
      integer(TagIntType), save :: tagEXCITINFO=0,tagEXCITSTORE=0,tagABCOUNTER=0,tagijorbs=0
      
      contains

      FUNCTION fMCPR3StarNodes(nI,Beta,i_P,nEl,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,nTay,RhoEps,ECore,dBeta,dLWdb)
      use HElem
      use neci_intfce
      TYPE(BasisFN) G1(*)
      INTEGER nEl,nI(nEl),nBasis,i_P,Brr(nBasis),nMsh
      INTEGER nMax,nTay(2),iMaxExcit,nExcitMemLen(1)
      INTEGER noij,noab,ierr,totexcits,nJ(nEl),Orbchange(4),noexcits
      INTEGER Height,TRIIND,INDX,i,ExcitInfoElems,j,exFlag
      INTEGER nStore(6),iExcit,invsbrr(nBasis),orbone,orbtwo,t
      INTEGER, ALLOCATABLE :: nExcit(:)
      complex(dp) fck(*)
      HElement_t(dp) UMat(*)
      real(dp) Beta,ALat(3),RhoEps,ECore,dBeta
      real(dp) dLWdB
      real(dp) fMCPR3StarNodes
      LOGICAL COMPIPATH
      character(*), parameter :: t_r='fMCPR3StarNodes'

      IF(HElement_t_size.GT.1) call stop_all(t_r, "NODEDIAG cannot function with complex orbitals currently")

!First need to allocate memory to hold the maximum number of double excitations
!No point to count them explicitly, since even if not connected to root, could still be connected to other double excitations.
      noij=nEl*(nEl-1)/2
      noab=(nBasis-nEl)*(nBasis-nEl-1)/2
      totexcits=noij*noab
      WRITE(6,*) "Total maximum vertices: ",totexcits
      WRITE(6,*) "Total maximum dets in each node: ",noab
      WRITE(6,*) "Total maximum nodes: ",noij
      ALLOCATE(EXCITSTORE(2,noab,noij),stat=ierr)
      CALL LogMemAlloc("EXCITSTORE",noij*noab,8,t_r,tagEXCITSTORE,iErr)
      EXCITSTORE(1:2,1:noij,1:noab)=0
      
!Allocate space for counting of ab excitations connected to each ij.
      ALLOCATE(ABCOUNTER(noij),stat=ierr)
      CALL LogMemAlloc("ABCOUNTER",noij/2,8,t_r,tagABCOUNTER,iErr)
      ABCOUNTER(1:noij)=0

!Create inversebrr for indexing purposes. Brr gives the orbitals in order of energy. If we create an inverse of this array, 
!then we can order the occupied orbitals in ascending order, and so create an indexing scheme for the {i,j} pairs, 
!which need to be stored.
      CALL CREATEINVSBRR(Brr,invsBrr,nBasis)

!Array to detail the ij orbs being excited from in each node
      ALLOCATE(ijorbs(2,noij),stat=ierr)
      CALL LogMemAlloc("IJORBS",noij,8,t_r,tagijorbs,iErr)
      IJORBS(1:2,1:noij)=0
      
!Fill array with ij orbs in each node
      do i=1,nEl
          do j=(i+1),nEl
              INDX=TRIIND(i,j,nEl)
              ijorbs(1,INDX)=Brr(i)
              ijorbs(2,INDX)=Brr(j)
!              WRITE(6,*) i,j,"Have index",INDX
          enddo
      enddo

!exFlag set to 2 to signify only wanting double excitations
      exFlag=2

!First call finds memory needed for excitation generator (nExcitMemLen)
!      nExcitMemLen=0
      nStore(1)=0
      CALL GenSymExcitIt2(nI,nEl,G1,nBasis,.TRUE.,nExcitMemLen,nJ,iMaxExcit,nStore,exFlag)
      Allocate(nExcit(nExcitMemLen(1)))

!Second call to calculate theoretical max number of excitations (iMaxExcit)
      nExcit(1)=0
      CALL GenSymExcitIt2(nI,nEl,G1,nBasis,.TRUE.,nExcit,nJ,iMaxExcit,nStore,exFlag)
      IF(iMaxExcit.gt.(noij*noab)) call stop_all(t_r, 'Incorrect calculation of number of excits')
      
!Go through all excitations and store them, even if have no direct connection to root.
      noexcits=0
      lp: do while(.true.)
          CALL GenSymExcitIt2(nI,nEl,G1,nBasis,.false.,nExcit,nJ,iExcit,nStore,exFlag)
          
!Shows that all double excitations have been accounted for
          IF(nJ(1).eq.0) exit lp
          
          IF(iExcit.ne.2) call stop_all(t_r, 'TUSEBRILLOUIN must be on for NODEDIAG')
          IF(COMPIPATH(nI,nJ,nEl)) THEN
              WRITE(6,*) "WARNING, nI and nJ the same!"
              WRITE(6,*) "nI = ", nI
              WRITE(6,*) "nJ = ", nJ
          ENDIF
          CALL GETEXCITSCHANGE(nI,nJ,nEl,Orbchange(:))
!The inverse of Brr is needed so that occupied orbitals are labelled in numerically ascending order.
          orbone=INVSBRR(Orbchange(1))
          orbtwo=INVSBRR(Orbchange(2))
!TRIIND stores each ij pair in a unique place
          INDX=TRIIND(orbone,orbtwo,nEl)
          Height=ABCOUNTER(INDX)+1
          EXCITSTORE(1:2,Height,INDX)=Orbchange(3:4)
          ABCOUNTER(INDX)=ABCOUNTER(INDX)+1
          noexcits=noexcits+1
      enddo lp
!      DO i=1,noij
!            WRITE(6,*) ABCOUNTER(i)
!      ENDDO

!Allocate Memory for ExcitInfo
      ALLOCATE(EXCITINFO(0:totexcits,0:2),stat=ierr)
      CALL LogMemAlloc("EXCITINFO",3*(totexcits+1),8,t_r,tagEXCITINFO,iErr)
      EXCITINFO=(0.0_dp)
      
!Calculate rho_ii and H_ii, and put into ExcitInfo. Again, we divide all rho elements through by rho_ii (Therefore rho_ii element=1)
      CALL CalcRho2(nI,nI,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rhii,nTay,0,ECore)
      ExcitInfo(0,2) = get_helement (nI, nI, 0)
      EXCITINFO(0,0)=1.0_dp
      EXCITINFO(0,1)=1.0_dp

!Go through all ij pairs, create the node matrix and diagonalise, then add elements to ExcitInfo.
      ExcitInfoElems=0
      totlinks=0.0_dp
      crosslinks=0
      do i=1,noij
          IF(ABCOUNTER(i).eq.0) CYCLE

          CALL CONSTRUCTNODE(ABCOUNTER(i),nEl,i,nI,Beta,i_P,G1,nBasis,nMsh,fck,nMax,ALat,UMat,nTay,ECore,RhoEps,ExcitInfoElems)

      enddo

      t=0
      do i=1,ExcitInfoElems
        IF(abs(EXCITINFO(i,1)).gt.RhoEps) t=t+1
      enddo
      WRITE(6,*) "Number of objects attached to resultant star = ",t
      WRITE(6,*) "Average rhoelement of links between excitations = ", totlinks/crosslinks
      
!Explicitly diagonalise resultant matrix - large scaling.
      CALL StarDiag(ExcitInfoElems+1,EXCITINFO,Totexcits+1,i_P,fmcpr3starnodes,dBeta,dLWdB)

      CALL LogMemDealloc(t_r,tagExcitInfo)
      DEALLOCATE(EXCITINFO)
      CALL LogMemDealloc(t_r,tagEXCITSTORE)
      DEALLOCATE(EXCITSTORE)
      CALL LogMemDealloc(t_r,tagABCOUNTER)
      DEALLOCATE(ABCOUNTER)
      CALL LogMemDealloc(t_r,tagijorbs)
      DEALLOCATE(ijorbs)
      DEALLOCATE(nExcit)
      
      RETURN
      END FUNCTION fMCPR3StarNodes
      

      
      
!From a given {i,j}, and a list of all {a,b}'s which result in possible double excitations from the HF, find 
!all the connections between them, and diagonalise the resulting matrix from this 'node'. 
!Finally, attach the resultant structures back to the HF in EXCITINFO star matrix.      
      SUBROUTINE CONSTRUCTNODE(novirt,nEl,node,nI,Beta,i_P,G1,nBasis,nMsh,fck,nMax,ALat,UMat,nTay,ECore,RhoEps,ExcitInfoElems)
        use MemoryManager, only: TagIntType
        IMPLICIT NONE
        Type(BasisFN) G1(*)
        complex(dp) fck(*)
        HElement_t(dp) UMat(*),rh,Hel
        INTEGER novirt,nEl,ierr,i,j,ijpair(2),node,nI(nEl),nJ(nEl),nK(nEl),i_P
        INTEGER nBasis,nMsh,nMax,nTay(2),WORKMEM
        INTEGER(int32) INFO
        INTEGER ExcitInfoElems,Orbchange(4),iExcit
        real(dp) Beta,ALat(3),RhoEps,ECore
        real(dp), ALLOCATABLE :: NODERHOMAT(:),WLIST(:)
        INTEGER, ALLOCATABLE :: FULLPATHS(:,:)
        real(dp), ALLOCATABLE :: WORK(:)
        integer(TagIntType), save :: tagNODERHOMAT=0,tagWLIST=0,tagFULLPATHS=0,tagWORK=0
        character(*),parameter :: t_r='CONSTRUCTNODE'

!iExcit should be the order of the excitation - parsed to HElement2
        iExcit=2

!Array to temporarily store full paths of excitations in node - needed for the calcrho2 routine to calculate rho elements.
        ALLOCATE(FULLPATHS(nEl,novirt),stat=ierr)
        CALL LogMemAlloc("FULLPATHS",novirt*nEl,8,t_r,tagFULLPATHS,iErr)
        FULLPATHS=0

        ijpair(:)=ijorbs(:,node)
            
!Array to store rho elements of node - this is the matrix to be diagonalised.
        ALLOCATE(NODERHOMAT(novirt*novirt),stat=ierr)
        CALL LogMemAlloc("NODERHOMAT",novirt*novirt,8,t_r,tagNODERHOMAT,iErr)
        NODERHOMAT=0.0_dp

!First deal with diagonal elements
        do i=1,novirt
!Store full paths as we calculate them
            Orbchange(1:2)=ijpair(:)
            Orbchange(3:4)=EXCITSTORE(1:2,i,node)
            CALL GETFULLPATH(nI,nEl,2,Orbchange,nJ)
            FULLPATHS(:,i)=nJ(:)
!            WRITE(68,*) nJ(:)
!Calculate diagonal rho elements and store them in the matrix            
            CALL CalcRho2(nJ,nJ,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
            NODERHOMAT((i-1)*novirt+i)=(rh/rhii)
        enddo

!Next, find offdiagonal elements
        do i=1,novirt
            do j=(i+1),novirt
                nJ(:)=FULLPATHS(:,i)
                nK(:)=FULLPATHS(:,j)
                IF(TDISCONODES) THEN
                    rh=(0.0_dp)
                ELSE
                    CALL CalcRho2(nJ,nK,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,2,ECore)
                ENDIF
!Fill matrix of determinants in node to diagonalise
                IF(abs(rh).gt.RhoEps) THEN
                    crosslinks=crosslinks+1
                    totlinks=totlinks+(rh)
                    NODERHOMAT(((i-1)*novirt)+j)=(rh/rhii)
                ENDIF
            enddo
        enddo

!OPTIONAL - remove determinants from matrix which are not attached to any other (TO DO - not finished)
!        do i=1,novirt
!            diag=(i-1)*novirt+i
!            KEEP=.false.
!            do j=1,novirt
!                elem=novirt*(j-1)+i
!                IF(elem.lt.diag) THEN
!                    IF(ABS(NODERHOMAT).gt.RhoEps) THEN
!                        KEEP=.true.
!                        EXIT
!                    ENDIF
!                ELSEIF(elem.gt.diag) THEN
!                    elem=novirt*(i-1)+j
!                    IF(ABS(NODERHOMAT).gt.RhoEps) THEN
!                        KEEP=.true.
!                        EXIT
!                    ENDIF
!                ELSEIF(elem.eq.diag) THEN
!                    CYCLE
!                ENDIF
!            enddo
!            IF(.NOT.KEEP) THEN
!                !REMOVE DETERMINANT FROM MATRIX
!            ENDIF
!        enddo

!Debugging - write out noderhomat
!        DO I=1,novirt
!            DO J=1,novirt
!                WRITE(68,"E14.6,$") NODERHOMAT((I-1)*novirt+j)
!            ENDDO
!            WRITE(68,*) ""
!            WRITE(68,*) ""
!        ENDDO
!        WRITE(68,*) "****************"
                

!Allocate Memory for diagonalisation
        ALLOCATE(WLIST(novirt),stat=ierr)
        CALL LogMemAlloc("WLIST",novirt,8,t_r,tagWLIST,ierr)
        WORKMEM=3*novirt
        ALLOCATE(WORK(WORKMEM),stat=ierr)
        CALL LogMemAlloc("WORK",3*novirt,8,t_r,tagWORK,ierr)

!Diagonalise
       CALL DSYEV('V','L',novirt,NODERHOMAT,novirt,WLIST,WORK,WORKMEM,INFO)
        IF(INFO.NE.0) THEN
            WRITE(6,*) 'DYSEV error: ',INFO
            call stop_all(t_r, "")
        ENDIF
        CALL LogMemDealloc(t_r,tagWORK)
        DEALLOCATE(WORK)

!Debugging - write out eigenvectors
!        DO I=1,novirt
!            DO J=1,novirt
!                WRITE(68,"E14.6,$") NODERHOMAT((I-1)*novirt+j)
!            ENDDO
!            WRITE(68,*) ""
!            WRITE(68,*) ""
!        ENDDO
!        WRITE(68,*) "*****************************************"

!Attached transformed objects back to root as a star graph...

!Diagonal elements of star graph are simply eigenvalues - add to ExcitInfo
        DO i=1,novirt
            IF(WLIST(i).gt.1.0e-9_dp) THEN
                ExcitInfoElems=ExcitInfoElems+1
                EXCITINFO(ExcitInfoElems,0)=WLIST(i)
                
!For offdiagonal elements sum over rho element of each determinant (k) to root, 
!times projection of that determinant onto eigenvector j.
! rho_ij = sum_k{ <D_0|rho|k><k|j> } where k are the original determinants in the 
!node, j are the eigenvectors from the diagonalisation, and D_0 is the root.

                DO j=1,novirt
                    nJ(:)=FULLPATHS(:,j)
                    CALL CalcRho2(nI,nJ,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,2,ECore)
!Eigenvectors are COLUMNS
                    EXCITINFO(ExcitInfoElems,1)=EXCITINFO(ExcitInfoElems,1)+(rh/rhii)*(NODERHOMAT((novirt*(j-1))+i))

!H Elements dealt with in the same way
                    Hel = get_helement(nI, nJ)
                    EXCITINFO(ExcitInfoElems,2)=EXCITINFO(ExcitInfoElems,2)+Hel*(NODERHOMAT((novirt*(j-1))+i))
                ENDDO
            ENDIF
        ENDDO

        CALL LogMemDealloc(t_r,tagWLIST)
        DEALLOCATE(WLIST)
        CALL LogMemDealloc(t_r,tagNODERHOMAT)
        DEALLOCATE(NODERHOMAT)
        CALL LogMemDealloc(t_r,tagFULLPATHS)
        DEALLOCATE(FULLPATHS)

        RETURN
      END SUBROUTINE CONSTRUCTNODE

      END MODULE NODEDIAG

      
      SUBROUTINE CREATEINVSBRR(Brr,INVSBRR,nBasis)
      IMPLICIT NONE
      INTEGER :: nBasis,Brr(nBasis),INVSBrr(nBasis),i
      INVSBRR(1:NBasis)=0
      DO i=1,nbasis
          INVSBrr(Brr(i))=i
      ENDDO
      RETURN
      END SUBROUTINE CREATEINVSBRR
      
      
!Returns the index in the third EXCITSTORE dimension which will hold the ij variable. Stored in triangular form, without diagonals.
      FUNCTION TRIIND(I,J,nEl)
      IMPLICIT NONE
      INTEGER TRIIND,I,J,nEl
      character(*), parameter :: this_routine = 'TRIIND'
      IF((I.gt.nEl).or.(J.gt.nEl)) THEN
          call stop_all(this_routine, 'PROBLEM WITH INDEXING')
      ENDIF
      IF(I.GT.J) THEN
          TRIIND=((I-1)*(I-2)/2)+J
          RETURN
      ELSEIF(J.GT.I) THEN
          TRIIND=((J-1)*(J-2)/2)+I
          RETURN
      ELSE
          call stop_all(this_routine, 'I should not equal J in TRIIND')
      ENDIF
      END FUNCTION TRIIND

!Compare two determinants and return true if they are the same, or false otherwise.
      LOGICAL FUNCTION COMPIPATH(nI,nJ,nEl)
      IMPLICIT NONE
      INTEGER :: nEl,nI(nEl),nJ(nEl),I

      COMPIPATH=.TRUE.
      DO I=1,nEl
          IF(nI(I).ne.nJ(I)) THEN
              COMPIPATH=.FALSE.
              RETURN
          ENDIF
      ENDDO
      RETURN
      END FUNCTION COMPIPATH
      
