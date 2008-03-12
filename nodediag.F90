!This function takes all {a,b} excitations from each i,j pair, and prediagonalises them, fully connecting all the determinants to each other. This is what is meant by a 'node' - a set of double excitations, all of which are excited from the same {i,j}. 
!The transformed objects from this diagonalisation of each node are then attached back to root in a star and solved.
!This should approximate CID, although O[N^2 M] connections are missing throughout the structure, correponding to the few connections there are between nodes.
!Current scaling is N^12 if the final star is explicitly diagonalised, or N^8 if it is done solving the polynomial (to do!).

    MODULE NODEDIAG
      USE HElem
      IMPLICIT NONE

!Stores the excitations by their {a,b} value, according to which {i,j} family they are under.
      INTEGER, ALLOCATABLE :: EXCITSTORE(:,:,:)

!Counts the {a,b} excitations in each {i,j}
      INTEGER, ALLOCATABLE :: ABCOUNTER(:)

!Gives the excited {i,j} orbitals in each node
      INTEGER, ALLOCATABLE :: ijorbs(:,:)

!Stores the final star matrix in the same way as in stardiag.F90
      TYPE(HElement), ALLOCATABLE :: EXCITINFO(:,:)

!The rho_ii rho matrix element
      TYPE(HElement) :: rhii
      REAL*8 :: totlinks
      INTEGER :: crosslinks
      
      contains

      FUNCTION fMCPR3StarNodes(nI,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,nTay,RhoEps,L,LT,nWHTay,iLogging,tSym,ECore,dBeta,dLWdb)

      IMPLICIT NONE
      INCLUDE 'basis.inc'
      TYPE(BasisFN) G1(*)
      INTEGER nI(nEl),nEl,i_P,nBasisMax(5,5),Brr(nBasis),nBasis,nMsh
      INTEGER nMax,nTay(2),L,LT,nWHTay,iLogging,iMaxExcit,nExcitMemLen
      INTEGER noij,noab,ierr,totexcits,nJ(nEl),Orbchange(4),noexcits
      INTEGER Height,TRIIND,INDX,i,ExcitInfoElems,j,exFlag
      INTEGER nStore(6),iExcit,invsbrr(nBasis),orbone,orbtwo,t
      INTEGER, ALLOCATABLE :: nExcit(:)
      COMPLEX*16 fck(*)
      TYPE(HElement) UMat(*)
      REAL*8 Beta,ALat(3),RhoEps,ECore,dBeta
      TYPE(HDElement) dLWdB
      TYPE(HDElement) fMCPR3StarNodes
      TYPE(HElement) HIJS(0:2)
      LOGICAL tSym,COMPIPATH

      IF(HElementSize.GT.1) STOP "NODEDIAG cannot function with complex orbitals currently"

!First need to allocate memory to hold the maximum number of double excitations
!No point to count them explicitly, since even if not connected to root, could still be connected to other double excitations.
      noij=nEl*(nEl-1)/2
      noab=(nBasis-nEl)*(nBasis-nEl-1)/2
      totexcits=noij*noab
      WRITE(6,*) "Total maximum vertices: ",totexcits
      WRITE(6,*) "Total maximum dets in each node: ",noab
      WRITE(6,*) "Total maximum nodes: ",noij
      ALLOCATE(EXCITSTORE(2,noab,noij),stat=ierr)
      CALL MemAlloc(iErr,EXCITSTORE,noij*noab,"EXCITSTORE")
      CALL IAZZERO(EXCITSTORE,2*noij*noab)
      
!Allocate space for counting of ab excitations connected to each ij.
      ALLOCATE(ABCOUNTER(noij),stat=ierr)
      CALL MemAlloc(iErr,ABCOUNTER,noij/2,"ABCOUNTER")
      CALL IAZZERO(ABCOUNTER,noij)

!Create inversebrr for indexing purposes. Brr gives the orbitals in order of energy. If we create an inverse of this array, then we can order the occupied orbitals in ascending order, and so create an indexing scheme for the {i,j} pairs, which need to be stored.
      CALL CREATEINVSBRR(Brr,invsBrr,nBasis)

!Array to detail the ij orbs being excited from in each node
      ALLOCATE(ijorbs(2,noij),stat=ierr)
      CALL MemAlloc(iErr,ijorbs,noij,"IJORBS")
      CALL IAZZERO(IJORBS,noij*2)
      
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
      CALL GenSymExcitIt2(nI,nEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
      Allocate(nExcit(nExcitMemLen))

!Second call to calculate theoretical max number of excitations (iMaxExcit)
      nExcit(1)=0
      CALL GenSymExcitIt2(nI,nEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
      IF(iMaxExcit.gt.(noij*noab)) STOP 'Incorrect calculation of number of excits'
      
!Go through all excitations and store them, even if have no direct connection to root.
      noexcits=0
      lp: do while(.true.)
          CALL GenSymExcitIt2(nI,nEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,exFlag)
          
!Shows that all double excitations have been accounted for
          IF(nJ(1).eq.0) exit lp
          
          IF(iExcit.ne.2) STOP 'TUSEBRILLOUIN must be on for NODEDIAG'
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
      CALL MemAlloc(iErr,EXCITINFO,3*(totexcits+1),"EXCITINFO")
      CALL AZZERO(EXCITINFO,3*(totexcits+1))
      
!Calculate rho_ii and H_ii, and put into ExcitInfo. Again, we divide all rho elements through by rho_ii (Therefore rho_ii element=1)
      CALL CalcRho2(nI,nI,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rhii,nTay,0,ECore)
      ExcitInfo(0,2)=GetHElement2(nI,nI,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,0,ECore)
      EXCITINFO(0,0)=1.D0
      EXCITINFO(0,1)=1.D0

!Go through all ij pairs, create the node matrix and diagonalise, then add elements to ExcitInfo.
      ExcitInfoElems=0
      totlinks=0.D0
      crosslinks=0
      do i=1,noij
          IF(ABCOUNTER(i).eq.0) CYCLE

          CALL CONSTRUCTNODE(ABCOUNTER(i),nEl,i,nI,Beta,i_P,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,nTay,ECore,RhoEps,ExcitInfoElems)

      enddo

      t=0
      do i=1,ExcitInfoElems
        IF(EXCITINFO(i,1).agt.RhoEps) t=t+1
      enddo
      WRITE(6,*) "Number of objects attached to resultant star = ",t
      WRITE(6,*) "Average rhoelement of links between excitations = ", totlinks/crosslinks
      
!Explicitly diagonalise resultant matrix - large scaling.
      CALL StarDiag(0,nEl,ExcitInfoElems+1,EXCITINFO,Totexcits+1,i_P,fmcpr3starnodes,dBeta,dLWdB)

      CALL MemDealloc(EXCITINFO)
      DEALLOCATE(EXCITINFO)
      CALL MemDealloc(EXCITSTORE)
      DEALLOCATE(EXCITSTORE)
      CALL MemDealloc(ABCOUNTER)
      DEALLOCATE(ABCOUNTER)
      CALL MemDealloc(ijorbs)
      DEALLOCATE(ijorbs)
      DEALLOCATE(nExcit)
      
      RETURN
      END FUNCTION fMCPR3StarNodes
      

      
      
!From a given {i,j}, and a list of all {a,b}'s which result in possible double excitations from the HF, find all the connections between them, and diagonalise the resulting matrix from this 'node'. 
!Finally, attach the resultant structures back to the HF in EXCITINFO star matrix.      
      SUBROUTINE CONSTRUCTNODE(novirt,nEl,node,nI,Beta,i_P,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,nTay,ECore,RhoEps,ExcitInfoElems)
        USE Integrals , only : TDISCONODES
        IMPLICIT NONE
        INCLUDE 'basis.inc'
        Type(BasisFN) G1(*)
        COMPLEX*16 fck(*)
        TYPE(HElement) UMat(*),rh,Hel
        INTEGER novirt,ierr,i,j,ijpair(2),node,nI(nEl),nJ(nEl),nK(nEl),i_P
        INTEGER nBasisMax(5,5),nBasis,Brr(nBasis),nMsh,nMax,nTay(2),WORKMEM,INFO
        INTEGER ExcitInfoElems,nEl,Orbchange(4),iExcit
        REAL*8 Beta,ALat(3),RhoEps,ECore
        REAL*8, ALLOCATABLE :: NODERHOMAT(:),WLIST(:)
        INTEGER, ALLOCATABLE :: FULLPATHS(:,:)
        REAL*8, ALLOCATABLE :: WORK(:)

!iExcit should be the order of the excitation - parsed to HElement2
        iExcit=2

!Array to temporarily store full paths of excitations in node - needed for the calcrho2 routine to calculate rho elements.
        ALLOCATE(FULLPATHS(nEl,novirt),stat=ierr)
        CALL MemAlloc(iErr,FULLPATHS,novirt*nEl,"FULLPATHS")
        CALL IAZZERO(FULLPATHS,novirt*nel)

        ijpair(:)=ijorbs(:,node)
            
!Array to store rho elements of node - this is the matrix to be diagonalised.
        ALLOCATE(NODERHOMAT(novirt*novirt),stat=ierr)
        CALL MemAlloc(iErr,NODERHOMAT,novirt*novirt,"NODERHOMAT")
        CALL AZZERO(NODERHOMAT,novirt*novirt)

!First deal with diagonal elements
        do i=1,novirt
!Store full paths as we calculate them
            Orbchange(1:2)=ijpair(:)
            Orbchange(3:4)=EXCITSTORE(1:2,i,node)
            CALL GETFULLPATH(nI,nEl,2,Orbchange,nJ)
            FULLPATHS(:,i)=nJ(:)
!            WRITE(68,*) nJ(:)
!Calculate diagonal rho elements and store them in the matrix            
            CALL CalcRho2(nJ,nJ,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
            NODERHOMAT((i-1)*novirt+i)=DREAL(rh/rhii)
        enddo

!Next, find offdiagonal elements
        do i=1,novirt
            do j=(i+1),novirt
                nJ(:)=FULLPATHS(:,i)
                nK(:)=FULLPATHS(:,j)
                IF(TDISCONODES) THEN
                    rh=HElement(0.D0)
                ELSE
                    CALL CalcRho2(nJ,nK,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,2,ECore)
                ENDIF
!Fill matrix of determinants in node to diagonalise
                IF(rh.agt.RhoEps) THEN
                    crosslinks=crosslinks+1
                    totlinks=totlinks+DREAL(rh)
                    NODERHOMAT(((i-1)*novirt)+j)=DREAL(rh/rhii)
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
        CALL MemAlloc(ierr,WLIST,novirt,"WLIST")
        WORKMEM=3*novirt
        ALLOCATE(WORK(WORKMEM),stat=ierr)
        CALL MemAlloc(ierr,WORK,3*novirt,"WORK")

!Diagonalise
       CALL DSYEV('V','L',novirt,NODERHOMAT,novirt,WLIST,WORK,WORKMEM,INFO)
        IF(INFO.NE.0) THEN
            WRITE(6,*) 'DYSEV error: ',INFO
            STOP
        ENDIF
        CALL MemDealloc(WORK)
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
            IF(WLIST(i).gt.1.D-09) THEN
                ExcitInfoElems=ExcitInfoElems+1
                EXCITINFO(ExcitInfoElems,0)=WLIST(i)
                
!For offdiagonal elements sum over rho element of each determinant (k) to root, times projection of that determinant onto eigenvector j.
! rho_ij = sum_k{ <D_0|rho|k><k|j> } where k are the original determinants in the node, j are the eigenvectors from the diagonalisation, and D_0 is the root.

                DO j=1,novirt
                    nJ(:)=FULLPATHS(:,j)
                    CALL CalcRho2(nI,nJ,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,2,ECore)
!Eigenvectors are COLUMNS
                    EXCITINFO(ExcitInfoElems,1)=EXCITINFO(ExcitInfoElems,1)+(rh/rhii)*HElement(NODERHOMAT((novirt*(j-1))+i))

!H Elements dealt with in the same way
                    Hel=GetHElement2(nI,nJ,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,iExcit,ECore)
                    EXCITINFO(ExcitInfoElems,2)=EXCITINFO(ExcitInfoElems,2)+Hel*HElement(NODERHOMAT((novirt*(j-1))+i))
                ENDDO
            ENDIF
        ENDDO

        CALL MemDealloc(WLIST)
        DEALLOCATE(WLIST)
        CALL MemDealloc(NODERHOMAT)
        DEALLOCATE(NODERHOMAT)
        CALL MemDealloc(FULLPATHS)
        DEALLOCATE(FULLPATHS)

        RETURN
      END SUBROUTINE CONSTRUCTNODE

      END MODULE NODEDIAG

      
      SUBROUTINE CREATEINVSBRR(Brr,INVSBRR,nBasis)
      IMPLICIT NONE
      INTEGER :: nBasis,Brr(nBasis),INVSBrr(nBasis),i
      CALL IAZZERO(INVSBRR,NBasis)
      DO i=1,nbasis
          INVSBrr(Brr(i))=i
      ENDDO
      RETURN
      END SUBROUTINE CREATEINVSBRR
      
      
!Returns the index in the third EXCITSTORE dimension which will hold the ij variable. Stored in triangular form, without diagonals.
      FUNCTION TRIIND(I,J,nEl)
      IMPLICIT NONE
      INTEGER TRIIND,I,J,nEl
      IF((I.gt.nEl).or.(J.gt.nEl)) THEN
          STOP 'PROBLEM WITH INDEXING'
      ENDIF
      IF(I.GT.J) THEN
          TRIIND=((I-1)*(I-2)/2)+J
          RETURN
      ELSEIF(J.GT.I) THEN
          TRIIND=((J-1)*(J-2)/2)+I
          RETURN
      ELSE
          STOP 'I should not equal J in TRIIND'
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
              EXIT
          ENDIF
      ENDDO
      RETURN
      END FUNCTION COMPIPATH
      
