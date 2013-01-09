    MODULE STARDIAGMOD
        use Determinants, only: get_helement
        use sort_mod
        use constants, only: dp,int32
        use MemoryManager, only: TagIntType
        use util_mod, only: NECI_ICOPY,get_unique_filename
        IMPLICIT NONE
      
!.. ExcitInfo will contain all the info needed to work out the value of the star
!.. ExcitInfo(0,...) corresponds to J=I
!.. ExcitInfo(J,0) = RHOJJ
!.. ExcitInfo(J,1) = RHOIJ
!.. ExcitInfo(J,2) = HIJ
         HElement_t, POINTER :: ExcitInfo(:,:)
         HElement_t, ALLOCATABLE, TARGET :: ExcitInfo2(:,:)
         integer(TagIntType) :: tagExcitInfo=0
         integer(TagIntType) :: tagExcitInfo2=0
         HElement_t, ALLOCATABLE :: temprhos(:,:)

!Offrho is used with TSumProd to store the original star off-diagonal elements when using TSumProd, which changes the values in EXCITINFO
         HElement_t, ALLOCATABLE :: Offrho(:)

!If we are keeping a list of excitations in the star, then allocate EXCITSTORE to hold the excitations, in the form of o o v v orbitals
         INTEGER, ALLOCATABLE :: EXCITSTORE(:,:)

!Prodpositions stores the position of the two excitations in EXCITSTORE which give rise to a real product excitation
         INTEGER, ALLOCATABLE :: ProdPositions(:,:)

! When calculating full products, the rho elements for these quadruple product excitations are given in OnDiagProdRho and OffDiagProdRho.
         real(dp), ALLOCATABLE :: OnDiagProdRho(:)
         real(dp), ALLOCATABLE :: OffDiagProdRho(:,:)

!Used with TSumProd, rhiiadd is the value to resum back into the root to account for quadruple excitations
         HElement_t :: rhiiadd

         integer(TagIntType), save :: tagOnDiagProdRho=0,tagOffDiagProdRho=0,tagEXCITSTORE=0,tagOffrho=0,tagtemprhos=0
         integer(TagIntType), save :: tagProdPositions=0

         contains
!  A function to generate all possible excitations of a determinant and link them together in a star
!    The rho_jj, rho_ij and H_ij values are stored for each connected determinant.
!   Based on a combined FMCPR3STAR and FMCPR3STAR2, this instead generates excitations on the fly.
   FUNCTION fMCPR3StarNewExcit(nI,Beta,i_P,nEl,nBasisMax,G1,nBasis,nMsh,fck,nMax,ALat,UMat,nTay, &
               RhoEps, L, LT,nWHTay, iLogging, ECore,dBeta,dLWdB,MP2E)
         use CalcData , only : TMPTHEORY,StarProd,TStarStars,TLanczos,TMCStar
         use SystemData , only : TSTOREASEXCITATIONS,BasisFN
         use IntegralsData , only : TCalcRhoProd,TSumProd,TCalcRealProd,TCalcExcitStar,TDiagStarStars,TLinRootChange
         use global_utilities
         use HElem
         IMPLICIT NONE
         Type(BasisFN) G1(*)
         INTEGER nEl,nI(nEl),i_P,nBasisMax(5,*),nBasis,nMsh
         INTEGER nMax,nTay(2),L,LT,nWHTay,iLogging
         complex(dp) fck(*)
         HElement_t UMat(*)
         real(dp) Beta, ALat(3),RhoEps,ECore,dBeta,MaxDiag
         real(dp) dLWdB
         real(dp) fMCPR3StarNewExcit
         HElement_t HIJS(0:2)
!         real(dp) LARGERHOJJ(10000)
         INTEGER iPath(nEl,0:2),UniqProd
         character(*), parameter :: this_routine='fMCPR3StarNewExcit'
!.. New lists are generated here
!.. This will contain all the info needed to work out the value of the
!.. star
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ

         INTEGER exFlag
         INTEGER, allocatable :: nExcit(:)
         INTEGER nExcitMemLen(1),nStore(6)
         INTEGER nJ(nEl),iExcit,iMaxExcit,excitcount,Prodnum
         INTEGER iErr
         INTEGER nRoots,i
         HElement_t rh,rhii,EHFDiff,Hii
         real(dp) MP2E(2:2)        
         LOGICAL tStarSingles,tCountExcits
         INTEGER nIExcitFormat(nEl)
!         real(dp) , ALLOCATABLE :: SortedHij(:)
!         real(dp) , Norm
         
!This needs to be removed, as it'll eventually be an input parameter

!         LARGERHOJJ(:)=0.0_dp
         hii = 0
         fMCPR3StarNewExcit=(0.0_dp)
         IF(tStoreAsExcitations) THEN
            nIExcitFormat(1)=-1
            nIExcitFormat(2)=0
         ELSE
            CALL NECI_ICOPY(NEL,nI,1,nIExcitFormat,1)
         ENDIF
         tStarSingles=BTEST(nWHTay,7)
         tCountExcits=BTEST(nWHTay,8)
         SELECT CASE (IAND(nWHTay,24))
         CASE(0)
!.. Allow both singles and doubles
            exFlag=3
         CASE(8)
!.. only singles
            WRITE(6,*) "FMCPR3Star Only singles."
            exFlag=1
         CASE(16)
            WRITE(6,*) "FMCPR3Star Only doubles."
            exFlag=2
         CASE(24)
            STOP "Invalid combination of flags in NWHTAY"
         END SELECT
!.. Count the excitations. - First call of GenSymExcitIt2 calculates memory needed for internal use in excitation generators

         nStore(1)=0
         CALL GenSymExcitIt2(nI,nEl,G1,nBasis,.TRUE.,nExcitMemLen,nJ,iMaxExcit,nStore,exFlag)
         Allocate(nExcit(nExcitMemLen(1)))
!Second call calculates size of arrays needed to store all symmetry allowed excitations - further calls will generate excitation on-the-fly(shown by the false in arg(6)
         nExcit(1)=0
         CALL GenSymExcitIt2(nI,nEl,G1,nBasis,.TRUE.,nExcit,nJ,iMaxExcit,nStore,exFlag)

!.. iMaxExcit now contains the number of excitations.
         !TCountExcits will run through all excitations possible, determine if they are connected, and then only store these.
         !Will be twice as expensive, as needs to run through all excitations twice - however, will only store memory needed.
         IF(tCountExcits) THEN
            Write(6,"(A,I10,A)") "Counting excitations - Running through all ",iMaxExcit, &
                " excitations to determine number connected"
            excitcount=0
            CALL CalcRho2(nI,nI,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rhii,nTay,0,ECore)
       lp2: do while(.true.)
                CALL GenSymExcitIt2(nI,nEl,G1,nBasis,.false.,nExcit,nJ,iExcit,nStore,exFlag)
                IF(nJ(1).eq.0) exit lp2
                CALL CalcRho2(nIExcitFormat,nJ,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,iExcit,ECore)
                IF(abs(rh).gt.RhoEps) excitcount=excitcount+1
            enddo lp2

            !Set number of excitations to number of connected determinants, and reset generator
            !This routine doesn't seem to work...would be good if didn't need to reinitialise excitation generator
!            CALL ResetExit2(nI,nEl,G1,nBasis,nBasisMax,nExcit,0)
            Deallocate(nExcit)
            nStore(1)=0
            CALL GenSymExcitIt2(nI,nEl,G1,nBasis,.TRUE.,nExcitMemLen,nJ,iMaxExcit,nStore,exFlag)
            Allocate(nExcit(nExcitMemLen(1)))
            nExcit(1)=0
            CALL GenSymExcitIt2(nI,nEl,G1,nBasis,.TRUE.,nExcit,nJ,iMaxExcit,nStore,exFlag)

            !set number of excitations to exact number there are
            iMaxExcit=excitcount
            
         ENDIF

!.. Allocate memory for the lists
         IF(tStarSingles) THEN
            iMaxExcit=iMaxExcit*5
            WRITE(6,*) "Adding StarSingles.  (5x storage space used)."
         ENDIF
         Write(6,*) "Allocating storage for ",iMaxExcit," excitations."
         Allocate(ExcitInfo(0:iMaxExcit,0:2),stat=iErr)
         call LogMemAlloc('ExcitInfo',3*(iMaxExcit+1),8*HElement_t_size,this_routine,tagExcitInfo,iErr)
!If we are keeping a list of excitations in the star, then allocate EXCITSTORE to hold the excitations, in the form of o o v v orbitals
         IF(TCalcRealProd.or.TSumProd.or.TCalcExcitStar) THEN
             ALLOCATE(EXCITSTORE(4,iMaxExcit),stat=iErr)
             CALL LogMemAlloc("EXCITSTORE",(iMaxExcit*2),4,this_routine,tagEXCITSTORE,iErr)
             ExcitStore(:,:)=0
         ENDIF

!.. This will contain all the info needed to work out the value of the
!.. star
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ
         i=0
         ExcitInfo(i,0)=1.0_dp
         ExcitInfo(i,1)=1.0_dp
         ExcitInfo(i,2) = get_helement (nI, nI, 0)
         IF(TMCStar) THEN
!Hii is HF energy - <D0|H|D0>, ExcitInfo(i,1) now is <Di|H|Di>-Hii. 
!All diagonal elements are therefore positive, and increasing down the leading diagonal.
!MaxDiag is the largest diagonal element.
             Hii=ExcitInfo(0,2)
             ExcitInfo(i,0)=0.0_dp
             MaxDiag=0.0_dp
         ENDIF
         ehfdiff = 0
         if(BTEST(nWhTay,5)) then
! We use the Zeroth order N-particle Hartree-Fock hamiltonian (as MP theory), but shifted by E_HF-E0.

!nMax has Arr hidden in it
            IF(nTay(2).eq.5) THEN
               call GetH0ElementDCCorr(nI,nI,nEl,G1,nBasis,nMax,ECore,rhii)
            ELSE
               call GetH0Element(nI,nEl,nMax,nBasis,ECore,rhii)
            ENDIF
            EHFDiff=ExcitInfo(i,2)-rhii
         endif
         CALL CalcRho2(nI,nI,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rhii,nTay,0,ECore)
!         write(6,*) "rhoii is", rhii
! Setup MP info
         CALL NECI_ICOPY(nEl,nI,1,iPath(1,0),1)
         CALL NECI_ICOPY(nEl,nI,1,iPath(1,2),1)
         HIJS(0)=ExcitInfo(0,2)
         HIJS(2)=ExcitInfo(0,2)

    lp:  do while(.true.)
            CALL GenSymExcitIt2(nI,nEl,G1,nBasis,.false.,nExcit,nJ,iExcit,nStore,exFlag)
            IF(nJ(1).eq.0) exit lp

!           Calculate rhoij element
            CALL CalcRho2(nIExcitFormat,nJ,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,iExcit,ECore)
            
            if(abs(rh ).gt. RhoEps) then
           
!               WRITE(33,*) nJ(:)
               i=i+1
!   Divide all elements though by rhoii
               ExcitInfo(i,1)=rh/rhii
               
               IF(TCalcRealProd.or.TSumProd.or.TCalcExcitStar) THEN
!   Stores all excitations as (occ,occ,vir,vir)
                   CALL GETEXCITSCHANGE(nI,nJ,nEl,EXCITSTORE(:,i))
               ENDIF
               if(btest(nwhtay,5)) then
                  call GetH0Element(nJ,nEl,nMax,nBasis,ECore,rh)
                  rh=rh+EHFDiff
                  rh=rh*(-Beta/I_P)
                  rh=exp(rh)
                  ExcitInfo(i,0)=rh/rhii
               else
                  !RHO_JJ elements
                   IF(TMCStar) THEN
!If we are solving star using MC, then we want the Hamiltonian matrix, rather than rho matrix elements for diagonal elements, and subtract the HF energy from them all
                       ExcitInfo(i,0) = get_helement (nJ, nJ, 0)
                       ExcitInfo(i,0)=ExcitInfo(i,0)-Hii
                       IF(abs(ExcitInfo(i,0)).gt.MaxDiag) MaxDiag=ExcitInfo(i,0)
                   ELSE
                       CALL CalcRho2(nJ,nJ,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
                       ExcitInfo(i,0)=rh/rhii
                   ENDIF
                  
!                  do j=1,10000
!                    IF(abs((rh).gt.LARGERHOJJ(J)).or.(LARGERHOJJ(J).eq.0.0_dp)) THEN
!                        LARGERHOJJ(J)=rh
!                        GOTO 765
!                    ENDIF
!                  ENDDO
!765               CONTINUE
               endif
               ExcitInfo(i,2) = get_helement (nIExcitFormat, nJ, iExcit)
!               write(75,*) rh,rh/rhii
!Now do MP2
               Hijs(1)=ExcitInfo(i,2)
               IF(tMPTheory) THEN
                  caLL NECI_ICOPY(nEl,nJ,1,iPath(1,1),1)
!nMax has Arr hidden in it
                  Call AddMP2E(Hijs,nMax,nBasis,iPath,nEl,BTEST(iLogging,0),MP2E)
               ENDIF
               IF(tStarSingles) Call StarAddSingles(nI,nJ,ExcitInfo,i,iMaxExcit,rhii,rhoeps,Beta,i_P,nEl,G1,nBasis, &
                    nMsh,fck,nMax,ALat,UMat,nTay,ECore)
            endif
         enddo lp
!Tell MCPATHS how many excitations there were and how many we are keeping
         L=i
         LT=iMaxExcit
         iExcit=i
         Deallocate(nExcit)

!Quick test to see range of Hij values
!         ALLOCATE(SortedHij(iExcit))
!         do l=1,iExcit
!            SortedHij(l)=ABS(REAL(ExcitInfo(l,2),KIND(0.0_dp)))
!            Norm=Norm+SortedHij(l)**5
!         enddo
!         CALL SORT(iExcit,SortedHij)
!         do l=1,iExcit
!            WRITE(19,*) l,SortedHij(l),(SortedHij(l)**5)/Norm
!         enddo
            

!.. we now have a list length NLCUR of dets in the star.
!.. Call a routine to generate the value of the star
         WRITE(6,*) iExcit," excited determinants in star"

!Routine for debugging/dev research to graph various eigenvalues and eigenvectors of varying star matrices - should generally be commented out.
!         CALL GraphRootChange(iMaxExcit,iExcit)

         
!If starprod is set, it means that one of a number of methods is used to attempt to indroduce quadruple excitations into the star in an approximate way to achieve size consistency for dissociation into two fragments.
         IF(StarProd) THEN
             CALL GetStarProds(iExcit,ProdNum,UniqProd,rhii,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,ECore)
         ENDIF
         
         IF(TStarStars) THEN
             IF(TCalcExcitStar) THEN

!CalcExcitStar explicitly calculates the excitations from each double excitation, and 
!forms a star graph out of these, prediagonalising them, and adding them to the original 
!star. This scales as N^8 M^8 as each excited star is fully diagonalised. The excitations 
!from these stars can be limited to quadruple excitations of the HF, or remove the double excitations.
                 CALL CalcExcitStar(iMaxExcit,iExcit,nI,rhii,Beta,i_p,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,nTay,ECore,RhoEps)
             
             ELSEIF(TDiagStarStars) THEN

!GetStarStars approximates excited stars as having the same connections as the original star, and so simply multiplies the diagonal elements by rho_jj and then diagonalises them.
                 CALL GetStarStars(iMaxExcit,iExcit,RhoEps)
             ELSEIF(TLinRootChange) THEN
 
!GetLinRootChangeStars creates excited stars, which are the same as the original star, but 
!the root changes. This form for the excited stars means that the eigenvectors are only 
!unity over disjoint ranges of the root value. This means that only one eigenvector is 
!attached per excited graph, and thus the final graph to diagonalise is the same size as 
!the original graph. It appears at the moment, that all the excited stars are within the 
!first eigenvector range. The eigenvalues within this range change linearly from the highest 
!eigenvalue of the original star. Once this gradient is found, the change in the diagonal 
!elements of the original matrix can be calculated and shifted. Offdiagonal elements, and 
!Hamiltonian elements are unchanged, since the eigenvectors are equal to unity.
                 CALL GetLinRootChangeStars(iMaxExcit,iExcit,nWHTay)
             ELSE

!GetLinStarStars assumes a linear variation of the star eigenvalues and vectors, and so uses this to approximate the diagonalised form for excited star matrices, given by multiplying the original star matrix by the rho_jj values.
                 CALL GetLinStarStars(iMaxExcit,iExcit,RhoEps)
                 
             ENDIF
         ENDIF

         IF(.NOT.BTEST(NWHTAY,0)) THEN
            WRITE(6,*) "Beginning Complete Star Diagonalization"
            IF(StarProd) THEN
                IF(TCalcRealProd) THEN
                    IF(TSumProd) THEN
                        WRITE(6,*) "Product Sum Star Diagonalisation"
                        CALL StarDiag(iExcit+1,ExcitInfo,iMaxExcit+1,i_p,fMCPR3StarNewExcit,dBeta,dLWdB)
                    ELSE
                        WRITE(6,*) "Real Product Star Diagonalisation"
                        CALL StarDiagRealProd(iExcit+1,ExcitInfo,iMaxExcit+1,i_P,fMCPR3StarNewExcit,dBeta, &
                            dLWdB,ProdNum,ProdPositions,OnDiagProdRho,OffDiagProdRho)
                    ENDIF
                ELSE
                    WRITE(6,*) "Complete Product Star Diagonalisation -  beware - large scaling!"
                    CALL StarDiagSC(iExcit+1,ExcitInfo,iMaxExcit+1,i_P,fMCPR3StarNewExcit,dBeta,dLWdB)
                ENDIF
            ELSEIF(TLanczos) THEN
!Lanczos diagonalisation

                WRITE(6,*) "Performing Lanczos diagonalisation of the star graph"
                CALL StarDiagLanc(iExcit+1,ExcitInfo,iMaxExcit+1,i_P,fMCPR3StarNewExcit,dBeta,dLWdB)
            ELSE
                CALL StarDiag(iExcit+1,ExcitInfo,iMaxExcit+1,i_P,fMCPR3StarNewExcit,dBeta,dLWdB)
            ENDIF

        ELSEIF(TMCStar) THEN

            WRITE(6,*) "Performing Monte Carlo diagonalization of the star graph"
            CALL StarDiagMC(iExcit+1,ExcitInfo,iMaxExcit+1,fMCPR3StarNewExcit,dLWdB,MaxDiag)
         
        ELSE

            WRITE(6,*) "Beginning Polynomial Star Diagonalization"
            nRoots=iExcit
            IF(BTEST(NWHTAY,1)) THEN
               nRoots=1
               WRITE(6,*) "Searching for 1 root"
            ELSEIF(BTEST(NWHTAY,2)) THEN
               WRITE(6,*) "Searching for enough roots to converge sum"
               nRoots=iExcit+1
            ELSEIF(BTEST(NWHTAY,6)) THEN
               WRITE(6,*) "Searching for enough roots to converge sum - Method 2"
               nRoots=iExcit+2
            ELSE
               WRITE(6,*) "Searching for all roots"
            ENDIF
            CALL StarDiag2(iExcit+1,ExcitInfo,iMaxExcit+1,i_P,fMCPR3StarNewExcit,dBeta,dLWdB,nRoots)
        ENDIF
         
         IF(ASSOCIATED(ExcitInfo)) THEN
            call LogMemDealloc(this_routine,tagExcitInfo)
            Deallocate(ExcitInfo)
         ENDIF
         IF(ALLOCATED(ExcitInfo2)) THEN
             call LogMemDealloc(this_routine,tagExcitInfo2)
             Deallocate(ExcitInfo2)
         ENDIF
         IF(ALLOCATED(EXCITSTORE)) THEN
            call LogMemDealloc(this_routine,tagEXCITSTORE)
            DEALLOCATE(EXCITSTORE)
         ENDIF
         IF(ALLOCATED(ProdPositions)) THEN
            call LogMemDealloc(this_routine,tagProdPositions)
            DEALLOCATE(ProdPositions)
         ENDIF
         IF(ALLOCATED(OnDiagProdRho)) THEN
            call LogMemDealloc(this_routine,tagOnDiagProdRho)
            DEALLOCATE(OnDiagProdRho)
            call LogMemDealloc(this_routine,tagOffDiagProdRho)
            DEALLOCATE(OffDiagProdRho)
         ENDIF
            
      END function

!This routine explicitly calculates all excited stars, and prediagonalises them.
!A star with its root at an excited Double excitation, can have excitations itself,
!which are single, double, triple and quadruple excitations of the HF determinant.
!To limit the excitations to ones which are Quadruple excitations, TJustQuads must be set.
!To limit the excitations to ones which are all but double excitations, i.e. no
!crosslinking, TNoDoubs must be set.
        SUBROUTINE CalcExcitStar(iMaxExcit,iExcit,nI,rhii,Beta,i_p,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,nTay,ECore,RhoEps) 
            use IntegralsData , only : TQuadValMax,TQuadVecMax,TJustQuads,TNoDoubs
            use SystemData, only: BasisFN
            use global_utilities
            use HElem
            IMPLICIT NONE
            TYPE(BasisFN) G1(*)
            HElement_t :: rhii,UMat(*),rh,rhij,Hij
            complex(dp) :: fck(*)
            real(dp) :: Beta,ALat(3),ECore,RhoEps
            INTEGER :: nEl
            INTEGER :: i,j,iExcit,nI(nEl),i_P,nBasis,nMsh
            INTEGER :: DoublePath(nEl),nStore2(6),exFlag2,iMaxExcit2,nJ(nEl),nExcitMemLen2(1)
            INTEGER :: QuadExcits,iExcit2,TotExcits,NextVertex,NoExcitsInStar(iExcit)
            type(timer), save :: proc_timer
            INTEGER :: nMax,nTay(2),iErr,ICMPDETS,iMaxExcit,temp,IGETEXCITLEVEL
            INTEGER(int32) Info
            INTEGER, ALLOCATABLE :: nExcit2(:)
            real(dp), ALLOCATABLE :: ExcitStarInfo(:,:),ExcitStarMat(:,:),WORK(:)
            real(dp), ALLOCATABLE :: ExcitStarVals(:),ExcitStarVecs(:)
            integer(TagIntType),save :: tagWORK=0,tagExcitStarVals=0,tagExcitStarMat=0
            integer(TagIntType),save :: tagExcitStarInfo=0
            LOGICAL :: HFFound
            character(*), parameter :: this_routine='CalcExcitStar'
            
            WRITE(6,*) "Explicitly calculating and prediagonalising all double excitations from the "&
                & //"original excitations of the star graph"
            
            IF(TJustQuads) THEN
                WRITE(6,*) "Excited stars are only allowed to contain excitations which are quadruple "&
                & //"excitations of the HF determinant"
            ENDIF
            
            proc_timer%timer_name='CalcExcitStar'
            call set_timer(proc_timer)
            IF(HElement_t_size.ne.1) STOP 'Only real orbitals allowed'
            
!Allow only double excitations
            exFlag2=2
            QuadExcits=0
            j=0
            NoExcitsInStar(1:iExcit)=0
            
!Initially, count all excitations
            do i=1,iMaxExcit
                
                IF(ExcitStore(1,i).eq.0) EXIT
                j=j+1
                CALL GetFullPath(nI,nEl,2,ExcitStore(:,i),DoublePath(:))
                nExcitMemLen2=0
                nStore2(1:6)=0
                CALL GenSymExcitIt2(DoublePath,nEl,g1,nBasis,.TRUE.,nExcitMemLen2,nJ,iMaxExcit2,nStore2,exFlag2)
                ALLOCATE(nExcit2(nExcitMemLen2(1)))
                nExcit2(1:nExcitMemLen2(1))=0
                CALL GenSymExcitIt2(DoublePath,nEl,G1,nBasis,.TRUE.,nExcit2,nJ,iMaxExcit2,nStore2,exFlag2)

                lpcount: do while(.true.)
                    CALL GenSymExcitIt2(DoublePath,nEl,G1,nBasis,.false.,nExcit2,nJ,iExcit2,nStore2,exFlag2)
                    IF(nJ(1).eq.0) exit lpcount

!If TNoDoubs is set, we disallow all double excitations of double excitations, which are themselves double excitations of the HF.
                    IF(TNoDoubs) THEN
                        IF(IGetExcitLevel(nI,nJ,nEl).eq.2) THEN
                            CYCLE
                        ENDIF
                    ENDIF
                    
!If TJustQuads is set, we only want to include excitations which are quadruple excitations of the HF in the excited stars.
                    IF(TJustQuads) THEN
!                        WRITE(6,*) IGetExcitLevel(nI,nJ,nEl)
                        IF(IGetExcitLevel(nI,nJ,nEl).ne.4) THEN
                            CYCLE
                        ENDIF
                    ENDIF
                    CALL CalcRho2(DoublePath,nJ,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,iExcit2,ECore)
                    IF(abs(rh).gt.RhoEps) THEN
                        QuadExcits=QuadExcits+1
                        NoExcitsInStar(i)=NoExcitsInStar(i)+1
                    ENDIF
                enddo lpcount
                
!Remove excitations back to HF (if not already removed by TJustQuads
                IF(.not.TJustQuads) THEN 
                    QuadExcits=QuadExcits-1
                    NoExcitsInStar(i)=NoExcitsInStar(i)-1
                ENDIF

!Uninitialise excitation generators
                Deallocate(nExcit2)

            enddo

            temp=0
            do i=1,iExcit
                temp=temp+NoExcitsInStar(i)
            enddo
            IF(temp.ne.QuadExcits) STOP 'Error when counting excitations here'
            IF(j.ne.iExcit) STOP 'Error when counting excitations'
            
!Total number of excitations is equal to the total excitations for the excited stars, + the original double excitations
            TotExcits=iExcit+QuadExcits

            WRITE(6,"(I10,A)") QuadExcits," extra excitations to attach to the HF-rooted star graph"

!Now go through all excitations, finding the excited star, and diagonalising, adding the eigenvalues and vectors to ExcitInfo2(0:TotExcits,0:2)
            ALLOCATE(ExcitInfo2(0:TotExcits,0:2),stat=iErr)
            call LogMemAlloc('ExcitInfo2',3*(TotExcits+1),8*HElement_t_size,this_routine,tagExcitInfo2,iErr)
            ExcitInfo2=(0.0_dp)

!Fill original star matrix - INCORRECT - do not want to include original excitations - these are already included in the prediagonalised elements
!            do i=0,iExcit
!
!                ExcitInfo2(i,0)=ExcitInfo(i,0)
!                ExcitInfo2(i,1)=ExcitInfo(i,1)
!                ExcitInfo2(i,2)=ExcitInfo(i,2)
!
!            enddo
!            NextVertex=iExcit+1

!All that is not included is the original i spoke
            ExcitInfo2(0,0)=(1.0_dp)
            ExcitInfo2(0,1)=(1.0_dp)
            ExcitInfo2(0,2)=ExcitInfo(0,2)
            
            NextVertex=1
            do i=1,iExcit
                
                IF(ExcitStore(1,i).eq.0) STOP 'Problem Here!'
                CALL GetFullPath(nI,nEl,2,ExcitStore(:,i),DoublePath(:))

!Calculate rhi_ij and H_ij from HF to excited star root
                CALL CalcRho2(nI,DoublePath,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,iExcit2,ECore)
                rhij=rh/rhii
                Hij = get_helement (nI, DoublePath, iExcit2)
                
!Reinitialise excitation generators
                HFFound=.false.
                nExcitMemLen2=0
                nStore2(1:6)=0
                CALL GenSymExcitIt2(DoublePath,nEl,G1,nBasis,.TRUE.,nExcitMemLen2,nJ,iMaxExcit2,nStore2,exFlag2)
                ALLOCATE(nExcit2(nExcitMemLen2(1)))
                nExcit2(1:nExcitMemLen2(1))=0
                CALL GenSymExcitIt2(DoublePath,nEl,G1,nBasis,.TRUE.,nExcit2,nJ,iMaxExcit2,nStore2,exFlag2)
                CALL CalcRho2(DoublePath,DoublePath,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
                
!Allocate Memory for excited star.
                ALLOCATE(ExcitStarInfo(0:NoExcitsInStar(i),0:1),stat=iErr)
                CALL LogMemAlloc("ExcitStarInfo",(NoExcitsInStar(i)+1)*2,8,this_routine,tagExcitStarInfo,iErr)
                ExcitStarInfo=0.0_dp
                
                ExcitStarInfo(0,0)=(rh/rhii)
                ExcitStarInfo(0,1)=(rh/rhii)
                j=0
                
                lp: do while(.true.)
                    CALL GenSymExcitIt2(DoublePath,nEl,G1,nBasis,.false.,nExcit2,nJ,iExcit2,nStore2,exFlag2)
                    IF(nJ(1).eq.0) exit lp

!If TNoDoubs is set, we disallow all double excitations of double excitations, which are themselves double excitations of the HF.
                    IF(TNoDoubs) THEN
                        IF(IGetExcitLevel(nI,nJ,nEl).eq.2) THEN
                            CYCLE
                        ENDIF
                    ENDIF
!If TJustQuads is set, we only want to include excitations which are quadruple excitations of the HF in the excited stars.
                    IF(TJustQuads) THEN
                        IF(IGetExcitLevel(nI,nJ,nEl).ne.4) THEN
                            CYCLE
                        ENDIF
                    ENDIF

!Remove excitation to HF generated in excited stars
                    IF((ICMPDETS(nJ,nI,nEl)).eq.0) THEN
                        IF(HFFound) THEN
                            WRITE(6,*) "Error-HF generated twice for Double excitation ",i
                            STOP
                        ENDIF
                        HFFound=.true.
                        CYCLE
                    ENDIF

!Calculate rho_jk for excited stars
                    CALL CalcRho2(DoublePath,nJ,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,iExcit2,ECore)
                    IF(abs(rh).gt.RhoEps) THEN
                        j=j+1
                        ExcitStarInfo(j,1)=(rh/rhii)

!Calculate rho_kk for quadruple excitations
                        CALL CalcRho2(nJ,nJ,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
                        ExcitStarInfo(j,0)=(rh/rhii)

                        IF(j.gt.NoExcitsInStar(i)) STOP 'Incorrect Counting here'
                    ENDIF

                enddo lp
                IF(j.ne.NoExcitsInStar(i)) STOP 'Incorrect Counting here 2'
                IF(.not.HFFound.and..not.TJustQuads) STOP 'HF excitation not generated in excited star'

!Now need to prepare to diagonalise excited star
                ALLOCATE(ExcitStarMat(j+1,j+1),stat=iErr)
                CALL LogMemAlloc("ExcitStarMat",(j+1)*(j+1),8,this_routine,tagExcitStarMat,iErr)
                ExcitStarMat=0.0_dp

                do j=1,(NoExcitsInStar(i)+1)
                    ExcitStarMat(j,j)=ExcitStarInfo(j-1,0)
                    ExcitStarMat(1,j)=ExcitStarInfo(j-1,1)
                enddo

                CALL LogMemDealloc(this_routine,tagExcitStarInfo)
                DEALLOCATE(ExcitStarInfo)

                ALLOCATE(WORK(3*(NoExcitsInStar(i)+1)),stat=iErr)
                CALL LogMemAlloc("WORK",3*(NoExcitsInStar(i)+1),8,this_routine,tagWORK,iErr)
                WORK=0.0_dp

                ALLOCATE(ExcitStarVals(NoExcitsInStar(i)+1),stat=iErr)
                CALL LogMemAlloc("ExcitStarVals",NoExcitsInStar(i)+1,8,this_routine,tagExcitStarVals,iErr)
                ExcitStarVals=0.0_dp

                ALLOCATE(ExcitStarVecs(NoExcitsInStar(i)+1),stat=iErr)
                CALL LogMemAlloc("ExcitStarVecs",NoExcitsInStar(i)+1,8,this_routine,tagExcitStarVals,iErr)
                ExcitStarVecs=0.0_dp

                CALL DSYEV('V','U',NoExcitsInStar(i)+1,ExcitStarMat,NoExcitsInStar(i)+1,ExcitStarVals, &
                    WORK,3*(NoExcitsInStar(i)+1),INFO)
                IF(INFO.ne.0) THEN
                    WRITE(6,*) "DSYEV error in CalcExcitStar: ",INFO
                    STOP
                ENDIF

                do j=1,(NoExcitsInStar(i)+1)
                    ExcitStarVecs(j)=ABS(ExcitStarMat(1,j))
                enddo

                CALL LogMemDealloc(this_routine,tagExcitStarMat)
                DEALLOCATE(ExcitStarMat)
                CALL LogMemDealloc(this_routine,tagWORK)
                DEALLOCATE(WORK)
                
!Now it is necessary to reattach the eigenvectors back to the original star matrix...
                Do j=1,(NoExcitsInStar(i)+1)
                    ExcitInfo2(NextVertex,0)=(ExcitStarVals(j))
                    ExcitInfo2(NextVertex,1)=rhij*(ExcitStarVecs(j))
                    ExcitInfo2(NextVertex,2)=Hij*(ExcitStarVecs(j))
                    NextVertex=NextVertex+1
                enddo

            enddo

            IF(NextVertex.ne.(TotExcits+1)) THEN
                WRITE(6,*) "Next Vertex is: ",NextVertex
                WRITE(6,*) "TotExcits is: ", TotExcits
                STOP 'Incorrect Counting Here 3'
            ENDIF

            iExcit=TotExcits
            iMaxExcit=TotExcits

            Call LogMemDealloc(this_routine,tagExcitInfo)
            DEALLOCATE(ExcitInfo)

            ExcitInfo => ExcitInfo2

            call halt_timer(proc_timer)

            RETURN

        END SUBROUTINE CalcExcitStar
       
!A testing routine which writes out eigenvectors and values from various modified star matrices
        SUBROUTINE GraphRootChange(iMaxExcit,iExcit)
            use global_utilities
            use HElem
            IMPLICIT NONE
            INTEGER :: iMaxExcit,iExcit,HalfiExcit,i,j,calcs,toprint
            INTEGER :: NoDegens,iErr,r
            INTEGER, ALLOCATABLE :: DegenPos(:)
            LOGICAL :: degen
            real(dp) :: gap,minimum
            real(dp), ALLOCATABLE :: Vals(:),Vecs(:),DiagRhos(:)
            HElement_t :: tmp(3)
            integer(TagIntType), save :: tagDegenPos=0
            character(*), parameter :: this_routine='GraphRootChange'

!            OPEN(47,FILE='Differences',STATUS='UNKNOWN')
            OPEN(48,FILE='FirstElemVecs',STATUS='UNKNOWN')
            OPEN(49,FILE='Vals',STATUS='UNKNOWN')
            
!First it is necessary to order the rho_jj elements, so that the range that the linear approximation needs to hold can be worked out.
!This routine sorts into ASCENDING order of rho_jj - therefore rho_jj max = ExcitInfo(iMaxExcit,0) = 1
            ! Sort elements 0:iMaxExcit.
            call sort (excitInfo(:,0), excitInfo(:,1), excitInfo(:,2))

!Reverse order of array ExcitInfo, as have coded up the other way round! - rho_jj max = ExcitInfo(0,0) = 1, and then in decending order.
            HalfiExcit=INT((iMaxExcit+1)/2)
            do i=1,HalfiExcit
                tmp(1)=ExcitInfo((iMaxExcit+1)-i,0)
                tmp(2)=ExcitInfo((iMaxExcit+1)-i,1)
                tmp(3)=ExcitInfo((iMaxExcit+1)-i,2)
                ExcitInfo((iMaxExcit+1)-i,0)=ExcitInfo(i-1,0)
                ExcitInfo((iMaxExcit+1)-i,1)=ExcitInfo(i-1,1)
                ExcitInfo((iMaxExcit+1)-i,2)=ExcitInfo(i-1,2)
                ExcitInfo(i-1,0)=tmp(1)
                ExcitInfo(i-1,1)=tmp(2)
                ExcitInfo(i-1,2)=tmp(3)
            enddo
            
            IF(.not.abs(ExcitInfo(iExcit,0)).gt.0.1_dp) STOP 'Reordering incorrect'
            WRITE(6,*) "Minimum rho_jj is :", ExcitInfo(iExcit,0)

! Ensures that all determinants are non-degenerate - useful for testing
!            do i=1,iExcit
!                ExcitInfo(i,0)=(exp(-0.0001*i))
!            enddo

!It is necessary to find the sets of degenerate rho_jj elements, in order to predict which are non-zero eigenvectors.
!First, it is necessary to count the number of degenerate sets.
            NoDegens=0
            i=0
            do while(i.lt.iExcit)
                i=i+1
                degen=.true.
                do while((REAL(ExcitInfo(i,0),dp).eq.REAL(ExcitInfo(i+1,0),dp)).and.(i.lt.iExcit))
                    i=i+1
                enddo
                NoDegens=NoDegens+1
            enddo
            WRITE(6,*) "Number of degenerate sets of excited determinants = ", NoDegens
            IF(NoDegens.gt.iExcit) THEN
                STOP 'Cannot have more degenerate sets that excitations!'
            ENDIF

            ALLOCATE(DegenPos(NoDegens),stat=iErr)
            CALL LogMemAlloc("DegenPos",NoDegens/2,4,this_routine,tagDegenPos,iErr)
            DegenPos(1:NoDegens)=0
!DegenPos shows the degeneracy structure of the excited determinants
!Each degenerate block extends from DegenPos(i-1)+1 --> DegenPos(i)
!DegenPos only refers to the excited determinants - the root is its own separate degenerate block
!DegenPos(NoDegens) should equal iExcit
            i=0
            j=0
            do while(i.lt.iExcit)
                i=i+1
                do while((REAL(ExcitInfo(i,0),dp).eq.REAL(ExcitInfo(i+1,0),dp)).and.(i.lt.iExcit))
                    i=i+1
                enddo
                j=j+1
                DegenPos(j)=i
!                WRITE(6,*) i
            enddo
!            WRITE(6,*) NoDegens,DegenPos(:)

!            do i=1,iExcit+1
!                WRITE(16,*) ExcitInfo(i-1,0),ExcitInfo(i-1,1)
!            enddo
            IF(DegenPos(NoDegens).ne.iExcit) THEN
                STOP 'Final element of DegenPos should equal iExcit to account for all degeneracies'
            ENDIF
            
            ALLOCATE(Vals(iExcit+1))
            ALLOCATE(Vecs(iExcit+1))
            ALLOCATE(DiagRhos(iExcit+1))
            Vals=0.0_dp
            Vecs=0.0_dp
            DiagRhos=0.0_dp

            calcs=100
            minimum=1.0_dp-((1.0_dp-REAL(ExcitInfo(iExcit,0),dp))*2)
            WRITE(6,*) "Minimum root value to search for is: ", minimum
            CALL neci_flush(6)
            gap=(1.0_dp-minimum)/(calcs-1)
            
!This is the number of non-degenerate eigenvectors to print out.            
            toprint=25
            IF(toprint.gt.(NoDegens+1)) THEN
                WRITE(6,*) 'Trying to print more eigenvectors than non-degenerate determinants'
                WRITE(6,*) "Resetting number to print out to ",NoDegens+1
                toprint=NoDegens+1
            ENDIF

!Diffs can be used to calculate the gradient of the eigenvector/root line quickly
!            ALLOCATE(Diffs(toprint,4))
!            Diffs=0.0_dp
!Store value for eigenvalue/vector at largest rhovalue in 1:2 - the lowest rho value in 3:4

!            ALLOCATE(TESTER(iExcit))
!            do i=1,iExcit
!                TESTER(i)=(0.0_dp)
!            enddo
            
            j=0
            do r=1,int(minimum),int(-gap)
                j=j+1
!Fill matrix as normal, but change root element from 1 -> little lower than rho_jj
                do i=2,iExcit+1
!Comment out as appropriate (all elements multiplied by r, or just root
!                    DiagRhos(i)=(ExcitInfo(i-1,0))
                    DiagRhos(i)=(ExcitInfo(i-1,0))*r
                enddo
                DiagRhos(1)=r

                WRITE(6,*) "Running calculation ",j, "out of ",calcs
                CALL neci_flush(6)
                CALL GetValsnVecs(iExcit+1,DiagRhos,ExcitInfo(1:iExcit,1),Vals,Vecs)
!                CALL GetValsnVecs(iExcit+1,DiagRhos,TESTER,Vals,Vecs)
                
                IF(r.eq.1.0_dp) THEN
                    WRITE(6,"(A,F15.11)") "For root equal to 1, highest eigenvalue is ", Vals(iExcit+1)
                    IF(j.ne.1) STOP 'Problem with counting'
                ENDIF

!                Write(6,*) Vals(1)
!Quick way of determining gradient, assuming linear relationship (look at difference over rho values
!                IF(j.eq.1) THEN
!                    Diffs(1,1)=Vals(iExcit+1)
!                    Diffs(1,2)=Vecs(iExcit+1)
!                    do i=2,toprint
!                        Diffs(i,1)=Vals(iExcit+1-DegenPos(i-1))
!                        Diffs(i,2)=Vecs(iExcit+1-DegenPos(i-1))
!                    enddo
!                ELSEIF(j.eq.calcs) THEN
!                    Diffs(1,3)=Vals(iExcit+1)
!                    Diffs(1,4)=Vecs(iExcit+1)
!                    do i=2,toprint
!                        Diffs(i,3)=Vals(iExcit+1-DegenPos(i-1))
!                        Diffs(i,4)=Vecs(iExcit+1-DegenPos(i-1))
!                    enddo
!                ENDIF
                        
!Just write out non-degenerate eigenvectors
                WRITE(48,"(F11.7)",advance='no') r
                WRITE(48,"(F13.9)",advance='no') Vecs(iExcit+1)
                do i=2,toprint
                    WRITE(48,"(F13.9)",advance='no') Vecs(iExcit+1-DegenPos(i-1))
                enddo
                WRITE(48,*) ""
                WRITE(49,"(F11.7)",advance='no') r
                WRITE(49,"(F13.9)",advance='no') Vals(iExcit+1)
                do i=2,toprint
                    WRITE(49,"(F13.9)",advance='no') Vals(iExcit+1-DegenPos(i-1))
                enddo
                WRITE(49,*) ""
                CALL neci_flush(48)
                CALL neci_flush(49)

!To write out all possible eigenvalues/vectors
!                WRITE(48,"(F11.7)",advance='no') r
!                do i=1,iExcit+1
!                    WRITE(48,"(F13.9)",advance='no') Vecs(i)
!                enddo
!                WRITE(48,*) ""
!                WRITE(49,"(F11.7)",advance='no') r
!                do i=1,iExcit+1
!                    WRITE(49,"(F13.9)",advance='no') Vals(i)
!                enddo
!                WRITE(49,*) ""
!                CALL neci_flush(48)
!                CALL neci_flush(49)

            enddo

!            WRITE(47,*) "1",Diffs(1,2)-Diffs(1,4),Diffs(1,1)-Diffs(1,3)
!            do i=2,toprint
!                WRITE(47,*) DegenPos(i-1)+1,Diffs(i,2)-Diffs(i,4),Diffs(i,1)-Diffs(i,3)
!            enddo

!Write out gnuscript for only non-degenerate eigenvectors
            OPEN(26,FILE='PlotDegen.gpi',STATUS='UNKNOWN')
            WRITE(26,*) "set key left"
            WRITE(26,"(A)",advance='no') "plot 'FirstElemVecs' u 1:(abs($2)) w lp t 'Vec 1', "
            do i=2,toprint-1
                WRITE(26,"(A,I3,A,I3,A)",advance='no') "'' u 1:(abs($",i+1,")) w lp t 'Vec",DegenPos(i-1)+1,"', "
            enddo
            WRITE(26,"(A,I3,A,I3,A)") "'' u 1:(abs($",toprint+1,")) w lp t 'Vec",DegenPos(toprint-1)+1,"'"
            WRITE(26,*) "pause -1"
            WRITE(26,"(A)",advance='no') "plot 'Vals' u 1:(abs($2)) w lp t 'Val 1', "
            do i=2,toprint-1
                WRITE(26,"(A,I3,A,I3,A)",advance='no') "'' u 1:(abs($",i+1,")) w lp t 'Val",DegenPos(i-1)+1,"',"
            enddo
            WRITE(26,"(A,I3,A,I3,A)") "'' u 1:(abs($",toprint+1,")) w lp t 'Val",DegenPos(toprint-1)+1,"'"
!            WRITE(26,*) "pause -1"
!            WRITE(26,*) "set xlabel 'Non-zero Eigenvector number'"
!            WRITE(26,*) "plot 'Differences' u :(abs($2)) w lp t 'Difference in Eigenvector', '' u :(abs($3)) w lp t 'Difference in Eigenvalue'"
            CLOSE(26)
!
!Write out gnuscript for all eigenvectors
!            OPEN(26,FILE='Plotall.gpi',STATUS='UNKNOWN')
!            WRITE(26,*) "set key left"
!            WRITE(26,"(A,I3,A)",advance='no') "plot 'FirstElemVecs' u 1:(abs($",iExcit+2,")) w lp t 'Vector 1',"
!            do i=iExcit+1,3,-1
!                WRITE(26,"(A,I3,A,I3,A)",advance='no') "'' u 1:(abs($",i,")) w lp t 'Vec",iExcit+3-i,"', "
!            enddo
!            WRITE(26,"(A,I3,A)") "'' u 1:(abs($2)) w lp t 'Vec",iExcit+1,"'"
!            WRITE(26,*) "pause -1"
!            WRITE(26,"(A,I3,A)",advance='no') "plot 'Vals' u 1:(abs($",iExcit+2,")) w lp t 'Val 1',"
!            do i=iExcit+1,3,-1
!                WRITE(26,"(A,I3,A,I3,A)",advance='no') "'' u 1:(abs($",i,")) w lp t 'Val",iExcit+3-i,"', "
!            enddo
!            WRITE(26,"(A,I3,A)") "'' u 1:(abs($2)) w lp t 'Val",iExcit+1,"'"
!            CLOSE(26)

            CALL LogMemDealloc(this_routine,tagDegenPos)
            DEALLOCATE(DegenPos)
            DEALLOCATE(Vals)
            DEALLOCATE(Vecs)
            DEALLOCATE(DiagRhos)

        END SUBROUTINE GraphRootChange

        
!GetStarStars approximates excited stars as having the same connections as the original star, and so simply multiplies the diagonal elements by rho_jj and then diagonalises them.
        SUBROUTINE GetStarStars(iMaxExcit,iExcit,RhoEps)
            use IntegralsData , only : TExcitStarsRootChange,TRmRootExcitStarsRootChange
            use global_utilities
            use HElem
            IMPLICIT NONE
            INTEGER :: NextVertex,i,j,iErr,iMaxExcit,iExcit
            type(timer), save :: proc_timer
            real(dp), ALLOCATABLE :: NewDiagRhos(:),Vals(:),Vecs(:)
            HElement_t, ALLOCATABLE :: NewOffDiagRhos(:)
            integer(TagIntType), save :: tagVals=0,tagVecs=0,tagNewDiagRhos=0,tagNewOffDiagRhos=0
            real(dp) :: RhoValue,RhoEps,OffRhoValue
            HElement_t :: Rhoia
            LOGICAL :: FoundRoot
            character(*), parameter :: this_routine='GetStarStars'

            proc_timer%timer_name='GetStarStars'
            call set_timer(proc_timer)
            IF(HElement_t_size.ne.1) STOP 'Only real orbitals allowed in GetStarStars'

            WRITE(6,*) "Explicitly diagonalising approximate excited stars from HF star template"
            IF(TExcitStarsRootChange) THEN
                WRITE(6,*) "Only changing root matrix element for excited stars"
            ENDIF

            ALLOCATE(NewDiagRhos(iExcit+1),stat=iErr)
            CALL LogMemAlloc("NewDiagRhos",iExcit+1,8,this_routine,tagNewDiagRhos,iErr)
            ALLOCATE(NewOffDiagRhos(iExcit),stat=iErr)
            CALL LogMemAlloc("NewOffDiagRhos",iExcit,8*HElement_t_size,this_routine,tagNewOffDiagRhos,iErr)
            NewOffDiagRhos=(0.0_dp)
            ALLOCATE(ExcitInfo2(0:iExcit*(iExcit+1),0:2),stat=iErr)
            call LogMemAlloc('ExcitInfo2',(iExcit*(iExcit+1)+1)*3,8*HElement_t_size,this_routine,tagExcitInfo2,iErr)
            ExcitInfo2=(0.0_dp)
            iMaxExcit=iExcit*(iExcit+1)
            
            ExcitInfo2(0,0)=(1.0_dp)
            ExcitInfo2(0,1)=(1.0_dp)
            ExcitInfo2(0,2)=ExcitInfo(0,2)
            NextVertex=1
            
            ALLOCATE(Vals(iExcit+1),stat=iErr)
            CALL LogMemAlloc("Vals",(iExcit+1)*iExcit,8,this_routine,tagVals,iErr)
            ALLOCATE(Vecs(iExcit+1),stat=iErr)
            CALL LogMemAlloc("Vecs",(iExcit+1)*iExcit,8,this_routine,tagVecs,iErr)

!Run through all excitations of original star
            do i=1,iExcit

                FoundRoot=.false.
!Refill offdiagonal elements the same each time
                do j=1,iExcit
                    NewOffDiagRhos(j)=ExcitInfo(j,1)
                enddo

                NewDiagRhos=0.0_dp
                Vals=0.0_dp
                Vecs=0.0_dp
                RhoValue=(ExcitInfo(i,0))
                OffRhoValue=(ExcitInfo(i,1))
                
!Fill matrix of excited star to prediagonalise
                NewDiagRhos(1)=RhoValue
                do j=1,iExcit
                    IF(TExcitStarsRootChange) THEN
!Only change the root element for the excited star matrix. 
                        NewDiagRhos(j+1)=(ExcitInfo(j,0))
                    
!Remove connection to itself in the excited star.
                    ELSEIF(TRmRootExcitStarsRootChange) THEN
                        IF((REAL(ExcitInfo(j,0),dp).eq.RhoValue).and.(OffRhoValue.eq.REAL(ExcitInfo(j,1),dp)) &
                            .and.(.not.FoundRoot)) THEN
                            NewDiagRhos(j+1)=0.0_dp
                            NewOffDiagRhos(j+1)=(0.0_dp)
                            FoundRoot=.true.
                        ELSE
                            NewDiagRhos(j+1)=(ExcitInfo(j,0))
                        ENDIF

!Multiply all diagonal elements by current rho_jj value
                    ELSE
                        NewDiagRhos(j+1)=(ExcitInfo(j,0))*RhoValue
                    ENDIF
                enddo

                IF(TRmRootExcitStarsRootChange.and..not.FoundRoot) THEN
                    STOP 'Could not find the root in excited star to remove'
                ENDIF

!Diagonalise
                CALL GetValsnVecs(iExcit+1,NewDiagRhos,NewOffDiagRhos,Vals,Vecs)
                
                do j=1,iExcit+1

                    Rhoia=ExcitInfo(i,1)*(Vecs(j))

                    IF(abs(Rhoia).gt.RhoEps) THEN
!Return excited star values to ExcitInfo2
                        ExcitInfo2(NextVertex,0)=(Vals(j))
                        ExcitInfo2(NextVertex,1)=Rhoia
                        ExcitInfo2(NextVertex,2)=ExcitInfo(i,2)*(Vecs(j))
                        NextVertex=NextVertex+1
                    ENDIF

                enddo
                
            enddo

            CALL LogMemDealloc(this_routine,tagNewDiagRhos)
            DEALLOCATE(NewDiagRhos)
            CALL LogMemDealloc(this_routine,tagVals)
            DEALLOCATE(Vals)
            CALL LogMemDealloc(this_routine,tagVecs)
            DEALLOCATE(Vecs)

            WRITE(6,"(I10,A)") NextVertex-1-iExcit, " extra vertices added to original star from excited stars"

            iExcit=NextVertex-1

            Call LogMemDealloc(this_routine,tagExcitInfo)
            DEALLOCATE(ExcitInfo)

            ExcitInfo => ExcitInfo2

            call halt_timer(proc_timer)

            RETURN

        END SUBROUTINE GetStarStars

        SUBROUTINE GetLinRootChangeStars(iMaxExcit,iExcit,nWHTay)
            use CalcData , only : LinePoints
            use global_utilities
            use HElem
            IMPLICIT NONE
            INTEGER :: i,j,iMaxExcit,iExcit,nWHTay,HalfiExcit,ierr
            type(timer), save :: proc_timer
            real(dp) :: LineRhoValues(LinePoints),RhoValue,Vals(LinePoints),meanx,RhoGap,EigenMax
            real(dp) :: MeanVal,Sxx,Sxy,Syy,GradVal,IncptVal,ExpctVal,Rsq,Vector,PreVec,lowrhojj
            LOGICAL :: ReachMax
            real(dp), ALLOCATABLE :: AllVals(:),AllVecs(:),DiagRhos(:)
            integer(TagIntType), save :: tagAllVals=0,tagAllVecs=0,tagDiagRhos=0
            HElement_t :: tmp(3)
            character(*), parameter :: this_routine='GetLinRootChangeStars'

            proc_timer%timer_name='GetLinRootChangeStars'
            call set_timer(proc_timer)
            IF(HElement_t_size.ne.1) STOP 'Only real orbitals allowed'
            WRITE(6,*) "Stars where only root changes to be included, using a linear approximation of eigenvalues"
            
!First it is necessary to order the rho_jj elements, so that the range that the linear approximation needs to hold can be worked out.
!This routine sorts into ASCENDING order of rho_jj - therefore rho_jj max = ExcitInfo(iMaxExcit,0) = 1
            call sort (excitInfo(:,0), excitInfo(:,1), excitInfo(:,2))

!Reverse order of array ExcitInfo, as have coded up the other way round! - rho_jj max = ExcitInfo(0,0) = 1 - rho_jj elements then decrease
            HalfiExcit=INT((iMaxExcit+1)/2)
            do i=1,HalfiExcit
                tmp(1)=ExcitInfo((iMaxExcit+1)-i,0)
                tmp(2)=ExcitInfo((iMaxExcit+1)-i,1)
                tmp(3)=ExcitInfo((iMaxExcit+1)-i,2)
                ExcitInfo((iMaxExcit+1)-i,0)=ExcitInfo(i-1,0)
                ExcitInfo((iMaxExcit+1)-i,1)=ExcitInfo(i-1,1)
                ExcitInfo((iMaxExcit+1)-i,2)=ExcitInfo(i-1,2)
                ExcitInfo(i-1,0)=tmp(1)
                ExcitInfo(i-1,1)=tmp(2)
                ExcitInfo(i-1,2)=tmp(3)
            enddo

            WRITE(6,*) "Total number of points from which to form linear approximation = ", LinePoints
            IF(LinePoints.lt.2) STOP 'LinePoints cannot be less than two'

!Take 'Linepoints' points along the change in rho_jj to calculate the gradient of the line for each eigenvalue, and the first element of the eigenvectors.
!Only two points are strictly needed, but 'Linepoints' will be taken so that the validity of the linear approximation can be calculated.
!Assign largest diagonal multiplicative constant to simply be the original star graph, i.e. rho_ii/rho_ii is the root
            IF((ABS(REAL(ExcitInfo(0,0),dp)-1.0_dp)).gt.1.0e-7_dp) THEN
                STOP 'First element of original star matrix should equal 1.0_dp'
            ENDIF
            LineRhoValues(1)=1.0_dp

!ExcitInfo(iExcit,0) is the smallest rho_jj value.
            WRITE(6,*) "Smallest rho_jj value is ", REAL(ExcitInfo(iExcit,0),dp)

!Calculate the number of eigenvectors which contribute to the excited star with the smallest root...
!First calculate the contribution from largest eigenvector...
            ALLOCATE(AllVals(iExcit+1),stat=ierr)
            CALL LogMemAlloc("AllVals",iExcit+1,8,this_routine,tagAllVals,iErr)
            ALLOCATE(AllVecs(iExcit+1),stat=ierr)
            CALL LogMemAlloc("AllVecs",iExcit+1,8,this_routine,tagAllVecs,iErr)
            ALLOCATE(DiagRhos(iExcit+1),stat=ierr)
            CALL LogMemAlloc("DiagRhos",iExcit+1,8,this_routine,tagDiagRhos,iErr)
            do j=2,iExcit+1
                DiagRhos(j)=REAL(ExcitInfo(j-1,0))
            enddo
            DiagRhos(1)=REAL(ExcitInfo(iExcit,0))
            CALL GetValsnVecs(iExcit+1,DiagRhos,ExcitInfo(1:iExcit,1),AllVals,AllVecs)
            Vector=AllVecs(iExcit+1)
            i=0
            reachmax = .false.
            do while((Vector.gt.0.1).or..not.Reachmax)
                i=i+1
                PreVec=Vector
                Vector=AllVecs(iExcit-(i-1))
                IF(PreVec.lt.Vector) THEN
!Largest contributing eigenvector has not yet been reached
                    ReachMax=.false.
                ELSE
                    ReachMax=.true.
                ENDIF
            enddo

            WRITE(6,*) i, "Eigenvectors needed to ensure complete contribution from smallest " &
                & //"root excited star, with eigenvector cutoff of 0.1"
            
!Try fitting not accros whole range of rho_jj values, but just the highest values - closer linear relationship
!            lowerrhojj=INT(iExcit/50)
!            RhoGap=(1.0_dp-REAL(ExcitInfo(lowerrhojj,0)))/(LinePoints-1)
!Instead of fitting to a certain distance down the list of rho_jjs, pick a minimum rootchange value directly from the spread.
             lowrhojj=1.0_dp-((1.0_dp-REAL(ExcitInfo(iExcit,0),dp))/50)
             RhoGap=(1.0_dp-lowrhojj)/(LinePoints-1)

!Calculate the spread of Rho_jj values which the linear approximation will be based around. Initially, this is just a linear spread.
            do i=2,LinePoints
                LineRhoValues(i)=LineRhoValues(i-1)-RhoGap
            enddo

!LineRhoValues(LinePoints) should be the same as ExcitInfo(iExcit,0)
!            IF(ABS(LineRhoValues(LinePoints)-(ExcitInfo(iExcit,0))).gt.1.0e-9_dp) THEN
!                STOP 'LineRhoValues(LinePoints) should be the same as the lowest rho_jj value'
!            ENDIF
            IF(.NOT.BTEST(NWHTAY,0)) THEN
                ALLOCATE(AllVals(iExcit+1),stat=ierr)
                CALL LogMemAlloc("AllVals",iExcit+1,8,this_routine,tagAllVals,iErr)
                ALLOCATE(AllVecs(iExcit+1),stat=ierr)
                CALL LogMemAlloc("AllVecs",iExcit+1,8,this_routine,tagAllVecs,iErr)
                ALLOCATE(DiagRhos(iExcit+1),stat=ierr)
                CALL LogMemAlloc("DiagRhos",iExcit+1,8,this_routine,tagDiagRhos,iErr)
            ENDIF
!Cycle through all possible roots
            do i=1,LinePoints
                
                RhoValue=LineRhoValues(i)

                write(6,*) ""
                write(6,*) "Rho value: ", RhoValue

!This indicates that a full diagonalisation should be done to calculate eigenvalues.
                IF(.NOT.BTEST(NWHTAY,0)) THEN

!First diagonal element is RhoValue. The other diagonal elements are given by ExcitInfo(1:iExcit,0)
                    do j=2,iExcit+1
                        DiagRhos(j)=(ExcitInfo(j-1,0))
                    enddo
                    DiagRhos(1)=RhoValue
                    
!Find the values for eigenvalues and eigenvectors of this matrix, and put them into the relevant AlllVals and AllVecs
                    CALL GetValsnVecs(iExcit+1,DiagRhos,ExcitInfo(1:iExcit,1),AllVals,AllVecs)

                    IF((RhoValue.eq.1.0_dp).and.(i.ne.1)) STOP 'Problem with assigning rho values for excited stars'

!Save largest eigenvalue
                    IF(RhoValue.eq.1.0_dp) THEN
                        EigenMax=AllVals(iExcit+1)
                        WRITE(6,*) "Largest eigenvalue of original star: ", EigenMax
                    ENDIF

!Initially, assume that all excited stars are only attached by their *first* eigenvectors
!Store the first eigenvalues for the excited stars in Vals. It is these that we assume linearly change.
                    Vals(i)=AllVals(iExcit+1)
                    IF(AllVecs(iExcit+1).lt.0.5) THEN
                        WRITE(6,*) "Eigenvector is equal to ",AllVecs(iExcit+1)," for Linepoint ",i," out of ",LinePoints
                        WRITE(6,*) "Assumption that only one eigenvector from each excited star attached is poor"
                    ENDIF

                    IF(i.eq.Linepoints) THEN
                        CALL LogMemDealloc(this_routine,tagAllVals)
                        DEALLOCATE(AllVals)
                        CALL LogMemDealloc(this_routine,tagAllVecs)
                        DEALLOCATE(AllVecs)
                        CALL LogMemDealloc(this_routine,tagDiagRhos)
                        DEALLOCATE(DiagRhos)
                    ENDIF

                ELSE
                    STOP 'Polynomial calculation of LinRootChangeStars not available yet'

                ENDIF

            enddo
            
!Use linear regression technique to find the best-fit gradient for the line for the eigenvalues

!Mean x value is simply the average of the roots chosen
            meanx=0.0_dp
            do i=1,LinePoints
                meanx=meanx+LineRhoValues(i)
            enddo
            meanx=meanx/LinePoints

!Now need to calculate average eigenvalues
            Meanval=0.0_dp
            do i=1,LinePoints
                Meanval=Meanval+Vals(i)
            enddo
            Meanval=Meanval/Linepoints

!Calculate Sxx
            Sxx=0.0_dp
            do i=1,LinePoints
                Sxx=Sxx+(LineRhoValues(i)-meanx)*(LineRhoValues(i)-meanx)
            enddo
    
!Calculate Sxy and Syy (to find R^2)
!Sxy = \sum{(x_i - X)(y_i - Y) where X and Y are the means of x and y respectivly

            Sxy=0.0_dp
            Syy=0.0_dp
            do i=1,LinePoints
                Sxy=Sxy+(LineRhoValues(i)-meanx)*(Vals(i)-Meanval)
                Syy=Syy+(Vals(i)-Meanval)*(Vals(i)-Meanval)
            enddo

!Best-fit gradients can be calculated for the first eigenvalues - Sxy/Sxx
!y-intercepts also calculated for use when calculating R^2 value.
            GradVal=Sxy/Sxx
            IncptVal=Meanval-GradVal*meanx

!To calculate R^2, also need \sum{(y_i - Y)^2} - the expected y value from the gradient calculation at each point
            ExpctVal=0.0_dp
            do i=1,LinePoints
                ExpctVal=ExpctVal+(Vals(i)-(IncptVal+GradVal*LineRhoValues(i)))**2
            enddo

!Calculate Rsq value
            Rsq=1.0_dp-ExpctVal/Syy

            IF(Rsq.lt.0.95) THEN
                WRITE(6,*) "Problem with linear approximation of eigenvalues, R^2 value: ", Rsq
                STOP
            ELSEIF(Rsq.gt.1.0_dp) THEN
                WRITE(6,*) "Fatal problem - R^2 value greater than 1!"
                STOP
            ENDIF

!Linearly change diagonal elements - rho_jj' = GradVal*(rho_jj - 1) + eigenmax
            ALLOCATE(ExcitInfo2(0:iExcit,0:2),stat=ierr)
            call LogMemAlloc('ExcitInfo2',3*(iExcit+1),8*HElement_t_size,this_routine,tagExcitInfo2,iErr)
            ExcitInfo2=(0.0_dp)
            
!Put HF determinant into element 0
            ExcitInfo2(0,0)=(1.0_dp)
            ExcitInfo2(0,1)=(1.0_dp)
            ExcitInfo2(0,2)=ExcitInfo(0,2)

!Diagonal elements and Hamiltonian elements do not change, since assume eigenvectors are 1 or 0...
            do i=1,iExcit
                ExcitInfo2(i,1)=ExcitInfo(i,1)
                ExcitInfo2(i,2)=ExcitInfo(i,2)
                ExcitInfo2(i,0)=EigenMax+(GradVal*((ExcitInfo(i,0))-1.0_dp))
            enddo

            CALL LogMemDealloc(this_routine,tagExcitInfo)
            DEALLOCATE(ExcitInfo)

            ExcitInfo => ExcitInfo2

            call halt_timer(proc_timer)

            RETURN
        END SUBROUTINE GetLinRootChangeStars

!GetLinStarStars uses the approximation that the change in the eigenvalues and the first element  
!of the eigenvectors, with respect to multiplying the diagonal elements by a constant, is linear. 

!If we try to prediagonalise a star, whose root is a double excitation determinant contained in 
!the original star graph, then all that is needed are the eigenvalues and first elements of the 
!eigenvectors of the original star, as well as how they change with diagonal constant (which will 
!simply be rho_jj/rho_ii for each excited star).

!Problems include the fact that if the set of excitations from the original star is applied to 
!excited stars, then at the level of double excitations, of the (N 2)(M 2) excitations, 
!(N-2 2)(M-2 2) excitations will be real, but the rest will be ficticious, applying annihilation 
!operators to empty orbitals, and creation operators to filled ones.

!Another problem is scaling: to find all the eigenvalues and eigenvectors of the original matrix 
!is an (N^2 M^2)^2 operation. Aside from this, the objects from the diagonalised excited stars 
!have to be attached back to the original star rooted at i, and solved - again potentially an 
!(N^2 M^2)^2 operation, though this could be reduced by only including the largest excitations 
!from the excited star, or MC.

!This could be extended to stars with roots in higher excitation space than quadruple excitations, 
!by simply multiplying the diagonal elements of the original star graph by products of rho_jj 
!elements. This however will make the final star graph even harder to solve. Resumming in the 
!contributions from excited stars using this linear approximation is the hope.

        SUBROUTINE GetLinStarStars(iMaxExcit,iExcit,RhoEps)
            use CalcData , only : LinePoints
            use IntegralsData , only : TQuadValMax,TQuadVecMax
            use global_utilities
            use HElem
            IMPLICIT NONE
            INTEGER :: i,j,iErr,CSE,NextVertex,iMaxExcit
            type(timer), save :: proc_timer
            INTEGER :: iExcit,TotExcits,HalfiExcit
            real(dp) :: tmp(3)
            LOGICAL :: TVal
            real(dp) :: RhoGap,RhoValue,RhoEps,Rhoia,StarEigens(iExcit+1,2)
            real(dp) :: meanx,Sxx,ValMax,VecMax,HMax,Rhoaa

!LineRhoValues gives the values of the diagonal elements multiplicative constant, over which the gradient of the linear approximation will be calculated.
            real(dp), ALLOCATABLE :: LineRhoValues(:)

!VecsDODMS is the array to hold the first elements of the eigenvectors from the diagonalisation of 
!the star matrix with different multiplicative factors down its diagonal. 'Vectors from 
!Diagonalisation Of Diagonally Multiplied Stars'. ValsDODMS is the same to hold the eigenvalues.
            real(dp), ALLOCATABLE :: VecsDODMS(:,:)
            real(dp), ALLOCATABLE :: ValsDODMS(:,:)

!NewDiagRhos holds the diagonal rho elements after they have been multiplied by the constant
            real(dp), ALLOCATABLE :: NewDiagRhos(:)

!These arrays hold various data for calculating the gradient, and R^2 value for the linear approximation for each eigenvalue & vector.
            real(dp), ALLOCATABLE :: RsqVals(:),RsqVecs(:),ExpctVals(:),ExpctVecs(:),IncptVals(:)
            real(dp), ALLOCATABLE :: IncptVecs(:),SxyVals(:),SxyVecs(:)
            real(dp), ALLOCATABLE :: SyyVals(:),SyyVecs(:),MeanVals(:),MeanVecs(:),GradVals(:),GradVecs(:)
            integer(TagIntType), save :: tagVecsDODMS=0
            integer(TagIntType), save :: tagValsDODMS=0
            integer(TagIntType), save :: tagRsqVals=0,tagRsqVecs=0,tagExpctVals=0,tagExpctVecs=0,tagIncptVals=0
            integer(TagIntType), save :: tagIncptVecs=0,tagSxyVals=0,tagSxyVecs=0
            integer(TagIntType), save :: tagSyyVecs=0,tagSyyVals=0,tagMeanVals=0,tagMeanVecs=0,tagGradVals=0,tagGradVecs=0
            integer(TagIntType), save :: tagNewDiagRhos=0
            character(*), parameter :: this_routine='GetLinStarStars'

            proc_timer%timer_name='GetLinStarStars'
            call set_timer(proc_timer)
            IF(HElement_t_size.ne.1) STOP 'Only real orbitals allowed'
            WRITE(6,*) "Stars of double excitations to be included in calculation using a linear approximation of eigenvalues"
            WRITE(6,*) iExcit*(iExcit+1)," possible extra excitations"
            
!First it is necessary to order the rho_jj elements, so that the range that the linear approximation needs to hold can be worked out.
!This routine sorts into ASCENDING order of rho_jj - therefore rho_jj max = ExcitInfo(iMaxExcit,0) = 1
            call sort (excitInfo(:,0), excitInfo(:,1), excitInfo(:,1))

!Reverse order of array ExcitInfo, as have coded up the other way round! - rho_jj max = ExcitInfo(0,0) = 1
            HalfiExcit=INT((iMaxExcit+1)/2)
            do i=1,HalfiExcit
                tmp(1)=ExcitInfo((iMaxExcit+1)-i,0)
                tmp(2)=ExcitInfo((iMaxExcit+1)-i,1)
                tmp(3)=ExcitInfo((iMaxExcit+1)-i,2)
                ExcitInfo((iMaxExcit+1)-i,0)=ExcitInfo(i-1,0)
                ExcitInfo((iMaxExcit+1)-i,1)=ExcitInfo(i-1,1)
                ExcitInfo((iMaxExcit+1)-i,2)=ExcitInfo(i-1,2)
                ExcitInfo(i-1,0)=tmp(1)
                ExcitInfo(i-1,1)=tmp(2)
                ExcitInfo(i-1,2)=tmp(3)
            enddo

            WRITE(6,*) "Total number of points from which to form linear approximation = ", LinePoints

!Take 'Linepoints' points along the change in rho_jj to calculate the gradient of the line for each eigenvalue, and the first element of the eigenvectors.
!Only two points are strictly needed, but 'Linepoints' will be taken so that the validity of the linear approximation can be calculated.
            ALLOCATE(LineRhoValues(LinePoints))

!Assign largest diagonal multiplicative constant to simply be the original star graph, i.e. rho_ii/rho_ii is the root
            IF((ABS(REAL(ExcitInfo(0,0),dp)-1.0_dp)).gt.1.0e-7_dp) THEN
                STOP 'First element of original star matrix should equal 1.0_dp'
            ENDIF
            LineRhoValues(1)=1.0_dp

!ExcitInfo(iExcit,0) is the smallest rho_jj value.
            WRITE(6,*) "Smallest rho_jj value is ", REAL(ExcitInfo(iExcit,0),dp)
            RhoGap=(1.0_dp-ExcitInfo(iExcit,0))/(LinePoints-1)

!Calculate the spread of Rho_jj values which the linear approximation will be based around. Initially, this is just a linear spread.
            do i=2,LinePoints
                LineRhoValues(i)=LineRhoValues(i-1)-RhoGap
            enddo

!LineRhoValues(LinePoints) should be the same as ExcitInfo(iExcit,0)
            IF(ABS(LineRhoValues(LinePoints)-REAL(ExcitInfo(iExcit,0),dp)).gt.1.0e-7_dp) THEN
                STOP 'LineRhoValues(LinePoints) should be the same as the lowest rho_jj value'
            ENDIF

!The values for the calculated eigenvectors (first elements) and eigenvalues for the 'LinePoints' points are put into ValsDODMS and VecsDODMS respectivly.
            ALLOCATE(ValsDODMS(LinePoints,(iExcit+1)),stat=ierr)
            CALL LogMemAlloc("ValsDODMS",(iExcit+1)*LinePoints,8,this_routine,tagValsDODMS,iErr)
            ALLOCATE(VecsDODMS(Linepoints,(iExcit+1)),stat=ierr)
            CALL LogMemAlloc("VecsDODMS",(iExcit+1)*LinePoints,8,this_routine,tagVecsDODMS,iErr)
            ALLOCATE(NewDiagRhos(iExcit+1),stat=ierr)
            CALL LogMemAlloc("NewDiagRhos",iExcit+1,8,this_routine,tagNewDiagRhos,iErr)

!Cycle through all the diagonal multiplicative factors to look at to find gradients
            do i=1,LinePoints

                RhoValue=LineRhoValues(i)

!                write(6,*) ""
!                write(6,*) "Rho value: ", RhoValue

!Multiply the diagonal elements by the value of rho_jj we want                
                do j=1,iExcit+1
                    IF(REAL(ExcitInfo(j-1,0),dp).lt.0.8) THEN
                        STOP 'rho_jj value too small, or incorrect for linear approximation'
                    ENDIF
                    NewDiagRhos(j)=(ExcitInfo(j-1,0))*RhoValue
                enddo

!Find the values for eigenvalues and eigenvectors of this matrix, and put them into the relevant ValsDODMS and VecsDODMS
                CALL GetValsnVecs(iExcit+1,NewDiagRhos,ExcitInfo(1:iExcit,1),ValsDODMS(i,:),VecsDODMS(i,:))

!                WRITE(6,*) "Eigenvalues: ", ValsDODMS(i,:)
!To check linear approximation...
!                WRITE(18,*) RhoValue, ValsDODMS(i,iExcit+1), VecsDODMS(i,iExcit+1)
                
            enddo

            CALL LogMemDealloc(this_routine,tagNewDiagRhos)
            DEALLOCATE(NewDiagRhos)

!Store the eigenvalues and first element of the eigenvectors of the original star matrix for later use.
            do i=1,iExcit+1

                StarEigens(i,1)=ValsDODMS(1,i)
                StarEigens(i,2)=VecsDODMS(1,i)

            enddo

!Use linear regression technique to find the best-fit gradient for the line for each eigenvalue & vector

!Mean x value is simply the average of the multiplicative constants, and is the same for all eigenvalues and vectors
            meanx=0.0_dp
            do i=1,LinePoints
                meanx=meanx+LineRhoValues(i)
            enddo
            meanx=meanx/LinePoints

!Need to calculate the average for each eigenvalue & vector
            ALLOCATE(MeanVecs(iExcit+1),stat=iErr)
            CALL LogMemAlloc("MeanVecs",iExcit+1,8,this_routine,tagMeanVecs,iErr)
            MeanVecs=0.0_dp
            ALLOCATE(MeanVals(iExcit+1),stat=iErr)
            CALL LogMemAlloc("MeanVals",iExcit+1,8,this_routine,tagMeanVals,iErr)
            MeanVals=0.0_dp

            do i=1,iExcit+1
                do j=1,LinePoints

                    MeanVals(i)=MeanVals(i)+ValsDODMS(j,i)
                    MeanVecs(i)=MeanVecs(i)+VecsDODMS(j,i)

                enddo
                MeanVals(i)=MeanVals(i)/LinePoints
                MeanVecs(i)=MeanVecs(i)/LinePoints
            enddo

!Calculate Sxx - this is the same for all eigenvector/values
            Sxx=0.0_dp
            do i=1,LinePoints
                Sxx=Sxx+(LineRhoValues(i)-meanx)*(LineRhoValues(i)-meanx)
            enddo

!Calculate Sxy and Syy (to calculate R^2) for each eigenvalue/vector.
!Sxy = \sum{(x_i - X)(y_i - Y) where X and Y are the means of x and y respectivly

            ALLOCATE(SxyVals(iExcit+1),stat=iErr)
            CALL LogMemAlloc("SxyVals",iExcit+1,8,this_routine,tagSxyVals,iErr)
            SxyVals=0.0_dp
            ALLOCATE(SxyVecs(iExcit+1),stat=iErr)
            CALL LogMemAlloc("SxyVecs",iExcit+1,8,this_routine,tagSxyVecs,iErr)
            SxyVecs=0.0_dp
            ALLOCATE(SyyVals(iExcit+1),stat=iErr)
            CALL LogMemAlloc("SyyVals",iExcit+1,8,this_routine,tagSyyVals,iErr)
            SyyVals=0.0_dp
            ALLOCATE(SyyVecs(iExcit+1),stat=iErr)
            CALL LogMemAlloc("SyyVecs",iExcit+1,8,this_routine,tagSyyVecs,iErr)
            SyyVecs=0.0_dp

            do i=1,iExcit+1
                do j=1,LinePoints

                    SxyVals(i)=SxyVals(i)+(LineRhoValues(j)-meanx)*(ValsDODMS(j,i)-MeanVals(i))
                    SxyVecs(i)=SxyVecs(i)+(LineRhoValues(j)-meanx)*(VecsDODMS(j,i)-MeanVecs(i))
                    SyyVals(i)=SyyVals(i)+(ValsDODMS(j,i)-MeanVals(i))*(ValsDODMS(j,i)-MeanVals(i))
                    SyyVecs(i)=SyyVecs(i)+(VecsDODMS(j,i)-MeanVecs(i))*(VecsDODMS(j,i)-MeanVecs(i))

                enddo
            enddo

!Now, the best-fit gradients can be calculated for each of the eigenvalues/vector - the gradient for each one is Sxy/Sxx.
!y-intercepts also calculated for use when calculating R^2 value.

            ALLOCATE(GradVals(iExcit+1),stat=iErr)
            CALL LogMemAlloc("GradVals",iExcit+1,8,this_routine,tagGradVals,iErr)
            GradVals=0.0_dp
            ALLOCATE(GradVecs(iExcit+1),stat=iErr)
            CALL LogMemAlloc("GradVecs",iExcit+1,8,this_routine,tagGradVecs,iErr)
            GradVecs=0.0_dp
            ALLOCATE(IncptVals(iExcit+1),stat=iErr)
            CALL LogMemAlloc("IncptVals",iExcit+1,8,this_routine,tagIncptVals,iErr)
            IncptVals=0.0_dp
            ALLOCATE(IncptVecs(iExcit+1),stat=iErr)
            CALL LogMemAlloc("IncptVecs",iExcit+1,8,this_routine,tagIncptVecs,iErr)
            IncptVecs=0.0_dp

            do i=1,iExcit+1
                
                GradVals(i)=SxyVals(i)/Sxx
                GradVecs(i)=SxyVecs(i)/Sxx
                IncptVals(i)=MeanVals(i)-GradVals(i)*meanx
                IncptVecs(i)=MeanVecs(i)-GradVecs(i)*meanx

            enddo

!To calculate the R^2 value for the linear regression, also need to calculate \sum{(y_i - Y)^2} - the expected y value from the gradient calculation at each point

            ALLOCATE(ExpctVals(iExcit+1),stat=iErr)
            CALL LogMemAlloc("ExpctVals",iExcit+1,8,this_routine,tagExpctVals,iErr)
            ExpctVals=0.0_dp
            ALLOCATE(ExpctVecs(iExcit+1),stat=iErr)
            CALL LogMemAlloc("ExpctVecs",iExcit+1,8,this_routine,tagExpctVecs,iErr)
            ExpctVecs=0.0_dp

            do i=1,iExcit+1
                do j=1,LinePoints
                    
                    ExpctVals(i)=ExpctVals(i)+(ValsDODMS(j,i)-(IncptVals(i)+GradVals(i)*LineRhoValues(j)))**2
                    ExpctVecs(i)=ExpctVals(i)+(ValsDODMS(j,i)-(IncptVals(i)+GradVals(i)*LineRhoValues(j)))**2

                enddo
            enddo

            !Now, calculate R^2 value at for each eigenvalue/vector. This has been quite a time-consuming operation as need to calculate Syy..,Incpt..,Expct.. however,does not increase scaling.

            ALLOCATE(RsqVals(iExcit+1),stat=iErr)
            CALL LogMemAlloc("RsqVals",iExcit+1,8,this_routine,tagRsqVals,iErr)
            RsqVals=0.0_dp
            ALLOCATE(RsqVecs(iExcit+1),stat=iErr)
            CALL LogMemAlloc("RsqVecs",iExcit+1,8,this_routine,tagRsqVecs,iErr)
            RsqVecs=0.0_dp

            do i=1,iExcit+1

                RsqVals(i)=1.0_dp-ExpctVals(i)/SyyVals(i)
                RsqVecs(i)=1.0_dp-ExpctVecs(i)/SyyVecs(i)

                IF((RsqVals(i).lt.0.95).or.(RsqVecs(i).lt.0.95)) THEN
                    WRITE(6,*) "Problem with linear approximation, R^2 value: ", RsqVals(i)," or, " &
                        ,RsqVecs(i)," for eigenvalue/vector : ", i
                ELSEIF((RsqVals(i).gt.1.0_dp).or.(RsqVecs(i).gt.1.0_dp)) THEN
                    WRITE(6,*) "Fatal problem in linear approximation, R^2 > 1 : ", RsqVals(i), &
                        " for eigenvalue/vector : ", i
                    STOP
                ENDIF

            enddo

!Deallocate all memory apart from gradients, which are still needed.
            CALL LogMemDealloc(this_routine,tagRsqVals)
            DEALLOCATE(RsqVals)
            CALL LogMemDealloc(this_routine,tagRsqVecs)
            DEALLOCATE(RsqVecs)
            CALL LogMemDealloc(this_routine,tagExpctVals)
            DEALLOCATE(ExpctVals)
            CALL LogMemDealloc(this_routine,tagExpctVecs)
            DEALLOCATE(ExpctVecs)
            CALL LogMemDealloc(this_routine,tagIncptVals)
            DEALLOCATE(IncptVals)
            CALL LogMemDealloc(this_routine,tagIncptVecs)
            DEALLOCATE(IncptVecs)
            CALL LogMemDealloc(this_routine,tagSxyVals)
            DEALLOCATE(SxyVals)
            CALL LogMemDealloc(this_routine,tagSxyVecs)
            DEALLOCATE(SxyVecs)
            CALL LogMemDealloc(this_routine,tagSyyVals)
            DEALLOCATE(SyyVals)
            CALL LogMemDealloc(this_routine,tagSyyVecs)
            DEALLOCATE(SyyVecs)
            CALL LogMemDealloc(this_routine,tagMeanVecs)
            DEALLOCATE(MeanVecs)
            CALL LogMemDealloc(this_routine,tagMeanVals)
            DEALLOCATE(MeanVals)
            CALL LogMemDealloc(this_routine,tagValsDODMS)
            DEALLOCATE(ValsDODMS)
            CALL LogMemDealloc(this_routine,tagVecsDODMS)
            DEALLOCATE(VecsDODMS)
            DEALLOCATE(LineRhoValues)

!Count all diagonalised eigenvectors from excited stars which are attached back to the root.
!The offdiagonal elements from the HF root to the spokes of each excited determinant, a,  with root at j, are given by rho_ij * <j|a>
!Since <j|a> is linearly related to <i|j> with the gradient calculated above, no more diagonalisations need to be performed.
!Cycle through all excitations of the root, and then all diagonalised eigenvectors of a star from each excitation.
!Taking the original star matrix eigenvectors, with rho_jj=1, and using the calculated gradient, the off diagonal element to any of these eigenvectors from the HF root can be calculated.

            CSE=0
!Cycle through excitations of the HF det
            do i=1,iExcit
            
!Cycle through all possible eigenvectors of each excited star
                do j=1,iExcit+1

                    rhoia=(ExcitInfo(i,1))*(GradVecs(j)*((ExcitInfo(i,0))-1.0_dp)+StarEigens(j,2))
                    
                    IF((ABS(rhoia).gt.RhoEps).and.(.not.TQuadValMax).and.(.not.TQuadVecMax)) THEN
                        !Include all quadruple excitations (ficticious and real)
                        CSE=CSE+1
                        
                    ELSEIF((ABS(rhoia).gt.RhoEps).and.(TQuadValMax.or.TQuadVecMax)) THEN
                        !Only want one eigenvector per double excitation to be included.
                        CSE=CSE+1
                        EXIT
                    ENDIF

                enddo
            enddo

            IF((CSE.gt.iExcit).and.(TQuadValMax.or.TQuadVecMax)) THEN
                WRITE(6,*) "Problem when considering only one quadruple eigenvector per double"
                STOP
            ENDIF

            WRITE(6,"(A,I9,A)") "There are ", CSE," extra excitations, from the inclusion of stars of all double excitations"

!Allocate memory to hold all excitations - but do not need to include original star

!            TotExcits=iExcit+CSE
            TotExcits=CSE
            ALLOCATE(ExcitInfo2(0:TotExcits,0:2),stat=iErr)
            CALL LogMemAlloc("ExcitInfo2",(TotExcits+1)*3*HElement_t_size,8,this_routine,tagExcitInfo2,iErr)
            call LogMemAlloc('ExcitInfo2',3*(TotExcits+1),8*HElement_t_size,this_routine,tagExcitInfo2,iErr)
            ExcitInfo2=(0.0_dp)

!Fill original star matrix - NO!! Do not want to double count i --> j excitations.
!Only need to include the original root, i
            
!            do i=0,iExcit
!
!                ExcitInfo2(i,0)=ExcitInfo(i,0)
!                ExcitInfo2(i,1)=ExcitInfo(i,1)
!                ExcitInfo2(i,2)=ExcitInfo(i,2)
!
!            enddo

             ExcitInfo2(0,0)=(1.0_dp)
             ExcitInfo2(0,1)=(1.0_dp)
             ExcitInfo2(0,2)=ExcitInfo(0,2)

!            WRITE(6,*) "Original Hij elements :"
!            WRITE(6,*) ExcitInfo(:,2)

!Add contributions from excited stars - offdiagonal elements as above - diagonal elements are linearly scaled eigenvalues, and the Hamiltonian elements are similarly scaled.

!            NextVertex=iExcit+1
            NextVertex=1

            hmax = 0
            do i=1,iExcit
                ValMax=0.0_dp
                VecMax=0.0_dp
                TVal=.false.
                
                do j=1,iExcit+1

                    rhoia=(ExcitInfo(i,1))*(GradVecs(j)*((ExcitInfo(i,0))-1.0_dp)+StarEigens(j,2))
                    rhoaa=GradVals(j)*(ExcitInfo(i,0)-1.0_dp)+StarEigens(j,1)
                    
!Include all quadruple excitations with large enough connection to root
                    IF((ABS(rhoia).gt.RhoEps).and.(.not.TQuadValMax).and.(.not.TQuadVecMax)) THEN
                        ExcitInfo2(NextVertex,1)=(rhoia)
                        ExcitInfo2(NextVertex,0)=(GradVals(j))*(ExcitInfo(i,0)-(1.0_dp))+(StarEigens(j,1))
                        ExcitInfo2(NextVertex,2)=ExcitInfo(i,2)*((GradVecs(j))*(ExcitInfo(i,0)-(1.0_dp))+(StarEigens(j,2)))
                        NextVertex=NextVertex+1
                    
                    ELSEIF((ABS(rhoia).gt.RhoEps).and.(TQuadValMax.or.TQuadVecMax)) THEN
!Only include from the quadruple excitations, the connection resulting from the largest eigenvector(TQuadVecMax),or the largest eigenvalue(TQuadValMax) - therefore only one quad contribution per double, and O[N^2M^2] scaling.
                        
                        IF((rhoaa.gt.ValMax).and.(TQuadValMax)) THEN
                            ValMax=rhoaa
                            VecMax=rhoia
                            HMax=(ExcitInfo(i,2))*(GradVecs(j)*((ExcitInfo(i,0))-1.0_dp)+StarEigens(j,2))
                            TVal=.true.
                        ELSEIF((ABS(rhoia).gt.VecMax).and.(TQuadVecMax)) THEN
                            ValMax=rhoaa
                            VecMax=rhoia
                            HMax=(ExcitInfo(i,2))*(GradVecs(j)*((ExcitInfo(i,0))-1.0_dp)+StarEigens(j,2))
                            TVal=.true.
                        ENDIF
                    ENDIF

                enddo

                IF(TVal.and.(ValMax.lt.1.0e-9_dp)) STOP 'Error in collecting maximum values'
                
!Put largest values into ExcitInfo2
                IF((TQuadValMax.or.TQuadVecMax).and.(ValMax.gt.1.0e-9_dp).and.TVal) THEN
                    ExcitInfo2(NextVertex,0)=(ValMax)
                    ExcitInfo2(NextVertex,1)=(VecMax)
                    ExcitInfo2(NextVertex,2)=(HMax)
                    NextVertex=NextVertex+1
                ENDIF
                
            enddo

            IF(NextVertex.ne.(TotExcits+1)) STOP 'Incorrect Counting in GetLinStarStars'

!            IF((TQuadValMax.or.TQuadVecMax).and.(NextVertex.gt.2*iExcit+1)) STOP 'Incorrect Counting in GetLinStarStars2'
            
!Return with the new information.
            iExcit=TotExcits
            iMaxExcit=TotExcits

            call LogMemDealloc(this_routine,tagExcitInfo)
            Deallocate(ExcitInfo)

            ! Resort again - so that root is in element 0, then ordered by 
            ! rho_jj in ascending order - probably not needed...
            call sort (excitInfo2(:,0))
            
            ExcitInfo => ExcitInfo2

!            WRITE(6,*) "Off-diag: "
!            WRITE(6,*) ExcitInfo2(:,1)
!            WRITE(6,*) "*****"
!            WRITE(6,*) "Diag: "
!            WRITE(6,*) ExcitInfo2(:,0)
!            WRITE(6,*) "*****"
!            WRITE(6,*) "Hia: "
!            WRITE(6,*) ExcitInfo2(:,2)
!            WRITE(6,*) ""
           
            call halt_timer(proc_timer)

            RETURN
        END SUBROUTINE GetLinStarStars
                
!This subroutine simply forms a star matrix, diagonalises it, and returns the eigenvalues and first elements of the eigenvectors.
        SUBROUTINE GetValsnVecs(Dimen,DiagRhos,OffDiagRhos,Vals,Vecs)
            use global_utilities
            use constants, only: int32
            IMPLICIT NONE
            INTEGER :: Dimen,i,iErr
            INTEGER(int32) INFO
            real(dp) :: DiagRhos(1:Dimen),Vals(Dimen),Vecs(Dimen)
            HElement_t :: OffDiagRhos(2:Dimen)
            real(dp), ALLOCATABLE :: StarMat(:,:),WORK(:)
            integer(TagIntType), save :: tagStarMat=0,tagWORK=0
            character(*), parameter :: this_routine='GetValsnVecs'

!Construct matrix - remember that the first element of OffDiagRhos now is the first offdiagonal element for the excitation, not root.
!For ease, this diagonalisation is to be done initially with DSYEV, though it is possible to reduce the scaling to N^2 rather than N^3

                ALLOCATE(StarMat(Dimen,Dimen),stat=iErr)
                CALL LogMemAlloc("StarMat",Dimen*Dimen,8,this_routine,tagStarMat,iErr)
                StarMat=0.0_dp
                Vals=0.0_dp
                Vecs=0.0_dp

                do i=2,Dimen
                    StarMat(i,i)=DiagRhos(i)
                    StarMat(1,i)=OffDiagRhos(i)
                    StarMat(i,1)=OffDiagRhos(i)
                enddo
                StarMat(1,1)=DiagRhos(1)

                ALLOCATE(WORK(3*Dimen),stat=iErr)
                CALL LogMemAlloc("WORK",Dimen*3,8,this_routine,tagWORK,iErr)
                WORK=0.0_dp

!                do i=1,
!                    do j=1,Dimen
!                        WRITE(6,"F20.14,$") StarMat(i,j)
!                    enddo
!                    write(6,*) ""
!                    write(6,*) ""
!                enddo
!
!                WRITE(6,*) "**************"
                
                CALL DSYEV('V','U',Dimen,StarMat,Dimen,Vals,WORK,3*Dimen,INFO)
                IF(INFO.ne.0) THEN
                    WRITE(6,*) "DYSEV error in GetValsnVecs: ",INFO
                    STOP
                ENDIF

!Store first elements of eigenvectors in 'Vecs'
!The Absolute value of the first element of the eigenvectors is taken - this *should* be ok, since the eigenvectors can always be multiplied by -1 without changing the validity of them.
                do i=1,Dimen
                    Vecs(i)=ABS(StarMat(1,i))
                enddo

                CALL LogMemDealloc(this_routine,tagStarMat)
                Deallocate(StarMat)
                CALL LogMemDealloc(this_routine,tagWORK)
                Deallocate(Work)

            RETURN
        END SUBROUTINE GetValsnVecs


         SUBROUTINE GetStarProds(iExcit,ProdNum,UniqProd,rhii,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,ECore)
            use global_utilities
            use IntegralsData , only : TCalcRealProd,TCalcRhoProd,TSumProd
            use SystemData, only: BasisFN
            use HElem
            IMPLICIT NONE
            INTEGER :: iExcit,ProdNum,Uniqprod,ProdOrbs(8),i_P,nEl,nBasis,nMsh,nMax,nTay(2),ierr,i,ni(nEl)
            INTEGER :: nj(nEl),nk(nEl),nl(nEl)
            complex(dp) fck(*)
            HElement_t UMat(*),rh,rhii
            TYPE(BasisFN) G1(*)
            real(dp) :: Beta,ALat(3),ECore
            character(*), parameter :: this_routine='GetStarProds'

!TCalcRealProd explicitly calculates 'real' quadruple product excitations i.e. neither of the two 'from' or 'to' orbitals are the same in the two excitations from which the product excitation is composed.
!Once found, these can be used in different ways - either they can be explicitly added to the original star graph, or they try to be resummed in to the original graph (TSumProd). 
             IF(TCalcRealProd.and.TSumProd) THEN
!TSumProd tries to find a way to resum the product excitations into the double excitations already available in the star graph.

                rhiiadd=0.0_dp
                ALLOCATE(Offrho(iExcit),stat=iErr)
                CALL LogMemAlloc("OffRho",iExcit,8*HElement_t_size,this_routine,tagOffrho,iErr)

!   Temporarily fill Offrho with the original off-diagonal rho elements, which can be used in countprodexcits, since the rhoelements in EXCITINFO are modified in the routine.
                DO I=1,iExcit
                    Offrho(I)=ExcitInfo(I,1)
                ENDDO

                !This routine finds all real quadruple excitations which are products of double excitations, add resums values from these.
                CALL CountProdExcits(ProdNum,.FALSE.,iExcit,UniqProd)
                WRITE(6,*) ProdNum, "product excitations found, and summed in..."
                Call LogMemDealloc(this_routine,tagOffrho)
                DEALLOCATE(Offrho)
                
                IF(ProdNum.gt.0) THEN
                    WRITE(6,*) "Sum of product diagonal rho elements is ", rhiiadd
                    WRITE(6,*) "Number of unique products is ", UniqProd
!                    rhiiadd=(rhiiadd+(UniqProd))/(2*UniqProd)
!   rhiiadd now in divided by the number of product excitations we are resumming in, to give an average value for the diagonal rho elements for the quadruple product excitations.
                    rhiiadd=rhiiadd/ProdNum
!                    ExcitInfo(0,0)=rhiiadd
!                    ExcitInfo(0,1)=rhiiadd
!   The resummed values for the diagonal elements of the product excitations are added to the root
                    ExcitInfo(0,0)=ExcitInfo(0,0)+rhiiadd
                    ExcitInfo(0,1)=ExcitInfo(0,1)+rhiiadd
                    WRITE(6,*) "New root is now ",ExcitInfo(0,0) 
                ELSE
                    IF(abs(rhiiadd).gt.(0.0_dp)) STOP 'rhiiadd should be zero as no products'
                ENDIF
            ENDIF

!This calculates all real product excitations, and attaches them explicitly to the double excitations from which they are derived
            IF(TCalcRealProd.and.(.NOT.TSumProd)) THEN

!Call twice - once to calculate number of products, so can allocate memory, then store products
                CALL CountProdExcits(ProdNum,.FALSE.,iExcit,UniqProd)
                WRITE(6,*) ProdNum, "product excitations found - allocating memory..."
                ALLOCATE(ProdPositions(2,ProdNum),stat=iErr)
                CALL LogMemAlloc("EXCITSTORE",ProdNum,4,this_routine,tagProdPositions,iErr)
!Prodpositions stores the position of the two excitations in EXCITSTORE which give rise to a real product excitation
                CALL CountProdExcits(ProdNum,.TRUE.,iExcit,UniqProd)
             
!Allocate memory for on diagonal rho elements of product excitations, and 2 x offdiagonal elements
                ALLOCATE(OnDiagProdRho(ProdNum),stat=ierr)
                CALL LogMemAlloc("OnDiagProdRho",ProdNum,8,this_routine,tagOnDiagProdRho,ierr)
                ALLOCATE(OffDiagProdRho(2,ProdNum),stat=ierr)
                CALL LogMemAlloc("OffDiagProdRho",2*ProdNum,8,this_routine,tagOffDiagProdRho,ierr)

                DO I=1,ProdNum
!Approximate on-diag elements as product of consituent excits (exact for FPLD)(remember, they are divided by rhoii^2)
                    OnDiagProdRho(I)=(ExcitInfo(ProdPositions(1,I),0)*ExcitInfo(ProdPositions(2,I),0))

                    IF(TCalcRhoProd) THEN
                        IF(I.eq.1) WRITE(6,*) "Calculating off-diagonal rho elements for product excitations exactly"
!This calculates the exact rho elements between the product excitation, and the constituent determinants
                        prodorbs(1:2)=EXCITSTORE(1:2,ProdPositions(1,I))
                        prodorbs(3:4)=EXCITSTORE(1:2,ProdPositions(2,I))
                        prodorbs(5:6)=EXCITSTORE(3:4,ProdPositions(1,I))
                        prodorbs(7:8)=EXCITSTORE(3:4,ProdPositions(2,I))
                        CALL GetFullPath(nI,nEl,4,prodorbs,nJ)
                        CALL GetFullPath(nI,nEl,2,EXCITSTORE(:,ProdPositions(1,I)),nK)
                        CALL GetFullPath(nI,nEl,2,EXCITSTORE(:,ProdPositions(2,I)),nL)
!nJ is now the product IPATH, and nK and nL the two consituent determinants
                        CALL CalcRho2(nK,nJ,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,2,ECore)
                        OffDiagProdRho(1,I)=(rh/rhii)
                        CALL CalcRho2(nL,nJ,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,2,ECore)
                        OffDiagProdRho(2,I)=(rh/rhii)
                    ELSE
                        IF(I.eq.1) WRITE(6,*) "Approximating off-diagonal elements for product excitations"
!Approximate the off-diag rho elements between product excitation and constituent determinant as being equal to the rho element between the root and other constituent determinant.
                        OffDiagProdRho(1,I)=ExcitInfo(ProdPositions(2,I),1)
                        OffDiagProdRho(2,I)=ExcitInfo(ProdPositions(1,I),1)
                    ENDIF
                ENDDO
            ENDIF

        RETURN

        END SUBROUTINE GetStarProds


!      ! This routine finds all "unique" product excitations of the star
!      SUBROUTINE COUNTUNIQPRODEXCITS(iMaxExcit,ProdNum,iExcit,nouniqprod,countprods)
!        IMPLICIT NONE
!        INTEGER iMaxExcit,ProdNum
!        INTEGER iExcit,I,J
!        INTEGER countprods,nouniqprod,tempprodfrom(4),tempprodto(4),ierr,k
!        LOGICAL UNIQPROD
!        INTEGER, ALLOCATABLE :: UNIQPRODS(:,:)
!!        INTEGER, POINTER :: tempprods(:,:)
!
!        !Space allocated to store the excitations of the unique products, so we can search through them to ensure they are not being duplicated
!        ALLOCATE(UNIQPRODS(8,(iExcit*(iExcit-1)/2)),stat=ierr)
!        CALL LogMemAlloc("UNIQPRODSALL",iExcit*(iExcit-1)*2,8,this_routine,tagUNIQPRODS,ierr)
!        
!        countprods=0
!        nouniqprod=0
!        DO I=1,iExcit
!            DO J=(I+1),iExcit
!                UNIQPROD=.TRUE.
!                ! Again, the criteria for a product excitation are that none of the "ij" or "ab" orbitals match in the constituent excitations
!                IF((EXCITSTORE(1,I).NE.EXCITSTORE(1,J)).AND.(EXCITSTORE(1,I).NE.EXCITSTORE(2,J))                   &
!                .AND.(EXCITSTORE(2,I).NE.EXCITSTORE(2,J)).AND.(EXCITSTORE(2,I).NE.EXCITSTORE(1,J))) THEN
!                    IF((EXCITSTORE(3,I).NE.EXCITSTORE(3,J)).AND.(EXCITSTORE(3,I).NE.EXCITSTORE(4,J))               &
!                    .AND.(EXCITSTORE(4,I).NE.EXCITSTORE(3,J)).AND.(EXCITSTORE(4,I).NE.EXCITSTORE(4,J))) THEN
!                
!                        countprods=countprods+1
!                        tempprodfrom(1:2)=EXCITSTORE(1:2,I)
!                        tempprodfrom(3:4)=EXCITSTORE(1:2,J)
!                        tempprodto(1:2)=EXCITSTORE(3:4,I)
!                        tempprodto(3:4)=EXCITSTORE(3:4,J)
!                        CALL NECI_SORTI(4,tempprodfrom(:))
!                        CALL NECI_SORTI(4,tempprodto(:))
!                        !Cycle through other stored products to see if unique
!                        DO K=1,nouniqprod
!                            !Test if unique
!                            IF((tempprodfrom(1).eq.UNIQPRODS(1,K)).AND.(tempprodfrom(2).eq.UNIQPRODS(2,K)).AND.     &
!                                (tempprodfrom(3).eq.UNIQPRODS(3,K)).AND.(tempprodfrom(4).eq.UNIQPRODS(4,K))) THEN
!                                IF((tempprodto(1).eq.UNIQPRODS(5,K)).AND.(tempprodto(2).eq.UNIQPRODS(6,K)).AND.     &
!                                    (tempprodto(3).eq.UNIQPRODS(7,K)).AND.(tempprodto(4).eq.UNIQPRODS(8,K))) THEN
!                                    UNIQPROD=.FALSE.
!                                    !Sum in off-diagonal contribution from another way to get same product
!                                    QUADRHOS(2,K)=QUADRHOS(2,K)+(ExcitInfo(I,1)*ExcitInfo(J,1))/(ExcitInfo(I,1)+ExcitInfo(J,1))
!                                    !If the excitation is found, then exit loop
!                                    EXIT
!                                ENDIF
!                            ENDIF
!                        ENDDO
!                        IF(UNIQPROD) THEN
!!                            WRITE(67,*) "Unique Product Excitation: ", I,J
!!                            WRITE(67,*) "Diagonal element is: ", ExcitInfo(I,0)*ExcitInfo(J,0)
!!                           If the product is unique, then count it, and store its information in QUADRHOS
!                            nouniqprod=nouniqprod+1
!                            UNIQPRODS(1:4,nouniqprod)=tempprodfrom(:)
!                            UNIQPRODS(5:8,nouniqprod)=tempprodto(:)
!                            QUADRHOS(1,nouniqprod)=ExcitInfo(I,0)*ExcitInfo(J,0)
!                            QUADRHOS(2,nouniqprod)=(ExcitInfo(I,1)*ExcitInfo(J,1))/(ExcitInfo(I,1)+ExcitInfo(J,1))
!                            !Let hamiltonian elements to these excitations = 0
!                            QUADRHOS(3,nouniqprod)=0.0_dp
!                        ENDIF
!                    ENDIF
!                ENDIF
!            ENDDO
!        ENDDO
!        CALL LogMemDealloc(this_routine,tagUNIQPRODS)
!        DEALLOCATE(UNIQPRODS)
!        
!        ProdNum=countprods
!      
!      END
                        

      SUBROUTINE CountProdExcits(ProdNum,setup,iExcit,UniqProd)
        use IntegralsData , only : TSumProd
        use HElem
        IMPLICIT NONE
        INTEGER ProdNum!,ProdPositions(2,length),EXCITSTORE(4,iMaxExcit),
        INTEGER iExcit,I,J
        INTEGER countprods,UniqProd
        LOGICAL setup,FIRST

        countprods=0
        UniqProd=0
        rhiiadd=0.0_dp
        !Cycle through all unique excitation products
        DO I=1,iExcit
            !Only allow one product excitation per original excitation attached to root
            FIRST=.TRUE.
            DO J=(I+1),iExcit
            !Test to see that none of the orbitals which are moving overlap - if not, then we have a real product excitation.
            IF((EXCITSTORE(1,I).NE.EXCITSTORE(1,J)).AND.(EXCITSTORE(1,I).NE.EXCITSTORE(2,J))               &
                .AND.(EXCITSTORE(2,I).NE.EXCITSTORE(2,J)).AND.(EXCITSTORE(2,I).NE.EXCITSTORE(1,J))) THEN
                IF((EXCITSTORE(3,I).NE.EXCITSTORE(3,J)).AND.(EXCITSTORE(3,I).NE.EXCITSTORE(4,J))               &
                    .AND.(EXCITSTORE(4,I).NE.EXCITSTORE(3,J)).AND.(EXCITSTORE(4,I).NE.EXCITSTORE(4,J))) THEN
                    countprods=countprods+1
     
                    IF(setup) THEN
                        !record the number of the two original excitations which create the product
                        ProdPositions(1,countprods)=I
                        ProdPositions(2,countprods)=J
                    ENDIF
                    IF(TSumProd) THEN
                        !Add the diagonal rho elements to the quadruple excitation back into the offdiagonal rho elements of the original matrix
                        ExcitInfo(I,1)=ExcitInfo(I,1)+Offrho(J)
                        ExcitInfo(J,1)=ExcitInfo(J,1)+Offrho(I)
                        !Only sum in the diagonal product rho element once per excitation
                        IF(FIRST) THEN
                            !Resum by adding in the product of the diagonal rhoelements from the constituent excitations
                            rhiiadd=rhiiadd+(ExcitInfo(J,0)*ExcitInfo(I,0))
                            !Count the products we are adding in
                            UniqProd=UniqProd+1
                            FIRST=.FALSE.
                        ENDIF
                    ENDIF
                ENDIF
            ENDIF
            ENDDO
        ENDDO
        IF(.NOT.setup) ProdNum=countprods
      END subroutine

      END MODULE
        
      
!.. This version creates a star of all 2-vertex terms.
!.. This sets up the excitation generators and the memory - using the old excitation generators
FUNCTION FMCPR3STAR(NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,NTAY, &
                     RHOEPS,LSTE,ICE,L,LT,NWHTAY,TSYM,ECORE,ILMAX,DBETA,DLWDB)
         use constants, only: dp      
         use SystemData, only: BasisFN
         use global_utilities
         use legacy_data, only: irat
         use HElem
         use MemoryManager, only: TagIntType
         use util_mod, only: neci_icopy
         IMPLICIT NONE
         INTEGER NEL,I_P,NBASISMAX(*),NBASIS,NMSH,NMAX
         INTEGER NTAY(2),NWHTAY,LT
         real(dp) ALAT(*),UMAT(*),ECORE
         complex(dp) FCK(*)
         INTEGER NI(NEL),ILMAX
         type(BasisFN) G1(*)
!.. These are old and no longer used
!.. LSTE is a list of excitations (which we will generate)
!.. ICE is the IC of each excitation (i.e. how much it differs from us (INODE)
         INTEGER LSTE(NEL,0:ILMAX)
         INTEGER ICE(0:ILMAX)
         real(dp) FMCPR3Star
!.. New lists are generated here
         INTEGER,ALLOCATABLE :: NLSTE(:,:)
         real(dp),ALLOCATABLE ::  NLIST(:)
         INTEGER,ALLOCATABLE :: NICE(:)
         INTEGER(TagIntType),SAVE :: tagNLSTE=0,tagNLIST=0,tagNICE=0
         INTEGER :: err

!.. This will contain all the info needed to work out the value of the
!.. star
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ

         real(dp) BETA,RHOEPS
         LOGICAL TSYM
         real(dp) DBETA,DLWDB
         INTEGER NLENLIST,L

         INTEGER NORDER,NMIN
         real(dp) FMCPR3STAR2
         character(*), parameter :: this_routine='FMCPR3STAR'
  
         IF(HElement_t_size.NE.1) STOP 'FMCPR3STAR cannot work with complex orbitals.' 
         SELECT CASE (IAND(NWHTAY,24))
         CASE(0)
!.. Allow both singles and doubles
            NORDER=2
            NMIN=1
         CASE(8)
!.. only singles
            NORDER=1
            NMIN=1
         CASE(16)
            NORDER=2
            NMIN=2
         CASE(24)
            STOP "Invalid combination of flags in NWHTAY"
         END SELECT
!.. Count the excitations.
         CALL GENEXCIT(NI,NORDER,NBASIS,NEL,LSTE,ICE,NLENLIST,NMIN,G1,TSYM,NBASISMAX,.TRUE.)
!.. Allocate memory for the lists
         ILMAX=NLENLIST
         allocate(NLIST(0:NLENLIST*3),stat=err)
         call LogMemAlloc('NLIST',(NLENLIST+1)*3,8,this_routine,tagNLIST)
         allocate(NLSTE(NEL,0:NLENLIST),stat=err)
         call LogMemAlloc('NLSTE',(NLENLIST+1)*NEL,4,this_routine,tagNLSTE)
         allocate(NICE(0:NLENLIST),stat=err)
         call LogMemAlloc('NICE',(NLENLIST+1),4,this_routine,tagNICE)

         CALL NECI_ICOPY(NEL,NI,1,NLSTE(1,0),1)
         NICE(0)=0
!.. Now generate the excitations
         CALL GENEXCIT(NI,NORDER,NBASIS,NEL,NLSTE(1,1),NICE(1),NLENLIST,NMIN,G1,TSYM,NBASISMAX,.FALSE.)

         FMCPR3STAR=FMCPR3STAR2(NI,BETA,I_P,NEL,G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,NTAY, &
           RHOEPS,NLSTE,NICE,L,LT,NWHTAY,ECORE,NLENLIST,DBETA,DLWDB)

         deallocate(NLIST,NICE,NLSTE)
         call LogMemDealloc(this_routine,tagNLIST)
         call LogMemDealloc(this_routine,tagNLSTE)
         call LogMemDealloc(this_routine,tagNICE)
         RETURN
      END

      
!.. This version creates a star of all 2-vertex terms.
!..   This is the heart of the function, called once the excitations are found.
      FUNCTION FMCPR3STAR2(NI,BETA,I_P,NEL,G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT,NTAY, &
           RHOEPS,LSTE,ICE,L,LT,NWHTAY,ECORE,ILMAX,DBETA,DLWDB)
         use constants, only: dp     
         use SystemData, only: BasisFN
         use Determinants, only: get_helement
         use util_mod, only: neci_icopy
         IMPLICIT NONE
         real(dp) FMCPR3Star2
         TYPE(BasisFN) G1(*)
         INTEGER NEL,I_P,NBASIS,NMSH,NMAX
         INTEGER NTAY(2),NWHTAY,LT
         real(dp) ALAT(*),ECORE
         HElement_t UMat(*)
         complex(dp) FCK(*)
         INTEGER NI(NEL),ILMAX
!.. LSTE is a list of excitations (which we will generate)
!.. ICE is the IC of each excitation (i.e. how much it differs from us (INODE)
         INTEGER LSTE(NEL,0:ILMAX)
         INTEGER ICE(0:ILMAX)
         INTEGER NROOTS
!.. This will contain all the info needed to work out the value of the
!.. star
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ

         HElement_t LIST(0:ILMAX,0:2)
         real(dp) BETA,RHOEPS
         real(dp) DBETA,DLWDB
         INTEGER NLIST,NLCUR,I,L

         HElement_t RHII,RH
         NLIST=ILMAX

         NLCUR=0
         rhii = 0
         DO I=0,NLIST
            IF(ICE(I).EQ.0.AND.I.GT.0) LSTE(1,I)=0
            IF(LSTE(1,I).NE.0) THEN
               CALL CALCRHO2(NI,LSTE(1,I),BETA,I_P,NEL,G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT, &
                  RH,NTAY,ICE(I),ECORE)
               IF(abs(RH ).gt. RHOEPS) THEN
                  IF(NLCUR.NE.I) CALL NECI_ICOPY(NEL,LSTE(1,I),1,LSTE(1,NLCUR),1)
                  ICE(NLCUR)=ICE(I)
                  IF(NLCUR.EQ.0) RHII=RH
                  LIST(NLCUR,1)=RH/RHII
                  CALL CALCRHO2(LSTE(1,I),LSTE(1,I),BETA,I_P,NEL,G1,NBASIS,NMSH,FCK,NMAX,ALAT,UMAT, &
                    RH,NTAY,0,ECORE)
                  LIST(NLCUR,0)=RH/RHII
                  LIST(nLCur, 2) = get_helement (nI, LSTE(:,I))
!                  CALL WRITEDET(6,LSTE(1,I),NEL,.FALSE.)
!                  WRITE(6,*) (LIST(NLCUR,J),J=0,2)
                  NLCUR=NLCUR+1
               ENDIF
            ENDIF
         ENDDO
!NLCUR is one after the last element, i.e. the total number of elements
         L=NLCUR
         LT=ILMAX+1
!.. we now have a list length NLCUR of dets in the star (from 0:NLCUR-1)
!.. Call a routine to generate the value of the star
         WRITE(6,*) NLCUR," determinants in star"
         IF(.NOT.BTEST(NWHTAY,0)) THEN
            WRITE(6,*) "Beginning Complete Star Diagonalization"
            CALL STARDIAG(NLCUR,LIST,ILMAX+1,I_P,FMCPR3STAR2,DBETA,DLWDB)
         ELSE
            WRITE(6,*) "Beginning Polynomial Star Diagonalization"
            NROOTS=NLCUR
            IF(BTEST(NWHTAY,1)) NROOTS=1
            CALL STARDIAG2(NLCUR,LIST,ILMAX+1,I_P,FMCPR3STAR2,DBETA,DLWDB,NROOTS)
         ENDIF
      END

      SUBROUTINE STARDIAGREALPROD(NLIST,LIST,ILMAX,I_P,SI,DBETA,DLWDB,ProdNum,ProdPositions,OnDiagProdRho,OffDiagProdRho)
         !NLIST is now no. original excitations (+ root) - ILMAX is max possible excitations +1
         use constants, only: dp, int32
         use HElem
         use global_utilities
        use MemoryManager, only: TagIntType
         IMPLICIT NONE
         INTEGER I_P
         INTEGER NLIST,ILMAX,ProdNum,TOTVERT
         real(dp) OnDiagProdRho(ProdNum),OffDiagProdRho(2,ProdNum)
         real(dp) LIST(ILMAX,0:2)
         real(dp),ALLOCATABLE ::  RIJMAT(:),WLIST(:),WORK(:)
         INTEGER(TagIntType), SAVE :: tagRIJMAT=0,tagWLIST=0,tagWORK=0

         type(timer), save :: proc_timer
         INTEGER WORKL,ProdPositions(2,ProdNum)
         INTEGER(int32) INFO
         real(dp) SI,DLWDB,DBETA
         INTEGER I,J,err
         character(*),parameter :: this_routine='STARDIAGREALPROD'
         
         IF(HElement_t_size.GT.1) STOP "STARDIAGREALPROD cannot function with complex orbitals."

         proc_timer%timer_name='STARDIAGRP'
         call set_timer(proc_timer)
         !Is there a need to sort the matrix? If there is, we have problems!
         !CALL SORT3RN(NLIST-1,LIST(2,0),LIST(2,1),LIST(2,2),HElement_t_size)
         
         TOTVERT=NLIST+ProdNum
         allocate(RIJMAT(TOTVERT**2),stat=err)
         call LogMemAlloc('RIJMAT',TOTVERT**2,8,this_routine,tagRIJMAT,err)
         RIJMAT=0.0_dp
        
         !Fill RIJMAT
         DO I=1,TOTVERT
            IF(I.LE.NLIST) THEN
                RIJMAT(I)=LIST(I,1)
                RIJMAT((I-1)*TOTVERT+I)=LIST(I,0)
            ELSE
                RIJMAT((ProdPositions(1,I-NLIST)*TOTVERT)+I)=OffDiagProdRho(1,I-NLIST)
                RIJMAT((ProdPositions(2,I-NLIST)*TOTVERT)+I)=OffDiagProdRho(2,I-NLIST)
                RIJMAT((I-1)*(TOTVERT)+I)=OnDiagProdRho(I-NLIST)
            ENDIF
         ENDDO

!.. Debug info         
!         WRITE(68,*) "ROOT RHOII, RHOIJ ", LIST(1,0),LIST(1,1)
!         DO I=2,NLIST
!            WRITE(68,"A,I3,A,2E14.6,4I4") "EXCITATION ",I-1," - RHOII, RHOIJ, IPATH: ",LIST(I,0),LIST(I,1),EXCITSTORE(:,I-1)
!         ENDDO
!         
!         DO I=1,ProdNum
!            WRITE(68,"A,I3,A,I3,I3") "PRODUCT EXCITATION ",I," FORMED FROM EXCITATIONS ", ProdPositions(1,I),ProdPositions(2,I)
!            WRITE(68,*) ""
!         ENDDO
!         
!         DO I=1,TOTVERT
!             DO J=1,TOTVERT
!                 WRITE(68,"E14.6,$") RIJMAT((I-1)*TOTVERT+J)
!             ENDDO
!             WRITE(68,*) ""
!             WRITE(68,*) ""
!         ENDDO
        
         allocate(WLIST(TOTVERT),stat=err)
         call LogMemAlloc('WLIST',TOTVERT,8,this_routine,tagWLIST,err)
         WORKL=3*TOTVERT
         allocate(WORK(WORKL),stat=err)
         call LogMemAlloc('WORK',WORKL,8,this_routine,tagWORK,err)
        
!.. Diagonalize
         CALL DSYEV('V','L',TOTVERT,RIJMAT,TOTVERT,WLIST,WORK,WORKL,INFO)
         IF(INFO.NE.0) THEN
            WRITE(6,*) 'DYSEV error: ',INFO
            STOP
         ENDIF
         deallocate(WORK)
         call LogMemDealloc(this_routine,tagWORK)
         WRITE(6,*)
         WRITE(6,*) "Highest root:",WLIST(TOTVERT)

         SI=0.0_dp
         DO I=1,TOTVERT
            SI=SI+RIJMAT(((I-1)*TOTVERT)+1)*RIJMAT(((I-1)*TOTVERT)+1)*(WLIST(I)**I_P)
            IF(DBETA.NE.0.0_dp) THEN
                DO J=1,NLIST
                    !only sum over vertices linked to i?
                    DLWDB=DLWDB+RIJMAT(((I-1)*TOTVERT)+1)*RIJMAT(((I-1)*TOTVERT)+J)*(WLIST(I)**I_P)*LIST(J,2)
                ENDDO
            ENDIF
         ENDDO
         WRITE(6,*) "Final SI= ", SI
         SI=SI-1.0_dp
         DLWDB=DLWDB-LIST(1,2)
         deallocate(WLIST,RIJMAT)
         call LogMemDealloc(this_routine,tagWLIST)
         call LogMemDealloc(this_routine,tagRIJMAT)
         call halt_timer(proc_timer)

         RETURN
      END SUBROUTINE

         
         
      SUBROUTINE STARDIAGSC(NLIST,LIST,ILMAX,I_P,SI,DBETA,DLWDB)
         use constants, only: dp,int32
         use global_utilities
         use sort_mod
         use helem, only: helement_t_size
         use MemoryManager, only: TagIntType
         IMPLICIT NONE
         INTEGER I_P
         INTEGER NLIST,ILMAX
         real(dp) LIST(ILMAX,0:2)
         real(dp),ALLOCATABLE ::  RIJMAT(:),WLIST(:),WORK(:)
         INTEGER(TagIntType), SAVE :: tagRIJMAT=0,tagWLIST=0,tagWORK=0
         real(dp), DIMENSION(:,:), POINTER :: AOFFDB
         real(dp), DIMENSION(:), POINTER :: AONDB
         integer(TagIntType), save :: tagAOFFDB,tagAONDB
         INTEGER IND,TOTVERT
         type(timer), save :: proc_timer
         INTEGER WORKL,PRODVERT,ierr
         INTEGER(int32) INFO
         real(dp) SI,DLWDB,DBETA
         INTEGER I,J,err
         character(*),parameter :: this_routine='STARDIAGSC'

         IF(HElement_t_size.GT.1) STOP "STARDIAGSC cannot function with complex orbitals."

         proc_timer%timer_name='STARDIAGSC'
         call set_timer(proc_timer)

         ! n.b. This ONLY works with helementsize = 1, i.e. real integrals.
         call sort (list(2:nList-1,0), list(2:nList-1,1), list(2:nList-1,2))
         
         PRODVERT=(NLIST-1)*(NLIST-2)/2
         TOTVERT=PRODVERT+NLIST
         
         ALLOCATE(AOFFDB(PRODVERT,NLIST-1),STAT=ierr)
         CALL LogMemAlloc('AOFFDB',(NLIST-1)*PRODVERT,8,this_routine,tagAOFFDB,ierr)
         AOFFDB=0.0_dp
         ALLOCATE(AONDB(PRODVERT),STAT=ierr)
         CALL LogMemAlloc('AONDB',PRODVERT,8,this_routine,tagAONDB,ierr)
         AONDB=0.0_dp

         IND=0
         DO I=2,NLIST
            DO J=(I+1),NLIST
                IND=IND+1
                IF(IND.gt.PRODVERT) STOP 'Error - IND larger than PRODVERT'
                AONDB(IND)=LIST(I,0)*LIST(J,0)
                AOFFDB(IND,(I-1))=LIST(J,1)
                AOFFDB(IND,(J-1))=LIST(I,1)
            ENDDO
        ENDDO
        WRITE(6,*) "VERTICES ADDED = ", IND
        IF(PRODVERT.NE.IND) THEN
            WRITE(6,*) "EXPECTED EXTRA VERTICES = ", PRODVERT
            WRITE(6,*) "VERTICES ADDED = ", IND
            CALL neci_flush(6)
            STOP 'WRONG NUMBER OF ADDED VERTICES'
        ENDIF

!        DO I=1,(NLIST-1)
!            DO J=1,PRODVERT
!                WRITE(68,"E15.6,$") AOFFDB(J,I)
!            ENDDO
!            WRITE(68,*) ""
!        ENDDO
!        WRITE(68,*) "*****************"
        
        allocate(RIJMAT(TOTVERT**2),stat=err)
        call LogMemAlloc('RIJMAT',TOTVERT**2,8,this_routine,tagRIJMAT,err)
        RIJMAT=0.0_dp
        
        DO I=1,TOTVERT
            IF(I.LE.NLIST) THEN
                RIJMAT(I)=LIST(I,1)
                RIJMAT((I-1)*TOTVERT+I)=LIST(I,0)
                IF(I.GT.1) THEN
                    DO J=1,PRODVERT
                        RIJMAT((I-1)*TOTVERT+NLIST+J)=AOFFDB(J,I-1)
                    ENDDO
                ENDIF
            ELSE
                RIJMAT((I-1)*(TOTVERT)+I)=AONDB(I-NLIST)
            ENDIF
        ENDDO

        CALL LogMemDealloc(this_routine,tagAONDB)
        Deallocate(AONDB)
        NULLIFY(AONDB)
        CALL LogMemDealloc(this_routine,tagAOFFDB)
        Deallocate(AOFFDB)
        NULLIFY(AOFFDB)

!        DO I=1,TOTVERT
!            DO J=1,TOTVERT
!                WRITE(68,"E14.6,$") RIJMAT((I-1)*TOTVERT+J)
!            ENDDO
!            WRITE(68,*) ""
!        ENDDO

        allocate(WLIST(TOTVERT),stat=err)
        call LogMemAlloc('WLIST',TOTVERT,8,this_routine,tagWLIST,err)
        WORKL=3*TOTVERT
        allocate(WORK(WORKL),stat=err)
        call LogMemAlloc('WORK',WORKL,8,this_routine,tagWORK,err)
        
!.. Diagonalize
         CALL DSYEV('V','L',TOTVERT,RIJMAT,TOTVERT,WLIST,WORK,WORKL,INFO)
         IF(INFO.NE.0) THEN
            WRITE(6,*) 'DYSEV error: ',INFO
            STOP
         ENDIF
         deallocate(WORK)
         call LogMemDealloc(this_routine,tagWORK)
         WRITE(6,*)
         WRITE(6,*) "Highest root:",WLIST(TOTVERT)

         SI=0.0_dp
         DO I=1,TOTVERT
            SI=SI+RIJMAT(((I-1)*TOTVERT)+1)*RIJMAT(((I-1)*TOTVERT)+1)*(WLIST(I)**I_P)
            IF(DBETA.NE.0.0_dp) THEN
!                OD=DLWDB
                DO J=1,NLIST
                    !Is this a correct formulation for the hamiltonian elements - only sum over vertices linked to i
                    DLWDB=DLWDB+RIJMAT(((I-1)*TOTVERT)+1)*RIJMAT(((I-1)*TOTVERT)+J)*(WLIST(I)**I_P)*LIST(J,2)
                ENDDO
            ENDIF
         ENDDO
         WRITE(6,*) "Final SI= ", SI
         SI=SI-1.0_dp
         DLWDB=DLWDB-LIST(1,2)
         deallocate(WLIST,RIJMAT)
         call LogMemDealloc(this_routine,tagRIJMAT)
         call LogMemDealloc(this_routine,tagWLIST)
         call halt_timer(proc_timer)
         
         RETURN
      END

!Use a MC sampling technique to find the eigenvector correesponding to the smallest eigenvalue.
!List now contains diagonal hamiltonian matrix elements-Hii in the ExcitInfo(i,0), rather than rho elements.
      SUBROUTINE StarDiagMC(NList,List,ILMax,SI,DLWDB,MaxDiag)
         use constants, only: dp,int32,sizeof_int
         use CalcData , only : InitWalkers,NMCyc,G_VMC_Seed,DiagSft,Tau,SftDamp,StepsSft
         use CalcData , only : TReadPops,ScaleWalkers,TBinCancel, iPopsFileNoRead, iPopsFileNoWrite
         USE Logging , only : TPopsFile,TCalcWavevector,WavevectorPrint, tIncrementPops
         USE global_utilities
         use sort_mod
         use StarDiagData, only: star_walker
         use helem, only: helement_t_size
         use MemoryManager, only: TagIntType
         use util_mod, only: get_unique_filename
         IMPLICIT NONE
         CHARACTER(len=*), PARAMETER :: this_routine='StarDiagMC'
         INTEGER :: i,j,NList,ILMax,ierr,WorkL,toprint,PreviousNMCyc
         INTEGER(int32) Info
         type(timer), save :: proc_timer
         INTEGER :: TotWalkers,Seed,VecSlot,TotWalkersNew,DetCurr,ReadWalkers
         INTEGER :: MaxWalkers,TotWalkersOld,NWalk,k,l,TotWalkersDet
         INTEGER(TagIntType) :: HMatTag=0,WListTag=0,WorkTag=0,WalkVecTag=0,WalkVec2Tag=0
         INTEGER :: SumNoatHF
         real(dp) :: List(ILMax,0:2),SI,DLWDB,MaxDiag,Ran2,Norm,GrowRate
         real(dp) :: rat,SumENum,ProjectionE
         LOGICAL :: DetSign
         real(dp) , ALLOCATABLE :: HMat(:,:),WList(:),Work(:)
         TYPE(star_walker) , POINTER :: WalkVec(:),WalkVec2(:)
         character(255) :: popsfile
         
         proc_timer%timer_name='StarDiagMC'
         call set_timer(proc_timer)

         IF(HElement_t_size.GT.1) THEN
             CALL Stop_All("StarDiagMC","StarDiagMC cannot function with complex orbitals.")
         ENDIF
         
         IF(Tau.eq.0.0_dp) THEN
            WRITE(6,*) "Choosing Tau so that Tau*Hii_max = 0.999 "
            Tau=0.999/MaxDiag
         ELSE
             IF(Tau.gt.0.999/MaxDiag) THEN
                 WRITE(6,*) "Illegal value of Tau chosen - resetting..."
                 WRITE(6,*) "Choosing Tau so that Tau*Hii_max = 0.999 "
                 Tau=0.999/MaxDiag
             ENDIF
         ENDIF
         WRITE(6,*) "Tau = ", Tau

         IF(TReadPops) THEN
             WRITE(6,*) "Reading in POPSFILE to restart calculation..."
             call get_unique_filename('POPSFILE',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
             OPEN(17,FILE=popsfile,Status='old')
             READ(17,*) InitWalkers
             READ(17,*) DiagSft
             READ(17,*) PreviousNMCyc
             WRITE(6,*) "Initial number of walkers read to be: ", InitWalkers
         ELSE
             WRITE(6,*) "Initial number of walkers chosen to be: ", InitWalkers
         ENDIF

         WRITE(6,*) "Damping parameter for Diag Shift set to: ", SftDamp

         toprint=MIN(NList,8)
!Initialise random number seed
         Seed=G_VMC_Seed
         SumNoatHF=0
         SumENum=0.0_dp
         OPEN(15,file='StarMCStats',status='unknown')
         IF(TCalcWavevector) THEN
             WRITE(6,*) "Writing first ",toprint," components of vector to file..."
             OPEN(14,file='VecConv',status='unknown')
         ENDIF

         WRITE(6,*) "Diagonal root value is: ", List(1,0)
         WRITE(6,*) "Largest Diagonal value is: ", MaxDiag
         IF(DiagSft.gt.0.0_dp) THEN
             CALL Stop_All("StarDiagMC","Intial value of DiagSft should be negative.")
         ELSE
             WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is: ", DiagSft
         ENDIF

         IF(TCalcWavevector) THEN
!If TCalcWavevector is on, then we perform a full diagonalisation of the system to compare, and output results
!In addition, we calculate the full wavefunction at each point of the MC, again for comparison

             ALLOCATE(HMat(NList,NList),stat=ierr)
             CALL LogMemAlloc('HMat',NList**2,8,this_routine,HMatTag,ierr)
             HMat=0.0_dp
             ALLOCATE(WList(NList),stat=ierr)
             CALL LogMemAlloc('WList',NList,8,this_routine,WListTag,ierr)
             WList=0.0_dp
             WorkL=3*NList
             ALLOCATE(Work(WorkL),stat=ierr)
             CALL LogMemAlloc('Work',WorkL,8,this_routine,WorkTag,ierr)
             

!.. Now we fill the HMat array
             DO i=1,NLIST
                 HMat(i,1)=List(i,2)
                 HMat(1,i)=List(i,2)
                 HMat(i,i)=List(i,0)
             ENDDO

!             WRITE(67,*) "Size of HMat is ", NLIST
!             WRITE(67,*) "********"
!             WRITE(67,*) "OFF DIAG LIST is "
!             DO I=1,ILMAX
!                 WRITE(67,"E14.6,$") LIST(I,2)
!             ENDDO
!             WRITE(67,*) ""
!             WRITE(67,*) "********"
!             WRITE(67,*) "ON DIAG LIST is "
!             DO I=1,ILMAX
!                 WRITE(67,"E14.6,$") LIST(I,0)
!             ENDDO
!             WRITE(67,*) ""
!             WRITE(67,*) "********"
!             WRITE(67,*) "HMAT is "
!             DO I=1,NLIST
!                 DO J=1,NLIST
!                     WRITE(67,"E14.6,$") HMat(i,j)
!                 ENDDO
!                 WRITE(67,*) ""
!             ENDDO
!             WRITE(67,*) "************"

!.. Diagonalize
             CALL DSYEV('V','L',NList,HMat,NList,WList,Work,WorkL,Info)
             IF(Info.NE.0) THEN
                 WRITE(6,*) 'DYSEV error: ',Info
                 STOP
             ENDIF

!Print out first few eigenvalues - want lowest ones
             WRITE(6,*) "Lowest eigenvalues are: "
             do i=1,toprint
                 WRITE(6,*) i,WList(i)
             enddo
!Print out ground state eigenvector
             OPEN(13,file='WAVEVECTOR',Status='unknown')
             do i=1,toprint
                 WRITE(13,*) i,HMat(i,1)
             enddo
             CLOSE(13)

             DEALLOCATE(HMat)
             CALL LogMemDealloc(this_routine,HMatTag)
             DEALLOCATE(WList)
             CALL LogMemDealloc(this_routine,WListTag)
             DEALLOCATE(Work)
             CALL LogMemDealloc(this_routine,WorkTag)
 
         ENDIF

!Set the maximum number of walkers allowed
         MaxWalkers=int(1000*InitWalkers,sizeof_int)

!Allocate memory to hold walkers
         ALLOCATE(WalkVec(MaxWalkers),stat=ierr)
         CALL LogMemAlloc('WalkVec',MaxWalkers,5,this_routine,WalkVecTag,ierr)
         ALLOCATE(WalkVec2(MaxWalkers),stat=ierr)
         CALL LogMemAlloc('WalkVec2',MaxWalkers,5,this_routine,WalkVec2Tag,ierr)

         IF(TReadPops) THEN
             IF((ABS(ScaleWalkers-1.0_dp)).lt.1.0e-8_dp) THEN
!Read in walker positions
                 do i=1,int(InitWalkers,sizeof_int)
                     READ(17,*) WalkVec(i)%Det,WalkVec(i)%WSign
                 enddo
             ELSE
!Read in walker positions - we will scale these later...
                 do i=1,int(InitWalkers,sizeof_int)
                     READ(17,*) WalkVec2(i)%Det,WalkVec2(i)%WSign
                 enddo
                 WRITE(6,*) "Scaling number of walkers by: ",ScaleWalkers
                 ReadWalkers=int(InitWalkers,sizeof_int)
                 InitWalkers=0
                 ! First, count the total number of initial walkers on each 
                 ! determinant - sort into list
                 call sort (WalkVec2(1:ReadWalkers))

                 j=1
!j is the counter over all read in walkers - it indicates when we have reached the end of the entire list
                 DetCurr=WalkVec2(j)%Det
!DetCurr is the current determinant
                 do while(j.le.ReadWalkers)
!Loop over all walkers
                     TotWalkersDet=0
                     do while ((WalkVec2(j)%Det.eq.DetCurr).and.(j.le.ReadWalkers))
!Loop over all walkers on DetCurr and count residual number after cancelling
                         IF(WalkVec2(j)%WSign) THEN
                             TotWalkersDet=TotWalkersDet+1
                         ELSE
                             TotWalkersDet=TotWalkersDet-1
                         ENDIF
                         j=j+1
                     enddo
                     DetCurr=WalkVec2(j)%Det
!Count total number of initial walkers
                     InitWalkers=InitWalkers+abs(nint((TotWalkersDet+0.0_dp)*ScaleWalkers))
                 enddo
                 WRITE(6,*) "Total number of walkers is now: ",InitWalkers
!Set the new maximum number of walkers allowed
                 MaxWalkers=int(100*InitWalkers,sizeof_int)

!Deallocate old memory block for WalkVec
                 DEALLOCATE(WalkVec)
                 CALL LogMemDealloc(this_routine,WalkVecTag)

!Allocate memory to hold new maximum number of walkers
                 ALLOCATE(WalkVec(MaxWalkers),stat=ierr)
                 CALL LogMemAlloc('WalkVec',MaxWalkers,5,this_routine,WalkVecTag,ierr)

!Now multiply them up...
                 j=1
                 VecSlot=1
!j is the counter over all read in walkers - it indicates when we have reached the end of the entire list
                 DetCurr=WalkVec2(j)%Det
!DetCurr is the current determinant
                 do while(j.le.ReadWalkers)
!Loop over all walkers
                     TotWalkersDet=0
                     do while ((WalkVec2(j)%Det.eq.DetCurr).and.(j.le.ReadWalkers))
!Loop over all walkers on DetCurr and count residual number after cancelling
                         IF(WalkVec2(j)%WSign) THEN
                             TotWalkersDet=TotWalkersDet+1
                         ELSE
                             TotWalkersDet=TotWalkersDet-1
                         ENDIF
                         j=j+1
                     enddo
!Now multiply up the number of walkers, and insert into WalkVec
                     TotWalkersDet=nint((TotWalkersDet+0.0_dp)*ScaleWalkers)
                     IF(TotWalkersDet.gt.0) THEN
                         do l=1,abs(TotWalkersDet)
                             WalkVec(VecSlot)%Det=DetCurr
                             WalkVec(VecSlot)%WSign=.true.
                             VecSlot=VecSlot+1
                         enddo
                     ELSE
                         do l=1,abs(TotWalkersDet)
                             WalkVec(VecSlot)%Det=DetCurr
                             WalkVec(VecSlot)%WSign=.false.
                             VecSlot=VecSlot+1
                         enddo
                     ENDIF
                     DetCurr=WalkVec2(j)%Det
                 enddo
                 IF((VecSlot-1).ne.InitWalkers) THEN
                     WRITE(6,*) "Problem scaling up walker number - exiting..."
                     STOP 'Problem scaling up walker number - exiting...'
                 ENDIF

!Now deallocate and reallocate WalkVec2 with correct number of total walkers
                 DEALLOCATE(WalkVec2)
                 CALL LogMemDealloc(this_routine,WalkVec2Tag)
                 ALLOCATE(WalkVec2(MaxWalkers),stat=ierr)
                 CALL LogMemAlloc('WalkVec2',MaxWalkers,5,this_routine,WalkVec2Tag,ierr)

             ENDIF

!End of reading in POPSFILE
             CLOSE(17)
         ENDIF

!Setup initial trial walker position - already they are all set to be at HF to start (with positive sign)
!TotWalkers contains the number of current walkers at each step
         TotWalkers=int(InitWalkers,sizeof_int)
         TotWalkersOld=int(InitWalkers,sizeof_int)

         IF(TCalcWavevector) THEN
!If we want to write out first few
             WRITE(14,*) 1,1.0_dp
             do i=2,toprint
                 WRITE(14,*) i,0.0_dp
             enddo
             WRITE(14,*) ""
         ENDIF
         
!Print out initial starting configurations
         WRITE(6,*) ""
         WRITE(6,*) "       Step  Shift  WalkerChange  GrowRate  TotWalkers"
         WRITE(15,*) "#       Step  Shift  WalkerChange  GrowRate  TotWalkers"
         IF(TReadPops) THEN
             WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7)") PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,0.0_dp
             WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7)") PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,0.0_dp
         ELSE
             WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,0.0_dp
             WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,0.0_dp
         ENDIF

!Start MC run - NMCyc indicates the number of times to run through all walkers
         do i=1,NMCyc

!VecSlot indicates the next free position in WalkVec2
             VecSlot=1

             do j=1,TotWalkers
!j runs through all current walkers

                 IF((WalkVec(j)%Det).eq.1) THEN

                     SumNoatHF=SumNoatHF+1
                     
!We are at HF - treat this walker slightly differently, since it is attached to all excits
!Run through all double excits and determine whether to create
!Change, so that now we are only selecting a single excitation at random - 
!if we do this, then we should increase prob by No.Excits
!                     r=Ran2(Seed)
!                     r=r*(NList-2)+2
                     do k=2,NList
!                     do k=int(r),int(r)
!Prob of creating a new walker on the excit is given by tau*abs(Hij)
                         rat=tau*abs(List(k,2))

                         IF(rat.gt.Ran2(Seed)) THEN
!Determine sign, and create new walker at k
                             WalkVec2(VecSlot)%Det=k
!If product of signs is +ve, create -ve walker, else create +ve walker
                             IF(WalkVec(j)%WSign) THEN
!Walker is positive
                                 IF(List(k,2).gt.0.0_dp) THEN
                                     WalkVec2(VecSlot)%WSign=.false. !-ve walker
                                 ELSE
                                     WalkVec2(VecSlot)%WSign=.true. !+ve walker
                                 ENDIF
                             ELSE
!Walker is negative
                                 IF(List(k,2).gt.0.0_dp) THEN
                                     WalkVec2(VecSlot)%WSign=.true. !+ve walker
                                 ELSE
                                     WalkVec2(VecSlot)%WSign=.false. !-ve walker
                                 ENDIF
                             ENDIF
!Increase walker number, and the next available slot
                             VecSlot=VecSlot+1
                         ENDIF

                     enddo

                 ELSE
!We are at an excitation - only possibility is to create walker at HF

                     DetCurr=WalkVec(j)%Det
                     
                     IF(WalkVec(j)%WSign) THEN
                         SumENum=SumENum+List(DetCurr,2)
                     ELSE
                         SumENum=SumENum-List(DetCurr,2)
                     ENDIF
                     rat=tau*abs(List(DetCurr,2))

                     IF(rat.gt.Ran2(Seed)) THEN
!Create new walker at HF
                         WalkVec2(VecSlot)%Det=1
!If product of signs is +ve, create -ve walker, else create +ve walker
                         IF(WalkVec(j)%WSign) THEN
!Walker is positive
                             IF(List(DetCurr,2).gt.0.0_dp) THEN
                                 WalkVec2(VecSlot)%WSign=.false. !-ve walker
                             ELSE
                                 WalkVec2(VecSlot)%WSign=.true. !+ve walker
                             ENDIF
                         ELSE
!Walker is negative
                             IF(List(DetCurr,2).gt.0.0_dp) THEN
                                 WalkVec2(VecSlot)%WSign=.true. !+ve walker
                             ELSE
                                 WalkVec2(VecSlot)%WSign=.false. !-ve walker
                             ENDIF
                         ENDIF
!Increase walker number, and the next available slot
                         VecSlot=VecSlot+1
                     ENDIF

!We have finished looking at excitations of walker determinant
                 ENDIF

!Next we have to decide if the walker wants to self-destruct or not
!Kill with prob (Hii-DiagSft)*tau - this should ALWAYS be positive
                 rat=tau*(List((WalkVec(j)%Det),0)-DiagSft)
                 IF(Ran2(Seed).gt.rat) THEN
!This walker is spared - copy him across to the new WalkVec - in the same position
!If the walker isn't copied across, it has self-destructed
                     WalkVec2(VecSlot)%Det=WalkVec(j)%Det
                     WalkVec2(VecSlot)%WSign=WalkVec(j)%WSign
                     VecSlot=VecSlot+1
                 ENDIF

!Finish cycling over walkers
             enddo

!Since VecSlot holds the next vacant slot in the array, TotWalkersNew will be one less than this.
             TotWalkersNew=VecSlot-1
             rat=(TotWalkersNew+0.0_dp)/(MaxWalkers+0.0_dp)
             IF(rat.gt.0.9) THEN
                 WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkers"
             ENDIF

!In this next section, we annihilate pairs of walkers on the same determinant with opposing signs
!This does not have to be performed every cycle - however, we have to at the moment,since it is the only way we transfer the determinants between arrays
!              IF(mod(i,1).eq.0) THEN
              IF(.true.) THEN

                  IF(.NOT.TBinCancel) THEN
                      ! This method of canceling down pairs of opposite sign 
                      ! walkers is based on ordering the list of determinants
                      ! Once ordered, each block of walkers on similar 
                      ! determinants can be analysed, and the residual walker
                      ! concentration moved to WalkVec First need to order 
                      ! the vector of new walkers according to determinant 
                      ! the walker is on - this is an NlogN scaling operation.
                      call sort (WalkVec2(1:TotWalkersNew))

                      j=1
!j is the counter over all TotWalkersNew - it indicates when we have reached the end of the entire list
                      DetCurr=WalkVec2(j)%Det
!DetCurr is the current determinant
                      VecSlot=1
                      do while(j.le.TotWalkersNew)
!Loop over all walkers
                          TotWalkersDet=0
                          do while ((WalkVec2(j)%Det.eq.DetCurr).and.(j.le.TotWalkersNew))
!Loop over all walkers on DetCurr and count residual number after cancelling
                              IF(WalkVec2(j)%WSign) THEN
                                  TotWalkersDet=TotWalkersDet+1
                              ELSE
                                  TotWalkersDet=TotWalkersDet-1
                              ENDIF
                              j=j+1
                          enddo
!Transfer residual population into VecSlot, along with residual sign
                          IF(TotWalkersDet.gt.0) THEN
                              do l=1,abs(TotWalkersDet)
                                  WalkVec(VecSlot)%Det=DetCurr
                                  WalkVec(VecSlot)%WSign=.true.
                                  VecSlot=VecSlot+1
                              enddo
                          ELSE
                              do l=1,abs(TotWalkersDet)
                                  WalkVec(VecSlot)%Det=DetCurr
                                  WalkVec(VecSlot)%WSign=.false.
                                  VecSlot=VecSlot+1
                              enddo
                          ENDIF
!Update new current determinant
                          DetCurr=WalkVec2(j)%Det
                      enddo
!The new number of cancelled determinants is given by one less than VecSlot again.
                      TotWalkers=VecSlot-1

                      IF(TCalcWavevector.and.(mod(i,WavevectorPrint).eq.0)) THEN
!We want to compare the evolution of wavefunction to exact wavefunction, so calculate it every WavevectorPrint
!Use the redundant List(i,1) vector to store wavefunctions, depicted by the residual concentration of walkers - zero it
                          do j=1,NList
                              List(j,1)=0.0_dp
                          enddo
                          do j=1,TotWalkers
!Run through all walkers, and bin walkers in the correct wavevector component
                              IF(WalkVec(j)%WSign) THEN
                                  List(WalkVec(j)%Det,1)=List(WalkVec(j)%Det,1)+1.0_dp
                              ELSE
                                  List(WalkVec(j)%Det,1)=List(WalkVec(j)%Det,1)-1.0_dp
                              ENDIF
                          enddo
!Now, normalise wavevector
                          Norm=0.0_dp
                          do j=1,NList
                              Norm=Norm+(List(j,1)**2)
                          enddo
                          Norm=SQRT(Norm)
                          do j=1,NList
                              List(j,1)=List(j,1)/Norm
                          enddo
                          do j=1,toprint
                              WRITE(14,*) j,List(j,1)
                          enddo
                          WRITE(14,*) ""

                      ENDIF

                  ELSE
!In this method of cancelling, we simply histogram the walkers according to the determinant they're on - scales with size of system - i.e. badly
!Use the redundant List(i,1) vector to store wavefunctions, depicted by the residual concentration of walkers - zero it
                      do j=1,NList
                          List(j,1)=0.0_dp
                      enddo
                
                      do j=1,TotWalkersNew
!Run through all walkers
                        
                          IF((WalkVec2(j)%Det).gt.NList) THEN
                              WRITE(6,*) "Serious problem here..."
                              STOP 'Serious problem here...'
                          ENDIF

                          IF(WalkVec2(j)%WSign) THEN
!Walker is positive - add to determinant contribution
                              List((WalkVec2(j)%Det),1)=List((WalkVec2(j)%Det),1)+1.0_dp
                          ELSE
!Walker is negative - subtract from determinant contribution
                              List((WalkVec2(j)%Det),1)=List((WalkVec2(j)%Det),1)-1.0_dp
                          ENDIF

                      enddo

                      VecSlot=1
                      Norm=0.0_dp
                      do j=1,NList
!Run through all determinants - normalise, find new number of walkers, and find new WalkVec
                          NWalk=nint(List(j,1))
                          IF(NWalk.ne.0) THEN
!Norm will be used to normalise the eigenvector
                              Norm=Norm+(List(j,1)**2)
                              IF((List(j,1)/abs(List(j,1))).lt.0.0_dp) THEN
!Component has overall negative sign
                                  DetSign=.false.
                              ELSE
                                  DetSign=.true.
                              ENDIF
                              do k=1,abs(NWalk)
                                  WalkVec(VecSlot)%Det=j
                                  WalkVec(VecSlot)%WSign=DetSign
                                  VecSlot=VecSlot+1
                              enddo

                          ENDIF

                      enddo

                      IF(TCalcWavevector.and.(mod(i,Wavevectorprint).eq.0)) THEN
!If we want to, we can renormalise the trial vector - only do this every 50 cycles
                          Norm=SQRT(Norm)
                          do j=1,NList
                              List(j,1)=List(j,1)/Norm
                          enddo
                          do j=1,toprint
                              WRITE(14,*) j,List(j,1)
                          enddo
                          WRITE(14,*) ""
                      ENDIF

!The new number of cancelled determinants is given by one less than VecSlot again.
                      TotWalkers=VecSlot-1

                  ENDIF

!End of walker cancellation
             ENDIF

!             WRITE(6,"(I9,G16.7,I9,G16.7,I9)") i,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                    
!Now we want to calculate the change in the shift of the diagonal elements...
!We only want to change every so often, so that walker numbers have a chance to acclimatise between shift change
!             StepsSft=100
             IF(mod(i,StepsSft).eq.0) THEN

                 ProjectionE=(SumENum/(SumNoatHF+0.0_dp))
                 
                 GrowRate=(TotWalkers+0.0_dp)/(TotWalkersOld+0.0_dp)
                 DiagSft=DiagSft-(log(GrowRate)*SftDamp)/(Tau*(StepsSft+0.0_dp))

!Write out MC cycle number, Shift, Change in Walker no, Growthrate, New Total Walkers
                 IF(TReadPops) THEN
                     WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7)") i+PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld, &
                        GrowRate,TotWalkers,ProjectionE
                     WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7)") i+PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld, &
                        GrowRate,TotWalkers,ProjectionE
                 ELSE
                     WRITE(15,"(I9,G16.7,I9,G16.7,I9,G16.7)") i,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE
                     WRITE(6,"(I9,G16.7,I9,G16.7,I9,G16.7)") i,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers,ProjectionE
                 ENDIF
                 CALL neci_flush(15)
                 CALL neci_flush(6)
                 IF((DiagSft).gt.0.0_dp) THEN
                     WRITE(6,*) "***WARNING*** - DiagSft trying to become positive...",DiagSft
                     STOP
                 ENDIF
                 TotWalkersOld=TotWalkers
             ENDIF

!Finish MC Cycle
         enddo

         IF(TPopsFile) THEN
!Print out current state of simulation, so it can be restarted if desired...
             call get_unique_filename('POPSFILE',tIncrementPops,.true.,iPopsFileNoWrite,popsfile)
             OPEN(16,file=popsfile,status='unknown')
             WRITE(16,*) TotWalkers, "   TOTWALKERS"
             WRITE(16,*) DiagSft, "   DIAG SHIFT"
             IF(TReadPops) THEN
                 WRITE(16,*) NMCyc+PreviousNMCyc, "   NO. CYCLES"
             ELSE
                 WRITE(16,*) NMCyc, "   MC CYCLES"
             ENDIF
             do i=1,TotWalkers
                 WRITE(16,*) WalkVec(i)%Det,WalkVec(i)%WSign
             enddo
             CLOSE(16)
         ENDIF

         SI=1.0_dp
         DLWDB=-2.0_dp*DiagSft

!Deallocate memory
         DEALLOCATE(WalkVec)
         CALL LogMemDealloc(this_routine,WalkVecTag)
         DEALLOCATE(WalkVec2)
         CALL LogMemDealloc(this_routine,WalkVec2Tag)

         IF(TCalcWavevector) CLOSE(14)
         CLOSE(15)

         call halt_timer(proc_timer)

         RETURN

      END SUBROUTINE StarDiagMC

!Use a Lanczos routine to find the first NEval eigenvectors of the rho matrix, and from this the energy of the star graph.
      SUBROUTINE StarDiagLanc(NList,List,ILMax,i_P,SI,DBeta,DLWDB)
         use constants, only: dp
         USE DetCalcData , only : B2L,nBlk,nKry,NEval
         USE global_utilities
         use HElem
        use MemoryManager, only: TagIntType
         IMPLICIT NONE
         CHARACTER(len=*), PARAMETER :: this_routine='StarDiagLanc'
         INTEGER :: i_P,NList,ILMax,i,j,LenMat,ICMax
         type(timer), save :: proc_timer
         real(dp) :: List(ILMax,0:2)
         real(dp) :: SI,DLWDB,DBeta
         LOGICAL :: TSeeded
         INTEGER :: NCycle,NBlock,NKry1,LScr,LIScr,ierr
         INTEGER(TagIntType) :: V2Tag=0,WTag=0,CKTag=0,CKNTag=0,AMTag=0,BMTag=0,Work2Tag=0
         INTEGER(TagIntType) :: TTag=0,SCRTag=0,ISCRTag=0,IndexTag=0,WHTag=0
         INTEGER(TagIntType) :: LabTag=0,ATag=0,MatTag=0,NRowTag=0,VTag=0,WTTag=0
         real(dp) , ALLOCATABLE :: Work2(:),WH(:),V2(:,:),W(:),A(:,:),V(:)
         real(dp) , ALLOCATABLE :: AM(:),BM(:),T(:),WT(:),SCR(:)
         real(dp) , ALLOCATABLE :: Mat(:),CK(:,:),CKN(:,:)
         INTEGER , ALLOCATABLE :: Lab(:),NRow(:),ISCR(:),Index(:)
         
         IF(HElement_t_size.GT.1) THEN
             CALL Stop_All("StarDiagLanc","StarDiagLanc cannot function with complex orbitals.")
         ENDIF

         proc_timer%timer_name='StarDiagLanc'
         call set_timer(proc_timer)

         LenMat=(NList*2)-1
         ICMax=NList

         ALLOCATE(Mat(LenMat),stat=ierr)
         CALL LogMemAlloc('Mat',LenMat,8,this_routine,MatTag,ierr)
         Mat=0.0_dp
         ALLOCATE(Lab(LenMat),stat=ierr)
         CALL LogMemAlloc('Lab',LenMat,4,this_routine,LabTag,ierr)
         Lab(1:LenMat)=0
         ALLOCATE(NRow(NList),stat=ierr)
         CALL LogMemAlloc('NRow',NList,4,this_routine,NRowTag,ierr)
         NRow(1:NList)=0

!Fill compressed rho-matrix
         IF((List(1,1).ne.List(1,0)).or.(List(1,1).ne.1.0_dp)) THEN
             CALL Stop_All("StarDiagLanc","Error with rho matrix elements")
         ENDIF
         NRow(1)=NList
         Mat(1)=List(1,1)
         Lab(1)=1
         do i=2,NList
             Mat(i)=List(i,1)
             Mat(NList+(i-1))=List(i,0)
             Lab(i)=i
             Lab(NList+(i-1))=i
             NRow(i)=1
         enddo

!Set up parameters for the diagonaliser.
!NEval indicates the number of eigenvalues we want to calculate
         IF(NEval.eq.0) THEN
!NEval is set to 0 by default, computing all eigenvectors
             WRITE(6,*) "Number of EIGENVALUES not specified. Computing all eigenvalues (slow)."
             NEval=NList
         ELSEIF(NEval.gt.NList) THEN
             WRITE(6,*) "Number of EIGENVALUES set to largest number than of connected doubles."
             WRITE(6,*) "Resetting number of EIGENVALUES to size of matrix."
             NEval=NList
         ELSE
             WRITE(6,"(A,I6,A)") "Calculating ",NEval," largest EIGENVALUES..."
         ENDIF
         NCycle=200
         NKry1=NKry+1
         NBlock=MIN(NEval,NBlk)
         LScr=MAX(NList*NEval,8*NBlock*NKry)
         LIScr=6*NBlock*NKry

!Allocate memory for diagonaliser
         ALLOCATE(A(NEval,NEval),stat=ierr)
         CALL LogMemAlloc('A',NEval*NEval,8,this_routine,ATag,ierr)
         A=0.0_dp
         ALLOCATE(V(NList*NBlock*NKry1),stat=ierr)
         CALL LogMemAlloc('V',NList*NBlock*NKry1,8,this_routine,VTag,ierr)
         V=0.0_dp
         ALLOCATE(AM(NBlock*NBlock*NKry1),stat=ierr)
         CALL LogMemAlloc('AM',NBlock*NBlock*NKry1,8,this_routine,AMTag,ierr)
         AM=0.0_dp
         ALLOCATE(BM(NBlock*NBlock*NKry),stat=ierr)
         CALL LogMemAlloc('BM',NBlock*NBlock*NKry,8,this_routine,BMTag,ierr)
         BM=0.0_dp
         ALLOCATE(T(3*NBlock*NKry*NBlock*NKry),stat=ierr)
         CALL LogMemAlloc('T',NBlock*NKry*NBlock*NKry*3,8,this_routine,TTag,ierr)
         T=0.0_dp
         ALLOCATE(WT(NBlock*NKry),stat=ierr)
         CALL LogMemAlloc('WT',NBlock*NKry,8,this_routine,WTTag,ierr)
         WT=0.0_dp
         ALLOCATE(SCR(LScr),stat=ierr)
         CALL LogMemAlloc('SCR',LScr,8,this_routine,SCRTag,ierr)
         SCR=0.0_dp
         ALLOCATE(ISCR(LIScr),stat=ierr)
         CALL LogMemAlloc('ISCR',LIScr,4,this_routine,ISCRTag,ierr)
         ISCR(1:LIScr)=0
         ALLOCATE(Index(NEval),stat=ierr)
         CALL LogMemAlloc('Index',NEval,4,this_routine,IndexTag,ierr)
         Index(1:NEval)=0
         ALLOCATE(WH(NList),stat=ierr)
         CALL LogMemAlloc('WH',NList,8,this_routine,WHTag,ierr)
         WH=0.0_dp
         ALLOCATE(Work2(3*NList),stat=ierr)
         CALL LogMemAlloc('Work2',3*NList,8,this_routine,Work2Tag,ierr)
         Work2=0.0_dp
         ALLOCATE(V2(NList,NEval),stat=ierr)
         CALL LogMemAlloc('V2',NList*NEval,8,this_routine,V2Tag,ierr)
         V2=0.0_dp
 
!W holds the eigenvalues 
         ALLOCATE(W(NEval),stat=ierr)
         CALL LogMemAlloc('W',NEval,8,this_routine,WTag,ierr)
         W=0.0_dp

!CK holds the eigenvectors
         ALLOCATE(CK(NList,NEval),stat=ierr)
         CALL LogMemAlloc('CK',NList*NEval,8,this_routine,CKTag,ierr)
!The initial trial wavefuntion is set to zero        
         CK=0.0_dp
         ALLOCATE(CKN(NList,NEval),stat=ierr)
         CALL LogMemAlloc('CKN',NList*NEval,8,this_routine,CKNTag,ierr)
         CKN=0.0_dp

         TSeeded=.false.

         CALL NECI_FRSBLKH(NList,ICMax,NEval,Mat,Lab,CK,CKN,NKry,NKry1,NBlock,NRow,LScr,LIScr,A,W,V,AM,BM,T,WT,SCR, &
            ISCR,Index,NCycle,B2L,.false.,.true.,TSeeded)
         
!Deallocate memory required by diagonaliser (including original matrix)
         DEALLOCATE(Mat)
         CALL LogMemDealloc(this_routine,MatTag)
         DEALLOCATE(Lab)
         CALL LogMemDealloc(this_routine,LabTag)
         DEALLOCATE(A)
         CALL LogMemDealloc(this_routine,ATag)
         DEALLOCATE(V)
         CALL LogMemDealloc(this_routine,VTag)
         DEALLOCATE(AM)
         CALL LogMemDealloc(this_routine,AMTag)
         DEALLOCATE(BM)
         CALL LogMemDealloc(this_routine,BMTag)
         DEALLOCATE(T)
         CALL LogMemDealloc(this_routine,TTag)
         DEALLOCATE(WT)
         CALL LogMemDealloc(this_routine,WTTag)
         DEALLOCATE(SCR)
         CALL LogMemDealloc(this_routine,ScrTag)
         DEALLOCATE(WH)
         CALL LogMemDealloc(this_routine,WHTag)
         DEALLOCATE(Work2)
         CALL LogMemDealloc(this_routine,Work2Tag)
         DEALLOCATE(V2)
         CALL LogMemDealloc(this_routine,V2Tag)
 
         WRITE(6,*) "Highest root: ",W(1)
         SI=0.0_dp

!         do i=1,NEval
!             WRITE(2,*) W(NEval-(i-1))
!         enddo

!         do i=1,NList
!This prints out the largest eigenvector
!             WRITE(3,"(F22.16)") CK(i,1)
!         enddo

!Eigenvectors are stored in CK(i,x) where i runs over the dimension of the matrix, and 
!x runs over all the eigenvectors (NEval), where 1 is the eigenvector corresponding to
!infinite temperature (i.e. largest eigenvalue)
         do i=1,NEval
             SI=SI+(CK(1,i)**2)*(W(i)**i_P)
             IF(DBeta.ne.0.0_dp) THEN
!Calculate <D|H exp(-b H)|D>/Rho_ii^P
                 do j=1,NList
                     DLWDB=DLWDB+List(j,2)*CK(j,i)*(W(i)**i_P)*CK(1,i)
                 enddo
             ENDIF
         enddo

         WRITE(6,*) "Final SI= ",SI
         SI=SI-1.0_dp
         DLWDB=DLWDB-LIST(1,2)


         call halt_timer(proc_timer)
         RETURN

      END SUBROUTINE StarDiagLanc
        
      
      SUBROUTINE STARDIAG(NLIST,LIST,ILMAX,I_P,SI,DBETA,DLWDB)
         use constants, only: dp,int32
         use IntegralsData , only : TCalcRealProd
         use global_utilities
         use sort_mod
         use helem, only: helement_t_size
        use MemoryManager, only: TagIntType
         IMPLICIT NONE
         INTEGER I_P
         INTEGER NLIST,ILMAX
         real(dp) LIST(ILMAX,0:2)
         real(dp),ALLOCATABLE ::  RIJMAT(:),WLIST(:),WORK(:)
         INTEGER(TagIntType), SAVE :: tagRIJMAT=0,tagWLIST=0,tagWORK=0
         type(timer), save :: proc_timer
         INTEGER WORKL
         INTEGER(int32) INFO
         real(dp) SI,DLWDB,DBETA,OD
         INTEGER I,J,err
         character(*),parameter :: this_routine='STARDIAG'
         
         IF(HElement_t_size.GT.1) THEN
             CALL Stop_All("StarDiag","STARDIAG cannot function with complex orbitals.")
         END IF

         proc_timer%timer_name='STARDIAG  '
         call set_timer(proc_timer)
         allocate(RIJMAT(NLIST**2),stat=err)
         call LogMemAlloc('RIJMAT',NLIST**2,8,this_routine,tagRIJMAT,err)
         allocate(WLIST(NLIST),stat=err)
         call LogMemAlloc('WLIST',NLIST,8,this_routine,tagWLIST,err)
         WORKL=3*NLIST
         allocate(WORK(WORKL),stat=err)
         call LogMemAlloc('WORK',WORKL,8,this_routine,tagWORK,err)

         ! Requires real integrals (not assuming helements)
         call sort (list(2:nList-1,0), list(2:nList-1,1), list(2:nList-1,2))

         RIJMAT=0.0_dp
!.. Now we fill the RIJ array
         DO I=0,NLIST-1
            RIJMAT(I*NLIST+I+1)=LIST(I+1,0)
            RIJMAT(I+1)=LIST(I+1,1)
         ENDDO

!         WRITE(67,*) "Size of rhomat is ", NLIST
!         WRITE(67,*) "********"
!         WRITE(67,*) "OFF DIAG LIST is "
!         DO I=1,ILMAX
!            WRITE(67,"E14.6,$") LIST(I,1)
!         ENDDO
!         WRITE(67,*) ""
!         WRITE(67,*) "********"
!         WRITE(67,*) "ON DIAG LIST is "
!         DO I=1,ILMAX
!            WRITE(67,"E14.6,$") LIST(I,0)
!         ENDDO
!         WRITE(67,*) ""
!         WRITE(67,*) "********"
!         WRITE(67,*) "RHOMAT is "
!         DO I=1,NLIST
!            DO J=1,NLIST
!                WRITE(67,"E14.6,$") RIJMAT(((I-1)*NLIST)+J)
!            ENDDO
!            WRITE(67,*) ""
!        ENDDO
!        WRITE(67,*) "************"

!.. Diagonalize
         CALL DSYEV('V','L',NLIST,RIJMAT,NLIST,WLIST,WORK,WORKL,INFO)
         IF(INFO.NE.0) THEN
            WRITE(6,*) 'DYSEV error: ',INFO
            STOP
         ENDIF
         deallocate(WORK)
         call LogMemDealloc(this_routine,tagWORK)
         WRITE(6,*)
         WRITE(6,*) "Highest root:",WLIST(NLIST)
!.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
         SI=0.0_dp
        
!         WRITE(67,*) "Eigenvalues are: "
!         DO I=1,NLIST
!            WRITE(67,"F22.16,$") WLIST(I)
!            WRITE(67,"(F22.16)") WLIST(I)
!         ENDDO

!         DO I=1,NLIST
!            WRITE(68,"(F22.16)") RIJMAT(((NLIST-1)*NLIST)+I)
!         ENDDO

!Divide through by largest eigenvalue to prevent blowing up in some cases
         IF(TCalcRealProd) THEN
            DO I=1,NLIST
                WLIST(I)=WLIST(I)/WLIST(NLIST)
            ENDDO
         ENDIF
         
!         DO I=0,NLIST-1
!            WRITE(6,*) WLIST(I+1)
!         ENDDO
!         DO I=NLIST-1,0,-1
!            WRITE(6,*) I+1,RIJMAT(I*NLIST+1)*RIJMAT(I*NLIST+1)*(WLIST(I+1)**I_P),RIJMAT(I*NLIST+1)*RIJMAT(I*NLIST+1)
!         ENDDO
!         DO I=NLIST,1,-1
!            RR=1.0_dp
!            IF(I.LT.NLIST) RR=(WLIST(I)-LIST(I+1,0))
!            write(6,"(I,3G)") I-1,WLIST(I),LIST(I+1,0),RR
!            IF(abs(RR).gt.1.0e-10_dp)  WRITE(6,*) 1/(RIJMAT((I-1)*NLIST+1)**2),RIJMAT((I-1)*NLIST+1)**2
!         ENDDO
         DO I=0,NLIST-1
            SI=SI+RIJMAT(I*NLIST+1)*RIJMAT(I*NLIST+1)*(WLIST(I+1)**I_P)
            IF(DBETA.NE.0.0_dp) THEN
!.. calculate <D|H exp(-b H)|D>/RHO_ii^P
               OD=DLWDB
               DO J=1,NLIST
                  DLWDB=DLWDB+LIST(J,2)*RIJMAT(I*NLIST+J)*RIJMAT(I*NLIST+1)*(WLIST(I+1)**I_P)
!                  WRITE(6,*) LIST(J,2),RIJMAT(I*NLIST+J)
               ENDDO
!               WRITE(6,*) WLIST(I+1)**I_P,DLWDB-OD
!               WRITE(6,*)
            ENDIF
         ENDDO
         WRITE(6,*) "Final SI=",SI
         SI=SI-1.0_dp
         DLWDB=DLWDB-LIST(1,2)
         deallocate(WLIST,RIJMAT)
         call LogMemDealloc(this_routine,tagWLIST)
         call LogMemDealloc(this_routine,tagRIJMAT)
         call halt_timer(proc_timer)
         RETURN
      END


!.. Use an iterative Order(N) root-finding method to diagonalize the
!.. star matrix.
!  LIST(0:NLIST-1,:) contains data
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ
      SUBROUTINE STARDIAG2(NLIST,LIST,ILMAX,I_P,SI,DBETA,DLWDB,NROOTS)
         use CalcData , only : STARCONV
         use IntegralsData , only : TQUADRHO,TEXPRHO
         use constants, only: dp
         use global_utilities
         use sort_mod
         IMPLICIT NONE
         INTEGER I_P
         INTEGER NLIST,ILMAX
         HElement_t LIST(0:ILMAX-1,0:2)
         type(timer), save :: proc_timer
         real(dp) SI,DLWDB,DBETA,NORM,E0
         HElement_t DLWDB2,RR
         INTEGER I,J,NROOTS
         real(dp) ROOTS(0:NROOTS),RPN
         INTEGER iEigv,iDegen
         LOGICAL lWarned
         real(dp) NORMCHECK,NORMROOTS

         proc_timer%timer_name='STARDIAG2 '
         call set_timer(proc_timer)
!.. we need to sort A and B (and the list of hamil values) into ascending A order
!         WRITE(6,*) (LIST(I,2),I=1,NLIST)
!         WRITE(6,*) (LIST(I,1),I=1,NLIST)
         call sort (list(1:nlist-1,0), list(1:nlist-1,1), list(1:nlist-1,2))
!         WRITE(6,*) (LIST(I,2),I=1,NLIST)
!         WRITE(6,*) (LIST(I,1),I=1,NLIST)

!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ
         if(nRoots.eq.nList) then
!  nRoots is one more than the number of excitations
!we've been asked to search to see how many roots we need for convergence.
            i=nRoots-1
            do while (i.gt.0.and.abs((nRoots-i)*(List(i,0)**I_P)).ge.STARCONV)
               i=i-1
            enddo
            nRoots=nRoots-1-i
            write(6,*) nRoots+1, " needed for convergence to", STARCONV
            nRoots=nRoots+1
         elseif(nRoots.eq.nList+1) then
!  nRoots is one more than the number of excitations
!we've been asked to search to see how many roots we need for convergence. - Method 2
!  Take into account the cumulative values as we go down the list, not just the absolute value of the rest
            i=nRoots-1
            si=1
            do while (i.gt.0.and.abs((nRoots-i)*(List(i,0)**I_P)/si).ge.1.0e-2_dp)
               Si=SI+List(i,0)**I_P
               i=i-1
            enddo
            nRoots=nRoots-1-i
            write(6,*) nRoots+1, " needed for method 2 convergence 1.0e-3_dp."
            nRoots=nRoots+1
         endif

!.. Find the eigenvalues
!  NLIST is the number of elements in the list, but we need to give the index of the last element, NLIST-1
!         CALL PLOTROOTSTAR(NLIST-1,LIST(0,0),LIST(0,1)) 
         CALL FINDROOTSTAR(NLIST-1,LIST(0,0),LIST(0,1),ROOTS,NROOTS)
         SI=0.0_dp
!         WRITE(6,*) "ROOTS ARE: " 
!         WRITE(6,*) ROOTS(:)
         WRITE(6,*) "Highest root:",ROOTS(NROOTS)
         E0=List(0,2)
         lWarned=.false.
!.. divide through by highest eigenvalue to stop things blowing up
!         DO I=1,NROOTS
!            WRITE(6,*) ROOTS(I)
!,list(NLIST-NROOTS+I-1,0)
!         ENDDO
         iEigv=0 
         iDegen=0
!NLIST is the length of list, and the max possible value of  NROOTS
         NORMCHECK=0
         DO I=NROOTS,1,-1
            iDegen=iDegen+1
            RR=1.0_dp
            IF(I.LT.NROOTS) RR=(ROOTS(I))-LIST(NLIST-NROOTS+I,0)
!            write(6,"(I,3G)") I,ROOTS(I),LIST(NLIST-NROOTS+I,0),RR
            IF(ROOTS(I).EQ.LIST(NLIST-NROOTS+I-1,0).OR..NOT.abs(RR).gt.1.0e-13_dp) THEN
!.. If we're in a degenerate set of eigenvectors, we calculate things a
!.. little differently
!.. k is the vertex which the degeneracies couple to
!.. and j is the current degenerate element
!.. <Di|H|Dj><Dj|Psi> L**P <Psi|Di>
!.. <Dk|Psi>=rho_ij/NORM
!.. <Dj|Psi>=-rho_ik/NORM
!.. NORM=rho_ij**2+rho_ik**2
!.. <Di|Psi>=0 so we have no contributions at all!


            ELSE
!.. We need to calculate the normalization of each eigenvector
               iEigv=iEigv+1
               NORM=1.0_dp
               if(iEigv.ge.1.and.iEigv.le.3) then
!!!                  write(6,*) iDegen-1
                  iDegen=1
               endif
!List(0,:) is the HF det.  We set its value in the eigenvector to 1.  The remaining NLIST-1 items are used to normalize.
               DO J=1,NLIST-1
                  RR=(ROOTS(I))-LIST(J,0)
                  IF(.NOT.abs(RR).gt.1d-13) THEN
!see comment below
                     WRITE(6,"(A,I6,A,G25.16,A,G25.16,A,I6)") "WARNING: Eigenvalue I=",I,":",ROOTS(I), &
                        " dangerously close to rhojj=",abs(LIST(J,0))," J=",J
                     WRITE(6,"(A,I6,2G25.16)") "POLE,NUMER",J,abs(LIST(J,0)),abs(LIST(J,1))
                     WRITE(6,"(A,I6,2G25.16)") "POLE,NUMER",J-1,abs(LIST(J-1,0)),abs(LIST(J-1,1))
                  ENDIF
                  NORM=NORM+abs(LIST(J,1)/RR)**2
               ENDDO
!.. We add in the first element of the eigenvector * lambda**P
! As a test for size consistency - use expansion of exp(rho)-1 instead of normal rho matrix - c.f QCISD
               IF(TQUADRHO) THEN
                   !In order to stop it blowing up, divide through by r_max + r_max^2/2
                   NORMROOTS=(ROOTS(I)+(ROOTS(I)**2)/2)/(ROOTS(NROOTS)+(ROOTS(NROOTS)**2)/2)
                   RPN=(NORMROOTS**I_P)*1.0_dp/NORM
               ELSEIF(TEXPRHO) THEN
                   NORMROOTS=(EXP(ROOTS(I))-1)/(EXP(ROOTS(NROOTS))-1)
                   RPN=(NORMROOTS**I_P)*1.0_dp/NORM
               ELSE
                   RPN=(ROOTS(I)**I_P)*1.0_dp/NORM
               ENDIF
               
!               write(6,*) NORM,RPN
               IF(.not.lWarned.and.RPN/SI.LT.1.0e-4_dp) then
                  lWarned=.true.
!!!                  WRITE(6,*) "Root ",NROOTS-I," has low contribution."
!!!                  WRITE(6,*) "SI=",SI
               ENDIF
               SI=SI+RPN
               NORMCHECK=NORMCHECK+1/NORM
!               WRITE(6,*) I,RPN,1/NORM
               IF(iEigv.le.2) then
!!!                  write(6,"(A,I,A,2G)",advance='no') "Eigenvalue ",iEigv," = ",roots(i),E0-(i_P/Beta)*log(roots(i))
!                  write(6,*) "***",E0
               endif
               IF(DBETA.NE.0.0_dp) THEN
                  DLWDB2=LIST(0,2)
!                  WRITE(6,*) LIST(1,2),SQRT(1/NORM)
                  DO J=1,NLIST-1
#if __CMPLX
                     DLWDB2=DLWDB2+LIST(J,2)*(conjg(LIST(J,1))/((ROOTS(I))-LIST(J,0)))
#else
                     DLWDB2=DLWDB2+LIST(J,2)*((LIST(J,1))/((ROOTS(I))-LIST(J,0)))
#endif
!                WRITE(6,*) LIST(J,2),
!     &            LIST(J,1)/((ROOTS(I+1)-LIST(J,0))*SQRT(NORM))
                  ENDDO
                  DLWDB=DLWDB+DLWDB2*RPN
               ENDIF
!               WRITE(6,*) ROOTS(I+1)**I_P,DLWDB2*(ROOTS(I+1)**I_P)/NORM
!               WRITE(6,*)
            ENDIF
         ENDDO
               if(iEigv.gt.1.and.iEigv.le.3) then
!                  write(6,*) iDegen-1
                  iDegen=1
               endif
         IF(TQUADRHO) WRITE(6,*) "QUADRATIC EXPANSION OF RHO MATRIX USED"
         IF(TEXPRHO) WRITE(6,*) "EXPONENTIAL EXPANSION OF RHO MATRIX USED"
         write(6,*)
         WRITE(6,*) "Final SI=",SI
         IF(NROOTS.eq.(NLIST-1)) THEN
!We are searching for all roots, therefore sum of squares of projection onto root should be unity
             WRITE(6,*) "Norm of i projection:", NORMCHECK
             IF(ABS(NORMCHECK-1).gt.0.01) WRITE(6,*)  "WARNING: Norm differs from 1 by more than 0.01.  " &
                & //"Convergence may not be reached."
         ENDIF
         SI=SI-1.0_dp
         DLWDB=DLWDB-LIST(0,2)
         call halt_timer(proc_timer)
         RETURN
      END
!Comment on numerical stability  AJWT 30/10/07
!It seems that we're going to have to be very careful about numerical stability for the 2-vertex star if we wish to use anything but the highest eigenvalue.
!
!The Newton-Raphson problems are now contained, and if a NR step takes us out of the search region, we revert to regular falsi.
!The best precision we can achieve is a convergence of about 1e-15 on the eigenvalue.

!This means that we must be exceedingly careful when calculating the eigenvector, and in particular normalising it.

!For det j, its component in the eigenvector corresponding to root r_a is
!v_aj = v_ai * rho_ji/(r_a-rho_jj)

!where v_ai is the component of the HF det in the eigenvector (which we set to 1).  From this we can normalize the eigenvector.
!However, if rho_ji is very small (say 1e-14), then its coupling with i is very small - i.e. the eigenvector will likely have a large component of j, and a small component of i.
!  This will manifest itself as a near-degeneracy of the eigenvalue r_a to the pole rho_jj of the polynomial.
!  However, as our accuracy in calculating the eigenvectors is limited to a little more than machine precision, we cannot calculate r_a-rho_jj very precisely at all, and so the eigenvector normalization gets very inaccurate.

!I've put a test in to check when r_a-rho_jj is <1e-13, and it will print a warning in this case.
!To remedy the solution is simple.  SInce the coupling to j is small, we can ignore it.  This is done by setting rhoepsilon to a value like 1e-12.

!Alex



!-----


!        ADDSINGLES specifies to add the singles which are en-route to each double to that double as spokes, and prediagonalize them.
!  i.e. if the double is (ij->ab), then singles (i->a),(i->b),(j->a) and (j->b) are created in a star with (ij->ab), the result diagonalized, and the eigenvalues and vectors used to create new spokes.  Only works with NEW
      SUBROUTINE StarAddSingles(nI,nJ,ExcitInfo,iExcit,iMaxExcit,rhii,rhoeps,Beta,i_P,nEl,G1, &
        nBasis,nMsh,fck,nMax,ALat,UMat,nTay,ECore)
         use constants, only: dp,int32      
         use SystemData, only: BasisFN
         use sort_mod
         use helem, only: helement_t_size
         use util_mod, only: neci_icopy
         IMPLICIT NONE
         Type(BasisFN) G1(*)
         INTEGER nEl,nI(nEl),i_P,nBasis,nMsh
         INTEGER nMax,nTay(2),iMaxExcit
         complex(dp) fck(*)
         HElement_t UMat(*)
         real(dp) Beta, ALat(3),RhoEps,ECore
         HElement_t  ExcitInfo(0:iMaxExcit,0:2)
!.. New lists are generated here
!.. This will contain all the info needed to work out the value of the
!.. star
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ

         INTEGER nJ(nEl),iEx(2,2),nK(nEl), iDummy
         HElement_t StarMat(5,5),rh,rhii
         integer i,iExcit,iExc

!Needed for diagonalizer
         real(dp) WLIST(5),WORK(3*5)         
         HElement_t NWORK(4*5)
         INTEGER(int32) INFO

         StarMat=(0.0_dp)
         iEx(1,1)=2
!.. Get the orbitals which are excited in going from I to J
!.. IEX(1,*) are in I, and IEX(2,*) are in J
         CALL GetExcitation(nI,nJ,nEl,iEx,iDummy)
         IF(iEx(1,1).GT.0.AND.iEx(1,2).GT.0) THEn
!  We've got a double excitation
            StarMat(1,1)=ExcitInfo(iExcit,0)
            iExc=1
!  Now generate all possible singles between nI and nJ, and put them into StarMat
!(i,a)
            CALL NECI_ICOPY(nEl,nI,1,nK,1)
      lp0:  DO i=1,nEl
               IF(nK(i).EQ.iEx(1,1)) THEN
                  nK(i)=iEx(2,1)
                  call sort (nK)
                  exit lp0
               endif
            end do lp0
            CALL CalcRho2(nJ,nK,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,1,ECore)
!            call writedet(6,nk,nel,.false.)
!            write(6,*) rh
            IF(abs(rh).ge.rhoeps) then
               iExc=iExc+1
               StarMat(1,iExc)=rh/rhii
               CALL CalcRho2(nK,nK,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
               StarMat(iExc,iExc)=rh/rhii
            ENDIF
!(i,b)
            CALL NECI_ICOPY(nEl,nI,1,nK,1)
        lp1:DO i=1,nEl
               IF(nK(i).EQ.iEx(1,1)) THEN
                  nK(i)=iEx(2,2)
                  call sort (nK)
                  exit lp1
               endif
            end do lp1
            CALL CalcRho2(nJ,nK,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,1,ECore)
!            call writedet(6,nk,nel,.false.)
!            write(6,*) rh
            IF(abs(rh).ge.rhoeps) then
               iExc=iExc+1
               StarMat(1,iExc)=rh/rhii
               CALL CalcRho2(nK,nK,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
               StarMat(iExc,iExc)=rh/rhii
            ENDIF
!(j,a)
            CALL NECI_ICOPY(nEl,nI,1,nK,1)
       lp2: DO i=1,nEl
               IF(nK(i).EQ.iEx(1,2)) THEN
                  nK(i)=iEx(2,1)
                  call sort (nK)
                  exit lp2
               endif
            end do lp2
            CALL CalcRho2(nJ,nK,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,1,ECore)
!            call writedet(6,nk,nel,.false.)
!            write(6,*) rh
            IF(abs(rh).ge.rhoeps) then
               iExc=iExc+1
               StarMat(1,iExc)=rh/rhii
               CALL CalcRho2(nK,nK,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
               StarMat(iExc,iExc)=rh/rhii
            ENDIF
!(j,b)
            CALL NECI_ICOPY(nEl,nI,1,nK,1)
   lp3:     DO i=1,nEl
               IF(nK(i).EQ.iEx(1,2)) THEN
                  nK(i)=iEx(2,2)
                  call sort (nK)
                  exit lp3
               endif
            end do lp3
            CALL CalcRho2(nJ,nK,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,1,ECore)
!            call writedet(6,nk,nel,.false.)
!            write(6,*) rh
            IF(abs(rh).ge.rhoeps) then
               iExc=iExc+1
               StarMat(1,iExc)=rh/rhii
               CALL CalcRho2(nK,nK,Beta,i_P,nEl,G1,nBasis,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
               StarMat(iExc,iExc)=rh/rhii
            ENDIF

! Now diagonalize.
            IF(HElement_t_size.EQ.1) THEN
!real case
               CALL DSYEV('V','U',iExc,StarMat,5,WLIST,WORK,3*iExc,INFO)
               IF(INFO.NE.0) THEN
                  WRITE(6,*) 'DYSEV error: ',INFO
                  STOP
               ENDIF
            ELSE
!.. The complex case
               CALL ZHEEV('V','U',iExc,StarMat,5,WLIST,NWORK,4*iExc,WORK,INFO)
               IF(INFO.NE.0) THEN
                  WRITE(6,*) 'ZHEEV error: ',INFO
                  STOP
               ENDIF
            ENDIF
!.. StarMat now contains the eigenvectors, and WLIST the eigenvalues         
            do i=iExc,1,-1
!               write(6,"(I,6F)") i,StarMat(1:iExc,i),wlist(i)
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ
!(rho__^P)_{jj}=sum_l v*_{lj} (r_j)^P v_{lj}
!  r_j -> rho_jj
!  |v_{lj}|^2 is the amount of the original j there is in each of the new eigenvectors l.
               ExcitInfo(iExcit+i-1,0)=WList(i)
#ifdef __CMPLX
               rh=sqrt((StarMat(1,i)*conjg(StarMat(1,I))))
#else
               rh=sqrt((StarMat(1,i)*(StarMat(1,I))))
#endif
               ExcitInfo(iExcit+i-1,1)=ExcitInfo(iExcit,1)*rh
               ExcitInfo(iExcit+i-1,2)=ExcitInfo(iExcit,2)*rh
            enddo
            iExcit=iExcit+iExc-1
         ENDIF
      END

      !Routine which takes the excitation form for the determinant, and calculates the full path.
      SUBROUTINE GETFULLPATH(nI,nEl,noexcits,Orbchange,nJ)
        use sort_mod
        IMPLICIT NONE
        INTEGER :: nEl,nI(nEl),nJ(nEl),noexcits,Orbchange(noexcits*noexcits)
        INTEGER :: I,J
        
        DO I=1,nEl
            nJ(I)=nI(I)
        ENDDO
        DO I=1,noexcits
            DO J=1,nEl
                IF(Orbchange(I).eq.nI(J)) THEN
                    nJ(J)=Orbchange(I+noexcits)
                    EXIT
                ENDIF
            ENDDO
        ENDDO
        call sort (nJ)
      END
      
     !Routine which takes a root determinant (nI), and a double excitation (nJ), and calculates the orbitals which have been excited.
     !This information is put into Orbchange(4), with the first two values being the excited from orbitals (ij), and the second two being the excited to orbitals (ab).
     SUBROUTINE GETEXCITSCHANGE(nI,nJ,nEl,Orbchange)
        IMPLICIT NONE
        INTEGER :: nEl,nI(nEl),nJ(nEl),Orbchange(4),q,I,J
        LOGICAL :: FOUND
        LOGICAL :: ROOT(nEl),EXCIT(nEl)
        ROOT(:)=.TRUE.
        EXCIT(:)=.TRUE.
        
        q=1
        DO I=1,nEl
            FOUND=.FALSE.
            DO J=1,nEl
                IF(nI(I).eq.nJ(J)) THEN
                    FOUND=.TRUE.
                    ROOT(I)=.FALSE.
                    EXCIT(J)=.FALSE.
                    EXIT
                ENDIF
            ENDDO
            IF(.NOT.FOUND) THEN
                IF(q.eq.1) THEN
                    Orbchange(1)=nI(I)
                    q=2
                ELSEIF(q.eq.2) THEN
                    Orbchange(2)=nI(I)
                    q=3
                ELSEIF(q.eq.3) THEN
                    STOP 'ERROR IN GETEXCITSCHANGE'
                ENDIF
            ENDIF
        ENDDO
        
        q=1
        DO I=1,nEl
            IF(EXCIT(I)) THEN
                IF(q.eq.1) THEN
                    Orbchange(3)=nJ(I)
                    q=2
                ELSEIF(q.eq.2) THEN
                    Orbchange(4)=nJ(I)
                    RETURN
                ELSE
                    STOP 'ERROR IN GETEXCITSCHANGE'
                ENDIF
            ENDIF
        ENDDO

        STOP 'ERROR IN GETEXCITSCHANGE'

        END
        
