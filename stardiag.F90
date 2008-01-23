    MODULE STARDIAGMOD
        USE HElement
        IMPLICIT NONE
      
         TYPE(HElement), ALLOCATABLE :: ExcitInfo(:,:)
         TYPE(HElement), ALLOCATABLE :: ExcitInfo2(:,:)
         TYPE(HElement), ALLOCATABLE :: QUADRHOS(:,:)
         TYPE(HElement), ALLOCATABLE :: temprhos(:,:)
         TYPE(HElement), ALLOCATABLE :: Offrho(:)
         INTEGER, ALLOCATABLE :: EXCITSTORE(:,:)
         INTEGER, ALLOCATABLE :: PRODPOSITIONS(:,:)
         REAL*8, ALLOCATABLE :: ONDIAGPRODRHO(:)
         REAL*8, ALLOCATABLE :: OFFDIAGPRODRHO(:,:)

         contains
!  A function to generate all possible excitations of a determinant and link them together in a star
!    The rho_jj, rho_ij and H_ij values are stored for each connected determinant.
!   Based on a combined FMCPR3STAR and FMCPR3STAR2, this instead generates excitations on the fly.
   FUNCTION fMCPR3StarNewExcit(nI,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,nTay, &
               RhoEps, L, LT,nWHTay, iLogging, tSym, ECore,dBeta,dLWdB,MP2E)
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         INCLUDE 'vmc.inc'
         INCLUDE 'uhfdet.inc'
         Type(BasisFN) G1(*)
         INTEGER nI(nEl),nEl,i_P,nBasisMax(*),Brr(nBasis),nBasis,nMsh
         INTEGER nMax,nTay(2),L,LT,nWHTay,iLogging
         COMPLEX*16 fck(*)
         TYPE(HElement) UMat(*)
         REAL*8 Beta, ALat(3),RhoEps,ECore,dBeta
         TYPE(HDElement) dLWdB
         TYPE(HDElement) fMCPR3StarNewExcit
         TYPE(HElement) rhiiadd
         TYPE(HElement) HIJS(0:2)
!         REAL*8 LARGERHOJJ(10000)
         INTEGER iPath(nEl,0:2),uniqprod,nouniqprod
         LOGICAL tSym
!.. New lists are generated here
!.. This will contain all the info needed to work out the value of the
!.. star
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ

         INTEGER exFlag,totvert
         INTEGER, allocatable :: nExcit(:)
         INTEGER nExcitMemLen,nStore(6),prodorbs(8)
         INTEGER nJ(nEl),iExcit,iMaxExcit,excitcount,Prodnum
         INTEGER iErr,nK(nEl),nL(nEl)
         INTEGER nRoots,i,j,rhoelem
         TYPE(HElement) rh,rhii,EHFDiff
         TYPE(HDElement) MP2E         
         LOGICAL tStarSingles,tCountExcits
         INTEGER nIExcitFormat(nEl)

!         LARGERHOJJ(:)=0.D0
         IF(tStoreAsExcitations) THEN
            nIExcitFormat(1)=-1
            nIExcitFormat(2)=0
         ELSE
            CALL ICOPY(NEL,nI,1,nIExcitFormat,1)
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
         CALL GenSymExcitIt2(nI,nEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
         Allocate(nExcit(nExcitMemLen))
!Second call calculates size of arrays needed to store all symmetry allowed excitations - further calls will generate excitation on-the-fly(shown by the false in arg(6)
         nExcit(1)=0
         CALL GenSymExcitIt2(nI,nEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
!.. iMaxExcit now contains the number of excitations.
         !TCountExcits will run through all excitations possible, determine if they are connected, and then only store these.
         !Will be twice as expensive, as needs to run through all excitations twice - however, will only store memory needed.
         IF(tCountExcits) THEN
            Write(6,"A,I10,A") "Counting excitations - Running through all ",iMaxExcit," excitations to determine number connected"
            excitcount=0
            CALL CalcRho2(nI,nI,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rhii,nTay,0,ECore)
            
       lp2: do while(.true.)
                CALL GenSymExcitIt2(nI,nEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,exFlag)
                IF(nJ(1).eq.0) exit lp2
                CALL CalcRho2(nIExcitFormat,nJ,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,iExcit,ECore)
                IF(rh.agt.RhoEps) excitcount=excitcount+1
            enddo lp2

            !Set number of excitations to number of connected determinants, and reset generator
            !This routine doesn't seem to work...would be good if didn't need to reinitialise excitation generator
!            CALL ResetExit2(nI,nEl,G1,nBasis,nBasisMax,nExcit,0)
            Deallocate(nExcit)
            nStore(1)=0
            CALL GenSymExcitIt2(nI,nEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
            Allocate(nExcit(nExcitMemLen))
            nExcit(1)=0
            CALL GenSymExcitIt2(nI,nEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)

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
         CALL MemAlloc(iErr,ExcitInfo,(iMaxExcit+1)*3*HElementSize,"ExcitInfo")
         
         !If we are keeping a list of excitations in the star, then allocate EXCITSTORE to hold the excitations, in the form of o o v v orbitals
         IF(TCALCREALPROD.or.TSUMPROD) THEN
             ALLOCATE(EXCITSTORE(4,iMaxExcit),stat=iErr)
             CALL MemAlloc(iErr,EXCITSTORE,(iMaxExcit*2),"EXCITSTORE")
         ENDIF

!.. This will contain all the info needed to work out the value of the
!.. star
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ
         i=0
         ExcitInfo(i,0)=1.D0
         ExcitInfo(i,1)=1.D0
         ExcitInfo(i,2)=GetHElement2(nI,nI,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,00,ECore)
         if(BTEST(nWhTay,5)) then
! We use the Zeroth order N-particle Hartree-Fock hamiltonian (as MP theory), but shifted by E_HF-E0.

!nMax has Arr hidden in it
            IF(nTay(2).eq.5) THEN
               call GetH0ElementDCCorr(nI,nI,nEl,nBasisMax,G1,nBasis,Brr,NMSH,FCK,nMax,ALAT,UMat,ECore,rhii)
            ELSE
               call GetH0Element(nI,nEl,nMax,nBasis,rhii)
            ENDIF
            EHFDiff=ExcitInfo(i,2)-rhii
         endif
         CALL CalcRho2(nI,nI,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rhii,nTay,0,ECore)
!         write(75,*) rhii
! Setup MP info
         CALL iCopy(nEl,nI,1,iPath(1,0),1)
         CALL iCopy(nEl,nI,1,iPath(1,2),1)
         HIJS(0)=ExcitInfo(0,2)
         HIJS(2)=ExcitInfo(0,2)

    lp:  do while(.true.)
            CALL GenSymExcitIt2(nI,nEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,exFlag)
            IF(nJ(1).eq.0) exit lp

!           Calculate rhoij element
            CALL CalcRho2(nIExcitFormat,nJ,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,iExcit,ECore)
            
            if(rh .agt. RhoEps) then
               i=i+1
!   Divide all elements though by rhoii
               ExcitInfo(i,1)=rh/rhii
               
               IF(TCALCREALPROD.or.TSUMPROD) THEN
!   Stores all excitations as (occ,occ,vir,vir)
                   CALL GETEXCITSCHANGE(nI,nJ,nEl,EXCITSTORE(:,i))
               ENDIF
               if(btest(nwhtay,5)) then
                  call GetH0Element(nJ,nEl,nMax,nBasis,rh)
                  rh=rh+EHFDiff
                  rh=rh*HElement(-Beta/I_P)
                  rh=exp(rh)
               else
                  !RHO_JJ elements
                   CALL CalcRho2(nJ,nJ,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
                  
!                  do j=1,10000
!                    IF((rh.agt.LARGERHOJJ(J)).or.(LARGERHOJJ(J).eq.0.D0)) THEN
!                        LARGERHOJJ(J)=rh%v
!                        GOTO 765
!                    ENDIF
!                  ENDDO
!765               CONTINUE
               endif
               ExcitInfo(i,0)=rh/rhii
               ExcitInfo(i,2)=GetHElement2(nIExcitFormat,nJ,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,iExcit,ECore)
!               write(75,*) rh,rh/rhii
!Now do MP2
               Hijs(1)=ExcitInfo(i,2)
               IF(tMPTheory) THEN
                  call iCopy(nEl,nJ,1,iPath(1,1),1)
!nMax has Arr hidden in it
                  Call AddMP2E(Hijs,nMax,nBasis,iPath,nEl,BTEST(iLogging,0),MP2E)
               ENDIF
               IF(tStarSingles) Call StarAddSingles(nI,nJ,ExcitInfo,i,iMaxExcit,rhii,rhoeps,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,nTay,ECore)
            endif
         enddo lp
!Tell MCPATHS how many excitations there were and how many we are keeping
         L=i
         LT=iMaxExcit
         iExcit=i
         Deallocate(nExcit)

         !FIND REAL PRODUCT EXCITATIONS i.e. neither of the two 'from' or 'two' orbitals are the same in the two excitations.
         IF(TCALCREALPROD) THEN
             
             !TSUMPROD tries to find a way to resum the 
             IF(TSUMPROD) THEN
                rhiiadd=0.D0
                rhoelem=iExcit
                ALLOCATE(Offrho(rhoelem),stat=iErr)
                CALL MemAlloc(iErr,Offrho,rhoelem*HElementSize,"OffRho")
!   Fill Offrho with the off-diagonal rho elements to pass to COUNTPRODEXCITS
                DO I=1,iExcit
                    Offrho(I)=ExcitInfo(I,1)
                ENDDO
                CALL COUNTPRODEXCITS(iMaxExcit,prodnum,.FALSE.,iExcit,rhoelem,rhiiadd,uniqprod)
                WRITE(6,*) prodnum, "product excitations found, and summed in..."
                Call MemDealloc(Offrho)
                DEALLOCATE(Offrho)
                
                IF(prodnum.gt.0) THEN
                    WRITE(6,*) "Sum of product diagonal rho elements is ", rhiiadd
                    WRITE(6,*) "Number of unique products is ", uniqprod
!                    rhiiadd=(rhiiadd%v+DREAL(uniqprod))/(2*uniqprod)
!   rhiiadd now in divided by the number of product excitations we are resumming in
                    rhiiadd=rhiiadd/prodnum
!                    ExcitInfo(0,0)=rhiiadd
!                    ExcitInfo(0,1)=rhiiadd
!   The resummed values for the diagonal elements of the product excitations are added to the root
                    ExcitInfo(0,0)=ExcitInfo(0,0)+rhiiadd
                    ExcitInfo(0,1)=ExcitInfo(0,1)+rhiiadd
                    WRITE(6,*) "New root is now ",ExcitInfo(0,0) 
                ELSE
                    IF(rhiiadd.agt.(0.D0)) STOP 'rhiiadd should be zero as no products'
                ENDIF

             !TQUASIEXCIT approximates the product excitations to a separate excitation attached to the root - maintains form of rho matrix
             ELSEIF(TQUASIEXCIT) THEN
                !Space is allocated for all possible product excitations (N choose 2), and elements(1) contain diagonal elements, (2)=offdiagonal (3)=Helement
                ALLOCATE(QUADRHOS(3,iExcit*(iExcit-1)/2),stat=ierr)
                CALL MemAlloc(ierr,QUADRHOS,iExcit*(iExcit-1)*HElementSize,"QUADRHOSALL")
                CALL COUNTUNIQPRODEXCITS(iMaxExcit,prodnum,iExcit,nouniqprod,prodnum)
                WRITE(6,*) "Theoretical max products is ", iExcit*(iExcit-1)/2
                WRITE(6,*) "Total number of products is ", prodnum
                WRITE(6,*) "Number of unique product excitations is ", nouniqprod
               
                CALL FLUSH(6)
                Totvert=nouniqprod+iExcit
                !Create new ExcitInfo, including the product excitations
                Allocate(ExcitInfo2(0:totvert,0:2),stat=iErr)
                CALL MemAlloc(iErr,ExcitInfo2,(nouniqprod+1)*3*HElementSize,"ExcitInfo2")
                ExcitInfo2(0:iExcit,0)=ExcitInfo(0:iExcit,0)
                ExcitInfo2(0:iExcit,1)=ExcitInfo(0:iExcit,1)
                ExcitInfo2(0:iExcit,2)=ExcitInfo(0:iExcit,2)
                ExcitInfo2(iExcit+1:iExcit+nouniqprod,0)=QUADRHOS(1,1:nouniqprod)
                ExcitInfo2(iExcit+1:iExcit+nouniqprod,1)=QUADRHOS(2,1:nouniqprod)
                ExcitInfo2(iExcit+1:iExcit+nouniqprod,2)=QUADRHOS(3,1:nouniqprod)

                !Deallocate information about the excitations
                CALL MemDealloc(QUADRHOS)
                DEALLOCATE(QUADRHOS)
                IF(ALLOCATED(ExcitInfo)) THEN
                    CALL MemDealloc(ExcitInfo)
                    DEALLOCATE(ExcitInfo)
                ENDIF
           
            !This calculates all product excitations
            ELSE
                !Call twice - once to calculate number of products, so can allocate memory, then store products
                CALL COUNTPRODEXCITS(iMaxExcit,prodnum,.FALSE.,iExcit,rhoelem,rhiiadd,uniqprod)
                WRITE(6,*) prodnum, "product excitations found - allocating memory..."
                ALLOCATE(prodpositions(2,prodnum),stat=iErr)
                CALL MemAlloc(iErr,prodpositions,prodnum,"EXCITSTORE")
                !Prodpositions stores the position of the two excitations in EXCITSTORE which give rise to a real product excitation
                CALL COUNTPRODEXCITS(iMaxExcit,prodnum,.TRUE.,iExcit,rhoelem,rhiiadd,uniqprod)
             
                !Allocate memory for on diagonal rho elements of product excitations, and 2 x offdiagonal elements
                ALLOCATE(ONDIAGPRODRHO(prodnum),stat=ierr)
                CALL MemAlloc(ierr,ONDIAGPRODRHO,prodnum,"ONDIAGPRODRHO")
                ALLOCATE(OFFDIAGPRODRHO(2,prodnum),stat=ierr)
                CALL MemAlloc(ierr,OFFDIAGPRODRHO,2*prodnum,"OFFDIAGPRODRHO")

                DO I=1,prodnum
                    !Approximate on-diag elements as product of consituent excits (exact for FPLD)(remember, they are divided by rhoii^2)
                    ONDIAGPRODRHO(I)=DREAL(ExcitInfo(prodpositions(1,I),0)*ExcitInfo(prodpositions(2,I),0))

                    IF(TCALCRHOPROD) THEN
                        IF(I.eq.1) WRITE(6,*) "Calculating off-diagonal rho elements for product excitations exactly"
                        !This calculates the exact rho elements between the product excitation, and the constituent determinants
                        prodorbs(1:2)=EXCITSTORE(1:2,prodpositions(1,I))
                        prodorbs(3:4)=EXCITSTORE(1:2,prodpositions(2,I))
                        prodorbs(5:6)=EXCITSTORE(3:4,prodpositions(1,I))
                        prodorbs(7:8)=EXCITSTORE(3:4,prodpositions(2,I))
                        CALL GetFullPath(nI,nEl,4,prodorbs,nJ)
                        CALL GetFullPath(nI,nEl,2,EXCITSTORE(:,prodpositions(1,I)),nK)
                        CALL GetFullPath(nI,nEl,2,EXCITSTORE(:,prodpositions(2,I)),nL)
                        !nJ is now the product IPATH, and nK and nL the two consituent determinants
                        CALL CalcRho2(nK,nJ,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,2,ECore)
                        OFFDIAGPRODRHO(1,I)=DREAL(rh/rhii)
                        CALL CalcRho2(nL,nJ,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,2,ECore)
                        OFFDIAGPRODRHO(2,I)=DREAL(rh/rhii)
                    ELSE
                        IF(I.eq.1) WRITE(6,*) "Approximating off-diagonal elements for product excitations"
                        !Approximate the off-diag rho elements between product excitation and constituent determinant as being equal to the rho element between the root and other constituent determinant.
                        OFFDIAGPRODRHO(1,I)=ExcitInfo(prodpositions(2,I),1)%v
                        OFFDIAGPRODRHO(2,I)=ExcitInfo(prodpositions(1,I),1)%v
                    ENDIF
                ENDDO
            ENDIF
        ENDIF
         
!         DO j=1,10000
!            WRITE(55,*) LARGERHOJJ(J)
!         ENDDO
!.. we now have a list length NLCUR of dets in the star.
!.. Call a routine to generate the value of the star
         WRITE(6,*) iExcit," excited determinants in star"
         IF(.NOT.BTEST(NWHTAY,0)) THEN
            WRITE(6,*) "Beginning Complete Star Diagonalization"
            IF(STARPROD) THEN
                IF(TCALCREALPROD) THEN
                    IF(TSUMPROD) THEN
                        WRITE(6,*) "Product Sum Star Diagonalisation"
                        CALL StarDiag(0,nEl,iExcit+1,ExcitInfo,iMaxExcit+1,i_p,fMCPR3StarNewExcit,dBeta,dLWdB)
                    ELSEIF(TQUASIEXCIT) THEN
                        WRITE(6,*) "Quadratic quasiexcitations Star Diagonalisation"
                        CALL StarDiag(0,nEl,totvert+1,ExcitInfo2,totvert+1,i_p,fMCPR3StarNewExcit,dBeta,dLWdB)
                    ELSE
                        WRITE(6,*) "Real Product Star Diagonalisation"
                        CALL StarDiagRealProd(nEl,iExcit+1,ExcitInfo,iMaxExcit+1,i_P,fMCPR3StarNewExcit,dBeta,dLWdB,prodnum,EXCITSTORE,prodpositions,ONDIAGPRODRHO,OFFDIAGPRODRHO)
                    ENDIF
                ELSE
                    WRITE(6,*) "Complete Product Star Diagonalisation -  beware - large scaling!"
                    CALL StarDiagSC(0,nEl,iExcit+1,ExcitInfo,iMaxExcit+1,i_P,fMCPR3StarNewExcit,dBeta,dLWdB)
                ENDIF
            ELSE
                CALL StarDiag(0,nEl,iExcit+1,ExcitInfo,iMaxExcit+1,i_P,fMCPR3StarNewExcit,dBeta,dLWdB)
            ENDIF
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
            CALL StarDiag2(0,nEl,iExcit+1,ExcitInfo,iMaxExcit+1,Beta,i_P,fMCPR3StarNewExcit,dBeta,dLWdB,nRoots,iLogging)
         ENDIF
         
         IF(ALLOCATED(ExcitInfo)) THEN
            call MemDealloc(ExcitInfo)
            Deallocate(ExcitInfo)
         ENDIF
         IF(ALLOCATED(ExcitInfo2)) THEN
             call MemDealloc(ExcitInfo2)
             Deallocate(ExcitInfo2)
         ENDIF
         IF(ALLOCATED(EXCITSTORE)) THEN
            call MemDealloc(EXCITSTORE)
            DEALLOCATE(EXCITSTORE)
         ENDIF
         IF(ALLOCATED(PRODPOSITIONS)) THEN
            call MemDealloc(PRODPOSITIONS)
            DEALLOCATE(PRODPOSITIONS)
         ENDIF
         IF(ALLOCATED(ONDIAGPRODRHO)) THEN
            call MemDealloc(ONDIAGPRODRHO)
            DEALLOCATE(ONDIAGPRODRHO)
            call MemDealloc(OFFDIAGPRODRHO)
            DEALLOCATE(OFFDIAGPRODRHO)
         ENDIF
            
      END

      ! This routine finds all "unique" product excitations of the star
      SUBROUTINE COUNTUNIQPRODEXCITS(iMaxExcit,prodnum,iExcit,nouniqprod,countprods)
        IMPLICIT NONE
        include 'vmc.inc'
        INTEGER iMaxExcit,prodnum
        INTEGER iExcit,I,J
        INTEGER countprods,nouniqprod,tempprodfrom(4),tempprodto(4),ierr,k
        LOGICAL UNIQPROD
        INTEGER, ALLOCATABLE :: UNIQPRODS(:,:)
!        INTEGER, POINTER :: tempprods(:,:)

        !Space allocated to store the excitations of the unique products, so we can search through them to ensure they are not being duplicated
        ALLOCATE(UNIQPRODS(8,(iExcit*(iExcit-1)/2)),stat=ierr)
        CALL MemAlloc(ierr,UNIQPRODS,iExcit*(iExcit-1)*2,"UNIQPRODSALL")
        
        countprods=0
        nouniqprod=0
        DO I=1,iExcit
            DO J=(I+1),iExcit
                UNIQPROD=.TRUE.
                ! Again, the criteria for a product excitation are that none of the "ij" or "ab" orbitals match in the constituent excitations
                IF((EXCITSTORE(1,I).NE.EXCITSTORE(1,J)).AND.(EXCITSTORE(1,I).NE.EXCITSTORE(2,J))                   &
                .AND.(EXCITSTORE(2,I).NE.EXCITSTORE(2,J)).AND.(EXCITSTORE(2,I).NE.EXCITSTORE(1,J))) THEN
                    IF((EXCITSTORE(3,I).NE.EXCITSTORE(3,J)).AND.(EXCITSTORE(3,I).NE.EXCITSTORE(4,J))               &
                    .AND.(EXCITSTORE(4,I).NE.EXCITSTORE(3,J)).AND.(EXCITSTORE(4,I).NE.EXCITSTORE(4,J))) THEN
                
                        countprods=countprods+1
                        tempprodfrom(1:2)=EXCITSTORE(1:2,I)
                        tempprodfrom(3:4)=EXCITSTORE(1:2,J)
                        tempprodto(1:2)=EXCITSTORE(3:4,I)
                        tempprodto(3:4)=EXCITSTORE(3:4,J)
                        CALL SORTI(4,tempprodfrom(:))
                        CALL SORTI(4,tempprodto(:))
                        !Cycle through other stored products to see if unique
                        DO K=1,nouniqprod
                            !Test if unique
                            IF((tempprodfrom(1).eq.UNIQPRODS(1,K)).AND.(tempprodfrom(2).eq.UNIQPRODS(2,K)).AND.     &
                                (tempprodfrom(3).eq.UNIQPRODS(3,K)).AND.(tempprodfrom(4).eq.UNIQPRODS(4,K))) THEN
                                IF((tempprodto(1).eq.UNIQPRODS(5,K)).AND.(tempprodto(2).eq.UNIQPRODS(6,K)).AND.     &
                                    (tempprodto(3).eq.UNIQPRODS(7,K)).AND.(tempprodto(4).eq.UNIQPRODS(8,K))) THEN
                                    UNIQPROD=.FALSE.
                                    !Sum in off-diagonal contribution from another way to get same product
                                    QUADRHOS(2,K)=QUADRHOS(2,K)+(ExcitInfo(I,1)*ExcitInfo(J,1))/(ExcitInfo(I,1)+ExcitInfo(J,1))
                                    !If the excitation is found, then exit loop
                                    EXIT
                                ENDIF
                            ENDIF
                        ENDDO
                        IF(UNIQPROD) THEN
!                            WRITE(67,*) "Unique Product Excitation: ", I,J
!                            WRITE(67,*) "Diagonal element is: ", ExcitInfo(I,0)*ExcitInfo(J,0)
!                           If the product is unique, then count it, and store its information in QUADRHOS
                            nouniqprod=nouniqprod+1
                            UNIQPRODS(1:4,nouniqprod)=tempprodfrom(:)
                            UNIQPRODS(5:8,nouniqprod)=tempprodto(:)
                            QUADRHOS(1,nouniqprod)=ExcitInfo(I,0)*ExcitInfo(J,0)
                            QUADRHOS(2,nouniqprod)=(ExcitInfo(I,1)*ExcitInfo(J,1))/(ExcitInfo(I,1)+ExcitInfo(J,1))
                            !Let hamiltonian elements to these excitations = 0
                            QUADRHOS(3,nouniqprod)=0.D0
                        ENDIF
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
        CALL MemDealloc(UNIQPRODS)
        DEALLOCATE(UNIQPRODS)
        
        prodnum=countprods
      
      END
                        

      SUBROUTINE COUNTPRODEXCITS(iMaxExcit,prodnum,setup,iExcit,rhoelem,rhiiadd,uniqprod)
        IMPLICIT NONE
        include 'vmc.inc'
        INTEGER iMaxExcit,prodnum!,prodpositions(2,length),EXCITSTORE(4,iMaxExcit),
        INTEGER iExcit,I,J
!        Type(HElement) :: ExcitInfo(0:iMaxExcit,0:2),OffRho(rhoelem),
        Type(Helement) rhiiadd
        INTEGER countprods,rhoelem,uniqprod!,length
        LOGICAL setup,FIRST

        countprods=0
        uniqprod=0
        rhiiadd=0.D0
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
                        prodpositions(1,countprods)=I
                        prodpositions(2,countprods)=J
                    ENDIF
                    IF(TSUMPROD) THEN
                        !Add the diagonal rho elements to the quadruple excitation back into the offdiagonal rho elements of the original matrix
                        ExcitInfo(I,1)=ExcitInfo(I,1)+Offrho(J)
                        ExcitInfo(J,1)=ExcitInfo(J,1)+Offrho(I)
                        !Only sum in the diagonal product rho element once per excitation
                        IF(FIRST) THEN
                            !Resum by adding in the product of the diagonal rhoelements from the constituent excitations
                            rhiiadd=rhiiadd+(ExcitInfo(J,0)*ExcitInfo(I,0))
                            !Count the products we are adding in
                            uniqprod=uniqprod+1
                            FIRST=.FALSE.
                        ENDIF
                    ENDIF
                ENDIF
            ENDIF
            ENDDO
        ENDDO
        IF(.NOT.setup) prodnum=countprods
      END

      END MODULE
        
      
!.. This version creates a star of all 2-vertex terms.
!.. This sets up the excitation generators and the memory - using the old excitation generators
FUNCTION FMCPR3STAR(NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY, &
                     RHOEPS,LSTE,ICE,LIST,L,LT,NWHTAY,ILOGGING,TSYM,ECORE,ILMAX,DBETA,DLWDB)
         USE HElement      
         IMPLICIT NONE
         INTEGER I_V,NEL,I_P,NBASISMAX(*),G1(*),NBASIS,BRR(*),NMSH,NMAX
         INTEGER NTAY,NWHTAY,ILOGGING,LT
         REAL*8 FCK(*), ALAT(*),UMAT(*),ECORE
         INTEGER NI(NEL),ILMAX
!.. These are old and no longer used
!.. LSTE is a list of excitations (which we will generate)
!.. ICE is the IC of each excitation (i.e. how much it differs from us (INODE)
         INTEGER LSTE(NEL,0:ILMAX)
         INTEGER ICE(0:ILMAX)
         REAL*8 LIST(0:ILMAX,0:2)
!.. New lists are generated here
         INTEGER NLSTE(NEL,0:*)
         REAL*8  NLIST(0:*)
         INTEGER NICE(0:*)
         TYPE(HDElement) FMCPR3Star
         POINTER (IP_NLIST,NLIST)
         POINTER (IP_NLSTE,NLSTE)
         POINTER (IP_NICE,NICE)

!.. This will contain all the info needed to work out the value of the
!.. star
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ

         REAL*8 BETA,RHOEPS
         LOGICAL TSYM
         REAL*8 DBETA,DLWDB
         INTEGER NLENLIST,NLCUR,I,J,L

         REAL*8 RH,RHII

         INTEGER NORDER,NMIN
         TYPE(HDElement) FMCPR3STAR2
         include 'irat.inc'
  
         IF(HElementSize.NE.1) STOP 'FMCPR3STAR cannot work with complex orbitals.' 
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
         CALL MEMORY(IP_NLIST,(NLENLIST+1)*3,'NLIST')
         CALL MEMORY(IP_NLSTE,((NLENLIST+1)*NEL)/IRAT+1,'NLSTE')
         CALL MEMORY(IP_NICE,(NLENLIST+1)/IRAT+1,'NICE')

         CALL ICOPY(NEL,NI,1,NLSTE(1,0),1)
         NICE(0)=0
!.. Now generate the excitations
         CALL GENEXCIT(NI,NORDER,NBASIS,NEL,NLSTE(1,1),NICE(1),NLENLIST,NMIN,G1,TSYM,NBASISMAX,.FALSE.)

         FMCPR3STAR=FMCPR3STAR2(NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY, &
           RHOEPS,NLSTE,NICE,NLIST,L,LT,NWHTAY,ILOGGING,TSYM,ECORE,NLENLIST,DBETA,DLWDB)

         CALL FREEM(IP_NLIST)
         CALL FREEM(IP_NLSTE)
         CALL FREEM(IP_NICE)
         RETURN
      END

      
!.. This version creates a star of all 2-vertex terms.
!..   This is the heart of the function, called once the excitations are found.
      FUNCTION FMCPR3STAR2(NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,NTAY, &
           RHOEPS,LSTE,ICE,LIST,L,LT,NWHTAY,ILOGGING,TSYM,ECORE,ILMAX,DBETA,DLWDB)
         USE HElement     
         IMPLICIT NONE
         TYPE(HDElement) FMCPR3Star2
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(*)
         INTEGER I_V,NEL,I_P,NBASISMAX(*),NBASIS,BRR(*),NMSH,NMAX
         INTEGER NTAY,NWHTAY,ILOGGING,LT
         REAL*8 ALAT(*),ECORE
         TYPE(HElement) UMat(*)
         COMPLEX*16 FCK(*)
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

         TYPE(HElement) LIST(0:ILMAX,0:2)
         REAL*8 BETA,RHOEPS
         LOGICAL TSYM
         REAL*8 DBETA,DLWDB
         INTEGER NLIST,NLCUR,I,J,L

         TYPE(HElement) RHII,RH
         NLIST=ILMAX

         NLCUR=0
         DO I=0,NLIST
            IF(ICE(I).EQ.0.AND.I.GT.0) LSTE(1,I)=0
            IF(LSTE(1,I).NE.0) THEN
               CALL CALCRHO2(NI,LSTE(1,I),BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT, &
                  RH,NTAY,ICE(I),ECORE)
               IF(RH .AGT. RHOEPS) THEN
                  IF(NLCUR.NE.I) CALL ICOPY(NEL,LSTE(1,I),1,LSTE(1,NLCUR),1)
                  ICE(NLCUR)=ICE(I)
                  IF(NLCUR.EQ.0) RHII=RH
                  LIST(NLCUR,1)=RH/RHII
                  CALL CALCRHO2(LSTE(1,I),LSTE(1,I),BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT, &
                    RH,NTAY,0,ECORE)
                  LIST(NLCUR,0)=RH/RHII
                  LIST(NLCUR,2)=GETHELEMENT2(NI,LSTE(1,I),NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,NMAX,ALAT,UMAT,-1,ECORE)
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
            CALL STARDIAG(LSTE,NEL,NLCUR,LIST,ILMAX+1,I_P,FMCPR3STAR2,DBETA,DLWDB)
         ELSE
            WRITE(6,*) "Beginning Polynomial Star Diagonalization"
            NROOTS=NLCUR
            IF(BTEST(NWHTAY,1)) NROOTS=1
            CALL STARDIAG2(LSTE,NEL,NLCUR,LIST,ILMAX+1,BETA,I_P,FMCPR3STAR2,DBETA,DLWDB,NROOTS,iLogging)
         ENDIF
      END

      SUBROUTINE STARDIAGREALPROD(NEL,NLIST,LIST,ILMAX,I_P,SI,DBETA,DLWDB,PRODNUM,EXCITSTORE,prodpositions,ONDIAGPRODRHO,OFFDIAGPRODRHO)
         !NLIST is now no. original excitations (+ root) - ILMAX is max possible excitations +1
         USE HElement
         IMPLICIT NONE
         INTEGER NEL,I_P
         INTEGER NLIST,ILMAX,PRODNUM,TOTVERT
         REAL*8 LIST(ILMAX,0:2)
         REAL*8 RIJMAT(*),WLIST(*)
         REAL*8 ONDIAGPRODRHO(PRODNUM),OFFDIAGPRODRHO(2,PRODNUM)
         POINTER (IP_RIJMAT,RIJMAT),(IP_WLIST,WLIST),(IP_WORK,WORK)
         INTEGER ISUB,EXCITSTORE(4,ILMAX-1)!contains all excitations (occ,occ,vir,vir)
         INTEGER WORKL,WORK(*),INFO,ierr,PRODPOSITIONS(2,PRODNUM)
         REAL*8 SI,DLWDB,DBETA
         INTEGER I,J
         
         IF(HElementSize.GT.1) STOP "STARDIAGREALPROD cannot function with complex orbitals."

         CALL TISET('STARDIAGRP',ISUB)
         !Is there a need to sort the matrix? If there is, we have problems!
         !CALL SORT3RN(NLIST-1,LIST(2,0),LIST(2,1),LIST(2,2),HElementSize)
         
         TOTVERT=NLIST+PRODNUM
         CALL MEMORY(IP_RIJMAT,TOTVERT*TOTVERT,"RIJMAT")
         CALL AZZERO(RIJMAT,TOTVERT*TOTVERT)
        
         !Fill RIJMAT
         DO I=1,TOTVERT
            IF(I.LE.NLIST) THEN
                RIJMAT(I)=LIST(I,1)
                RIJMAT((I-1)*TOTVERT+I)=LIST(I,0)
            ELSE
                RIJMAT((PRODPOSITIONS(1,I-NLIST)*TOTVERT)+I)=OFFDIAGPRODRHO(1,I-NLIST)
                RIJMAT((PRODPOSITIONS(2,I-NLIST)*TOTVERT)+I)=OFFDIAGPRODRHO(2,I-NLIST)
                RIJMAT((I-1)*(TOTVERT)+I)=ONDIAGPRODRHO(I-NLIST)
            ENDIF
         ENDDO

!.. Debug info         
!         WRITE(68,*) "ROOT RHOII, RHOIJ ", LIST(1,0),LIST(1,1)
!         DO I=2,NLIST
!            WRITE(68,"A,I3,A,2E14.6,4I4") "EXCITATION ",I-1," - RHOII, RHOIJ, IPATH: ",LIST(I,0),LIST(I,1),EXCITSTORE(:,I-1)
!         ENDDO
!         
!         DO I=1,PRODNUM
!            WRITE(68,"A,I3,A,I3,I3") "PRODUCT EXCITATION ",I," FORMED FROM EXCITATIONS ", PRODPOSITIONS(1,I),PRODPOSITIONS(2,I)
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
        
         CALL MEMORY(IP_WLIST,TOTVERT,"WLIST")
         WORKL=3*TOTVERT
         CALL MEMORY(IP_WORK,WORKL,"WORK")
        
!.. Diagonalize
         CALL DSYEV('V','L',TOTVERT,RIJMAT,TOTVERT,WLIST,WORK,WORKL,INFO)
         IF(INFO.NE.0) THEN
            WRITE(6,*) 'DYSEV error: ',INFO
            STOP
         ENDIF
         CALL FREEM(IP_WORK)
         WRITE(6,*)
         WRITE(6,*) "Highest root:",WLIST(TOTVERT)

         SI=0.D0
         DO I=1,TOTVERT
            SI=SI+RIJMAT(((I-1)*TOTVERT)+1)*RIJMAT(((I-1)*TOTVERT)+1)*(WLIST(I)**I_P)
            IF(DBETA.NE.0.D0) THEN
                DO J=1,NLIST
                    !only sum over vertices linked to i?
                    DLWDB=DLWDB+RIJMAT(((I-1)*TOTVERT)+1)*RIJMAT(((I-1)*TOTVERT)+J)*(WLIST(I)**I_P)*LIST(J,2)
                ENDDO
            ENDIF
         ENDDO
         WRITE(6,*) "Final SI= ", SI
         SI=SI-1.D0
         DLWDB=DLWDB-LIST(1,2)
         CALL FREEM(IP_WLIST)
         CALL FREEM(IP_RIJMAT)
         CALL TIHALT("STARDIAGRP",ISUB)

         RETURN
      END SUBROUTINE

         
         
      SUBROUTINE STARDIAGSC(LSTE,NEL,NLIST,LIST,ILMAX,I_P,SI,DBETA,DLWDB)
         USE HElement
         IMPLICIT NONE
         INTEGER NEL,I_P
         INTEGER LSTE(NEL,NLIST),NLIST,ILMAX
         REAL*8 LIST(ILMAX,0:2)
         REAL*8 RIJMAT(*),WLIST(*)
         REAL*8, DIMENSION(:,:), POINTER :: AOFFDB
         REAL*8, DIMENSION(:), POINTER :: AONDB
         POINTER (IP_RIJMAT,RIJMAT),(IP_WLIST,WLIST),(IP_WORK,WORK)
         INTEGER ISUB,IND,TOTVERT
         INTEGER WORKL,WORK(*),INFO,PRODVERT,ierr
         REAL*8 SI,DLWDB,DBETA
         INTEGER I,J
         TYPE(HElement) RR

         IF(HElementSize.GT.1) STOP "STARDIAGSC cannot function with complex orbitals."

         CALL TISET('STARDIAGSC',ISUB)

         CALL SORT3RN(NLIST-1,LIST(2,0),LIST(2,1),LIST(2,2),HElementSize)
         
         PRODVERT=(NLIST-1)*(NLIST-2)/2
         TOTVERT=PRODVERT+NLIST
         
         ALLOCATE(AOFFDB(PRODVERT,NLIST-1),STAT=ierr)
         CALL MemAlloc(ierr,AOFFDB,(NLIST-1)*PRODVERT,'AOFFDB')
         CALL AZZERO(AOFFDB,(NLIST-1)*PRODVERT)
         ALLOCATE(AONDB(PRODVERT),STAT=ierr)
         CALL MemAlloc(ierr,AONDB,PRODVERT,'AONDB')
         CALL AZZERO(AONDB,PRODVERT)

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
            CALL FLUSH(6)
            STOP 'WRONG NUMBER OF ADDED VERTICES'
        ENDIF

!        DO I=1,(NLIST-1)
!            DO J=1,PRODVERT
!                WRITE(68,"E15.6,$") AOFFDB(J,I)
!            ENDDO
!            WRITE(68,*) ""
!        ENDDO
!        WRITE(68,*) "*****************"
        
        CALL MEMORY(IP_RIJMAT,TOTVERT*TOTVERT,"RIJMAT")
        CALL AZZERO(RIJMAT,TOTVERT*TOTVERT)
        
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

        CALL MemDealloc(AONDB)
        Deallocate(AONDB)
        NULLIFY(AONDB)
        CALL MemDealloc(AOFFDB)
        Deallocate(AOFFDB)
        NULLIFY(AOFFDB)

!        DO I=1,TOTVERT
!            DO J=1,TOTVERT
!                WRITE(68,"E14.6,$") RIJMAT((I-1)*TOTVERT+J)
!            ENDDO
!            WRITE(68,*) ""
!        ENDDO

        CALL MEMORY(IP_WLIST,TOTVERT,"WLIST")
        WORKL=3*TOTVERT
        CALL MEMORY(IP_WORK,WORKL,"WORK")
        
!.. Diagonalize
         CALL DSYEV('V','L',TOTVERT,RIJMAT,TOTVERT,WLIST,WORK,WORKL,INFO)
         IF(INFO.NE.0) THEN
            WRITE(6,*) 'DYSEV error: ',INFO
            STOP
         ENDIF
         CALL FREEM(IP_WORK)
         WRITE(6,*)
         WRITE(6,*) "Highest root:",WLIST(TOTVERT)

         SI=0.D0
         DO I=1,TOTVERT
            SI=SI+RIJMAT(((I-1)*TOTVERT)+1)*RIJMAT(((I-1)*TOTVERT)+1)*(WLIST(I)**I_P)
            IF(DBETA.NE.0.D0) THEN
!                OD=DLWDB
                DO J=1,NLIST
                    !Is this a correct formulation for the hamiltonian elements - only sum over vertices linked to i
                    DLWDB=DLWDB+RIJMAT(((I-1)*TOTVERT)+1)*RIJMAT(((I-1)*TOTVERT)+J)*(WLIST(I)**I_P)*LIST(J,2)
                ENDDO
            ENDIF
         ENDDO
         WRITE(6,*) "Final SI= ", SI
         SI=SI-1.D0
         DLWDB=DLWDB-LIST(1,2)
         CALL FREEM(IP_WLIST)
         CALL FREEM(IP_RIJMAT)
         CALL TIHALT("STARDIAGSC",ISUB)
         
         RETURN
      END
                
      
      SUBROUTINE STARDIAG(LSTE,NEL,NLIST,LIST,ILMAX,I_P,SI,DBETA,DLWDB)
         USE HElement
         IMPLICIT NONE
         INTEGER NEL,I_P
         INTEGER LSTE(NEL,NLIST),NLIST,ILMAX
         REAL*8 LIST(ILMAX,0:2)
         REAL*8 RIJMAT(*),WLIST(*)
         POINTER (IP_RIJMAT,RIJMAT),(IP_WLIST,WLIST),(IP_WORK,WORK)
         INTEGER ISUB
         INTEGER WORKL,WORK(*),INFO
         REAL*8 SI,DLWDB,DBETA,OD
         INTEGER I,J
         TYPE(HElement) RR
         INCLUDE 'vmc.inc'
         
         IF(HElementSize.GT.1) STOP "STARDIAG cannot function with complex orbitals."

         CALL TISET('STARDIAG  ',ISUB)
         CALL MEMORY(IP_RIJMAT,NLIST*NLIST,"RIJMAT")
         CALL MEMORY(IP_WLIST,NLIST,"WLIST")
         WORKL=3*NLIST
         CALL MEMORY(IP_WORK,WORKL,"WORK")

         CALL SORT3RN(NLIST-1,LIST(2,0),LIST(2,1),LIST(2,2),HElementSize)

         CALL AZZERO(RIJMAT,NLIST*NLIST)
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
         CALL FREEM(IP_WORK)
         WRITE(6,*)
         WRITE(6,*) "Highest root:",WLIST(NLIST)
!.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
         SI=0.D0
        
!         WRITE(67,*) "Eigenvalues are: "
!         DO I=1,NLIST
!            WRITE(67,"F12.6,$") WLIST(I)
!         ENDDO

!Divide through by largest eigenvalue to prevent blowing up in some cases
         IF(TCALCREALPROD) THEN
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
!            RR=1.d0
!            IF(I.LT.NLIST) RR=HElement(WLIST(I)-LIST(I+1,0))
!            write(6,"(I,3G)") I-1,WLIST(I),LIST(I+1,0),RR
!            IF(RR.AGT.1.d-10)  WRITE(6,*) 1/(RIJMAT((I-1)*NLIST+1)**2),RIJMAT((I-1)*NLIST+1)**2
!         ENDDO
         DO I=0,NLIST-1
            SI=SI+RIJMAT(I*NLIST+1)*RIJMAT(I*NLIST+1)*(WLIST(I+1)**I_P)
            IF(DBETA.NE.0.D0) THEN
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
         SI=SI-1.D0
         DLWDB=DLWDB-LIST(1,2)
         CALL FREEM(IP_WLIST)
         CALL FREEM(IP_RIJMAT)
         CALL TIHALT("STARDIAG  ",ISUB)
         RETURN
      END


!.. Use an iterative Order(N) root-finding method to diagonalize the
!.. star matrix.
!  LIST(0:NLIST-1,:) contains data
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ
      SUBROUTINE STARDIAG2(LSTE,NEL,NLIST,LIST,ILMAX,BETA,I_P,SI,DBETA,DLWDB,NROOTS,iLogging)
         USE HElement
         IMPLICIT NONE
         INTEGER NEL,I_P
         INTEGER LSTE(NEL,NLIST),NLIST,ILMAX
         TYPE(HElement) LIST(0:ILMAX-1,0:2)
         INTEGER ISUB
         REAL*8 SI,DLWDB,DBETA,NORM,E0,BETA,osi
         TYPE(HElement) DLWDB2,RR
         INTEGER I,J,NROOTS
         REAL*8 ROOTS(0:NROOTS),RPN,R
         INTEGER iLogging
         INTEGER iEigv,iDegen
         LOGICAL lWarned
         REAL*8 NORMCHECK,NORMROOTS
         include 'vmc.inc'
         include 'uhfdet.inc'
         CALL TISET('STARDIAG2 ',ISUB)
!.. we need to sort A and B (and the list of hamil values) into ascending A order
!         WRITE(6,*) (LIST(I,2),I=1,NLIST)
!         WRITE(6,*) (LIST(I,1),I=1,NLIST)
         CALL SORT3RN(NLIST-1,LIST(1,0),LIST(1,1),LIST(1,2),HElementSize)
!         WRITE(6,*) (LIST(I,2),I=1,NLIST)
!         WRITE(6,*) (LIST(I,1),I=1,NLIST)

!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ
         if(nRoots.eq.nList) then
!  nRoots is one more than the number of excitations
!we've been asked to search to see how many roots we need for convergence.
            i=nRoots-1
            do while (i.gt.0.and.abs((nRoots-i)*(List(i,0)%v**I_P)).ge.STARCONV)
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
            do while (i.gt.0.and.abs((nRoots-i)*(List(i,0)%v**I_P)/si).ge.1.d-2)
               Si=SI+List(i,0)%v**I_P
               i=i-1
            enddo
            nRoots=nRoots-1-i
            write(6,*) nRoots+1, " needed for method 2 convergence 1.d-3."
            nRoots=nRoots+1
         endif

!.. Find the eigenvalues
!  NLIST is the number of elements in the list, but we need to give the index of the last element, NLIST-1
!         CALL PLOTROOTSTAR(NLIST-1,LIST(0,0),LIST(0,1),ROOTS,NROOTS) 
         CALL FINDROOTSTAR(NLIST-1,LIST(0,0),LIST(0,1),ROOTS,NROOTS)
         SI=0.D0
         WRITE(6,*) "Highest root:",ROOTS(NROOTS)
         E0=List(0,2)%v
         lWarned=.false.
!.. divide through by highest eigenvalue to stop things blowing up
!         DO I=1,NROOTS
!            WRITE(6,*) ROOTS(I)
!,list(NLIST-NROOTS+I-1,0)%v
!         ENDDO
         iEigv=0 
         iDegen=0
!NLIST is the length of list, and the max possible value of  NROOTS
         NORMCHECK=0
         DO I=NROOTS,1,-1
            iDegen=iDegen+1
            RR=1.d0
            IF(I.LT.NROOTS) RR=HElement(ROOTS(I))-LIST(NLIST-NROOTS+I,0)
!            write(6,"(I,3G)") I,ROOTS(I),LIST(NLIST-NROOTS+I,0),RR
            IF(ROOTS(I).EQ.LIST(NLIST-NROOTS+I-1,0)%v.OR..NOT.(RR.AGT.1.d-13)) THEN
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
               NORM=1.D0
               if(iEigv.ge.1.and.iEigv.le.3) then
!!!                  write(6,*) iDegen-1
                  iDegen=1
               endif
!List(0,:) is the HF det.  We set its value in the eigenvector to 1.  The remaining NLIST-1 items are used to normalize.
               DO J=1,NLIST-1
                  RR=HElement(ROOTS(I))-LIST(J,0)
                  IF(.NOT.(RR.AGT.1e-13)) THEN
!see comment below
                     WRITE(6,"(A,I,A,G,A,G,A,I)") "WARNING: Eigenvalue I=",I,":",ROOTS(I), " dangerously close to rhojj=",LIST(J,0)," J=",J
                     WRITE(6,"(A,I,2G)") "POLE,NUMER",J,LIST(J,0),LIST(J,1)
                     WRITE(6,"(A,I,2G)") "POLE,NUMER",J-1,LIST(J-1,0),LIST(J-1,1)
                  ENDIF
                  NORM=NORM+SQ(LIST(J,1)/RR)
               ENDDO
!.. We add in the first element of the eigenvector * lambda**P
! As a test for size consistency - use expansion of exp(rho)-1 instead of normal rho matrix - c.f QCISD
               IF(TQUADRHO) THEN
                   !In order to stop it blowing up, divide through by r_max + r_max^2/2
                   NORMROOTS=(ROOTS(I)+(ROOTS(I)**2)/2)/(ROOTS(NROOTS)+(ROOTS(NROOTS)**2)/2)
                   RPN=(NORMROOTS**I_P)*1.D0/NORM
               ELSEIF(TEXPRHO) THEN
                   NORMROOTS=(EXP(ROOTS(I))-1)/(EXP(ROOTS(NROOTS))-1)
                   RPN=(NORMROOTS**I_P)*1.D0/NORM
               ELSE
                   RPN=(ROOTS(I)**I_P)*1.D0/NORM
               ENDIF
               
!               write(6,*) NORM,RPN
               IF(.not.lWarned.and.RPN/SI.LT.1.d-4) then
                  lWarned=.true.
!!!                  WRITE(6,*) "Root ",NROOTS-I," has low contribution."
!!!                  WRITE(6,*) "SI=",SI
               ENDIF
               SI=SI+RPN
               NORMCHECK=NORMCHECK+1/NORM
!               WRITE(6,*) I,RPN,1/NORM
               IF(iEigv.le.2) then
!!!                  write(6,"(A,I,A,2G,$)") "Eigenvalue ",iEigv," = ",roots(i),E0-(i_P/Beta)*log(roots(i))
!                  write(6,*) "***",E0
               endif
               IF(DBETA.NE.0.D0) THEN
                  DLWDB2=LIST(0,2)
!                  WRITE(6,*) LIST(1,2),SQRT(1/NORM)
                  DO J=1,NLIST
                     DLWDB2=DLWDB2+LIST(J,2)*(DCONJG(LIST(J,1))/(HElement(ROOTS(I))-LIST(J,0)))
!                WRITE(6,*) LIST(J,2),
!     &            LIST(J,1)/((ROOTS(I+1)-LIST(J,0))*SQRT(NORM))
                  ENDDO
                  DLWDB=DLWDB+DLWDB2%v*RPN
               ENDIF
!               WRITE(6,*) ROOTS(I+1)**I_P,DLWDB2*(ROOTS(I+1)**I_P)/NORM
!               WRITE(6,*)
            ENDIF
         ENDDO
               if(iEigv.gt.1.and.iEigv.le.3) then
                  write(6,*) iDegen-1
                  iDegen=1
               endif
         IF(TQUADRHO) WRITE(6,*) "QUADRATIC EXPANSION OF RHO MATRIX USED"
         IF(TEXPRHO) WRITE(6,*) "EXPONENTIAL EXPANSION OF RHO MATRIX USED"
         write(6,*)
         WRITE(6,*) "Final SI=",SI
         WRITE(6,*) "Norm of i projection:", NORMCHECK
         IF(ABS(NORMCHECK-1).gt.0.01) WRITE(6,*)  "WARNING: Norm differs from 1 by more than 0.01.  Convergence may not be reached."
         SI=SI-1.D0
         DLWDB=DLWDB-LIST(0,2)%v
         CALL TIHALT("STARDIAG2 ",ISUB)
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
      SUBROUTINE StarAddSingles(nI,nJ,ExcitInfo,iExcit,iMaxExcit,rhii,rhoeps,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,nTay,ECore)
         USE HElement      
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         Type(BasisFN) G1(*)
         INTEGER nI(nEl),nEl,i_P,nBasisMax(*),Brr(nBasis),nBasis,nMsh
         INTEGER nMax,nTay(2),L,LT,nWHTay,iLogging
         COMPLEX*16 fck(*)
         TYPE(HElement) UMat(*)
         REAL*8 Beta, ALat(3),RhoEps,ECore,dBeta
         TYPE(HElement)  ExcitInfo(0:iMaxExcit,0:2)
!.. New lists are generated here
!.. This will contain all the info needed to work out the value of the
!.. star
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ

         INTEGER nJ(nEl),iExcit,iEx(2,2),nK(nEl)
         LOGICAL tDummy
         type(HElement) StarMat(5,5),rh,rhii
         integer i,iExc,iMaxExcit

!Needed for diagonalizer
         REAL*8 WLIST(5),WORK(3*5)         
         TYPE(HElement) NWORK(4*5)
         INTEGER INFO

         call azzero(StarMat,25*HElementSize)
         iEx(1,1)=2
!.. Get the orbitals which are excited in going from I to J
!.. IEX(1,*) are in I, and IEX(2,*) are in J
         CALL GetExcitation(nI,nJ,nEl,iEx,tDummy)
         IF(iEx(1,1).GT.0.AND.iEx(1,2).GT.0) THEn
!  We've got a double excitation
            StarMat(1,1)=ExcitInfo(iExcit,0)
            iExc=1
!  Now generate all possible singles between nI and nJ, and put them into StarMat
!(i,a)
            CALL ICOPY(nEl,nI,1,nK,1)
      lp0:  DO i=1,nEl
               IF(nK(i).EQ.iEx(1,1)) THEN
                  nK(i)=iEx(2,1)
                  CALL SORTI(nEl,nK)
                  exit lp0
               endif
            end do lp0
            CALL CalcRho2(nJ,nK,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,1,ECore)
!            call writedet(6,nk,nel,.false.)
!            write(6,*) rh
            IF(rh.age.rhoeps) then
               iExc=iExc+1
               StarMat(1,iExc)=rh/rhii
               CALL CalcRho2(nK,nK,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
               StarMat(iExc,iExc)=rh/rhii
            ENDIF
!(i,b)
            CALL ICOPY(nEl,nI,1,nK,1)
        lp1:DO i=1,nEl
               IF(nK(i).EQ.iEx(1,1)) THEN
                  nK(i)=iEx(2,2)
                  CALL SORTI(nEl,nK)
                  exit lp1
               endif
            end do lp1
            CALL CalcRho2(nJ,nK,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,1,ECore)
!            call writedet(6,nk,nel,.false.)
!            write(6,*) rh
            IF(rh.age.rhoeps) then
               iExc=iExc+1
               StarMat(1,iExc)=rh/rhii
               CALL CalcRho2(nK,nK,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
               StarMat(iExc,iExc)=rh/rhii
            ENDIF
!(j,a)
            CALL ICOPY(nEl,nI,1,nK,1)
       lp2: DO i=1,nEl
               IF(nK(i).EQ.iEx(1,2)) THEN
                  nK(i)=iEx(2,1)
                  CALL SORTI(nEl,nK)
                  exit lp2
               endif
            end do lp2
            CALL CalcRho2(nJ,nK,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,1,ECore)
!            call writedet(6,nk,nel,.false.)
!            write(6,*) rh
            IF(rh.age.rhoeps) then
               iExc=iExc+1
               StarMat(1,iExc)=rh/rhii
               CALL CalcRho2(nK,nK,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
               StarMat(iExc,iExc)=rh/rhii
            ENDIF
!(j,b)
            CALL ICOPY(nEl,nI,1,nK,1)
   lp3:     DO i=1,nEl
               IF(nK(i).EQ.iEx(1,2)) THEN
                  nK(i)=iEx(2,2)
                  CALL SORTI(nEl,nK)
                  exit lp3
               endif
            end do lp3
            CALL CalcRho2(nJ,nK,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,1,ECore)
!            call writedet(6,nk,nel,.false.)
!            write(6,*) rh
            IF(rh.age.rhoeps) then
               iExc=iExc+1
               StarMat(1,iExc)=rh/rhii
               CALL CalcRho2(nK,nK,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,0,ECore)
               StarMat(iExc,iExc)=rh/rhii
            ENDIF

! Now diagonalize.
            IF(HElementSize.EQ.1) THEN
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
               rh=sqrt(dreal(StarMat(1,i)*DCONJG(StarMat(1,I))))
               ExcitInfo(iExcit+i-1,1)=ExcitInfo(iExcit,1)*rh
               ExcitInfo(iExcit+i-1,2)=ExcitInfo(iExcit,2)*rh
            enddo
            iExcit=iExcit+iExc-1
         ENDIF
      END

      !Routine which takes the excitation form for the determinant, and calculates the full path.
      SUBROUTINE GETFULLPATH(nI,nEl,noexcits,Orbchange,nJ)
        IMPLICIT NONE
        INTEGER :: nI(nEl),nJ(nEl),nEl,noexcits,Orbchange(noexcits*noexcits)
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
        CALL SORTI(nEl,nJ)
      END
      
     !Routine which takes a root determinant (nI), and a double excitation (nJ), and calculates the orbitals which have been excited.
     !This information is put into Orbchange(4), with the first two values being the excited from orbitals (ij), and the second two being the excited to orbitals (ab).
     SUBROUTINE GETEXCITSCHANGE(nI,nJ,nEl,Orbchange)
        IMPLICIT NONE
        INTEGER :: nI(nEl),nJ(nEl),nEl,Orbchange(4,1),q,I,J
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
                    Orbchange(1,1)=nI(I)
                    q=2
                ELSEIF(q.eq.2) THEN
                    Orbchange(2,1)=nI(I)
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
                    Orbchange(3,1)=nJ(I)
                    q=2
                ELSEIF(q.eq.2) THEN
                    Orbchange(4,1)=nJ(I)
                    RETURN
                ELSE
                    STOP 'ERROR IN GETEXCITSCHANGE'
                ENDIF
            ENDIF
        ENDDO

        STOP 'ERROR IN GETEXCITSCHANGE'

        END
        
