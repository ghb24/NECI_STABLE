!  A function to generate all possible excitations of a determinant and link them together in a star
!    The rho_jj, rho_ij and H_ij values are stored for each connected determinant.
!   Based on a combined FMCPR3STAR and FMCPR3STAR2, this instead generates excitations on the fly.
   FUNCTION fMCPR3StarNewExcit(nI,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,TMat,nMax,ALat,UMat,nTay, &
               RhoEps, L, LT,nWHTay, iLogging, tSym, ECore,dBeta,dLWdB)
         USE HElement      
         IMPLICIT NONE
         INCLUDE 'basis.inc'
         Type(BasisFN) G1(*)
         INTEGER nI(nEl),nEl,i_P,nBasisMax(*),Brr(nBasis),nBasis,nMsh
         INTEGER nMax,nTay(2),L,LT,nWHTay,iLogging
         COMPLEX*16 fck(*)
         TYPE(HElement) TMat(*),UMat(*)
         REAL*8 Beta, ALat(3),RhoEps,ECore,dBeta
         TYPE(HDElement) dLWdB
         TYPE(HDElement) fMCPR3StarNewExcit
         TYPE(HElement), Allocatable :: ExcitInfo(:,:)
         LOGICAL tSym
!.. New lists are generated here
!.. This will contain all the info needed to work out the value of the
!.. star
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ

         INTEGER exFlag
         INTEGER, allocatable :: nExcit(:)
         INTEGER nExcitMemLen,nStore(6)
         INTEGER nJ(nEl),iExcit,iMaxExcit
         INTEGER iErr
         INTEGER nRoots,i
         TYPE(HElement) rh,rhii,EHFDiff
         
   
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
!.. Count the excitations.
         nStore(1)=0
         CALL GenSymExcitIt2(nI,nEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
         Allocate(nExcit(nExcitMemLen))
         nExcit(1)=0
         CALL GenSymExcitIt2(nI,nEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)
!.. iC now contains the number of excitations.
!.. Allocate memory for the lists
         Allocate(ExcitInfo(0:iMaxExcit,0:2),stat=iErr)
         CALL MemAlloc(iErr,ExcitInfo,(iMaxExcit+1)*3*HElementSize,"ExcitInfo")
!.. This will contain all the info needed to work out the value of the
!.. star
!.. LIST(0,...) corresponds to J=I
!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ
         i=0
         ExcitInfo(i,0)=1.D0
         ExcitInfo(i,1)=1.D0
         ExcitInfo(i,2)=GetHElement2(nI,nI,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,TMat,nMax,ALat,UMat,00,ECore)
         if(BTEST(nWhTay,5)) then
! We use the Zeroth order N-particle Hartree-Fock hamiltonian (as MP theory), but shifted by E_HF-E0.

!nMax has Arr hidden in it
            IF(nTay(2).eq.5) THEN
               call GetH0ElementDCCorr(nI,nI,nEl,nBasisMax,G1,nBasis,Brr,NMSH,FCK,TMat,nMax,ALAT,UMat,ECore,rhii)
            ELSE
               call GetH0Element(nI,nEl,nMax,nBasis,rhii)
            ENDIF
            EHFDiff=ExcitInfo(i,2)-rhii
         endif
         CALL CalcRho2(nI,nI,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,TMat,nMax,ALat,UMat,rhii,nTay,0,ECore)
!         write(75,*) rhii
    lp:  do while(.true.)
            CALL GenSymExcitIt2(nI,nEl,G1,nBasis,nBasisMax,.false.,nExcit,nJ,iExcit,0,nStore,exFlag)
            IF(nJ(1).eq.0) exit lp
            CALL CalcRho2(nI,nJ,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,TMat,nMax,ALat,UMat,rh,nTay,-1,ECore)
            if(rh .agt. RhoEps) then
               i=i+1
               ExcitInfo(i,1)=rh/rhii
               if(btest(nwhtay,5)) then
                  call GetH0Element(nJ,nEl,nMax,nBasis,rh)
                  rh=rh+EHFDiff
                  rh=rh*HElement(-Beta/I_P)
                  rh=exp(rh)
               else
                  CALL CalcRho2(nJ,nJ,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,TMat,nMax,ALat,UMat,rh,nTay,0,ECore)
               endif
               ExcitInfo(i,0)=rh/rhii
               ExcitInfo(i,2)=GetHElement2(nI,nJ,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,TMat,nMax,ALat,UMat,-1,ECore)
!               write(75,*) rh,rh/rhii
            endif
         enddo lp
!Tell MCPATHS how many excitations there were and how many we are keeping
         L=i
         LT=iMaxExcit
         iExcit=i
         Deallocate(nExcit)
!.. we now have a list length NLCUR of dets in the star.
!.. Call a routine to generate the value of the star
         WRITE(6,*) iExcit," excited determinants in star"
         IF(.NOT.BTEST(NWHTAY,0)) THEN
            WRITE(6,*) "Beginning Complete Star Diagonalization"
            CALL StarDiag(0,nEl,iExcit+1,ExcitInfo,iMaxExcit+1,i_P,fMCPR3StarNewExcit,dBeta,dLWdB)
         ELSE
            WRITE(6,*) "Beginning Polynomial Star Diagonalization"
            nRoots=iExcit
            IF(BTEST(NWHTAY,1)) THEN
               nRoots=1
               WRITE(6,*) "Searching for 1 root"
            ELSEIF(BTEST(NWHTAY,2)) THEN
               WRITE(6,*) "Searching for enough roots to converge sum"
               nRoots=iExcit+1
            ELSE
               WRITE(6,*) "Searching for all roots"
            ENDIF
            CALL StarDiag2(0,nEl,iExcit+1,ExcitInfo,iMaxExcit+1,i_P,fMCPR3StarNewExcit,dBeta,dLWdB,nRoots)
         ENDIF
         call MemDealloc(ExcitInfo)
         Deallocate(ExcitInfo)
      END
!.. This version creates a star of all 2-vertex terms.
!.. This sets up the excitation generators and the memory - using the old excitation generators
FUNCTION FMCPR3STAR(NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY, &
                     RHOEPS,LSTE,ICE,LIST,L,LT,NWHTAY,ILOGGING,TSYM,ECORE,ILMAX,DBETA,DLWDB)
         USE HElement      
         IMPLICIT NONE
         INTEGER I_V,NEL,I_P,NBASISMAX(*),G1(*),NBASIS,BRR(*),NMSH,NMAX
         INTEGER NTAY,NWHTAY,ILOGGING,LT
         REAL*8 FCK(*), TMat(*),ALAT(*),UMAT(*),ECORE
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

         FMCPR3STAR=FMCPR3STAR2(NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY, &
           RHOEPS,NLSTE,NICE,NLIST,L,LT,NWHTAY,ILOGGING,TSYM,ECORE,NLENLIST,DBETA,DLWDB)

         CALL FREEM(IP_NLIST)
         CALL FREEM(IP_NLSTE)
         CALL FREEM(IP_NICE)
         RETURN
      END
!.. This version creates a star of all 2-vertex terms.
!..   This is the heart of the function, called once the excitations are found.
      FUNCTION FMCPR3STAR2(NI,BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,NTAY, &
           RHOEPS,LSTE,ICE,LIST,L,LT,NWHTAY,ILOGGING,TSYM,ECORE,ILMAX,DBETA,DLWDB)
         USE HElement     
         IMPLICIT NONE
         TYPE(HDElement) FMCPR3Star2
         INCLUDE 'basis.inc'
         TYPE(BasisFN) G1(*)
         INTEGER I_V,NEL,I_P,NBASISMAX(*),NBASIS,BRR(*),NMSH,NMAX
         INTEGER NTAY,NWHTAY,ILOGGING,LT
         REAL*8 ALAT(*),ECORE
         TYPE(HElement) TMat(*),UMat(*)
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
               CALL CALCRHO2(NI,LSTE(1,I),BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT, &
                  RH,NTAY,ICE(I),ECORE)
               IF(RH .AGT. RHOEPS) THEN
                  IF(NLCUR.NE.I) CALL ICOPY(NEL,LSTE(1,I),1,LSTE(1,NLCUR),1)
                  ICE(NLCUR)=ICE(I)
                  IF(NLCUR.EQ.0) RHII=RH
                  LIST(NLCUR,1)=RH/RHII
                  CALL CALCRHO2(LSTE(1,I),LSTE(1,I),BETA,I_P,NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT, &
                    RH,NTAY,0,ECORE)
                  LIST(NLCUR,0)=RH/RHII
                  LIST(NLCUR,2)=GETHELEMENT2(NI,LSTE(1,I),NEL,NBASISMAX,G1,NBASIS,BRR,NMSH,FCK,TMat,NMAX,ALAT,UMAT,-1,ECORE)
!                  CALL WRITEDET(6,LSTE(1,I),NEL,.FALSE.)
!                  WRITE(6,*) (LIST(NLCUR,J),J=0,2)
                  NLCUR=NLCUR+1
               ENDIF
            ENDIF
         ENDDO
         L=NLCUR
         LT=ILMAX+1
!.. we now have a list length NLCUR of dets in the star.
!.. Call a routine to generate the value of the star
         WRITE(6,*) NLCUR," determinants in star"
         IF(.NOT.BTEST(NWHTAY,0)) THEN
            WRITE(6,*) "Beginning Complete Star Diagonalization"
            CALL STARDIAG(LSTE,NEL,NLCUR,LIST,ILMAX+1,I_P,FMCPR3STAR2,DBETA,DLWDB)
         ELSE
            WRITE(6,*) "Beginning Polynomial Star Diagonalization"
            NROOTS=NLCUR
            IF(BTEST(NWHTAY,1)) NROOTS=1
            CALL STARDIAG2(LSTE,NEL,NLCUR,LIST,ILMAX+1,I_P,FMCPR3STAR2,DBETA,DLWDB,NROOTS)
         ENDIF
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
         IF(HElementSize.GT.1) STOP "STARDIAG cannot function with complex orbitals."

         CALL TISET('STARDIAG  ',ISUB)
         CALL MEMORY(IP_RIJMAT,NLIST*NLIST,"RIJMAT")
         CALL MEMORY(IP_WLIST,NLIST,"WLIST")
         WORKL=3*NLIST
         CALL MEMORY(IP_WORK,WORKL,"WORK")


         CALL AZZERO(RIJMAT,NLIST*NLIST)
!.. Now we fill the RIJ array
         DO I=0,NLIST-1
            RIJMAT(I*NLIST+I+1)=LIST(I+1,0)
            RIJMAT(I+1)=LIST(I+1,1)
         ENDDO

!.. Diagonalize
         CALL DSYEV('V','L',NLIST,RIJMAT,NLIST,WLIST,WORK,WORKL,INFO)
         IF(INFO.NE.0) THEN
            WRITE(6,*) 'DYSEV error: ',INFO
            STOP
         ENDIF
         CALL FREEM(IP_WORK)
!.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
         SI=0.D0
         DO I=0,NLIST-1
!            WRITE(6,*) WLIST(I+1)
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
         SI=SI-1.D0
         DLWDB=DLWDB-LIST(1,2)
         CALL FREEM(IP_WLIST)
         CALL FREEM(IP_RIJMAT)
         CALL TIHALT("STARDIAG  ",ISUB)
         RETURN
      END


!.. Use an iterative Order(N) root-finding method to diagonalize the
!.. star matrix.
      SUBROUTINE STARDIAG2(LSTE,NEL,NLIST,LIST,ILMAX,I_P,SI,DBETA,DLWDB,NROOTS)
         USE HElement
         IMPLICIT NONE
         INTEGER NEL,I_P
         INTEGER LSTE(NEL,NLIST),NLIST,ILMAX
         TYPE(HElement) LIST(ILMAX,0:2)
         INTEGER ISUB
         REAL*8 SI,DLWDB,DBETA,NORM
         TYPE(HElement) DLWDB2
         INTEGER I,J,NROOTS
         REAL*8 ROOTS(1:NROOTS+1),RPN,R
         CALL TISET('STARDIAG2 ',ISUB)
      
!.. we need to sort A and B (and the list of hamil values) into ascending A order
!         WRITE(6,*) (LIST(I,2),I=1,NLIST)
!         WRITE(6,*) (LIST(I,1),I=1,NLIST)
         CALL SORT3R(NLIST-1,LIST(2,0),LIST(2,1),LIST(2,2))
!         WRITE(6,*) (LIST(I,2),I=1,NLIST)
!         WRITE(6,*) (LIST(I,1),I=1,NLIST)

!.. LIST(J,0) = RHOJJ
!.. LIST(J,1) = RHOIJ
!.. LIST(J,2) = HIJ
         if(nRoots.eq.nList) then
!  nRoots is one more than the number of excitations
!we've been asked to search to see how many roots we need for convergence.
            i=nRoots-1
            do while (i.gt.0.and.((nRoots-i)*(List(i,0)%v**I_P)).ge.1.d-5)
               i=i-1
            enddo
            nRoots=nRoots-1-i
            write(6,*) nRoots+1, " needed for convergence 1.d-5."
         endif

!.. Find the eigenvalues
         CALL FINDROOTSTAR(NLIST-1,LIST(1,0),LIST(1,1),ROOTS,NROOTS)
         SI=0.D0
         WRITE(6,*) "Highest root:",ROOTS(NROOTS+1)
!.. divide through by highest eigenvalue to stop things blowing up
         DO I=NROOTS,0,-1
            IF(ROOTS(I+1).EQ.LIST(NLIST-NROOTS+I,0)%v) THEN
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
               NORM=1.D0
               DO J=2,NLIST
                  NORM=NORM+SQ(LIST(J,1)/(HElement(ROOTS(I+1))-LIST(J,0)))
               ENDDO
!.. We add in the first element of the eigenvector * lambda**P
               RPN=(ROOTS(I+1)**I_P)*1.D0/NORM
               SI=SI+RPN
               IF(DBETA.NE.0.D0) THEN
                  DLWDB2=LIST(1,2)
!                  WRITE(6,*) LIST(1,2),SQRT(1/NORM)
                  DO J=2,NLIST
                     DLWDB2=DLWDB2+LIST(J,2)*(DCONJG(LIST(J,1))/(HElement(ROOTS(I+1))-LIST(J,0)))
!                WRITE(6,*) LIST(J,2),
!     &            LIST(J,1)/((ROOTS(I+1)-LIST(J,0))*SQRT(NORM))
                  ENDDO
                  DLWDB=DLWDB+DLWDB2%v*RPN
               ENDIF
!               WRITE(6,*) ROOTS(I+1)**I_P,DLWDB2*(ROOTS(I+1)**I_P)/NORM
!               WRITE(6,*)
            ENDIF
         ENDDO
         SI=SI-1.D0
         DLWDB=DLWDB-LIST(1,2)%v
         CALL TIHALT("STARDIAG2 ",ISUB)
         RETURN
      END

