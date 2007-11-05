!  A function to generate all possible excitations of a determinant and link them together in a star
!    The rho_jj, rho_ij and H_ij values are stored for each connected determinant.
!   Based on a combined FMCPR3STAR and FMCPR3STAR2, this instead generates excitations on the fly.
   FUNCTION fMCPR3StarNewExcit(nI,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,nTay, &
               RhoEps, L, LT,nWHTay, iLogging, tSym, ECore,dBeta,dLWdB,MP2E)
         USE HElement      
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
         TYPE(HElement), Allocatable :: ExcitInfo(:,:)
         TYPE(HElement) HIJS(0:2)
!         REAL*8 LARGERHOJJ(10000)
         INTEGER iPath(nEl,0:2)
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
         INTEGER nRoots,i,j
         TYPE(HElement) rh,rhii,EHFDiff

         TYPE(HDElement) MP2E         
         LOGICAL tStarSingles
         INTEGER nIExcitFormat(nEl)
!         LARGERHOJJ(:)=0.D0
         IF(tStoreAsExcitations) THEN
            nIExcitFormat(1)=-1
            nIExcitFormat(2)=0
         ELSE
            CALL ICOPY(NEL,nI,1,nIExcitFormat,1)
         ENDIF
         tStarSingles=BTEST(nWHTay,7)
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
!.. iMaxExcit now contains the number of excitations.
!.. Allocate memory for the lists
         IF(tStarSingles) THEN
            iMaxExcit=iMaxExcit*5
            WRITE(6,*) "Adding StarSingles.  (5x storage space used)."
         ENDIF
         Write(6,*) "Allocating storage for ",iMaxExcit," excitations."
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
            CALL CalcRho2(nIExcitFormat,nJ,Beta,i_P,nEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,nMax,ALat,UMat,rh,nTay,iExcit,ECore)
            if(rh .agt. RhoEps) then
               i=i+1
               ExcitInfo(i,1)=rh/rhii
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
!         DO j=1,10000
!            WRITE(55,*) LARGERHOJJ(J)
!         ENDDO
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
            ELSEIF(BTEST(NWHTAY,6)) THEN
               WRITE(6,*) "Searching for enough roots to converge sum - Method 2"
               nRoots=iExcit+2
            ELSE
               WRITE(6,*) "Searching for all roots"
            ENDIF
            CALL StarDiag2(0,nEl,iExcit+1,ExcitInfo,iMaxExcit+1,Beta,i_P,fMCPR3StarNewExcit,dBeta,dLWdB,nRoots,iLogging)
         ENDIF
         call MemDealloc(ExcitInfo)
         Deallocate(ExcitInfo)
      END
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
         WRITE(6,*)
         WRITE(6,*) "Highest root:",WLIST(NLIST)
!.. RIJMAT now contains the eigenvectors, and WLIST the eigenvalues         
         SI=0.D0
!         DO I=0,NLIST-1
!            WRITE(6,*) WLIST(I+1)
!         ENDDO
!         DO I=NLIST-1,0,-1
!            WRITE(6,*) I+1,RIJMAT(I*NLIST+1)*RIJMAT(I*NLIST+1)*(WLIST(I+1)**I_P),RIJMAT(I*NLIST+1)*RIJMAT(I*NLIST+1)
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
         include 'vmc.inc'
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
         DO I=NROOTS,1,-1
            iDegen=iDegen+1
            IF(ROOTS(I).EQ.LIST(NLIST-NROOTS+I-1,0)%v) THEN
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
                  write(6,*) iDegen-1
                  iDegen=1
               endif
!List(1,:) is the HF det.  We set its value in the eigenvector to 1.
               DO J=1,NLIST
                  RR=HElement(ROOTS(I))-LIST(J,0)
                  IF(.NOT.(RR.AGT.1e-15)) THEN
!see comment below
                     WRITE(6,*) "WARNING: Eigenvalue I=",I,":",ROOTS(I), " dangerously close to rhojj=",LIST(J,0)," J=",J
                  ENDIF
                  NORM=NORM+SQ(LIST(J,1)/RR)
               ENDDO
!.. We add in the first element of the eigenvector * lambda**P
!               write(6,*) ROOTS(i+1),NORM
               RPN=(ROOTS(I)**I_P)*1.D0/NORM
               IF(.not.lWarned.and.RPN/SI.LT.1.d-4) then
                  lWarned=.true.
                  WRITE(6,*) "Root ",NROOTS-I," has low contribution."
                  WRITE(6,*) "SI=",SI
               ENDIF
               SI=SI+RPN
!               WRITE(6,*) I,RPN,1/NORM
               IF(iEigv.le.2) then
                  write(6,"(A,I,A,2G,$)") "Eigenvalue ",iEigv," = ",roots(i),E0-(i_P/Beta)*log(roots(i))
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
         write(6,*)
         WRITE(6,*) "Final SI=",SI
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
